#!/usr/bin/env python3
"""Build the real original-layout interactive report with corrected loaders.

This uses the maintained PerturbScope HTML template, not the old saved HTML and
not the simplified static fallback.  Fixes applied:
  - embed interactive_data.json in a DOM JSON block, so file:// works;
  - replace fetch() loader with JSON.parse() from that block;
  - keep the original sidebar, tabs, and Plotly interactions;
  - keep Plotly local via plotly.min.js;
  - add visible global JS error reporting.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).parent
sys.path.insert(0, str(REPO))

from perturbscope.interactive import _copy_plotly_bundle, _html_template

RESULTS = REPO / "results_replogle_weissman_k562_essential_mini20"

FETCH_LOAD = """\
function loadReportData(){
  if(_dataReady) return _dataReady;
  const hdr=document.getElementById('hdr-summary');
  if(hdr) hdr.textContent='Loading results data...';
  _dataReady=fetch('interactive_data.json', {cache:'no-store'})
    .then(r=>{
      if(!r.ok) throw new Error('Could not load interactive_data.json: HTTP '+r.status);
      return r.json();
    })
    .then(data=>{
      D=Object.assign(emptyDataStub(), data||{});
      if(hdr) hdr.textContent='Results data loaded';
      return D;
    })
    .catch(err=>{
      showStartupError(new Error(String(err && err.message || err) + '. If you opened the HTML directly with file://, serve the results folder with a small web server, for example: python -m http.server 8899'));
      throw err;
    });
  return _dataReady;
}"""

DOM_LOAD = """\
function loadReportData(){
  if(_dataReady) return _dataReady;
  const hdr=document.getElementById('hdr-summary');
  if(hdr) hdr.textContent='Loading embedded results data...';
  _dataReady=Promise.resolve().then(()=>{
    const el=document.getElementById('interactive-data-json');
    if(!el) throw new Error('Embedded data block #interactive-data-json is missing');
    const data=JSON.parse(el.textContent || '{}');
    D=Object.assign(emptyDataStub(), data||{});
    if(hdr) hdr.textContent='Results data loaded';
    return D;
  }).catch(err=>{showStartupError(err);throw err;});
  return _dataReady;
}"""

GLOBAL_ERROR = """\
window.addEventListener('error', ev=>{
  try{showStartupError(new Error((ev.message||'JavaScript error')+'\\n'+(ev.filename||'')+':'+(ev.lineno||'')+':'+(ev.colno||'')));}catch(e){}
});
window.addEventListener('unhandledrejection', ev=>{
  try{showStartupError(ev.reason instanceof Error ? ev.reason : new Error(String(ev.reason||'Unhandled promise rejection')));}catch(e){}
});
"""

HOME_SHELL = """\
/* HOME */
function renderHomeShell(){
  const s = D.summary || {n_cells:0,n_genes:0,n_perturbations:0};
  document.getElementById('hdr-summary').textContent =
    s.n_cells.toLocaleString()+' cells \\u00b7 '+s.n_genes.toLocaleString()+' genes \\u00b7 '+s.n_perturbations.toLocaleString()+' perturbations';
  const cards = [{l:'Cells',v:s.n_cells.toLocaleString(),c:''},{l:'Genes',v:s.n_genes.toLocaleString(),c:''},
    {l:'Perturbations',v:s.n_perturbations.toLocaleString(),c:''}];
  document.getElementById('home-cards').innerHTML = cards.map(c=>
    '<div class="card"><div class="label">'+c.l+'</div><div class="value '+c.c+'">'+c.v+'</div></div>').join('');
  document.getElementById('home-analysis-text').textContent =
    'This dataset contains '+s.n_cells.toLocaleString()+' single cells with '+s.n_perturbations.toLocaleString()+
    ' perturbations and '+s.n_genes.toLocaleString()+' measured genes. '+
    (D.deg_summary.length?'DEG analysis was performed on '+D.deg_summary.length+' perturbations. ':'')+
    'Use the tabs above to explore QC, UMAP, heatmaps, per-perturbation DEGs, and gene expression.';
  const ALL_STEPS=['qc','preprocess','eda','score','effects','trajectory','programs','interaction','state_enrich','deg','genenet','regulatory','report','bundle'];
  const done = new Set(D.config._completed_steps||[]);
  document.getElementById('home-steps').innerHTML = ALL_STEPS.map(st=>
    '<span class="step-chip '+(done.has(st)?'done':'')+'">'+( done.has(st)?'\\u2713 ':'')+st+'</span>').join('');
  ['home-top-perts-chart','home-deg-bar','home-effect-scatter'].forEach(id=>{
    const el=document.getElementById(id);
    if(el) el.innerHTML='<p style="padding:18px;color:var(--muted)">Loading interactive chart...</p>';
  });
  const gl=document.getElementById('umap-gene-list');
  if(gl) gl.innerHTML=(D.umap_genes||D.genes||[]).map(g=>'<option value="'+g+'"></option>').join('');
  renderStaticUmap();
}

function renderHomePlotlyCharts(){
  setTimeout(()=>{try{renderHomeTopPertChart();}catch(e){showStartupError(e);}}, 50);
  setTimeout(()=>{try{renderHomeDegChart();}catch(e){showStartupError(e);}}, 120);
  setTimeout(()=>{try{renderHomeEffectScatter();}catch(e){showStartupError(e);}}, 190);
}

function renderHome(){
  renderHomeShell();
  renderHomePlotlyCharts();
}
"""

TEMPLATE_HOME_BLOCK = """\
/* HOME */
function renderHome(){
  const s = D.summary;
  document.getElementById('hdr-summary').textContent =
    s.n_cells.toLocaleString()+' cells \\u00b7 '+s.n_genes.toLocaleString()+' genes \\u00b7 '+s.n_perturbations.toLocaleString()+' perturbations';
  const cards = [{l:'Cells',v:s.n_cells.toLocaleString(),c:''},{l:'Genes',v:s.n_genes.toLocaleString(),c:''},
    {l:'Perturbations',v:s.n_perturbations.toLocaleString(),c:''}];
  document.getElementById('home-cards').innerHTML = cards.map(c=>
    '<div class="card"><div class="label">'+c.l+'</div><div class="value '+c.c+'">'+c.v+'</div></div>').join('');
  document.getElementById('home-analysis-text').textContent =
    'This dataset contains '+s.n_cells.toLocaleString()+' single cells with '+s.n_perturbations.toLocaleString()+
    ' perturbations and '+s.n_genes.toLocaleString()+' measured genes. '+
    (D.deg_summary.length?'DEG analysis was performed on '+D.deg_summary.length+' perturbations. ':'')+
    'Use the tabs above to explore QC, UMAP, heatmaps, per-perturbation DEGs, and gene expression.';
  const ALL_STEPS=['qc','preprocess','eda','score','effects','trajectory','programs','interaction','state_enrich','deg','genenet','regulatory','report','bundle'];
  const done = new Set(D.config._completed_steps||[]);
  document.getElementById('home-steps').innerHTML = ALL_STEPS.map(st=>
    '<span class="step-chip '+(done.has(st)?'done':'')+'">'+( done.has(st)?'\\u2713 ':'')+st+'</span>').join('');
  ['home-top-perts-chart','home-deg-bar','home-effect-scatter'].forEach(id=>{
    const el=document.getElementById(id);
    if(el) el.innerHTML='<p style="padding:18px;color:var(--muted)">Loading chart…</p>';
  });
  const gl=document.getElementById('umap-gene-list');
  if(gl) gl.innerHTML=(D.umap_genes||D.genes||[]).map(g=>'<option value="'+g+'"></option>').join('');
  renderStaticUmap();
  setTimeout(renderHomeTopPertChart, 50);
  setTimeout(renderHomeDegChart, 120);
  setTimeout(renderHomeEffectScatter, 190);
}
"""


def safe_json(data: object) -> str:
    return json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def main() -> int:
    results = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else RESULTS
    data_path = results / "interactive_data.json"
    data = json.loads(data_path.read_text(encoding="utf-8"))

    html = _html_template()
    marker = "<script>\nlet D = {DATA_STUB};"
    if marker not in html:
        raise SystemExit("ERROR: template data marker not found")
    html = html.replace(
        marker,
        '<script type="application/json" id="interactive-data-json">\n'
        + safe_json(data)
        + '\n</script>\n<script>\nlet D = {};',
        1,
    )

    if FETCH_LOAD not in html:
        raise SystemExit("ERROR: fetch loadReportData block not found")
    html = html.replace(FETCH_LOAD, DOM_LOAD, 1)
    html = html.replace(
        "function showTab(name, btnEl, src) {",
        GLOBAL_ERROR + "\nfunction showTab(name, btnEl, src) {",
        1,
    )
    html = html.replace(
        "  if (n==='home')       renderHome();",
        "  if (n==='home')       ensurePlotly(renderHome);",
        1,
    )
    html = html.replace(
        "      requestAnimationFrame(()=>{try{renderHome();}catch(e){showStartupError(e);}});",
        "      requestAnimationFrame(()=>{try{ensurePlotly(renderHome);}catch(e){showStartupError(e);}});",
        1,
    )
    if TEMPLATE_HOME_BLOCK not in html:
        raise SystemExit("ERROR: template renderHome block not found")
    html = html.replace(TEMPLATE_HOME_BLOCK, HOME_SHELL, 1)
    html = html.replace(
        "      requestAnimationFrame(()=>{try{ensurePlotly(renderHome);}catch(e){showStartupError(e);}});",
        "      renderHomeShell();\n      requestAnimationFrame(()=>{try{ensurePlotly(renderHomePlotlyCharts);}catch(e){showStartupError(e);}});",
        1,
    )
    html = html.replace(
        "function ensurePlotly(cb){\n  if(typeof Plotly!=='undefined'){_plotlyLoaded=true;cb();return;}",
        "function ensurePlotly(cb){\n  if(typeof Plotly!=='undefined'){_plotlyLoaded=true;cb();return;}\n  const hdr=document.getElementById('hdr-summary');\n  if(hdr) hdr.textContent='Loading local Plotly bundle...';",
        1,
    )
    html = html.replace(
        "  s.onload=()=>{\n    _plotlyLoaded=true;",
        "  s.onload=()=>{\n    _plotlyLoaded=true;\n    const hdr=document.getElementById('hdr-summary');\n    if(hdr) hdr.textContent='Plotly loaded, rendering report...';",
        1,
    )
    html = html.replace(
        "<title>PerturbFlow-Analyzer — Interactive Report</title>",
        "<title>PerturbScope Corrected Original Interactive Report</title>",
        1,
    )

    out = results / "interactive_report.html"
    debug = results / "interactive_report_CORRECTED_ORIGINAL_OPEN_THIS.html"
    backup = results / "interactive_report_before_corrected_original.html"
    if out.exists() and not backup.exists():
        shutil.copy2(out, backup)
    out.write_text(html, encoding="utf-8")
    debug.write_text(html, encoding="utf-8")
    _copy_plotly_bundle(results)

    print(f"Wrote {out} size={out.stat().st_size}")
    print(f"Wrote {debug} size={debug.stat().st_size}")
    print("Source: maintained perturbscope.interactive._html_template()")
    print("Layout: original sidebar/tabs/interactive Plotly")
    print("Data: embedded DOM JSON, no fetch()")
    print("Plotly: local plotly.min.js")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
