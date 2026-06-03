#!/usr/bin/env python3
"""Build a file-openable interactive_report.html using a JSON script block.

This is an alternate standalone format for browsers that fail to finish the
large `let D = JSON.parse("...")` startup assignment.  It embeds the report data
as inert JSON text and parses it from the DOM after DOMContentLoaded.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).parent
sys.path.insert(0, str(REPO))

from perturbscope.interactive import _copy_plotly_bundle, _html_template

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

DOMJSON_LOAD = """\
function loadReportData(){
  if(_dataReady) return _dataReady;
  const hdr=document.getElementById('hdr-summary');
  if(hdr) hdr.textContent='Loading embedded results data...';
  _dataReady=Promise.resolve().then(()=>{
    startupNote('Looking for embedded data block');
    const el=document.getElementById('interactive-data-json');
    if(!el) throw new Error('Embedded interactive data block is missing');
    startupNote('Embedded data block found: '+(el.textContent||'').length+' characters');
    const data=JSON.parse(el.textContent || '{}');
    startupNote('Embedded JSON parsed: '+((data.summary&&data.summary.n_cells)||0)+' cells, '+((data.perturbations&&data.perturbations.length)||0)+' perturbations');
    D=Object.assign(emptyDataStub(), data||{});
    if(hdr) hdr.textContent='Results data loaded';
    return D;
  }).catch(err=>{showStartupError(err);throw err;});
  return _dataReady;
}"""

DIAGNOSTIC_JS = """\
function startupNote(msg){
  try{
    const box=document.getElementById('startup-diagnostics');
    if(!box) return;
    const line=document.createElement('div');
    line.textContent=new Date().toLocaleTimeString()+' - '+msg;
    box.appendChild(line);
  }catch(e){}
}"""

DIAGNOSTIC_HTML = """\
<div id="startup-diagnostics" style="margin:0 0 16px;padding:10px 12px;border:1px solid #c0d1de;border-radius:6px;background:#fffdf3;color:#3a4a58;font:12px/1.5 monospace;">
  <div>Interactive report startup diagnostics</div>
</div>"""

TOP_BANNER_HTML = """\
<div id="debug-top-banner" style="position:relative;z-index:999999;margin:0;padding:14px 18px;background:#fff3cd;border:4px solid #b45309;color:#111827;font:16px/1.45 Arial,sans-serif;">
  <strong>DEBUG BUILD ACTIVE:</strong>
  This is the rebuilt DOM-JSON interactive report. If you can see this banner, the browser is opening the newest rebuilt file.
</div>"""


def _safe_json_text(data: object) -> str:
    text = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    return text.replace("</", "<\\/")


def build(results_dir: Path) -> None:
    results_dir = results_dir.resolve()
    json_path = results_dir / "interactive_data.json"
    out_path = results_dir / "interactive_report.html"
    debug_path = results_dir / "interactive_report_DEBUG_OPEN_THIS.html"

    if not json_path.exists():
        raise SystemExit(f"ERROR: missing {json_path}")

    print(f"Reading {json_path}")
    data = json.loads(json_path.read_text(encoding="utf-8"))
    json_text = _safe_json_text(data)

    html = _html_template()
    html = html.replace(
        "<title>PerturbScope Interactive Report</title>",
        "<title>PerturbScope Interactive Report DEBUG DOMJSON</title>",
        1,
    )
    html = html.replace("<body>", "<body>\n" + TOP_BANNER_HTML, 1)
    html = html.replace('<div class="content">', '<div class="content">\n' + DIAGNOSTIC_HTML, 1)
    marker = "<script>\nlet D = {DATA_STUB};"
    replacement = (
        '<script type="application/json" id="interactive-data-json">\n'
        f"{json_text}\n"
        "</script>\n"
        "<script>\n"
        f"{DIAGNOSTIC_JS}\n"
        "let D = {};"
    )
    if marker not in html:
        raise SystemExit("ERROR: could not find template data marker")
    html = html.replace(marker, replacement, 1)

    if FETCH_LOAD not in html:
        raise SystemExit("ERROR: could not find fetch-based loadReportData block")
    html = html.replace(FETCH_LOAD, DOMJSON_LOAD, 1)

    backup = results_dir / "interactive_report_before_domjson_fix.html"
    if out_path.exists() and not backup.exists():
        shutil.copy2(out_path, backup)
        print(f"Backed up old HTML -> {backup.name}")

    out_path.write_text(html, encoding="utf-8")
    debug_path.write_text(html, encoding="utf-8")
    _copy_plotly_bundle(results_dir)

    print(f"Wrote {out_path}")
    print(f"Wrote {debug_path}")
    print(f"Size: {out_path.stat().st_size} bytes")
    print("Mode: DOM JSON block; no fetch; no huge JSON.parse string assignment")


if __name__ == "__main__":
    targets = [Path(a) for a in sys.argv[1:]] or [
        REPO / "results_replogle_weissman_k562_essential_mini20"
    ]
    for target in targets:
        build(target)
