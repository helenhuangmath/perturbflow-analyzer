#!/usr/bin/env python3
"""Write a robust standalone interactive report.

This intentionally avoids the large original report JavaScript bundle.  The
browser has proven it can parse the embedded JSON, so this viewer keeps the JS
small and defensive: clickable tabs, summary tables, perturbation tables, and a
plot browser over the existing `plots/` assets.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path


def _collect_plot_paths(obj, out: set[str]) -> None:
    if isinstance(obj, dict):
        for value in obj.values():
            _collect_plot_paths(value, out)
    elif isinstance(obj, list):
        for value in obj:
            _collect_plot_paths(value, out)
    elif isinstance(obj, str) and obj.startswith("plots/"):
        out.add(obj)


def _safe_json(data: object) -> str:
    return json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def build(results_dir: Path) -> None:
    results_dir = results_dir.resolve()
    data_path = results_dir / "interactive_data.json"
    out_path = results_dir / "interactive_report.html"
    debug_path = results_dir / "interactive_report_ROBUST_OPEN_THIS.html"

    data = json.loads(data_path.read_text(encoding="utf-8"))
    plot_paths: set[str] = set()
    _collect_plot_paths(data, plot_paths)
    plot_paths = {p for p in plot_paths if (results_dir / p).exists()}
    data["_robust_plot_paths"] = sorted(plot_paths)

    payload = _safe_json(data)
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>PerturbScope Robust Interactive Report</title>
  <style>
    :root {{
      --bg:#f6f8fb; --panel:#ffffff; --text:#172331; --muted:#5d6b7a;
      --border:#cbd7e3; --accent:#1f5f8b; --accent2:#1f7a4d; --warn:#9a5b00;
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; font-family:Arial,Helvetica,sans-serif; background:var(--bg); color:var(--text); }}
    header {{ padding:18px 24px; background:#12344f; color:white; }}
    header h1 {{ margin:0 0 6px; font-size:24px; }}
    header .sub {{ color:#dce8f4; font-size:14px; }}
    #status {{ margin:14px 24px 0; padding:10px 12px; border:2px solid #b45309; background:#fff8e1; color:#332400; font:13px/1.45 monospace; }}
    nav {{ display:flex; gap:6px; padding:16px 24px 0; border-bottom:1px solid var(--border); background:white; flex-wrap:wrap; }}
    nav button {{ border:0; background:transparent; border-bottom:3px solid transparent; padding:10px 14px; cursor:pointer; color:var(--muted); font-weight:700; font-size:14px; }}
    nav button.active {{ color:var(--accent); border-bottom-color:var(--accent); }}
    main {{ padding:20px 24px 40px; }}
    .tab {{ display:none; }}
    .tab.active {{ display:block; }}
    .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:12px; margin-bottom:18px; }}
    .card {{ background:var(--panel); border:1px solid var(--border); border-radius:8px; padding:14px; }}
    .label {{ color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.04em; }}
    .value {{ font-size:26px; font-weight:800; margin-top:5px; }}
    .toolbar {{ display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin:10px 0 14px; }}
    input, select {{ border:1px solid var(--border); border-radius:6px; padding:8px 10px; font-size:14px; background:white; }}
    table {{ width:100%; border-collapse:collapse; background:white; border:1px solid var(--border); }}
    th, td {{ padding:8px 10px; border-bottom:1px solid var(--border); text-align:left; font-size:13px; vertical-align:top; }}
    th {{ background:#eaf1f7; color:#33485d; position:sticky; top:0; }}
    tr:hover td {{ background:#f3f8fc; }}
    .table-wrap {{ max-height:520px; overflow:auto; border:1px solid var(--border); background:white; }}
    .plot-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(260px,1fr)); gap:14px; }}
    .plot-card {{ background:white; border:1px solid var(--border); border-radius:8px; padding:10px; }}
    .plot-card img {{ width:100%; height:auto; display:block; border:1px solid #e5edf4; background:white; }}
    .plot-card a {{ display:block; margin-bottom:8px; color:var(--accent); font-size:13px; overflow-wrap:anywhere; }}
    .note {{ color:var(--muted); font-size:13px; margin:8px 0 14px; }}
    .ok {{ color:var(--accent2); font-weight:700; }}
    .err {{ color:#b91c1c; font-weight:700; white-space:pre-wrap; }}
  </style>
</head>
<body>
  <header>
    <h1>PerturbScope Interactive Report</h1>
    <div class="sub">Robust standalone viewer generated from interactive_data.json</div>
  </header>
  <div id="status">Starting JavaScript...</div>
  <nav>
    <button data-tab="summary" class="active">Summary</button>
    <button data-tab="perturbations">Perturbations</button>
    <button data-tab="deg">DEG</button>
    <button data-tab="qc">QC</button>
    <button data-tab="plots">Plots</button>
    <button data-tab="raw">Data Keys</button>
  </nav>
  <main>
    <section id="tab-summary" class="tab active"></section>
    <section id="tab-perturbations" class="tab"></section>
    <section id="tab-deg" class="tab"></section>
    <section id="tab-qc" class="tab"></section>
    <section id="tab-plots" class="tab"></section>
    <section id="tab-raw" class="tab"></section>
  </main>
  <script type="application/json" id="interactive-data-json">
{payload}
  </script>
  <script>
  const statusEl = document.getElementById('status');
  function status(msg, cls) {{
    statusEl.className = cls || '';
    statusEl.textContent = msg;
  }}
  function esc(v) {{
    return String(v == null ? '' : v).replace(/[&<>"']/g, s => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[s]));
  }}
  function fmt(v) {{
    if (v == null) return '';
    if (typeof v === 'number') return Number.isInteger(v) ? String(v) : v.toPrecision(4);
    return String(v);
  }}
  function table(rows, columns, limit=250) {{
    rows = Array.isArray(rows) ? rows.slice(0, limit) : [];
    if (!columns) {{
      const keys = new Set();
      rows.slice(0, 50).forEach(r => Object.keys(r || {{}}).forEach(k => keys.add(k)));
      columns = Array.from(keys);
    }}
    if (!rows.length) return '<p class="note">No rows available.</p>';
    return '<div class="table-wrap"><table><thead><tr>' +
      columns.map(c => '<th>'+esc(c)+'</th>').join('') +
      '</tr></thead><tbody>' +
      rows.map(r => '<tr>'+columns.map(c => '<td>'+esc(fmt((r||{{}})[c]))+'</td>').join('')+'</tr>').join('') +
      '</tbody></table></div>' +
      (rows.length === limit ? '<p class="note">Showing first '+limit+' rows.</p>' : '');
  }}
  function keysTable(obj) {{
    return table(Object.keys(obj).sort().map(k => {{
      const v = obj[k];
      return {{key:k, type:Array.isArray(v) ? 'array' : typeof v, size:Array.isArray(v) ? v.length : (v && typeof v === 'object' ? Object.keys(v).length : '')}};
    }}), ['key','type','size'], 1000);
  }}
  function showTab(name) {{
    document.querySelectorAll('nav button').forEach(b => b.classList.toggle('active', b.dataset.tab === name));
    document.querySelectorAll('.tab').forEach(t => t.classList.toggle('active', t.id === 'tab-'+name));
  }}
  document.querySelectorAll('nav button').forEach(b => b.addEventListener('click', () => showTab(b.dataset.tab)));

  try {{
    const text = document.getElementById('interactive-data-json').textContent;
    status('Embedded data block found: '+text.length+' characters');
    const D = JSON.parse(text);
    window.D = D;
    status('Data loaded: '+D.summary.n_cells+' cells, '+D.summary.n_genes+' genes, '+D.perturbations.length+' perturbations', 'ok');

    document.getElementById('tab-summary').innerHTML =
      '<div class="grid">' +
      '<div class="card"><div class="label">Cells</div><div class="value">'+esc(D.summary.n_cells)+'</div></div>' +
      '<div class="card"><div class="label">Genes</div><div class="value">'+esc(D.summary.n_genes)+'</div></div>' +
      '<div class="card"><div class="label">Perturbations</div><div class="value">'+esc(D.summary.n_perturbations)+'</div></div>' +
      '<div class="card"><div class="label">Plot Assets</div><div class="value">'+esc((D._robust_plot_paths||[]).length)+'</div></div>' +
      '</div><h2>Cells Per Perturbation</h2>' + table(D.cells_per_pert, ['perturbation','n_cells','is_control'], 1000);

    document.getElementById('tab-perturbations').innerHTML =
      '<h2>Perturbations</h2><div class="toolbar"><input id="pert-filter" placeholder="Filter perturbations"></div><div id="pert-list"></div>';
    function renderPertList() {{
      const q = document.getElementById('pert-filter').value.toLowerCase();
      const rows = (D.cells_per_pert || []).filter(r => String(r.perturbation).toLowerCase().includes(q));
      document.getElementById('pert-list').innerHTML = table(rows, ['perturbation','n_cells','is_control'], 1000);
    }}
    document.getElementById('pert-filter').addEventListener('input', renderPertList);
    renderPertList();

    document.getElementById('tab-deg').innerHTML =
      '<h2>DEG Summary</h2>' + table(D.deg_summary, null, 500) +
      '<h2>Available DEG Tables</h2>' + keysTable(D.deg || {{}});

    document.getElementById('tab-qc').innerHTML =
      '<h2>QC Cell Sample</h2><p class="note">Showing embedded QC rows from the report data.</p>' + table(D.qc_cells, null, 300);

    const plots = D._robust_plot_paths || [];
    document.getElementById('tab-plots').innerHTML =
      '<h2>Plot Browser</h2><div class="toolbar"><input id="plot-filter" placeholder="Filter plots"></div><div id="plot-list" class="plot-grid"></div>';
    function renderPlots() {{
      const q = document.getElementById('plot-filter').value.toLowerCase();
      const subset = plots.filter(p => p.toLowerCase().includes(q)).slice(0, 120);
      document.getElementById('plot-list').innerHTML = subset.map(p =>
        '<div class="plot-card"><a href="'+esc(p)+'" target="_blank">'+esc(p)+'</a><img src="'+esc(p)+'" loading="lazy" onerror="this.style.display=\\'none\\';"></div>'
      ).join('') || '<p class="note">No matching plots.</p>';
    }}
    document.getElementById('plot-filter').addEventListener('input', renderPlots);
    renderPlots();

    document.getElementById('tab-raw').innerHTML = '<h2>Embedded Data Keys</h2>' + keysTable(D);
  }} catch (err) {{
    status('ERROR: '+(err && (err.stack || err.message) || String(err)), 'err');
  }}
  </script>
</body>
</html>
"""

    backup = results_dir / "interactive_report_before_robust_viewer.html"
    if out_path.exists() and not backup.exists():
        shutil.copy2(out_path, backup)
    out_path.write_text(html, encoding="utf-8")
    debug_path.write_text(html, encoding="utf-8")
    print(f"Wrote {out_path} size={out_path.stat().st_size}")
    print(f"Wrote {debug_path} size={debug_path.stat().st_size}")


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results_replogle_weissman_k562_essential_mini20"
    build(target)
