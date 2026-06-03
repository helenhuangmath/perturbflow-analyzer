#!/usr/bin/env python3
from __future__ import annotations

import html
import json
import shutil
import sys
from pathlib import Path


def esc(value) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def table(rows, columns=None, limit=300) -> str:
    if not rows:
        return '<p class="note">No rows available.</p>'
    rows = list(rows)[:limit]
    if columns is None:
        keys = []
        seen = set()
        for row in rows[:50]:
            if isinstance(row, dict):
                for key in row:
                    if key not in seen:
                        seen.add(key)
                        keys.append(key)
        columns = keys
    head = "".join(f"<th>{esc(c)}</th>" for c in columns)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{esc(fmt((row or {}).get(c)))}</td>" for c in columns) + "</tr>")
    note = f'<p class="note">Showing first {limit} rows.</p>' if len(rows) == limit else ""
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body_rows)}</tbody></table></div>{note}'


def collect_plot_paths(obj, out: set[str]) -> None:
    if isinstance(obj, dict):
        for value in obj.values():
            collect_plot_paths(value, out)
    elif isinstance(obj, list):
        for value in obj:
            collect_plot_paths(value, out)
    elif isinstance(obj, str) and obj.startswith("plots/"):
        out.add(obj)


def keys_table(data: dict) -> str:
    rows = []
    for key in sorted(data):
        value = data[key]
        if isinstance(value, list):
            typ, size = "array", len(value)
        elif isinstance(value, dict):
            typ, size = "object", len(value)
        else:
            typ, size = type(value).__name__, ""
        rows.append({"key": key, "type": typ, "size": size})
    return table(rows, ["key", "type", "size"], limit=1000)


def safe_json(data) -> str:
    return json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def build(results_dir: Path) -> None:
    results_dir = results_dir.resolve()
    data = json.loads((results_dir / "interactive_data.json").read_text(encoding="utf-8"))
    summary = data.get("summary", {})

    plots: set[str] = set()
    collect_plot_paths(data, plots)
    plots = sorted(p for p in plots if (results_dir / p).exists())

    plot_cards = []
    for path in plots[:120]:
        plot_cards.append(
            '<div class="plot-card">'
            f'<a href="{esc(path)}" target="_blank">{esc(path)}</a>'
            f'<img src="{esc(path)}" loading="lazy">'
            "</div>"
        )

    config_rows = []
    config = data.get("config") or {}
    for key in sorted(config):
        value = config[key]
        if isinstance(value, (dict, list)):
            value = json.dumps(value, ensure_ascii=False)[:500]
        config_rows.append({"parameter": key, "value": value})

    perts = data.get("cells_per_pert", [])
    sidebar_buttons = "".join(
        f'<a class="side-pert" href="#tab-perturbations" data-pert="{esc(row.get("perturbation"))}">{esc(row.get("perturbation"))}<span>{esc(row.get("n_cells"))}</span></a>'
        for row in perts
    )
    chart_payload = safe_json(
        {
            "cells_per_pert": perts,
            "deg_summary": data.get("deg_summary", []),
            "umap": data.get("umap", []),
        }
    )

    html_text = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>PerturbScope Static V2 Interactive Report</title>
  <style>
    :root{{--bg:#f6f8fb;--panel:#fff;--text:#172331;--muted:#5d6b7a;--border:#cbd7e3;--accent:#1f5f8b}}
    *{{box-sizing:border-box}} body{{margin:0;font-family:Arial,Helvetica,sans-serif;background:var(--bg);color:var(--text)}}
    header{{padding:18px 24px;background:#12344f;color:white}} header h1{{margin:0 0 6px;font-size:24px}} header .sub{{color:#dce8f4;font-size:14px}}
    #status{{margin:14px 24px 0;padding:10px 12px;border:2px solid #1f7a4d;background:#ecfff4;color:#123;font:13px/1.45 monospace}}
    nav{{display:flex;gap:6px;padding:16px 24px 0;border-bottom:1px solid var(--border);background:white;flex-wrap:wrap}}
    nav a{{border:0;background:transparent;border-bottom:3px solid transparent;padding:10px 14px;cursor:pointer;color:var(--muted);font-weight:700;font-size:14px;text-decoration:none}}
    nav a:hover{{color:var(--accent)}} nav a.active{{color:var(--accent);border-bottom-color:var(--accent)}}
    .layout{{display:grid;grid-template-columns:260px minmax(0,1fr);gap:0;align-items:start}}
    aside{{position:sticky;top:0;max-height:calc(100vh - 1px);overflow:auto;background:#fff;border-right:1px solid var(--border);padding:16px}}
    aside h3{{font-size:13px;margin:0 0 10px;color:#33485d;text-transform:uppercase}}
    .side-pert{{width:100%;display:flex;justify-content:space-between;gap:8px;border:1px solid var(--border);background:#f8fbfe;color:#172331;border-radius:6px;padding:7px 9px;margin-bottom:5px;cursor:pointer;text-align:left;font-size:13px;text-decoration:none}}
    .side-pert:hover,.side-pert.active{{border-color:var(--accent);background:#eaf4fc}}
    .side-pert span{{color:var(--muted);font-variant-numeric:tabular-nums}}
    main{{padding:20px 24px 40px;min-width:0}}
    .tab{{display:block;margin-bottom:28px;scroll-margin-top:12px}}
    .tab h2:first-child{{border-top:1px solid var(--border);padding-top:18px}}
    .grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin-bottom:18px}}
    .card{{background:var(--panel);border:1px solid var(--border);border-radius:8px;padding:14px}} .label{{color:var(--muted);font-size:12px;text-transform:uppercase}} .value{{font-size:26px;font-weight:800;margin-top:5px}}
    .chart-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:14px;margin:14px 0 20px}}
    .chart{{min-height:330px;background:white;border:1px solid var(--border);border-radius:8px;padding:8px}}
    .table-wrap{{max-height:540px;overflow:auto;border:1px solid var(--border);background:white}} table{{width:100%;border-collapse:collapse;background:white}} th,td{{padding:8px 10px;border-bottom:1px solid var(--border);text-align:left;font-size:13px;vertical-align:top}} th{{background:#eaf1f7;color:#33485d;position:sticky;top:0}}
    .plot-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}} .plot-card{{background:white;border:1px solid var(--border);border-radius:8px;padding:10px}} .plot-card img{{width:100%;height:auto;display:block;border:1px solid #e5edf4;background:white}} .plot-card a{{display:block;margin-bottom:8px;color:var(--accent);font-size:13px;overflow-wrap:anywhere}}
    .note{{color:var(--muted);font-size:13px;margin:8px 0 14px}}
    @media(max-width:900px){{.layout{{grid-template-columns:1fr}} aside{{position:relative;max-height:280px;border-right:0;border-bottom:1px solid var(--border)}}}}
  </style>
</head>
<body>
  <header><h1>PerturbScope Interactive Report</h1><div class="sub">Static-rendered robust report; no data-loading JavaScript required</div></header>
  <div id="status">Data loaded statically: {esc(summary.get("n_cells"))} cells, {esc(summary.get("n_genes"))} genes, {esc(summary.get("n_perturbations"))} perturbations</div>
  <nav>
    <a class="active" href="#tab-summary">Summary</a>
    <a href="#tab-perturbations">Perturbations</a>
    <a href="#tab-deg">DEG</a>
    <a href="#tab-qc">QC</a>
    <a href="#tab-plots">Plots</a>
    <a href="#tab-params">Parameters</a>
    <a href="#tab-keys">Data Keys</a>
  </nav>
  <div class="layout">
  <aside>
    <h3>Perturbations</h3>
    {sidebar_buttons}
  </aside>
  <main>
    <section id="tab-summary" class="tab">
      <div class="grid">
        <div class="card"><div class="label">Cells</div><div class="value">{esc(summary.get("n_cells"))}</div></div>
        <div class="card"><div class="label">Genes</div><div class="value">{esc(summary.get("n_genes"))}</div></div>
        <div class="card"><div class="label">Perturbations</div><div class="value">{esc(summary.get("n_perturbations"))}</div></div>
        <div class="card"><div class="label">Plot Assets</div><div class="value">{len(plots)}</div></div>
      </div>
      <div class="chart-grid">
        <div id="cells-chart" class="chart"></div>
        <div id="deg-chart" class="chart"></div>
      </div>
      <h2>Cells Per Perturbation</h2>{table(data.get("cells_per_pert", []), ["perturbation", "n_cells", "is_control"], 1000)}
    </section>
    <section id="tab-perturbations" class="tab"><h2>Perturbations</h2><div id="selected-pert-card" class="card"><div class="label">Selected perturbation</div><div class="value">All</div></div>{table(data.get("cells_per_pert", []), ["perturbation", "n_cells", "is_control"], 1000)}</section>
    <section id="tab-deg" class="tab"><h2>DEG Summary</h2>{table(data.get("deg_summary", []), None, 1000)}<h2>Available DEG Tables</h2>{keys_table(data.get("deg") or {})}</section>
    <section id="tab-qc" class="tab"><h2>QC Cell Sample</h2><p class="note">Showing embedded QC rows.</p>{table(data.get("qc_cells", []), None, 300)}</section>
    <section id="tab-plots" class="tab"><h2>Plot Browser</h2><p class="note">Showing up to 120 referenced plot assets.</p><div class="plot-grid">{"".join(plot_cards)}</div></section>
    <section id="tab-params" class="tab"><h2>Configuration</h2>{table(config_rows, ["parameter", "value"], 1000)}</section>
    <section id="tab-keys" class="tab"><h2>Embedded Data Keys</h2>{keys_table(data)}</section>
  </main>
  </div>
  <script type="application/json" id="chart-data">{chart_payload}</script>
  <script src="plotly.min.js"></script>
  <script>
    document.querySelectorAll('nav a').forEach(a => a.addEventListener('click', () => {{
      document.querySelectorAll('nav a').forEach(x => x.classList.toggle('active', x === a));
    }}));
    document.querySelectorAll('.side-pert').forEach(btn => btn.addEventListener('click', () => {{
      document.querySelectorAll('.side-pert').forEach(b => b.classList.toggle('active', b === btn));
      const card = document.getElementById('selected-pert-card');
      if (card) card.innerHTML = '<div class="label">Selected perturbation</div><div class="value">' + btn.dataset.pert + '</div>';
    }}));
    try {{
      const C = JSON.parse(document.getElementById('chart-data').textContent || '{{}}');
      if (window.Plotly) {{
        const cells = C.cells_per_pert || [];
        Plotly.newPlot('cells-chart', [{{type:'bar', x:cells.map(r=>r.perturbation), y:cells.map(r=>r.n_cells), marker:{{color:'#1f5f8b'}}}}],
          {{title:'Cells per perturbation', margin:{{l:50,r:20,t:45,b:110}}, xaxis:{{tickangle:-45}}, yaxis:{{title:'Cells'}}}},
          {{responsive:true, displaylogo:false}});
        const deg = C.deg_summary || [];
        const numericKeys = deg.length ? Object.keys(deg[0]).filter(k => typeof deg[0][k] === 'number') : [];
        const yKey = numericKeys.find(k => /sig|deg|n_|count|up|down/i.test(k)) || numericKeys[0];
        const xKey = deg.length && ('perturbation' in deg[0] ? 'perturbation' : Object.keys(deg[0])[0]);
        if (yKey && xKey) {{
          Plotly.newPlot('deg-chart', [{{type:'bar', x:deg.map(r=>r[xKey]), y:deg.map(r=>r[yKey]), marker:{{color:'#1f7a4d'}}}}],
            {{title:'DEG summary: '+yKey, margin:{{l:50,r:20,t:45,b:110}}, xaxis:{{tickangle:-45}}}},
            {{responsive:true, displaylogo:false}});
        }} else {{
          document.getElementById('deg-chart').innerHTML = '<p class="note">No numeric DEG summary column available for charting.</p>';
        }}
      }}
    }} catch (err) {{
      document.getElementById('cells-chart').innerHTML = '<p class="note">Interactive charts unavailable: '+String(err.message || err)+'</p>';
    }}
  </script>
</body>
</html>
"""

    out = results_dir / "interactive_report.html"
    debug = results_dir / "interactive_report_STATIC_V2_OPEN_THIS.html"
    backup = results_dir / "interactive_report_before_static_viewer.html"
    if out.exists() and not backup.exists():
        shutil.copy2(out, backup)
    out.write_text(html_text, encoding="utf-8")
    debug.write_text(html_text, encoding="utf-8")
    print(f"Wrote {out} size={out.stat().st_size}")
    print(f"Wrote {debug} size={debug.stat().st_size}")
    print(f"plots={len(plots)}")


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results_replogle_weissman_k562_essential_mini20"
    build(target)
