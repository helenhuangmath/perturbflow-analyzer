#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> int:
    results_dir = Path(sys.argv[1]).resolve()
    data = json.loads((results_dir / "interactive_data.json").read_text(encoding="utf-8"))
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    out = results_dir / "OPEN_THIS_DATA_PARSE_TEST.html"
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Interactive Data Parse Test</title>
  <style>
    body{{font-family:Arial,sans-serif;margin:40px;background:#f8fafc;color:#102033}}
    .box{{border:4px solid #1f4e79;background:#eef7ff;padding:24px;border-radius:8px;max-width:1000px}}
    pre{{background:#fff;border:1px solid #cbd5e1;padding:12px;white-space:pre-wrap}}
  </style>
</head>
<body>
  <div class="box">
    <h1>Interactive Data Parse Test</h1>
    <p>This page uses the same embedded JSON data as the interactive report, but no report UI.</p>
    <pre id="out">JavaScript has not run yet.</pre>
  </div>
  <script type="application/json" id="interactive-data-json">
{payload}
  </script>
  <script>
  const out = document.getElementById('out');
  function log(msg) {{ out.textContent += "\\n" + msg; }}
  out.textContent = "JavaScript started";
  try {{
    const el = document.getElementById('interactive-data-json');
    log("data block chars: " + (el.textContent || "").length);
    const data = JSON.parse(el.textContent || "{{}}");
    log("JSON.parse OK");
    log("n_cells: " + data.summary.n_cells);
    log("n_genes: " + data.summary.n_genes);
    log("perturbations: " + data.perturbations.length);
    log("first perturbation: " + data.perturbations[0]);
  }} catch (err) {{
    log("ERROR: " + (err && (err.stack || err.message) || String(err)));
  }}
  </script>
</body>
</html>
"""
    out.write_text(html, encoding="utf-8")
    print(f"Wrote {out}")
    print(f"Size: {out.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
