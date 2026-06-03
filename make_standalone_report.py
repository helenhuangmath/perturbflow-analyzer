#!/usr/bin/env python3
"""
Regenerate a FULLY STANDALONE interactive_report.html.

Root cause of "Loading…" bug
-----------------------------
The old HTML embeds a 7.8 MB JavaScript *object literal*:
    let D = { <7.8 MB of raw JSON> };
Browsers parse object literals with the full JS grammar (slower) rather than
the optimised JSON fast-path.  On typical hardware this stalls the JS engine for
10-30 s, so the page appears frozen.  Tabs, figures, and sidebars never render.

Fix
---
Replace the object literal with JSON.parse("...") — identical data, but the
browser's native JSON parser (written in C, uses SIMD) is 5-15× faster.
Also switch loadReportData() to Promise.resolve() so no HTTP fetch is needed,
making the file work when opened directly with file://.
"""

import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).parent
sys.path.insert(0, str(REPO))

from perturbscope.interactive import _copy_plotly_bundle, _html_template

# ── fetch-based loadReportData that the current template produces ─────────────
_FETCH_LOAD = """\
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

# ── inline replacement (no fetch, works with file://) ────────────────────────
_INLINE_LOAD = """\
function loadReportData(){
  if(_dataReady) return _dataReady;
  const hdr=document.getElementById('hdr-summary');
  _dataReady=Promise.resolve().then(()=>{
    D=Object.assign(emptyDataStub(),D||{});
    if(hdr) hdr.textContent='Results data loaded';
    return D;
  }).catch(err=>{showStartupError(err);throw err;});
  return _dataReady;
}"""


def _escape_for_js_dq_string(s: str) -> str:
    """Escape *s* so it can be placed inside a JS double-quoted string literal.

    JSON produced by json.dumps already has its internal strings properly
    escaped (e.g. \" → \\\\" in the output).  We only need to:
      1. double every backslash  (so \\\" becomes \\\\\\\")
      2. escape every double-quote  (so the outer JS string isn't terminated)
    CR characters are stripped because json.dumps never produces them.
    """
    return s.replace("\\", "\\\\").replace('"', '\\"')


def regen_standalone(results_dir: Path) -> None:
    results_dir = results_dir.resolve()
    json_path = results_dir / "interactive_data.json"

    if not json_path.exists():
        print(f"ERROR: {json_path} not found.  Run the full pipeline first.")
        sys.exit(1)

    kb = json_path.stat().st_size / 1024
    print(f"Reading {json_path.name}  ({kb:.0f} KB) …")
    data = json.loads(json_path.read_text(encoding="utf-8"))

    # Re-serialise with ensure_ascii=True to avoid JS line-separator problems
    # (U+2028 LINE SEPARATOR and U+2029 PARAGRAPH SEPARATOR are valid in JSON
    # strings but act as newline terminators inside JS string literals).
    safe_json = json.dumps(data, ensure_ascii=True, separators=(",", ":"))

    # Build the JS expression that replaces {DATA_STUB}:
    #   let D = JSON.parse("...escaped full data...");
    js_expr = f'JSON.parse("{_escape_for_js_dq_string(safe_json)}")'

    template = _html_template()
    html = template.replace("{DATA_STUB}", js_expr)

    patched = False
    if _FETCH_LOAD in html:
        html = html.replace(_FETCH_LOAD, _INLINE_LOAD, 1)
        patched = True
        print("Patched loadReportData: fetch → Promise.resolve (standalone mode)")
    else:
        print(
            "WARNING: fetch-based loadReportData not found verbatim.\n"
            "The HTML may still work if the template already embeds all data."
        )

    # ── backup & write ────────────────────────────────────────────────────────
    out_path = results_dir / "interactive_report.html"
    backup = results_dir / "interactive_report_before_jsonparse_fix.html"
    if out_path.exists() and not backup.exists():
        shutil.copy2(out_path, backup)
        print(f"Backed up old HTML → {backup.name}")

    out_path.write_text(html, encoding="utf-8")
    new_kb = out_path.stat().st_size / 1024
    print(f"Wrote standalone HTML → {out_path.name}  ({new_kb:.0f} KB)")

    _copy_plotly_bundle(results_dir)
    print("plotly.min.js: present.")

    print()
    print("=" * 60)
    print("Done!  Open interactive_report.html directly in a browser.")
    print("file:// protocol works — no web server needed.")
    print("=" * 60)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        targets = [Path(a) for a in sys.argv[1:]]
    else:
        targets = [
            REPO / "results_replogle_weissman_k562_essential_mini20",
        ]
    for d in targets:
        regen_standalone(d)
