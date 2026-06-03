#!/usr/bin/env python3
"""Restore the original interactive report layout with local Plotly.

The simplified/static viewers were useful for debugging, but this restores the
old report UI that has the sidebar, original tabs, and full interactive Plotly
rendering.  The only functional patch is replacing the CDN Plotly script with
the local plotly.min.js next to the report, plus a small visible error reporter
so browser-side failures do not silently leave the page at Loading.
"""

from __future__ import annotations

import shutil
from pathlib import Path


RESULTS = Path("/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20")
SRC = RESULTS / "interactive_report old.html"
OUT = RESULTS / "interactive_report.html"
DEBUG = RESULTS / "interactive_report_ORIGINAL_LAYOUT_OPEN_THIS.html"


ERROR_HELPER = """\
<script>
window.addEventListener('error', function(ev) {
  try {
    var box = document.getElementById('browser-startup-error');
    if (!box) {
      box = document.createElement('div');
      box.id = 'browser-startup-error';
      box.style.cssText = 'position:relative;z-index:999999;margin:10px 18px;padding:12px;border:3px solid #b91c1c;background:#fff1f2;color:#111;font:13px/1.45 monospace;white-space:pre-wrap;';
      document.body.insertBefore(box, document.body.firstChild);
    }
    box.textContent = 'Browser JavaScript error:\\n' + (ev.message || '') + '\\n' + (ev.filename || '') + ':' + (ev.lineno || '') + ':' + (ev.colno || '');
  } catch (e) {}
});
</script>"""


def main() -> int:
    if not SRC.exists():
        raise SystemExit(f"Missing source old-layout report: {SRC}")
    if not (RESULTS / "plotly.min.js").exists():
        raise SystemExit(f"Missing local plotly.min.js next to report: {RESULTS / 'plotly.min.js'}")

    text = SRC.read_text(encoding="utf-8", errors="replace")
    text = text.replace(
        '<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>',
        '<script src="plotly.min.js"></script>\n' + ERROR_HELPER,
        1,
    )
    text = text.replace(
        "<title>PerturbFlow-Analyzer — Interactive Report</title>",
        "<title>PerturbScope Original Layout Interactive Report</title>",
        1,
    )

    backup = RESULTS / "interactive_report_before_original_layout_restore.html"
    if OUT.exists() and not backup.exists():
        shutil.copy2(OUT, backup)

    OUT.write_text(text, encoding="utf-8")
    DEBUG.write_text(text, encoding="utf-8")
    print(f"Wrote {OUT} size={OUT.stat().st_size}")
    print(f"Wrote {DEBUG} size={DEBUG.stat().st_size}")
    print("Patched Plotly: CDN -> local plotly.min.js")
    print("Original sidebar/tabs/interactions preserved from interactive_report old.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
