#!/usr/bin/env python3
from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path


def main() -> int:
    report = Path(sys.argv[1]).resolve()
    text = report.read_text(encoding="utf-8", errors="replace")
    print(f"report={report}")
    print(f"size={report.stat().st_size}")
    checks = {
        "original_sidebar": '<aside class="sidebar">' in text,
        "original_tabs": 'class="tab-btn"' in text and "showTab('pert'" in text,
        "dom_json": 'id="interactive-data-json"' in text,
        "no_fetch_call": "fetch('interactive_data.json'" not in text,
        "local_plotly_loader": "plotly.min.js" in text,
        "no_cdn_plotly": "cdn.plot.ly" not in text,
        "startup_error_reporter": "unhandledrejection" in text,
    }
    for key, value in checks.items():
        print(f"{key}={value}")

    match = re.search(
        r'<script type="application/json" id="interactive-data-json">\n(?P<data>.*?)\n</script>',
        text,
        flags=re.S,
    )
    if not match:
        print("ERROR: DOM JSON block missing")
        return 2
    data = json.loads(html.unescape(match.group("data")).replace("<\\/", "</"))
    print(f"summary={data.get('summary')}")
    print(f"perturbations={len(data.get('perturbations', []))}")
    print(f"deg_summary={len(data.get('deg_summary', []))}")
    print(f"umap={len(data.get('umap', []))}")
    if not all(checks.values()):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
