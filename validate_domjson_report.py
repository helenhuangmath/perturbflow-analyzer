#!/usr/bin/env python3
"""Validate the DOM-JSON standalone report format without a browser."""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path


def _collect_plot_paths(obj, out):
    if isinstance(obj, dict):
        for value in obj.values():
            _collect_plot_paths(value, out)
    elif isinstance(obj, list):
        for value in obj:
            _collect_plot_paths(value, out)
    elif isinstance(obj, str) and obj.startswith("plots/"):
        out.append(obj)


def main() -> int:
    report = Path(sys.argv[1]).resolve()
    root = report.parent
    text = report.read_text(encoding="utf-8", errors="replace")

    print(f"report={report}")
    print(f"size={report.stat().st_size}")
    print(f"fetch_calls={text.count('fetch(')}")
    print(f"json_parse_string_assignment={'let D = JSON.parse(' in text}")
    has_dom_json = 'id="interactive-data-json"' in text
    print(f"dom_json_block={has_dom_json}")
    print(f"local_plotly_ref={'plotly.min.js' in text}")

    match = re.search(
        r'<script type="application/json" id="interactive-data-json">\n(?P<data>.*?)\n</script>',
        text,
        flags=re.S,
    )
    if not match:
        print("ERROR: embedded DOM JSON block not found")
        return 2

    payload = html.unescape(match.group("data")).replace("<\\/", "</")
    data = json.loads(payload)
    print(f"summary={data.get('summary')}")
    print(f"keys={len(data)} perturbations={len(data.get('perturbations', []))}")

    plotly = root / "plotly.min.js"
    print(f"plotly_exists={plotly.exists()} size={plotly.stat().st_size if plotly.exists() else 0}")

    plot_paths = []
    _collect_plot_paths(data, plot_paths)
    missing = sorted({p for p in plot_paths if not (root / p).exists()})
    print(f"referenced_plot_paths={len(set(plot_paths))} missing_plot_paths={len(missing)}")
    for path in missing[:20]:
        print(f"MISSING: {path}")

    return 1 if missing or not plotly.exists() else 0


if __name__ == "__main__":
    raise SystemExit(main())
