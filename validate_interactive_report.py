#!/usr/bin/env python
"""Validate the generated interactive_report.html without a browser."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


def _collect_paths(obj, out):
    if isinstance(obj, dict):
        for value in obj.values():
            _collect_paths(value, out)
    elif isinstance(obj, list):
        for value in obj:
            _collect_paths(value, out)
    elif isinstance(obj, str):
        if obj.startswith("plots/") or obj.startswith("data:image/"):
            out.append(obj)


def main() -> int:
    report = Path(sys.argv[1]).resolve()
    root = report.parent
    text = report.read_text(encoding="utf-8", errors="replace")

    print(f"report={report}")
    print(f"size={report.stat().st_size}")

    if report.name.startswith("._"):
        print("ERROR: this is an AppleDouble metadata sidecar, not the real HTML report")
        return 2

    print(f"script_tags={text.count('<script')}, closing_script_tags={text.count('</script>')}")
    print(f"cdn_plotly={'cdn.plot.ly' in text}, local_plotly_ref={'plotly.min.js' in text}")

    plotly = root / "plotly.min.js"
    print(f"plotly_exists={plotly.exists()} size={plotly.stat().st_size if plotly.exists() else 0}")

    m = re.search(r'let\s+D\s*=\s*JSON\.parse\("(?P<payload>.*)"\);\n', text)
    if not m:
        print("ERROR: could not find standalone JSON.parse data assignment")
        return 2

    payload = m.group("payload")
    try:
        decoded = json.loads(f'"{payload}"')
        data = json.loads(decoded)
    except Exception as exc:
        print(f"ERROR: embedded JSON.parse payload is invalid: {exc}")
        return 2

    print(f"summary={data.get('summary')}")
    print(f"keys={len(data)} perturbations={len(data.get('perturbations', []))}")

    paths = []
    _collect_paths(data, paths)
    plot_paths = sorted({p for p in paths if p.startswith("plots/")})
    missing = [p for p in plot_paths if not (root / p).exists()]
    print(f"referenced_plot_paths={len(plot_paths)} missing_plot_paths={len(missing)}")
    for p in missing[:20]:
        print(f"MISSING: {p}")

    apple = sorted(root.glob("._interactive_report*.html"))
    print(f"appledouble_sidecars={len(apple)}")
    for p in apple[:10]:
        print(f"SIDECAR: {p.name} size={p.stat().st_size}")

    if missing or not plotly.exists():
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
