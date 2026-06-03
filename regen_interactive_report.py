#!/usr/bin/env python3
"""
Regenerate interactive_report.html from the existing interactive_data.json
WITHOUT re-running the full pipeline or loading the h5ad file.

The old HTML embedded 7.8 MB of JSON as a JavaScript literal (let D = {...}),
causing the browser's JS engine to stall or fail to parse it.
The current generator uses fetch('interactive_data.json') + a tiny stub,
producing a ~50 KB HTML that loads the data lazily.

Usage:
    python regen_interactive_report.py <results_dir>

If no results_dir is given, defaults to the mini20 results folder.
"""

import json
import shutil
import sys
from pathlib import Path

# Allow running from anywhere in the repo
REPO = Path(__file__).parent
sys.path.insert(0, str(REPO))

from perturbscope.interactive import _html_template, _copy_plotly_bundle


def regen(results_dir: Path) -> None:
    results_dir = results_dir.resolve()
    json_path = results_dir / "interactive_data.json"

    if not json_path.exists():
        print(f"ERROR: {json_path} not found. Run the report step first.")
        sys.exit(1)

    print(f"Reading {json_path.name}  ({json_path.stat().st_size / 1024:.0f} KB) ...")
    data = json.loads(json_path.read_text(encoding="utf-8"))

    # Build the minimal stub that gets embedded in the HTML
    stub = {
        "summary": data.get("summary", {}),
        "config":  data.get("config",  {}),
    }
    stub_str = json.dumps(stub, ensure_ascii=False, separators=(",", ":"))

    template = _html_template()
    html = template.replace("{DATA_STUB}", stub_str)

    out_path = results_dir / "interactive_report.html"

    # Back up the old file once
    backup = results_dir / "interactive_report_old_embedded.html"
    if out_path.exists() and not backup.exists():
        shutil.copy2(out_path, backup)
        print(f"Backed up old HTML → {backup.name}")

    out_path.write_text(html, encoding="utf-8")
    new_kb = out_path.stat().st_size / 1024
    print(f"Wrote new HTML → {out_path.name}  ({new_kb:.0f} KB)")

    # Make sure plotly bundle is present
    _copy_plotly_bundle(results_dir)

    print()
    print("=" * 60)
    print("The new HTML uses fetch('interactive_data.json').")
    print("It must be served over HTTP — open via file:// will NOT work")
    print("in Chrome (CORS). Run ONE of the following:")
    print()
    print(f"  cd {results_dir}")
    print("  python -m http.server 8899")
    print()
    print("Then open in a browser:")
    print("  http://localhost:8899/interactive_report.html")
    print()
    print("If on a remote cluster with VS Code, forward port 8899:")
    print("  VS Code → Ports panel → Forward Port 8899")
    print("=" * 60)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        target = Path(sys.argv[1])
    else:
        target = Path(
            "/vast/projects/wherry/foundation-models-immuno/hhua/"
            "sc_perturbation/PerturbScope_v1/"
            "results_replogle_weissman_k562_essential_mini20"
        )
    regen(target)
