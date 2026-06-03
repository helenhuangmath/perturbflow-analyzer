#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path

report = Path(sys.argv[1]).resolve()
text = report.read_text(encoding="utf-8", errors="replace")
print(f"report={report}")
print(f"size={report.stat().st_size}")
print(f"title_loading={'<title>Loading' in text or '<title>Loading...' in text}")
print(f"fetch_calls={text.count('fetch(')}")
print(f"jsonparse_assignment={'let D = JSON.parse(' in text}")
has_dom_json = 'id="interactive-data-json"' in text
print(f"dom_json_block={has_dom_json}")
print(f"diagnostics_panel={'startup-diagnostics' in text}")
print(f"plotly_ref={'plotly.min.js' in text}")
print(f"script_tags={text.count('<script')} closing_script_tags={text.count('</script>')}")
print(f"appledouble={report.name.startswith('._')}")

match = re.search(
    r'<script type="application/json" id="interactive-data-json">\n(?P<data>.*?)\n</script>',
    text,
    flags=re.S,
)
if match:
    data = match.group("data")
    print(f"dom_json_chars={len(data)}")
    print(f"dom_json_prefix={data[:80]}")
else:
    print("dom_json_chars=0")
