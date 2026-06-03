#!/usr/bin/env python3
from pathlib import Path

results_dir = Path("/vast/projects/wherry/foundation-models-immuno/hhua/sc_perturbation/PerturbScope_v1/results_replogle_weissman_k562_essential_mini20")
out = results_dir / "OPEN_THIS_PLOTLY_TEST.html"
out.write_text(
    """<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Plotly Test</title>
  <style>
    body{font-family:Arial,sans-serif;margin:32px;background:#f8fafc;color:#102033}
    #status{padding:10px 12px;border:2px solid #1f4e79;background:#eef7ff;margin-bottom:16px;font:13px/1.45 monospace}
    #chart{height:420px;border:1px solid #cbd5e1;background:white}
  </style>
</head>
<body>
  <h1>Plotly Test</h1>
  <div id="status">Starting...</div>
  <div id="chart"></div>
  <script src="plotly.min.js"></script>
  <script>
    const status = document.getElementById('status');
    try {
      status.textContent = 'Plotly type: ' + typeof Plotly;
      if (typeof Plotly === 'undefined') throw new Error('Plotly is undefined');
      Plotly.newPlot('chart', [{type:'bar', x:['A','B','C'], y:[1,3,2]}],
        {title:'If you see bars, local Plotly works'}, {responsive:true, displaylogo:false});
      status.textContent += '\\nPlotly.newPlot OK';
    } catch (err) {
      status.textContent += '\\nERROR: ' + (err && (err.stack || err.message) || String(err));
    }
  </script>
</body>
</html>
""",
    encoding="utf-8",
)
print(f"Wrote {out} size={out.stat().st_size}")
