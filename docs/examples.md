# Example Notebooks

The `examples/` folder contains runnable notebook templates for common PerturbFlow workflows.

| Notebook | Purpose |
| --- | --- |
| `examples/01_prepare_and_run.ipynb` | Prepare a user `.h5ad`, run PerturbFlow, and locate reports. |
| `examples/02_step_rerun_and_config.ipynb` | Tune config values and rerun selected analysis steps. |
| `examples/03_interpret_with_agents.ipynb` | Export structured interpretation context from completed results. |
| `examples/04_explore_outputs.ipynb` | Inspect summary JSON, CSV tables, reports, plots, and bundle files. |

Each notebook starts with path variables near the top. Update those paths before running. Heavy analysis cells are guarded by flags such as `RUN_ANALYSIS = False` so opening a notebook does not start a pipeline run by accident.
