# AI And Agents

PerturbFlow is designed to connect analysis outputs to LLMs and agents without making hidden network calls.

## Export Interpretation Context

```bash
perturbflow interpret --results results/my_run --project-name "My screen"
```

This writes:

```text
agent_handoff/
├── agent_manifest.json
├── agent_prompt.md
├── interpretation_context.md
└── machine_context.json
```

## Suggested Agent Roles

- QC reviewer.
- Perturbation prioritizer.
- Pathway interpreter.
- Network rewiring interpreter.
- Report writer.

## Privacy

The handoff exporter summarizes derived tables and report artifacts. It does not include raw count matrices. Review the files before sending them to an external model provider.
