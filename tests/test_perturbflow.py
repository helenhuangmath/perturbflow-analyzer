from __future__ import annotations

import json
from pathlib import Path

from perturbflow.ai import write_agent_handoff
from perturbflow.cli import build_parser


def test_cli_exposes_public_commands() -> None:
    parser = build_parser()
    help_text = parser.format_help()
    assert "prepare" in help_text
    assert "run" in help_text
    assert "interpret" in help_text


def test_agent_handoff_writes_expected_files(tmp_path: Path) -> None:
    results = tmp_path / "results"
    csv_dir = results / "csv"
    csv_dir.mkdir(parents=True)
    (results / "summary.json").write_text(
        json.dumps({"n_cells": 100, "n_genes": 2000, "n_perturbations": 3}),
        encoding="utf-8",
    )
    (results / "checkpoint.json").write_text(
        json.dumps({"completed_steps": ["qc", "deg", "report", "bundle"]}),
        encoding="utf-8",
    )
    (csv_dir / "deg_summary.csv").write_text(
        "perturbation,n_de_total,n_de_up,n_de_down,top_up_gene\nA,12,7,5,GENE1\n",
        encoding="utf-8",
    )

    files = write_agent_handoff(results, project_name="Test")

    assert Path(files["interpretation_context"]).exists()
    assert Path(files["agent_prompt"]).exists()
    assert Path(files["machine_context"]).exists()
    assert Path(files["agent_manifest"]).exists()
    text = Path(files["interpretation_context"]).read_text(encoding="utf-8")
    assert "Test interpretation context" in text
    assert "GENE1" in text
