"""Structured interpretation exports for PerturbFlow result folders."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _read_csv_rows(path: Path, limit: int = 20) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [row for _, row in zip(range(limit), reader)]


def _numeric(row: dict[str, Any], key: str) -> float:
    try:
        return float(row.get(key, 0) or 0)
    except (TypeError, ValueError):
        return 0.0


def _top_rows(rows: list[dict[str, str]], keys: list[str], limit: int) -> list[dict[str, str]]:
    for key in keys:
        if rows and key in rows[0]:
            return sorted(rows, key=lambda row: abs(_numeric(row, key)), reverse=True)[:limit]
    return rows[:limit]


def _markdown_table(rows: list[dict[str, str]], preferred_cols: list[str]) -> str:
    if not rows:
        return "No table was found.\n"
    cols = [col for col in preferred_cols if col in rows[0]]
    if not cols:
        cols = list(rows[0].keys())[:6]
    header = "| " + " | ".join(cols) + " |"
    divider = "| " + " | ".join(["---"] * len(cols)) + " |"
    body = []
    for row in rows:
        values = [str(row.get(col, "")).replace("\n", " ")[:120] for col in cols]
        body.append("| " + " | ".join(values) + " |")
    return "\n".join([header, divider, *body]) + "\n"


def collect_result_context(results_dir: str | Path, max_rows: int = 10) -> dict[str, Any]:
    """Collect compact, non-expression-matrix context from a completed run."""
    root = Path(results_dir)
    csv_dir = root / "csv"
    summary = _read_json(root / "summary.json")
    checkpoint = _read_json(root / "checkpoint.json")
    bundle_manifest = _read_json(root / "bundle" / "manifest.json")

    deg_rows = _top_rows(
        _read_csv_rows(csv_dir / "deg_summary.csv", limit=200),
        ["n_de_total", "n_de_up", "n_de_down"],
        max_rows,
    )
    cscore_rows = _top_rows(
        _read_csv_rows(csv_dir / "cscore_summary.csv", limit=200),
        ["c_total", "c_shift", "c_gain", "c_loss"],
        max_rows,
    )
    effect_rows = _top_rows(
        _read_csv_rows(csv_dir / "effect_decomposition.csv", limit=200),
        ["total_effect", "transcriptional_effect", "state_shift_effect", "effect_size"],
        max_rows,
    )
    trajectory_rows = _top_rows(
        _read_csv_rows(csv_dir / "trajectory_effects.csv", limit=200),
        ["trajectory_shift", "pseudotime_shift", "effect_size"],
        max_rows,
    )

    return {
        "results_dir": str(root),
        "summary": summary,
        "checkpoint": checkpoint,
        "bundle_manifest": bundle_manifest,
        "deg_top": deg_rows,
        "cscore_top": cscore_rows,
        "effect_top": effect_rows,
        "trajectory_top": trajectory_rows,
    }


def build_interpretation_markdown(context: dict[str, Any], project_name: str = "PerturbFlow run") -> str:
    """Render a concise interpretation context document from result summaries."""
    summary = context.get("summary") or {}
    checkpoint = context.get("checkpoint") or {}
    bundle_manifest = context.get("bundle_manifest") or {}

    lines = [
        f"# {project_name} agent interpretation context",
        "",
        "This file is designed for a human analyst to interpret a completed PerturbFlow run.",
        "It summarizes derived tables and report artifacts; it does not include raw count matrices.",
        "",
        "## Dataset overview",
        "",
        f"- Results directory: `{context.get('results_dir', '')}`",
        f"- Cells: {summary.get('n_cells', 'unknown')}",
        f"- Genes: {summary.get('n_genes', 'unknown')}",
        f"- Perturbations: {summary.get('n_perturbations', 'unknown')}",
        f"- Completed steps: {', '.join(checkpoint.get('completed_steps', [])) or 'unknown'}",
        f"- Bundle schema: {bundle_manifest.get('schema_version', 'unknown')}",
        "",
        "## Top differential-expression perturbations",
        "",
        _markdown_table(
            context.get("deg_top", []),
            ["perturbation", "n_de_total", "n_de_up", "n_de_down", "top_up_gene", "top_down_gene"],
        ),
        "## Strongest connectivity rewiring signals",
        "",
        _markdown_table(
            context.get("cscore_top", []),
            ["perturbation", "c_total", "c_gain", "c_loss", "c_shift", "n_edges_control", "n_edges_perturbed"],
        ),
        "## Largest effect-decomposition signals",
        "",
        _markdown_table(
            context.get("effect_top", []),
            ["perturbation", "total_effect", "transcriptional_effect", "state_shift_effect", "dominant_effect"],
        ),
        "## Trajectory effects",
        "",
        _markdown_table(
            context.get("trajectory_top", []),
            ["perturbation", "trajectory_shift", "pseudotime_shift", "direction", "effect_size"],
        ),
        "## Suggested analyst questions",
        "",
        "1. Which perturbations have concordant DEG, state-shift, trajectory, and network-rewiring evidence?",
        "2. Which perturbations are strong only in one modality, and could that indicate a specific mechanism or artifact?",
        "3. Are top genes/pathways consistent with the expected biology of the perturbation targets?",
        "4. Which perturbations should be prioritized for validation, and what evidence supports each priority?",
        "5. What QC caveats, cell-state composition shifts, or low-cell-count cases should limit interpretation?",
        "",
    ]
    return "\n".join(lines)


def build_agent_prompt(project_name: str = "PerturbFlow run") -> str:
    """Return a reusable analysis prompt."""
    return f"""You are a careful single-cell perturbation analysis assistant.

Use the attached `{project_name} agent interpretation context` and the linked PerturbFlow result artifacts to write a scientist-facing interpretation.

Requirements:
- Separate strong findings from hypotheses.
- Cite the exact table or figure artifact that supports each claim.
- Mention QC limitations, low-cell-count perturbations, and multiple-testing caveats.
- Prioritize perturbations for follow-up with a short rationale.
- Do not infer causal biology from UMAP position alone.
- Do not request or expose raw expression matrices unless the user explicitly provides them.
"""


def write_agent_handoff(
    results_dir: str | Path,
    output_dir: str | Path | None = None,
    project_name: str = "PerturbFlow run",
    max_rows: int = 10,
) -> dict[str, str]:
    """Write structured Markdown and JSON files for a completed run."""
    root = Path(results_dir)
    out = Path(output_dir) if output_dir else root / "agent_handoff"
    out.mkdir(parents=True, exist_ok=True)

    context = collect_result_context(root, max_rows=max_rows)
    context["created_at"] = datetime.now(timezone.utc).isoformat()
    context["project_name"] = project_name

    context_md = build_interpretation_markdown(context, project_name=project_name)
    prompt_md = build_agent_prompt(project_name=project_name)
    manifest = {
        "project_name": project_name,
        "created_at": context["created_at"],
        "results_dir": str(root),
        "files": {
            "interpretation_context": "interpretation_context.md",
            "agent_prompt": "agent_prompt.md",
            "machine_context": "machine_context.json",
        },
        "recommended_review_roles": [
            "qc_reviewer",
            "perturbation_prioritizer",
            "pathway_interpreter",
            "network_rewiring_interpreter",
            "report_writer",
        ],
        "privacy_note": "Review files before sharing them outside your analysis environment. Raw count matrices are not included by this exporter.",
    }

    (out / "interpretation_context.md").write_text(context_md, encoding="utf-8")
    (out / "agent_prompt.md").write_text(prompt_md, encoding="utf-8")
    (out / "machine_context.json").write_text(json.dumps(context, indent=2), encoding="utf-8")
    (out / "agent_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    return {key: str(out / value) for key, value in manifest["files"].items()} | {
        "agent_manifest": str(out / "agent_manifest.json")
    }
