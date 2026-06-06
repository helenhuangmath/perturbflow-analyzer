# Contributing to PerturbFlow

PerturbFlow is open infrastructure for perturbation biology, and it improves
through community contributions. Thank you for considering one.

## Ways to contribute

- **Report bugs or rough edges** — open an issue with a minimal reproduction
  (ideally a small `.h5ad` or synthetic AnnData) and the command you ran.
- **Improve documentation** — clarifications, examples, and tutorials are
  high-value and low-risk.
- **Propose features** — open an issue first, describing the use case. Changes
  that advance one of the three aims in [`ROADMAP.md`](ROADMAP.md) are
  especially welcome.
- **Contribute datasets or benchmarks** — reference dataset loaders (Aim 1) and
  evaluation metrics (Aim 3) are actively sought.

## Development setup

```bash
git clone https://github.com/helenhuangmath/PerturbFlow.git
cd PerturbFlow
python -m pip install -e ".[dev,bundle]"
pytest
```

## Pull request guidelines

1. **Branch** from `main` and keep each PR focused on one logical change.
2. **Preserve interoperability.** PerturbFlow is AnnData-native and complements
   the scverse ecosystem; new analysis steps should read and write AnnData, not
   bespoke formats.
3. **Keep pipeline contracts stable.** Add new pipeline steps as additional
   terminal steps rather than changing existing module inputs/outputs (see
   [`DESIGN.md`](DESIGN.md)).
4. **Test what you add.** Include or extend tests under `tests/`; `pytest` must
   pass.
5. **Document it.** Update the relevant page under `docs/`, the `README.md`, and
   `ROADMAP.md` status if your change advances an aim.
6. **Privacy by default.** PerturbFlow never uploads data to an external service
   automatically; AI integrations must keep the analyst in control of when and
   where data leaves their machine.

## Code style

Match the surrounding code: type hints with `from __future__ import
annotations`, clear function-level docstrings, and small composable modules.

## Community standards

Participation is governed by the [Code of Conduct](CODE_OF_CONDUCT.md). Be
constructive and assume good faith.

## License

By contributing, you agree that your contributions are licensed under the
project's [MIT License](LICENSE).
