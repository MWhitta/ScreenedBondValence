# ScreenedBondValence

`ScreenedBondValence` fits bond-valence parameters `R0` and `B` for cation-anion pairs from crystal structures.

The repo is now package-first:

- one Python package under `screened_bond_valence/`
- one CLI entrypoint, `screened-bond-valence`
- one test suite under `tests/`

Companion project docs:

- [AGENTS.md](AGENTS.md)
- [ARCHITECTURE.md](ARCHITECTURE.md)
- [docs/adr/](docs/adr/)
- [docs/integrations/critical_mineral_project.md](docs/integrations/critical_mineral_project.md)

The fitted bond-valence model uses:

`Sij = exp((R0 - Rij) / B)`

where `Sij` is theoretical bond valence and `Rij` is bond length.

## Installation

```bash
git clone https://github.com/MWhitta/ScreenedBondValence.git
cd ScreenedBondValence
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

You need a Materials Project API key for workflows that call MP directly.

## Python API

```python
from screened_bond_valence import (
    BondValenceProcessor,
    ScreenedBondValenceService,
    build_material_from_cif,
    build_summary_payload,
)
```

### Fit a single CIF

If oxidation states are missing from the CIF, the package asks `pymatgen`
to infer them with `BVAnalyzer`, then falls back to
`Structure.add_oxidation_state_by_guess()` if needed.

```python
from screened_bond_valence import ScreenedBondValenceService, build_material_from_cif

material = build_material_from_cif(
    "1011191_aspod.cif",
    cation="Li",
    anion="O",
)

service = ScreenedBondValenceService(algorithms=("shgo", "diff"))
result = service.fit_material(material)
summary = result.aggregate()

print(summary.to_dict())
```

### Fit Materials Project structures

```python
from screened_bond_valence import BondValenceProcessor

processor = BondValenceProcessor(
    api_key="your_api_key",
    algos=["shgo", "brute", "diff", "dual_annealing", "direct"],
    cations=["Li"],
    anions=["O"],
)

processor.process_cation_system("Li", "O")
```

The batch processor resumes previous runs by skipping materials that already
have fitted outputs for every configured algorithm, while preserving the
existing `dict_sijs.json`, `dict_charges.json`, and `no_solu/` records.

### Export a downstream JSON payload

`build_summary_payload` accepts:

- `classifier`: decide which output bucket each fitted material belongs to
- `serializer`: control the JSON schema for each summary

```python
from screened_bond_valence import build_summary_payload

payload = build_summary_payload(
    results,
    classifier=lambda summary: "hydroxides" if "OH" in (summary.formula_pretty or "") else "oxides",
    serializer=lambda summary: summary.to_compact_dict(),
)
```

This matches the minimal shape currently consumed by `CriticalMineralProject`.

## CLI

The installed CLI keeps the batch and web workflows available without
separate top-level scripts.

### Batch fitting

```bash
screened-bond-valence batch \
  --api-key "$MP_API_KEY" \
  --cations Li Na \
  --anions O \
  --algorithms shgo diff \
  --output-dir res
```

### Web app

```bash
screened-bond-valence web --server-name 127.0.0.1 --server-port 7860
```

You can also launch the CLI with:

```bash
python -m screened_bond_valence web
```

## Outputs

```text
res/
└── LiO/
    ├── params/
    │   └── dict_matID_possible_species.json
    ├── R0Bs/
    ├── no_solu/
    ├── dict_sijs.json
    └── dict_charges.json
```

## Development

```bash
PYTHONPATH=. python -m unittest discover -s tests -v
```

Golden regression coverage for representative sample materials lives in
`tests/golden/` and `tests/test_golden_outputs.py`.

## References

- Brown, I. D. (2009). Recent Developments in the Methods and Applications of the Bond Valence Model. *Chemical Reviews*, 109(12), 6858-6919.
- Materials Project: https://materialsproject.org/
