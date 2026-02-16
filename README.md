# ScreenedBondValence

A Python package for fitting bond valence parameters $R_0$ and $B$ for cation-anion pairs using crystal structure data from the [Materials Project](https://materialsproject.org/).

Bond valence analysis relates bond lengths to bond strengths through the equation:

$$S_{ij} = \exp\left(\frac{R_0 - R_{ij}}{B}\right)$$

where $S_{ij}$ is the bond valence, $R_{ij}$ is the bond length, and $R_0$ and $B$ are empirical parameters specific to each cation-anion pair. This package provides data-driven refinement of these parameters by:

1. Fetching crystal structure data from the Materials Project API
2. Computing theoretical bond valences using network equations (valence sum rule and Kirchhoff's laws)
3. Optimizing $R_0$ and $B$ to match computed and empirical values using multiple global optimization algorithms

## Prerequisites

- **Python 3.9+**
- A **Materials Project API key** (free). Register at [materialsproject.org](https://materialsproject.org/) and get your key from [your dashboard](https://materialsproject.org/api#api-key).

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/MWhitta/ScreenedBondValence.git
cd ScreenedBondValence
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv SBV
source SBV/bin/activate    # macOS/Linux
# SBV\Scripts\activate     # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

## Usage

### As a Python script

```python
from bond_valence_processor import BondValenceProcessor

cations = ['Li']  # list of cation species to process
anions = ['O']    # list of anion species to process
my_api_key = "your_api_key"  # your Materials Project API key

# optimization algorithms to use
algos = ['shgo', 'brute', 'diff', 'dual_annealing', 'direct']

processor = BondValenceProcessor(my_api_key, algos, cations, anions)

for cation in cations:
    for anion in anions:
        processor.process_cation_system(cation, anion)
```

### As a Jupyter Notebook

See [BVparams_fit_example.ipynb](BVparams_fit_example.ipynb) for an interactive walkthrough. To run it:

```bash
jupyter notebook BVparams_fit_example.ipynb
```

**Important:** Replace `"your_api_key"` in the notebook with your actual Materials Project API key before running.

## Optimization Algorithms

The package supports five global optimization algorithms from `scipy.optimize`:

| Algorithm | Key | Description |
|---|---|---|
| SHGO | `shgo` | Simplicial Homology Global Optimization |
| Brute Force | `brute` | Grid search over parameter space |
| Differential Evolution | `diff` | Stochastic population-based optimizer |
| Dual Annealing | `dual_annealing` | Combines classical simulated annealing with fast simulated annealing |
| DIRECT | `direct` | DIviding RECTangles algorithm |

## Output

Results are saved in a `res/` directory organized by cation-anion system:

```
res/
└── LiO/
    ├── params/
    │   └── dict_matID_possible_species.json
    ├── R0Bs/
    │   ├── shgo/
    │   ├── brute/
    │   ├── diff/
    │   ├── dual_annealing/
    │   └── direct/
    ├── no_solu/
    ├── dict_sijs.json
    └── dict_charges.json
```

- **`R0Bs/<algorithm>/`** - Optimized $R_0$ and $B$ values per material, per algorithm
- **`no_solu/`** - Materials where parameter fitting failed
- **`dict_sijs.json`** - Computed theoretical bond valences
- **`dict_charges.json`** - Element charge assignments

## Project Structure

```
ScreenedBondValence/
├── bond_valence_processor.py    # Main processor: API calls, pipeline orchestration
├── BVparams_search.py           # Core: theoretical bond valence solver & parameter optimizer
├── BVparams_fit_example.ipynb   # Example Jupyter notebook
├── element2charge.json          # Default element-to-oxidation-state mapping
├── 1011191_aspod.cif            # Sample CIF file (Spodumene, LiAlSi2O6)
├── requirements.txt             # Python dependencies
└── README.md
```

## References

- Brown, I. D. (2009). Recent Developments in the Methods and Applications of the Bond Valence Model. *Chemical Reviews*, 109(12), 6858-6919.
- Materials Project: https://materialsproject.org/
