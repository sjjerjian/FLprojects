# pyDots3DMP

Python codes for dots3DMP experiments modelling and analysis.

### Basic Usage

#### Environment setup 🔧

- Option A — **recommended**: install `uv` and use `uv sync` to install dependencies directly from `pyproject.toml`:
  1. Install `uv` if needed:
     ```bash
     curl -LsSf https://astral.sh/uv/install.sh | sh  # macOS / Linux
     ```
     Visit [Getting Started](https://docs.astral.sh/uv/getting-started/installation/) for more information.
  2. From the project root, sync and install dependencies - uv will automatically create a venv in the project workspace.
     ```bash
     uv sync
     ```

- Option B — Standard venv creation and pip installation:
    Ensure you are using Python >= 3.12
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # macOS / Linux
    pip install -r requirements.txt
    pip install -e .  # to include the package as an editable install
    ```

---

#### Project structure 📁

Package code lives under `src/` (installed as `ddm`, `behavior`, and `neural`):

- `src/ddm/` — DDM and accumulator code: `Accumulator.py`, `selfmotionddm.py`, `moi.py` (method of images), `np_cache.py`
- `src/behavior/` — behavioral helpers and descriptive analysis: `utils.py`, `descriptive.py`
- `src/neural/` — neural data loading, tuning, and decoding modules. Not updated since late 2023.

At project root:

- `scripts/` — example and analysis scripts (e.g. `ddm_testing_script.py`, `run_selfmotion_ddm.py`, `ddm_demo.ipynb`)
- `archive/` — deprecated or legacy code

---

#### Using `SelfMotionDDM` (behavior/selfmotionddm.py) ✅

`SelfMotionDDM` implements the 3DMP accumulator model with optional confidence/wager readouts. Example workflow:

1. Instantiate the model:
```python
import numpy as np
from ddm import SelfMotionDDM

# set diffusion grid resolution
grid_vec = np.arange(-3, 0, 0.05)
time_vec = np.arange(0, 2, 0.01)

# define initial parameters
init_params = {
    'kmult': [0.6, 0.6],
    'bound': [1.0, 1.0, 1.0],
    'non_dec_time': [0.3],
    'wager_thr': [1, 1, 1],
    'wager_alpha': [0.05],
}
ddm = SelfMotionDDM(
  grid_vec=grid_vec,
  tvec=time_vec,
  **init_params,
  stim_scaling=False,
  return_wager=False
  )
```

2. Prepare data (use provided helpers in `behavior.preprocessing`):
- `data = data_cleanup(filepath)`
- `data = format_onetargconf(data, remove_one_targ=True)`

Data shapes expected:
- X DataFrame: columns `['modality', 'coherence', 'delta', 'heading']`
- y DataFrame: columns `['choice', 'PDW', 'RT']` (RT is used for RT likelihood computations)

3. Fit the model:
```python
# optionally specify parameters to hold fixed during fitting
accum.fit(X, y, fixed_params=['kmult'])
```
- `fixed_params` is an optional list of parameter names to keep constant during optimization - these do not get passed to the optimization call but are used by the `predict` method.
- By default, optimization uses `pybads.BADS` (requires `pybads`); the code also contains a hook to use `scipy.optimize.minimize`.

4. Generate predictions:
```python
y_pred, y_pred_samp = ddm.predict(X, n_samples=1, cache_accumulators=True, seed=1)
```
- `predict` returns a DataFrame with predicted likelihoods for `choice`, `PDW`, and `RT` and an optional sampled predictions DataFrame when `n_samples > 0`.

Key parameters (see `behavior/selfmotionddm.py` for defaults):
- `kmult` — k multipliers per modality (e.g., `[k_ves, k_vis]`)
- `bound` — bounds per modality (e.g., `[ves, vis, comb]`)
- `non_dec_time` — non-decision times per modality
- `wager_thr` — wager threshold(s)
- `wager_alpha` — wager mapping alpha(s)
- `return_wager` — whether to compute wager/confidence predictions
- `stim_scaling` — whether to compute stimulus-driven urgency signals (or pass tuple of urgency arrays)

---

#### Notes & roadmap

- The `Accumulator` class and helpers (`behavior/Accumulator.py`, `behavior/moi.py`) provides the low-level method-of-images computations (CDF/PDF, RT distributions, and log posterior odds) that `SelfMotionDDM` uses.
- The older `ddm_2d` codebase is deprecated; you may still find legacy code under `archive/` but active analyses should use `SelfMotionDDM` and `Accumulator`.

---

## TO DO

1. Some experimental features (cue-combination strategies, different confidence mappings) are skeletons and not fully implemented or tested.
3. add unit tests.
4. Improve and expand documentation.
5. Improve and expand diagnostic visualizations of accumulators.
6. Improve logging/saving of json results and iteration history for checks and resuming optimization runs.
