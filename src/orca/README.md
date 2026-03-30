# 🐋 ORCA — Optimized aRbitrary-ensemble Compression Algorithm

> **ORCA** compresses a trained Random Forest using a Genetic Algorithm that searches for the optimal `BitVector` mask — minimising the number of active bits (smaller forest) while preserving classification accuracy.

---

## How It Works

A trained Random Forest of `N` trees is represented by a `BitVector` of length `3 × N`, split into three segments:

```
[ b_1 … b_N | b_(N+1) … b_(2N) | b_(2N+1) … b_(3N) ]
  presence     depth/pruning       alphabet
```

The Genetic Algorithm evolves this bitvector by minimising a combined fitness function:

```
fitness = error_rate + penalty_weight × (n_ones / n_bits)
```

where `error_rate` is measured on a held-out **validation set** and `n_ones / n_bits` is the compression ratio. The result is a compressed forest that trades a small accuracy loss for a significant reduction in model size.

### Pipeline Stages

Once the best bitvector is found, it is applied through three sequential compression stages:

| Stage         | Bits used       | What it does                                                      |
|---------------|-----------------|-------------------------------------------------------------------|
| `first_part`  | bits `1..N`     | **Tree selection** — drops trees with bit = 0                     |
| `second_part` | bits `N+1..2N`  | **Tree pruning** — cuts trees exceeding the reference depth       |
| `third_part`  | bits `2N+1..3N` | **Alphabet unification** — aligns feature thresholds across trees |

---
## Usage

### Given a forest, find the optimal compression mask

```
compression(original_f, mode, f_val, l_val)
```

### Custom parameters

```julia
compression(
    population_size = 100,     # GA population — more = better exploration, slower
    n_generations   = 200,     # GA iterations — more = better convergence, slower
    penalty_weight  = 0.3      # Trade-off weight: higher = more compression, lower = more accuracy
)
```

### Quick test run

```julia
compression(population_size = 10, n_generations = 10)
```

---

## Parameters

### `compression` keyword arguments

| Parameter         | Default | Description                                       |
|-------------------|---------|---------------------------------------------------|
| `population_size` | `50`    | Number of bitvector individuals per GA generation |
| `n_generations`   | `100`   | Number of GA iterations                           |
| `penalty_weight`  | `0.3`   | Compression vs accuracy trade-off (see below)     |

### Tuning `penalty_weight`

`penalty_weight` controls the balance between compression and accuracy:

```
penalty_weight = 0.0  →  ignore compression, maximise accuracy only
penalty_weight = 0.3  →  balanced (default)
penalty_weight = 1.0  →  aggressive compression, accept accuracy loss
```

To explore the full trade-off frontier, run multiple times with different values:

```julia
for w in [0.1, 0.2, 0.3, 0.5, 0.8]
    println("=== penalty_weight = $w ===")
    compression(penalty_weight = w)
end
```
---

## 📦 Dependencies

```julia
using DecisionTree    # RF training and raw tree manipulation
using SoleModels      # DecisionEnsemble representation and apply()
using SoleData        # Dataset utilities
using Evolutionary    # Genetic Algorithm (GA, tournament, UX)
using DataFrames      # DataFrame conversion for SoleModels.apply
using StatsBase       # mode() for majority class voting
```

---

## ⚠️ Notes

- **Validation vs Test**: the GA optimises on the validation set (`15%` of data). The test set (`10%`) is never seen during optimisation and is the honest accuracy measure.
- **Scalarised objective**: ORCA uses a weighted-sum fitness — it is not true multi-objective optimisation. Each run yields one point on the Pareto frontier. To explore the full frontier, vary `penalty_weight`.
- **Reproducibility**: training uses a fixed seed (`2025`). GA results may vary across runs due to random initialisation of the population.

---

## 📄 License

MIT