# sourmash-plugin-repeat-robust-mutation-rate-estimators

[sourmash](https://sourmash.bio) is a tool for biological sequence analysis and comparisons.

This plugin implements repeat-robust substitution rate estimators r_pp, r_pc, and r_cc based on FracMinHash sketches, as described in:

> Wu, H. and Medvedev, P. (2026). Repeat-robust estimation of substitution rates from k-mer sketches. *bioRxiv*. https://www.biorxiv.org/content/10.64898/2026.04.01.715966v1

## Installation

Install sourmash, then install this plugin:

```
# Option 1:
conda install -c conda-forge -c bioconda sourmash
pip install sourmash-plugin-repeat-robust-mutation-rate-estimators

# Option 2:
pip install sourmash
pip install sourmash-plugin-repeat-robust-mutation-rate-estimators
```

Verify the plugin is recognized:

```
sourmash scripts
```

You should see `sketch` and `mutation_rate` listed under available plugin commands.

## Usage

### Background

The three estimators treat the two input sequences **asymmetrically**: we assume string t is mutated from string s.

If unsure which is s and which is t, use the longer sequence as s.

Each estimator requires a specific sketch mode:

| Estimator | s sketch mode | t sketch mode |
|-----------|--------------|--------------|
| `r_pp`    | `standard`   | `standard`   |
| `r_pc`    | `standard`   | `multiplicity` |
| `r_cc`    | `extended`   | `multiplicity` |

In general, estimators that use more information achieve higher accuracy.

### Step 1: Sketch your sequences

```
# For r_pp
sourmash scripts sketch s.fa --sketch-mode standard    -o s.sig -k 21 --scaled 1000
sourmash scripts sketch t.fa --sketch-mode standard    -o t.sig -k 21 --scaled 1000

# For r_pc
sourmash scripts sketch s.fa --sketch-mode standard    -o s.sig -k 21 --scaled 1000
sourmash scripts sketch t.fa --sketch-mode multiplicity -o t.sig -k 21 --scaled 1000

# For r_cc
sourmash scripts sketch s.fa --sketch-mode extended    -o s.sig -k 21 --scaled 1000
sourmash scripts sketch t.fa --sketch-mode multiplicity -o t.sig -k 21 --scaled 1000
```

Sketch modes:
- `standard`: stores distinct k-mer hashes and L, where L = |x| - k + 1 is the total number of k-mers in string x. Use as s or t for r_pp.
- `multiplicity`: stores k-mer hashes with per-hash counts and L. Use as t for r_pc and r_cc.
- `extended`: stores distinct k-mer hashes, L, and a precomputed correction constant `sum_occ_h1`. Use as s for r_cc. Note: computing `sum_occ_h1` requires reading the full sequence and may take longer for large genomes.

### Step 2: Estimate mutation rate

```
sourmash scripts mutation_rate --estimator r_pp --s-sig s.sig --t-sig t.sig
sourmash scripts mutation_rate --estimator r_pc --s-sig s.sig --t-sig t.sig
sourmash scripts mutation_rate --estimator r_cc --s-sig s.sig --t-sig t.sig
```

Example output:
```
Estimator : r_cc
k         : 21
scaled    : 1000
L_s       : 4800000
Estimated mutation rate : 0.012345
```

## Support

Please file issues at https://github.com/Wu-Haonan/sourmash-plugin-repeat-robust-mutation-rate-estimators/issues

## Dev docs

`sourmash-plugin-repeat-robust-mutation-rate-estimators` is developed at https://github.com/Wu-Haonan/sourmash-plugin-repeat-robust-mutation-rate-estimators.


## Citation

If you use this plugin, please cite:

```
Wu, H. and Medvedev, P. (2026). Repeat-robust estimation of substitution rates
from k-mer sketches. bioRxiv.
https://www.biorxiv.org/content/10.64898/2026.04.01.715966v1
```

## License

MIT License. See [LICENSE](LICENSE) for details.
