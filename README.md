# sourmash-plugin-repeat-robust-mutation-rate-estimators

This is a [sourmash](https://sourmash.bio/) plugin that estimates the **mutation rate** between two DNA sequences using repeat-robust $k$-mer-based estimators. It is particularly well-suited for highly repetitive sequences such as **centromeres**.

The method is described in:

> Wu, H. and Medvedev, P. (2026). The gift of novelty: repeat-robust *k*-mer-based estimators of mutation rates. *bioRxiv*. https://www.biorxiv.org/content/10.64898/2026.04.01.715966v1

---

## Installation

```
# Option 1: conda + pip (recommended)
conda install -c conda-forge -c bioconda sourmash
pip install sourmash-plugin-repeat-robust-mutation-rate-estimators

# Option 2: pip only
pip install sourmash
pip install sourmash-plugin-repeat-robust-mutation-rate-estimators
```

Verify the plugin is recognized:

```
sourmash scripts
```

You should see `sketch` and `mutation_rate` listed under available plugin commands.

---

## Supported Input

This tool accepts **FASTA files** (`.fa`, `.fasta`) as input. Each FASTA file may contain multiple records. This tool is designed for **full-length assembled sequences** only, such as whole genome vs. whole genome or full centromere sequence vs. full centromere sequence. We currently **Do NOT** supported:

- Sequencing reads vs. genome
- Sequencing reads vs. sequencing reads

---

## Mutation Model and Sequence Assignment ($s$ and $t$)

We consider the following substitution model, parameterized by a rate $0 \le r \le 1$. Given a string $s$, the character at each position independently mutates to one of the three other nucleotides with probability $r/3$. We denote the mutated string as $t$.

**The roles of $s$ and $t$ are not symmetric.** Swapping them may give different results, especially for repetitive sequences.

- **If you know the biological direction** (e.g., one sequence is ancestral), assign the ancestral sequence as $s$ and the derived sequence as $t$.
- **If you are unsure**, use the longer sequence as $s$ and the shorter as $t$. We cannot guarantee accuracy if the assignment is unclear. If you are comparing many sequences (e.g., building a phylogenetic tree), be consistent: always assign the longer sequence as $s$ across all pairs.

---

## Estimators and Usage

This plugin provides three estimators. Each estimator requires a specific **sketch mode** for $s$ and $t$. Choose your estimator first, then sketch accordingly.

| Estimator | Sketch mode for s | Sketch mode for t | Accuracy |
|-----------|------------------|------------------|----------|
| `r_pp`    | `standard`       | `standard`       | Lower    |
| `r_pc`    | `standard`       | `multiplicity`   | Medium   |
| `r_cc`    | `extended`       | `multiplicity`   | Highest  |

- **r_pp** (presence-to-presence): uses distinct $k$-mers of both $s$ and $t$.
- **r_pc** (presence-to-count): uses distinct $k$-mers of $s$ and $k$-mer counts of $t$.
- **r_cc** (count-to-count): uses $k$-mer counts of both $s$ and $t$, with a bias correction term precomputed from $s$. Most accurate for repetitive sequences. Note: sketching $s$ with `extended` mode may take longer for large genomes.

### Step 1: Sketch your sequences

Choose the commands that match your chosen estimator:

**For r_pp** — sketch both sequences with `standard`:
```
sourmash scripts sketch s.fa --sketch-mode standard -o s.sig -k 21 --scaled 1000
sourmash scripts sketch t.fa --sketch-mode standard -o t.sig -k 21 --scaled 1000
```

**For r_pc** — sketch $s$ with `standard`, sketch $t$ with `multiplicity`:
```
sourmash scripts sketch s.fa --sketch-mode standard    -o s.sig -k 21 --scaled 1000
sourmash scripts sketch t.fa --sketch-mode multiplicity -o t.sig -k 21 --scaled 1000
```

**For r_cc** — sketch $s$ with `extended`, sketch $t$ with `multiplicity`:
```
sourmash scripts sketch s.fa --sketch-mode extended    -o s.sig -k 21 --scaled 1000
sourmash scripts sketch t.fa --sketch-mode multiplicity -o t.sig -k 21 --scaled 1000
```

Parameters:
- `-k`: $k$-mer size.
- `--scaled`: subsampling rate. Larger value leads to smaller sketch, faster but less precise.

### Step 2: Estimate the mutation rate

```
sourmash scripts mutation_rate --estimator r_pp --s-sig s.sig --t-sig t.sig
sourmash scripts mutation_rate --estimator r_pc --s-sig s.sig --t-sig t.sig
sourmash scripts mutation_rate --estimator r_cc --s-sig s.sig --t-sig t.sig
```

If the sketch modes do not match the estimator requirements, the tool will report an error and tell you which modes are needed.

Example output:
```
Estimator : r_cc
k         : 21
scaled    : 1000
L_s       : 4800000
Estimated mutation rate : 0.012345
```

---

## Support

Please file issues at https://github.com/Wu-Haonan/sourmash-plugin-repeat-robust-mutation-rate-estimators/issues

---

## Citation

If you use this plugin, please cite:

> Wu, H. and Medvedev, P. (2026). The gift of novelty: repeat-robust *k*-mer-based estimators of mutation rates. *bioRxiv*. https://www.biorxiv.org/content/10.64898/2026.04.01.715966v1

## License

MIT License. See [LICENSE](LICENSE) for details.