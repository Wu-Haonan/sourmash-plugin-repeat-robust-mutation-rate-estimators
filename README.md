# Repeat-robust estimator for mutation rate

This repo describes a [sourmash](https://sourmash.bio/) plugin that estimates the mutation rate between two DNA sequences.
It implements the estimators from this [ISMB paper](https://doi.org/10.64898/2026.04.01.715966).
These estimators are well-suited for highly repetitive sequences such as centromeres.

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
## Quick start

To start, you need two FASTA files (`.fa` or `.fasta`) as input. Each FASTA file may contain multiple records. This tool is designed for full-length assembled sequences. We currently do not support sequencing reads as input.

You first need to sketch each of the FASTA files. There are three sketch modes: standard, multipliticy, and extended. We will describe which one to use later. As an example, if your two files are s.fa and t.fa and you want to sketch them using standard and multiplicity modes respectively, then you can run:

```
sourmash scripts sketch s.fa --sketch-mode standard    -o s.sig -k 21 --scaled 1000
sourmash scripts sketch t.fa --sketch-mode multiplicity -o t.sig -k 21 --scaled 1000
```

This uses a k-mer size of 21 and a subsampling rate of 1/1000. A smaller subsampling rate leads to a smaller sketch size and faster downstream estimation, but introduces more variance in downstream estimation.

Once you've sketched the files, you are ready to compare them. There various estimators which we describe below. For example, to compare them using the `r_pc` estimator use

```
sourmash scripts mutation_rate --estimator r_pc --s-sig s.sig --t-sig t.sig
```

This will generate output that looks like

```
Estimator : r_pc
k         : 21
scaled    : 1000
L_s       : 4800000
Estimated mutation rate : 0.012345
```


---

## Choice of $s$ and $t$ 

Our estimators are not symmetric and your decision on which fasta file to label as $s$ and which as $t$ may give different results. 
Our underlying model assumes that there is a starting string $s$ where each position has a substitution at rate $r$, resulting in the mutated string $t$.
To choose which of your FASTA files to use as $s$ and which as $t$, we give the following guidelines:
- If you know the biological direction (e.g., one sequence is ancestral), assign the ancestral sequence as $s$ and the derived sequence as $t$.
- If you are unsure, use the longer sequence as $s$ and the shorter as $t$.
- If you are comparing many sequences (e.g., building a phylogenetic tree), be consistent in the way you assign $s$ and $t$, e.g. always assign the longer sequence as $s$ across all pairs.
- If the biological direction of the sequences is unclear, we recommend running the estimator twice, swapping the roles of $s$ and $t$ between the runs. That way, you can compare the difference in the estimates to the uncertainty you are willing to tolerate. 
---

## Choice of estimator 

We provide three estimators, and each estimator requires a specific sketch mode for $s$ and $t$. 

| Estimator | Sketch mode for s | Sketch mode for t | Accuracy |
|-----------|------------------|------------------|----------|
| `r_pp`    | `standard`       | `standard`       | Lower    |
| `r_pc`    | `standard`       | `multiplicity`   | Medium   |
| `r_cc`    | `extended`       | `multiplicity`   | Highest  |

- `r_pp` (presence-to-presence): ignores the multiplicity of $k$-mers of both $s$ and $t$. Use this estimator if you do not trust the multiplicity information in your input sequences, e.g. you expect collapsed repeats or the input is unitigs.
- `r_pc` (presence-to-count): ignores the multiplicity of $k$-mers in $s$ but uses the multiplicity of $k$-mers in $t$. 
- `r_cc` (count-to-count): uses $k$-mer multiplicites of both $s$ and $t$, with a bias correction term precomputed from $s$. This is the most accurate estimator, especially for highly repetitive sequences. However, sketching $s$ in `extended` mode will take longer.

To estimate the mutation rate, you can use one of the following commands:

```
sourmash scripts mutation_rate --estimator r_pp --s-sig s.sig --t-sig t.sig
sourmash scripts mutation_rate --estimator r_pc --s-sig s.sig --t-sig t.sig
sourmash scripts mutation_rate --estimator r_cc --s-sig s.sig --t-sig t.sig
```

If the estimator finds that the sketch modes do not match its requirements, it will report an error and tell you which modes are needed. 

---

## Support

Please file issues at https://github.com/Wu-Haonan/sourmash-plugin-repeat-robust-mutation-rate-estimators/issues

---

## Citation

If you use this plugin, please cite:

> Wu, H. and Medvedev, P. (2026). The gift of novelty: repeat-robust *k*-mer-based estimators of mutation rates. Proceedings of ISMB 2026, to appear in *Bionformatics*.

## License

MIT License. See [LICENSE](LICENSE) for details.
