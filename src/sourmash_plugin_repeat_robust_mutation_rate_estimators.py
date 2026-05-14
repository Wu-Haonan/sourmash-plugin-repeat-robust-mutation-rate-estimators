"""
Repeat-robust mutation rate estimators for sourmash.

Three estimators: r_pp, r_pc, r_cc.
Each requires different sketch modes for s (source) and t (query).

Sketch modes:
  standard     - distinct k-mer hashes + L
                 use as s or t for r_pp
  multiplicity - k-mer hashes with per-hash count + L
                 use as t for r_pc and r_cc
  extended     - distinct k-mer hashes + sum_occ_h1 + L
                 use as s for r_cc

Valid combinations:
  r_pp: s=standard,    t=standard
  r_pc: s=standard,    t=multiplicity
  r_cc: s=extended,    t=multiplicity
"""

usage_sketch = """
   sourmash scripts sketch input.fa --sketch-mode <standard|multiplicity|extended> -o output.sig
"""

usage_mutation_rate = """
   sourmash scripts mutation_rate --estimator <r_pp|r_pc|r_cc> --s-sig S.sig --t-sig T.sig
"""

epilog = """
Note: our three repeat-robust estimators treat the roles of the two input sequences asymmetrically.
In general, estimators that use more information from s and t achieve higher accuracy.

Sketch modes and valid combinations:
  r_pp : s sketched with --sketch-mode standard,      t sketched with --sketch-mode standard
  r_pc : s sketched with --sketch-mode standard,      t sketched with --sketch-mode multiplicity
  r_cc : s sketched with --sketch-mode extended,      t sketched with --sketch-mode multiplicity

If unsure which sequence should be s or t:
  s = the longer sequence
  t = the shorter sequence

To learn more about our estimators, see:
https://www.biorxiv.org/content/10.64898/2026.04.01.715966v1.full

Need help?
http://github.com/sourmash-bio/sourmash/issues
"""

import sys
import json
import argparse
from collections import defaultdict

import sourmash
from sourmash.plugins import CommandLinePlugin
from sourmash.logging import notify, error


# ---------------------------------------------------------------------------
# Math
# ---------------------------------------------------------------------------

def _r_from_q(q, k):
    q = max(0.0, min(float(q), 1.0))
    if q == 0.0:
        return 0.0
    return 1.0 - (1.0 - q) ** (1.0 / k)


def estimate_r_pp(diff, L_t, k):
    return _r_from_q(diff / L_t, k)


def estimate_r_pc(diff_multi, L_t, k):
    return _r_from_q(diff_multi / L_t, k)


def estimate_r_cc(diff_multi, L_t, k, sum_occ_h1_normalized):
    r_pc = estimate_r_pc(diff_multi, L_t, k)
    correction = sum_occ_h1_normalized * ((1.0 - r_pc) ** (k - 1)) * r_pc / 3.0
    q_cc = diff_multi / L_t + correction
    return _r_from_q(q_cc, k)


# ---------------------------------------------------------------------------
# h1 computation
# ---------------------------------------------------------------------------

def compute_h1_values(kmer_set, k):
    mask_to_kmers = defaultdict(set)
    for kmer in kmer_set:
        for pos in range(k):
            masked = kmer[:pos] + 'N' + kmer[pos + 1:]
            mask_to_kmers[masked].add(kmer)
    h1_map = {}
    for kmer in kmer_set:
        count = 0
        for pos in range(k):
            masked = kmer[:pos] + 'N' + kmer[pos + 1:]
            count += len(mask_to_kmers[masked]) - 1
        h1_map[kmer] = count
    return h1_map


def compute_sum_occ_h1(sequence, k):
    occ = defaultdict(int)
    seq = sequence.upper()
    for i in range(len(seq) - k + 1):
        occ[seq[i:i + k]] += 1
    h1_map = compute_h1_values(set(occ.keys()), k)
    return sum(occ[kmer] * h1_map[kmer] for kmer in occ)


# ---------------------------------------------------------------------------
# Sketch functions
# ---------------------------------------------------------------------------

def _scaled_threshold(scaled):
    MAXHASH = (1 << 64) - 1
    return MAXHASH // scaled


def _hash_kmer(kmer):
    return sourmash.minhash.hash_murmur(kmer)


def _iter_kmers(sequence, k):
    seq = sequence.upper()
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]


def sketch_standard(sequence, k, scaled):
    threshold = _scaled_threshold(scaled)
    hashes = set()
    for kmer in _iter_kmers(sequence, k):
        h = _hash_kmer(kmer)
        if h <= threshold:
            hashes.add(h)
    return {
        'mode': 'standard',
        'k': k, 'scaled': scaled,
        'L': len(sequence) - k + 1,
        'hashes': hashes,
    }


def sketch_multiplicity(sequence, k, scaled):
    threshold = _scaled_threshold(scaled)
    counts = defaultdict(int)
    for kmer in _iter_kmers(sequence, k):
        h = _hash_kmer(kmer)
        if h <= threshold:
            counts[h] += 1
    return {
        'mode': 'multiplicity',
        'k': k, 'scaled': scaled,
        'L': len(sequence) - k + 1,
        'counts': dict(counts),
    }


def sketch_extended(sequence, k, scaled):
    std = sketch_standard(sequence, k, scaled)
    L_s = std['L']
    sum_occ_h1 = compute_sum_occ_h1(sequence, k)
    return {
        'mode': 'extended',
        'k': k, 'scaled': scaled,
        'L': L_s,
        'hashes': std['hashes'],
        'sum_occ_h1': sum_occ_h1 / L_s, 
    }


SKETCH_MODE_HELP = {
    'standard':
        'store distinct k-mer hashes and L '
        '[used as s or t for r_pp], \n'
        'where L = |x| - k + 1 is the total number of k-mers in string x',

    'multiplicity':
        'store k-mer hashes, counts, and L '
        '[used as t for r_pc and r_cc]',

    'extended':
        'store distinct k-mer hashes, sum_occ_h1, and L '
        '[used as s for r_cc], \n'
        'where sum_occ_h1 is a constant computed from string s for bias correction',
}


# ---------------------------------------------------------------------------
# Signature save/load
# ---------------------------------------------------------------------------

def save_sketch_sig(sketch_dict, output_path):
    d = dict(sketch_dict)
    if 'hashes' in d:
        d['hashes'] = sorted(d['hashes'])
    if 'counts' in d:
        d['counts'] = {str(h): v for h, v in d['counts'].items()}

    mh = sourmash.MinHash(n=0, ksize=d['k'], scaled=d['scaled'])
    if 'hashes' in d:
        mh.add_many(d['hashes'])
    elif 'counts' in d:
        mh.add_many([int(h) for h in d['counts'].keys()])

    ss = sourmash.SourmashSignature(mh)
    ss._name = json.dumps({
        'mode': d['mode'],
        'L': d['L'],
        'sum_occ_h1': d.get('sum_occ_h1', None),
        'counts': d.get('counts', None),
    })

    with sourmash.sourmash_args.SaveSignaturesToLocation(output_path) as save_sig:
        save_sig.add(ss)


def load_sketch_sig(sig_path):
    siglist = list(sourmash.load_file_as_signatures(sig_path))
    if not siglist:
        raise ValueError(f"No signatures found in {sig_path}")
    ss = siglist[0]
    meta = json.loads(ss.name)
    mh = ss.minhash

    d = {
        'mode': meta['mode'],
        'k': mh.ksize,
        'scaled': mh.scaled,
        'L': meta['L'],
    }

    if meta['mode'] in ('standard', 'extended'):
        d['hashes'] = set(mh.hashes.keys())
        if meta.get('sum_occ_h1') is not None:
            d['sum_occ_h1'] = meta['sum_occ_h1']
    elif meta['mode'] == 'multiplicity':
        d['hashes'] = set(int(h) for h in meta['counts'].keys())
        d['counts'] = {int(h): c for h, c in meta['counts'].items()}

    return d


# ---------------------------------------------------------------------------
# Diff computation
# ---------------------------------------------------------------------------

def compute_diffs(sk_s, sk_t):
    scaled = sk_s['scaled']
    L_t = sk_t['L']          
    s_hashes = sk_s['hashes']
    t_counts = sk_t['counts']

    diff = 0
    diff_multi = 0
    for h, count in t_counts.items():
        if h not in s_hashes:
            diff += 1
            diff_multi += count

    return min(diff * scaled, L_t), min(diff_multi * scaled, L_t)


# ---------------------------------------------------------------------------
# CLI: sourmash scripts sketch
# ---------------------------------------------------------------------------

class Command_Sketch(CommandLinePlugin):
    command = 'sketch'
    description = 'Sketch a sequence for mutation rate estimation'
    usage = usage_sketch
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        subparser.add_argument('input_fasta',
                               help='input FASTA file')
        subparser.add_argument('-o', '--output', required=True,
                               help='output signature file (.sig)')
        subparser.add_argument('--sketch-mode',
                               choices=['standard', 'multiplicity', 'extended'],
                               required=True,
                               help='\n'.join(f'  {m}: {h}'
                                              for m, h in SKETCH_MODE_HELP.items()))
        subparser.add_argument('-k', '--ksize', type=int, default=21,
                               help='k-mer size (default: 21)')
        subparser.add_argument('--scaled', type=int, default=1000,
                               help='scaled factor for FracMinHash (default: 1000)')

    def main(self, args):
        super().main(args)
        import screed

        notify(f"Reading sequence from {args.input_fasta}")
        sequence = ''
        with screed.open(args.input_fasta) as f:
            for record in f:
                sequence += record.sequence.upper()

        notify(f"Sketching with --sketch-mode={args.sketch_mode}, k={args.ksize}, scaled={args.scaled}")

        if args.sketch_mode == 'standard':
            sk = sketch_standard(sequence, args.ksize, args.scaled)
        elif args.sketch_mode == 'multiplicity':
            sk = sketch_multiplicity(sequence, args.ksize, args.scaled)
        elif args.sketch_mode == 'extended':
            notify("Computing sum_occ_h1 (may take a moment for large sequences)...")
            sk = sketch_extended(sequence, args.ksize, args.scaled)

        notify(f"Saving signature to {args.output}")
        save_sketch_sig(sk, args.output)
        notify("Done.")


# ---------------------------------------------------------------------------
# CLI: sourmash scripts mutation_rate
# ---------------------------------------------------------------------------

VALID_COMBINATIONS = {
    'r_pp': ('standard',    'standard'),
    'r_pc': ('standard',    'multiplicity'),
    'r_cc': ('extended',    'multiplicity'),
}


class Command_MutationRate(CommandLinePlugin):
    command = 'mutation_rate'
    description = 'Estimate mutation rate between two sequences using r_pp, r_pc, or r_cc'
    usage = usage_mutation_rate
    epilog = epilog
    formatter_class = argparse.RawTextHelpFormatter

    def __init__(self, subparser):
        super().__init__(subparser)
        subparser.add_argument('--estimator', required=True,
                               choices=['r_pp', 'r_pc', 'r_cc'],
                               help='estimator to use:\n'
                                    '  r_pp: requires s=standard,    t=standard\n'
                                    '  r_pc: requires s=standard,    t=multiplicity\n'
                                    '  r_cc: requires s=extended,    t=multiplicity')
        subparser.add_argument('--s-sig', required=True,
                               help='signature for s (source/reference sequence)')
        subparser.add_argument('--t-sig', required=True,
                               help='signature for t (query sequence)')

    def main(self, args):
        super().main(args)

        notify(f"Loading s signature from {args.s_sig}")
        sk_s = load_sketch_sig(args.s_sig)
        notify(f"Loading t signature from {args.t_sig}")
        sk_t = load_sketch_sig(args.t_sig)

        # Validate modes
        s_mode = sk_s['mode']
        t_mode = sk_t['mode']
        required_s, required_t = VALID_COMBINATIONS[args.estimator]

        if s_mode != required_s or t_mode != required_t:
            error(f"ERROR: estimator {args.estimator} requires "
                  f"s sketched with --sketch-mode {required_s}, "
                  f"t sketched with --sketch-mode {required_t}")
            error(f"  but got s mode='{s_mode}', t mode='{t_mode}'")
            error(f"  Please re-sketch with the correct --sketch-mode.")
            sys.exit(1)

        # Validate k and scaled match
        if sk_s['k'] != sk_t['k']:
            error(f"ERROR: k mismatch: s has k={sk_s['k']}, t has k={sk_t['k']}")
            sys.exit(1)
        if sk_s['scaled'] != sk_t['scaled']:
            error(f"ERROR: scaled mismatch: s scaled={sk_s['scaled']}, t scaled={sk_t['scaled']}")
            sys.exit(1)

        k = sk_s['k']
        L_s = sk_s['L']
        L_t = sk_t['L']

        # r_pp: t is standard (no counts), treat each hash as count=1
        if args.estimator == 'r_pp':
            sk_t['counts'] = {h: 1 for h in sk_t['hashes']}

        diff, diff_multi = compute_diffs(sk_s, sk_t)

        if args.estimator == 'r_pp':
            result = estimate_r_pp(diff, L_t, k)
        elif args.estimator == 'r_pc':
            result = estimate_r_pc(diff_multi, L_t, k)
        elif args.estimator == 'r_cc':
            result = estimate_r_cc(diff_multi, L_t, k, sk_s['sum_occ_h1'])

        notify(f"\nEstimator : {args.estimator}")
        notify(f"k         : {k}")
        notify(f"scaled    : {sk_s['scaled']}")
        notify(f"L_s       : {L_s}")
        notify(f"Estimated mutation rate : {result:.6f}")