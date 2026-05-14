"""
Tests for sourmash-plugin-repeat-robust-mutation-rate-estimators.
"""
import os
import random
import pytest

import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_fasta(path, name, sequence):
    with open(path, 'w') as f:
        f.write(f'>{name}\n{sequence}\n')


def random_dna(length, seed=42):
    random.seed(seed)
    return ''.join(random.choice('ACGT') for _ in range(length))


def mutate(sequence, r, seed=1):
    random.seed(seed)
    bases = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
    return ''.join(
        c if random.random() > r else random.choice(bases.get(c, 'ACG'))
        for c in sequence
    )


def sketch_pair(runtmp, k=21, scaled=10, L=10000):
    """Create s and t FASTA files and sketch them in all modes."""
    s_seq = random_dna(L + k - 1, seed=42)
    t_seq = mutate(s_seq, 0.05, seed=7)

    s_fa = runtmp.output('s.fa')
    t_fa = runtmp.output('t.fa')
    make_fasta(s_fa, 's', s_seq)
    make_fasta(t_fa, 't', t_seq)

    s_std   = runtmp.output('s_standard.sig')
    s_ext   = runtmp.output('s_extended.sig')
    t_std   = runtmp.output('t_standard.sig')
    t_multi = runtmp.output('t_multiplicity.sig')

    runtmp.sourmash('scripts', 'sketch', s_fa, '--sketch-mode', 'standard',
                    '-o', s_std, '-k', str(k), '--scaled', str(scaled))
    runtmp.sourmash('scripts', 'sketch', s_fa, '--sketch-mode', 'extended',
                    '-o', s_ext, '-k', str(k), '--scaled', str(scaled))
    runtmp.sourmash('scripts', 'sketch', t_fa, '--sketch-mode', 'standard',
                    '-o', t_std, '-k', str(k), '--scaled', str(scaled))
    runtmp.sourmash('scripts', 'sketch', t_fa, '--sketch-mode', 'multiplicity',
                    '-o', t_multi, '-k', str(k), '--scaled', str(scaled))

    return s_std, s_ext, t_std, t_multi


# ---------------------------------------------------------------------------
# Basic CLI tests
# ---------------------------------------------------------------------------

def test_sketch_help(runtmp):
    """sketch command shows help."""
    runtmp.sourmash('scripts', 'sketch', '--help')
    out = runtmp.last_result.out + runtmp.last_result.err
    assert '--sketch-mode' in out
    assert 'standard' in out
    assert 'multiplicity' in out
    assert 'extended' in out


def test_mutation_rate_help(runtmp):
    """mutation_rate command shows help."""
    runtmp.sourmash('scripts', 'mutation_rate', '--help')
    out = runtmp.last_result.out + runtmp.last_result.err
    assert '--estimator' in out
    assert 'r_pp' in out
    assert 'r_pc' in out
    assert 'r_cc' in out


# ---------------------------------------------------------------------------
# Sketch tests
# ---------------------------------------------------------------------------

def test_sketch_standard(runtmp):
    """sketch --sketch-mode standard produces a .sig file."""
    s_fa = runtmp.output('s.fa')
    s_sig = runtmp.output('s.sig')
    make_fasta(s_fa, 's', random_dna(10000))

    runtmp.sourmash('scripts', 'sketch', s_fa,
                    '--sketch-mode', 'standard',
                    '-o', s_sig, '-k', '21', '--scaled', '10')

    assert os.path.exists(s_sig)
    assert runtmp.last_result.status == 0


def test_sketch_multiplicity(runtmp):
    """sketch --sketch-mode multiplicity produces a .sig file."""
    t_fa = runtmp.output('t.fa')
    t_sig = runtmp.output('t.sig')
    make_fasta(t_fa, 't', random_dna(10000, seed=99))

    runtmp.sourmash('scripts', 'sketch', t_fa,
                    '--sketch-mode', 'multiplicity',
                    '-o', t_sig, '-k', '21', '--scaled', '10')

    assert os.path.exists(t_sig)
    assert runtmp.last_result.status == 0


def test_sketch_extended(runtmp):
    """sketch --sketch-mode extended produces a .sig file."""
    s_fa = runtmp.output('s.fa')
    s_sig = runtmp.output('s.sig')
    make_fasta(s_fa, 's', random_dna(10000))

    runtmp.sourmash('scripts', 'sketch', s_fa,
                    '--sketch-mode', 'extended',
                    '-o', s_sig, '-k', '21', '--scaled', '10')

    assert os.path.exists(s_sig)
    assert runtmp.last_result.status == 0


def test_sketch_missing_mode(runtmp):
    """sketch without --sketch-mode fails."""
    s_fa = runtmp.output('s.fa')
    make_fasta(s_fa, 's', random_dna(10000))

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('scripts', 'sketch', s_fa,
                        '-o', runtmp.output('s.sig'), '-k', '21', '--scaled', '10')


# ---------------------------------------------------------------------------
# Estimator tests (just verify commands run successfully)
# ---------------------------------------------------------------------------

def test_r_pp(runtmp):
    """r_pp estimator runs successfully."""
    s_std, s_ext, t_std, t_multi = sketch_pair(runtmp)

    runtmp.sourmash('scripts', 'mutation_rate',
                    '--estimator', 'r_pp',
                    '--s-sig', s_std, '--t-sig', t_std)

    assert runtmp.last_result.status == 0
    out = runtmp.last_result.out + runtmp.last_result.err
    assert 'Estimated mutation rate' in out


def test_r_pc(runtmp):
    """r_pc estimator runs successfully."""
    s_std, s_ext, t_std, t_multi = sketch_pair(runtmp)

    runtmp.sourmash('scripts', 'mutation_rate',
                    '--estimator', 'r_pc',
                    '--s-sig', s_std, '--t-sig', t_multi)

    assert runtmp.last_result.status == 0
    out = runtmp.last_result.out + runtmp.last_result.err
    assert 'Estimated mutation rate' in out


def test_r_cc(runtmp):
    """r_cc estimator runs successfully."""
    s_std, s_ext, t_std, t_multi = sketch_pair(runtmp)

    runtmp.sourmash('scripts', 'mutation_rate',
                    '--estimator', 'r_cc',
                    '--s-sig', s_ext, '--t-sig', t_multi)

    assert runtmp.last_result.status == 0
    out = runtmp.last_result.out + runtmp.last_result.err
    assert 'Estimated mutation rate' in out


# ---------------------------------------------------------------------------
# Error handling tests
# ---------------------------------------------------------------------------

def test_mode_mismatch_r_cc(runtmp):
    """r_cc with wrong s sketch mode gives a clear error."""
    s_std, s_ext, t_std, t_multi = sketch_pair(runtmp)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('scripts', 'mutation_rate',
                        '--estimator', 'r_cc',
                        '--s-sig', s_std, '--t-sig', t_multi)

    assert 'ERROR' in runtmp.last_result.err


def test_mode_mismatch_r_pc(runtmp):
    """r_pc with wrong t sketch mode gives a clear error."""
    s_std, s_ext, t_std, t_multi = sketch_pair(runtmp)

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('scripts', 'mutation_rate',
                        '--estimator', 'r_pc',
                        '--s-sig', s_std, '--t-sig', t_std)

    assert 'ERROR' in runtmp.last_result.err


def test_k_mismatch(runtmp):
    """Mismatched k between s and t gives a clear error."""
    s_fa = runtmp.output('s.fa')
    t_fa = runtmp.output('t.fa')
    make_fasta(s_fa, 's', random_dna(10000))
    make_fasta(t_fa, 't', random_dna(10000, seed=99))

    s_sig = runtmp.output('s.sig')
    t_sig = runtmp.output('t.sig')

    runtmp.sourmash('scripts', 'sketch', s_fa, '--sketch-mode', 'standard',
                    '-o', s_sig, '-k', '21', '--scaled', '10')
    runtmp.sourmash('scripts', 'sketch', t_fa, '--sketch-mode', 'standard',
                    '-o', t_sig, '-k', '31', '--scaled', '10')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('scripts', 'mutation_rate',
                        '--estimator', 'r_pp',
                        '--s-sig', s_sig, '--t-sig', t_sig)

    assert 'ERROR' in runtmp.last_result.err