# Importing Dependencies
import pytest
from align import read_fasta, NeedlemanWunsch
import numpy as np

def test_nw_alignment():
    """
    Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    nw.align(seq1, seq2)

    expected_align = np.array([[  0, -10, -11, -12], 
                               [-10,   5,  -5,  -6], 
                               [-11,  -5,   4,  -6], 
                               [-12,  -6,   0,   5], 
                               [-13,  -7,  -5,   5]
    ])

    expected_gapA = np.array([[  0, -np.inf, -np.inf, -np.inf],
                              [-10, -20, -21, -22],
                              [-11,  -5, -15, -16],
                              [-12,  -6,  -6, -16],
                              [-13,  -7,  -7,  -5]
    ])

    expected_gapB = np.array([[  0, -10, -11, -12],
                              [-np.inf, -20,  -5,  -6],
                              [-np.inf, -21, -15,  -6],
                              [-np.inf, -22, -16, -10],
                              [-np.inf, -23, -17, -15]
    ])

    assert np.allclose(nw._align_matrix, expected_align)
    assert np.allclose(nw._gapA_matrix, expected_gapA)
    assert np.allclose(nw._gapB_matrix, expected_gapB)
    

def test_nw_backtrace():
    """
    Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, alignA, alignB = nw.align(seq3, seq4)

    expected_score = 18.0
    expected_alignA = "MAVHQLIRRP"
    expected_alignB = "M---QLIRHP"

    assert score == expected_score
    assert alignA == expected_alignA
    assert alignB == expected_alignB



