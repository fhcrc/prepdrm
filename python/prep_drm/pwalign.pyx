# distutils: language = c++
# distutils: include_dirs = seqan-library-1.4.1/include prep_drm
# distutils: sources = prep_drm/simplepwalign.cpp
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdlib cimport free, malloc
from libc.string cimport strdup
from libc.stdint cimport int32_t, uint16_t, uint8_t

cdef extern from "simplepwalign.hpp" namespace "simplepwalign":

    cdef cppclass PairwiseAlignment:
        string aligned_ref
        string aligned_qry
        int score

    PairwiseAlignment pairwise_global(const string& ref, const string& qry,
                                      const int match_score, const int mismatch_score,
                                      const int gap_open, const int gap_extend,
                                      const bool boundary_gaps_qry)


cdef int BAM_CMATCH = 0
cdef int BAM_CINS = 1
cdef int BAM_CDEL = 2

CIGAR_OPS = ['M', 'I', 'D']

def pw_global(bytes ref, bytes qry, int match_score=4, int mismatch_score=-2, int gap_open=-4, int gap_extend=-2,
              bool boundary_gaps_qry=False):
    cdef string r = ref, q = qry
    cdef PairwiseAlignment result = pairwise_global(r, q, match_score,
                                                    mismatch_score, gap_open,
                                                    gap_extend,
                                                    boundary_gaps_qry)
    return result.aligned_ref, result.aligned_qry, result.score

def calculate_cigar(bytes aligned_ref not None, bytes aligned_query not None):
    assert len(aligned_ref) == len(aligned_query)
    if len(aligned_ref) == 0:
        return []

    cdef int last_op = -1, new_op = -1, n_ops = 0
    cdef size_t i
    cdef list cigar_ops = []
    for i in xrange(len(aligned_ref)):
        if aligned_ref[i] == '-':
            assert aligned_query[i] != '-'
            new_op = BAM_CINS
        elif aligned_query[i] == '-':
            assert aligned_ref[i] != '-'
            new_op = BAM_CDEL
        else:
            new_op = BAM_CMATCH
        if last_op == -1 or new_op != last_op:
            if last_op != -1:
                cigar_ops.append((last_op, n_ops))
            last_op = new_op
            n_ops = 1
        else:
            n_ops += 1
    cigar_ops.append((last_op, n_ops))
    return cigar_ops
