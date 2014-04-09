import collections


cdef char* OPNAMES = "MIDNSH"
cdef enum cigar_op_t:
    BAM_CMATCH = 0
    BAM_CINS = 1
    BAM_CDEL = 2
    BAM_CSKIP = 3
    BAM_CSOFTCLIP = 4
    BAM_CHARDCLIP = 5

AlignedPair = collections.namedtuple('AlignedPair', ['qpos', 'rpos',
                                                     'qbase', 'rbase',
                                                     'qual', 'cigar_op'])

cdef int consumes_query(int cigar_op):
    return cigar_op == BAM_CMATCH or cigar_op == BAM_CINS or cigar_op == BAM_CSOFTCLIP

cdef int consumes_ref(int cigar_op):
    return cigar_op == BAM_CMATCH or cigar_op == BAM_CDEL or cigar_op == BAM_CSKIP

def all_pairs_iter(read not None, bytes reference not None):
    cdef int op, cq, cr
    cdef size_t oplength, qi = 0, ri = read.pos, i
    cdef bytes seq = read.seq
    cdef bytes qual = read.qual

    for op, oplength in read.cigar:
        cq = consumes_query(op)
        cr = consumes_ref(op)
        for i in xrange(oplength):
            yield AlignedPair(qpos=qi if cq else None,
                              rpos=ri if cr else None,
                              qbase=seq[qi] if cq else None,
                              rbase=reference[ri] if cr else None,
                              qual=qual[qi] if cq else None,
                              cigar_op=op)
            if cq:
                qi += 1
            if cr:
                ri += 1
    assert qi == len(seq)

def cigar_to_pwalign(read, bytes reference not None):
    cdef int op, cq, cr
    cdef size_t oplength, qi = 0, ri = read.pos, i
    cdef bytes seq = read.seq
    cdef list qaln = [], raln = []

    for op, oplength in read.cigar:
        cq = consumes_query(op)
        cr = consumes_ref(op)
        if cq:
            qaln.append(seq[qi:qi+oplength])
            qi += oplength
        else:
            qaln.append('-' * oplength)
        if cr:
            raln.append(reference[ri:ri+oplength])
            ri += oplength
        else:
            raln.append('-' * oplength)
        if op == BAM_CSOFTCLIP:
            qaln[-1] = qaln[-1].lower()
    return ''.join(raln), ''.join(qaln)

