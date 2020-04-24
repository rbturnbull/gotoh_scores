# cython: boundscheck=False, wraparound=False
# cython: cdivision=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# cython: language_level=3
import numpy as np
cimport numpy as np

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
cdef extern from "Python.h":
    ctypedef void PyObject
    PyObject *PyString_FromStringAndSize(char *, size_t)
    int _PyString_Resize(PyObject **, size_t)
    char * PyUnicode_AsUTF8(PyObject *)


cimport cython

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.int8_t DTYPE_BOOL

cdef size_t UP = 1, LEFT = 2, DIAG = 3, NONE = 4


cpdef (float, int, int, int, int) scores(str seqi, str seqj,
                              float match = 0,
                              float mismatch = -1,
                              float gap_open = -1,
                              float gap_extend = -1):
    """
    perform a global sequence alignment (needleman-wunsch) on seq and and 2. using
    the matrix for nucleotide transition from matrix if available.
    where matrix is of the format provided in the ncbi/data directory.

    >>> from nwalign import global_align
    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')
    
    # See https://bitbucket.org/brentp/biostuff/src/default/

    """

    cdef bint flip = 0
    
    cdef size_t max_j = len(seqj) 
    cdef size_t max_i = len(seqi) 
    if max_i == max_j == 0:
        return 0, 0, 0, 0, 0


    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i

    # need to initialize j for the case when it's a zero-length string.
    cdef size_t i = 0, j = 0, p
    cdef float diag_score, up_score, left_score, tscore

    cdef char ci, cj
    cdef int zero=0, one=1
    cdef int match_count=0, mismatch_count=0, gap_count=0, extension_count=0

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_i = np.ones((max_i + 1,), dtype=np.int8)
    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_j = np.ones((max_j + 1,), dtype=np.int8)

    cdef np.ndarray[float, ndim=2] score = np.empty((max_i + 1, max_j + 1), dtype=np.float32)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.empty((max_i + 1, max_j + 1), dtype=np.uint)
  
    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.float32)
    score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.float32)

    agap_i[0] = zero

    for i in range(1, max_i + 1):
        ci = seqi[<size_t>(i - 1)]

        agap_j[0] = zero
        for j in range(1, max_j + 1):
            agap_j[j] = one
            cj = seqj[<size_t>(j - 1)]
            tscore = match if ci == cj else mismatch

            diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore
            up_score   = score[<size_t>(i - 1), j] + (gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend)
            left_score = score[i, <size_t>(j - 1)] + (gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend)

            """
            fix cases where scores are tied.
            choose diagonal when not at the ends of either string.
            and choose up/left (gap) when at the end of the string.
            this forces gaps to the ends of the aligments.
            """
            if diag_score == left_score:
                # so here. if we're at the end we choose a gap.
                # otherwise, we choose diag
                if i == max_i or i == 1:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            elif diag_score == up_score:
                if j == max_j or j == 1:
                    score[i, j] = up_score
                    pointer[i, j] = UP
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            # end of ambiguous score checks.

            elif up_score > diag_score: #checked.
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
            elif diag_score > up_score:
                if left_score > diag_score:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero

    cdef float score_max = score[-1,-1]

    p = pointer[i, j]
    cdef size_t previous_pointer = NONE
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            if seqj[j] == seqi[i]:
                match_count += 1
            else:
                mismatch_count += 1
        elif p == LEFT:
            j -= 1
            if p == previous_pointer:
                extension_count += 1
            else:
                gap_count += 1
            
        elif p == UP:
            i -= 1
            if p == previous_pointer:
                extension_count += 1
            else:
                gap_count += 1
        else:
            raise Exception('Bad Pointer!:pointer: %i', p)
        previous_pointer = p
        p = pointer[i, j]
        
        
    #print(score_max, match_count, mismatch_count, gap_count, extension_count)
    return score_max, match_count, mismatch_count, gap_count, extension_count
    
    
@cython.boundscheck(False)
@cython.nonecheck(False)
cpdef align(str seqi, str seqj, int match, int mismatch, int gap_open, int gap_extend):
    """
    perform a global sequence alignment (needleman-wunsch) on seq and and 2. using
    the matrix for nucleotide transition from matrix if available.
    where matrix is of the format provided in the ncbi/data directory.

    >>> from nwalign import global_align
    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')

    """
    cdef bint flip = 0
    
    cdef size_t max_j = len(seqj) 
    cdef size_t max_i = len(seqi) 
    if max_i == max_j == 0:
        return "", ""


    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i

    # need to initialize j for the case when it's a zero-length string.
    cdef size_t i = 0, j = 0, seqlen, align_counter = 0, p
    cdef int diag_score, up_score, left_score, tscore

    cdef char *align_j, *align_i
    cdef char ci, cj
    cdef str ai, aj
    cdef int zero=0, one=1

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_i = np.ones((max_i + 1,), dtype=np.int8)
    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_j = np.ones((max_j + 1,), dtype=np.int8)

    cdef np.ndarray[DTYPE_INT, ndim=2] score = np.empty((max_i + 1, max_j + 1), dtype=np.int)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.empty((max_i + 1, max_j + 1), dtype=np.uint)

  
    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.int)
    score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.int)

    agap_i[0] = zero

    for i in range(1, max_i + 1):
        ci = seqi[<size_t>(i - 1)]

        agap_j[0] = zero
        for j in range(1, max_j + 1):
            agap_j[j] = one
            cj = seqj[<size_t>(j - 1)]
            tscore = match if ci == cj else mismatch

            diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore
            up_score   = score[<size_t>(i - 1), j] + (gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend)
            left_score = score[i, <size_t>(j - 1)] + (gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend)

            """
            fix cases where scores are tied.
            choose diagonal when not at the ends of either string.
            and choose up/left (gap) when at the end of the string.
            this forces gaps to the ends of the aligments.
            """
            if diag_score == left_score:
                # so here. if we're at the end we choose a gap.
                # otherwise, we choose diag
                if i == max_i or i == 1:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            elif diag_score == up_score:
                if j == max_j or j == 1:
                    score[i, j] = up_score
                    pointer[i, j] = UP
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            # end of ambiguous score checks.

            elif up_score > diag_score: #checked.
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
            elif diag_score > up_score:
                if left_score > diag_score:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero

    seqlen = max_i + max_j
    ai = ""
    aj = ""

    cdef int score_max #, = score[:, -1].max()


    p = pointer[i, j]
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            aj += seqj[j]
            ai += seqi[i]
        elif p == LEFT:
            j -= 1
            aj += seqj[j]
            ai += "-"
        elif p == UP:
            i -= 1
            aj += "-"
            ai += seqi[i]
        else:
            raise Exception('wtf!:pointer: %i', p)
        align_counter += 1
        p = pointer[i, j]

    #_PyString_Resize(&aj, align_counter)
    #_PyString_Resize(&ai, align_counter)
    if flip:
        return (<object>ai)[::-1], (<object>aj)[::-1] #, score.max()
    else:
        return (<object>aj)[::-1], (<object>ai)[::-1] #, score.max()