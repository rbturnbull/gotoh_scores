# cython: boundscheck=False, wraparound=False
# cython: cdivision=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# cython: language_level=3

from libc cimport limits
from libc.stdlib cimport malloc, free

cpdef (float, int) maximum(float a, float b):
    if a <= b:
        return b, 1
    return a, 0


cpdef (float, int, int, int, int) scores(str string_a, str string_b,
                              float matchWeight = 0,
                              float mismatchWeight = -1,
                              float gapWeight = -1,
                              float spaceWeight = -1):
    """
    Calculate the affine gap distance between two strings 
    
    Default weights are from Alvaro Monge and Charles Elkan, 1996, 
    "The field matching problem: Algorithms and applications" 
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.23.9685
    """

    cdef int length1 = len(string_a)
    cdef int length2 = len(string_b)

    if (string_a == string_b and
        matchWeight == max(matchWeight,
                           mismatchWeight,
                           gapWeight)):
        return matchWeight * length1, length1, 0, 0, 0

    if length1 < length2 :
        string_a, string_b = string_b, string_a
        length1, length2 = length2, length1

    # Initialize C Arrays
    cdef int memory_size = sizeof(float) * (length1+1)
    cdef float *D = <float*> malloc(memory_size)
    cdef float *V_current = <float*> malloc(memory_size)
    cdef int *V_current_matches    = <int*> malloc(memory_size)
    cdef int *V_current_mismatches = <int*> malloc(memory_size)
    cdef int *V_current_gaps       = <int*> malloc(memory_size)
    cdef int *V_current_extensions = <int*> malloc(memory_size)
    
    cdef float *V_previous = <float*> malloc(memory_size)
    cdef int *V_previous_matches    = <int*> malloc(memory_size)
    cdef int *V_previous_mismatches = <int*> malloc(memory_size)
    cdef int *V_previous_gaps       = <int*> malloc(memory_size)
    cdef int *V_previous_extensions = <int*> malloc(memory_size)

    cdef int i, j
    cdef float distance
    cdef int matches, mismatches, gaps, extensions
    
    cdef int M_choice, I_choice, D_choice

    # Set up Recurrence relations
    #
    # Base conditions
    # V(0,0) = 0
    # V(0,j) = gapWeight + spaceWeight * i
    # D(0,j) = Infinity
    V_current[0] = 0
    for j in range(1, length1 + 1) :
        V_current[j] = gapWeight + spaceWeight * j
        D[j] = -limits.INT_MAX # Should this be FLOAT_MAX??

    for i in range(1, length2 +1) :
        char2 = string_b[i-1]
        # V_previous = V_current
        for _ in range(0, length1 + 1) :
            V_previous[_] = V_current[_]
            V_previous_matches[_] = V_current_matches[_]
            V_previous_mismatches[_] = V_current_mismatches[_]
            V_previous_gaps[_] = V_current_gaps[_]
            V_previous_extensions[_] = V_current_extensions[_]

        # Base conditions    
        # V(i,0) = gapWeight + spaceWeight * i
        # I(i,0) = Infinity 
        V_current[0] = gapWeight + spaceWeight * (i-1)
        I = -limits.INT_MAX
    
        for j in range(1, length1+1) :
            char1 = string_a[j-1]

            # I(i,j) is the edit distance if the jth character of string 1
            # was inserted into string 2.
            #
            # I(i,j) = min(I(i,j-1), V(i,j-1) + gapWeight) + spaceWeight
            I, I_choice = maximum(I + spaceWeight, V_current[j-1] + gapWeight) 
        
            # D(i,j) is the edit distance if the ith character of string 2
            # was deleted from string 1
            #
            # D(i,j) = min(D(i-1,j), V(i-1,j) + gapWeight) + spaceWeight
            D[j], D_choice = maximum(D[j] + spaceWeight, V_previous[j] + gapWeight)
                    
            # M(i,j) is the edit distance if the ith and jth characters
            # match or mismatch
            #
            # M(i,j) = V(i-1,j-1) + (matchWeight | misMatchWeight)    
            if char2 == char1 :
                M = V_previous[j-1] + matchWeight
                M_choice = 1
            else:
                M = V_previous[j-1] + mismatchWeight
                M_choice = 0
            
            # V(i,j) is the minimum edit distance 
            #    
            # V(i,j) = min(E(i,j), F(i,j), G(i,j))
            if I >= D[j] and I >= M:
                V_current[j] = I
                V_current_matches[j]    = V_current_matches[j-1]
                V_current_mismatches[j] = V_current_mismatches[j-1]
                V_current_gaps[j]       = V_current_gaps[j-1] + I_choice
                V_current_extensions[j] = V_current_extensions[j-1] + 1 - I_choice
            elif D[j] >= I and D[j] >= M:
                V_current[j]            = D[j]
                V_current_matches[j]    = V_previous_matches[j]
                V_current_mismatches[j] = V_previous_mismatches[j]
                V_current_gaps[j]       = V_previous_gaps[j] + D_choice
                V_current_extensions[j] = V_previous_extensions[j] + 1 - D_choice
            else:
                V_current[j]            = M
                V_current_matches[j]    = V_previous_matches[j-1] + M_choice
                V_current_mismatches[j] = V_previous_mismatches[j-1] + 1 - M_choice
                V_current_gaps[j]       = V_previous_gaps[j-1]
                V_current_extensions[j] = V_previous_extensions[j-1]
            

    distance   = V_current[length1]
    matches    = V_current_matches[length1]
    mismatches = V_current_mismatches[length1]
    gaps       = V_current_gaps[length1]
    extensions = V_current_extensions[length1]

    free(D)
    free(V_current)
    free(V_current_matches)
    free(V_current_mismatches)
    free(V_current_gaps)
    free(V_current_extensions)
    
    free(V_previous)
    free(V_previous_matches)
    free(V_previous_mismatches)
    free(V_previous_gaps)
    free(V_previous_extensions)
    
    return distance, matches, mismatches, gaps, extensions