from libc.math cimport sqrt, pow, fabs
import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

INT = np.int
ctypedef np.int_t INT_t

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
def distance_matrix(np.ndarray[DTYPE_t, ndim=2] X):

    cdef:
        int L = X.shape[0]
        int N = X.shape[1]
        int i, j
        double dist
        np.ndarray[DTYPE_t, ndim=2] result = np.empty((L, L), dtype=DTYPE)

    for i in xrange(L):
        for j in xrange(i + 1, L):
            dist = 0
            for k in xrange(N):
                dist += pow(X[i, k] - X[j, k], 2)
            dist = sqrt(dist)
            result[i, i] = 0.0
            result[i, j] = dist
            result[j, i] = dist
    return result

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def similarity_matrix(np.ndarray[DTYPE_t, ndim=2] X,
                      np.ndarray[DTYPE_t, ndim=2] Y):

    cdef:
        double winsize = 8.0
        double sumsize = (winsize - 1.0) * (winsize - 2.0) / 2.0
        int i, j, row, col
        int Lx = X.shape[0]
        int Ly = Y.shape[0]
        double score
        np.ndarray[DTYPE_t, ndim=2] S = np.empty((Lx, Ly), dtype=DTYPE)

    for i in range(Lx):
        for j in range(Ly):
            S[i, j] = -1.0
            if i > Lx - winsize or j > Ly - winsize:
                continue

            score = 0.0

            for row in range(<long> winsize - 2):
                for col in range(<long> winsize):
                    score += fabs(X[i + row, i + col] - Y[j + row, j + col])

            S[i, j] = score / sumsize
    return S


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def find_path(np.ndarray[DTYPE_t, ndim=2] S,
              np.ndarray[DTYPE_t, ndim=2] dA,
              np.ndarray[DTYPE_t, ndim=2] dB):
    cdef:
        double winsize = 8.0
        double D0 = 3.0
        double D1 = 4.0
        int MAX_KEPT = 20
        int GAP_MAX = 30

        int iA, iB, i, g, s, k, jA, jB, jgap, gA, gB
        int lenA = S.shape[0]
        int lenB = S.shape[1]
        double best_path_score = 1e6
        int best_path_len = 0
        int smaller = S.shape[0] if S.shape[0] < S.shape[1] else S.shape[1]
        int winsum = 21

        np.ndarray[INT_t, ndim=2] bestpath = -1 * np.ones((smaller, 2), dtype=INT)

        int bufferindex = 0
        int buffersize = 0

        np.ndarray[INT_t, ndim=1] lenbuffer = np.zeros(MAX_KEPT, dtype=INT)
        np.ndarray[DTYPE_t, ndim=1] scorebuffer = 1e6 * np.ones(MAX_KEPT, dtype=DTYPE)
        np.ndarray[INT_t, ndim=3] pathbuffer = np.empty((MAX_KEPT, smaller, 2), dtype=INT)
        np.ndarray[INT_t, ndim=1] wincache = np.array(
            [(i * 1) * i * winsize / 2 + (i + 1) * winsum for i in range(smaller)],
            dtype=INT
        )
        np.ndarray[DTYPE_t, ndim=2] all_score_buffer = 1e6 * np.ones(
            (smaller, 2 * (GAP_MAX) + 1), dtype=DTYPE
        )
        np.ndarray[INT_t, ndim=1] tindex = np.empty(smaller, dtype=INT)
        np.ndarray[INT_t, ndim=2] cur_path
        int gap_best_index = -1
        int cur_path_len, done
        double cur_total_score, gap_best_score, cur_score, score1, score2
        np.ndarray[INT_t, ndim=3] results
        np.ndarray[INT_t, ndim=2] p
        np.ndarray[INT_t, ndim=1] r

    for iA in range(lenA):
        if iA > lenA - winsize * (best_path_len - 1):
            break
        for iB in range(lenB):
            if S[iA, iB] >= D0:
                continue
            if S[iA, iB] == -1:
                continue
            if iB > lenB - winsize * (best_path_len - 1):
                break

            cur_path = -1 * np.ones((smaller, 2), dtype=INT)

            cur_path[0, 0] = iA
            cur_path[0, 1] = iB
            cur_path_len = 1
            tindex[0] = 0
            cur_total_score = 0.0

            done = 0

            while not done:
                gap_best_score = 1e6
                gap_best_index = -1

                for g in range(2 * GAP_MAX + 1):
                    jA = cur_path[cur_path_len - 1, 0] + <int> winsize
                    jB = cur_path[cur_path_len - 1, 1] + <int> winsize

                    if (g + 1) % 2 == 0:
                        jA += (g + 1) / 2
                    else:
                        jB += (g + 1) / 2

                    if jA > lenA - winsize - 1 or jB > lenB - winsize - 1:
                        continue
                    if S[jA, jB] > D0:
                        continue
                    if S[jA, jB] == -1.0:
                        continue

                    cur_score = 0.0
                    for s in range(cur_path_len):
                        cur_score += fabs(
                            dA[cur_path[s, 0], jA] - dB[cur_path[s, 1], jB]
                        )
                        cur_score += fabs(
                            dA[cur_path[s, 0] + (<int> winsize - 1), jA + (<int> winsize - 1)] -
                            dB[cur_path[s, 1] + (<int> winsize - 1), jB + (<int> winsize - 1)]
                        )
                        for k in range(1, <long> winsize - 1):
                            cur_score += fabs(
                                dA[cur_path[s, 0] + k, jA + (<int> winsize - 1) - k] -
                                dB[cur_path[s, 1] + k, jB + (<int> winsize - 1) - k]
                            )
                    cur_score /= winsize * <double> cur_path_len

                    if cur_score >= D1:
                        continue

                    if cur_score < gap_best_score:
                        cur_path[cur_path_len, 0] = jA
                        cur_path[cur_path_len, 1] = jB
                        gap_best_score = cur_score
                        gap_best_index = g
                        all_score_buffer[cur_path_len - 1, g] = cur_score

                cur_total_score = 0.0

                score1 = 0.0
                score2 = 0.0

                if gap_best_index != -1:
                    jgap = (gap_best_index + 1) / 2
                    if (gap_best_index + 1) % 2 == 0:
                        gA = cur_path[cur_path_len - 1, 0] + <int> winsize + jgap
                        gB = cur_path[cur_path_len - 1, 1] + <int> winsize
                    else:
                        gA = cur_path[cur_path_len - 1, 0] + <int> winsize
                        gB = cur_path[cur_path_len - 1, 1] + <int> winsize + jgap

                    score1 = (
                        all_score_buffer[cur_path_len - 1, gap_best_index] * winsize *
                        cur_path_len + S[gA, gB] * winsum
                    ) / (winsize * cur_path_len + winsum)

                    if cur_path_len > 1:
                        score2 = all_score_buffer[cur_path_len - 2, tindex[cur_path_len - 1]]
                    else:
                        score2 = S[iA, iB]
                    score2 *= (
                        wincache[cur_path_len - 1] + score1 * (wincache[cur_path_len] -
                                                               wincache[cur_path_len - 1])
                    )
                    score2 /= wincache[cur_path_len]
                    cur_total_score = score2
                    if cur_total_score > D1:
                        done = 1
                        gap_best_index = -1
                        break
                    else:
                        all_score_buffer[cur_path_len - 1, gap_best_index] = cur_total_score
                        tindex[cur_path_len] = gap_best_index
                        cur_path_len += 1
                else:
                    done = 1
                    cur_path_len -= 1
                    break

                if cur_path_len > best_path_len or (
                    cur_path_len == best_path_len and cur_total_score < best_path_score
                ):
                    best_path_len = cur_path_len
                    best_path_score = cur_total_score
                    bestpath = np.copy(cur_path)

            if best_path_len > lenbuffer[bufferindex] or (
                best_path_len == lenbuffer[bufferindex] and best_path_score < scorebuffer[bufferindex]
            ):
                bufferindex = 0 if bufferindex == MAX_KEPT - 1 else bufferindex + 1
                buffersize = buffersize + 1 if buffersize < MAX_KEPT else MAX_KEPT

                if bufferindex == 0 and buffersize == MAX_KEPT:
                    pathbuffer[MAX_KEPT - 1] = np.copy(bestpath)
                    scorebuffer[MAX_KEPT - 1] = best_path_score
                    lenbuffer[MAX_KEPT - 1] = best_path_len

                else:
                    pathbuffer[bufferindex - 1] = np.copy(bestpath)
                    scorebuffer[bufferindex - 1] = best_path_score
                    lenbuffer[bufferindex - 1] = best_path_len

    results = np.array([[r for r in p if r[0] != -1] for p in pathbuffer])
    return results
