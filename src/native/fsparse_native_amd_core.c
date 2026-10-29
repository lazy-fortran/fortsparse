/* MIT: project-owned approximate minimum degree ordering.
 * Algorithm: P. R. Amestoy, T. A. Davis and I. S. Duff, "An approximate
 * minimum degree ordering algorithm", SIAM J. Matrix Anal. Appl. 17(4),
 * 886-905, 1996: quotient graph, approximate external degrees, element and
 * aggressive absorption, mass elimination, supervariable detection and dense
 * row deferral. Independent implementation; it uses only caller arrays so the
 * bounded adapter accounts for every word without SuiteSparse. */
#include "fsparse_native_amd_core.h"
#include <limits.h>

#define EMPTY (-1)
/* Encodes a node as a value <= -2: absorbed parents and elimination steps. */
#define FLIP(x) (-(x) - 2)

enum { CORE_OK = 0, CORE_INVALID = 1, CORE_BROKEN = 2 };

struct buckets { int32_t *head, *next, *last, *degree; };

static void bucket_insert(const struct buckets *b, int32_t i, int32_t deg) {
    int32_t first = b->head[deg];
    b->next[i] = first;
    b->last[i] = EMPTY;
    if (first != EMPTY) b->last[first] = i;
    b->head[deg] = i;
    b->degree[i] = deg;
}

static void bucket_remove(const struct buckets *b, int32_t i) {
    int32_t after = b->next[i], before = b->last[i];
    if (after != EMPTY) b->last[after] = before;
    if (before != EMPTY) b->next[before] = after;
    else b->head[b->degree[i]] = after;
}

static int32_t reset_flags(int32_t n, int32_t *w) {
    for (int32_t x = 0; x < n; ++x)
        if (w[x] != 0) w[x] = 1;
    return 2;
}

static int32_t root_floor(int32_t n) {
    int32_t r = 0;
    while ((int64_t)(r + 1) * (r + 1) <= n) ++r;
    return r;
}

/* Moves every live list (pe >= 0, len > 0) to the front of iw, followed by
 * the partial new element iw[pme1..pfree). The first word of each list is
 * replaced by its flipped owner so the sweep needs no extra storage. */
static int32_t compact(int32_t n, int32_t *pe, int32_t *iw, const int32_t *len,
                       int32_t *pme1, int32_t pfree) {
    for (int32_t j = 0; j < n; ++j) {
        int32_t p = pe[j];
        if (p >= 0 && len[j] > 0) {
            pe[j] = iw[p];
            iw[p] = FLIP(j);
        }
    }
    int32_t src = 0, dst = 0;
    while (src < *pme1) {
        int32_t x = iw[src++];
        if (x >= 0) continue;
        int32_t j = FLIP(x);
        iw[dst] = pe[j];
        pe[j] = dst++;
        for (int32_t k = 1; k < len[j]; ++k) iw[dst++] = iw[src++];
    }
    int32_t start = dst;
    for (src = *pme1; src < pfree; ++src) iw[dst++] = iw[src];
    *pme1 = start;
    return dst;
}

/* Heap sort of the new element's variables by their stored list hash. */
static void sift(int32_t *a, const int32_t *key, int32_t root, int32_t count) {
    int32_t value = a[root];
    for (;;) {
        int32_t child = 2 * root + 1;
        if (child >= count) break;
        if (child + 1 < count && key[a[child + 1]] > key[a[child]]) ++child;
        if (key[a[child]] <= key[value]) break;
        a[root] = a[child];
        root = child;
    }
    a[root] = value;
}

static void sort_by_hash(int32_t *a, int32_t count, const int32_t *key) {
    for (int32_t root = count / 2 - 1; root >= 0; --root) sift(a, key, root, count);
    for (int32_t end = count - 1; end > 0; --end) {
        int32_t top = a[0];
        a[0] = a[end];
        a[end] = top;
        sift(a, key, 0, end);
    }
}

/* Variables with identical element and variable lists are merged into one
 * supervariable; equal hashes only nominate candidates. */
static void detect_supervariables(int32_t *iw, int32_t pme1, int32_t pme2,
                                  int32_t *pe, const int32_t *len, int32_t *elen,
                                  int32_t *nv, const int32_t *hash, int32_t *w,
                                  int32_t *wflg) {
    sort_by_hash(iw + pme1, pme2 - pme1, hash);
    for (int32_t g = pme1, gend; g < pme2; g = gend) {
        gend = g + 1;
        while (gend < pme2 && hash[iw[gend]] == hash[iw[g]]) ++gend;
        if (gend - g < 2) continue;
        for (int32_t a = g; a < gend - 1; ++a) {
            int32_t i = iw[a];
            if (nv[i] == 0) continue;
            for (int32_t p = pe[i]; p < pe[i] + len[i]; ++p) w[iw[p]] = *wflg;
            for (int32_t b = a + 1; b < gend; ++b) {
                int32_t j = iw[b], same = 1;
                if (nv[j] == 0 || len[j] != len[i] || elen[j] != elen[i]) continue;
                for (int32_t p = pe[j]; p < pe[j] + len[j] && same; ++p)
                    same = w[iw[p]] == *wflg;
                if (!same) continue;
                nv[i] += nv[j];
                nv[j] = 0;
                elen[j] = EMPTY;
                pe[j] = FLIP(i);
            }
            ++*wflg;
        }
    }
}

int fsparse_amd_core(int32_t n, int32_t *pe, int32_t *iw, int32_t *len,
                     int32_t iwlen, int32_t pfree, int32_t *nv, int32_t *next,
                     int32_t *last, int32_t *head, int32_t *elen,
                     int32_t *degree, int32_t *w) {
    if (n <= 0 || n >= INT32_MAX / 4 || pfree < 0 ||
        (int64_t)pfree + n > iwlen)
        return CORE_INVALID;
    struct buckets b = { head, next, last, degree };
    /* Rows denser than this are ordered last, as recommended by the paper. */
    int32_t dense = 10 * root_floor(n);
    if (dense < 16) dense = 16;
    const int32_t wbig = INT32_MAX - 2 * n - 2;
    int32_t wflg = 2, mindeg = 0, nel = 0, ndense = 0, nstep = 0;

    for (int32_t i = 0; i < n; ++i) {
        head[i] = EMPTY;
        nv[i] = 1;
        w[i] = 1;
        elen[i] = 0;
    }
    for (int32_t i = 0; i < n; ++i) {
        if (len[i] < 0 || pe[i] < 0 || (int64_t)pe[i] + len[i] > pfree)
            return CORE_INVALID;
        if (len[i] > dense && len[i] > 0) {
            nv[i] = 0;
            elen[i] = EMPTY;
            pe[i] = EMPTY;
            ++ndense;
        } else {
            bucket_insert(&b, i, len[i] < n ? len[i] : n - 1);
        }
    }

    while (nel < n - ndense) {
        int32_t deg = mindeg;
        while (deg < n && head[deg] == EMPTY) ++deg;
        if (deg >= n) return CORE_BROKEN;
        mindeg = deg;
        int32_t me = head[deg];
        bucket_remove(&b, me);
        int32_t elenme = elen[me], nvpiv = nv[me];
        nel += nvpiv;
        nv[me] = -nvpiv;
        int32_t degme = 0, pme1, pme2;

        /* Construct the new element Lme from me's variables and elements. */
        if (elenme == 0) {
            int32_t q = pe[me];
            for (int32_t p = pe[me]; p < pe[me] + len[me]; ++p) {
                int32_t i = iw[p];
                if (nv[i] <= 0) continue;
                degme += nv[i];
                nv[i] = -nv[i];
                iw[q++] = i;
                bucket_remove(&b, i);
            }
            pme1 = pe[me];
            pme2 = q;
        } else {
            int32_t mp = pe[me], mrem = len[me], pending = elenme;
            pme1 = pfree;
            for (;;) {
                int32_t e, sp, srem;
                if (pending > 0) {
                    e = iw[mp++];
                    --mrem;
                    --pending;
                    if (w[e] == 0) continue;
                    sp = pe[e];
                    srem = len[e];
                } else {
                    e = me;
                    sp = mp;
                    srem = mrem;
                }
                for (; srem > 0; ++sp, --srem) {
                    int32_t i = iw[sp];
                    if (nv[i] <= 0) continue;
                    if (pfree >= iwlen) {
                        if (e == me) {
                            mp = sp;
                            mrem = srem;
                        } else {
                            pe[e] = sp;
                            len[e] = srem;
                        }
                        pe[me] = mp;
                        len[me] = mrem;
                        pfree = compact(n, pe, iw, len, &pme1, pfree);
                        if (pfree >= iwlen) return CORE_BROKEN;
                        mp = pe[me];
                        sp = e == me ? mp : pe[e];
                    }
                    degme += nv[i];
                    nv[i] = -nv[i];
                    iw[pfree++] = i;
                    bucket_remove(&b, i);
                }
                if (e == me) break;
                pe[e] = FLIP(me);
                w[e] = 0;
            }
            pme2 = pfree;
        }
        elen[me] = FLIP(nstep);
        ++nstep;
        if (wflg > wbig) wflg = reset_flags(n, w);

        /* w[e]-wflg becomes |Le \ Lme| for every element next to Lme. */
        int32_t lemax = 0;
        for (int32_t p = pme1; p < pme2; ++p) {
            int32_t i = iw[p], nvi = -nv[i];
            for (int32_t q = pe[i]; q < pe[i] + elen[i]; ++q) {
                int32_t e = iw[q], we = w[e];
                if (we >= wflg) {
                    we -= nvi;
                } else if (we != 0) {
                    we = degree[e] + wflg - nvi;
                    if (degree[e] > lemax) lemax = degree[e];
                }
                w[e] = we;
            }
        }

        /* Prune each list, bound its external degree, absorb covered elements. */
        for (int32_t p = pme1; p < pme2; ++p) {
            int32_t i = iw[p];
            if (nv[i] >= 0) continue;
            int32_t nvi = -nv[i], p1 = pe[i], p2 = p1 + elen[i], pend = p1 + len[i];
            int32_t pn = p1, ext = 0;
            uint64_t hash = 0;
            for (int32_t q = p1; q < p2; ++q) {
                int32_t e = iw[q], we = w[e];
                if (we == 0) continue;
                if (we > wflg) {
                    ext += we - wflg;
                    iw[pn++] = e;
                    hash += (uint64_t)e;
                } else {
                    pe[e] = FLIP(me);
                    w[e] = 0;
                }
            }
            int32_t p3 = pn;
            for (int32_t q = p2; q < pend; ++q) {
                int32_t j = iw[q];
                if (nv[j] <= 0) continue;
                ext += nv[j];
                iw[pn++] = j;
                hash += (uint64_t)j;
            }
            if (p3 == p1 && pn == p1) {
                /* Only adjacent to me: eliminate together with the pivot. */
                pe[i] = FLIP(me);
                nvpiv += nvi;
                nel += nvi;
                degme -= nvi;
                nv[i] = 0;
                elen[i] = EMPTY;
                continue;
            }
            if (pn >= pend) return CORE_BROKEN;
            if (ext < degree[i]) degree[i] = ext;
            iw[pn] = iw[p3];
            iw[p3] = iw[p1];
            iw[p1] = me;
            elen[i] = p3 - p1 + 1;
            len[i] = pn - p1 + 1;
            last[i] = (int32_t)(hash % (uint64_t)n);
        }
        wflg += lemax + 1;

        int32_t q = pme1;
        for (int32_t p = pme1; p < pme2; ++p)
            if (nv[iw[p]] < 0) iw[q++] = iw[p];
        pme2 = q;
        detect_supervariables(iw, pme1, pme2, pe, len, elen, nv, last, w, &wflg);

        /* Final approximate degrees; Lme keeps only principal variables. */
        q = pme1;
        for (int32_t p = pme1; p < pme2; ++p) {
            int32_t i = iw[p], nvi = -nv[i];
            if (nvi <= 0) continue;
            nv[i] = nvi;
            int32_t d = degree[i] + degme - nvi, cap = n - nel - nvi;
            if (d > cap) d = cap;
            if (d < 0) d = 0;
            bucket_insert(&b, i, d);
            if (d < mindeg) mindeg = d;
            iw[q++] = i;
        }
        nv[me] = nvpiv;
        degree[me] = degme;
        len[me] = q - pme1;
        if (len[me] == 0) {
            pe[me] = EMPTY;
            w[me] = 0;
        } else {
            pe[me] = pme1;
            w[me] = 1;
        }
        if (elenme != 0) pfree = q;
    }

    /* Pivot blocks in elimination order, then the deferred dense rows. */
    for (int32_t s = 0; s < n; ++s) head[s] = EMPTY;
    for (int32_t x = 0; x < n; ++x)
        if (elen[x] <= -2) head[FLIP(elen[x])] = x;
    int32_t pos = 0;
    for (int32_t s = 0; s < nstep; ++s) {
        int32_t e = head[s];
        if (e == EMPTY) return CORE_BROKEN;
        degree[e] = pos;
        pos += nv[e];
    }
    if (pos != n - ndense) return CORE_BROKEN;
    for (int32_t j = 0; j < n; ++j) {
        int32_t x = j;
        while (elen[x] == EMPTY && pe[x] != EMPTY) x = FLIP(pe[x]);
        int32_t k;
        if (elen[x] == EMPTY) {
            k = pos++;
        } else {
            k = degree[x]++;
            for (int32_t y = j; y != x;) {
                int32_t up = FLIP(pe[y]);
                pe[y] = FLIP(x);
                y = up;
            }
        }
        if (k < 0 || k >= n) return CORE_BROKEN;
        last[k] = j;
        next[j] = k;
    }
    return CORE_OK;
}
