// permpoly: distribution of S(p) = sum_i E[i][p(i)] over all permutations p of [n], modulo a prime q.
// Part of permsum-oeis, (c) 2026 Victor Nidal Bassem Nagib, MIT licence (see LICENSE).
// Build: gcc -O3 -march=native -fopenmp -o permpoly src/permpoly.c
// Computes P(x) = per( x^E[i][j] ) at the N-th roots of unity with Glynn's formula (Gray code order),
// then recovers the coefficients c_S = #{p : S(p) = S} for Smin <= S <= Smax (needs Smax - Smin < N).
//
// usage: permpoly q g N Smin Smax sym T [eval:lo:hi | combine:FILE] < matrix   (matrix: n, then n*n exponents)
//   eval:lo:hi    evaluate only the points m in [lo, hi) and print "m P(omega^m)" (for splitting long runs)
//   combine:FILE  read all "m value" lines from FILE (points 0..M-1) and print the coefficients
//   q: prime with N | q-1, q < 2^62;  g: primitive root mod q;  N: power of two
//   sym=1: the distribution is palindromic, c_S = c_{T-S}; only half of the points are evaluated
// output: lines "S c_S" (c_S mod q)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
static u64 Q, QINV, R2;

static inline u64 mmul(u64 a, u64 b) {            // Montgomery product a*b/2^64 mod Q
    u128 t = (u128)a * b;
    u64 m = (u64)t * QINV;
    u64 u = (u64)((t + (u128)m * Q) >> 64);
    return u >= Q ? u - Q : u;
}
static inline u64 madd(u64 a, u64 b) { u64 s = a + b; return s >= Q ? s - Q : s; }
static inline u64 msub(u64 a, u64 b) { return a >= b ? a - b : a + Q - b; }
static u64 tomont(u64 a) { return mmul(a % Q, R2); }
static u64 frommont(u64 a) { return mmul(a, 1); }
static u64 mpow(u64 a, u64 e) {                  // a in Montgomery form
    u64 r = tomont(1);
    while (e) { if (e & 1) r = mmul(r, a); a = mmul(a, a); e >>= 1; }
    return r;
}

#define PB 8      // evaluation points processed together (independent multiplication chains)
#define MAXN 32

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: permpoly q g N Smin Smax sym T < matrix\n"); return 1; }
    Q = strtoull(argv[1], 0, 10);
    u64 g = strtoull(argv[2], 0, 10);
    long N = atol(argv[3]);
    long Smin = atol(argv[4]), Smax = atol(argv[5]);
    int sym = atoi(argv[6]);
    long T = atol(argv[7]);
    long lo = -1, hi = -1; const char *cfile = 0;
    if (argc > 8 && !strncmp(argv[8], "eval:", 5)) sscanf(argv[8] + 5, "%ld:%ld", &lo, &hi);
    if (argc > 8 && !strncmp(argv[8], "combine:", 8)) cfile = argv[8] + 8;
    if (Smax - Smin >= N) { fprintf(stderr, "N too small\n"); return 1; }
    // Montgomery constants
    u64 inv = 1;
    for (int i = 0; i < 7; i++) inv *= 2 - Q * inv;      // inverse of Q mod 2^64 (Newton)
    QINV = (u64)0 - inv;
    { u128 r = ((u128)1 << 64) % Q; R2 = (u64)((r * r) % Q); }

    int n; if (scanf("%d", &n) != 1 || n < 1 || n > MAXN) { fprintf(stderr, "bad n\n"); return 1; }
    static long E[MAXN][MAXN];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (scanf("%ld", &E[i][j]) != 1) { fprintf(stderr, "bad matrix\n"); return 1; }

    // omega = g^((Q-1)/N), powers table in Montgomery form
    u64 gm = tomont(g);
    u64 w = mpow(gm, (Q - 1) / N);
    u64 *pw = malloc(sizeof(u64) * N);
    pw[0] = tomont(1);
    for (long t = 1; t < N; t++) pw[t] = mmul(pw[t - 1], w);
    if (frommont(mmul(pw[N - 1], w)) != 1) { fprintf(stderr, "omega has wrong order\n"); return 1; }

    long M = sym ? N / 2 + 1 : N;                     // points to evaluate: m = 0 .. M-1
    u64 *val = calloc(N, sizeof(u64));                 // P(omega^m), Montgomery form
    long first = 0, last = M;
    if (lo >= 0) { first = lo < M ? lo : M; last = hi < M ? hi : M; }
    if (cfile) { first = last = 0;
        FILE *f = fopen(cfile, "r"); long m; unsigned long long v; long got = 0;
        if (!f) { fprintf(stderr, "cannot open %s\n", cfile); return 1; }
        while (fscanf(f, "%ld %llu", &m, &v) == 2) { if (m < 0 || m >= M) { fprintf(stderr, "bad m\n"); return 1; } val[m] = tomont(v); got++; }
        fclose(f); if (got != M) { fprintf(stderr, "combine: got %ld of %ld points\n", got, M); return 1; } }
    long nblk = (last - first + PB - 1) / PB;
    // 1/2^(n-1) in Montgomery form
    u64 two = tomont(2), inv2 = mpow(two, Q - 2), invh = mpow(inv2, n - 1);
    u64 steps = (u64)1 << (n - 1);

    #pragma omp parallel for schedule(dynamic, 1)
    for (long b = 0; b < nblk; b++) {
        u64 Mt[MAXN][MAXN][PB], M2[MAXN][MAXN][PB], s[MAXN][PB], acc[PB], prod[PB];
        long m0 = first + b * PB;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int p = 0; p < PB; p++) {
            long m = m0 + p; if (m >= last) m = first;
            long e = (long)(((__int128)m * E[i][j]) % N);
            Mt[i][j][p] = pw[e]; M2[i][j][p] = madd(pw[e], pw[e]);
        }
        for (int i = 0; i < n; i++) for (int p = 0; p < PB; p++) {
            u64 t = 0; for (int j = 0; j < n; j++) t = madd(t, Mt[i][j][p]); s[i][p] = t;
        }
        for (int p = 0; p < PB; p++) { u64 pr = s[0][p]; for (int i = 1; i < n; i++) pr = mmul(pr, s[i][p]); acc[p] = pr; }
        u64 delta = 0; int neg = 0;                   // bit j set: column j currently has sign -1 (column 0 fixed +1)
        for (u64 k = 1; k < steps; k++) {
            int j = __builtin_ctzll(k) + 1;
            if ((delta >> j) & 1) {
                for (int i = 0; i < n; i++) for (int p = 0; p < PB; p++) s[i][p] = madd(s[i][p], M2[i][j][p]);
                neg--;
            } else {
                for (int i = 0; i < n; i++) for (int p = 0; p < PB; p++) s[i][p] = msub(s[i][p], M2[i][j][p]);
                neg++;
            }
            delta ^= (u64)1 << j;
            for (int p = 0; p < PB; p++) prod[p] = s[0][p];
            for (int i = 1; i < n; i++) for (int p = 0; p < PB; p++) prod[p] = mmul(prod[p], s[i][p]);
            if (neg & 1) for (int p = 0; p < PB; p++) acc[p] = msub(acc[p], prod[p]);
            else         for (int p = 0; p < PB; p++) acc[p] = madd(acc[p], prod[p]);
        }
        for (int p = 0; p < PB; p++) { long m = m0 + p; if (m < last) val[m] = mmul(acc[p], invh); }
    }
    if (lo >= 0) { for (long m = first; m < last; m++) printf("%ld %llu\n", m, (unsigned long long)frommont(val[m])); return 0; }
    if (sym) {   // P(omega^(N-m)) = omega^(-m*T) * P(omega^m)
        for (long m = 1; m < N / 2; m++) {
            long e = (long)((((__int128)(N - m)) * (T % N)) % N);   // omega^(-m*T) = omega^((N-m)*T)
            val[N - m] = mmul(val[m], pw[e]);
        }
    }
    // c_S = (1/N) sum_m val[m] * omega^(-m*S)
    u64 invN = mpow(tomont(N % Q), Q - 2);
    for (long S = Smin; S <= Smax; S++) {
        long sm = ((S % N) + N) % N;
        u64 c = 0;
        for (long m = 0; m < N; m++) {
            long e = (long)((((__int128)m) * (N - sm)) % N);
            c = madd(c, mmul(val[m], pw[e]));
        }
        printf("%ld %llu\n", S, (unsigned long long)frommont(mmul(c, invN)));
    }
    return 0;
}
