#!/usr/bin/env python3
"""permsum: exact distributions of S(p) = sum_k f(k, p(k)) over all n! permutations p of [n].

Kinds   kp   f = k*p(k)         A073114, and through D = sum (k - p(k))^2 = 2*(n(n+1)(2n+1)/6 - S) also
                                 A180111, A180112, A179263, A179265, A179264, A179266 (Spearman's rho)
        gcd  f = gcd(k, p(k))    A073343
        k2p  f = k^2*p(k)        A189043

Usage   python3 permsum.py compute KIND N [THREADS]       one run; writes results/distributions/KIND_N.json
        python3 permsum.py chunk SECONDS KIND:N[:PIECES]...  long runs in resumable pieces (stops starting new
                                                              pieces after SECONDS; run again to continue)
        python3 permsum.py finalize KIND:N[:PIECES] ...      combine the pieces into the distribution
        python3 permsum.py report                            all sequence terms, compared with the known OEIS terms
        python3 permsum.py bfiles                            OEIS b-files (known + new terms) in results/bfiles
        python3 permsum.py newterms                          markdown summary of the new terms (results/new_terms.md)

Method  P_n(x) = sum_p x^S(p) is the permanent of the n x n matrix (x^f(i,j)). permpoly (src/permpoly.c) evaluates
        it at the N-th roots of unity modulo a prime q = c*2^20 + 1 with Glynn's formula in Gray-code order, and
        recovers the coefficients (the counts) by an inverse transform; two primes and the Chinese remainder
        theorem give the exact counts (exact while n! < q1*q2, i.e. n <= 33).
Checks  every distribution must have integer counts in [0, n!], sum to n!, and have exactly the mean and
        variance of Hoeffding's formula for sum_i a(i, p(i)) over a uniformly random permutation.

(c) 2026 Victor Nidal Bassem Nagib. Code: MIT licence. Data: CC BY 4.0. See README.md.
"""
import json, math, os, subprocess, sys, time
from fractions import Fraction

ROOT = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(ROOT, "permpoly")
DIST = os.path.join(ROOT, "results", "distributions")
CHUNKS = os.path.join(ROOT, "results", "chunks")
PRIMES = [(4611686018405367809, 3), (4611686018326724609, 3)]     # q = c*2^20 + 1 < 2^62, primitive root g


# ---------- the statistic ----------
def matrix(kind, n):
    f = {"kp": lambda i, j: i * j, "gcd": math.gcd, "k2p": lambda i, j: i * i * j}[kind]
    return [[f(i, j) for j in range(1, n + 1)] for i in range(1, n + 1)]


def span(kind, n):
    """(Smin, Smax, T): range of S, and T with c_S = c_(T-S) when the distribution is palindromic (else None)."""
    s2, s3 = n * (n + 1) * (2 * n + 1) // 6, (n * (n + 1) // 2) ** 2
    if kind == "kp":
        return n * (n + 1) * (n + 2) // 6, s2, n * (n + 1) * (n + 1) // 2
    if kind == "gcd":
        return (n if n > 1 else 1), n * (n + 1) // 2, None
    if kind == "k2p":
        return (n + 1) * s2 - s3, s3, (n + 1) * s2
    raise ValueError(kind)


def plan(kind, n):
    """(Smin, Smax, T, N, M): N = transform size (power of 2 > Smax - Smin), M = points evaluated."""
    smin, smax, T = span(kind, n)
    N = 1
    while N < smax - smin + 1:
        N *= 2
    return smin, smax, T, N, (N // 2 + 1 if T is not None else N)


# ---------- running the engine ----------
def engine(kind, n, q, g, mode=None, threads=2):
    smin, smax, T, N, M = plan(kind, n)
    inp = f"{n}\n" + "\n".join(" ".join(map(str, row)) for row in matrix(kind, n)) + "\n"
    args = [BIN, str(q), str(g), str(N), str(smin), str(smax), "1" if T is not None else "0", str(T or 0)]
    if mode:
        args.append(mode)
    out = subprocess.run(args, input=inp, capture_output=True, text=True, check=True,
                         env=dict(os.environ, OMP_NUM_THREADS=str(threads))).stdout
    return {int(a): int(b) for a, b in (line.split() for line in out.splitlines() if line.strip())}


def crt(residues, primes):
    x, m = 0, 1
    for r, q in zip(residues, primes):
        x, m = x + m * (((r - x) * pow(m, -1, q)) % q), m * q
    return x


def combine(residue_lists):
    qs = [q for q, g in PRIMES]
    return {S: crt([r[S] for r in residue_lists], qs) for S in residue_lists[0]}


# ---------- checks ----------
def exact_moments(a):
    """Mean and variance of sum_i a[i][p(i)] for a uniformly random permutation p (Hoeffding 1951)."""
    n = len(a)
    tot = sum(map(sum, a))
    if n == 1:
        return Fraction(tot), Fraction(0)
    r = [Fraction(sum(row), n) for row in a]
    c = [Fraction(sum(a[i][j] for i in range(n)), n) for j in range(n)]
    g = Fraction(tot, n * n)
    return Fraction(tot, n), sum((a[i][j] - r[i] - c[j] + g) ** 2 for i in range(n) for j in range(n)) / (n - 1)


def check(kind, n, counts):
    f = math.factorial(n)
    if not all(0 <= v <= f for v in counts.values()): return "count outside [0, n!]"
    if sum(counts.values()) != f: return "counts do not sum to n!"
    mean, var = exact_moments(matrix(kind, n))
    m1 = Fraction(sum(S * v for S, v in counts.items()), f)
    m2 = Fraction(sum(S * S * v for S, v in counts.items()), f)
    if m1 != mean: return "mean differs from the exact value"
    if m2 - m1 * m1 != var: return "variance differs from the exact value"
    return None


def save(kind, n, counts, seconds):
    os.makedirs(DIST, exist_ok=True)
    problem = check(kind, n, counts)
    json.dump({"kind": kind, "n": n, "seconds": round(seconds, 1), "check": problem or "passed",
               "counts": {str(S): c for S, c in sorted(counts.items())}},
              open(os.path.join(DIST, f"{kind}_{n}.json"), "w"))
    print(f"{kind} n={n}: {'checks passed' if not problem else 'CHECK FAILED: ' + problem} [{seconds:.0f} s]", flush=True)
    return problem


def load(kind, n):
    d = json.load(open(os.path.join(DIST, f"{kind}_{n}.json")))
    return {int(S): c for S, c in d["counts"].items()}


# ---------- sequences ----------
def is_prime(m):
    if m < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if m % p == 0: return m == p
    d, s = m - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):     # deterministic for m < 3.1e23
        x = pow(a, d, m)
        if x in (1, m - 1): continue
        for _ in range(s - 1):
            x = x * x % m
            if x == m - 1: break
        else: return False
    return True


def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def derived(kind, n, counts):
    if kind == "kp":
        sq = n * (n + 1) * (2 * n + 1) // 6
        D = {2 * (sq - S): c for S, c in counts.items()}
        t1, t2, t3 = (n - 1) * n * (n + 1) // 6, n * (n - 1) // 2, n * (n + 1) // 2
        return {"A073114": sum(c for S, c in counts.items() if is_prime(S)),
                "A180111": sum(c for d, c in D.items() if d <= t1), "A180112": sum(c for d, c in D.items() if d < t1),
                "A179263": sum(c for d, c in D.items() if d <= t2), "A179265": sum(c for d, c in D.items() if d < t2),
                "A179264": sum(c for d, c in D.items() if d <= t3), "A179266": sum(c for d, c in D.items() if d < t3)}
    if kind == "gcd":
        return {"A073343": counts.get(sigma(n), 0)}
    return {"A189043": sum(1 for c in counts.values() if c)}


def all_terms():
    rows = {}
    for f in sorted(os.listdir(DIST)) if os.path.isdir(DIST) else []:
        if not f.endswith(".json"): continue
        kind, n = f[:-5].rsplit("_", 1)
        for a, v in derived(kind, int(n), load(kind, int(n))).items():
            rows.setdefault(a, {})[int(n)] = v
    return rows


def known():
    return json.load(open(os.path.join(ROOT, "data", "known_terms.json")))["sequences"]


# ---------- command line ----------
def pieces(kind, n, k):
    smin, smax, T, N, M = plan(kind, n)
    step = -(-M // k)
    return [(qi, q, g, lo, min(M, lo + step)) for qi, (q, g) in enumerate(PRIMES) for lo in range(0, M, step)]


def main(argv):
    cmd = argv[0] if argv else "help"
    if cmd == "compute":
        kind, n = argv[1], int(argv[2]); threads = int(argv[3]) if len(argv) > 3 else 2
        t0 = time.time()
        counts = combine([engine(kind, n, q, g, threads=threads) for q, g in PRIMES])
        return 1 if save(kind, n, counts, time.time() - t0) else 0
    if cmd == "chunk":
        budget, t0 = float(argv[1]), time.time()
        os.makedirs(CHUNKS, exist_ok=True)
        for spec in argv[2:]:
            p = spec.split(":"); kind, n, k = p[0], int(p[1]), int(p[2]) if len(p) > 2 else 8
            for qi, q, g, lo, hi in pieces(kind, n, k):
                path = os.path.join(CHUNKS, f"{kind}_{n}_{qi}_{lo}_{hi}.txt")
                if os.path.exists(path): continue
                if time.time() - t0 > budget:
                    print("time budget used; run the same command again to continue"); return 0
                t1 = time.time()
                vals = engine(kind, n, q, g, mode=f"eval:{lo}:{hi}")
                with open(path + ".tmp", "w") as fh:
                    fh.write("".join(f"{m} {v}\n" for m, v in sorted(vals.items())))
                os.replace(path + ".tmp", path)
                print(f"{kind} n={n} prime {qi} points {lo}-{hi}: {time.time() - t1:.0f} s", flush=True)
        print("all pieces done"); return 0
    if cmd == "finalize":
        for spec in argv[1:]:
            p = spec.split(":"); kind, n, k = p[0], int(p[1]), int(p[2]) if len(p) > 2 else 8
            lists = []
            for qi, (q, g) in enumerate(PRIMES):
                files = [os.path.join(CHUNKS, f"{kind}_{n}_{i}_{lo}_{hi}.txt") for i, qq, gg, lo, hi in pieces(kind, n, k) if i == qi]
                missing = [f for f in files if not os.path.exists(f)]
                if missing: print(f"{kind} n={n}: {len(missing)} pieces missing"); break
                allpts = os.path.join(CHUNKS, f"{kind}_{n}_{qi}_all.txt")
                with open(allpts, "w") as fh:
                    for f in files: fh.write(open(f).read())
                lists.append(engine(kind, n, q, g, mode=f"combine:{allpts}"))
            else:
                save(kind, n, combine(lists), 0.0)
        return 0
    if cmd == "report":
        rows, kn, bad = all_terms(), known(), 0
        for a in sorted(rows):
            off, ks = kn[a]["offset"], kn[a]["terms"]
            last = off + len(ks) - 1
            same = sum(1 for n, v in rows[a].items() if off <= n <= last and ks[n - off] == v)
            diff = [n for n, v in rows[a].items() if off <= n <= last and ks[n - off] != v]
            bad += len(diff)
            new = sorted(n for n in rows[a] if n > last)
            print(f"{a}: {same} known terms reproduced" + (f", DIFFERENT at n = {diff}" if diff else "")
                  + (f"; new terms for n = {new[0]}..{new[-1]}" if new else ""))
            for n in new: print(f"    a({n}) = {rows[a][n]}")
        return 1 if bad else 0
    if cmd == "bfiles":
        rows, kn = all_terms(), known()
        out = os.path.join(ROOT, "results", "bfiles"); os.makedirs(out, exist_ok=True)
        for a in sorted(rows):
            off, ks = kn[a]["offset"], kn[a]["terms"]
            last = off + len(ks) - 1
            new = sorted(n for n in rows[a] if n > last)
            assert new == list(range(last + 1, last + 1 + len(new))), (a, new)
            terms = ks + [rows[a][n] for n in new]
            with open(os.path.join(out, f"b{a[1:]}.txt"), "w") as fh:
                fh.writelines(f"{off + i} {t}\n" for i, t in enumerate(terms))
            print(f"b{a[1:]}.txt: n = {off}..{off + len(terms) - 1}")
        return 0
    if cmd == "newterms":
        rows, kn = all_terms(), known()
        order = ["A073114", "A180111", "A180112", "A179263", "A179265", "A179264", "A179266", "A073343", "A189043"]
        lines = ["# New terms", "",
                 "Terms beyond those in the OEIS (export of 2026-10-29), computed and checked as described in README.md.", "",
                 "| Sequence | Counts | Known up to | New terms |", "| --- | --- | --- | --- |"]
        detail = []
        for a in order:
            if a not in rows: continue
            off, ks = kn[a]["offset"], kn[a]["terms"]; last = off + len(ks) - 1
            new = sorted(n for n in rows[a] if n > last)
            if not new: continue
            lines.append(f"| [{a}](https://oeis.org/{a}) | {kn[a]['name']} | n = {last} | n = {new[0]}..{new[-1]} ({len(new)}) |")
            detail += ["", f"## {a}", "", "| n | a(n) |", "| --- | --- |"] + [f"| {n} | {rows[a][n]} |" for n in new]
        print("\n".join(lines + detail)); return 0
    print(__doc__); return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
