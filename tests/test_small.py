"""Brute-force test: for n <= 8 every distribution must equal a direct count over all n! permutations,
and the split (chunk + finalize) path must give the same result as a one-shot run."""
import itertools, os, shutil, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import permsum as P

for kind in ("kp", "gcd", "k2p"):
    for n in range(1, 9):
        E = P.matrix(kind, n)
        brute = {}
        for p in itertools.permutations(range(n)):
            S = sum(E[i][p[i]] for i in range(n))
            brute[S] = brute.get(S, 0) + 1
        got = P.combine([P.engine(kind, n, q, g, threads=1) for q, g in P.PRIMES])
        assert {S: c for S, c in got.items() if c} == brute, (kind, n)
        assert P.check(kind, n, got) is None, (kind, n)
    print(f"{kind}: n = 1..8 match brute force and pass the checks")

# split path
P.CHUNKS = os.path.join(ROOT, "results", "test_chunks")
shutil.rmtree(P.CHUNKS, ignore_errors=True)
for kind, n in (("kp", 11), ("gcd", 11), ("k2p", 9)):
    P.main(["chunk", "1000", f"{kind}:{n}:3"])
    lists = []
    for qi, (q, g) in enumerate(P.PRIMES):
        files = [os.path.join(P.CHUNKS, f"{kind}_{n}_{i}_{lo}_{hi}.txt") for i, qq, gg, lo, hi in P.pieces(kind, n, 3) if i == qi]
        allpts = os.path.join(P.CHUNKS, f"all_{qi}.txt")
        open(allpts, "w").write("".join(open(f).read() for f in files))
        lists.append(P.engine(kind, n, q, g, mode=f"combine:{allpts}"))
    assert P.combine(lists) == P.combine([P.engine(kind, n, q, g) for q, g in P.PRIMES]), (kind, n)
shutil.rmtree(P.CHUNKS, ignore_errors=True)
print("split runs equal one-shot runs")
