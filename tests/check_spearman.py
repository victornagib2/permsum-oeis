"""Compare results/distributions/kp_n.json with the independent exact distribution of Spearman's rho for n <= 22
in the R package pspearman (P. Savicky; values also checked against M. A. van de Wiel's program).
    git clone https://github.com/cran/pspearman
    python3 tests/check_spearman.py pspearman/R/spearman.data.R
The package lists, for n = 2..22, the counts of v = sum (i - p(i))^2 / 2 = n(n+1)(2n+1)/6 - S for the first half
of the (symmetric) range."""
import os, re, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import permsum as P

txt = open(sys.argv[1]).read()
vecs = re.findall(r"c\(([^)]*)\)", txt[txt.index("list(") + 5:])
total = 0
for idx, vec in enumerate(vecs):
    n = idx + 1
    if n < 2: continue
    theirs = [int(x) for x in vec.replace("\n", "").split(",") if x.strip()]
    c = P.load("kp", n)
    sq = n * (n + 1) * (2 * n + 1) // 6
    ours = [c.get(sq - v, 0) for v in range(len(theirs))]
    assert ours == theirs, f"n={n} differs"
    total += len(theirs)
print(f"identical for n = 2..{len(vecs)}: {total} counts compared")
