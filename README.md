# permsum-oeis

Exact distributions of permutation sums S(p) = Σ f(k, p(k)) over all n! permutations p of {1, …, n}, and the new
OEIS terms that follow from them.

| Statistic | Sequences extended |
| --- | --- |
| Σ k·p(k) (equivalently Spearman's Σ (k − p(k))²) | [A073114](https://oeis.org/A073114), [A180111](https://oeis.org/A180111), [A180112](https://oeis.org/A180112), [A179263](https://oeis.org/A179263), [A179265](https://oeis.org/A179265), [A179264](https://oeis.org/A179264), [A179266](https://oeis.org/A179266) |
| Σ gcd(k, p(k)) | [A073343](https://oeis.org/A073343) |
| Σ k²·p(k) | [A189043](https://oeis.org/A189043) |

The new terms are listed in [results/new_terms.md](results/new_terms.md), the complete term lists in OEIS b-file
format in [results/bfiles](results/bfiles), and the full exact distributions (every value of S with its count) in
[results/distributions](results/distributions).

## Results

| Sequence | Counts | Known up to | New terms |
| --- | --- | --- | --- |
| [A073114](https://oeis.org/A073114) | Number of permutations p of [n] such that 1*p(1) + 2*p(2) + 3*p(3) + ... + n*p(n) is prime. | n = 18 | n = 19..26 (8) |
| [A180111](https://oeis.org/A180111) | Number of permutations of 1..n with the Sum_{i=1..n} of (i-p(i))^2 <= (n-1)*n*(n+1)/6. | n = 22 | n = 23..26 (4) |
| [A180112](https://oeis.org/A180112) | Number of permutations of 1..n with the Sum_{i=1..n} of (i-p(i))^2 < (n-1)*n*(n+1)/6. | n = 22 | n = 23..26 (4) |
| [A179263](https://oeis.org/A179263) | Number of permutations of 1..n with sum (i-p(i))^2 <= n*(n-1)/2. | n = 22 | n = 23..26 (4) |
| [A179265](https://oeis.org/A179265) | Number of permutations of 1..n with sum (i-p(i))^2 < n*(n-1)/2. | n = 22 | n = 23..26 (4) |
| [A179264](https://oeis.org/A179264) | Number of permutations p of 1..n with sum (i-p(i))^2 <= (n+1)*n/2. | n = 25 | n = 26..26 (1) |
| [A179266](https://oeis.org/A179266) | Number of permutations p of 1..n with sum (i-p(i))^2 < (n+1)*n/2. | n = 25 | n = 26..26 (1) |
| [A073343](https://oeis.org/A073343) | Number of permutations p of (1,2,3,...,n) such that Sum_{k=1..n} gcd(k,p(k)) = sigma(n). | n = 20 | n = 21..27 (7) |
| [A189043](https://oeis.org/A189043) | For all permutations of [1..n]: number of distinct values taken by sum(k=1..n, k^2 * pi(k) ). | n = 14 | n = 15..22 (8) |

A189043: for 8 <= n <= 22 every value between the smallest and largest possible sum is attained except S_max - d and S_min + d for d in {1, 2, 4, 6, 8}, so a(n) = (n-1)n(n+1)^2/6 - 9 there. We conjecture this holds for all n >= 8.

## How to cite

If you use this code or data, please cite it (GitHub's "Cite this repository" button uses
[CITATION.cff](CITATION.cff)):

> Victor Nidal Bassem Nagib, *permsum-oeis: exact distributions of permutation sums and new OEIS terms*,
> 2026, https://github.com/victornagib2/permsum-oeis

The code is MIT licensed and the computed data is CC BY 4.0; both require credit when reused (see
[LICENSE](LICENSE) and [DATA-LICENSE.md](DATA-LICENSE.md)).

## Method

The generating polynomial of the statistic is a permanent:

    P_n(x) = Σ_p x^S(p) = per[ x^f(i,j) ]_{i,j=1..n},

so the count of permutations with S(p) = s is the coefficient of x^s. `src/permpoly.c` evaluates P_n at the
N-th roots of unity modulo a prime q = c·2^20 + 1 < 2^62 (N a power of two larger than the range of S), using
Glynn's formula

    per(M) = 2^(1−n) Σ_{δ ∈ {±1}^n, δ_1 = 1} (Π_k δ_k) Π_i Σ_j δ_j M_ij

with the sign vectors in Gray-code order (one column update and n multiplications per step), Montgomery
arithmetic, eight evaluation points per inner loop and OpenMP over points. An inverse transform recovers the
counts modulo q; two primes and the Chinese remainder theorem give the exact counts (valid while n! < q1·q2,
i.e. n ≤ 33). When the distribution is symmetric (Σ k·p(k) and Σ k²·p(k), via p → n + 1 − p) only half of the
points are evaluated. The work is about N · 2^(n−1) · n modular multiplications per prime.

`permsum.py` drives the engine, splits long runs into resumable pieces, checks every distribution and derives
the sequence terms.

## How the results were checked

1. **Brute force** for n ≤ 8, all three statistics (`tests/test_small.py`).
2. **Every known OEIS term** of the nine sequences is reproduced (`python3 permsum.py report`); these were
   computed between 2002 and 2025 by several contributors with independent programs.
3. **Spearman's rho tables**: for n ≤ 22 the full distribution of Σ k·p(k) is identical, count for count
   (5,331 values), to the exact null distribution of Spearman's rho of
   [van de Wiel and Di Bucchianico (1998)](https://research.tue.nl/en/publications/fast-computation-of-the-exact-null-distribution-of-spearmans-rho-/)
   as packaged in the R package [pspearman](https://cran.r-project.org/package=pspearman) by P. Savicky
   (`tests/check_spearman.py`).
4. **Independent recent terms**: for n = 23, 24 and 25 the derived values of A179264 and A179266 equal those
   added to the OEIS by Alois P. Heinz in December 2025.
5. **Exact moments**: for every n, the counts sum to n! and have exactly the mean and variance given by
   Hoeffding's (1951) formula for Σ a(i, p(i)) over a uniformly random permutation.

The terms of A073114 for n = 19..22 could in principle have been derived from the 1998 Spearman tables, but they
were not in the OEIS.

## Reproducing

    gcc -O3 -march=native -fopenmp -o permpoly src/permpoly.c
    python3 tests/test_small.py
    python3 permsum.py compute kp 20          # one exact distribution (about 20 s on 2 cores)
    python3 permsum.py chunk 6000 kp:27:8     # long runs in resumable pieces; repeat until done
    python3 permsum.py finalize kp:27:8
    python3 permsum.py report                 # all terms, compared with data/known_terms.json
    python3 permsum.py bfiles

Measured times on 2 CPU cores (both primes): Σ k·p(k) n = 24: 10 min, n = 25: 21 min, n = 26: 44 min;
Σ gcd(k, p(k)) n = 26: 11 min. Each further n costs about 2.3 to 2.4 times as much.

## Credits

- Sequence authors and earlier extensions: Benoit Cloitre (A073114, A073343), R. H. Hardin (A179263–A179266,
  A180111, A180112), Joerg Arndt (A189043), with earlier terms by Hagen von Eitzen, Donovan Johnson,
  Mohammed Bouayoun, Rick L. Shepherd, Ryan Propper, Robert Gerbicz and Alois P. Heinz. Data in
  `data/known_terms.json` is from the OEIS (CC BY-SA 4.0).
- D. G. Glynn, The permanent of a square matrix, European J. Combin. 31 (2010) 1887–1891.
- M. A. van de Wiel and A. Di Bucchianico, Fast computation of the exact null distribution of Spearman's rho
  and Page's L statistic, Memorandum COSOR 98-17, TU Eindhoven (1998).
- W. Hoeffding, A combinatorial central limit theorem, Ann. Math. Statist. 22 (1951) 558–566.
- Background for A189043: R. Carpenter, C. Defant, N. Kravitz,
  [On the number of permutation-twisted dot products](https://arxiv.org/abs/2601.15276) (2026), and
  C. Pohoata, [Distinct permutation dot products](https://arxiv.org/abs/2601.12445) (2026).

The code and computations were produced with the help of an AI assistant (Claude, by Anthropic), working
under the direction of the author; every result was checked as described above.
