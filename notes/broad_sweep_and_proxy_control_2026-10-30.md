# Broader literature sweep and proxy-control calculation — 30 October 2026

**Status:** Plausible narrow gap, not a verified novelty claim. There is NO physical QES, replica saddle, or fine-grained Page crossing yet. Our new results below are necessary *instantaneous* control diagnostics of **Easson's thermodynamic proxy** only.

## Broad prior-art comparison

- Easson, PRD **113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk — full anomaly flux and thermodynamic Page proxy, defers full microscopic QES. APS cited-by (checked Oct 30) lists Shi–Turiaci–Wu and Easson's RN-dS extension, neither a neutral UdS fine-grained QES phase diagram.
- Goswami & Narayan, JHEP10(2022)031 / JHEP05(2024)016 — small/frozen SdS islands, future-boundary collector; no anomaly-backreacted global phase diagram.
- Yadav & Joshi, PRD **107**, 026009 (2023) — two-horizon SdS islands using opaque membranes and frozen opposite patch.
- **Aalsma, Aguilar-Gutierrez & Sybesma**, JHEP01(2023)129, https://doi.org/10.1007/JHEP01(2023)129 — **controlled** information recovery is already demonstrated in a different 2D de Sitter–Rindler wall model with a *finite-duration* departure from equilibrium. Do not claim generic controlled de Sitter Page recovery is new.
- Aalsma & Sybesma, JHEP05(2021)291; Kames-King, Verheijden & Verlinde, JHEP03(2022)040 — collection-induced backreaction and obstructions in other de Sitter setups, not neutral SdS universal no-go.
- Balasubramanian, Kar & Ugajin, JHEP02(2021)072 — backreacted dS JT islands plus auxiliary universe.
- Yu & Ge, Chinese Phys C **49** 045107 (2025), https://doi.org/10.1088/1674-1137/ada34e — static near-horizon island criterion, with f''<0 sufficient in their construction; classical SdS already satisfies it.
- Matsuo, JHEP03(2025)068 — universal Unruh-vacuum interior QES placement.
- Ageev et al., JHEP05(2024)308 — de Sitter purity/complementarity issues; QES existence alone need not resolve entropy.
- Shi, Turiaci & Wu, JHEP08(2026)191 — smooth state evaporating dS black holes, different observer and state; no collector Page diagram.
- Easson, JHEP07(2026)212 — RN-dS discharge and anomaly flux, only conservative islands/QES estimates.

Search scope: targeted publisher pages, APS cited-by, preprint/INSPIRE-oriented queries and varied combinations of SdS, QES, Polyakov, island, backreaction, observer, bath, state, validity and 2025–Oct 2026. No exact **full neutral UdS, state-specific, anomaly-backreacted existence + dominance + control** theorem was located. Bibliographic absence is NOT a proof of novelty.

## New exact algebraic diagnostic: thermodynamic-proxy mass and adiabaticity

Use dimensionless x=sqrt(Lambda)*r_b in (0,1), y=sqrt(Lambda)*r_c=(sqrt(12-3*x^2)-x)/2, mu=G*M*sqrt(Lambda)=x*(1-x^2/3)/2, and q=N*G*Lambda. Then x^2+x*y+y^2=3.

From Easson's **proxy** crossing S_c(M_P)=S_b(M_0)+S_c(M_0), setting p=x_0*y_0:

    mu_P=p*sqrt(3-p)/6
    x_P=mu^{-1}(mu_P) on 0<x<1.

The temperature-scale adiabaticity diagnostics, under Easson's declared mass-loss convention, are

    eps_b=|dot(kappa_b)|/kappa_b^2=q*A_b(x)/(48*pi)
    A_b=(1+x^2)/(x^2*(1-x^2)) * [1-(kappa_c/kappa_b)^2]
    eps_c=|dot(kappa_c)|/kappa_c^2=q*A_c(x)/(48*pi)
    A_c=(1+y^2)/(y^2*(y^2-1)) * [(kappa_b/kappa_c)^2-1]
    kappa_b/sqrt(Lambda)=(1-x^2)/(2*x)
    kappa_c/sqrt(Lambda)=(y^2-1)/(2*y).

SymPy proves lim_(x->1) A_b=A_c=8/3 (apparent Nariai singularity cancels in these two specific ratios); at x->0, x^2*A_b->1 and x^2*A_c->1/2. The **near-Nariai coordinate and boundary-layer perturbation can still fail**: small adiabaticity is not sufficient for global regularity.

Scan along x∈[x_P,x_0] and also monitor gamma=q/x^2=N*G/r_b^2. At q=10^-6:

| x0 | mu_P/mu0 | x_P | max(eps_b,eps_c) | max gamma |
|---:|---:|---:|---:|---:|
| .05 | .97195138 | .04859532 | 2.8130e-6 | 4.2346e-4 |
| .20 | .89704423 | .17892614 | 2.1292e-7 | 3.1236e-5 |
| .50 | .78467529 | .37758733 | 5.3281e-8 | 7.0140e-6 |
| .80 | .72041827 | .49342853 | 3.4526e-8 | 4.1073e-6 |
| .99 | .70714202 | .51757391 | 3.2156e-8 | 3.7330e-6 |

These weak instantaneous bounds show ONLY that a thermodynamic-proxy crossover can occur while both instantaneous horizon adiabaticity measures and cG/r_b^2 are small. They do **not** show a physical island, controlled collector, Page crossing, horizon-boundary-layer validity, or a new theorem.

## Reproduction and next physical task

Run `python src/proxy_viability.py` and `python -m pytest -q`. Local complete snapshot: **60 passed** (prior 54 + six new tests, including symbolic Nariai and independent Vieta/finite-difference checks). Remote GitHub Actions result unconfirmed.

Next: specify a finite-duration free-Dirac observer–bath coupling and purification, derive the transparent interface's energy conservation and gravitational boundary/junction conditions, and solve the behind-horizon two-null QES using a regular chart with the SAME radiation algebra for both competing saddles. Compare to Matsuo 2025, Goswami–Narayan, and Aalsma et al. 2023 before claiming either a new allowed domain or obstruction. The 2D-to-4D flux normalization and proper counterterms remain audit items.
