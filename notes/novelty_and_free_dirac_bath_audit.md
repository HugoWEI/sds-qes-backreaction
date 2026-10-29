# Novelty audit and exact free-fermion bath benchmark — 29 October 2026

## Research status

**No established novelty claim. No physical fine-grained SdS QES or Page-time shift has been derived.** The source-matched two-horizon EF geometry and the conditional free-fermion calculations are reproducible infrastructure. A new result would need a specified time-dependent gravity–bath boundary-value problem and replica saddles, with a correction not already present in the prior literature.

## Known prior art (check full papers before any novelty claim)

- Easson (2026), Phys Rev D **113** 084014, https://doi.org/10.1103/trlf-d6sk : anomaly-backreacted nonequilibrium UdS flux and **thermodynamic Page proxy**, explicitly defers full QES extremization.
- *Cosmological and black hole islands in multi-event horizon spacetimes*, Phys Rev D **107** 026009 (2023), https://doi.org/10.1103/PhysRevD.107.026009 : two-horizon SdS islands with opaque membranes.
- Goswami & Narayan, arXiv:2207.10724 (2022) and 2312.05904 (2024): frozen-background SdS islands and future-boundary radiation collection.
- *Universal structure of islands in evaporating black holes*, JHEP03(2025)068: general QES placement in an evaporating Unruh state.
- Đorđević & Radovanović, Phys Rev D **113** 105020 (2026): Polyakov-corrected *ordinary Schwarzschild* fine-grained Page curve.
- Balasubramanian, Kar, Ugajin, JHEP02(2021)072: CFT bath and gravitational replicas in backreacted de Sitter JT.
- Casini & Huerta, arXiv:0812.1773, 0903.5284: **exact** free Dirac multi-interval entropy is already known.
- Calabrese, Cardy, Tonni (2009) and nonselfdual compact-boson calculations: multi-interval entropy depends on four-point twist data, not just central charge. Headrick, Lawrence & Roberts (2012), arXiv:1209.2428, show important spin-structure/self-dual caveats.

## Exact free Dirac CFT result (established prior art)

For two disjoint spacelike intervals with endpoints 1<2<3<4 ordered along one slice and vacuum conformal coordinates U,V:

```
S(A union R) = c/6*(L(U)+L(V)+sum_i rho_i - sum_i log(a_i)) + const
L(z)=log|z12|+log|z34|+log|z23|+log|z14|-log|z13|-log|z24|
x_z=(z12*z34)/(z13*z24)
I(A:R)=-(c/6)*(log(1-x_U)+log(1-x_V)).
```

This formula assumes the fixed-spin-structure Gaussian massless Dirac vacuum. The mutual information is Weyl- and cutoff-independent. For static SdS (G=1, Lambda=3e-6, M=100, c=12), radii (300,350) and (450,500) on the same static slice give S_union=35.44159035500496 and MI=0.5002485802218247. Independent lattice tight-binding correlator with two intervals length 80/gap 40 gives MI=0.19590884984357437 versus continuum 0.19592888830070632. The prescribed first-order chart/geometry continuation has dS_union/d epsilon=1.7323804001387145e-6, verified by finite differences, **not** a physical QES shift.

## An explicit *auxiliary* nongravitating bath

Attach the static SdS exterior at r=b=500 with tortoise y=rstar(r)-rstar(b) ≤0 and prescribe a conformally flat bath at y≥0. Smooth conformal factor:

```
rho_b=0.5log f(b), A=0.5f'(b), D=0.5f(b) f''(b), ell=50
h1(y)=ell*(1-exp(-y/ell))
h2(y)=ell^2*(1-(1+y/ell)*exp(-y/ell))
rho_B(y)=rho_b + A*h1(y) + (D+A/ell)*h2(y).
```

This matches the metric through second derivative at the seam, tends to a constant Weyl factor at bath infinity, and permits a **chosen** conformal Dirac state U=-exp(-kb(t-y)), V=exp(kc(t+y)) with transparent field matching. **It is not Easson's model and its gravitational junction conditions have not been solved.**

Example for gravitational exterior interval (320,400) and bath interval (20,100) at t=0: S_union=34.877923218822886, MI=0.38294516896267033. Changing collar ell=50 to ell=10 changes S_union to 34.89219861689867, leaves MI unchanged; this is normal Weyl/cutoff behavior, not a claim about physically renormalized Page times.

The large spatial bath interval entropy slope tends to c*(kb+kc)/12=0.0028263436238952105, **not** outgoing-only chiral entropy slope c*kb/12=0.0020769186253272033 nor Easson's net thermodynamic entropy proxy.

## Two-null *negative control* with free-Dirac bath

Artificially fix one gravitational endpoint at r=400; vary another gravitational exterior endpoint (a,t_a) and use bath interval [20,100]. Extremize S_trial=pi*a²/G+S_Dirac([a,400] union [20,100]). Exact analytic gradients match independent finite differences. Formal stationary point a=399.9984084483509, t_a~0, radial gap=0.0015915516491, proper gap=0.0027294899974, far **below UV cutoff 1**. The leading gap cG/(6pi*400)=0.0015915494309 explains this UV-adjacent artifact. **Reject this saddle. It is NOT an island or Page transition and does not explore the black-hole interior.**

## Criteria before claiming novelty

1. Derive a globally consistent, time-dependent gravitating SdS plus nongravitating collector junction (including energy flux and boundary terms), beyond a toy smooth fixed bath.
2. Specify purification, free Dirac state/spin structure and a single radiation algebra on a spacelike slice, identical in the island/no-island branches.
3. Solve the matched backreacted geometry through the black-hole horizon; outgoing and ingoing EF outer expansions are nonuniform at opposite horizons.
4. Derive the replica channels and physical island endpoint boundary conditions, then solve both null QES equations and compare on-shell entropies.
5. Numerically check the nonperturbative QES roots, first-order displacement and Page crossing, including UV and near-Nariai control.
6. Compare the resulting term specifically with 2023 two-horizon SdS islands, 2025 universal Unruh QES and 2026 Polyakov Schwarzschild.

**Bottom line:** validated free-field entropy machinery, but no novel physics demonstrated yet.