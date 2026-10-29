# Research update 4: outgoing radiation algebra, UV-saddle rejection, and conditional peeling response

**Date:** 29 October 2026. **Status:** 8 new local tests passed; **37 locally tested cases counting the prior 29**. GitHub Actions remote result not independently verified. No physical SdS island, QES displacement or fine-grained Page transition has been computed.

## Literature and scope

- Easson, Phys. Rev. D **113**, 084014 (2026), [DOI](https://doi.org/10.1103/trlf-d6sk): two-horizon anomaly flux and thermodynamic Page **proxy**. Full QES extremization deferred.
- Goswami & Narayan, JHEP **10** (2022) 031: frozen SdS islands, small black hole.
- Goswami & Narayan, JHEP **05** (2024) 016: SdS future-boundary radiation region.
- Aalsma, Parikh & van der Schaar, JHEP **11** (2019) 136: two-dimensional Unruh-de Sitter stress/backreaction.
- Calabrese & Cardy, J. Phys. A **42** (2009) 504005: universal 2D CFT interval entropy.
- Đorđević & Radovanović, Phys. Rev. D **113** (2026) 105020: evaporating Schwarzschild islands (no cosmological horizon).

These ingredients are known; their combination here is a reproducibility exercise, **not an originality claim**.

## 1. Precisely chosen radiation **subalgebra**, not yet a bath

Place an ideal *mode reader* at fixed radius `b`, `rb<b<rc`. Let `R_[u0,u1]` be the outgoing chiral subalgebra of conformal modes passing the reader between retarded times `u0,u1` (not a region on a single spacelike slice). Assume it has the Gaussian chiral vacuum `U=-exp(-kb*u)`. Its regulated chiral interval entropy, for fixed observer proper cutoff `a_uv`, is

```
S_out(Delta_u) = (c/6) log[(2 sqrt(f(b)) / (kb*a_uv)) * sinh(kb*Delta_u/2)] + C
dS_out/du = c*kb/12 * coth(kb*Delta_u/2)
lim_{Delta_u -> large} dS_out/du = c*kb/12.
```

This formula follows from the chiral conformal transformation law for interval entropy and is **not** the entropy of all radiation emitted to a nongravitating collector. The constant C, the algebraic UV regulator, and any memory/storage/bath coupling have not been fixed. In particular, **the formula alone does not justify applying the island replica formula**.

Compare Easson's steady entropy-production proxy:

```
J = c*(kb**2-kc**2)/(48*pi)
s_proxy = J*(2*pi/kc-2*pi/kb)
        = c*(kb-kc)**2*(kb+kc)/(24*kb*kc)
```

At `G=1, Lambda=3e-6, M=100, c=12, b=500`, the late chiral mode-interval entropy growth rate is `0.0020769186253272033` versus Easson's net thermodynamic proxy `0.0015999713235596304`. These are **different observables**, not conflicting calculations. Do not extrapolate the stationary chiral expression to an evaporation-time-scale Page crossing.

## 2. A two-null exterior trial extremum, explicitly rejected

Stress-test the **conditional** generalized entropy `S_trial(a,t_a)=pi*a²/G+S_CFT((a,t_a),(b,t_b))`. This is an ordinary single-interval CFT entropy between two **exterior** spacelike endpoints, not a replica entropy of `R union I`; b is fixed. The formula is not justified for a point behind the black-hole horizon.

Let `D=rstar(b)-rstar(a)>0`, `dt=t_a-t_b`, so `Delta_u=dt+D`, `Delta_v=dt-D`, and `|dt|<D`. For `f=f0(a)`, the two exact extremality equations are

```
0 = (c/12) [kb*coth(kb*(dt+D)/2) - kc*coth(kc*(D-dt)/2)]
0 = (2*pi*a/G) + (c/(12*f))[
        f'(a)-kb*coth(kb*(dt+D)/2)-kc*coth(kc*(D-dt)/2)].
```

At each `a<b`, the time-gradient is strictly decreasing with dt, and its signs at the two spacelike endpoints are opposite, so it has **exactly one** time-stationary solution. Substitute that into the radial equation.

For the reference parameters, the only exterior radial root found on the logarithmic/adaptive scan is

```
a = 499.99872675744456,   dt ≈ 0,
b-a = 0.00127324255544,
near-endpoint estimate c G/(6 pi b) = 0.00127323954474.
```

The leading near-endpoint formula follows from `S_CFT~(c/3)log(ell/a_uv)`, `ell~(b-a)/sqrt(f(b))`, giving `d_a S_trial~2*pi*b/G-c/[3(b-a)]`.

The proper separation `ell~0.002152172726` is **below the explicit unit UV cutoff**. The local Hessian in (a,dt) has negative eigenvalues near `-2.456e6` and `-3.009e5`. Negative Hessian alone does not disqualify a Lorentzian QES, but the below-cutoff distance and absence of a legitimate radiation replica problem **do**. No macroscopic black-hole-horizon island arises from this particular exterior interval ansatz for these parameters. This is a parameter-specific negative control, not a theorem excluding interior islands.

## 3. First-order adiabatic **peeling ansatz** for a chiral entropy

A time-dependent quantum state is not specified by the instantaneous stress alone. Define an *additional* approximation:

```
M(u) = M0 - epsilon J0*u,
kappa(u) = kb0 + epsilon*alpha*u,
alpha = -J0 * (d kb/dM)_M0,
d kb/dM = -G*(1/rb²+Lambda)/(1-Lambda*rb²).
```

Assume the outgoing vacuum coordinate obeys the **peeling prescription** `U'(u)=exp(-int_0^u kappa(s)ds)`. This is an ansatz, not derived from the complete EF geometry and not a global state solution. With fixed coordinate UV cutoffs (so no detector redshift contribution), the regulated chiral interval entropy for `[0,T]` is

```
K(u)=kb0*u+epsilon*alpha*u²/2
I(epsilon)=int_0^T exp[-K(u)]du
S(epsilon)=(c/6)[log(I(epsilon)/a_u)+K(T)/2].
```

Differentiate analytically at `epsilon=0`:

```
I0 = (1-exp(-kb0*T))/kb0
I2 = 2/kb0³ * [1-exp(-kb0*T)*(1+kb0*T+(kb0*T)²/2)]
dS/depsilon|0 = (c*alpha/6) * [T²/4-I2/(2*I0)].
```

At `T=20000` in the reference parameters, `alpha=8.887557839243327e-12`, `J0*T/M0=5.9714188971398796e-5`, `alpha*T/kb0=8.558407374138848e-5`. The **coordinate-cutoff** entropy is `S0=53.8921123341493`, with first-order response `S1=0.0017733908443548453`. Numerical quadrature for epsilon=1 gives `53.89388572501488`; an independently integrated ODE for `(U',int U')` gives `53.893885725014876`; the linear prediction differs from the finite-epsilon change by `2.12e-11`. Central finite differences confirm the derivative.

The baseline entropy normalization here differs from the earlier fixed-*proper*-cutoff formula by a known stationary regulator term. The first-order response excludes metric/detector redshift, state-dependent counterterms, higher-order in mass drift, and all replica-wormhole and island effects. Thus **it is not a physical backreaction Page-time shift**.

## 4. Why the Page transition remains underdetermined

To evaluate `S_matter(R union I)` one needs the density matrix and a four-point twist correlator of the relevant interval endpoints, with appropriate boundary channel/cosmic-horizon continuation:

```
Tr rho_(R union I)^n ~ <Twist_n(x1) AntiTwist_n(x2) Twist_n(x3) AntiTwist_n(x4)>.
```

The correlation function generally includes a theory-dependent function of conformal cross ratios, not determined solely by the central charge or Polyakov one-point stress. The gravitational replica saddles, initial state/purifier, boundary conditions, UV counterterms, and a common no-island radiation algebra also need specification. Easson's flux law and our coordinate-patch metric solver alone are insufficient to fix a unique fine-grained Page time.

**Next defensible research extension:** choose free massless Dirac (or other explicit unitary CFT), a genuine transparent coupling to a nongravitating collector with boundary conditions, and a fixed initial state/purification. Then compute Gaussian multi-interval entropy or explicit twist correlators, solve both null QES equations, and compare on-shell saddles; perform horizon boundary-layer matching before importing the first-order metric into a horizon-scale QES. Until then no novelty claim.

Run `python src/radiation_algebra.py`, `python src/static_patch_candidate.py`, `python src/peeling_response.py`, and `python -m pytest -q`.
