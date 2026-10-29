# Model audit and first-order derivations (2026-10-29)

**Do not equate a thermodynamic proxy with the fine-grained entropy of a specified radiation subsystem.**

## Source equations and assumptions

Easson [arXiv:2511.11873v3](https://arxiv.org/html/2511.11873), PRD 113 084014 (2026), Sections II, III, V, VIII, IX, Appendices A/B:

- Spherical gravity: `S2=(1/(4G4)) integral sqrt(-g)[XR+(grad X)^2/(2X)+2(1-Lambda X)]`, `X=r^2`.
- SdS `f0=1-2G4 M/r-Lambda*r^2/3`; `T_h=kappa_h/(2pi)`.
- Polyakov `S_P=-(N/(96pi)) integral sqrt(-g) R box^{-1} R`. The state-dependent conformal stress is `T_uu=-N/(12pi)[(partial_u rho)^2-partial_u^2 rho-t_u(u)]`; likewise for `vv`; `T_uv=-N/(12pi)partial_u partial_v rho`.
- Effective outgoing flux `J=N*(kappa_b^2-kappa_c^2)/(48pi)`, slow mass loss `dot M=-J`.
- The published QES prescription is `partial_+ S_gen=partial_- S_gen=0`, `S_gen=X/(4G2)+S_matter`, **but not an evaluated QES**.
- EF gauge: `ds^2=-f(v,r) dv^2+2 dv dr`, with `X=r^2`, auxiliary field `2 partial_v psi_prime+(f psi_prime)prime=-f_double_prime`.

## Nonnegotiable conceptual checks

1. **Quantum state:** UdS regularity of the stress on appropriate future horizons and constant leading state functions do not alone determine the global density matrix, all interval/twist entropies, the purification and cross-ratio channels. A conformal vacuum prescription must be stated explicitly.
2. **Radiation algebra:** stationary in-patch entropy is assigned `dot S_out=0` in Section VIII, while Section IX defines an *increasing thermodynamic entropy-production proxy*. They are not the same observable. Fix `R` before comparing island and no-island saddles.
3. **Geometry:** no explicit unique first-order flux-carrying `f1(v,r)` and `X1(v,r)` is tabulated for a physical QES calculation. Static first integrals refer to a zero-flux sector. The time-dependent metric response must be solved with boundary/gauge conditions.
4. **Limits:** the paper's slow-time criterion `|partial_v f| << |kappa f|` is non-uniform as `f -> 0`. A QES near a *moving* horizon requires a regular matched/EF expansion.
5. **Effective theory:** an anomaly model for minimally coupled conformal 2D matter is not identical to fully reduced 4D matter; the dimensional-reduction anomaly is omitted (Easson Appendix A.5).
6. **Dominance:** `kappa_b>kappa_c` determines flux direction, not which replica saddle dominates a chosen entropy.

## General QES first-order response

Assume a consistent perturbative semiclassical state and solution have been supplied:

```
g_e = g0 + e*g1,   X_e = X0 + e*X1,
S(q,e) = S0(q) + e*S1(q) + O(e**2),
q = (u_I,v_I).
S1(q) = X1(q)/(4 G2) + S_matter,1(R union I; q).
```

At an isolated base QES `grad S0(q0)=0`, define `H_ab=partial_a partial_b S0(q0)`. Provided `det H != 0`, differentiating both null equations yields

```
q1^a = -(H^{-1})^{ab} partial_b S1(q0).
```

A singular Hessian marks a saddle merger; an ordinary linear expansion need not be uniform there. At first order the on-shell island entropy shifts by `S1(q0)`, independent of `q1` (envelope theorem).

For two entropies defining the **same observable** set `D(t,e)=S_no(t,e)-S_is,on_shell(t,e)`. For a simple crossing `D(t0,0)=0`, `partial_t D(t0,0)!=0`,

```
t_Page,1 = -partial_e D(t0,0)/partial_t D(t0,0).
```

This is *not* a physically computed Page time without the matter entropies and a flux-compatible metric.

## Horizon and proxy perturbations (conditional)

If `f=f0+e*h` at fixed Casimir mass, the shift of any simple horizon is

```
r_h,1 = -h(r_h,0) / f0_prime(r_h,0).
kappa_h,1 = sigma_h*(h_prime(r_h)+r_h,1*f0_double_prime(r_h))/2,
sigma_b=+1, sigma_c=-1.
```

With `X=r^2+e*x`, the fixed-mass area contribution shifts by `S_h,1=pi*(2*r_h*r_h,1+x(r_h))/G4`. Neither `h` nor `x` is known until the semiclassical equations are solved.

If we additionally posit a perturbed version of **Easson's thermodynamic Page proxy**,

```
D(M,e)=integral_M^M0 (1/Tc(m,e)-1/Tb(m,e)) dm - Sb(M,e),
M_Page,1=-D1(M_Page,0)/D0_prime(M_Page,0),
D1(M)=integral_M^M0 (-Tc1/Tc0**2 + Tb1/Tb0**2) dm-Sb1(M).
```

For `J=J0+e*J1`, fixed initial mass:

```
t_Page,1 = -M_Page,1/J0(M_Page,0) -
           integral_M_Page,0^M0 J1(m)/J0(m)**2 dm.
```

One must recheck the gravitational first laws when using perturbed temperatures and entropies; do not assume this proxy is a physical fine-grained QES calculation.

**Illustrative response only:** `J_e(M)=J0(M)[1+e M/M0]` with unchanged proxy entropy implies `M_Page,1=0` and `t_Page,1= -integral (M/M0)/J0 dM`. For `G=1, Lambda=3e-6, N=12, M0=100`, `t_Page,0=43899849.00766`, `dt/de=-40511130.26773`; at `e=0.01` quadrature gives `43498452.01144` versus independent time-domain ODE `43498452.01891`.

The one-coordinate manufactured test in `tests/test_baseline.py` uses `S_toy(r)=pi*r^2/G+alpha*log(r_b-r)`, giving `r_toy=(r_b+sqrt(r_b**2-2alpha G/pi))/2`. **This is not a regular global UdS entropy or a two-null QES**, just a derivative unit test.

## Prior art / no-novelty claim

- Easson, PRD 113 084014 (2026): mass evolution and thermodynamic proxy; full QES extremization explicitly deferred.
- Goswami & Narayan, JHEP10(2022)031: frozen SdS, small-black-hole regime.
- Goswami & Narayan, JHEP05(2024)016: future-boundary radiation, distinct state/region prescription.
- Đorđević & Radovanović, PRD 113 105020 (2026): evaporating Schwarzschild *without* a cosmological horizon, Polyakov QES.
- Baek & Choi, JHEP05(2023)098: backreacted islands in de Sitter JT-style theory.

**Literature search is preliminary; novelty remains unverified.**