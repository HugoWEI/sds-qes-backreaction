# Dual EF patches and conditional CFT interval entropy — 29 October 2026

**Status:** local first-order overlap calculation (not a global horizon-crossing solution, QES saddle or fine-grained Page time). Reproduce with `python -m pytest -q`: 29 local tests.

## Source and literature audit

- D. A. Easson, *Fate of Schwarzschild–de Sitter black holes: Nonequilibrium evaporation*, PRD **113**, 084014 (2026), DOI [10.1103/trlf-d6sk](https://doi.org/10.1103/trlf-d6sk): anomaly-generated steady flux and thermodynamic Page proxy, but not explicit full time-dependent QES.
- L. Aalsma, M. Parikh, J. P. van der Schaar, *Back(reaction) to the future in the Unruh-de Sitter state*, JHEP **11** (2019) 136, DOI [10.1007/JHEP11(2019)136](https://doi.org/10.1007/JHEP11(2019)136): existing 2D de Sitter state/backreaction work; the two-horizon regularity idea is **not** claimed new.
- K. Goswami and K. Narayan, JHEP **10** (2022) 031, DOI [10.1007/JHEP10(2022)031](https://doi.org/10.1007/JHEP10(2022)031): frozen SdS background QES.
- Same authors, JHEP **05** (2024) 016, DOI [10.1007/JHEP05(2024)016](https://doi.org/10.1007/JHEP05(2024)016): future-boundary SdS islands.

## Explicit chosen state and basic stress

Instantaneous `f(r)=1-2GM/r-Lambda*r²/3`, `rb<rc`, `kb=f'(rb)/2`, `kc=-f'(rc)/2` and `drstar/dr=1/f`. Choose the *additional* conformal vacuum prescription `U_b=-exp(-kb*u)`, `V_c=-exp(-kc*v)` with `u=t-rstar`, `v=t+rstar`. With chiral stress functions `t_u=kb²/4`, `t_v=kc²/4`, the Polyakov state has `T_uu-T_vv=J=N(kb²-kc²)/(48pi)`. **The stress tensor alone does not determine a global density matrix or entropy; this chiral-vacuum assumption is extra input.**

## Both EF metric patches

Black-hole future horizon: `ds²=-exp(2sigma_plus)*F dv²+2exp(sigma_plus)dv dr`. Cosmological future horizon: `ds²=-exp(2sigma_minus)*F du²-2exp(sigma_minus)du dr`. Here `F=(grad r)²` in the areal-radius gauge and the outgoing stress components obey

```text
T_rr(-)=4*T_vv/f²
T_ur(-)=2*(T_vv+T_uv)/f
sigma_minus,1' = G*T_rr(-)/r
(r*F1_minus)'=-2*G*T_ur(-)
partial_u F1_minus=2*G*J/r
```

The outgoing Einstein-dilaton metric equations were independently derived with SymPy, and the radial quadratures and independent ODE solver agree. With `rstar(r0)=0` and the same mass/reference integration constant in both charts,

```text
F1_plus(v,r)=F1_minus(u=v-2*rstar(r),r).
```

This exact **first-order scalar matching identity** follows from `T_uu-T_vv=J`. This does not provide a physical global horizon-matched metric, or fix the homogeneous mass/time gauge parameters.

## Null-coordinate matching including the redshift factors

Write `u=v-2rstar(r)+epsilon*chi(v,r)+O(epsilon²)`. The null condition and metric cross term imply

```text
chi_v = sigma_plus,1(r) - sigma_minus,1(r)
chi_r = 2*sigma_minus,1(r)/f + 2*F1_plus(v,r)/f²
```

The integrability condition holds identically:

```text
d_r(sigma_plus,1-sigma_minus,1)
 = 4G*J/(r*f²)
 = d_v(2*sigma_minus,1/f+2*F1_plus/f²).
```

A solution satisfying `chi(v,r0)=0` is

```text
chi(v,r)=v*(sigma_plus,1(r)-sigma_minus,1(r))
       +integral_(r0)^r [2*sigma_minus,1(s)/f(s)
                         +2*F1_plus(0,s)/f(s)²] ds.
```

SymPy and numerical transformation checks validate the two EF metric forms to first order. Example `G=1, Lambda=3e-6, M=100, N=12, r0=400, r=500, v=0`:

```text
sigma_plus,1=+1.8264080066815465e-6
sigma_minus,1=-3.474260169511522e-7
F1_plus=F1_minus=-3.679751368641078e-7
```

The opposite chart develops coordinate divergences near the other future horizon, so the overlap expansion is non-uniform at the endpoints. No uniform near-horizon, fully dynamical solution is claimed.

## A **specified** finite CFT interval; not collected Hawking radiation

For N central-charge-one conformal fields in the chosen chiral vacuum, the standard two-dimensional interval formula for a spacelike interval `p=(u1,v1)`, `q=(u2,v2)`, physical endpoint cutoffs `a1,a2`, and `ds²=-e^(2rho)du dv` is

```text
S_interval = N/6 * log[
 (4*exp(rho(p)+rho(q))/(kb*kc*a1*a2))
 * abs(sinh(kb*(u1-u2)/2)*sinh(kc*(v1-v2)/2))].
```

For endpoints at **equal static time** `Delta u=-Delta rstar`, `Delta v=+Delta rstar`, this stationary background entropy is time-independent. It therefore must **not** be equated with Easson's growing thermodynamic Page proxy.

At fixed reference advanced-coordinate endpoints `(v_i,r_i)`, and **conditional on continuing the same chiral-vacuum prescription** in perturbed null coordinates, first-order conformal geometry yields

```text
rho1(v,r)=0.5*(F1_plus(v,r)/f+sigma_plus,1(r)+sigma_minus,1(r))
d_epsilon S_interval|0 = N/6*[
  (kb/2)*coth(kb*(u1,0-u2,0)/2)*(chi1-chi2) + rho1(p)+rho1(q)].
```

For background simultaneous-time endpoints at `r1=400,r2=500,t=0`, fixed cutoffs `a1=a2=1` and the example above: `S_interval,0=20.520454221120325`, `dS_interval/d_epsilon=-1.1063082186803846e-7`. Independent central finite differences reproduce the derivative at ~1e-14 absolute error. **Neither number is a QES entropy or Page-time shift.** The response depends on the endpoint/worldline and reference-gauge prescription; the quantities have not been made diffeomorphism-invariant detector observables.

## Blockers before an actual QES paper

1. State an explicit radiation algebra, purifier, endpoint worldlines and renormalization (for example a nongravitating collector with specified matching).
2. Match solutions through a controlled horizon boundary layer and specify boundary/gauge conditions; do not extrapolate the outer EF series through `f=0`.
3. Construct the fine-grained entropy of `R union I` and the no-island radiation entropy using *the same* radiation subsystem.
4. Extremize both null coordinates, compare on-shell saddles, and verify a genuine first-order Page transition against full numerical roots.
5. Compare the candidate result thoroughly with prior work. **Novelty not established.**

Code: [src/two_patch.py](../src/two_patch.py), [src/interval_entropy.py](../src/interval_entropy.py), [tests/test_two_patch.py](../tests/test_two_patch.py), [tests/test_interval_entropy.py](../tests/test_interval_entropy.py).
