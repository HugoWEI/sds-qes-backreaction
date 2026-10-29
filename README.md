# SdS QES backreaction — research notebook

**Started:** 29 October 2026. **Status:** literature/model audit, analytic thermodynamic-proxy reproduction, and conditional first-order perturbation theory. **A physical QES displacement or fine-grained Page-time correction has NOT yet been calculated.**

## Research question

Can we construct explicit quantum extremal surfaces and their saddle transitions for a *defined fine-grained radiation subsystem* in the anomaly-backreacted, nonequilibrium Schwarzschild–de Sitter (SdS) model of Damien A. Easson, [Phys. Rev. D 113, 084014 (2026)](https://doi.org/10.1103/trlf-d6sk), [arXiv:2511.11873v3](https://arxiv.org/html/2511.11873)?

This is distinct from (i) a coarse-grained thermodynamic Page-time **proxy**, (ii) frozen-background SdS islands, and (iii) evaporating ordinary Schwarzschild with Polyakov backreaction.

## First exact baseline: the **proxy**, not microscopic QES

For the classical instantaneous SdS metric
```
f(r) = 1 - 2*G*M/r - Lambda*r*r/3
S_b = pi*r_b*r_b/G
S_c = pi*r_c*r_c/G
J(M) = N*(kappa_b**2-kappa_c**2)/(48*pi)
```
Easson uses `dS_proxy/dt = J*(1/T_c-1/T_b)`, `dM/dt=-J`. The geometric first laws imply `dS_b/dM=1/T_b` and `dS_c/dM=-1/T_c`, so the integral reduces exactly to

```
S_proxy(M) = S_b(M)+S_c(M)-S_b(M0)-S_c(M0)
S_proxy(M_Page) = S_b(M_Page)
=> S_c(M_Page) = S_b(M0)+S_c(M0).
```

Thus `r_c(M_Page) = hypot(r_b(M0),r_c(M0))` and `M_Page = r_c*(1-Lambda*r_c**2/3)/(2*G)`, **within the thermodynamic proxy only**. No claim of new QES physics is made.

Illustration: `G=1, Lambda=3e-6, N=12, M0=100`; `r_b=209.14884844`, `r_c=878.88506625`, proxy `M_Page=83.0330759993`, proxy `t_Page=43899849.00766` in model units.

## Conditional first-order formulas

Let `q=(u_I,v_I)` and `S_gen(q,epsilon)=S0(q)+epsilon*S1(q)+O(epsilon**2)`. For a nondegenerate QES `grad S0(q0)=0`, with invertible `H=Hessian(S0)(q0)`,

```
q1 = - inverse(H) @ grad(S1)(q0).
```

At the **on-shell** island saddle the displacement term cancels at first order (envelope theorem). If `D(t,epsilon)=S_no_island(t,epsilon)-S_island_on_shell(t,epsilon)` has a simple zero at `t0`,

```
t_Page,1 = - partial_epsilon D(t0,0) / partial_t D(t0,0).
```

These are universal implicit-function formulas, **not physically evaluated SdS QES corrections**. We must derive the anomaly-corrected metric/dilaton, a well-defined UdS-state matter entropy and the radiation algebra first.

For a manufactured sensitivity experiment `J_epsilon(M)=J0(M)*(1+epsilon*M/M0)`, holding the entropy **proxy** fixed, `dM_Page/depsilon=0` and

```
dt_Page/depsilon = -integral(M_Page, M0) (M/M0)/J0(M) dM.
```

This sensitivity is independently verified against central differences and time-domain ODE integration.

## Reproduce

```bash
python -m pip install -r requirements.txt
python src/baseline.py
python -m pytest -q
```

## Physical blockers / next objectives

1. Fix the radiation subsystem and algebra: detector interval inside the static patch, collected flux beyond the cosmological horizon, or explicit nongravitating collector. These are not the same observable.
2. Specify a globally regular UdS-state conformal/twist-correlator prescription, including the two-horizon chiral sectors, interval endpoints and UV subtraction. Polyakov stress functions alone cannot fix all entropies.
3. Solve the first-order *time-dependent* Einstein–dilaton–Polyakov equations in EF gauge and check gauge fixing, horizon regularity and adiabatic uniformity near QES.
4. Extremize **both** null directions; check every candidate saddle and compare entropies in the **same** radiation problem.
5. Validate perturbation formulas against numerical roots of the nonlinear QES equations; analyze singular Hessians and near-Nariai/small-black-hole limits.

## Literature comparison (non-exhaustive)

- Easson (2026), [PRD 113 084014](https://doi.org/10.1103/trlf-d6sk): anomaly flux, mass loss, GSL, **thermodynamic Page proxy**; leaves full QES extremization open.
- Goswami & Narayan (2022), [JHEP10(2022)031](https://doi.org/10.1007/JHEP10(2022)031): SdS islands for small black holes on a **frozen** de Sitter background.
- Goswami & Narayan (2024), [JHEP05(2024)016](https://doi.org/10.1007/JHEP05(2024)016): future-boundary radiation region and SdS islands; a different entropy setup.
- Đorđević & Radovanović (2026), [PRD 113 105020](https://doi.org/10.1103/nc2z-hx6s): Polyakov-corrected evaporating **ordinary Schwarzschild** and fine-grained Page curve.
- Baek & Choi (2023), [JHEP05(2023)098](https://doi.org/10.1007/JHEP05(2023)098): backreacted islands in **JT-like two-dimensional de Sitter**, not full SdS.

**Novelty is not established.** See `notes/model_audit.md` for assumptions and derivations.
## Update: local first-order EF metric (29 October 2026)

- [Detailed derivation and physical limitations](notes/linearized_solution.md)
- [First-order anomaly-sourced geometry solver](src/linearized_geometry.py)
- [Symbolic/independent numerical tests](tests/test_linearized_geometry.py)

The exact metric equations in areal-radius EF gauge yield, at first order on classical SdS:
```
sigma1(r)=G*integral_(r0)^r T_rr(s)/s ds
F1(v,r)=[2G*integral_(r0)^r T_vr(s) ds +2G*J*(v-v0)]/r
```
The remaining dilaton equation follows from an independently proven Polyakov stress-conservation identity. At the representative test point (G=1, Lambda=3e-6, M=100, N=12, r0=400, r=500), sigma1=1.8264080066815465e-6 and F1_radial=-3.679751368641078e-7. The first-order tensor field equations agree with the 2D anomaly stress.

**Boundary of validity:** reference gauge/mass constants are not fixed by global boundary data, and the ingoing EF expansion is nonuniform at the cosmological future horizon. A two-patch matched geometry and defined radiation entropy remain outstanding; no physical QES position or fine-grained Page-time shift is claimed. See the derivation note for source normalization and prior-art checks.

Run `python src/linearized_geometry.py` and `python -m pytest -q` to reproduce all 16 local test cases.
