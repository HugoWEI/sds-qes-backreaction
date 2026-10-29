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

### Additional validation: conditional black-hole horizon root

A one-sided regular-horizon extrapolation predicts `delta r_b/eta=+2.58385855e-4` for the **chosen** r0=400 integration constant. Direct nonlinear roots of `f0+eta*F1` agree within `8.2e-6` relative error even at eta=1, with smaller errors for smaller eta. See [regression](tests/test_horizon_response.py). This root is not a QES and is not independent of the global mass/gauge matching convention.

**Latest local verification: 17 tests passed.** GitHub Actions has been configured; its remote completion state has not been separately confirmed.

## Update 3: dual horizon-adapted patches and conditional CFT entropy (29 October 2026)

**29 local tests passed**, including exact outgoing-EF metric identities, independent ODE/quadrature, null-coordinate integrability, same-event matching, stationary CFT interval entropy and its first-order finite-difference response.

- [Derivation and limitations](notes/two_patch_entropy_audit.md)
- [Two-patch solver](src/two_patch.py)
- [Conditional interval-entropy response](src/interval_entropy.py)
- [Two-patch tests](tests/test_two_patch.py)
- [Interval tests](tests/test_interval_entropy.py)

**Results:** Using the same Unruh–de Sitter chiral stress and homogeneous reference normalization, the two EF charts satisfy the first-order overlap identity `F1_+(v,r)=F1_-(v-2rstar(r),r)`. The null-coordinate correction obeys `chi_v=sigma_+-sigma_-`, `chi_r=2sigma_-/f+2F1_+/f^2`; its integrability follows exactly from `T_uu-T_vv=J`. A separately specified two-chiral vacuum yields a **stationary** equal-static-time finite-interval entropy, rather than Easson's growing thermodynamic Page proxy; its first-order geometric response is verified against finite differences.

**Not done:** physical global two-horizon matching, diffeomorphism-invariant radiation algebra and entropy, QES extremization, or a fine-grained Page transition. Novelty is **not established**. None of the test responses should be presented as a physical QES displacement or a Page-time correction.

Run `python src/two_patch.py`, `python src/interval_entropy.py`, and `python -m pytest -q`.

## Update 4: explicit radiation subalgebra, two-null UV negative control and peeling entropy response

**Local validation:** 8 new tests passed in this update. Earlier 29-test result was reported for the preceding research snapshot and has not been rerun here as a full clone; GitHub Actions result is unconfirmed.

- [Detailed scientific audit and derivations](notes/radiation_algebra_and_negative_control.md)
- [Outgoing chiral subalgebra entropy](src/radiation_algebra.py)
- [Both-null exterior-interval trial extremization](src/static_patch_candidate.py)
- [Conditional time-dependent peeling-map entropy response](src/peeling_response.py)
- [Eight new local independent tests](tests/test_radiation_qes_negative_control.py)

**Conditional numerical results** (`G=1,Lambda=3e-6,M=100,c=12,b=500`):
- Outgoing chiral late-time entropy rate: `0.0020769186253272033`, distinct from Easson's **thermodynamic proxy** `0.0015999713235596304`. Different entropy observables; no Page-time comparison.
- Simultaneous two-null stationary point for a one-interval exterior **trial** functional: `a=499.99872675744`, only `0.00127324` from the artificial boundary at `b=500`. Proper separation `0.00215217` is **smaller than cutoff=1**, so this is NOT a controlled island. The parameter-specific scan found no other exterior root.
- With the separately assumed peeling map `U'(u)=exp[-int_0^u kappa_b(M(s)) ds]`, `M(u)=M0-epsilon J0 u`, `T=20000`, the outgoing **coordinate-cutoff** interval entropy has analytic first-order response `dS/depsilon=0.00177339084435485`. Direct quadrature and an independent ODE agree at `epsilon=1` to numerical precision. This is **not a QES displacement or Page-time shift**.

The missing genuine collector/bath coupling, global quantum state, multi-interval twist correlator and replica channel must be specified before comparing island and no-island branches. **Novelty has not been demonstrated.**

## Update 5: targeted novelty audit and exact free Dirac bath benchmark

**Scientific status:** No original SdS QES theorem or physical fine-grained Page-time shift is established. Existing literature has two-horizon SdS islands (2023), frozen SdS future-boundary islands (2022/24), universal evaporating-Unruh QES results (2025), and Polyakov-corrected Schwarzschild Page curves (2026). Easson (2026) explicitly defers the fully backreacted QES, but that alone is not a novelty proof.

- [Prior-art audit, formulas, numerical values and blockers](notes/novelty_and_free_dirac_bath_audit.md)
- [Exact free Dirac multi-interval and independent lattice correlation matrix](src/free_dirac_multiinterval.py)
- [Auxiliary nongravitating conformal bath, smoothly matched to static SdS](src/transparent_bath_benchmark.py)
- [Both-null trial extremization and UV-artifact rejection](src/bath_trial_extremum.py)
- [Automated free-field, bath and root checks](tests/test_free_dirac_bath.py)

**Reference values** (G=1, Lambda=3e-6, M=100, c=12): free Dirac MI between static intervals [300,350] and [450,500] = 0.5002485802218247; exact two-interval entropy geometric first-order response in the *prescribed* state = 1.7323804001387145e-6 (not an island correction). In the auxiliary bath, MI between gravity [320,400] and bath [20,100] = 0.38294516896267033. The toy two-null extremum a=399.9984084483509 is only 0.00272949 proper units from an artificial endpoint, far **below** the UV cutoff of 1, and is rejected as a non-island artifact.

The complete local research snapshot passed **54 tests** including 17 new tests (7 exact multi-interval, 7 smooth-bath, 3 two-null control). A compact GitHub regression suite is also committed, but its GitHub Actions outcome is not independently verified. The snapshot archive is available as a conversation artifact.

Run `python src/free_dirac_multiinterval.py`, `python src/transparent_bath_benchmark.py`, `python src/bath_trial_extremum.py` and `python -m pytest -q`.

**Next physical requirement:** a genuine gravity–bath interface with conserved stress and specified boundary terms; a globally defined quantum state and purifier; the behind-horizon geometry and replica channel; two-null QES and same-algebra Page-saddle comparison. Only then assess novelty of the *resulting correction*.

## Refreshed literature novelty check — 29 October 2026

[Read focused 2026 novelty audit, primary source comparison, and falsifiable next milestones](notes/novelty_audit_2026-10-29.md).

The original SdS **island** phenomenon, free-fermion multi-interval formula, universal near-horizon QES rule, generic first-order response, and standard EF field constraints are all **known or not yet original**. Easson (2026) explicitly defers full nonequilibrium two-horizon QES extremization. The narrow conditional opportunity is a **physical, UV-controlled, state-specific and finite-cosmological-horizon-dependent QES saddle/transition correction**, or a sharply delimited **obstruction** to that transition due to breakdown of semiclassical control. Prior de Sitter work (Aalsma/Sybesma 2021; Kames-King/Verheijden/Verlinde 2022) found severe backreaction in related collection setups, so this must be confronted and not glossed over. No physical new QES or Page shift is yet established.

## Refined novelty check: existence, dominance and semiclassical validity (29 October 2026)

See [targeted literature re-check](notes/formulation_literature_check_2026-10-29.md). The static near-horizon criterion of Yu & Ge (Chinese Phys. C 2025) already covers classical neutral SdS because f''(r)=-4GM/r³-2Lambda/3<0. Ageev et al. (JHEP 2024) explicitly show that formal island existence is different from island **dominance** in de Sitter setups. Akhmedov/Bazarov (PRD 2023) find strong horizon backreaction for certain thermal states; Shi/Turiaci/Wu (JHEP 2026) study smooth evaporating nucleated dS black holes under different state conditions. **We therefore seek a controlled, state-specific phase classification of QES existence, saddle dominance, and semiclassical validity—not merely a yes/no island.**

No such complete classification in Easson's particular time-dependent UdS Polyakov neutral SdS model was located in the targeted search, but originality is **unproven**, and our code still has no physical behind-horizon Page transition.

## Broad literature sweep and Phase-0 necessary-control scan (30 October 2026)

- [Broad novelty and state/collector protocol literature audit](notes/broad_sweep_and_proxy_control_2026-10-30.md)
- [Exact dimensionless SdS **thermodynamic-proxy** adiabaticity scan](src/proxy_viability.py)
- [Six new independent/symbolic regression tests](tests/test_proxy_viability.py)

**Complete local test result: 60 passed.** The new scan compares the thermodynamic **proxy** crossover with *necessary only* small instantaneous adiabaticity at **both** horizons, and weak (NG/r_b^2) coupling. In the illustrative (q=NG\Lambda=10^{-6}) scan, all five sampled initial conditions satisfy these necessary checks at the *proxy* crossing. This **does not** establish a physical Page curve or an actual controlled island: full collector backreaction, horizon boundary-layer validity, quantum state, replica channels, and two-null behind-horizon QES remain unsolved.

Broader prior art includes the important **controlled finite-duration** 2023 dS–Rindler collector of Aalsma, Aguilar-Gutierrez and Sybesma. Thus neither generic controlled de Sitter information recovery nor a universal de Sitter no-go is a new claim. A **model-specific finite-\Lambda** QES existence/dominance/control result in Easson's neutral time-dependent UdS background remains a *plausible but unverified* opportunity.

## Phase 1 — exact finite-window Gaussian collector and purification (30 October 2026)

**Scientific status: physically incomplete collector, NOT a QES or Page curve.** [Detailed model, results, prior art and limitations](notes/finite_collector_phase1_2026-10-30.md).

- [Gaussian free-fermion mixing and exact external-work energy balance](src/finite_collector.py)
- [Canonical purifier and collector–hot reference mutual information](src/purified_collector.py)
- [Energy conservation, switching, finite-size, and information tests](tests/test_finite_collector.py) · [purification tests](tests/test_purified_collector.py)

A finite-support mixing profile requires external switching work: `Delta E_hot + Delta E_collector = W_switch`, verified to machine precision. For the example reference SdS temperature ratio `kappa_b/kappa_c=2.77135` and a 40-bin half-transmitting Gaussian window, the one-species lattice collector gains `1.3557780` energy units, the switch supplies `0.0378271`, and mutual information with a **chosen** hot-channel Gaussian purifier is `13.0379724` nats. The collector's mean energy is **not** enough to infer its correlations: at equal channel temperature, a uniform rotation gives zero energy transfer and zero work but positive hot-purifier MI.

**Combined local snapshot: 77 tests passed** (previous 60 plus 17 new tests; remote GitHub Actions not separately confirmed). This is a regulated one-shot/time-bin model, not yet a causal continuum defect or a gravitational interface; the ancilla is NOT the black-hole interior. A physical behind-horizon QES, replica saddle comparison and Page-time shift remain unsolved. The result is not claimed as novel fundamental physics.

## Phase 2: causal continuum chiral collector and switching-energy audit (30 October 2026)

**Status:** analytically derived and independently verified *free-field* continuum scattering, **not a gravitationally coupled SdS collector, fine-grained Page transition, QES or novelty claim**.

- [Exact derivations, circle-mode checks, numerical results and missing physics](notes/continuum_collector_phase2_2026-10-30.md)
- [Continuum finite-duration free-Dirac scattering and independent Fourier sideband solver](src/continuum_collector.py)
- [Nine analytical and independent spectral tests](tests/test_continuum_collector.py)

For smooth compact-supported mixing `R(theta(u))`, the per-complex-Dirac-chirality thermal current `J_beta=pi/(12 beta^2)` transforms with switching-energy addition `theta_dot^2/(4pi)` in **each** output channel. Total external work is `W=integral theta_dot^2/(2pi) du`, and `Delta E_C=integral sin^2(theta)(J_H-J_C) du+W/2`. We independently recover the result from Fourier-sideband scattering on an antiperiodic circle, correctly treating the Dirac sea (a naively truncated negative-energy band creates fake divergent work).

The thermal ratio `beta_C/beta_H=kappa_b/kappa_c` reproduces Easson's **2D** flux `N*(kappa_b^2-kappa_c^2)/(48pi)`. For a declared time-window example, switching work is about **28.8% of transferred thermal energy**, so a physical gravity–bath interface cannot ignore the external source.

**Full local regression:** 86 tests passed (77 earlier + nine Phase-2 tests). A causal chiral boundary-scattering prescription is not a complete counterpropagating SdS field theory; curved-space contact action, defect Ward stress, gravitational junction, purification, *behind-horizon* island and saddle comparison remain unsolved. Remote GitHub Actions result has not been confirmed.
