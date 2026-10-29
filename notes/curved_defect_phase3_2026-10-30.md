# Phase 3 — SdS worldtube contact, localized Killing Ward identity and conditional dilaton shell

**30 October 2026. Scientific status:** model-defined free-field contact and necessary gravitational conservation constraints, **NOT** a solved gravity–bath junction, horizon-regular backreacted solution, QES, or Page crossing. No novelty claim.

## Primary literature and limits

Easson, *Phys. Rev. D* **113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk , computes the 2D Polyakov two-horizon energy flux and a **thermodynamic Page proxy**. Review of gravitational bath coupling and islands: *Solvable models of quantum black holes: a review on Jackiw–Teitelboim gravity*, Living Rev. Relativity **26** (2023), https://pmc.ncbi.nlm.nih.gov/articles/PMC10390458/. The conformal-interface energy-transport framework predates this calculation, e.g. Meineri–Penedones–Rousset, JHEP02(2020)138, https://doi.org/10.1007/JHEP02(2020)138. The external work source has to be included when a contact is switched. The present thin-shell equations are conditional consequences of established 2D dilaton-gravity field equations, not original theorems.

## 1. Three channels, rather than misidentifying bath and cosmological radiation

Choose static background `f(r)=1-2GM/r-Lambda*r²/3` with `rb<b<rc`, `ds²=f(-dt²+drstar²)`. Three *distinct* chiral inputs have Killing currents

```
J_H=N*kappa_b²/(48pi)     BH outgoing
J_C=N*kappa_c²/(48pi)     cosmological incoming; spectator
J_B=N*pi/(12*beta_B^2)   independent nongravitating bath input.
```

Easson's intrinsic 2D net flux is `J_H-J_C`. We choose `beta_B/beta_H=kappa_b/kappa_c` only as an illustrative thermal-matched input; **B and C remain distinct Hilbert spaces**. A true curved-space bath boundary action and state preparation remain unspecified.

After folding the BH-outgoing and bath-input chiral legs, model a chiral point contact in a locally matched conformal frame by `S_def=-int dt Psi^dagger V(t) delta(x-b) Psi`, where `V(t)=-g(t) sigma_y`. Symmetric contact averaging gives

```
Psi(b+)-Psi(b-)=-i*V(t)/2*[Psi(b+)+Psi(b-)]
S=(I+i V/2)^(-1)(I-i V/2)=R(theta),
theta(t)=2atan(g(t)/2).
```

The cosmological input passes unchanged (3-port `diag[R(theta),1]`). The contact is causal at the level of unfolded massless chiral fields and unitary, but an unfolded scattering prescription is not the whole curved external-bath theory.

## 2. Killing-energy Ward identity

Set `T=sin²theta`, `P_K=N*(dtheta/dt)²/(2pi)`. The chiral output currents obey

```
J_patch,BEFORE=J_H-J_C
J_patch,AFTER =J_patch,BEFORE-T*(J_H-J_B)+P_K/2
Delta J_bath = T*(J_H-J_B)+P_K/2
J_patch,AFTER-J_patch,BEFORE+Delta J_bath = P_K.
```

This is a **flux ledger across branches**. For the folded network, the source is distributional `partial_x J_network=P_K delta(x-x_b)` with branch orientations. In static coordinates `sqrt(-g)=1` for `ds²=-f dt²+dr²/f`. A covariant Killing-energy current `j_K^a=-T^a{}_b K^b` then has a localized source, `nabla_a j_K^a=P_K delta(r-b)`, understood as the integrated branch current with appropriate orientations. It is invalid to continue asserting *one conserved radially constant flux* through the contact while omitting the source.

Proper lapse at b is `A=sqrt(f(b))`: `d tau=A dt`, `E_local=E_K/A`, `J_local=J_K/A²`, `P_local=P_K/A²`. Hence fixed Killing input has b-independent integrated energy in this **fixed-geometry model**, while the detector's locally measured power, energy and acceleration grow as either horizon is approached. The static proper acceleration is `|f'(b)|/(2sqrt(f(b)))`. These are necessary control diagnostics, not curvature bounds in the full corrected model.

## 3. Consistency obstruction: an external controller is indispensable

The matter contact receives positive `P_K`. In a fully diffeomorphism-invariant treatment, including gravitational equations, the combined stress must be covariantly conserved. Therefore the time-dependent prescribed coupling must be sourced by an apparatus/controller with work supply `dE_battery^K/dt=-P_K` (or an explicitly named external reservoir). A bare prescribed source with no compensating stress is **not** a complete solution of the Einstein–dilaton–Polyakov equations. A fully dynamical gravitational apparatus could involve momentum/displacement forces and geometric counterterms not fixed by this scalar energy equation.

The repository contains an exact energy-ledger *illustrative battery* starting with twice the total switch work, integrating `E_K(t)=E_K(0)-int^t P_K ds`. No physical gravitational battery is derived here.

## 4. Optional INTERNAL timelike shell: controlled jump coefficients

This is a separate diagnostic, **not** the external gravity–bath interface. In Gaussian normal coordinates `ds²=dell²-A(ell,t)² dt²`, assuming continuous induced metric and X and no delta-function dilaton coupling, the bulk field equations for spherically reduced gravity give

```
E_tt = -A² * partial_ell² X + nonsingular terms
E_X  = -2 * (partial_ell² A)/A -
       (partial_ell² X)/X + nonsingular terms.
```

Let `epsilon_shell` be the line energy density so `T_tt=A² epsilon_shell delta(ell)`. Using `E_ab=2G T_ab` and `E_X=0` for a shell with no dilaton charge yields

```
[partial_ell X]=-2G epsilon_shell,
[K]=G epsilon_shell/X,       K=(partial_ell A)/A.
```

These coefficients are independently derived symbolically, and checked by smoothing the derivative jumps with tanh profiles and integrating the singular second derivatives. A shell action depending on X changes the second equation. A bath with gravity on **only one** side has boundary rather than two-sided internal-shell conditions; **do not** apply these jumps to it without constructing its boundary action/GHY terms.

With a static-lapse illustrative controller `epsilon_shell(t)=E_battery^K(t)/A`, the jump coefficients evolve as its energy is consumed. This only demonstrates that the controller's gravitational contribution cannot simply be omitted; it is not a time-dependent solved geometry.

## 5. Numerical results and control

Parameters `G=1,Lambda=3e-6,M=100,N=12,b=500`, Phase-2 smooth pulse:

```
Q_K (thermal collector transfer)       +0.0023456322760689844
W_K (external switching work)          +0.0006757117012267977
Delta E_patch, K                       -0.0020077764254555854
Delta E_bath, K                        +0.0026834881266823833
total balance residual                  1.1e-19
redshift sqrt(f(b))                     0.5916079783099616
W_local                                 0.0011421612385233446
proper pulse duration                   21477.10080653789
max local switching power               5.193793030118504e-7
W_K/M                                   6.757117012267977e-6
bath energy/M                           2.6834881266823833e-5
unperturbed mass drift/M                1.0838998992061009e-4
```

These last ratios indicate *small integrated background energy fractions*, not control of delta-source curvature, local proper frame at horizons, the bath junction, or fine-grained QES.

A fixed Killing pulse at b=250,300,500,700,850 gives equal Q_K and W_K, but different W_local and local proper acceleration. This is exactly expected redshift behavior and also demonstrates a limitation of inferring local junction control from global mass drift.

## 6. Reproduction and the next genuine physical task

```bash
python src/curved_defect.py
python src/defect_junction_constraints.py
python -m pytest -q
```

**105 local tests passed** in the complete extracted Phase-2-plus-Phase-3 snapshot (86 preceding + 13 worldtube/contact tests + 6 shell/battery tests). The remote GitHub Actions outcome is not independently verified.

Next step requiring original physical assumptions: write a *single diffeomorphism-invariant* gravity–bath + controller boundary action, including GHY/dilaton boundary and proper-time interface coupling, derive all tangential and normal Ward identities and gravitational boundary equations, and solve for a globally regular time-dependent state. Only then is a behind-horizon QES and Page-transition calculation warranted.

**Novelty:** no publishable original SdS QES or gravitational island phase boundary established. The energy/source bookkeeping is mandatory consistency, not a new island theorem.
