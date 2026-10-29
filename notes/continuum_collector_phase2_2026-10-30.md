# Phase 2 — continuum chiral Dirac collector, verified against Fourier modes

**30 October 2026.** This is a finite-duration, local-in-retarded-time two-channel free-field scattering boundary condition. **It is not a physically solved SdS gravity–bath junction, global UdS density matrix, behind-horizon QES or Page transition. Novelty is not claimed.**

## Connection to literature

- Easson, *Phys. Rev. D* **113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk : anomaly flux and thermodynamic Page proxy, without a microscopic collector.
- *Quantum transport phenomena induced by time-dependent fields*, *Rivista del Nuovo Cimento* (2026), https://doi.org/10.1007/s40766-025-00074-3 : time-dependent scattering mixes energy sidebands and requires accounting for external driving.
- Meineri, Penedones, Rousset, JHEP **02** (2020) 138, https://doi.org/10.1007/JHEP02(2020)138 : transmission through conformal interfaces (stationary framework). Our time-dependent point contact is *not* claimed to remain a stationary conformal defect.

## Local scattering condition and continuum energy identity

For two **unfolded co-propagating** massless complex chiral Dirac fields, impose
```text
(Psi_H^out, Psi_C^out)^T = R[theta(u)] (Psi_H^in,Psi_C^in)^T,
R(theta) = [[cos(theta),sin(theta)],[-sin(theta),cos(theta)]].
```
The mixing angle is smooth with compact temporal support. This is a local causal chiral boundary **scattering condition**; a point-contact Hamiltonian and regularization implementing it in curved SdS space have not yet been fixed.

For a Gibbs state, the per-species incoming currents are `J_i=pi/(12 beta_i^2)`. The output two-point function is `C_ab^out(u,u')=sum_i R_ai(u) C_i(u,u') R_bi(u')^*`. Point-splitting of the universal short-distance singularity gives
```text
J_H_out = cos^2(theta)*J_H+sin^2(theta)*J_C+theta_dot^2/(4pi)
J_C_out = sin^2(theta)*J_H+cos^2(theta)*J_C+theta_dot^2/(4pi).
Q       = integral sin^2(theta)*(J_H-J_C) du
W       = integral theta_dot^2/(2pi) du
Delta_E_H = -Q+W/2
Delta_E_C = +Q+W/2
Delta_E_H+Delta_E_C = W.
```
For N independent identically coupled species multiply Q, W, and both changes by N. This is standard free-fermion stress/transport physics, **not** an island theorem.

## Independent Fourier-mode and Dirac-sea regularization check

On an antiperiodic spatial circle L, frequencies are `omega_m=2pi(m+1/2)/L`. If `c_q,s_q` are Fourier coefficients of `cos(theta),sin(theta)`, then
```text
Delta_f_H(m) = sum_q [
    |c_q|^2*(f_H(omega_(m-q))-f_H(omega_m))
  + |s_q|^2*(f_C(omega_(m-q))-f_H(omega_m))]
Delta_E_H=sum_m omega_m*Delta_f_H(m).
```
Similarly for C, swapping H and C. **Do not** sharply truncate both incoming and outgoing filled Dirac seas to the same spectral band: that generates a spurious UV-edge switching-energy divergence. Instead, sum every incoming sideband m-q and increase the outgoing cutoff. At zero temperature the identity `sum_m omega_m*(f0(omega_(m-q))-f0(omega_m))=pi*q^2/L` and Parseval give
```text
W = (2pi/L)*sum_q q^2*(|c_q|^2+|s_q|^2)
  = (1/(2pi))*integral theta_dot^2 du.
```
Use the exact **finite-circle** NS thermal currents `J_beta^(L)=(1/L)sum_m omega_m*(f_beta(omega_m)-Theta(-omega_m))` for mode benchmarks, not infinite-line currents. The common NS Casimir term cancels in the temperature **difference**.

## Reproducible results

Pulse period L=24, support [4,20], smooth C-infinity ramps of width 3, peak angle 0.7, beta_H=2 and beta_C=7:
```text
Infinite line: Q=0.3121814286751057, W=0.08517469464265984
Delta_E_H=-0.26959408135377577
Delta_E_C=+0.3547687759964356

Finite-circle corrected prediction:
Delta_E_H=-0.2696080641774304
Delta_E_C=+0.3547827588200902
W=0.08517469464265984

Independent Fourier sideband result (m<=160, q<=120, grid 4096):
Delta_E_H=-0.2696080642192624
Delta_E_C=+0.35478275877821763
W=0.08517469455895521
Absolute disagreements < 1e-9.
```

A Cauchy–Schwarz switch-work bound for reaching angle theta0 in ramp durations tau_on/off is `W>=(theta0^2/(2pi))*(1/tau_on+1/tau_off)`. Under time dilation of a fixed pulse, the transported heat Q scales as T, while W scales as 1/T and W/Q scales as T^-2. This is a generic bound, not a new physics theorem.

## SdS temperature-ratio match and unresolved junction

For G=1, Lambda=3e-6, M=100, N=12:
```text
kappa_b=0.0020769186253272033, kappa_c=0.0007494249985680072
J_Easson = N*(kappa_b^2-kappa_c^2)/(48pi)=2.98570944856994e-7.
```
Choosing dimensionless x=(kappa_b/pi)*u gives beta_H,x=2, beta_C,x=2*kappa_b/kappa_c. Then the free-Dirac thermal current difference exactly matches Easson's **2D** current after multiplying by N*(kappa_b/pi)^2. For the sample pulse in x:
```text
Killing-time window = 36302.926251757504
thermal transferred energy Q = 0.0023456322760689852
external switching work W = 0.0006757117012267976
collector total energy gain = 0.0026834881266823837
W/Q = 0.28807230703664005.
```
The *unmodified* Easson mass-loss law during this window would give Delta M/M~1.0839e-4, a necessary slow-evolution scale check only.

**Crucial limitations:** The actual SdS outgoing BH and incoming cosmological streams have opposite propagation directions before unfolding; no physical covariant folding, contact action, defect Ward identity, dilaton-gravity junction, input purification, or replica saddle is solved. A detector/work source changes stress near the interface; transmitted energy is **not** an extra term in the intrinsic BH mass loss without a coupled field-equation derivation. The local Killing-current identification also needs 2D/4D reduction and proper-time redshift treatment.

## Next physics target

Derive a counterpropagating curved-space defect action and stress-energy junction with switch-work source, then solve the anomaly-backreacted two-patch geometry and calculate a *defined* radiation algebra's fine-grained entropy. Only with a UV-separated behind-horizon two-null QES and same-algebra on-shell competing saddles can the existence/dominance/validity phase diagram or a Page-time shift be claimed.

Local entire-repository regression: **86 tests passed** (77 earlier, nine new) on the working snapshot. Remote GitHub Actions not yet confirmed.

Run `python src/continuum_collector.py` and `python -m pytest -q`.
