# Phase 6 — Transmitted free-Dirac boundary state and purification audit (30 Oct 2026)

**Status:** 155 **local** tests pass (137 before this update + 18 new). This is an explicitly prescribed, frozen-background Gaussian matter **boundary** state. It is NOT a solved curved-space gravitational bath state, black-hole interior purifier, replica saddle, QES or Page time. No original theorem is claimed.

## Primary prior art

- Easson, PRD **113** 084014 (2026), https://doi.org/10.1103/trlf-d6sk — reference SdS two-horizon Polyakov flux; no microscopic radiation-collector QES.
- Gutperle & Miller, *Entanglement entropy at CFT junctions*, PRD **95**, 106008 (2017), https://doi.org/10.1103/PhysRevD.95.106008 — entanglement for transmitting free-boson/fermion junctions.
- Casini & Huerta, *Entanglement entropy in free quantum field theory*, J Phys A **42**, 504007 (2009), https://doi.org/10.1088/1751-8113/42/50/504007 — correlation-matrix free-field entropy.
- Lombardo & Mazzitelli, *Influence functional in two-dimensional dilaton gravity*, PRD **58** 024009 (1998), https://doi.org/10.1103/PhysRevD.58.024009 — **closed-time-path/in-in** treatment, essential for real causal semiclassical equations.
- Potaux, Sarkar & Solodukhin, PRD **105** 025015 (2022), https://doi.org/10.1103/PhysRevD.105.025015 — auxiliary-field state choice affects backreacted geometry.
- Akhmedov, Anempodistov & Bazarov, PRD **109** 065026 (2024), https://doi.org/10.1103/PhysRevD.109.065026 — self-consistent two-horizon dilaton gravity backreaction, important separate precedent.
- Barenboim, Frolov & Kunstatter, PRD **114** 024067 (2026), https://doi.org/10.1103/q7hb-xjtf — quantum-corrected Brown–York/Bondi mass and Polyakov+GHY boundary Hamiltonian.

## 1. Unequal-time boundary Wightman state

The independent hot H and nongravitating bath B chiral input fields have thermal kernels (for real separated s=u-v, in a chosen convention)

```
G_beta(s) = i/[2 beta sinh(pi*s/beta)]
          = i/(2pi*s) * r_beta(s),
r_beta(s)=(pi*s/beta)/sinh(pi*s/beta).
```

With a prescribed smooth finite pulse mixing angle theta(u), the *exact* Gaussian two-time boundary covariance is

```
G_out(u,v) = R(theta(u)) * diag(G_beta_H(u-v),G_beta_B(u-v))
                          * R(theta(v))^T,
R(theta)=[[cos(theta),sin(theta)],[-sin(theta),cos(theta)]].
```

This is **not** the retarded Green function for a globally joined gravity/bath spacetime. The original incoming SdS cosmological C chirality is a *separate spectator*. A switched output is generically nonstationary/non-KMS: it cannot be assumed to have the unchanged steady UdS Polyakov state.

The stable point-split stress for each output a is

```
J_a(t;s) = [1-sum_i R_ai(t+s/2)R_ai(t-s/2)r_beta_i(s)]/(2pi*s²).
```

At s->0 this independently reproduces the Phase-2 result

```
J_Hout = cos²theta J_H+sin²theta J_B+theta_dot²/(4pi),
J_Bout = sin²theta J_H+cos²theta J_B+theta_dot²/(4pi),
J_beta = pi/(12 beta²).
```

The cross-channel Wightman kernel is
`G_HB(u,v)=-cos(theta(u))*sin(theta(v))*G_H(u-v)+sin(theta(u))*cos(theta(v))*G_B(u-v)`.
Its short-distance coincident limit is **i*theta_dot(t)/(2pi)**. An event while the switch is on remains correlated with an event before the switch; an instantaneous local flux value cannot determine this correlation. Correlation across times is consistent with *causal* local scattering; it is not an acausal signal.

## 2. Regulated NS-circle density matrix and entropy

For even n, circle circumference L, u_j=j L/n and frequencies omega_m=2pi(m+1/2)/L for m=-n/2,...,n/2-1,

```
F_jm = exp(-i omega_m*u_j)/sqrt(n)
C_beta = F diag[(1+exp(beta*omega_m))^-1] F^dagger.
```

The full 2n-mode covariance is rotated by an exact bin-diagonal Gaussian SO(2) unitary. The output bath covariance is `C_B'=s*C_H*s+c*C_B*c`, and cross-channel block `-c*C_H*s+s*C_B*c`, with diagonal c,s. The entropy is exactly `S(C)=-Tr[C log C+(1-C)log(1-C)]`, with total H+B entropy preserved.

Introduce a **declared** hot-channel Gaussian canonical purifier R_H with covariance `C_R=I-C_H` and cross `D=sqrt(C_H(I-C_H))`. Its bath-output cross correlation is `C_(B,R_H)=-s D`; the resulting `I(B:R_H)` is not automatically a physical black-hole island entropy.

Reference pulse L=24, beta_H=2, beta_B=7, support [4,20], ramp 3, peak theta=0.7:

| n | bath entropy increase | I(bath:hot reference) | I(hot:bath) |
|---:|---:|---:|---:|
| 32 | 1.9681756955 | 2.6507000702 | 1.8894184903 |
| 64 | 1.9818710389 | 2.6505593960 | 1.9178026713 |
| 96 | 1.9818873017 | 2.6505469661 | 1.9178691771 |
| 128 | 1.9818810265 | 2.6505472308 | 1.9178573510 |

These are nats for **one** Dirac species in a *finite-band antiperiodic circle*, not a horizon Page curve. The global 2-channel entropy residual is order 1e-13; point-split energy agrees with analytic local currents.

## 3. A particularly sharp purifier ambiguity (known quantum information)

Prepare *identical* thermal H and B states with covariance `C_H=C_B=C`. A global constant channel rotation commutes with the physical H+B product thermal state: `diag(C,C)` is unchanged as a full 2n x 2n covariance. Thus **all H+B physical observables and stress components are unchanged** by the mixing. But the correlations with a *fixed labeled* external hot purifier R_H do change: `C_(B,R_H)=-sin(theta)*sqrt(C(1-C))`. With n=32, beta=2 and theta=pi/4, `I(B:R_H)=5.148248758839879` nats versus exactly zero at theta=0.

This is **not novel** as an abstract result. It proves that the Polyakov stress (indeed even the physical H+B density matrix *without the labeled external purifier*) does not specify the desired *purifier-tagged* information flow. A fine-grained SdS island problem must independently fix the global state and the precise interior/radiation algebras.

## 4. Implications and next genuine physics

The Phase-5 auxiliary scalar chi can encode the anomaly and **stress one-point functions** with state-dependent homogeneous solutions; it does **not** by itself uniquely fix the transmitted unequal-time matrix kernel or its external reference correlations. Do **not** infer that a single chi can never reproduce an arbitrary stress transient; the stronger no-go is about **state/purifier determination from stress alone**.

Next:
1. Define a *globally causal*, two-sided, curved-space Dirac boundary initial-value problem and a physical collapse/initial state.
2. Derive a **Schwinger–Keldysh** influence action for the interface, including retarded response and quantum noise; do not identify a naive in-out inverse box with a physical causal bath Green function.
3. Dynamically include collector/controller stress, radial support, the Polyakov boundary Hamiltonian and renormalization in gravity–bath junction equations.
4. Locate regular behind-horizon QESs with **the same physical purifier and radiation algebra** in both replica branches, or prove a state-specific validity obstruction.
5. Redo novelty audit against the established CFT junction and 2D two-horizon backreaction literature before a public claim.

Run `python src/transparent_state_kernel.py`, `python src/gaussian_boundary_state.py`, and `python -m pytest -q`. Full locally assembled repository: **155 passed**. The new GitHub files are committed; remote Actions status not independently confirmed.
