# Phase 7 — Exact finite-mode closed-time-path collector response

**Date:** 30 October 2026. **Local regression:** 168/168 passed: 155 previously established tests plus **13 new** CTP/contact tests. **Status:** a physically explicit *finite-dimensional quantum contact Hamiltonian* with exact causal influence functional and driving-work accounting. NOT a curved SdS–bath influence functional, dynamical dilaton solution, physically identified purifier, QES, or Page curve.

## Literature and novelty

- I. Klich, *Full counting statistics: an elementary derivation of Levitov's formula* (2002), arXiv:cond-mat/0209642, https://arxiv.org/abs/cond-mat/0209642 : the fermionic Fock-space trace / one-particle determinant identity is standard.
- F. Lombardo and F. Mazzitelli, *Influence functional in two-dimensional dilaton gravity*, PRD **58**, 024009 (1998), https://doi.org/10.1103/PhysRevD.58.024009 : in-in/CTP dilaton-gravity influence action and state/boundary effects. Our exact finite-contact determinant **does not** reproduce or replace their geometric effective action.
- D. Easson, *Fate of Schwarzschild–de Sitter black holes*, PRD **113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk, supplies the desired eventual background flux, not the finite-chain response.
- Casini–Huerta 2009 and Gaussian transport literature: finite-fermion covariance, Wick contractions and thermal/KMS response are established techniques.

## 1. Microscopic quantum Hamiltonian instead of an instantaneous gate

There are two independently prepared open finite hopping chains H and B, each of n sites. A Hermitian local contact at site j has one-particle Hamiltonian
\[
h(t)=h_0+g(t)V,\quad h_0=h_H\oplus h_B,\quad
V=\begin{pmatrix}0&iP_j\\-iP_j&0\end{pmatrix},\quad P_j=|j\rangle\langle j|.
\]
Initial covariance
\[
C_0=f_{\beta_H}(h_H)\oplus f_{\beta_B}(h_B),
\quad f_\beta(E)=(1+e^{\beta E})^{-1}.
\]
The one-particle evolution is the chronological product of \(\exp[-i(h_0+g_k V)\Delta t]\). At nonzero local contact, the Hamiltonian does not commute with the free chain: this provides a genuine retarded response and an exact external-work identity. The model is **not a chiral continuum, nor a physical 4D SdS worldtube**. The two input temperatures are independent parameters; bath energy offset distinguishes channel spectra.

## 2. Exact two-branch influence functional

For controller histories \(g_+,g_-\), define
\[
\boxed{\mathcal F[g_+,g_-]
= {\rm Tr}_{\rm Fock}(U_+\rho_0U_-^\dagger)
=\det\big[I-C_0+C_0\,u_-^\dagger u_+\big],}
\]
with \(u_\pm\) the single-particle time-ordered unitaries. This is exactly Klich's known formula. Tests verify:
- \(\mathcal F[g,g]=1\) (CTP normalization);
- \(\mathcal F[g_+,g_-]^*=\mathcal F[g_-,g_+]\);
- \(|\mathcal F|\le1\) within floating error;
- a **separately implemented 16-dimensional full fermionic Fock-space trace** agrees with the determinant to 1e-12.
A forward history with a different return history produces a nontrivial overlap modulus and phase. This captures decoherence of *controller histories* after tracing the finite matter channels; it is NOT the fully renormalized gravitational influence functional.

## 3. Causal linear response and positive noise

At zero drive, \(V_I(t)=e^{+i h_0t}V e^{-ih_0t}\). Connected free-fermion Wick contraction yields
\[
N(t,s) = \operatorname{Re}\operatorname{Tr}
[C_0 V_I(t)(I-C_0)V_I(s)].
\]
The Kubo retarded susceptibility is
\[
\chi_R(t,s)= -i\Theta(t-s)\operatorname{Tr}
\left[C_0\left[V_I(t),V_I(s)\right]\right].
\]
The noise kernel is symmetric, real, and positive semidefinite on sampled grids (up to floating eigenvalues approximately -1.9e-16), and the retarded matrix vanishes before its source times. A separate finite-difference calculation of the exact dynamically evolved final contact force agrees with the independently calculated Kubo integral (relative error <3e-5 for the tested perturbation). Exact analytic integration of interaction-picture force over each time step avoids midpoint integration bias.

For a small single-history probe \(g(t)=\epsilon p(t)\) with the other CTP history zero,
\[
-\frac{2}{\epsilon^2}\log|\mathcal F[\epsilon p,0]|
=\operatorname{Var}\left(\int p(t)\hat V_I(t)\,dt\right)+O(\epsilon),
\]
which is independently evaluated as
\(\operatorname{Tr}[C_0 A(I-C_0)A]\).
For n=4 chains, 28 steps of 0.15 and a unit sin-squared profile, the independently computed variance is \(\approx0.58052429475\). The determinant gives \(\approx0.58052400953\) at \(\epsilon=0.002\), with convergence as epsilon decreases. This is the standard Gaussian influence/noise relation—not new physics.

## 4. Exact switching-work conservation

In a **piecewise constant** coupling history, energy stays constant under each segment's Hamiltonian \(h_0+g_j V\). At an instantaneous coupling jump \(g_{\rm old}\to g_{\rm new}\), the state is unchanged and switching work is \(\Delta W=(g_{\rm new}-g_{\rm old})\langle \hat V\rangle\). Including initial and final jumps to/from zero gives **exactly**
\[
\boxed{W_{\rm controller}
=\Delta E_H+\Delta E_B
=\Delta\operatorname{Tr}(h_0 C).}
\]
For n=5, hopping 0.7, beta_H=2, beta_B=5, bath on-site offset 0.18, 48 steps of 0.12, and peak g=0.45:
- External switch work: **0.06937972616684293** lattice units.
- Hot-channel energy change: **−0.00531914817325560**.
- Bath energy change: **+0.07469887434009767**.
- Energy-balance residual: **8.74e-16**.
These results are for one finite number-conserving fermion species and **cannot be identified with a 2D/4D SdS mass change** without the curved-space contact Hamiltonian and gravitational junction.

## 5. Equilibrium KMS versus nonequilibrium driving

The positive-frequency unsymmetrized force-noise weights for a pair of energy eigenmodes i,j (with \(\omega=E_j-E_i>0\)) are proportional to
\[
w_+ = f_i(1-f_j)|V_{ij}|^2,\quad
w_- = f_j(1-f_i)|V_{ij}|^2.
\]
At common equilibrium inverse temperature beta, \(w_+/w_-=e^{\beta\omega}\) exactly: the tests verify this to <2e-12 in logarithms. When beta_H=2, beta_B=5, there is generally **no single KMS temperature** for this contact force, as checked numerically. The nonequilibrium kernel must not be replaced by a thermal Langevin/noise kernel by fiat.

## 6. Scientific limitations and the true next work

The exact influence functional is for a **finite local Hermitian contact controlled by an external c-number history**, not for metric/dilaton histories \(g_{ab}^{\pm},X^\pm\). Its reservoir is finite and has revivals, not irreversible thermodynamic friction or Markovian noise. It is not a curved, covariantly regulated stress-energy influence functional. The presence of separate hot and bath thermal density matrices still does not identify a physical purifier with a SdS interior. Thus:
1. Derive a covariant, counterpropagating curved-space contact action whose weak-field limit matches a regulated contact Hamiltonian and respects SdS Killing/proper-time normalization.
2. Include controller, normal-support stress, gravitational boundary counterterms and Polyakov state data in **both** Keldysh branches.
3. Construct initial collapse/field state and explicit physically justified purifier/radiation algebra, and compute its causal propagation through the horizon.
4. Solve retarded semiclassical dilaton equations and check quantum stress/boundary noise, without inferring full backreaction from the contact's finite energy alone.
5. Only after that compute two-null **behind-horizon** QES saddles and the same-algebra fine-grained Page transition (or a controlled obstruction).

**Conclusion:** robust causal contact-response infrastructure and known determinant/KMS physics, not an original SdS QES or Page-curve discovery.
