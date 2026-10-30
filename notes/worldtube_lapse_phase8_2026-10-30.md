# Phase 8: SdS proper-time worldtube and two-branch lapse response

**30 October 2026.** Exact *anchored finite-mode worldline* CTP benchmark on the **frozen classical SdS static patch**. **No physical gravitating–nongravitating-bath junction, globally horizon-regular state, full covariant metric/dilaton influence functional, QES or Page transition.** The ingredients are standard general-relativistic redshift, defect Ward identities and Gaussian/Klich quantum dynamics; **originality not established**.

## Literature check

- Easson, PRD **113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk: SdS two-horizon anomaly flux, quasistationary rather than an exact static nonzero-flux solution.
- Lombardo & Mazzitelli, PRD **58**, 024009 (1998), https://doi.org/10.1103/PhysRevD.58.024009: closed-time-path influence action for 2D dilaton gravity, including boundary contributions.
- Billò, Gonçalves, Lauria & Meineri, JHEP **04** (2016) 091, https://doi.org/10.1007/JHEP04(2016)091: displacement operators and normal defect Ward identities.
- Major, PRD **105**, 104050 (2022), https://doi.org/10.1103/PhysRevD.105.104050: energetics of accelerated observers in SdS.
- Klich, arXiv:cond-mat/0209642: exact fermion Fock trace/determinant.

## 1. Proper-time anchored action

For 
\[
ds^2=-f(r)\,dt^2+f(r)^{-1}dr^2,\quad
f(r)=1-\frac{2GM}{r}-\frac{\Lambda r^2}{3},\quad
A(b)=\sqrt{f(b)},\quad d\tau=A(b)\,dt,
\]
model a *finite set of detector/contact modes* by
\[
S_{\rm contact}=\int d\tau\, c^\dagger[i\partial_\tau-h_0-g(\tau)V]c .
\]
At fixed radius, \(H_K(t)=A(b)[h_0+g(A(b)t)V]\); a Gibbs state has inverse Killing temperature \(\beta_K=\beta_{\rm proper}/A\). For a properly transformed drive, **the CTP influence determinant is identical in both time coordinates**. Killing energies and work scale as \(A\), and the Killing contact-force noise and retarded kernels scale as \(A^2\) at corresponding times.

The exact finite-mode determinant under independent **lapse histories** is
\[
\boxed{\mathcal F[A_+,A_-;g_+,g_-]
=\det\!\left[1-C_0+C_0U_-^\dagger U_+\right]},
\quad
U_\pm=\mathcal T\exp[-i\int dt\,A_\pm(t)(h_0+g_\pm(t)V)] .
\]
This is only the **clock/lapse restriction** of a metric influence functional, not the full \(\Gamma[g_{ab}^\pm,X^\pm]\).

For no contact drive, \([C_0,h_0]=0\) and \(\Delta\Theta=\int dt(A_+-A_-)\),
\[
\mathcal F=\prod_j[1-f_j+f_j e^{-i\epsilon_j\Delta\Theta}],\quad
\partial_{\Delta\Theta}\log\mathcal F|_0=-i\langle h_0\rangle,\quad
-\partial^2_{\Delta\Theta}\log\mathcal F|_0=\mathrm{Tr}[C_0h_0(1-C_0)h_0].
\]
For a driven interval \(j\) with \(H_j=h_0+g_jV\) and initial segment state \(\rho_j=U_{j-1}\rho_0U_{j-1}^\dagger\),
\[
\boxed{\left.\partial_{A_{+,j}}\log\mathcal F\right|_{A_+=A_-}
=-i\Delta t\,\mathrm{Tr}(\rho_j H_j).}
\]
The source is the **interacting** instantaneous local energy, not simply \(h_0\). This has been verified by independent finite differences at multiple pulse times.

## 2. Both components of the localized worldline Ward source

A static worldline with unit velocity \(u^a=(1/A,0)\) and local energy \(E(t)\) contributes
\[
T^{ab}_\Sigma=E(t)\,u^a u^b\,A\,\delta(r-b).
\]
In static coordinates \(\sqrt{-g}=1\), the delta-supported divergence is
\[
\boxed{\nabla_aT^{at}_\Sigma=(\partial_tE/A)\delta(r-b),\qquad
\nabla_aT^{ar}_\Sigma=(E A f'(b)/2)\delta(r-b).}
\]
The local signed radial support force is
\[
\boxed{F_{\rm support}^{\hat r}
=E\frac{f'(b)}{2A}
=\left.\frac{\partial}{\partial b}[A(b)E]\right|_{E\,\mathrm{fixed}}.}
\]
It vanishes at the familiar classical SdS static radius \(b_*=(3GM/\Lambda)^{1/3}\). This is a **point-worldline** stress identity: propagating radiation momentum, full interface displacement operator and dilaton surface terms remain uncomputed.

For a **closed** finite fermion contact+battery, switch work increases contact excitation energy by \(W\) while the battery energy decreases by \(W\). Total local energy and, for a fixed worldline, total static support force remain unchanged. Omitting the battery would create a spurious tangent source; claiming an additional total weight without a momentum channel would be incorrect.

## 3. Numerical results, proper-time pulse fixed

For \(G=1,\Lambda=3\times10^{-6},M=100\), the 4+4-mode 24-step detector benchmark with local pulse peak \(g=0.35\) and \(d\tau=0.12\) gives \(W_{\rm proper}=0.0664365959603\).

| radius b | lapse A | signed local acceleration | Killing work |
|---:|---:|---:|---:|
| 250 | 0.3708099244 | +0.0036406793 | 0.0246353491 |
| 350 | 0.5532372263 | +0.0008429052 | 0.0367551981 |
| 500 | 0.5916079783 | -0.0001690309 | 0.0393044202 |
| 700 | 0.4735881273 | -0.0010471512 | 0.0314635831 |
| 850 | 0.2054407028 | -0.0034637328 | 0.0136487810 |

The CTP equality under reparameterization holds to ~\(6\times10^{-16}\). An 18-step driven-lapse numerical derivative \(+0.2926825183840i\) independently matches \(+0.2926825183919i\), absolute discrepancy ~\(1.4\times10^{-11}\). A static example has energy mean \(-2.2485793601\) and variance \(0.3342520943\) in *arbitrary lattice energy units*. The absolute energy depends on its reference and is not ADM/Brown–York mass.

The radius comparison fixes the **proper-time pulse and proper apparatus Hamiltonian**, not the Killing-time pulse. Holding the latter fixed instead represents a different experiment.

## 4. Scope and follow-up

The committed [compact worldtube solver](../src/worldtube_phase8.py) and [independent regression tests](../tests/test_worldtube_phase8.py) cover the CTP lapse response, static radial force, redshift and energy cumulants. A more extensive downloadable local snapshot also contains separate full worldtube and lapse modules plus 23 detailed tests; the complete expanded **local** snapshot passes **197 tests**, including six compact integration tests, while the base Phase-7 snapshot had 168. The remote CI run has not been verified.

Unsolved: vary the spatial metric, dilaton and worldline embedding in a **propagating curved-space free-Dirac action**; derive the normal displacement/field momentum source and its quantum noise, including controller/strut and state-dependent Polyakov GHY terms; solve the global gravity–bath junction and retarded semiclassical geometry; then seek UV-controlled behind-horizon two-null QES and a saddle-dominance Page crossing.

**Not a QES discovery or a new geometric theorem.**
