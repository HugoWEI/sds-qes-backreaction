# Phase 4 — Classical boundary variation and the missing SdS shell support force

**30 October 2026. Status:** classical GHY/Brown–York and **two-sided internal-shell** constraints only; **not** a solved external nongravitating bath, semiclassical Polyakov boundary problem, replica island, or novelty claim. Complete local snapshot: **125 tests passed**, 20 additional relative to Phase 3.

## Literature audit

The gravitational boundary variation and quasilocal energy here are standard physics: Brown & York, *Phys Rev D* **47**, 1407 (1993), https://doi.org/10.1103/PhysRevD.47.1407 ; canonical 1+1 gravity boundary energy, *Ann Phys* **258**, 37 (1997), https://doi.org/10.1006/aphy.1997.5687 ; general 2D dilaton GHY/Hamilton-Jacobi counterterms (Grumiller & McNees). A **particularly relevant 2026 predecessor** is Barenboim, Frolov & Kunstatter, *Quantum-corrected Bondi mass for 2D Hawking radiation*, *Phys Rev D* **114**, 024067, https://doi.org/10.1103/q7hb-xjtf, which includes Polyakov and GHY boundary contributions. **The present classical results are not novel and omit their quantum corrections.**

## 1. Classical variational action and radial boundary momenta

In our areal-dilaton convention X=r², use one chosen outward-normal sign convention:
```text
S_g = (1/(4G)) integral_M sqrt(-g) [X R +(grad X)^2/(2X)+2(1-Lambda X)]
    +(1/(2G)) integral_boundary sqrt(-h) X K
    -(1/(2G)) integral_boundary sqrt(-h) W(X),
W(X)=2 sqrt(X)*sqrt(1-Lambda X/3).
```
W is a **finite de Sitter reference subtraction** on its valid static patch; another reference changes the absolute energy, not dE/dM.

In Gaussian normal ds²=dell²-A(t,ell)²dt², with outward n=+partial_ell at the right boundary, K=A_ell/A and R=-2A_ellell/A. Integration by parts exactly cancels the GHY second-normal-derivative term, leaving
```text
L_eff=(1/(4G))[2 X_ell A_ell +A X_ell²/(2X)-X_t²/(2X A)+2A(1-Lambda X)].
Pi_A=(X_n-W)/(2G),
Pi_X=A_n/(2G)+A X_n/(4G X)-A W_X/(2G),
delta S_on_shell=int dt [Pi_A delta A+Pi_X delta X] + other field terms.
```
The counterterm convention gives classical E_BY=-Pi_A at a static b:
```text
E_BY(b,M)=(b/G)[sqrt(1-Lambda b²/3)-sqrt(1-2GM/b-Lambda b²/3)]
d E_BY/dM =1/sqrt(f(b,M)).
```
Thus first-order Killing-mass change delta M=sqrt(f)*delta E_BY, as required for the static-worldtube energy normalization. At G=1,Lambda=3e-6,M=100,b=500, A=0.5916079783099616 and E_BY=137.2087127372385.

**Boundary-condition warning:** fixing BOTH induced lapse A and dilaton X as Dirichlet data makes their boundary variations vanish; it does **not** supply a dynamical gravity–nongravitating-bath junction. A coupled problem must allow appropriate matching variations and include the bath, controller, contact, Polyakov boundary functional and gravitational counterterms. Under tangential diffeo xi(t), delta A=(A xi)_t, delta X=xi X_t; the gravitational boundary variation contributes xi[-A(Pi_A)_t+Pi_X X_t] after integration by parts. The FULL matter+gravity+controller system must satisfy the Ward identity, including incoming/outgoing flux, not the isolated gravity boundary term.

## 2. Exact *internal* shell and second normal-force condition

For two **classical SdS gravity** regions separated by a static timelike shell at b, with shell action S_wall=-integral d tau m(X), X=b², and the SAME Lambda on both sides, choose one continuous increasing-r normal. This is not the external nongravitating bath construction.

Denote A_±=sqrt(f_±(b)), K_±=f_±'(b)/(2A_±). Distributional metric and dilaton equations require:
```text
[X_n]=2b(A_+-A_-)=-2G m,
-2[K]-[X_n]/X=4G m_X,  m_X=d m/d X.
```
The first yields an exact nonlinear gravitational mass dressing:
```text
A_+=A_--Gm/b,
Delta M=M_+-M_-=A_- m-Gm²/(2b).
```
The **independent** dilaton/normal-force equation requires
```text
m_X^required=m/(4b²) [1-(1-Lambda b²)/(A_+ A_-)].
```
A source-free/dilaton-neutral m_X=0 static internal shell **does not exist at a generic prescribed b**, even if the first metric-jump equation holds. Independent radial-equilibrium differentiation with fixed M_± yields exactly the same required force:
```text
m_X= -[m A_-'(b)+Gm²/(2b²)]/(2b A_+).
```
SymPy verifies both forms agree; numeric checks verify the full jump equations.

The force-free shell can exist at a special classical static balance radius. At m->0 this tends to the **known** SdS static radius b*=(3GM_-/Lambda)^(1/3); small m shifts it by Gm A_-(b*)/(2 Lambda b*²). For reference M_-=100,Lambda=3e-6,m=0.1: b*=464.1588833613, b_equilibrium=464.2048847807. At b=500, m_X_required=+2.8547273141e-8; neglecting support leaves a dilaton delta-equation residual +1.1418909256e-7. At b=300, the required dilaton charge has opposite sign. This is classical force-balance, **not a quantum island**.

## 3. Correcting the fixed-background battery ledger

Phase-3 battery Killing energy E_K=A_- m and switch work W=E_K(initial)-E_K(final) are correct **within the fixed-background collector benchmark**. A static self-gravitating internal shell instead has
```text
Delta M=E_K-G E_K²/(2b A_-²),
d Delta M / d E_K = A_+/A_-,
W-[Delta M(initial)-Delta M(final)]
  = G [E_K(initial)²-E_K(final)²]/(2b A_-²).
```
In the reference pulse with initial battery energy 2W, W=0.000675711701226798 and static-shell mass loss=0.000675707787629914. The binding-energy gap=3.91359688435554e-9. It is not a conservation-law violation: a real dynamical solution must include the evolution of binding energy and the supporting force. **Comparing static snapshots is not solving a time-dependent shell or coupled collector.**

## 4. Reproduction and next physical step

Code: [boundary momenta](../src/boundary_variation.py), [exact shell/force and battery dressing](../src/static_shell_support.py), [GHY/BY regression tests](../tests/test_boundary_variation.py), [shell/force regression tests](../tests/test_static_shell_support.py). Complete local regression: **125 passed**.

The outstanding *physical* step is a **single covariant action with the Polyakov boundary/auxiliary field**, the nongravitating bath (not an internal gravity shell), proper-time/embedding matching, independent Dirac fields, and a dynamically supported controller. Derive both tangential and normal Ward identities, solve the horizon-regular time-dependent geometry, and only then attempt a behind-horizon QES / saddle dominance comparison. These classical results do not establish novelty or a Page transition.
