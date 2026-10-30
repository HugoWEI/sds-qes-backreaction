# Phase 5: Polyakov anomaly boundary condition in Unruh–de Sitter

**30 October 2026.** This is a fixed-background state/variational consistency calculation, NOT a physical gravity–bath junction, QES, Page curve, or original theorem. The full local suite passes **137 tests** (125 previous plus 12 new).

## Primary sources and novelty warning

- Barenboim, Frolov & Kunstatter, *Quantum-corrected Bondi mass for 2D Hawking radiation*, PRD **114**, 024067 (2026), https://doi.org/10.1103/q7hb-xjtf, gives the complete semiclassical 2D Polyakov+GHY Hamiltonian framework. Our scheme-dependent bare energy is **not** their generalized Bondi mass.
- Fitkevich, Levkov & Zenkevich, JHEP **06** (2020) 184, https://doi.org/10.1007/JHEP06(2020)184, derive Polyakov auxiliary boundary conditions in RST/JT **reflective** 2D models; these are not a transparent SdS bath.
- Easson, PRD **113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk, fixes our anomaly-flux convention but does not supply a microscopic collection model.

## Localized action, additional GHY term and momenta

In areal-dilaton gauge X=r², let

\[
F(X,\chi)=\frac{X}{4G}-\frac{N\chi}{48\pi},
\qquad W(X)=2\sqrt{X}\sqrt{1-\frac{\Lambda X}{3}}.
\]

We use the **specified** local auxiliary action

\[
S=\int_{\mathcal M}\sqrt{-g}\left[
FR+\frac{(\nabla X)^2}{8GX}
+\frac{1-\Lambda X}{2G}-\frac{N}{96\pi}(\nabla\chi)^2\right]
+2\int_{\partial\mathcal M}\sqrt{|h|}FK
-\frac1{2G}\int_{\partial\mathcal M}\sqrt{|h|}W(X).
\]

The required additional Polyakov Dirichlet-metric boundary term is
\(-N/(24\pi)\int\sqrt{|h|}\chi K\). The finite W subtraction is classical, not a completed quantum state-dependent counterterm. Variation gives bulk equation \(\Box\chi=R\). The auxiliary-field Green function and boundary conditions encode the state; \(\chi\) is not an independent collectible radiation degree of freedom.

For \(ds^2=d\ell^2-A^2dt^2\) in Gaussian normal coordinates, \(R=-2A_{\ell\ell}/A\).
The GHY term exactly removes the second normal derivatives and yields

\[
\mathcal L_{\rm GN}=
2F_\ell A_\ell+\frac{AX_\ell^2}{8GX}-\frac{X_t^2}{8GXA}
+\frac{A(1-\Lambda X)}{2G}
-\frac{NA\chi_\ell^2}{96\pi}+\frac{N\chi_t^2}{96\pi A}.
\]

With increasing-r outward normal, the radial momenta are

\[
\boxed{
\Pi_A=\frac{X_n-W}{2G}-\frac{N\chi_n}{24\pi},
\quad
\Pi_X=\frac{A_n}{2G}+\frac{AX_n}{4GX}-\frac{AW_X}{2G},
\quad
\Pi_\chi=-\frac{N}{48\pi}(2A_n+A\chi_n).
}
\]

Independent symbolic differentiation of the complete improved Lagrangian validates all three. If \(\chi\) is **freely varied** at a boundary with no additional field action, \(\Pi_\chi=0\), equivalently \(\chi_n=-2K\). **Dirichlet** \(\delta\chi=0\) does not imply this Neumann condition. Neither can be declared a physically transparent collector without a microscopic state construction.

## SdS state and exact boundary-flux identity

For frozen classical \(f(r)=1-2GM/r-\Lambda r^2/3\), \(u=t-r_*\), \(v=t+r_*\), choose the two-future-horizon-adapted auxiliary solution

\[
\boxed{\chi=-\log f-\kappa_bu-\kappa_cv}.
\]

It gives \(\Box\chi=R=-f''\) and the previous Polyakov state stress components:

\[
T_{uu}=\frac{N}{192\pi}(2ff''-f'^2+4\kappa_b^2),\quad
T_{vv}=\frac{N}{192\pi}(2ff''-f'^2+4\kappa_c^2),\quad
T_{uv}=\frac{Nff''}{96\pi}.
\]

The trace is \(NR/(24\pi)\). The state is regular in the *appropriate future* EF charts from the static patch:
\((\partial_r\chi)_{v,BH}=(2\kappa_b-f')/f\to-f''(r_b)/f'(r_b)\) and
\((\partial_r\chi)_{u,cosmo}=(-2\kappa_c-f')/f\to-f''(r_c)/f'(r_c)\).
It is NOT asserted regular at the past or bifurcate horizons.

At any \(r_b<b<r_c\):

\[
A=\sqrt f,\quad A_n=f'/2,\quad X_n=2b\sqrt f,\quad
\chi_t=-(\kappa_b+\kappa_c),\quad
\chi_n=(\kappa_b-\kappa_c-f')/\sqrt f.
\]

Substitution yields the exact fixed-background formulas

\[
\boxed{\Pi_\chi=-\frac{N(\kappa_b-\kappa_c)}{48\pi}},
\qquad
\boxed{\Pi_\chi\,\dot\chi=
\frac{N(\kappa_b^2-\kappa_c^2)}{48\pi}=\mathcal J_{\rm Easson}.}
\]

The first is independent of the observation radius; the second is a **kinematic boundary-state flux factorization, not the full backreacted gravity+bath Ward identity**. It is algebraically implied by established Polyakov theory.

**Conditional incompatibility:** generic neutral nondegenerate SdS has \(\kappa_b\ne\kappa_c\), so this chosen UdS state cannot satisfy the free-auxiliary condition \(\Pi_\chi=0\) without an additional source. It can instead have appropriate Dirichlet state data or a derived bath boundary coupling. A merely formal source \(\int dt\,j_\chi\chi\) would require \(j_\chi=-\Pi_\chi=N(\kappa_b-\kappa_c)/(48\pi)\), contributing \(j_\chi\dot\chi=-\mathcal J\). This source is **not a physical, covariant gravity–bath interface**: its metric, embedding and state-dependent variations have not been derived.

Under boundary time reparametrization
\(\delta_\xi A=(A\xi)_t,\delta_\xi X=\xi X_t,\delta_\xi\chi=\xi\chi_t\),
the on-shell boundary variation has kinematic coefficient
\(-A\dot\Pi_A+\Pi_X\dot X+\Pi_\chi\dot\chi\), up to total derivatives.
On the **frozen** state it equals \(\mathcal J\). The frozen background with nonzero quantum flux does not solve the full semiclassical equations; the actual Ward balance needs dynamical metric, bath and controller contributions.

## Numerical example and renormalization caveat

At \(G=1,\Lambda=3\times10^{-6},M=100,N=12,b=500\):

| Object | Computed value |
|---|---:|
| Auxiliary boundary momentum \(\Pi_\chi\) | -0.000105638586310857 |
| Flux factor \(\Pi_\chi\dot\chi\) | +0.000000298570944856994 |
| Formal source \(j_\chi=-\Pi_\chi\) | +0.000105638586310857 |
| Bare anomalous shift of \(-\Pi_A\) | +0.000410927793663 |

The bare quantity \(E_{\rm bare}=-\Pi_A=E_{\rm BY,classical}+N\chi_n/(24\pi)\) depends on radius, state and the reference/counterterm convention, and must NOT be presented as the unique fully quantum-corrected Brown–York mass. Sign choices for the homogeneous auxiliary-field chiral slopes can leave stress one-point functions unchanged but alter \(\Pi_\chi\); horizon-state selection is necessary.

## Reproducibility and next physics

Implementation: src/polyakov_boundary.py. Twelve independent symbolic/numerical checks: tests/test_polyakov_boundary.py. Run the repository Python tests; **137 passed locally**.

Next: construct the boundary Green function for the **actual transmitted free-Dirac state** on a joined SdS/nongravitating-bath manifold; derive the covariant contact/controller action, including metric and auxiliary-field boundary variations, and solve the dynamical bulk constraints. Only then can a behind-horizon QES and common-algebra island/no-island Page crossing be tested.

**Originality: NOT established.**
