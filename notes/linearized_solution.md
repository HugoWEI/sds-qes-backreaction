# Audited local first-order EF geometry (29 October 2026)

**Status: model-specific local leading-order metric solution, not a completed global Unruh–de Sitter QES or fine-grained Page curve.** Reference: Easson, *Phys. Rev. D* **113**, 084014 (2026), [DOI](https://doi.org/10.1103/trlf-d6sk), [arXiv:2511.11873](https://arxiv.org/abs/2511.11873). The precise two-dimensional versus area-weighted flux convention in the paper still requires reconciliation.

In spherical reduction the metric action is
```
S_g=(1/(4G)) integral sqrt(-g) [X R + (grad X)^2/(2X) + 2(1-Lambda X)].
```
With `T_ab=-2 δS_m/(sqrt(-g) δg^ab)`, the metric equation is `E_ab=2G T_ab`. Use areal-radius gauge `X=r²` and `ds²=-e^(2σ)F dv² +2 e^σ dv dr`. SymPy independently derives **the exact equations**
```
E_rr=2r σ_r,
E_vr=e^σ (Lambda r²+r F_r+F-1),
E_vv=-e^(2σ)F(Lambda r²+r F_r+F-1)-r e^σ F_v.
```

On classical `f=1-2GM/r-Lambda r²/3`, with `u=v-2r_*`, `dr_*/dr=1/f`, `t_u=kappa_b²/4`, `t_v=kappa_c²/4`, the zeroth-order Polyakov stress is:
```
T_uu=-N/(192π)(f'²-2ff''-4κb²)
T_vv(conformal)=-N/(192π)(f'²-2ff''-4κc²)
T_uv=N ff''/(96π)
T_rr(EF)=4T_uu/f²
T_vr(EF)=-2(T_uu+T_uv)/f
T_vv(EF)=T_uu+2T_uv+T_vv(conformal).
```
This fixes `T_vv(EF)+f T_vr=-J`, `J=N(κb²-κc²)/(48π)>0`. This is the **two-dimensional** conserved flux in the declared convention. The consistency of Easson's area-weighted flux notation remains an audit item.

At linear order, evaluated with the classical stress, `σ1_r=G T_rr/r`, `(r F1)_r=2G T_vr`, `F1_v=2G J/r`. The local solution, with reference normalization `σ1(r0)=0`, `F1(v0,r0)=0`, is
```
σ1(r) = G ∫_(r0)^r T_rr(s)/s ds,
F1(v,r) = [2G ∫_(r0)^r T_vr(s) ds + 2G J(v-v0)]/r.
```
The reference normalization fixes the homogeneous mass/redshift integration constants, **not** physical global boundary data; the secular mass shift `δM=-J(v-v0)` requires `|δM| << M`.

Independently, the dilaton equation `E_X=R + (grad X)^2/(2X²)-Box X/X-2Lambda=0` yields, for stationary `σ1`, the linear residual `E_X1=-(G/r)[2(T_vr)' +3f'T_rr+2f(T_rr)']`. **SymPy verifies the bracket vanishes exactly**, an additional constraint beyond the three integrated metric equations.

**Horizon patch behavior:** at the future black-hole horizon in ingoing EF, `lim T_vr=-N f''(rb)/(48π)` and `lim T_rr=N f'''(rb)/(48π f'(rb))`, both finite. Near the future cosmological horizon `δ=rc-r→0+`, the same ingoing EF chart is singular and `T_uu→J`, `T_vr≈-J/(κc δ)`, `T_rr≈J/(κc² δ²)`. Hence `σ1≈GJ/(rc κc² δ)` and `F1≈(2GJ/(κc rc))log(δ)` up to less singular terms. **These are coordinate-patch/perturbative expressions, not evidence of a physical curvature singularity.** A horizon-adapted outgoing EF/Kruskal patch and matching are necessary before defining a cosmological-horizon displacement.

**Numerical regression (G=1, Lambda=3e-6, M=100, N=12, r0=400, v=v0, r=500):** `J=2.98570944856994e-7`, `σ1=1.8264080066815465e-6`, `F1=-3.679751368641078e-7`; `E_rr=2G T_rr=1.7636139712627618e-5`, `E_vr=2G T_vr=-2.2587187456319774e-6`, `E_vv=2G T_vv=1.9340967125720405e-7`. Numerical quadrature and an independently solved radial IVP agree under the regression tolerances.

**Novelty unverified.** Important prior art: Goswami/Narayan (2022, 2024) derive frozen-background SdS islands; Easson (2026) gives a thermodynamic proxy, not an evaluated full-QES transition; Đorđević/Radovanović (2026) calculate evaporating Schwarzschild QES without a cosmological horizon; Baek/Choi (2023) treat JT-like 2D de Sitter backreaction. Next: reconcile source action/flux normalizations; solve two horizon-adapted patches and match; define the radiation algebra, purification and UV subtractions; find simultaneous two-null island saddles and compare on-shell entropies in the same setup.
