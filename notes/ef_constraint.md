# Explicit leading-order EF radial constraint: initial audit

**Provisional calculation in the 2D Polyakov convention.** This does not prove an inconsistency in Easson (2026), does not solve the full field equations, and does not determine a QES.

Starting with the spherically reduced action as written in Easson's Eq. (5), variation with respect to the metric gives, for the 2D conformal stress defined as `T_ab=-2 delta S_m/(sqrt(-g) delta g^{ab})`,

```
E_ab = g_ab Box X - nabla_a nabla_b X
       + (partial_a X partial_b X)/(2X)
       - g_ab (nabla X)^2/(4X)
       - g_ab (1-Lambda X)
     = 2 G4 T_ab.
```

Use areal-radius gauge `X=r^2` with a **general** Eddington–Finkelstein ansatz

```
ds^2=-exp(2sigma(v,r))*F(v,r)*dv^2+2exp(sigma(v,r))*dv*dr.
```

Symbolic evaluation gives the exact reduced-gravity components

```
E_rr = 2*r*partial_r sigma,
E_vr = exp(sigma)*(Lambda*r^2+r*partial_r F+F-1).
```

Thus the radial constraint is

```
partial_r sigma = G4*T_rr/r.
```

The simpler ansatz `sigma=0`, `X=r^2` would require `T_rr=0`. Outgoing Polyakov radiation generally produces nonzero `T_rr` in an ingoing EF chart. This is a **convention/gauge consistency check that must be reconciled** with the paper's adiabatic ansatz, not grounds for claiming the paper is wrong without a full review.

On an instantaneous static SdS background let `u=v-2r_*(r)`, `dr_*/dr=1/f0(r)`, and `rho=log(f0)/2`. The regular black-hole outgoing chiral choice is `t_u=kappa_b^2/4`; the incoming cosmological sector has `t_v=kappa_c^2/4`. Easson's stress formula yields

```
T_uu=-N/(12pi) [ f0'^2/16 - f0*f0''/8 - kappa_b^2/4 ].
T_rr (in ingoing EF) = 4*T_uu/f0^2.
sigma_1'(r) = G4*T_rr^(0)(r)/r.
```

This is a **local first-order redshift constraint** at fixed instantaneous mass; integration determines `sigma_1(r)-sigma_1(r_ref)`, not the gauge-dependent constant. It is *not* the full metric response: the `vr`, `vv` and dilaton/auxiliary-field equations, state regularity on both horizons and evolution must also be solved. The ingoing EF chart and naive adiabatic expansion are nonuniform at the cosmological horizon; a second regular patch or matched expansion will likely be necessary.

For `G4=1,Lambda=3e-6,M=100,N=12`, at `r=500`:

```
f0=0.35
T_uu=2.7005338935e-7
T_rr=8.8180698563e-6
sigma((rb+rc)/2)-sigma(500) = 7.84977926e-7
```

Direct quadrature and a separately integrated radial ODE for `sigma(500)-sigma((rb+rc)/2)` agree to approximately 1e-15 absolute. SymPy independently reproduces the exact `E_rr` identity.

## Stress/flux normalization issue to check

With a **minimally coupled 2D** Polyakov stress, the stationary conservation law uses the 2D volume density. In static coordinates `T^r_v=T_vv-T_uu=N(t_v-t_u)/(12pi)=-J`, independent of radius. The paper also writes `J=-2r^2 T^r_v` (its Eq. (17)), suggesting that its `T^r_v` might instead denote a suitably reduced/4D-normalized density in that particular expression. We must explicitly reconcile these definitions before importing its flux into a new metric solver. This is an *unresolved normalization audit*, not a confirmed contradiction.

## Next

Derive and solve **all** linearized Einstein–dilaton–Polyakov equations with a consistently normalized 2D stress in regular two-horizon patches. Only then compute conformal-twist entropies of a declared radiation algebra and QES corrections.