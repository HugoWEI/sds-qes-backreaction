import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
import sympy as s
from baseline import horizons,kappa,flux
from linearized_geometry import (f0, f0_prime, f0_double_prime, stress,
                                 radial_solution, radial_solution_ode,
                                 local_metric_check)


def test_symbolic_ef_metric_components():
    v,r,L=s.symbols("v r L",real=True)
    F=s.Function("F")(v,r);q=s.Function("q")(v,r);z=s.exp(q)
    g=s.Matrix([[-z*z*F,z],[z,0]]);gi=s.simplify(g.inv());coords=(v,r);X=r*r
    Gamma=[[[s.simplify(sum(gi[i,k]*(s.diff(g[k,a],coords[b])+
        s.diff(g[k,b],coords[a])-s.diff(g[a,b],coords[k]))/2
        for k in range(2))) for b in range(2)] for a in range(2)] for i in range(2)]
    def nab(a,b):
        return s.diff(X,coords[a],coords[b])-sum(Gamma[k][a][b]*s.diff(X,coords[k]) for k in range(2))
    bx=sum(gi[a,b]*nab(a,b) for a in range(2) for b in range(2))
    grad=sum(gi[a,b]*s.diff(X,coords[a])*s.diff(X,coords[b]) for a in range(2) for b in range(2))
    def e(a,b):
        return s.simplify(g[a,b]*bx-nab(a,b)+s.diff(X,coords[a])*s.diff(X,coords[b])/(2*X)-g[a,b]*grad/(4*X)-g[a,b]*(1-L*X))
    P=L*r*r+r*s.diff(F,r)+F-1
    assert s.simplify(e(1,1)-2*r*s.diff(q,r))==0
    assert s.simplify(e(0,1)-z*P)==0
    assert s.simplify(e(0,0)+z*z*F*P+r*z*s.diff(F,v))==0


def test_independent_dilaton_conservation_identity():
    r,N,k=s.symbols("r N k",real=True)
    f=s.Function("f")(r)
    tuu=-N/(192*s.pi)*(s.diff(f,r)**2-2*f*s.diff(f,r,2)-4*k**2)
    tuv=N*f*s.diff(f,r,2)/(96*s.pi)
    trr=4*tuu/f**2;tvr=-2*(tuu+tuv)/f
    assert s.simplify(2*s.diff(tvr,r)+3*s.diff(f,r)*trr+2*f*s.diff(trr,r))==0


def test_flux_and_trace():
    for r in [350.,400.,500.,600.]:
        t=stress(r);j=flux(100.,3e-6,1.,12)
        assert np.isclose(t["Tvv_ef"]+f0(r)*t["Tvr"],-j,rtol=1e-12,atol=1e-16)
        R0=4*100/r**3+2*3e-6/3
        assert np.isclose(-4*t["Tuv"]/f0(r),12*R0/(24*np.pi),rtol=1e-12)


def test_independent_radial_integrations():
    for r in [300.,500.,600.]:
        assert np.allclose(radial_solution(r),radial_solution_ode(r),
                           rtol=1e-9,atol=1e-13)


def test_all_first_order_metric_equations():
    for r in [300.,400.,500.,600.]:
        q=local_metric_check(r)
        for component in ["rr","vr","vv"]:
            assert np.isclose(q["E_"+component],q["T_"+component],rtol=1e-10,atol=1e-13)


def test_independent_radial_finite_difference():
    r=500.;h=.005
    s0,F0=radial_solution(r)
    sp,Fp=radial_solution(r+h)
    sm,Fm=radial_solution(r-h)
    q=local_metric_check(r)
    assert np.isclose((sp-sm)/(2*h),q["T_rr"]/(2*r),rtol=3e-7)
    assert np.isclose(((r+h)*Fp-(r-h)*Fm)/(2*h),q["T_vr"],rtol=3e-7)


def test_black_hole_regular_limits():
    rb,_=horizons(100.)
    q=stress(rb+.01)
    assert np.isclose(q["Tvr"],-12*f0_double_prime(rb)/(48*np.pi),rtol=3e-4)
    assert np.isclose(q["Trr"],12*(12*100/rb**4)/(48*np.pi*f0_prime(rb)),rtol=2e-4)


def test_cosmological_patch_asymptotics():
    _,rc=horizons(100.);_,kc=kappa(100.)
    j=flux(100.,3e-6,1.,12);delta=.1
    q=stress(rc-delta)
    assert np.isclose(q["Tuu"],j,rtol=1e-7)
    assert np.isclose(q["Tvr"]*delta,-j/kc,rtol=5e-4)
    assert np.isclose(q["Trr"]*delta**2,j/kc**2,rtol=5e-4)
