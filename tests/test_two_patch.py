"""State and two-chart compatibility: do not mistake this for a QES test."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import sympy as s
from baseline import horizons,kappa,flux
from linearized_geometry import f0,stress,radial_solution
from two_patch import (outgoing_stress,outgoing_radial_solution,
    outgoing_radial_solution_ode,outgoing_field_residuals,coordinate_match,
    tortoise)


def test_exact_outgoing_ef_components():
    x,r,L,sign=s.symbols('x r L sign',real=True)
    F=s.Function('F')(x,r);q=s.Function('q')(x,r);z=s.exp(q)
    g=s.Matrix([[-z*z*F,-z],[-z,0]]);gi=s.simplify(g.inv());co=(x,r);X=r*r
    ch=[[[s.simplify(sum(gi[i,k]*(s.diff(g[k,a],co[b])+s.diff(g[k,b],co[a])-s.diff(g[a,b],co[k]))/2 for k in range(2))) for b in range(2)] for a in range(2)] for i in range(2)]
    def nab(a,b):return s.diff(X,co[a],co[b])-sum(ch[k][a][b]*s.diff(X,co[k]) for k in range(2))
    bx=sum(gi[a,b]*nab(a,b) for a in range(2) for b in range(2))
    gr=sum(gi[a,b]*s.diff(X,co[a])*s.diff(X,co[b]) for a in range(2) for b in range(2))
    def E(a,b):return s.simplify(g[a,b]*bx-nab(a,b)+s.diff(X,co[a])*s.diff(X,co[b])/(2*X)-g[a,b]*gr/(4*X)-g[a,b]*(1-L*X))
    P=L*r*r+r*s.diff(F,r)+F-1
    assert s.simplify(E(1,1)-2*r*s.diff(q,r))==0
    assert s.simplify(E(0,1)+z*P)==0
    assert s.simplify(E(0,0)+z*z*F*P-r*z*s.diff(F,x))==0


def test_state_has_regular_future_patch_components():
    rb,rc=horizons(100.)
    # r=rb+.01 near future BH horizon, outgoing EF regular at rc-.01
    assert np.isfinite(stress(rb+.01)['Trr'])
    assert np.isfinite(outgoing_stress(rc-.01)['Trr'])
    # sanity of leading finite asymptotic values without cancellation
    for d in [0.1,0.05]:
        v=outgoing_stress(rc-d)['Trr']
        assert abs(v)<1.


def test_independent_outgoing_ode():
    for r in [250.,300.,500.,650.,800.]:
        assert np.allclose(outgoing_radial_solution(r),
                           outgoing_radial_solution_ode(r),rtol=1e-9,atol=1e-12)


def test_outgoing_field_equation_sources():
    for r in [250.,300.,500.,650.,800.]:
        a,b,c,d=outgoing_field_residuals(r)
        assert np.isclose(a,b,rtol=1e-12,atol=1e-14)
        assert np.isclose(c,d,rtol=1e-12,atol=1e-14)


def test_scalar_F1_matching_two_charts():
    # gauge-invariant F=(grad r)^2, compared at identical physical events
    # to O(epsilon); rstar(400)=0 and matching homogeneous constants.
    for r in [240.,260.,300.,400.,500.,650.,800.,840.]:
        for v in [0.,100.,1000.]:
            ins,outs,u=coordinate_match(v,r)
            assert np.isclose(ins,outs,rtol=1e-8,atol=1e-11)


def test_tortoise_derivative():
    r=500.;h=1e-3
    dr=(tortoise(r+h)-tortoise(r-h))/(2*h)
    assert np.isclose(dr,1/f0(r),rtol=1e-9)


def test_sigma_integrability_identity():
    from two_patch import integrability_identity
    for r in [240.,300.,500.,700.,850.]:
        a,b=integrability_identity(r)
        assert np.isclose(a,b,rtol=2e-12,atol=1e-14)


def test_metric_overlap_coordinate_transform():
    from two_patch import coordinate_correction
    v=300.;r=500.;h=2e-3
    sigp,F1=radial_solution(r)
    from two_patch import outgoing_radial_solution
    sigm,_=outgoing_radial_solution(r)
    chi_v=sigp-sigm
    chi_r=(coordinate_correction(v,r+h)-coordinate_correction(v,r-h))/(2*h)
    pred=2*sigm/f0(r)+2*(F1+2*flux(100.,3e-6,1.,12)*v/r)/f0(r)**2
    assert np.isclose(chi_r,pred,rtol=1e-6,atol=1e-11)
    assert np.isclose((coordinate_correction(v+h,r)-coordinate_correction(v-h,r))/(2*h),chi_v,rtol=1e-8,atol=1e-11)
    # Numerical tensor transformation at fixed (v,r). No truncation of
    # rational Jacobian, verify O(eps^2) covariance against outgoing ansatz.
    sigp,F1r=radial_solution(r)
    F1=F1r+2*flux(100.,3e-6,1.,12)*v/r
    sigm=outgoing_radial_solution(r)[0]
    for eps in [0.25,0.125,0.0625]:
        u_v=1+eps*chi_v
        u_r=-2/f0(r)+eps*chi_r
        v_u=1/u_v
        v_r=-u_r/u_v
        gvv=-(1+eps*sigp)**2*(f0(r)+eps*F1)
        gvr=(1+eps*sigp)
        transformed_rr=gvv*v_r*v_r+2*gvr*v_r
        transformed_ur=gvv*v_u*v_r+gvr*v_u
        transformed_uu=gvv*v_u*v_u
        assert abs(transformed_rr)<1e-8
        assert abs(transformed_ur+(1+eps*sigm))<2e-8
        assert abs(transformed_uu+(1+eps*sigm)**2*(f0(r)+eps*F1))<1e-8
