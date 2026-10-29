import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from two_patch import tortoise
from interval_entropy import (interval_entropy, interval_entropy_response,
                              background_uv,linearized_rho)


def test_stationary_interval_entropy_not_radiation_page_curve():
    r1,r2=400.,500.;rs=tortoise(r2)
    values=[interval_entropy(t,r1,t+rs,r2) for t in [0.,100.,500.,1500.]]
    assert np.allclose(values,values[0],rtol=0,atol=1e-12)


def test_symmetry_endpoint_exchange():
    rs=tortoise(500.)
    assert np.isclose(interval_entropy(0.,400.,rs,500.),
                      interval_entropy(rs,500.,0.,400.))


def test_interval_linear_response_vs_finite_difference():
    rs=tortoise(500.)
    pred=interval_entropy_response(0.,400.,rs,500.)
    for h in [.5,.1,.02]:
        fp=interval_entropy(0.,400.,rs,500.,epsilon=h)
        fm=interval_entropy(0.,400.,rs,500.,epsilon=-h)
        assert np.isclose((fp-fm)/(2*h),pred,rtol=2e-6,atol=2e-9)


def test_conformal_factor_linearization():
    # In double-null coordinates: e^(2rho) = -2 exp(sigma_plus)/u_r.
    from linearized_geometry import f0,radial_solution
    from two_patch import outgoing_radial_solution
    r,v=500.,40.
    G=1.
    sigp,F1_rad=radial_solution(r)
    sigm,_=outgoing_radial_solution(r)
    from baseline import flux
    F1=F1_rad+2*G*flux(100.,3e-6,G,12)*v/r
    rho0,rho1=linearized_rho(v,r)
    assert np.isclose(2*rho1,F1/f0(r)+sigp+sigm)
    h=.1
    def rho_exact(e):
        ur=-2/f0(r)+e*(2*sigm/f0(r)+2*F1/f0(r)**2)
        return .5*np.log(-2*np.exp(e*sigp)/ur)
    diff=(rho_exact(h)-rho_exact(-h))/(2*h)
    assert np.isclose(diff,rho1,rtol=1e-6)
