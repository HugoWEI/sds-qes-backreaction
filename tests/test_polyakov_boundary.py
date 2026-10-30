"""Independent symbolic and numerical tests of localized Polyakov GHY boundary."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import sympy as s
from baseline import horizons,kappa,flux
from linearized_geometry import f0,f0_prime,stress
from boundary_variation import classical_sds_boundary
from polyakov_boundary import (coefficients,state_chi,state_derivatives,local_stress,
    canonical_momenta,sds_boundary_data,horizon_regular_derivatives)


def test_symbolic_complete_GHY_integration_by_parts():
    t,l,G,N,L=s.symbols('t ell G N Lambda',positive=True)
    A=s.Function('A')(t,l);X=s.Function('X')(t,l);chi=s.Function('chi')(t,l)
    F=X/(4*G)-N*chi/(48*s.pi)
    bulk=-2*F*s.diff(A,l,2)
    improved=s.expand(bulk+2*s.diff(F*s.diff(A,l),l))
    expected=s.expand(2*s.diff(F,l)*s.diff(A,l))
    assert s.simplify(improved-expected)==0


def test_symbolic_Polyakov_canonical_momenta():
    t,l,G,N,L=s.symbols('t ell G N Lambda',positive=True)
    A=s.Function('A')(t,l);X=s.Function('X')(t,l);chi=s.Function('chi')(t,l)
    W=s.Function('W')(X)
    F=X/(4*G)-N*chi/(48*s.pi)
    Lag=(2*s.diff(F,l)*s.diff(A,l)
        + A*s.diff(X,l)**2/(8*G*X)-s.diff(X,t)**2/(8*G*X*A)
        + A*(1-L*X)/(2*G)
        -N*A*s.diff(chi,l)**2/(96*s.pi)
        +N*s.diff(chi,t)**2/(96*s.pi*A))
    pA=s.diff(Lag,s.diff(A,l))-W/(2*G)
    pX=s.diff(Lag,s.diff(X,l))-A*s.diff(W,X)/(2*G)
    pchi=s.diff(Lag,s.diff(chi,l))
    assert s.simplify(pA-(s.diff(X,l)-W)/(2*G)+N*s.diff(chi,l)/(24*s.pi))==0
    assert s.simplify(pX-s.diff(A,l)/(2*G)-A*s.diff(X,l)/(4*G*X)+A*s.diff(W,X)/(2*G))==0
    assert s.simplify(pchi+N*s.diff(A,l)/(24*s.pi)+N*A*s.diff(chi,l)/(48*s.pi))==0


def test_symbolic_chi_gauge_stress_and_eom():
    f,fp,fpp,kb,kc,N=s.symbols('f fp fpp kb kc N',positive=True)
    rho_u=-fp/4;rho_v=fp/4
    chi_u=fp/2-kb;chi_v=-fp/2-kc
    chi_uu=-f*fpp/4;chi_vv=-f*fpp/4
    Tuu=N/(48*s.pi)*(chi_u**2-2*chi_uu+4*rho_u*chi_u)
    Tvv=N/(48*s.pi)*(chi_v**2-2*chi_vv+4*rho_v*chi_v)
    expected_u=N*(2*f*fpp-fp**2+4*kb**2)/(192*s.pi)
    expected_v=N*(2*f*fpp-fp**2+4*kc**2)/(192*s.pi)
    assert s.simplify(Tuu-expected_u)==0
    assert s.simplify(Tvv-expected_v)==0
    chi_uv=f*fpp/4
    assert s.simplify(-4*chi_uv/f+fpp)==0
    Tuv=N*chi_uv/(24*s.pi)
    assert s.simplify(Tuv-N*f*fpp/(96*s.pi))==0
    assert s.simplify(-4*Tuv/f+N*fpp/(24*s.pi))==0


def test_local_auxiliary_stress_matches_prior_anomaly_stress():
    for r in [240.,300.,400.,500.,700.,850.]:
        a=local_stress(r);b=stress(r)
        for ch in ('Tuu','Tvv','Tuv'):
            assert np.isclose(a[ch],b[ch],rtol=2e-12,atol=1e-17)
        assert np.isclose(a['trace'],a['anomaly'],rtol=1e-12,atol=1e-17)
        assert np.isclose(a['net_chiral_flux'],flux(100.),rtol=1e-12)


def test_auxiliary_satisfies_full_static_polyakov_equation():
    r=500.;h=.01
    chp=state_chi(0,r+h);chm=state_chi(0,r-h);ch0=state_chi(0,r)
    chip=(chp-chm)/(2*h);chipp=(chp-2*ch0+chm)/h**2
    from linearized_geometry import f0_double_prime
    box=f0(r)*chipp+f0_prime(r)*chip
    assert np.isclose(box,-f0_double_prime(r),rtol=4e-7)
    assert np.isclose((state_chi(.1,r)-state_chi(-.1,r))/.2,
                      state_derivatives(r)['chi_t'],rtol=3e-11)


def test_state_horizon_regular_future_EF_radial_derivatives():
    rb,rc=horizons(100.)
    for which,rh,direction in [('bh',rb,1),('cosmo',rc,-1)]:
        limit=horizon_regular_derivatives(which)
        for d in (.05,.005):
            der=state_derivatives(rh+direction*d)
            value=der['future_bh_ingoing_r_slope' if which=='bh' else
                      'future_cosmo_outgoing_r_slope']
            assert np.isclose(value,limit,rtol=1e-3)


def test_radius_independent_auxiliary_boundary_momentum():
    kb,kc=kappa(100.)
    expected=-12*(kb-kc)/(48*np.pi)
    for r in [220.,250.,350.,400.,500.,700.,850.,875.]:
        a=sds_boundary_data(r)
        assert np.isclose(a['Pi_chi'],expected,rtol=2e-12,atol=2e-16)
        assert np.isclose(a['chi_ward_term'],a['state_killing_flux'],rtol=2e-12)
        assert np.isclose(a['Pi_A'],-a['bare_total_energy'],atol=1e-13)
        assert np.isclose(a['bare_total_energy']-a['classical_energy'],
                          a['bare_polyakov_energy_correction'],rtol=1e-9,atol=1e-12)


def test_natural_neumann_requires_boundary_state_source():
    z=sds_boundary_data(500.)
    assert z['Pi_chi']!=0
    assert np.isclose(z['required_coordinate_chi_source'],-z['Pi_chi'])
    assert np.isclose(z['required_proper_chi_source']*z['lapse'],
                      -z['Pi_chi'])
    assert np.isclose(z['required_coordinate_chi_source']*z['chi_t'],
                      -z['state_killing_flux'])


def test_conformal_state_onepoint_not_unique_for_auxiliary_boundary():
    kb,kc=kappa(100.)
    values=[]
    for su,sv in [(1,1),(1,-1),(-1,1),(-1,-1)]:
        hu=su*kb;hv=sv*kc
        values.append(12*(hu-hv)/(48*np.pi))
        assert np.isclose(hu**2,kb**2)
        assert np.isclose(hv**2,kc**2)
    assert len(set(np.round(values,14)))==4
    assert np.isclose(values[-1],sds_boundary_data(500.)['Pi_chi'])


def test_tangential_ward_kinematics_including_auxiliary():
    t=s.symbols('t',real=True)
    A=s.Function('A')(t);X=s.Function('X')(t);chi=s.Function('chi')(t)
    pa=s.Function('pa')(t);px=s.Function('px')(t);pc=s.Function('pc')(t)
    xi=s.Function('xi')(t)
    before=pa*s.diff(A*xi,t)+px*xi*s.diff(X,t)+pc*xi*s.diff(chi,t)
    after=xi*(-A*s.diff(pa,t)+px*s.diff(X,t)+pc*s.diff(chi,t))
    assert s.simplify(before-after-s.diff(pa*A*xi,t))==0


def test_momentum_finite_scheme_dependence_clarity():
    z=sds_boundary_data(500.)
    assert not np.isclose(z['bare_polyakov_energy_correction'],0,atol=1e-12)
    assert np.isclose(z['Pi_X'],classical_sds_boundary(500.)['Pi_X'])


def test_parameter_and_worldtube_failures():
    import pytest
    for r in (50.,1000.):
        with pytest.raises(ValueError):sds_boundary_data(r)
    with pytest.raises(ValueError):coefficients(N=0)
    with pytest.raises(ValueError):horizon_regular_derivatives('bad')
