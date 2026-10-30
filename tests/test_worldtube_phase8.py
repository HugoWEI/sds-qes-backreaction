"""Compact independent source/force/regression tests for Phase 8."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from keldysh_contact import Contact
from worldtube_phase8 import (StaticWorldtube,lapse_ctp,exact_lapse_source,
                              static_energy_variance,clock_covariance_check)


def test_local_clock_ctp_reparametrization():
    c=Contact(sites=3);p=np.array([0.,.2,.4,.2,0.]);dt=.17
    for b in (250.,500.,850.):
        z=clock_covariance_check(StaticWorldtube(b),c,p,dt)
        assert z['overlap_difference']<1e-12
        assert np.allclose(z['killing_noise'],StaticWorldtube(b).lapse**2*z['local_noise'])
        assert np.allclose(z['killing_response'],StaticWorldtube(b).lapse**2*z['local_response'])


def test_lapse_perturbation_source_matches_independent_finite_difference():
    c=Contact(sites=3);g=np.array([.1,.2,.5,.3,.1]);a=np.full(5,.61);dt=.2
    src=exact_lapse_source(c,a,g,dt)
    for j in (0,2,4):
        eps=1e-5;up=a.copy();dn=a.copy();up[j]+=eps;dn[j]-=eps
        fp=lapse_ctp(c,up,a,g,g,dt);fm=lapse_ctp(c,dn,a,g,g,dt)
        fd=(np.log(fp)-np.log(fm))/(2*eps)
        assert np.isclose(fd,-1j*dt*src[j],atol=1e-8)


def test_static_lapse_second_cumulant_gives_energy_variance():
    c=Contact(sites=3);E,var=static_energy_variance(c)
    a=np.array([.5,.6]);g=np.zeros(2);dt=.1;eps=1e-3
    def F(s):return lapse_ctp(c,a+s/(len(a)*dt),a,g,g,dt)
    first=(np.log(F(eps))-np.log(F(-eps)))/(2*eps)
    second=(np.log(F(eps))+np.log(F(-eps)))/(eps*eps)
    assert np.isclose(first,-1j*E,atol=1e-6)
    assert np.isclose(second,-var,atol=1e-6)
    assert var>0


def test_radial_support_and_distributional_ward():
    w=StaticWorldtube(250.)
    t,r=w.stress_divergence_delta_coefficients(6.,.02)
    assert np.isclose(t,.02/w.lapse)
    assert np.isclose(r,w.lapse**2*w.support_force(6.))
    assert w.support_force(6.)>0 and StaticWorldtube(850.).support_force(6.)<0


def test_force_vanishes_at_static_balance_radius():
    b=(3*100./3e-6)**(1/3)
    assert abs(StaticWorldtube(b).support_force(6.))<1e-14


def test_shared_branches_have_unit_ctp_normalization():
    c=Contact(sites=3);a=np.array([.5,.7,.9]);g=np.array([.1,.2,.3])
    assert abs(lapse_ctp(c,a,a,g,g,.12)-1)<1e-12
