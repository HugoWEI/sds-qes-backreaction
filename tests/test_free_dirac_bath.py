"""Free Dirac CFT and toy-bath checks; NO physical island claims."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from dataclasses import replace
from baseline import kappa
from free_dirac_multiinterval import (
    entropy_one_interval,entropy_two_intervals,mutual_information,
    same_static_time_points,first_order_entropy_response,
    lattice_mutual_information,continuum_mutual_information)
from transparent_bath_benchmark import SmoothBath
from bath_trial_extremum import analytic_gradient,trial_entropy,find_root,proper_gap

def test_dirac_MI_and_cutoff_invariance():
    p=same_static_time_points((300,350,450,500))
    S=entropy_two_intervals(p)
    MI=entropy_one_interval(*p[:2])+entropy_one_interval(*p[2:])-S
    assert np.isclose(MI,mutual_information(p),atol=1e-12)
    assert np.isclose(S-entropy_two_intervals(p,cutoffs=(2,2,2,2)),8*np.log(2),atol=1e-12)

def test_lattice_vs_known_continuum_cross_ratio():
    q=continuum_mutual_information(80,40,80)
    a=lattice_mutual_information(20,10,20)
    b=lattice_mutual_information(40,20,40)
    c=lattice_mutual_information(80,40,80)
    assert abs(q-c)<abs(q-b)<abs(q-a)
    assert abs(q-c)<3e-5

def test_dirac_geom_response_independent_finite_difference():
    p=same_static_time_points((300,350,450,500))
    ana=first_order_entropy_response(p)
    h=.1
    fd=(entropy_two_intervals(p,eps=h)-entropy_two_intervals(p,eps=-h))/(2*h)
    assert np.isclose(ana,fd,rtol=5e-8,atol=2e-11)

def test_bath_C2_matching_and_state():
    m=SmoothBath();rho,d1,d2=m.collar_derivatives()
    assert np.isclose(rho,m.bath_rho(0),atol=1e-14)
    h=.01
    fd1=(-3*m.bath_rho(0)+4*m.bath_rho(h)-m.bath_rho(2*h))/(2*h)
    fd2=(2*m.bath_rho(0)-5*m.bath_rho(h)+4*m.bath_rho(2*h)-m.bath_rho(3*h))/h**2
    assert np.isclose(d1,fd1,rtol=1e-5,atol=1e-8)
    assert np.isclose(d2,fd2,rtol=5e-3,atol=5e-10)
    assert np.allclose(m.endpoint('grav',77,m.b),m.endpoint('bath',77,0),atol=1e-12)

def _p():
    return [('grav',0,320),('grav',0,400),
            ('bath',0,20),('bath',0,100)]

def test_bath_MI_and_weyl_collars():
    m=SmoothBath();p=_p()
    MI=m.entropy_single(*p[:2])+m.entropy_single(*p[2:])-m.entropy_two(p)
    assert np.isclose(MI,m.MI(p),atol=1e-12)
    m2=replace(m,cutoff=.03,collar=10)
    assert np.isclose(m.MI(p),m2.MI(p),atol=1e-12)
    assert not np.isclose(m.entropy_two(p),replace(m,collar=10).entropy_two(p))

def test_bath_static_and_spatial_rate():
    m=SmoothBath();p=_p()
    shifted=[(side,t+1000,pos) for side,t,pos in p]
    assert np.isclose(m.entropy_two(p),m.entropy_two(shifted),atol=1e-10)
    kb,kc=kappa(m.mass,m.L,m.G)
    S=lambda t:m.entropy_single(('bath',t,0),('bath',t,t))
    assert np.isclose((S(40000.)-S(35000.))/5000.,m.N*(kb+kc)/12,rtol=1e-6)

def test_trial_two_null_gradient_numerically():
    m=SmoothBath();a=350.;t=2.;h=.01
    dr,dt=analytic_gradient(a,t,m)
    fr=(trial_entropy(a+h,t,m)-trial_entropy(a-h,t,m))/(2*h)
    ft=(trial_entropy(a,t+h,m)-trial_entropy(a,t-h,m))/(2*h)
    assert np.isclose(dr,fr,rtol=1e-7,atol=1e-6)
    assert np.isclose(dt,ft,rtol=1e-7,atol=1e-8)

def test_artificial_root_is_below_continuum_cutoff():
    m=SmoothBath();a,t,grad=find_root(m)
    assert np.linalg.norm(grad)<1e-4
    assert np.isclose(400-a,m.N*m.G/(6*np.pi*400),rtol=2e-5)
    assert proper_gap(a,m=m)<m.cutoff/300
    assert abs(t)<1e-6
