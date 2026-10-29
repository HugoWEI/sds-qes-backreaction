"""Independent analytical, root, ultraviolet and dynamical-response tests."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from dataclasses import replace
from scipy.integrate import quad
from static_patch_candidate import Model
from radiation_algebra import OutgoingAlgebra
from peeling_response import (kappa_mass_derivative,alpha_slope,stationary_entropy,
    entropy_quadrature,entropy_independent_ode,analytic_response)


def test_chiral_rate_and_distinct_thermodynamic_proxy():
    m=Model();a=OutgoingAlgebra(m)
    assert np.isclose(a.rate(20000.),a.asymptotic_rate(),rtol=2e-9)
    assert np.isclose(a.rate(20000.),
        (a.entropy(20000.01)-a.entropy(19999.99))/.02,rtol=1e-9)
    assert a.asymptotic_rate()>a.thermodynamic_proxy_rate()


def test_two_null_closed_gradients_against_finite_differences():
    m=Model();r=350.;t=m.time_stationary(r);h=.001
    S=m.generalized_entropy
    num_t=(S(r,t+h)-S(r,t-h))/(2*h)
    num_r=(S(r+h,t)-S(r-h,t))/(2*h)
    assert abs(num_t-m.time_gradient(r,t))<1e-7
    assert abs(num_r-m.radial_gradient(r,t))<1e-6


def test_time_extremum_unique_in_spacelike_range():
    m=Model()
    for r in [220.,250.,350.,450.,490.]:
        D=m.effective_D(r);t=m.time_stationary(r)
        assert -D<t<D
        assert abs(m.time_gradient(r,t))<1e-9
        assert m.time_gradient(r,-.9*D)>0
        assert m.time_gradient(r,.9*D)<0


def test_uv_near_detector_root_not_controlled_island():
    m=Model();a,t=m.find_near_endpoint_root()
    gap=m.b-a;leading=m.c*m.G/(6*np.pi*m.b)
    assert np.isclose(gap,leading,rtol=3e-6)
    assert abs(m.time_gradient(a,t))<1e-7
    assert abs(m.radial_gradient(a,t))<1e-5
    assert m.uv_scale(a)/m.uv_cutoff<0.003
    assert max(np.linalg.eigvalsh(m.hessian(a,t)))<0


def test_no_macroscopic_exterior_root_in_reference_scan():
    m=Model();rb,_=m.horizons()
    gaps=np.geomspace(.01,m.b-rb-.00001,160)
    assert min(m.reduced_gradient(m.b-d) for d in gaps)>0


def test_tortoise_partial_fractions_against_independent_quadrature():
    m=Model()
    for r in [220.,300.,400.,450.]:
        direct=quad(lambda x:1/m.f(x),m.b,r,epsrel=1e-11)[0]
        assert np.isclose(m.rstar(r),direct,rtol=1e-10,atol=1e-9)


def test_kappa_M_derivative_finite_difference():
    m=Model();h=.0001
    fd=(replace(m,M=m.M+h).kappas()[0]-
        replace(m,M=m.M-h).kappas()[0])/(2*h)
    assert np.isclose(fd,kappa_mass_derivative(m),rtol=2e-8)
    assert alpha_slope(m)>0


def test_peeling_analytic_response_vs_quadrature_and_independent_ode():
    m=Model();T=20000.;h=1e-3
    s0=stationary_entropy(T,m)
    assert np.isclose(s0,entropy_quadrature(T,0,m),rtol=1e-12)
    fd=(entropy_quadrature(T,h,m)-
        entropy_quadrature(T,-h,m))/(2*h)
    s1=analytic_response(T,m)
    assert np.isclose(fd,s1,rtol=2e-8)
    assert np.isclose(entropy_quadrature(T,1,m),
                      entropy_independent_ode(T,1,m),atol=1e-10)
    assert np.isclose(entropy_quadrature(T,1,m)-s0,s1,rtol=3e-4)
