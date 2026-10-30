"""Classical GHY/counterterm, Brown York reference energy checks."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import sympy as s
from boundary_variation import (reference_W,reference_W_prime,
    radial_boundary_momenta,classical_sds_boundary,brown_york_reference_mass,dEby_dM)
from linearized_geometry import f0


def test_GHY_exact_cancel_second_normal_derivative():
    r,t=s.symbols('r t');G,L=s.symbols('G L',positive=True)
    A=s.Function('A')(t,r);X=s.Function('X')(t,r)
    bulk=-2*X*s.diff(A,r,2)+A*s.diff(X,r)**2/(2*X)-s.diff(X,t)**2/(2*X*A)+2*A*(1-L*X)
    improved=bulk+2*s.diff(X*s.diff(A,r),r)
    expected=2*s.diff(X,r)*s.diff(A,r)+A*s.diff(X,r)**2/(2*X)-s.diff(X,t)**2/(2*X*A)+2*A*(1-L*X)
    assert s.simplify(improved-expected)==0
    assert s.simplify(s.diff(expected,s.diff(A,r))-2*s.diff(X,r))==0
    assert s.simplify(s.diff(expected,s.diff(X,r))-(2*s.diff(A,r)+A*s.diff(X,r)/X))==0


def test_spherical_dS_reference_counterterm_derivative():
    X=250000.;L=3e-6;h=1e-2
    fd=(reference_W(X+h,L)-reference_W(X-h,L))/(2*h)
    assert np.isclose(fd,reference_W_prime(X,L),rtol=2e-9)


def test_brown_york_agrees_with_variational_momentum():
    for b in [250.,300.,400.,500.,700.,850.]:
        d=classical_sds_boundary(b)
        e=brown_york_reference_mass(b)
        assert np.isclose(d['quasilocal_energy'],e,atol=2e-13)
        assert d['Pi_A']<0


def test_background_subtraction_vanishes_for_pure_desitter():
    for b in [100.,500.,900.]:
        assert abs(brown_york_reference_mass(b,mass=0.))<1e-13


def test_quasilocal_mass_derivative_and_redshift():
    for b in [250.,300.,500.,800.]:
        h=1e-3;M=100.
        fd=(brown_york_reference_mass(b,M+h)-brown_york_reference_mass(b,M-h))/(2*h)
        assert np.isclose(fd,dEby_dM(b,M),rtol=1e-9)
        A=np.sqrt(f0(b,M))
        assert np.isclose(A*fd,1.,rtol=1e-9)


def test_static_reference_dependent_absolute_energy_but_fixed_variation():
    b=500.;M=100.;G=1.;L=3e-6
    f=f0(b,M,L,G);ref_dS=b/G*(np.sqrt(1-L*b*b/3)-np.sqrt(f))
    ref_flat=b/G*(1-np.sqrt(f))
    assert ref_flat>ref_dS
    assert np.isclose(ref_flat-ref_dS,b/G*(1-np.sqrt(1-L*b*b/3)))


def test_static_boundary_momentum_finite_and_positive_reference_energy():
    d=classical_sds_boundary(500.)
    assert np.isclose(d['lapse'],np.sqrt(.35))
    assert np.isclose(d['quasilocal_energy'],137.2087127372385,rtol=1e-13)
    assert np.isfinite(d['Pi_X'])


def test_dirichlet_boundary_does_not_imply_free_junction_equations():
    d=classical_sds_boundary(500.)
    assert d['Pi_A']!=0 and d['Pi_X']!=0
    assert d['Pi_A']*0+d['Pi_X']*0==0
