"""Independent NS-circle Gaussian entropy and purifier regression tests."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from continuum_collector import Pulse
from gaussian_boundary_state import (
    circle_thermal_covariance,transmitted_state,purification_cross,
    entropy,state_observables)


def test_ns_fourier_covariance_is_hermitian_and_physical():
    for beta in (2.,7.):
        C=circle_thermal_covariance(32,beta=beta)
        assert np.allclose(C,C.conj().T,atol=1e-14)
        eig=np.linalg.eigvalsh(C)
        assert np.all((eig>0)&(eig<1))
        assert np.isclose(np.trace(C),16.,atol=1e-12)


def test_exact_gaussian_scattering_matches_explicit_unitary():
    z=transmitted_state(n=32)
    CH,CB=z["initial_hot"],z["initial_bath"]
    c=np.diag(np.cos(z["angle"]));s=np.diag(np.sin(z["angle"]))
    U=np.block([[c,s],[-s,c]])
    old=np.block([[CH,np.zeros_like(CH)],[np.zeros_like(CB),CB]])
    new=np.block([[z["hot"],z["cross"]],
                  [z["cross"].conj().T,z["bath"]]])
    assert np.allclose(U@old@U.conj().T,new,atol=2e-14)
    assert np.allclose(U.conj().T@U,np.eye(64),atol=2e-14)


def test_covariance_positivity_and_entropy_unitarity():
    z=state_observables(n=32)
    assert 0<z["covariance_eigenvalue_bounds"][0]
    assert z["covariance_eigenvalue_bounds"][1]<1
    assert abs(z["global_entropy_invariant_residual"])<1e-10
    assert z["mutual_information_channels"]>0


def test_hot_purifier_projection_and_mutual_information():
    H=circle_thermal_covariance(32,beta=2.);D=purification_cross(H)
    P=np.block([[H,D],[D.conj().T,np.eye(32)-H]])
    assert np.allclose(P@P,P,atol=2e-13)
    assert entropy(P)<1e-10
    a=state_observables(n=32)
    assert 0<a["mutual_information_bath_hot_reference"]<=2*entropy(H)+1e-12


def test_no_contact_gives_zero_transfer_and_reference_information():
    z=state_observables(n=32,pulse=Pulse(angle=0.))
    assert abs(z["bath_entropy_change"])<1e-11
    assert abs(z["mutual_information_bath_hot_reference"])<1e-10
    assert abs(z["mutual_information_channels"])<1e-10


def test_equal_temperature_uniform_rotation_has_invariant_physical_marginals():
    # Uniform SO(2) mixing identical Gibbs fields leaves ALL H+B
    # observables unchanged, while changing MI with labeled R_H.
    H=circle_thermal_covariance(32,beta=2.)
    th=np.pi/4;c,s=np.cos(th),np.sin(th)
    Ho=(c*c+s*s)*H;Bo=(s*s+c*c)*H
    X=(-c*s+s*c)*H
    assert np.allclose(Ho,H,atol=1e-14)
    assert np.allclose(Bo,H,atol=1e-14)
    assert np.allclose(X,0.,atol=1e-14)
    D=purification_cross(H)
    BR=np.block([[Bo,-s*D],[-s*D.conj().T,np.eye(32)-H]])
    MI=entropy(Bo)+entropy(np.eye(32)-H)-entropy(BR)
    assert MI>1.


def test_finite_mode_entropy_convergence_under_smooth_switch():
    z=[state_observables(n=n) for n in (32,64,96)]
    for a in z:
        assert a["bath_entropy_change"]>1.9
        assert a["mutual_information_bath_hot_reference"]>2.6
    assert abs(z[1]["bath_entropy_change"]-z[2]["bath_entropy_change"])<1e-3
    assert abs(z[1]["mutual_information_bath_hot_reference"]-
               z[2]["mutual_information_bath_hot_reference"])<1e-3


def test_cross_channel_information_is_not_collected_entropy():
    z=state_observables(n=64)
    assert abs(z["mutual_information_bath_hot_reference"]-
               z["bath_entropy_change"])>.5
    assert not np.isclose(z["mutual_information_channels"],
                          z["mutual_information_bath_hot_reference"])
