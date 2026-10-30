"""Independent local boundary Wightman and point-split stress tests."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
import pytest
from transparent_state_kernel import (
    Pulse,thermal_kernel,thermal_regularizer,outgoing_kernel,
    point_split_currents,local_currents,offdiag_point_split,
    cross_coincident,smooth_pulse_derivative,theta_at,
    state_thermal_translation_violation)
from continuum_collector import thermal_current


def test_thermal_hadamard_and_hermiticity():
    for beta in (1.,2.,7.):
        for s in (.01,.1,1.):
            z=thermal_kernel(s,beta)
            assert np.isclose(z.conjugate(),thermal_kernel(-s,beta),atol=1e-14)
            assert np.isclose(thermal_regularizer(s,beta),thermal_regularizer(-s,beta))
            assert np.isclose(thermal_regularizer(s,beta),np.pi*s/beta/np.sinh(np.pi*s/beta))
    assert thermal_regularizer(0,2.)==1


def test_thermal_point_split_from_known_current():
    p=Pulse(angle=0.)
    for beta in (2.,7.):
        assert np.allclose(point_split_currents(7.,.001,beta,beta,p),
                           thermal_current(beta),rtol=2e-7)


def test_output_kernel_hermitian_matrix_reversing_times():
    for u,v in [(5.7,11.2),(10.,2.),(20.,7.)]:
        assert np.allclose(outgoing_kernel(u,v),
              outgoing_kernel(v,u).conj().T,atol=1e-14)


def test_output_point_split_stress_including_switch_work():
    for t in (4.8,5.,6.,8.,12.,18.,19.):
        prediction=local_currents(t)
        a=point_split_currents(t,.001)
        b=point_split_currents(t,.0005)
        assert np.allclose(a,prediction,rtol=3e-6,atol=1e-7)
        assert np.max(abs(b-prediction))<1.1*np.max(abs(a-prediction))+2e-8


def test_coincident_offdiagonal_state_coherence():
    for t in (5.,6.,18.,19.):
        expected=cross_coincident(t)
        a=offdiag_point_split(t,1e-3);b=offdiag_point_split(t,5e-4)
        assert abs(b-expected)<abs(a-expected)
        assert abs(b-expected)<2e-5
        assert np.isclose(expected.imag,smooth_pulse_derivative(t)/(2*np.pi),
                          rtol=1e-10)


def test_pulse_cross_channel_memory_at_separated_times():
    z=outgoing_kernel(9.,2.)
    assert theta_at(9.)>0 and theta_at(2.)==0
    assert abs(z[0,1])>1e-4
    assert abs(outgoing_kernel(9.,2.,pulse=Pulse(angle=0.))[0,1])<1e-14


def test_nonstationary_bath_wightman_and_KMS_failure():
    u,v=1.,2.
    assert abs(outgoing_kernel(u+8.,v+8.)[1,1]-outgoing_kernel(u,v)[1,1])>1e-3
    assert abs(state_thermal_translation_violation(u,v,8.))>1e-3
    assert abs(state_thermal_translation_violation(u,v,8.,pulse=Pulse(angle=0.)))<1e-13


def test_equal_temperature_rotation_still_has_switching_energy():
    t=5.
    actual=point_split_currents(t,1e-3,3.,3.)
    expected=thermal_current(3.)+smooth_pulse_derivative(t)**2/(4*np.pi)
    assert np.allclose(actual,[expected,expected],rtol=3e-6)
    assert expected>thermal_current(3.)


def test_two_point_kernel_requires_two_time_and_chirality_labels():
    z=outgoing_kernel(9.,2.)
    assert z.shape==(2,2)
    assert abs(z[0,1])>1e-4 and abs(z[1,0])>1e-6
    assert not np.allclose(z,np.eye(2)*np.trace(z)/2)


def test_input_validation():
    with pytest.raises(ValueError):thermal_regularizer(1.,0.)
    with pytest.raises(ValueError):thermal_kernel(0.,2.)
    with pytest.raises(ValueError):outgoing_kernel(1.,1.)
    with pytest.raises(ValueError):point_split_currents(5.,0.)
