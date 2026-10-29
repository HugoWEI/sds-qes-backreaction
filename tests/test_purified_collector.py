"""Independent canonical Gaussian purification checks."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from finite_collector import hopping,thermal_covariance,switching_profile,entropy
from purified_collector import (purified_channel,
    collector_hot_reference_information,collector_cold_reference_information,
    matching_temperature_ratio,sds_kappas)


def sample(n=32):
    h=hopping(n)
    return thermal_covariance(h,2.),thermal_covariance(h,7.)


def test_purification_is_projector():
    H,C=sample()
    for x in (H,C):
        P=purified_channel(x)
        np.testing.assert_allclose(P@P,P,rtol=0,atol=2e-13)
        assert abs(entropy(P))<1e-10


def test_no_pulse_no_hot_reference_mi():
    H,C=sample()
    assert abs(collector_hot_reference_information(H,C,np.zeros(32)))<1e-10


def test_full_swap_transfers_purification():
    H,C=sample();t=np.full(32,np.pi/2)
    I=collector_hot_reference_information(H,C,t)
    assert abs(I-2*entropy(H))<1e-10
    assert abs(collector_cold_reference_information(H,C,t))<1e-10


def test_finite_window_transfers_reference_information():
    H,C=sample()
    theta=switching_profile(32,5,27,5,np.pi/4)
    reg=np.where(theta>0)[0]
    I=collector_hot_reference_information(H,C,theta,reg)
    Ic=collector_cold_reference_information(H,C,theta,reg)
    assert 0<I<2*entropy(H)
    assert 0<Ic<2*entropy(C)


def test_sds_temperature_ratio():
    kb,kc=sds_kappas()
    assert np.isclose(matching_temperature_ratio()/2,kb/kc)


def test_equal_temperature_information_with_zero_work():
    # Uniform channel mixing commutes with the equal-temperature
    # free Hamiltonian. Information can transfer with zero heat current.
    from finite_collector import detailed_budget
    h=hopping(32);H=thermal_covariance(h,2.)
    theta=np.full(32,np.pi/4)
    r=detailed_budget(h,H,H,theta)
    assert abs(r['thermal_contrast_to_collector'])<1e-12
    assert abs(r['net_external_switch_work'])<1e-12
    assert abs(r['delta_energy_collector'])<1e-12
    assert collector_hot_reference_information(H,H,theta)>1


def test_reference_information_stable_with_lattice_size():
    from finite_collector import detailed_budget
    values=[]
    for n in (80,96):
        h=hopping(n)
        H=thermal_covariance(h,2.)
        C=thermal_covariance(h,matching_temperature_ratio())
        start=(n-40)//2
        t=switching_profile(n,start,start+40,10,np.pi/4)
        reg=np.where(t>0)[0]
        r=detailed_budget(h,H,C,t)
        mi=collector_hot_reference_information(H,C,t,reg)
        values.append((r['net_external_switch_work'],
                       r['delta_energy_collector'],mi))
    np.testing.assert_allclose(values[0],values[1],rtol=0,atol=1e-9)
