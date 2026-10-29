"""Test energy conservation, exact unitary and finite-window work."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from finite_collector import (hopping,thermal_covariance,channel_covariances,
    switching_profile,full_unitary,entropy,energy,work_formula,
    detailed_budget)


def sample(n=48):
    h=hopping(n)
    return h,thermal_covariance(h,2.),thermal_covariance(h,8.)


def test_exact_unitarity():
    t=switching_profile(48,8,40,8,.71)
    U=full_unitary(t)
    np.testing.assert_allclose(U@U.T,np.eye(96),rtol=0,atol=2e-15)


def test_full_covariance_vs_explicit_unitary():
    h,a,b=sample();t=switching_profile(48,8,40,8)
    ap,bp,x=channel_covariances(a,b,t)
    U=full_unitary(t)
    old=np.block([[a,np.zeros_like(a)],[np.zeros_like(b),b]])
    np.testing.assert_allclose(U@old@U.T,np.block([[ap,x],[x.T,bp]]),atol=3e-15)


def test_number_and_global_entropy_invariants():
    h,a,b=sample();t=switching_profile(48,8,40,8)
    ap,bp,x=channel_covariances(a,b,t)
    C=np.block([[ap,x],[x.T,bp]])
    assert abs(np.trace(C)-np.trace(a)-np.trace(b))<1e-12
    assert abs(entropy(C)-entropy(a)-entropy(b))<1e-10


def test_exact_energy_and_work_decomposition():
    h,a,b=sample();t=switching_profile(48,8,40,8)
    r=detailed_budget(h,a,b,t)
    assert abs(r['delta_energy_hot']+r['delta_energy_collector']-
               r['net_external_switch_work'])<2e-12
    assert abs(r['delta_energy_collector']-r['thermal_contrast_to_collector']-
               r['collector_switch_work'])<2e-12
    assert abs(r['net_external_switch_work']-r['collector_switch_work']-
               r['hot_switch_work'])<2e-12
    assert r['net_external_switch_work']>0


def test_zero_switch_is_identity():
    h,a,b=sample()
    ap,bp,x=channel_covariances(a,b,np.zeros(48))
    np.testing.assert_allclose(ap,a,atol=1e-15)
    np.testing.assert_allclose(bp,b,atol=1e-15)
    assert np.max(abs(x))==0


def test_uniform_switch_has_no_work():
    h,a,b=sample()
    r=detailed_budget(h,a,b,np.full(48,np.pi/4))
    assert abs(r['net_external_switch_work'])<1e-12
    contrast=.5*(energy(h,a)-energy(h,b))
    assert abs(r['delta_energy_collector']-contrast)<1e-12
    assert abs(r['delta_energy_hot']+contrast)<1e-12


def test_temperature_gradient_transfers_energy():
    h,a,b=sample();r=detailed_budget(h,a,b,switching_profile(48,8,40,8))
    assert r['thermal_contrast_to_collector']>0
    assert r['delta_energy_collector']>0
    assert r['total_channel_mutual_information']>0


def test_support_and_work():
    h,a,b=sample();t=switching_profile(48,8,40,8)
    assert t[:9].sum()==0 and t[40:].sum()==0
    assert np.isclose(t.max(),np.pi/4)
    assert work_formula(h,a,b,t)>0


def test_small_angle_work_is_quadratic():
    h,a,b=sample();t=switching_profile(48,8,40,8,peak_angle=.1)
    r=work_formula(h,a,b,2*t)/work_formula(h,a,b,t)
    assert 3.95<r<4.01


def test_smoother_window_reduces_work():
    h,a,b=sample()
    W4=work_formula(h,a,b,switching_profile(48,8,40,4))
    W8=work_formula(h,a,b,switching_profile(48,8,40,8))
    assert W8<W4
