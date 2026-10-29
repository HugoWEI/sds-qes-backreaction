"""Causal free-Dirac continuum stress and spectral regularization regressions."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from continuum_collector import (
    Pulse,smooth_step,thermal_current,spectral_derivative,
    analytic_energy_budget,circle_current,finite_circle_prediction,
    spectral_energy_budget,sds_collector_scale)
from baseline import kappa


def test_compact_smooth_switch_and_local_unitarity():
    p=Pulse()
    u,t=p.grid(2048)
    assert t[0]==0 and np.all(t[u<p.start]==0)
    assert np.all(t[u>p.stop]==0)
    assert np.all((t>=0)&(t<=p.angle))
    assert smooth_step(np.array([-1.,0.,1.,2.])).tolist()==[0.,0.,1.,1.]
    for z in [0.01,.1,.2,.5,.8,.99]:
        assert np.isclose(smooth_step(z)+smooth_step(1-z),1.,atol=2e-14)
    assert np.allclose(np.cos(t)**2+np.sin(t)**2,1.,atol=1e-15)


def test_continuum_conservation_work_positive_and_nontrivial():
    a=analytic_energy_budget()
    assert a["zero_outside"]
    assert a["heat_transfer"]>0 and a["external_work"]>0
    assert a["hot_delta"]<0 and a["collector_delta"]>0
    assert abs(a["hot_delta"]+a["collector_delta"]-a["external_work"])<1e-14
    assert np.isclose(a["external_work"],.08517469464266,atol=4e-12)
    assert np.isclose(a["heat_transfer"],.31218142867511,atol=4e-12)


def test_spectral_derivative_matches_finite_difference():
    u,t=Pulse().grid(4096)
    derivative=spectral_derivative(t,24)
    central=(np.roll(t,-1)-np.roll(t,1))/(2*(u[1]-u[0]))
    assert np.max(abs(derivative-central))<4e-5
    assert abs(derivative[0])<1e-11


def test_fourier_sea_properly_converges_to_circle_prediction():
    p=Pulse()
    ref=finite_circle_prediction(p)
    errors=[]
    for cutoff in (20,40,80,160):
        a=spectral_energy_budget(p,n_modes=cutoff)
        errors.append(max(abs(a["hot_delta"]-ref["hot_delta"]),
                          abs(a["collector_delta"]-ref["collector_delta"])))
        assert abs(a["coefficient_unitarity"]-1)<1e-12
    assert errors[-1]<1e-9
    assert errors[0]>errors[1]>errors[2]>errors[3]
    assert abs(a["fermion_number_delta"])<1e-10
    assert np.isclose(a["external_work"],ref["external_work"],atol=1e-9)


def test_thermal_circle_is_not_infinite_line_at_finite_period():
    assert abs(circle_current(7,24)-thermal_current(7))>1e-6
    assert np.isclose(circle_current(2,24)-thermal_current(2),
                      np.pi/(12*24**2),atol=2e-9)
    a=analytic_energy_budget()
    b=finite_circle_prediction()
    assert abs(a["external_work"]-b["external_work"])<1e-15
    assert abs(a["heat_transfer"]-b["heat_transfer"])>1e-6


def test_zero_pulse_and_equal_temperature_limits():
    zero=analytic_energy_budget(Pulse(angle=0.))
    assert abs(zero["external_work"])<1e-12
    assert abs(zero["heat_transfer"])<1e-12
    same=analytic_energy_budget(beta_hot=3.,beta_cold=3.)
    assert same["heat_transfer"]==0 and same["external_work"]>0
    assert np.isclose(same["hot_delta"],same["collector_delta"])
    beta=analytic_energy_budget(Pulse(angle=0.),beta_hot=2.,beta_cold=7.)
    assert abs(beta["collector_delta"])<1e-12


def test_dilation_tradeoff_work_vs_heat():
    p=Pulse()
    q=Pulse(period=48.,start=8.,stop=40.,ramp=6.,angle=.7)
    a=analytic_energy_budget(p)
    b=analytic_energy_budget(q)
    assert np.isclose(b["heat_transfer"],2*a["heat_transfer"],rtol=1e-13)
    assert np.isclose(b["external_work"],a["external_work"]/2,rtol=1e-13)


def test_cauchy_schwarz_lower_bound_for_switch_work():
    p=Pulse()
    bound=p.angle**2/(2*np.pi)*(1/p.ramp+1/p.ramp)
    assert analytic_energy_budget(p)["external_work"]>bound


def test_sds_temperature_ratio_and_killing_units():
    d=sds_collector_scale()
    kb,kc=kappa(100.,3e-6,1.)
    assert np.isclose(d["thermal_ratio"],kb/kc)
    assert np.isclose(d["sdS_net_flux"],d["physical_flux_difference"],rtol=1e-14)
    assert np.isclose(d["physical_collector_energy_gain"],
                      d["physical_collected_heat"]+d["physical_switching_work"]/2,
                      rtol=1e-14)
    assert np.isclose(d["physical_switching_work"]/d["physical_collected_heat"],
                      .28807230703664,rtol=1e-11)
