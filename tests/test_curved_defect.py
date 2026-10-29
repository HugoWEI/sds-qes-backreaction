"""Three-channel static-SdS Killing-flux Ward regression."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
import pytest
from baseline import kappa,flux,horizons
from continuum_collector import Pulse,analytic_energy_budget
from curved_defect import (Interface,cayley_scattering,rotation,
                            killing_curvature_scale,observer_acceleration)

def test_contact_cayley_and_unitarity():
    for angle in (-1.2,-.4,0.,.3,.7,1.3):
        S=cayley_scattering(angle)
        assert np.allclose(S,rotation(angle),atol=5e-16)
        assert np.allclose(S@S.conj().T,np.eye(2),atol=5e-16)

def test_full_3port_unitarity_cosmological_spectator():
    for angle in (-.7,0.,.7):
        S=np.eye(3);S[:2,:2]=rotation(angle)
        assert np.allclose(S.T@S,np.eye(3),atol=3e-16)
        assert np.allclose(S[2],(0,0,1))

def test_bath_does_not_double_count_cosmological_chiral_stream():
    z=Interface();JH,JC,JB=z.input_currents()
    assert np.isclose(JH-JC,z.integrated_ward()["easson_flux"],rtol=4e-15)
    assert np.isclose(JB,JC,atol=1e-19,rtol=1e-14)
    assert JH>JC>0

def test_local_killing_ward_including_work():
    w=Interface().local_ward()
    assert np.max(abs(w["local_residual"]))<2e-22
    assert np.any(w["switching_power"]>0)
    assert np.all(w["collector_gain"][w["angle"]==0]>=0)

def test_independent_integral_against_phase2():
    z=Interface(); d=z.integrated_ward()
    kb,kc=kappa(z.mass,z.L,z.G);scale=kb/np.pi
    ref=analytic_energy_budget(z.pulse,2,2*kb/kc,z.grid)
    assert np.isclose(d["heat"],z.N*scale*ref["heat_transfer"],rtol=3e-13)
    assert np.isclose(d["work"],z.N*scale*ref["external_work"],rtol=3e-13)
    assert abs(d["mismatch"])<1e-14
    assert np.isclose(d["bath_gain"],d["heat"]+d["work"]/2,atol=1e-14)
    assert np.isclose(d["patch_change"],-d["heat"]+d["work"]/2,atol=1e-14)

def test_proper_time_local_redshift_and_conservation():
    z=Interface();d=z.integrated_ward();w=z.local_proper_ward()
    assert np.isclose(d["window_proper"],d["redshift"]*d["window_killing"])
    assert np.isclose(d["work_local"],d["work"]/d["redshift"])
    assert np.isclose(np.max(w["proper_power"]),d["max_switch_power"]/d["redshift"]**2)
    assert np.max(abs(w["local_residual"]))<1e-21

def test_b_invariance_of_killing_flux_and_proper_growth():
    ds=[Interface(b=b).integrated_ward() for b in (250.,350.,500.,750.,850.)]
    for d in ds[1:]:
        for key in ("heat","work","bath_gain","patch_change"):
            assert np.isclose(d[key],ds[0][key],rtol=1e-14)
    assert len({round(d["work_local"],8) for d in ds})>1

def test_static_bath_observer_horizon_limit():
    rb,rc=horizons(100.)
    for rh in (rb,rc):
        far=rh+1. if rh==rb else rh-1.
        near=rh+.001 if rh==rb else rh-.001
        assert killing_curvature_scale(near)<killing_curvature_scale(far)
        assert observer_acceleration(near)>observer_acceleration(far)
        assert Interface(b=near).integrated_ward()["work_local"]>Interface(b=far).integrated_ward()["work_local"]

def test_work_scales_quadratically_under_small_angle():
    a=Interface(pulse=Pulse(angle=.1)).integrated_ward()["work"]
    b=Interface(pulse=Pulse(angle=.2)).integrated_ward()["work"]
    assert np.isclose(b/a,4.,rtol=4e-12)

def test_zero_switch_and_equal_bath_temperature():
    d=Interface(pulse=Pulse(angle=0.)).integrated_ward()
    assert abs(d["bath_gain"])<1e-14 and abs(d["work"])<1e-14
    h=Interface(beta_bath_ratio=1.).integrated_ward()
    assert abs(h["heat"])<1e-14 and h["work"]>0
    assert np.isclose(h["bath_gain"],h["work"]/2)

def test_two_species_scaling():
    one=Interface(N=1).integrated_ward()
    many=Interface(N=12).integrated_ward()
    for key in ("heat","work","bath_gain","patch_change"):
        assert np.isclose(many[key],12*one[key],rtol=2e-14)

def test_background_necessary_controls_not_full_junction():
    a=Interface().weak_backreaction_indicators()
    assert a["background_mass_drift_over_M"]<1e-3
    assert a["work_over_M"]<1e-4
    assert a["work_over_heat"]>0

def test_invalid_worldtubes_rejected():
    rb,rc=horizons(100.)
    for b in (rb,rc,100.,1000.):
        with pytest.raises(ValueError):Interface(b=b)
