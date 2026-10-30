"""Classical two-sided fixed-radius SdS shell: energy AND force condition."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import sympy as s
from scipy.optimize import brentq
from static_shell_support import Shell,force_free_static_radius,controller_shell_benchmark
from linearized_geometry import f0
from baseline import horizons


def test_exact_shell_mass_jump_from_sqrt_f():
    for e in [1e-3,.05,.1,1.]:
        a=Shell(proper_energy=e).data()
        assert np.isclose(a['mass_jump'],a['inner_lapse']*e-e**2/(2*500),rtol=1e-12)
        assert np.isclose(a['outer_lapse']**2,f0(500,a['outer_mass']),atol=3e-15)
        assert abs(a['metric_jump_residual'])<1e-12
        assert a['mass_jump']>0


def test_symbolic_dilaton_charge_equals_radial_force_equilibrium():
    b,G,L,M,m=s.symbols('b G L M m',positive=True)
    f=1-2*G*M/b-L*b*b/3
    Ai=s.sqrt(f);Ao=Ai-G*m/b
    dM=Ai*m-G*m*m/(2*b)
    mx=s.symbols('mx')
    radial=s.diff(dM,b)+s.diff(dM,m)*2*b*mx
    solve=s.solve(radial,mx)[0]
    closed=m/(4*b*b)*(1-(1-L*b*b)/(Ao*Ai))
    assert s.simplify(s.together(solve-closed))==0


def test_numeric_dilaton_equation_and_independent_force():
    for b in [300.,400.,500.,700.,850.]:
        d=Shell(b=b,proper_energy=.01).data()
        assert np.isclose(d['required_dilaton_charge'],
                          d['required_charge_from_radial_equilibrium'],rtol=8e-7,atol=2e-15)
        assert np.isclose(d['required_dilaton_charge'],
                          d['required_charge_closed'],rtol=8e-7,atol=2e-15)
        assert abs(d['full_dilaton_residual'])<1e-20


def test_no_dilaton_charge_generically_does_not_hold_shell_static():
    d=Shell(b=500.,proper_energy=.1).data()
    assert abs(d['static_dilaton_residual_if_no_support'])>1e-8
    assert abs(d['metric_jump_residual'])<1e-12


def test_static_force_balance_radius_has_shift():
    for e in [.001,.1,1.]:
        r,r0=force_free_static_radius(proper_energy=e)
        assert r>r0
        d=Shell(b=r,proper_energy=e).data()
        assert abs(d['required_charge_closed'])<1e-14
        A=np.sqrt(f0(r0,100.,3e-6,1.))
        pred=e*A/(2*3e-6*r0*r0)
        assert np.isclose(r-r0,pred,rtol=.002)


def test_charge_changes_sign_on_sides_of_force_balance():
    r0=(3*100./3e-6)**(1/3)
    a=Shell(b=r0-10,proper_energy=.001).data()['required_charge_closed']
    b=Shell(b=r0+10,proper_energy=.001).data()['required_charge_closed']
    assert a<0<b


def test_controller_battery_needs_adjusted_support_force():
    a=controller_shell_benchmark()
    assert a['initial']['proper_energy']>a['final']['proper_energy']
    assert a['initial']['required_dilaton_charge']>a['final']['required_dilaton_charge']>0
    assert a['needed_external_support_variation']<0
    assert np.isclose(a['initial']['mass_jump']-a['final']['mass_jump'],
                      a['battery_killing_work'],rtol=3e-6)


def test_nonrelativistic_shell_mass_binding_term():
    b=500.;e=.1
    a=Shell(b=b,proper_energy=e).data()
    assert np.isclose(a['mass_jump']+a['binding_correction'],a['inner_lapse']*e,atol=1e-14)
    assert a['binding_correction']>0


def test_mass_difference_matches_brown_york_to_first_order():
    m=1e-4;d=Shell(proper_energy=m).data()
    assert np.isclose(d['energy_BY_outer']-d['energy_BY_inner'],m,rtol=5e-7)


def test_battery_gravitational_binding_gap_exact():
    from static_shell_support import battery_gravitational_dressing
    z=battery_gravitational_dressing()
    assert z['gap_between_work_and_static_mass_loss']>0
    assert np.isclose(z['switching_work']-z['mass_loss_in_static_shell'],
        z['gap_between_work_and_static_mass_loss'],rtol=1e-10,atol=2e-15)
    assert np.isclose(z['mass_jump_initial'],
        z['initial_battery_K']-z['binding_energy_initial'],rtol=1e-14,atol=1e-14)
    assert 0<z['instantaneous_mass_to_battery_response']<1


def test_static_dressed_mass_response_finite_difference():
    from static_shell_support import battery_gravitational_dressing
    z=battery_gravitational_dressing()
    e=z['initial_battery_K'];h=1e-7
    b=500.;G=1.;A=np.sqrt(f0(b,100.,3e-6,G))
    fun=lambda E:E-G*E**2/(2*b*A*A)
    fd=(fun(e+h)-fun(e-h))/(2*h)
    assert np.isclose(fd,z['instantaneous_mass_to_battery_response'],rtol=1e-12,atol=1e-12)


def test_tangential_boundary_reparametrization_ward_kinematics():
    t=s.symbols('t',real=True)
    A=s.Function('A')(t);X=s.Function('X')(t)
    pA=s.Function('Pi_A')(t);pX=s.Function('Pi_X')(t)
    xi=s.Function('xi')(t)
    varied=pA*s.diff(A*xi,t)+pX*xi*s.diff(X,t)
    integrated=-xi*A*s.diff(pA,t)+xi*pX*s.diff(X,t)
    assert s.simplify(varied-integrated-s.diff(xi*A*pA,t))==0
