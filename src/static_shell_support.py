"""Exact classical static SdS internal-shell consistency, no Polyakov stress.

A two-sided shell at fixed b joins two SdS vacua M_in and M_out with the
same cosmological constant, with continuous increasing-r normal.
X=r²; S_wall=-int d tau m(X). Distributional equations:
[n X]=-2G m,  -2[K]-[n X]/X=4G m_X.
They are NOT the boundary conditions between SdS gravity and external bath.
Time-dependent battery evolution is NOT solved by this static model.
"""
from dataclasses import dataclass
import numpy as np
from linearized_geometry import f0,f0_prime
from boundary_variation import brown_york_reference_mass


@dataclass(frozen=True)
class Shell:
    b:float=500.
    inner_mass:float=100.
    proper_energy:float=0.0011421612385233446
    G:float=1.
    L:float=3e-6

    def __post_init__(self):
        if self.G<=0 or self.proper_energy<0 or self.b<=0:
            raise ValueError('invalid shell parameters')
        if f0(self.b,self.inner_mass,self.L,self.G)<=0:
            raise ValueError('inner surface outside static patch')
        if self.G*self.proper_energy/self.b >= np.sqrt(f0(self.b,self.inner_mass,self.L,self.G)):
            raise ValueError('shell would cross outer static-patch horizon')

    def data(self):
        b,G,L=self.b,self.G,self.L
        m=self.proper_energy
        Ai=np.sqrt(f0(b,self.inner_mass,L,G))
        Ao=Ai-G*m/b
        dM=Ai*m-G*m*m/(2*b)
        Mo=self.inner_mass+dM
        Ko=f0_prime(b,Mo,L,G)/(2*Ao)
        Ki=f0_prime(b,self.inner_mass,L,G)/(2*Ai)
        dK=Ko-Ki
        dXn=2*b*(Ao-Ai)
        X=b*b
        dilaton_delta=-2*dK-dXn/X
        m_X_required=dilaton_delta/(4*G)
        Ai_prime=f0_prime(b,self.inner_mass,L,G)/(2*Ai)
        m_X_from_equilibrium=-(m*Ai_prime+G*m*m/(2*b*b))/(2*b*Ao)
        m_X_closed=m/(4*b*b)*(1-(1-L*b*b)/(Ao*Ai))
        Ei=brown_york_reference_mass(b,self.inner_mass,L,G)
        Eo=brown_york_reference_mass(b,Mo,L,G)
        return dict(inner_lapse=Ai,outer_lapse=Ao,inner_mass=self.inner_mass,
            outer_mass=Mo,mass_jump=dM,proper_energy=m,
            binding_correction=G*m*m/(2*b),
            jump_X_n=dXn,jump_K=dK,dilaton_delta=dilaton_delta,
            required_dilaton_charge=m_X_required,
            required_charge_from_radial_equilibrium=m_X_from_equilibrium,
            required_charge_closed=m_X_closed,
            energy_BY_inner=Ei,energy_BY_outer=Eo,
            static_dilaton_residual_if_no_support=dilaton_delta,
            metric_jump_residual=dXn+2*G*m,
            full_dilaton_residual=dilaton_delta-4*G*m_X_required)


def controller_shell_benchmark(interface=None,initial_multiplier=2.):
    """Quasistatic snapshots of a battery at b; NOT a dynamical solution."""
    if interface is None:
        from curved_defect import Interface
        interface=Interface()
    from defect_junction_constraints import ShellBattery
    battery=ShellBattery(interface,initial_multiplier).trajectory()
    initial=Shell(b=interface.b,inner_mass=interface.mass,
                  proper_energy=float(battery['epsilon_proper'][0]),
                  G=interface.G,L=interface.L).data()
    final=Shell(b=interface.b,inner_mass=interface.mass,
                proper_energy=float(battery['epsilon_proper'][-1]),
                G=interface.G,L=interface.L).data()
    return dict(initial=initial,final=final,
        needed_external_support_variation=final['required_dilaton_charge']-initial['required_dilaton_charge'],
        battery_killing_work=interface.integrated_ward()['work'])


def force_free_static_radius(inner_mass=100.,proper_energy=.001,G=1.,L=3e-6):
    """Two-vacuum classical dust shell's static equilibrium radius, if any."""
    from scipy.optimize import brentq
    from baseline import horizons
    rb,rc=horizons(inner_mass,L,G)
    def f(b):
        try:return Shell(b,inner_mass,proper_energy,G,L).data()['required_charge_closed']
        except ValueError:return np.nan
    lo=rb+.02*(rc-rb);hi=rc-.02*(rc-rb)
    if not(np.isfinite(f(lo)) and np.isfinite(f(hi)) and f(lo)*f(hi)<0):
        raise ValueError('no bracket for force-free shell within chosen static patch')
    root=brentq(f,lo,hi,xtol=1e-11)
    return root,(3*G*inner_mass/L)**(1/3)


def battery_gravitational_dressing(interface=None,initial_multiplier=2.):
    """Static-shell comparison only. E_K=A_in*m and Delta M=E_K-G E_K²/(2b A_in²).

    Finite work loss need not equal the static shell mass loss because the
    binding contribution changes. A real dynamic wall needs a support
    and gravitational energy ledger; it is NOT solved here.
    """
    from curved_defect import Interface
    from defect_junction_constraints import ShellBattery
    if interface is None:interface=Interface()
    traj=ShellBattery(interface,initial_multiplier).trajectory()
    Ek0=float(traj['energy_K'][0]);Ek1=float(traj['energy_K'][-1])
    A=np.sqrt(f0(interface.b,interface.mass,interface.L,interface.G))
    G=interface.G;b=interface.b
    binding=lambda E:G*E**2/(2*b*A*A)
    dM=lambda E:E-binding(E)
    return dict(initial_battery_K=Ek0,final_battery_K=Ek1,
      switching_work=Ek0-Ek1,
      mass_jump_initial=dM(Ek0),mass_jump_final=dM(Ek1),
      mass_loss_in_static_shell=dM(Ek0)-dM(Ek1),
      binding_energy_initial=binding(Ek0),
      binding_energy_final=binding(Ek1),
      gap_between_work_and_static_mass_loss=binding(Ek0)-binding(Ek1),
      instantaneous_mass_to_battery_response=1-G*Ek0/(b*A*A))


if __name__=='__main__':
    for e in (1e-3,.01,.1,1.0):
        d=Shell(proper_energy=e).data()
        print('proper shell energy',e,'mass shift',d['mass_jump'],
              'required m_X',d['required_dilaton_charge'],
              'no-support residual',d['static_dilaton_residual_if_no_support'])
    print('Battery snapshots:',controller_shell_benchmark())
