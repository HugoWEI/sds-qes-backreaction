"""Curved SdS worldtube: exact 3-channel chiral contact energy ledger.

Two unfolded right movers H (outgoing black hole) and B (independent bath)
mix; incoming cosmological C is spectator. The background conformal metric
is ds^2=f(r)(-dt^2+drstar^2). This is a prescribed scattering action,
NOT a solved gravitating-to-nongravitating junction or QES.
"""
from dataclasses import dataclass
import numpy as np
from baseline import kappa,flux,horizons
from linearized_geometry import f0,f0_prime
from continuum_collector import Pulse,spectral_derivative

SIGMA_Y=np.array([[0.,-1.j],[1.j,0.]])

def cayley_scattering(theta):
    """Symmetrically regulated contact V=-g sigma_y, theta=2atan(g/2)."""
    g=2*np.tan(float(theta)/2)
    V=-g*SIGMA_Y;I=np.eye(2,dtype=complex)
    return np.linalg.solve(I+.5j*V,I-.5j*V)

def rotation(theta):
    c,s=np.cos(theta),np.sin(theta)
    return np.array([[c,s],[-s,c]])

def killing_curvature_scale(r,mass=100.,L=3e-6,G=1.):
    f=f0(r,mass,L,G)
    if f<=0:raise ValueError("outside static patch")
    return np.sqrt(f)

def observer_acceleration(r,mass=100.,L=3e-6,G=1.):
    return abs(f0_prime(r,mass,L,G))/(2*killing_curvature_scale(r,mass,L,G))

@dataclass(frozen=True)
class Interface:
    mass:float=100.
    L:float=3e-6
    G:float=1.
    N:int=12
    b:float=500.
    beta_bath_ratio:float|None=None
    pulse:Pulse=Pulse()
    grid:int=4096

    def __post_init__(self):
        rb,rc=horizons(self.mass,self.L,self.G)
        if not rb<self.b<rc:raise ValueError("b must lie between horizons")
        if self.N<=0 or self.grid<128:raise ValueError("invalid N/grid")

    def input_currents(self):
        kb,kc=kappa(self.mass,self.L,self.G)
        JH=self.N*kb**2/(48*np.pi);JC=self.N*kc**2/(48*np.pi)
        ratio=kb/kc if self.beta_bath_ratio is None else self.beta_bath_ratio
        if ratio<=0:raise ValueError("nonpositive bath beta ratio")
        return JH,JC,JH/ratio**2

    def time_data(self):
        kb,_=kappa(self.mass,self.L,self.G);scale=kb/np.pi
        x,theta=self.pulse.grid(self.grid)
        td=spectral_derivative(theta,self.pulse.period)
        return x,theta,td,scale

    def local_ward(self):
        x,theta,td,scale=self.time_data()
        JH,JC,JB=self.input_currents()
        trans=np.sin(theta)**2
        P=self.N*scale**2*td**2/(2*np.pi)
        Q=trans*(JH-JB)
        before=np.full(self.grid,JH-JC)
        after=before-Q+P/2
        bath_gain=Q+P/2
        return dict(x=x,angle=theta,transmission=trans,
            switching_power=P,heat_flux=Q,patch_before=before,
            patch_after=after,collector_gain=bath_gain,
            local_residual=(after-before)+bath_gain-P,
            dK_dt=self.pulse.period/(self.grid*scale),
            scale=scale,J_hot=JH,J_cosmological=JC,J_bath=JB)

    def integrated_ward(self):
        w=self.local_ward();dt=w["dK_dt"]
        Q=float(np.sum(w["heat_flux"])*dt)
        W=float(np.sum(w["switching_power"])*dt)
        patch=float(np.sum(w["patch_after"]-w["patch_before"])*dt)
        bath=float(np.sum(w["collector_gain"])*dt)
        lapse=killing_curvature_scale(self.b,self.mass,self.L,self.G)
        return dict(heat=Q,work=W,patch_change=patch,bath_gain=bath,
            mismatch=patch+bath-W,window_killing=self.pulse.period/w["scale"],
            window_proper=lapse*self.pulse.period/w["scale"],
            b=self.b,redshift=lapse,
            input_net=w["J_hot"]-w["J_cosmological"],
            easson_flux=flux(self.mass,self.L,self.G,self.N),
            max_switch_power=float(np.max(w["switching_power"])),
            max_local_power=float(np.max(w["switching_power"]))/lapse**2,
            work_local=W/lapse,
            local_acceleration=observer_acceleration(self.b,self.mass,self.L,self.G))

    def local_proper_ward(self):
        w=self.local_ward()
        f=f0(self.b,self.mass,self.L,self.G)
        return dict(proper_power=w["switching_power"]/f,
                    proper_heat_flux=w["heat_flux"]/f,
                    proper_bath_gain=w["collector_gain"]/f,
                    proper_patch_change=(w["patch_after"]-w["patch_before"])/f,
                    local_residual=w["local_residual"]/f)

    def weak_backreaction_indicators(self):
        d=self.integrated_ward()
        return dict(work_over_M=d["work"]/self.mass,
            collected_energy_over_M=d["bath_gain"]/self.mass,
            background_mass_drift_over_M=d["easson_flux"]*d["window_killing"]/self.mass,
            proper_acceleration=d["local_acceleration"],
            work_over_heat=d["work"]/d["heat"] if d["heat"]>0 else np.inf)

if __name__=="__main__":
    for b in (250.,300.,500.,700.,850.):
        m=Interface(b=b)
        print("b=",b,"ledger=",m.integrated_ward())
