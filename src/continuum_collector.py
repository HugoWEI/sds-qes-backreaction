"""Causal chiral free-Dirac collector: continuum point-splitting versus modes.

Two unfolded co-propagating complex chiral fields:
Psi_out(u)=R(theta(u))*Psi_in(u).
This is a postulated local scattering boundary condition, NOT an SdS
gravity-bath action, globally matched state, or island calculation.

Per chiral complex Dirac species, J_beta=pi/(12 beta^2).
At a smooth SO(2) interface:
JHout=c^2 JH+s^2 JC+theta_dot^2/(4 pi),
JCout=s^2 JH+c^2 JC+theta_dot^2/(4 pi).
External work density =theta_dot^2/(2 pi).
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import expit
from baseline import kappa, flux


def smooth_step(z):
    z=np.asarray(z,float)
    y=np.zeros_like(z)
    active=(z>0)&(z<1)
    zz=z[active]
    y[active]=expit(1/(1-zz)-1/zz)
    y[z>=1]=1.
    return y


@dataclass(frozen=True)
class Pulse:
    period:float=24.
    start:float=4.
    stop:float=20.
    ramp:float=3.
    angle:float=0.7

    def __post_init__(self):
        if not(self.period>0 and self.ramp>0 and
               0<=self.start<self.start+self.ramp<=self.stop-self.ramp<self.stop<=self.period):
            raise ValueError("Invalid smooth compact pulse")

    def values(self,u):
        u=np.asarray(u,float)
        return self.angle*smooth_step((u-self.start)/self.ramp)*smooth_step((self.stop-u)/self.ramp)

    def grid(self,n):
        u=np.arange(n)*self.period/n
        return u,self.values(u)


def thermal_current(beta):
    return np.pi/(12*beta**2)


def spectral_derivative(values,period):
    n=len(values)
    momenta=2*np.pi*np.fft.fftfreq(n,d=period/n)
    return np.fft.ifft(1j*momenta*np.fft.fft(values)).real


def analytic_energy_budget(pulse=Pulse(),beta_hot=2.,beta_cold=7.,n_grid=4096):
    if beta_hot<=0 or beta_cold<=0:raise ValueError("beta must be positive")
    u,t=pulse.grid(n_grid)
    dot=spectral_derivative(t,pulse.period)
    du=pulse.period/n_grid
    trans=np.sin(t)**2
    jh=thermal_current(beta_hot);jc=thermal_current(beta_cold)
    q=float(trans.sum()*du*(jh-jc))
    work=float(np.sum(dot**2)*du/(2*np.pi))
    return dict(heat_transfer=q,external_work=work,
                hot_delta=-q+work/2,collector_delta=q+work/2,
                zero_outside=bool(np.all(t[(u<pulse.start)|(u>pulse.stop)]==0)),
                max_transmission=float(trans.max()),
                energy_balance_residual=float(-q+work/2+q+work/2-work),
                thermal_hot_current=jh,thermal_cold_current=jc)


def circle_current(beta,period,n_max=350):
    """Antiperiodic circle thermal energy density relative to filled NS sea."""
    n=np.arange(-n_max,n_max+1)
    w=2*np.pi*(n+.5)/period
    return float(np.sum(w*(expit(-beta*w)-(w<0)))/period)


def spectral_energy_budget(pulse=Pulse(),beta_hot=2.,beta_cold=7.,
                           n_grid=4096,n_modes=80,n_sidebands=120):
    """Independent Fourier scattering with untruncated incoming Dirac sea.

    For outgoing mode m sum incoming modes n=m-q. Never sharply truncate both
    incoming and outgoing seas to the same finite band: that creates fake
    negative-energy-edge work diverging with the mode cutoff.
    """
    if n_sidebands>=n_grid//2:raise ValueError("Nyquist")
    _,theta=pulse.grid(n_grid)
    c=np.fft.fft(np.cos(theta))/n_grid
    s=np.fft.fft(np.sin(theta))/n_grid
    q=np.arange(-n_sidebands,n_sidebands+1)
    wc=np.abs(c[q%n_grid])**2
    ws=np.abs(s[q%n_grid])**2
    m=np.arange(-n_modes,n_modes+1)[:,None]
    w=2*np.pi*(m+.5)/pulse.period
    shifted=2*np.pi*(m-q[None,:]+.5)/pulse.period
    fh=expit(-beta_hot*w);fc=expit(-beta_cold*w)
    dh=(wc*(expit(-beta_hot*shifted)-fh)+
        ws*(expit(-beta_cold*shifted)-fh)).sum(axis=1)
    dc=(wc*(expit(-beta_cold*shifted)-fc)+
        ws*(expit(-beta_hot*shifted)-fc)).sum(axis=1)
    eh=float(np.sum(w[:,0]*dh));ec=float(np.sum(w[:,0]*dc))
    return dict(hot_delta=eh,collector_delta=ec,external_work=eh+ec,
                fermion_number_delta=float(np.sum(dh+dc)),
                coefficient_unitarity=float(np.sum(wc+ws)),
                mode_limit=n_modes,sideband_limit=n_sidebands)


def finite_circle_prediction(pulse=Pulse(),beta_hot=2.,beta_cold=7.,n_grid=4096):
    """Circle thermal currents plus continuum local switching work."""
    a=analytic_energy_budget(pulse,beta_hot,beta_cold,n_grid)
    theta=pulse.grid(n_grid)[1]
    Q=(circle_current(beta_hot,pulse.period)-
       circle_current(beta_cold,pulse.period))*np.sum(np.sin(theta)**2)*pulse.period/n_grid
    W=a["external_work"]
    return dict(heat_transfer=float(Q),hot_delta=float(-Q+W/2),
                collector_delta=float(Q+W/2),external_work=W)


def sds_collector_scale(m=100.,L=3e-6,G=1.,N=12,pulse=Pulse(),n_grid=4096):
    """Match SdS *thermal ratio and Killing current*, not physical junction."""
    kb,kc=kappa(m,L,G)
    scale=kb/np.pi
    bc=2*kb/kc
    d=analytic_energy_budget(pulse,2.,bc,n_grid)
    return dict(surface_gravity_black_hole=kb,
                surface_gravity_cosmological=kc,
                physical_window_duration=pulse.period/scale,
                thermal_ratio=kb/kc,sdS_net_flux=flux(m,L,G,N),
                physical_flux_difference=N*scale**2*
                    (thermal_current(2)-thermal_current(bc)),
                physical_collected_heat=N*scale*d["heat_transfer"],
                physical_switching_work=N*scale*d["external_work"],
                physical_collector_energy_gain=N*scale*d["collector_delta"])


if __name__=="__main__":
    p=Pulse()
    print("Infinite line:",analytic_energy_budget(p))
    print("Antiperiodic circle:",finite_circle_prediction(p))
    for cutoff in (10,20,40,80,160):
        print("modes",cutoff,spectral_energy_budget(p,n_modes=cutoff))
    print("SdS scale",sds_collector_scale(pulse=p))
