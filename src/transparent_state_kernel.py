"""Exact 2-channel free-Dirac Wightman BOUNDARY KERNEL at a switched contact.

G_beta(u-v)=i/(2 beta sinh(pi*(u-v)/beta));
G_out(u,v)=R(theta(u)) diag(G_H(u-v),G_B(u-v)) R(theta(v))^T.

Frozen unfolded chiral model, NOT global curved-space retarded bath Green
function, black-hole purifier, or QES. The two-time correlator determines
the specified Gaussian interface state's matter correlations.
"""
import numpy as np
from continuum_collector import Pulse,thermal_current


def theta_at(u,pulse=Pulse()):
    return float(pulse.values(np.array([u]))[0])


def smooth_pulse_derivative(t,pulse=Pulse(),h=1e-4):
    return (theta_at(t-2*h,pulse)-8*theta_at(t-h,pulse)
            +8*theta_at(t+h,pulse)-theta_at(t+2*h,pulse))/(12*h)


def thermal_regularizer(s,beta):
    if beta<=0: raise ValueError("positive inverse temperature required")
    x=np.pi*np.asarray(s,dtype=float)/beta
    z=np.ones_like(x)
    sel=np.abs(x)>1e-4
    z[sel]=x[sel]/np.sinh(x[sel])
    xs=x[~sel]
    z[~sel]=1-xs**2/6+7*xs**4/360-31*xs**6/15120
    return z.item() if z.ndim==0 else z


def thermal_kernel(s,beta):
    if s==0: raise ValueError("Wightman singular at coincidence; point split")
    return 1j/(2*np.pi*s)*thermal_regularizer(s,beta)


def rotate(angle):
    c=np.cos(angle);z=np.sin(angle)
    return np.array([[c,z],[-z,c]])


def outgoing_kernel(u,v,beta_H=2.,beta_B=7.,pulse=Pulse()):
    if u==v: raise ValueError("use analytic coincident limit")
    a=theta_at(u,pulse);b=theta_at(v,pulse)
    C=np.diag([thermal_kernel(u-v,beta_H),thermal_kernel(u-v,beta_B)])
    return rotate(a)@C@rotate(b).T


def point_split_currents(t,split,beta_H=2.,beta_B=7.,pulse=Pulse()):
    """Renormalized point-split stress: (1-overlap)/(2pi split²).
    Converges to mixed thermal currents + theta_dot²/(4pi) per channel.
    """
    if split==0: raise ValueError("nonzero point splitting required")
    ra=rotate(theta_at(t+split/2,pulse))
    rb=rotate(theta_at(t-split/2,pulse))
    ratio=np.array([thermal_regularizer(split,beta_H),
                    thermal_regularizer(split,beta_B)])
    return (1-(ra*rb)@ratio)/(2*np.pi*split*split)


def local_currents(t,beta_H=2.,beta_B=7.,pulse=Pulse(),h=1e-4):
    ang=theta_at(t,pulse);dot=smooth_pulse_derivative(t,pulse,h)
    jh=thermal_current(beta_H);jb=thermal_current(beta_B)
    c2=np.cos(ang)**2;s2=np.sin(ang)**2
    w=dot*dot/(4*np.pi)
    return np.array([c2*jh+s2*jb+w,s2*jh+c2*jb+w])


def cross_coincident(t,pulse=Pulse(),h=1e-4):
    """Exact limit of G_out[H,B](t+s/2,t-s/2) = i theta_dot/(2pi)."""
    return 1j*smooth_pulse_derivative(t,pulse,h)/(2*np.pi)


def offdiag_point_split(t,split,beta_H=2.,beta_B=7.,pulse=Pulse()):
    if split==0:raise ValueError("nonzero split required")
    return outgoing_kernel(t+split/2,t-split/2,beta_H,beta_B,pulse)[0,1]


def state_thermal_translation_violation(u,v,shift=1.,pulse=Pulse(),
                                        beta_H=2.,beta_B=7.):
    return (outgoing_kernel(u+shift,v+shift,beta_H,beta_B,pulse)[1,1]
            -outgoing_kernel(u,v,beta_H,beta_B,pulse)[1,1])


if __name__=="__main__":
    for t in (5.,6.,8.,12.,18.):
        print(t,local_currents(t),point_split_currents(t,1e-3),
              offdiag_point_split(t,1e-3),cross_coincident(t))
