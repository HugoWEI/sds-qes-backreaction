"""Gaussian free-Dirac state of switched contact in a finite NS circle basis.

The n-mode antiperiodic Fourier covariance C=F f_beta(omega) F^dagger
is rotated by an exact local SO(2) contact on time bins. A declared
canonical purifier of the hot input allows exact mutual information;
it is NOT a physical SdS black-hole interior or gravitational QES.
"""
import numpy as np
from scipy.special import expit
from continuum_collector import Pulse


def circle_thermal_covariance(n=64,L=24.,beta=2.):
    if n%2 or n<8 or beta<=0:raise ValueError("even n>=8, beta>0 required")
    j=np.arange(n);m=np.arange(-n//2,n//2)
    omega=2*np.pi*(m+.5)/L
    F=np.exp(-2j*np.pi*j[:,None]*(m[None,:]+.5)/n)/np.sqrt(n)
    return (F*expit(-beta*omega)[None,:])@F.conj().T


def entropy(C):
    eig=np.clip(np.linalg.eigvalsh((C+C.conj().T)/2).real,1e-14,1-1e-14)
    return float(-np.sum(eig*np.log(eig)+(1-eig)*np.log1p(-eig)))


def purification_cross(C):
    eig,V=np.linalg.eigh((C+C.conj().T)/2)
    return (V*np.sqrt(np.maximum(eig*(1-eig),0))[None,:])@V.conj().T


def transmitted_state(n=64,L=24.,beta_H=2.,beta_B=7.,pulse=Pulse()):
    if not np.isclose(L,pulse.period):raise ValueError("L must match pulse period")
    CH=circle_thermal_covariance(n,L,beta_H)
    CB=circle_thermal_covariance(n,L,beta_B)
    times=np.arange(n)*L/n;angle=pulse.values(times)
    co=np.cos(angle);si=np.sin(angle)
    Ho=np.outer(co,co)*CH+np.outer(si,si)*CB
    Bo=np.outer(si,si)*CH+np.outer(co,co)*CB
    cross=-np.outer(co,si)*CH+np.outer(si,co)*CB
    return dict(hot=Ho,bath=Bo,cross=cross,initial_hot=CH,initial_bath=CB,
                time=times,angle=angle,transmission=si*si)


def state_observables(n=64,L=24.,beta_H=2.,beta_B=7.,pulse=Pulse()):
    z=transmitted_state(n,L,beta_H,beta_B,pulse)
    H=z["hot"];B=z["bath"];X=z["cross"]
    CH=z["initial_hot"];CB=z["initial_bath"]
    R=np.eye(n)-CH;D=purification_cross(CH)
    BR=-np.sin(z["angle"])[:,None]*D
    joint_BR=np.block([[B,BR],[BR.conj().T,R]])
    joint_HB=np.block([[H,X],[X.conj().T,B]])
    return dict(n=n,period=L,hot_entropy_before=entropy(CH),
      hot_entropy_after=entropy(H),bath_entropy_before=entropy(CB),
      bath_entropy_after=entropy(B),bath_entropy_change=entropy(B)-entropy(CB),
      mutual_information_channels=entropy(H)+entropy(B)-entropy(joint_HB),
      mutual_information_bath_hot_reference=entropy(B)+entropy(R)-entropy(joint_BR),
      global_entropy_invariant_residual=entropy(joint_HB)-entropy(CH)-entropy(CB),
      correlation_norm=float(np.linalg.norm(X)),
      local_max_transmission=float(np.max(z["transmission"])),
      covariance_eigenvalue_bounds=(float(np.min(np.linalg.eigvalsh(joint_HB))),
                                  float(np.max(np.linalg.eigvalsh(joint_HB)))))


if __name__=="__main__":
    for n in (32,64,96,128):
        a=state_observables(n=n)
        print(n,a["bath_entropy_change"],
              a["mutual_information_bath_hot_reference"],
              a["mutual_information_channels"],
              a["global_entropy_invariant_residual"])
