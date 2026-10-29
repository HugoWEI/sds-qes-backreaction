"""Finite-duration two-channel free-fermion collector benchmark.

Exact Gaussian number-conserving local unitary, NOT a solved SdS
gravity-bath junction or a continuum dynamical conformal defect.
A finite spatial/window profile models regularized time bins.
External switching work is explicitly accounted for.

C_ij=<c_j^dagger c_i>, E=Tr(hC).
U_j=[[cos(theta_j),sin(theta_j)],[-sin(theta_j),cos(theta_j)]].
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import expit


@dataclass(frozen=True)
class CollectorParameters:
    n: int = 96
    beta_hot: float = 2.0
    beta_cold: float = 8.0
    peak_angle: float = np.pi/4
    start: int = 28
    stop: int = 68
    taper: int = 10


def hopping(n):
    if n < 4: raise ValueError("n >= 4")
    h=np.zeros((n,n))
    x=np.arange(n-1)
    h[x,x+1]=h[x+1,x]=-0.5
    return h


def thermal_covariance(h,beta):
    if beta<=0:raise ValueError("beta must be positive")
    energies,vectors=np.linalg.eigh(h)
    return (vectors*expit(-beta*energies)[None,:])@vectors.T


def switching_profile(n=96,start=28,stop=68,taper=10,peak_angle=np.pi/4):
    if not(0<=start<stop<=n and 1<=taper<=(stop-start)/2):
        raise ValueError("invalid finite window")
    x=np.arange(n)
    t=np.zeros(n)
    a=(x>=start)&(x<start+taper)
    b=(x>=start+taper)&(x<stop-taper)
    c=(x>=stop-taper)&(x<stop)
    t[a]=np.sin(np.pi*(x[a]-start)/(2*taper))**2
    t[b]=1
    t[c]=np.sin(np.pi*(stop-x[c])/(2*taper))**2
    return peak_angle*t


def channel_covariances(C_hot,C_cold,theta):
    c=np.cos(theta);s=np.sin(theta)
    hot=np.outer(c,c)*C_hot+np.outer(s,s)*C_cold
    bath=np.outer(s,s)*C_hot+np.outer(c,c)*C_cold
    cross=-np.outer(c,s)*C_hot+np.outer(s,c)*C_cold
    return hot,bath,cross


def full_unitary(theta):
    c=np.diag(np.cos(theta));s=np.diag(np.sin(theta))
    return np.block([[c,s],[-s,c]])


def entropy(C):
    z=np.clip(np.linalg.eigvalsh(C),1e-15,1-1e-15)
    return float(-np.sum(z*np.log(z)+(1-z)*np.log1p(-z)))


def energy(h,C):
    return float(np.trace(h@C))


def work_formula(h,C_hot,C_cold,theta):
    diff=theta[:,None]-theta[None,:]
    return float(np.sum(h*(np.cos(diff)-1)*(C_hot+C_cold).T))


def detailed_budget(h,C_hot,C_cold,theta):
    A,B,cross=channel_covariances(C_hot,C_cold,theta)
    dA=energy(h,A)-energy(h,C_hot)
    dB=energy(h,B)-energy(h,C_cold)
    c=np.cos(theta);s=np.sin(theta)
    contrast=energy(h,np.outer(s,s)*(C_hot-C_cold))
    edge=np.cos(theta[:,None]-theta[None,:])-1
    bath_work=energy(h,edge*C_cold)
    hot_work=energy(h,edge*C_hot)
    work=work_formula(h,C_hot,C_cold,theta)
    region=np.flatnonzero(theta!=0)
    return dict(delta_energy_hot=dA,delta_energy_collector=dB,
        net_external_switch_work=work,
        thermal_contrast_to_collector=contrast,
        collector_switch_work=bath_work,hot_switch_work=hot_work,
        total_channel_mutual_information=entropy(A)+entropy(B)-entropy(C_hot)-entropy(C_cold),
        global_entropy_before=entropy(C_hot)+entropy(C_cold),
        global_entropy_after=entropy(np.block([[A,cross],[cross.T,B]])),
        window_entropy_before=entropy(C_cold[np.ix_(region,region)]),
        window_entropy_after=entropy(B[np.ix_(region,region)]),
        transmission_max=float(np.max(s*s)),
        support_bins=len(region))


def run(p=CollectorParameters()):
    h=hopping(p.n)
    H=thermal_covariance(h,p.beta_hot)
    C=thermal_covariance(h,p.beta_cold)
    t=switching_profile(p.n,p.start,p.stop,p.taper,p.peak_angle)
    return detailed_budget(h,H,C,t)


if __name__=="__main__":
    for taper in (2,5,10,16):
        r=run(CollectorParameters(taper=taper))
        print("taper",taper,r)
        print("energy residual",r["delta_energy_hot"]+r["delta_energy_collector"]-
              r["net_external_switch_work"])
