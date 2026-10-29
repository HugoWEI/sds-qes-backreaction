"""Gaussian thermofield purification for finite-window collector.

An exact free-fermion ancilla purification, NOT the SdS black-hole
interior or a gravitational island. Shows collector information is
not determined by mean energy flux alone.
"""
import numpy as np
from scipy.optimize import brentq
from finite_collector import (hopping,thermal_covariance,switching_profile,
    channel_covariances,entropy,CollectorParameters)


def sds_kappas(m=100.,L=3e-6,G=1.):
    f=lambda r:1-2*G*m/r-L*r*r/3
    turn=1/np.sqrt(L)
    rb=brentq(f,G*m,turn);rc=brentq(f,turn,np.sqrt(3/L))
    return (1-L*rb*rb)/(2*rb),(L*rc*rc-1)/(2*rc)


def canonical_purification_cross(C):
    w,V=np.linalg.eigh(C)
    return (V*np.sqrt(np.clip(w*(1-w),0,None))[None,:])@V.T


def purified_channel(C):
    D=canonical_purification_cross(C)
    return np.block([[C,D],[D,np.eye(len(C))-C]])


def collector_hot_reference_information(C_hot,C_cold,theta,region=None):
    n=len(theta)
    region=np.arange(n) if region is None else np.asarray(region,int)
    _,B,_=channel_covariances(C_hot,C_cold,theta)
    D=canonical_purification_cross(C_hot)
    cross=-np.sin(theta[region,None])*D[region,:]
    Breg=B[np.ix_(region,region)]
    ref=np.eye(n)-C_hot
    joint=np.block([[Breg,cross],[cross.T,ref]])
    return entropy(Breg)+entropy(ref)-entropy(joint)


def collector_cold_reference_information(C_hot,C_cold,theta,region=None):
    n=len(theta)
    region=np.arange(n) if region is None else np.asarray(region,int)
    _,B,_=channel_covariances(C_hot,C_cold,theta)
    D=canonical_purification_cross(C_cold)
    cross=np.cos(theta[region,None])*D[region,:]
    Breg=B[np.ix_(region,region)]
    ref=np.eye(n)-C_cold
    joint=np.block([[Breg,cross],[cross.T,ref]])
    return entropy(Breg)+entropy(ref)-entropy(joint)


def matching_temperature_ratio(m=100.,L=3e-6,G=1.,beta_hot=2.):
    kb,kc=sds_kappas(m,L,G)
    return beta_hot*kb/kc


def run(p=CollectorParameters()):
    h=hopping(p.n)
    H=thermal_covariance(h,p.beta_hot)
    C=thermal_covariance(h,p.beta_cold)
    t=switching_profile(p.n,p.start,p.stop,p.taper,p.peak_angle)
    region=np.flatnonzero(t>1e-12)
    return dict(beta_hot=p.beta_hot,beta_cold=p.beta_cold,
        window_bins=len(region),
        I_collector_window_hot_reference=collector_hot_reference_information(H,C,t,region),
        I_collector_window_cold_reference=collector_cold_reference_information(H,C,t,region),
        I_entire_collector_hot_reference=collector_hot_reference_information(H,C,t))


if __name__=="__main__":
    for a in (0,np.pi/8,np.pi/4,np.pi/2):
        p=CollectorParameters(peak_angle=a,beta_cold=matching_temperature_ratio())
        print("angle",a,run(p))
