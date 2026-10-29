import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from scipy.integrate import quad
from sympy import symbols,pi,diff,simplify
from baseline import (horizons,kappa,entropies,page_mass_proxy,
                      page_time_proxy,page_time_derivative,check_ode)


def test_symbolic_first_laws():
    r,G,L=symbols("r G L",positive=True)
    m=r*(1-L*r**2/3)/(2*G)
    S=pi*r**2/G
    Tsigned=(1-L*r**2)/(4*pi*r)
    assert simplify(diff(S,r)/diff(m,r)-1/Tsigned)==0


def test_proxy_entropy_identity():
    M0=100.; Mp=page_mass_proxy(M0)
    integral=quad(lambda m: 2*np.pi*(1/kappa(m)[1]-1/kappa(m)[0]),Mp,M0)[0]
    sb0,sc0=entropies(M0); sb,sc=entropies(Mp)
    assert np.isclose(integral,sb+sc-sb0-sc0,rtol=1e-10)
    assert np.isclose(sc,sb0+sc0,rtol=1e-10)


def test_proxy_time_linear_response():
    M0=100.; h=1e-4
    fd=(page_time_proxy(M0,eps=h)-page_time_proxy(M0,eps=-h))/(2*h)
    assert np.isclose(fd,page_time_derivative(M0),rtol=2e-8)


def test_independent_ode():
    tq,to=check_ode()
    assert np.isclose(tq,to,rtol=1e-8)


def test_manufactured_one_dimensional_qes():
    # NUMERICAL PIPELINE ONLY: not global UdS entropy or a two-null QES.
    M=100.; L=3e-6; G=1.; alpha=2.; mu=2.
    def toy(m):
        rb,_=horizons(m,L,G)
        return (rb+np.sqrt(rb**2-2*alpha*G/np.pi))/2
    rb,_=horizons(M,L,G)
    drb=2*G/(1-L*rb**2)
    dtdrb=(1+rb/np.sqrt(rb**2-2*alpha*G/np.pi))/2
    predicted=mu*drb*dtdrb
    h=1e-3
    observed=(toy(M+h*mu)-toy(M-h*mu))/(2*h)
    assert np.isclose(predicted,observed,rtol=2e-8)
