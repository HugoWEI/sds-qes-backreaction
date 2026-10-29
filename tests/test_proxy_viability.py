import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import sympy as s
from scipy.optimize import brentq
from proxy_viability import (cosmological_y,mass_mu,kappas_scaled,
    proxy_page,adiabatic_coefficient,cosmological_adiabatic_coefficient,ProxyDiagnostic)


def test_horizon_vieta_and_mass():
    for x in [.05,.2,.5,.8,.99]:
        y=cosmological_y(x)
        assert np.isclose(x*x+x*y+y*y,3,rtol=0,atol=1e-14)
        assert np.isclose(mass_mu(x),y*(1-y*y/3)/2,rtol=1e-11)


def test_proxy_page_identity_independent_root():
    for x in [.05,.2,.5,.8,.99]:
        z=proxy_page(x)
        xp=z['x_page']
        yp=cosmological_y(xp)
        yi=cosmological_y(x)
        assert np.isclose(yp*yp,yi*yi+x*x,rtol=2e-12)
        direct_y=np.sqrt(yi*yi+x*x)
        mu_from_y=direct_y*(1-direct_y*direct_y/3)/2
        assert np.isclose(mu_from_y,z['mu_page'],rtol=2e-12)
        assert 0<z['mass_ratio']<1


def test_adiabatic_formula_direct_finite_difference():
    # G=Lambda=1; masses scale as mu. d kappa/dmu = -(1/x^2+1)/(1-x^2).
    for x in [.05,.2,.5,.8,.99]:
        mu=mass_mu(x)
        h=1e-7*min(mu,1/3-mu)
        def kb(mu0):
            xx=brentq(lambda r:mass_mu(r)-mu0,1e-12,1-1e-12)
            return kappas_scaled(xx)[0]
        dkb=(kb(mu+h)-kb(mu-h))/(2*h)
        direct=abs(dkb)*(kappas_scaled(x)[0]**2-kappas_scaled(x)[1]**2)/kappas_scaled(x)[0]**2
        assert np.isclose(direct,adiabatic_coefficient(x),rtol=4e-6)


def test_nariai_limit_symbolic():
    x=s.symbols('x',positive=True)
    y=(s.sqrt(12-3*x*x)-x)/2
    kb=(1-x*x)/(2*x)
    kc=(y*y-1)/(2*y)
    a=(1+x*x)/(x*x*(1-x*x))*(1-(kc/kb)**2)
    assert s.limit(a,x,1,dir='-')==s.Rational(8,3)
    assert s.limit(x*x*a,x,0,dir='+')==1
    for t in [.9999,.9999999]:
        assert np.isclose(adiabatic_coefficient(t),8/3,rtol=1e-6)


def test_cosmological_horizon_adiabatic_limits():
    x=s.symbols('x',positive=True)
    y=(s.sqrt(12-3*x*x)-x)/2
    kb=(1-x*x)/(2*x)
    kc=(y*y-1)/(2*y)
    a=(1+y*y)/(y*y*(y*y-1))*((kb/kc)**2-1)
    assert s.limit(a,x,1,dir='-')==s.Rational(8,3)
    assert s.limit(x*x*a,x,0,dir='+')==s.Rational(1,2)
    for xx in [.05,.2,.5,.8,.99]:
        assert cosmological_adiabatic_coefficient(xx)>0


def test_proxy_necessary_regime_diagnostics():
    for x in [.05,.2,.5,.8,.99]:
        d=ProxyDiagnostic(x).evaluate(grid=101)
        assert d['max_adiabatic']>0
        assert d['max_adiabatic']>=d['max_bh_adiabatic']
        assert d['max_adiabatic']>=d['max_cosmological_adiabatic']
        assert d['max_coupling']>0
        assert np.isclose(d['max_coupling'],1e-6/d['x_page']**2)
        assert d['max_adiabatic']<.01
        assert d['max_coupling']<.01
        assert np.isclose(d['critical_q_for_one_percent']/d['q'] * d['max_adiabatic'],.01)
