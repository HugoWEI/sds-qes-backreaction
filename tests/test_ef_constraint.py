import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from sympy import symbols,Function,Matrix,diff,exp,simplify
from ef_constraint import (blackening,static_chiral_tuu,ef_trr,sigma_prime,
                           sigma_difference,sigma_difference_ode)


def test_symbolic_rr_equation():
    v,r=symbols("v r",real=True)
    F=Function("F")(v,r); sig=Function("sig")(v,r)
    g=Matrix([[-exp(2*sig)*F,exp(sig)],[exp(sig),0]])
    inv=simplify(g.inv()); coords=(v,r); X=r*r
    gamma_r_rr=simplify(sum(inv[1,k]*(2*diff(g[k,1],r)
                       -diff(g[1,1],coords[k]))/2 for k in range(2)))
    e_rr=simplify(-(diff(X,r,2)-gamma_r_rr*diff(X,r))
                  +diff(X,r)**2/(2*X))
    assert simplify(e_rr-2*r*diff(sig,r))==0


def test_polyakov_stress_nonzero():
    r=500.
    assert blackening(r)>0
    assert abs(static_chiral_tuu(r))>1e-12
    assert np.isclose(ef_trr(r),4*static_chiral_tuu(r)/blackening(r)**2)
    assert np.isclose(2*r*sigma_prime(r),2*ef_trr(r))


def test_quadrature_vs_ode():
    a,b=520.,500.
    assert np.isclose(sigma_difference(a,b),sigma_difference_ode(a,b),
                      rtol=1e-9,atol=1e-14)
