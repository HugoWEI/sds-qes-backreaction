"""Root-sensitivity validation only: NOT a QES or physical horizon prediction."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from scipy.optimize import brentq
import numpy as np
from baseline import horizons
from linearized_geometry import f0,f0_prime,stress,radial_solution

def test_perturbative_black_hole_root_vs_nonlinear():
    rb,_=horizons(100.)
    d=.005
    F=radial_solution(rb+d)[1]
    Fp=(2*stress(rb+d)["Tvr"]-F)/(rb+d)
    Fb=F-d*Fp
    delta=-Fb/f0_prime(rb)
    assert delta>0  # Only in our r_reference=400 normalization.
    for eta in [.1,.2,.5,1.]:
        root=brentq(lambda r:f0(r)+eta*radial_solution(r)[1],
                    rb+1e-8,rb+.01,xtol=1e-12)
        assert np.isclose((root-rb)/eta,delta,rtol=2e-5,atol=1e-10)
