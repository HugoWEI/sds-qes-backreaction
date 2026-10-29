"""Symbolic and independent smooth-shell distributional checks."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
import sympy as sp
from scipy.integrate import quad
from curved_defect import Interface
from defect_junction_constraints import ShellBattery,shell_distributional_coefficients

def test_exact_symbolic_gravity_shell_coefficients():
    l,L=sp.symbols("ell Lambda",real=True)
    A=sp.Function("A")(l);X=sp.Function("X")(l)
    BX=sp.diff(X,l,2)+sp.diff(A,l)/A*sp.diff(X,l)
    E_tt=sp.simplify(-A**2*BX+A*sp.diff(A,l)*sp.diff(X,l)+
        A**2*sp.diff(X,l)**2/(4*X)+A**2*(1-L*X))
    E_X=-2*sp.diff(A,l,2)/A+sp.diff(X,l)**2/(2*X**2)-BX/X-2*L
    assert sp.simplify(sp.diff(E_tt,sp.diff(X,l,2))+A**2)==0
    assert sp.simplify(sp.diff(E_X,sp.diff(X,l,2))+1/X)==0
    assert sp.simplify(sp.diff(E_X,sp.diff(A,l,2))+2/A)==0

def test_internal_shell_dilaton_junction_only_with_no_charge():
    A=.6;X=250000.;G=1.;epsilon=.0001
    dX=-2*G*epsilon;dA=A*G*epsilon/X
    C=shell_distributional_coefficients(A,X,dX,dA)
    assert np.isclose(C["Ett"],2*G*(A*A*epsilon),rtol=2e-15)
    assert abs(C["EX"])<1e-20

def test_shell_battery_energy_exhaustion_and_recoil():
    a=ShellBattery().trajectory()
    assert np.max(abs(a["integrated_conservation"]))<1e-13
    assert np.all(a["energy_K"]>0)
    assert np.all(np.diff(a["used_work"])>=-1e-17)
    assert np.isclose(a["energy_K"][-1],.5*a["energy_K"][0],rtol=1e-12)
    assert abs(a["jump_X_normal_derivative"][-1])<abs(a["jump_X_normal_derivative"][0])
    assert abs(a["jump_extrinsic_curvature"][-1])<abs(a["jump_extrinsic_curvature"][0])

def test_local_battery_conservation_derivative():
    a=ShellBattery().trajectory();t=a["killing_time"];E=a["energy_K"];P=a["power_K"]
    num=(E[2:]-E[:-2])/(t[2:]-t[:-2])
    assert np.max(abs(num+P[1:-1]))<1e-10

def test_bianchi_external_work_no_source_inconsistent():
    w=Interface().local_ward()
    assert np.max(w["switching_power"])>0
    assert np.max(abs(w["collector_gain"]+w["patch_after"]-w["patch_before"]-
                      w["switching_power"]))<1e-20
    b=ShellBattery().trajectory()
    assert np.allclose(np.diff(b["energy_K"]),-np.diff(b["used_work"]),atol=1e-14)

def test_delta_approximation_smoothing():
    A0=.6;X0=250000.;dX=-.0002;dA=A0*.0001/X0
    ref=shell_distributional_coefficients(A0,X0,dX,dA)
    w=1e-3
    ix=quad(lambda x:dX/(2*w)/np.cosh(x/w)**2,-.02,.02,
            epsabs=1e-18,epsrel=1e-13)[0]
    ia=quad(lambda x:dA/(2*w)/np.cosh(x/w)**2,-.02,.02,
            epsabs=1e-22,epsrel=1e-13)[0]
    assert np.isclose(-A0**2*ix,ref["Ett"],rtol=1e-12)
    assert np.isclose(-2*ia/A0-ix/X0,ref["EX"],atol=1e-20)
