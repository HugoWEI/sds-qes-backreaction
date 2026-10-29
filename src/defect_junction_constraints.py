"""Conditional 2D dilaton internal-shell junction and battery stress audit.

A finite-duration chiral contact requires an energy supplier; omitting the
controller breaks total stress conservation. In Gaussian normal coordinates
ds²=dell²-A(ell,t)²dt² with continuous X, the delta parts of the spherical
2D gravitational equations are E_tt=-A0²[X_ell] delta(ell) and
E_X=(-2[K]-[X_ell]/X0) delta(ell). For an INTERNAL timelike shell with
T_tt=A0² epsilon delta(ell), no X-dependent surface charge:
[X_ell]=-2G epsilon, [K]=G epsilon/X.
These are NOT boundary conditions for the real external nongravitating bath.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import cumulative_trapezoid
from curved_defect import Interface,killing_curvature_scale

@dataclass(frozen=True)
class ShellBattery:
    interface:Interface=Interface()
    initial_work_multiple:float=2.

    def trajectory(self):
        if self.initial_work_multiple<1:
            raise ValueError("initial battery less than planned work")
        w=self.interface.local_ward()
        W=self.interface.integrated_ward()["work"]
        t=w["x"]/w["scale"]
        used=cumulative_trapezoid(w["switching_power"],t,initial=0.)
        lapse=killing_curvature_scale(self.interface.b,self.interface.mass,
                                      self.interface.L,self.interface.G)
        E0=self.initial_work_multiple*W
        E=E0-used
        eps=E/lapse
        G=self.interface.G;X=self.interface.b**2
        return dict(killing_time=t,power_K=w["switching_power"],
           used_work=used,energy_K=E,epsilon_proper=eps,
           jump_X_normal_derivative=-2*G*eps,
           jump_extrinsic_curvature=G*eps/X,
           integrated_conservation=E+used-E0,
           initial_killing_energy=E0,redshift=lapse,
           jump_dilaton_residual=-2*G*eps+2*G*eps,
           jump_scalar_residual=-2*G*eps/X+2*G*eps/X)

def shell_distributional_coefficients(A0,X0,delta_X_deriv,delta_A_deriv):
    if A0<=0 or X0<=0: raise ValueError("positive lapse/dilaton required")
    return dict(Ett=-A0**2*delta_X_deriv,
                EX=-2*delta_A_deriv/A0-delta_X_deriv/X0)

if __name__=="__main__":
    a=ShellBattery().trajectory()
    print("controller beginning/end Killing energy",a["energy_K"][0],a["energy_K"][-1])
    print("conservation residual",np.max(abs(a["integrated_conservation"])))
