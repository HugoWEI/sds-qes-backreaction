"""Compact, tested Phase-8 static SdS worldline/Keldysh source benchmark.

Only an anchored finite-mode contact on a frozen static SdS geometry.
Neither the full curved-space fermion influence functional nor a QES.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm
from baseline import horizons
from linearized_geometry import f0, f0_prime
from keldysh_contact import Contact, influence, exact_energy_ledger, noise_and_retarded


@dataclass(frozen=True)
class StaticWorldtube:
    radius: float = 500.
    mass: float = 100.
    cosmological_constant: float = 3e-6
    G: float = 1.

    def __post_init__(self):
        rb, rc = horizons(self.mass, self.cosmological_constant, self.G)
        if not rb < self.radius < rc:
            raise ValueError('worldtube outside static patch')

    @property
    def lapse(self):
        return np.sqrt(f0(self.radius, self.mass, self.cosmological_constant, self.G))

    @property
    def signed_acceleration(self):
        return f0_prime(self.radius, self.mass, self.cosmological_constant, self.G)/(2*self.lapse)

    def support_force(self, local_energy):
        return local_energy*self.signed_acceleration

    def stress_divergence_delta_coefficients(self, energy, dE_dt=0.):
        """Coefficients of delta(r-radius) in (nabla T)^t and (nabla T)^r."""
        fprime=f0_prime(self.radius, self.mass, self.cosmological_constant, self.G)
        return dE_dt/self.lapse, energy*self.lapse*fprime/2


def lapse_ctp(contact, Aplus, Aminus, gplus, gminus, dt):
    """Exact CTP determinant, H_t,± = A_±(t) (h0+g_±(t)V)."""
    ap=np.asarray(Aplus,float);am=np.asarray(Aminus,float)
    gp=np.asarray(gplus,float);gm=np.asarray(gminus,float)
    if dt<=0 or ap.ndim!=1 or not(ap.shape==am.shape==gp.shape==gm.shape):
        raise ValueError('invalid branch histories')
    if np.any(ap<=0) or np.any(am<=0):
        raise ValueError('lapse must be positive')
    h,V,C=contact.matrices()
    def U(a,g):
        u=np.eye(len(C),dtype=complex)
        for ai,gi in zip(a,g):
            u=expm(-1j*ai*(h+gi*V)*dt)@u
        return u
    up=U(ap,gp);um=U(am,gm)
    return np.linalg.det(np.eye(len(C))-C+C@um.conj().T@up)


def exact_lapse_source(contact, lapse, drive, dt):
    """Instantaneous expectation Tr[rho_before (h+g_j V)] for each bin."""
    a=np.asarray(lapse,float);g=np.asarray(drive,float)
    if a.shape!=g.shape or dt<=0 or np.any(a<=0):raise ValueError('invalid drive')
    h,V,C=contact.matrices();u=np.eye(len(C),dtype=complex)
    source=[]
    for ai,gi in zip(a,g):
        hj=h+gi*V
        source.append(np.trace(u@C@u.conj().T@hj).real)
        u=expm(-1j*ai*hj*dt)@u
    return np.asarray(source)


def static_energy_variance(contact):
    h,V,C=contact.matrices()
    return float(np.trace(C@h).real), float(np.trace(C@h@(np.eye(len(C))-C)@h).real)


def clock_covariance_check(worldtube, contact, drive, dtau):
    """Killing/proper clock comparison on equal physical histories."""
    a=worldtube.lapse
    zero=np.zeros_like(drive)
    proper=influence(contact,drive,zero,dtau)
    killing=lapse_ctp(contact,np.full(len(drive),a),np.full(len(drive),a),drive,zero,dtau/a)
    local_work=exact_energy_ledger(contact,drive,dtau)['external_work']
    tau=np.linspace(0,2,6)
    N,R=noise_and_retarded(contact,tau)
    return dict(overlap_difference=abs(proper-killing),
                local_work=local_work,killing_work=a*local_work,
                local_noise=N,killing_noise=a*a*N,
                local_response=R,killing_response=a*a*R)
