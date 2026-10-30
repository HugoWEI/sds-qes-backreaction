"""Exact finite-mode *real-time* Gaussian fermion influence functional.

Two coupled tight-binding channels, one local Hermitian contact V:
  h(t) = h0 + g(t) V.
For two controller histories g_+,g_-, the closed-time-path overlap is
  F[g_+,g_-] = det(I - C0 + C0 U_-^dagger U_+),
where C0 is the initial Gaussian number-conserving covariance and U_±
are independently time ordered one-particle unitaries. This is the
Klich determinant, an established formula. No continuum, curvature,
Polyakov subtraction, bulk gravity, or QES is claimed.

All units and energy are those of finite open-boundary lattice chains.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm
from scipy.special import expit


def thermal_covariance(h, beta):
    if beta <= 0: raise ValueError('beta must be positive')
    energies, basis = np.linalg.eigh(h)
    return (basis * expit(-beta * energies)[None, :]) @ basis.conj().T


@dataclass(frozen=True)
class Contact:
    sites: int = 5
    hopping: float = 0.7
    beta_hot: float = 2.0
    beta_bath: float = 5.0
    bath_offset: float = 0.18
    contact_site: int | None = None

    def matrices(self):
        n = self.sites
        if n < 2 or self.hopping <= 0: raise ValueError('invalid chain')
        h = np.zeros((n, n), dtype=complex)
        for j in range(n-1):
            h[j, j+1] = h[j+1,j] = -self.hopping
        h += np.diag([0.08 * (j-n/2)/n for j in range(n)])
        hb = h + self.bath_offset * np.eye(n)
        h0 = np.block([[h, np.zeros_like(h)], [np.zeros_like(hb), hb]])
        j = n//2 if self.contact_site is None else self.contact_site
        if not (0 <= j < n): raise ValueError('contact outside lattice')
        p = np.zeros((n,n), dtype=complex);p[j,j] = 1
        v = np.block([[np.zeros_like(h), 1j*p], [-1j*p, np.zeros_like(h)]])
        C0 = np.block([[thermal_covariance(h,self.beta_hot),np.zeros_like(h)],
                       [np.zeros_like(hb),thermal_covariance(hb,self.beta_bath)]])
        return h0,v,C0


def evolve(h0,V,history,dt,collect=False):
    """Time-ordered piecewise constant drive at each interval midpoint."""
    if dt<=0: raise ValueError('dt positive')
    U=np.eye(h0.shape[0],dtype=complex)
    intermediates=[]
    for g in history:
        U=expm(-1j*(h0+float(g)*V)*dt)@U
        if collect: intermediates.append(U)
    return (U,intermediates) if collect else U


def influence(contact,plus,minus,dt):
    """F[g+,g-]=Tr(U+ rho0 U-^dagger) as a one-particle determinant."""
    plus=np.asarray(plus,float);minus=np.asarray(minus,float)
    if plus.shape!=minus.shape or plus.ndim!=1:raise ValueError('mismatched histories')
    h0,V,C0=contact.matrices()
    Up=evolve(h0,V,plus,dt)
    Um=evolve(h0,V,minus,dt)
    M=np.eye(len(C0))-C0+C0@Um.conj().T@Up
    return np.linalg.det(M)


def force_onebody(h0,V,t):
    """Unperturbed interaction-picture contact force V_I(t)."""
    u=expm(-1j*h0*t)
    return u.conj().T@V@u


def noise_and_retarded(contact,times):
    """Symmetrized connected force noise and retarded Kubo response.

    N(t,s) = Re Tr[C V_I(t) (1-C) V_I(s)].
    R(t,s) = -i theta(t-s) Tr[C [V_I(t),V_I(s)]].
    Purely finite quantum-system correlators; no Markov/thermal-bath limit.
    """
    h,V,C=contact.matrices()
    times=np.asarray(times,float)
    Vs=[force_onebody(h,V,t) for t in times]
    I=np.eye(len(C))
    N=np.zeros((len(times),len(times)),float)
    R=np.zeros_like(N)
    for i,A in enumerate(Vs):
        for j,B in enumerate(Vs):
            N[i,j]=np.trace(C@A@(I-C)@B).real
            if times[i]>times[j]:
                R[i,j]=(-1j*np.trace(C@(A@B-B@A))).real
    return N,R


def integrated_force(contact,history,dt):
    """Exact first-order drive generator A=integral g(t) V_I(t) dt.

    Stable analytic integral of frequencies, not midpoint quadrature.
    """
    h,V,C=contact.matrices()
    energies,P=np.linalg.eigh(h)
    W=P.conj().T@V@P
    dw=energies[:,None]-energies[None,:]
    t=(np.arange(len(history))+.5)*dt
    A=np.zeros_like(W)
    for center,amp in zip(t,history):
        z=dw*dt/2
        A+=float(amp)*dt*np.sinc(z/np.pi)*np.exp(1j*dw*center)*W
    return P@A@P.conj().T


def linear_noise_variance(contact,history,dt):
    h,V,C=contact.matrices()
    A=integrated_force(contact,history,dt)
    return np.trace(C@A@(np.eye(len(C))-C)@A).real


def retarded_integral_final(contact,history,dt):
    h,V,C=contact.matrices()
    A=integrated_force(contact,history,dt)
    VT=force_onebody(h,V,len(history)*dt)
    return (-1j*np.trace(C@(VT@A-A@VT))).real


def driven_force_at_end(contact,history,dt):
    h,V,C=contact.matrices()
    U=evolve(h,V,history,dt)
    return np.trace(U@C@U.conj().T@V).real


def exact_energy_ledger(contact,history,dt):
    """Exact time-resolved switching work for stepwise g_j(t).

    At each instantaneous jump g_previous -> g_new, the state is fixed and
    work is Delta g * <V>. Between jumps energy of h0+gV is conserved.
    Final jump returns g to zero. Thus W = Delta Tr(h0 C) exactly.
    """
    history=np.asarray(history,float)
    h,V,C=contact.matrices()
    n=contact.sites
    U=np.eye(len(C),dtype=complex)
    old_g=0.
    work=0.
    for g in history:
        state=U@C@U.conj().T
        work+=(g-old_g)*np.trace(state@V).real
        U=expm(-1j*(h+g*V)*dt)@U
        old_g=g
    state=U@C@U.conj().T
    work-=old_g*np.trace(state@V).real
    deltaE=np.trace(h@(state-C)).real
    dh=np.trace(h[:n,:n]@(state[:n,:n]-C[:n,:n])).real
    db=np.trace(h[n:,n:]@(state[n:,n:]-C[n:,n:])).real
    return dict(external_work=work,total_energy_change=deltaE,
                hot_energy_change=dh,bath_energy_change=db,
                residual=work-deltaE,final_covariance=state)


def demo():
    m=Contact()
    nsteps=48;dt=.12
    x=(np.arange(nsteps)+.5)*dt
    profile=0.45*np.sin(np.pi*x/(nsteps*dt))**2
    led=exact_energy_ledger(m,profile,dt)
    N,R=noise_and_retarded(m,np.linspace(0,nsteps*dt,24))
    return dict(ledger={k:v for k,v in led.items() if k!='final_covariance'},
        influence_equal=influence(m,profile,profile,dt),
        influence_forward=influence(m,profile,np.zeros_like(profile),dt),
        noise_min_eigenvalue=float(np.linalg.eigvalsh(N).min()),
        response_earliest_to_later=float(R[-1,0]),
        force_linear=retarded_integral_final(m,profile,dt))


if __name__=='__main__':
    print(demo())


def spectral_detailed_balance(contact):
    """Unsymmetrized force-noise weights for positive energy differences.

    Records (omega, w_plus, w_minus); w_plus/w_minus=exp(beta*omega)
    at a common equilibrium beta, even for channels with different offsets.
    At different temperatures there is generally NO single KMS beta.
    """
    h,V,C=contact.matrices()
    E,P=np.linalg.eigh(h)
    Y=P.conj().T@V@P
    f=np.real(np.diag(P.conj().T@C@P))
    records=[]
    for i in range(len(E)):
        for j in range(len(E)):
            omega=E[j]-E[i]
            strength=abs(Y[i,j])**2
            if omega>.05 and strength>1e-12:
                records.append((omega,f[i]*(1-f[j])*strength,
                                f[j]*(1-f[i])*strength))
    return np.asarray(records,float)
