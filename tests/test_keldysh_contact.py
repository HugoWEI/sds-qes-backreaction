"""CTP fermion determinant, retarded Kubo kernel, noise, exact work audit."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pytest
from scipy.linalg import expm
from keldysh_contact import (Contact,thermal_covariance,evolve,influence,
    force_onebody,noise_and_retarded,integrated_force,linear_noise_variance,
    retarded_integral_final,driven_force_at_end,exact_energy_ledger)


def pulse(n=28,dt=.15):
    t=(np.arange(n)+.5)*dt
    return np.sin(np.pi*t/(n*dt))**2


def test_hamiltonians_hermitian_gaussian_state_positive():
    h,V,C=Contact().matrices()
    np.testing.assert_allclose(h,h.conj().T,atol=1e-15)
    np.testing.assert_allclose(V,V.conj().T,atol=1e-15)
    assert np.linalg.eigvalsh(C).min()>0
    assert np.linalg.eigvalsh(C).max()<1
    assert np.linalg.norm(h@C-C@h)<1e-11


def test_closed_time_path_normalization_and_branch_conjugation():
    model=Contact(sites=3);p=pulse(14,.19)
    plus=.37*p;minus=-.13*p
    assert abs(influence(model,plus,plus,.19)-1)<1e-12
    assert abs(influence(model,np.zeros_like(p),np.zeros_like(p),.19)-1)<1e-12
    assert np.isclose(influence(model,plus,minus,.19),
                      np.conj(influence(model,minus,plus,.19)),atol=1e-12)
    assert abs(influence(model,plus,minus,.19))<=1+1e-12


def fock_annihilation(mode,number):
    """Independent Jordan-Wigner/Fock-space construction."""
    size=1<<number
    a=np.zeros((size,size),complex)
    for bit in range(size):
        if bit&(1<<mode):
            flipped=bit^(1<<mode)
            a[flipped,bit]=(-1)**((bit&((1<<mode)-1)).bit_count())
    return a


def second_quantize(a,ops):
    return sum((a[i,j]*ops[i].conj().T@ops[j]
                for i in range(len(ops)) for j in range(len(ops))),
               np.zeros_like(ops[0]))


def test_klich_determinant_against_independent_full_fock_trace():
    model=Contact(sites=2)
    h,V,C=model.matrices()
    ops=[fock_annihilation(j,len(h)) for j in range(len(h))]
    H=second_quantize(h,ops);W=second_quantize(V,ops)
    n=model.sites
    Kh=np.zeros_like(h);Kb=np.zeros_like(h)
    Kh[:n,:n]=h[:n,:n];Kb[n:,n:]=h[n:,n:]
    rho=expm(-model.beta_hot*second_quantize(Kh,ops)-
             model.beta_bath*second_quantize(Kb,ops))
    rho/=np.trace(rho)
    plus=np.array([.2,.4,-.2]);minus=np.array([-.1,.25,.15]);dt=.23
    up=evolve(H,W,plus,dt);um=evolve(H,W,minus,dt)
    FF=np.trace(up@rho@um.conj().T)
    FD=influence(model,plus,minus,dt)
    assert np.isclose(FF,FD,atol=1e-12)


def test_exact_switch_work_equals_both_channel_energy_changes():
    model=Contact(sites=5)
    ledger=exact_energy_ledger(model,.45*pulse(38,.1),.1)
    assert abs(ledger['residual'])<1e-12
    assert np.isclose(ledger['external_work'],ledger['hot_energy_change']+
                      ledger['bath_energy_change'],atol=1e-12)
    assert ledger['external_work']>0
    assert ledger['bath_energy_change']>0


def test_zero_pulse_has_zero_energy_and_identity():
    p=np.zeros(10);m=Contact(sites=3)
    a=exact_energy_ledger(m,p,.2)
    assert abs(a['external_work'])<1e-12
    assert abs(a['total_energy_change'])<1e-12
    assert abs(influence(m,p,p,.2)-1)<1e-12


def test_symmetrized_noise_is_real_symmetric_psd():
    m=Contact(sites=4)
    times=np.linspace(0,4,32)
    N,R=noise_and_retarded(m,times)
    np.testing.assert_allclose(N,N.T,atol=2e-13)
    assert np.linalg.eigvalsh(N).min()>-1e-11
    assert np.all(np.diag(N)>=-1e-14)
    assert np.allclose(R[np.triu_indices(len(R))],0,atol=1e-14)


def test_retarded_response_independent_finite_difference():
    m=Contact(sites=4);dt=.15;p=pulse(28,dt)
    pred=retarded_integral_final(m,p,dt)
    eps=.0005
    fd=(driven_force_at_end(m,eps*p,dt)-
        driven_force_at_end(m,-eps*p,dt))/(2*eps)
    assert abs(pred)>1e-4
    assert np.isclose(fd,pred,rtol=3e-5,atol=1e-8)


def test_decoherence_quadratic_term_matches_gaussian_noise():
    m=Contact(sites=4);dt=.15;p=pulse(28,dt)
    var=linear_noise_variance(m,p,dt)
    assert var>0
    for eps in (.005,.002):
        f=influence(m,eps*p,np.zeros_like(p),dt)
        sample=-2*np.log(abs(f))/eps**2
        assert np.isclose(sample,var,rtol=1e-4)


def test_influence_magnitude_only_state_and_history_dependent():
    m=Contact(sites=3)
    p=pulse(18,.15)
    a=influence(m,.3*p,np.zeros_like(p),.15)
    b=influence(m,.3*p,-.3*p,.15)
    assert not np.isclose(abs(a),abs(b),atol=1e-3)


def test_quadratic_noise_integral_exact_generator_is_hermitian():
    m=Contact(sites=4)
    A=integrated_force(m,pulse(18,.12),.12)
    np.testing.assert_allclose(A,A.conj().T,atol=2e-14)
    assert np.linalg.norm(A)>0


def test_future_coupling_cannot_affect_past_measured_force():
    model=Contact(sites=4);h,V,C=model.matrices();dt=.1
    x=.4*pulse(20,dt)
    a=evolve(h,V,x[:8],dt)
    y=x.copy();y[8:]+=5
    b=evolve(h,V,y[:8],dt)
    np.testing.assert_allclose(a,b,atol=1e-14)


def test_bad_inputs_rejected():
    model=Contact(sites=3)
    with pytest.raises(ValueError):evolve(*model.matrices()[:2],[1.],-1)
    with pytest.raises(ValueError):influence(model,[.1],[.1,.2],.1)
    with pytest.raises(ValueError):Contact(sites=1).matrices()
    with pytest.raises(ValueError):thermal_covariance(np.eye(2),0)


def test_equilibrium_kms_detailed_balance_and_nonequilibrium_violation():
    from keldysh_contact import spectral_detailed_balance
    eq=spectral_detailed_balance(Contact(sites=5,beta_hot=3.,beta_bath=3.))
    assert len(eq)>3
    assert np.max(abs(np.log(eq[:,1]/eq[:,2])-3*eq[:,0]))<2e-12
    neq=spectral_detailed_balance(Contact(sites=5,beta_hot=2.,beta_bath=5.))
    assert np.std(np.log(neq[:,1]/neq[:,2])/neq[:,0])>.1
