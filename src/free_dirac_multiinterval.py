"""Known Casini-Huerta two-interval entropy of c free Dirac fields.
This is a *conditional field-theory benchmark*, not an island or Page curve.
Uses the declared conformal vacuum U=-exp(-kb*u), V=exp(kc*v)
and the earlier two-patch first-order metric response.
"""
import numpy as np
from baseline import kappa
from two_patch import tortoise,coordinate_correction
from interval_entropy import linearized_rho

def _state(v,r,eps=0.,r_reference=400.,**kw):
    kb,kc=kappa(kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.))
    u=v-2*tortoise(r,r_reference,**kw)
    chi=coordinate_correction(v,r,r_reference,**kw) if eps!=0 else 0.
    u+=eps*chi
    U=-np.exp(-kb*u);V=np.exp(kc*v)
    rho0,rho1=linearized_rho(v,r,r_reference,**kw)
    rho=rho0+eps*rho1-.5*(np.log(kb*kc)-kb*u+kc*v)
    return U,V,rho

def _L(z,i,j):
    d=abs(z[i]-z[j])
    if d<=0:raise ValueError('coincident endpoints')
    return np.log(d)

def _twointerval_log(z):
    return (_L(z,0,1)+_L(z,2,3)+_L(z,1,2)+_L(z,0,3)-
            _L(z,0,2)-_L(z,1,3))

def _cross_ratio(z):
    return (z[1]-z[0])*(z[3]-z[2])/((z[2]-z[0])*(z[3]-z[1]))

def _endpoint_state(points,eps=0.,**kw):
    if len(points)!=4:raise ValueError('four endpoints needed')
    vals=np.array([_state(v,r,eps=eps,**kw) for v,r in points],float)
    for k in (0,1):
        d=np.diff(vals[:,k])
        if not(np.all(d>0) or np.all(d<0)):
            raise ValueError('points not consistently spacelike ordered')
    return vals

def entropy_two_intervals(points,eps=0.,cutoffs=(1.,1.,1.,1.),**kw):
    z=_endpoint_state(points,eps,**kw)
    a=np.asarray(cutoffs,float)
    if np.any(a<=0):raise ValueError('nonpositive cutoff')
    return kw.get('N',12)/6*(_twointerval_log(z[:,0])+
            _twointerval_log(z[:,1])+z[:,2].sum()-np.log(a).sum())

def entropy_one_interval(p,q,eps=0.,cutoffs=(1.,1.),**kw):
    a=np.array(_state(*p,eps=eps,**kw))
    b=np.array(_state(*q,eps=eps,**kw))
    return kw.get('N',12)/6*(np.log(abs((a[0]-b[0])*(a[1]-b[1])))+
                             a[2]+b[2]-np.log(cutoffs[0]*cutoffs[1]))

def mutual_information(points,eps=0.,**kw):
    z=_endpoint_state(points,eps,**kw)
    xu=_cross_ratio(z[:,0]);xv=_cross_ratio(z[:,1])
    if not(0<xu<1 and 0<xv<1):raise ValueError('bad cross ratio')
    return -kw.get('N',12)/6*(np.log1p(-xu)+np.log1p(-xv))

def same_static_time_points(radii,t=0.,r_reference=400.,**kw):
    return [(t+tortoise(r,r_reference,**kw),r) for r in radii]

def first_order_entropy_response(points,**kw):
    """Derivative for prescribed geometry/state continuation, NOT full entropy."""
    kb,_=kappa(kw.get('mass',100.),kw.get('L',3e-6),kw.get('G',1.))
    z=_endpoint_state(points,**kw)
    rref=kw.get('r_reference',400.)
    sk={k:v for k,v in kw.items() if k!='r_reference'}
    chi=np.array([coordinate_correction(v,r,rref,**sk) for v,r in points])
    U=z[:,0];dU=-kb*U*chi
    terms=[(0,1,1),(2,3,1),(1,2,1),(0,3,1),(0,2,-1),(1,3,-1)]
    bnd=sum(sign*(dU[i]-dU[j])/(U[i]-U[j]) for i,j,sign in terms)
    drho=sum(linearized_rho(v,r,rref,**sk)[1]+kb*ch/2
              for (v,r),ch in zip(points,chi))
    return kw.get('N',12)/6*(bnd+drho)

def lattice_free_dirac_entropy(indices):
    """Independent half-filled infinite-chain fermion correlator."""
    x=np.asarray(sorted(indices),int)
    D=x[:,None]-x[None,:]
    C=np.empty(D.shape,float);C[D==0]=.5
    sel=D!=0;C[sel]=np.sin(np.pi*D[sel]/2)/(np.pi*D[sel])
    lam=np.clip(np.linalg.eigvalsh(C),1e-15,1-1e-15)
    return float(-np.sum(lam*np.log(lam)+(1-lam)*np.log1p(-lam)))

def lattice_mutual_information(L1,gap,L2):
    a=np.arange(L1);b=np.arange(L1+gap,L1+gap+L2)
    return (lattice_free_dirac_entropy(a)+lattice_free_dirac_entropy(b)-
            lattice_free_dirac_entropy(np.r_[a,b]))

def continuum_mutual_information(L1,gap,L2,c=1.):
    x=L1*L2/((L1+gap)*(L2+gap))
    return -c/3*np.log1p(-x)

if __name__=='__main__':
    p=same_static_time_points((300.,350.,450.,500.))
    print('two interval S',entropy_two_intervals(p),'MI',mutual_information(p))
    print('first-order derivative',first_order_entropy_response(p))
    for n in [10,20,40,80]:
        print('lattice',n,lattice_mutual_information(n,n//2,n),
              'continuum',continuum_mutual_information(n,n//2,n))
