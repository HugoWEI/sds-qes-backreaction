# Phase 1: finite-window Gaussian collector with exact work and purification

**Date: 30 October 2026.** Reproducible calculation; **not a new SdS QES result**. The complete local Phase-0 + Phase-1 snapshot passes **77 tests** (60 previous + 17 new). A current remote GitHub Actions result has not been independently checked.

## Source/prior-art scope

- Easson, Phys. Rev. D **113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk : classical SdS mass loss and Polyakov energy current, not a specified microscopic collector.
- Aalsma, Aguilar-Gutierrez, Sybesma, JHEP **01** (2023) 129, https://doi.org/10.1007/JHEP01(2023)129 : finite-duration de Sitter departure from equilibrium with controlled backreaction **in a different gravity model**; not ours.
- Meineri, Penedones, Rousset, JHEP **02** (2020) 138, https://doi.org/10.1007/JHEP02(2020)138 : energy transmission/reflection across conformal interfaces and nonequilibrium currents.
- Standard Gaussian fermion correlation-matrix and canonical-purification techniques (Casini–Huerta / Peschel and subsequent work): **not a new quantum-information theorem**.

## Exactly defined discrete quantum model

Two independent chains, H (hot) and C (cool collector input), each with an open-boundary identical free hopping Hamiltonian `h_(j,j+1)=h_(j+1,j)=-1/2`; no on-site potential. Initially independent Gaussian Gibbs states

```
C_H=[I+exp(beta_H*h)]^-1,    C_C=[I+exp(beta_C*h)]^-1,
beta_H < beta_C.
```

Define a local pulse-angle array `theta_j` with zero outside a compact *bin window*, sine-squared ramps on its edges and a plateau in between. In a one-shot exact canonical fermion operation,

```
[H_out]   [ diag(cos theta)    diag(sin theta) ] [H_in]
[C_out] = [-diag(sin theta)    diag(cos theta) ] [C_in].
```

At each bin `j`, the transmission is `sin² theta_j`. This is a **regulated time-bin/chain model**. It is not an actual continuum time-dependent conformal defect Hamiltonian, has no physical worldtube action, and is not yet coupled dynamically to gravity.

The output reduced covariances are

```
C_H'=c_i c_j C_H + s_i s_j C_C,
C_C'=s_i s_j C_H + c_i c_j C_C,
C_HC'=-c_i s_j C_H + s_i c_j C_C.
```

Exact `E=Tr(hC)` bookkeeping yields, for a temperature contrast term `Q`,

```
Delta E_C = Q + W_C,
Delta E_H = -Q + W_H,
W_switch = W_H + W_C = Delta E_H + Delta E_C,
Q = Tr[h * (s_i s_j (C_H-C_C)_ij)],
W_K = Tr[h * ((cos(theta_i-theta_j)-1) C_K,ij)].
```

Thus the externally applied *switching operation* supplies work. For the nearest-neighbor real hopping chosen here,

```
W_switch=sum_(j=0)^(n-2) (C_H[j,j+1]+C_C[j,j+1])
         *(1-cos(theta[j+1]-theta[j])) >= 0
```

for these positive-temperature half-filled chains, which have positive nearest-neighbor correlators. A spatially uniform rotation costs no work; a finite ramp costs work. The small-angle work is quadratic in theta gradients.

The total fermion number and **global Gaussian von Neumann entropy** are exactly preserved by the full two-chain unitary.

## Exact Gaussian purifier and information observable

For each thermal input covariance C, a chosen canonical fermion purification on doubled modes is

```
P = [[C, sqrt(C(I-C))],
     [sqrt(C(I-C)), I-C]],  P²=P.
```

Its ancilla (R_H) purifies the hot input. Under a collector rotation, the off-diagonal collector-to-hot-reference covariance is

```
C_(C_out,R_H) = -diag(sin theta)*sqrt(C_H(I-C_H)).
```

For a chosen finite subset of collector bins K, compute exactly

```
I(C_K:R_H) = S(C_K)+S(R_H)-S(C_K union R_H)
S(A) = -Tr[C_A log C_A+(I-C_A)log(I-C_A)].
```

The hot reference is **an arbitrary canonical Gaussian ancilla, NOT an identified SdS black-hole interior/island**, and the MI is not a fine-grained Page entropy by itself.

### Main numerical baseline

`n=96, beta_H=2, beta_C=8, start=28,stop=68,taper=10, peak_angle=pi/4`:

```
thermal contrast energy Q           1.4835351111991573
collector delta energy              1.506875329676788
hot delta energy                   -1.468233842817409
external switch work                0.038641486859380
energy conservation residual        ~3.6e-15
```

For a *temperature ratio* chosen to match Easson's reference classical SdS `G=1,Lambda=3e-6,M=100`, with `beta_C=beta_H*kappa_b/kappa_c=5.542699080750591`, the same pulse yields

```
Q                                1.3332521679960507
W_switch                         0.03782713668218591
Delta E_C                        1.355778036296492
I(collector_window : R_H)       13.037972391975138
```

This temperature-ratio match is only a **dimensionless analogy**. No physical lattice spacing, 2D-to-4D energy conversion, detector worldtube redshift or coupling to the gravity boundary conditions is claimed. For independent N identical Dirac species, extensive energy and entropy terms multiply by N.

**Finite-size control:** At fixed 40-bin pulse window centered in chains of n=80,96,128,160, the above Q, W and I agree numerically to ~1e-12–1e-13 absolute, due to exponential thermal locality.

## Energy flux does NOT determine the fine-grained purifier information

Take **equal** input temperatures `C_H=C_C=C`, and a **uniform** nontrivial rotation `theta=pi/4`. Then exactly

```
Delta E_C=0, Delta E_H=0, W_switch=0,
C_C'=C_C=C, but I(C_out : R_H)>0.
```

For `n=32,beta=2` the collector–hot reference mutual information is `13.374493083232661` nats at `theta=pi/4`, vs 0 with `theta=0`. The two situations have identical collector marginal density matrix and mean energy, but different correlations with a labeled reference. This is an **elementary known general lesson**, made explicit in the free-fermion benchmark: a stress tensor / conserved mean energy flux does not fix the radiation density matrix or its correlations.

## Validation

- Exact full one-particle unitary compared with compact covariance formulas;
- Total number and global entropy invariance;
- Zero pulse, global uniform mixing and full-swap purifier limits;
- Work identity and thermal-versus-switching split to ~1e-12 or better;
- Positive temperature-gradient transfer, small-angle quadratic work and reduced work for smoother ramps;
- Exact `P²=P` purifier, nonnegative mutual information, hot-reference transfer;
- Equal-temperature zero-net-heat positive-hot-reference-information example;
- Finite-size reference-information convergence.

**Local combined regression:** `77 passed`. Run `python -m pytest -q` for entire repo; additional direct scripts are `python src/finite_collector.py`, `python src/purified_collector.py`.

## Unsolved physical tasks — no novelty claims yet

1. Replace the one-shot *time-bin unitary* by a **causal real-time conformal defect/worldtube Hamiltonian**, including derivative/contact energy production for switching and a properly normalized stress tensor.
2. Specify a full 2D gravity–bath interface action, Israel/dilaton boundary conditions, energy balance in the same Killing normalization, and gravitational response to any reflected channel.
3. Choose a globally well-defined UdS state with a **physically justified purifier**. Canonical thermal ancilla is not a black-hole interior.
4. Propagate the free-field covariance on dynamical SdS, with valid horizon patches, renormalized UV and boundary counterterms.
5. Compare **the same radiation algebra** on no-island and behind-horizon two-null island replica saddles; test existence, dominance, and semiclassical control.
6. Search published and unpublished literature for the exact resulting model-specific finite-Lambda result before claiming novelty.

**Conclusion:** This is a verified microscopic collection/energy-information toy model and a useful falsification tool for careless entropy-flux arguments, but it is NOT a physically completed finite-duration SdS collector and proves no QES or Page time.
