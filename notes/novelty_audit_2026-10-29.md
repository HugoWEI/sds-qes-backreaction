# Targeted novelty audit — 29 October 2026

## Verdict (do not overstate)

**The specific question is plausibly open but originality of a result is NOT established.** No physical two-null QES in a consistent time-dependent Polyakov-backreacted neutral Schwarzschild–de Sitter (SdS) spacetime has been computed in this repository. Exact model-specific QES locations and a fine-grained Page-time correction would be potentially original only after checking both the prior universal result and the validity of a fully defined gravitating-radiation system.

### What Easson actually leaves open

Damien A. Easson, *Fate of Schwarzschild–de Sitter black holes: Nonequilibrium evaporation*, **PRD 113**, 084014 (2026), https://doi.org/10.1103/trlf-d6sk ; arXiv:2511.11873.

Section I explicitly says *“A detailed QES extremization in the fully backreacted SdS geometry is left for future work.”* Section IX gives the two-null QES prescription and a **thermo-controlled Page proxy**, but it does not solve explicit replica saddles for a fully specified radiation algebra. This is strong evidence of a gap in **that paper**, not a proof the wider literature leaves it untouched. His model is adiabatic; one must not claim a globally exact time-dependent metric based only on the instantaneous mass-loss law.

## Direct prior art and overlap

| Reference | Result already in the literature | Our required differentiation |
| --- | --- | --- |
| Goswami & Narayan, JHEP **10** (2022) 031, https://doi.org/10.1007/JHEP10(2022)031 | Frozen SdS static patch, growing radiation entropy, late island and Page curve for small BH | Add explicitly derived, controlled anomaly-sourced dynamical geometry; no frozen-background rephrasing |
| Yadav & Joshi, PRD **107** (2023) 026009, https://doi.org/10.1103/PhysRevD.107.026009 | Black-hole- and cosmological-side SdS islands and Page curves, using thermal opaque membranes and frozen opposite patch | No claim that two horizon islands alone are new |
| Goswami & Narayan, JHEP **05** (2024) 016, https://doi.org/10.1007/JHEP05(2024)016 | Future-boundary radiation collection and BH islands just inside horizon | Specify a distinct and operationally physical algebra / time-dependent correction |
| Balasubramanian, Kar & Ugajin, JHEP **02** (2021) 072, https://doi.org/10.1007/JHEP02(2021)072 | Gravitating de Sitter JT coupled to non-gravitating degrees of freedom, replica wormholes, both horizons, backreaction | Full spherical SdS versus de Sitter JT, state/flux and finite SdS horizon physics |
| Aalsma & Sybesma, JHEP **05** (2021) 291, https://doi.org/10.1007/JHEP05(2021)291 | Collecting de Sitter radiation drives backreaction/singularity as an island appears | Need a controlled backreaction-size/validity diagnostic, not merely a formal QES |
| Kames-King, Verheijden & Verlinde, JHEP **03** (2022) 040, https://doi.org/10.1007/JHEP03(2022)040 | Unruh-dS observer collection and auxiliary bath, found catastrophic backreaction by Page time in their 3D/partial-reduction setup | Such an obstruction might occur, but CANNOT be assumed to hold automatically in 4D-spherical SdS |
| Matsuo, JHEP **03** (2025) 068, https://doi.org/10.1007/JHEP03(2025)068 | General evaporating Unruh-state QES inside horizon vs equilibrium Hartle-Hawking outside, independent of BH details | Leading near-horizon location/sign alone is not new; isolate any novel two-horizon finite-Lambda correction |
| Đorđević & Radovanović, PRD **113** (2026) 105020, https://doi.org/10.1103/nc2z-hx6s | Polyakov-Liouville backreaction + fine-grained Page curve for evaporating ordinary Schwarzschild | Demonstrate finite cosmological-horizon dependence absent at Lambda=0 |
| Casini & Huerta, JHEP **03** (2009) 048, https://doi.org/10.1088/1126-6708/2009/03/048 ; Class Quantum Grav **26** (2009) 185005, https://doi.org/10.1088/0264-9381/26/18/185005 | Exact free 1+1 massless Dirac multi-interval entropy/modular density matrix | Free-Dirac formulas and lattice MI verification are useful known benchmarks, not original results |

**Search coverage:** APS, Springer/JHEP, arXiv-related records and title/keyword searches for 2025–October 2026 SdS Polyakov/QES/backreaction/islands, plus specific older state/bath work. Search is not a systematic bibliographic proof of absence; unpublished, very recent and differently titled results could exist.

## Audit of existing project artifacts

- Exact SdS thermodynamic-proxy Page-mass identity: algebraic simplification of Easson’s first laws; **not a new fine-grained Page result**.
- EF radial constraints, two EF charts and coordinate integrability: model-specific implementation of established Einstein-dilaton, stress-conservation and coordinate-transform identities; not a stand-alone novel theorem.
- Chiral entropy and adiabatic peeling derivative: standard 2D CFT interval transformation in an explicitly postulated evolving conformal state; not a derived physical quantum state.
- Free Dirac two-interval entropy, smooth conformal bath collar, lattice correlation-matrix comparison: known Gaussian free-field machinery and an **auxiliary** bath model, not a gravity–bath junction solution.
- Exterior near-boundary roots below the UV cutoff: negative controls; **not physical islands**. No behind-horizon two-null saddle or Page crossing has been established.
- Generic first-order root or Page shift via implicit-function/envelope theorem: standard perturbation theory. The novelty, if any, must be in the explicitly evaluated, physically well-defined **coefficient**, its sign/parameter dependence/phase diagram, or a nonexistence/validity theorem.

## Proposed narrow testable claims (NOT YET PROVED)

### Option A — Quantitative two-horizon correction

Specify one matter theory (e.g. free Dirac with fixed spin structure), initial purification, physical detector/bath coupling, gravitational junction conditions and renormalization scheme. Construct both semiclassical horizon patches including correct mass/redshift integration constants, match/continue across horizons, compute the same radiation algebra for both saddle branches and solve both null QES conditions. Isolate the finite-`Lambda` anomaly-induced change to the **on-shell** saddle difference, e.g.
```
D(t,epsilon)=S_no-island(R;t,epsilon)-
             S_island(R;t,epsilon)|QES,
delta t_Page=-epsilon*(partial_epsilon D)/(partial_t D) at crossing.
```
The identity is universal; only its evaluated, cutoff/gauge-independent finite-`Lambda` term might be new. Check against Matsuo 2025 and Schwarzschild Polyakov 2026 in the appropriate limit.

### Option B — Validity/obstruction theorem (potentially more robust)

For the same radiation algebra and physically valid geometry, compare candidate island-saddle formation and saddle-crossing times to the onset of order-one backreaction/curvature or nonuniform horizon boundary-layer effects. A proof that no **controlled** fine-grained Page transition exists in a precisely delimited Easson-model parameter regime could be important, but comparable de Sitter no-Page-curve arguments exist (Aalsma/Sybesma 2021, Kames-King et al. 2022). Therefore any result must state precisely what SdS-specific assumptions/parameter region are new.

## Minimum claim checklist

- [ ] Matter state and reduced density matrix defined; not inferred solely from stress flux
- [ ] Physical radiation algebra and purifier defined on a Cauchy slice
- [ ] Junction energy conservation and gravitational boundary/counterterms fixed
- [ ] Both horizon boundary layers controlled; no use of singular EF overlap beyond validity
- [ ] At least one **behind-horizon**, UV-separated saddle solved in **both** null directions, or a rigorous nonexistence theorem
- [ ] Same-algebra no-island and island entropies computed with renormalized constants
- [ ] Full numerical equations compared with first-order expansion
- [ ] Cutoff, gauge, matter-species, and small-Lambda sensitivity checked
- [ ] Compare directly with Matsuo 2025 and the de Sitter backreaction obstructions
- [ ] Verify specific originality claim with current bibliography again before public preprint

**Safe external statement:** “Developed and numerically validated a reproducible semiclassical SdS Polyakov model and free-field entropy benchmarks; investigating whether a controlled two-horizon backreaction correction to QESs or an SdS-specific obstruction exists.” **Not safe:** “Discovered a new QES,” “proved a novel Page-time shift,” or “first SdS islands.”