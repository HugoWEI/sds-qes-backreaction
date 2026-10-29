# Literature re-check of the controlled-Page-transition formulation (2026-10-29)

**Question tested:** For a specified radiation algebra in the neutral Unruh–de Sitter (UdS), Polyakov-anomaly-backreacted, spherically reduced Schwarzschild–de Sitter model, when does a UV-separated behind-horizon QES *exist, dominate and remain semiclassically valid through a fine-grained Page transition*?

**Verdict:** The individual notions and several partial answers are old, but this exact conjunction for Easson's time-dependent neutral SdS model is **not answered in the papers located**. This is **not proof of novelty**, and no new result has been obtained by merely posing it.

## Crucial directly relevant papers

1. **Easson**, PRD 113, 084014 (2026), https://doi.org/10.1103/trlf-d6sk — two-horizon anomaly stress/flux and thermodynamic Page *proxy*, microscopic full QES left open. No defined collector/purifier/renormalized competing saddles in the proxy.
2. **Yu & Ge**, *Geometric constraints via Page curves: insights from island rule and quantum focusing conjecture*, Chinese Phys. C 49, 045107 (2025), https://doi.org/10.1088/1674-1137/ada34e — near-horizon geometric existence condition for **static** spherically symmetric spacetimes, including (A)dS. Their sufficient criterion is `f''(rh)<0`, which automatically holds for classical neutral SdS: `f''(r)=-4GM/r^3-2Lambda/3<0` at positive radius. **This does NOT establish the far-from-equilibrium, observable Page transition or validity of the full time-dependent model**. Do not claim that static curvature implies UdS Page dominance.
3. **Ageev, Aref'eva, Belokon, Pushkarev & Rusalev**, *Entanglement entropy in de Sitter: no pure states for conformal matter*, JHEP05(2024)308, https://doi.org/10.1007/JHEP05(2024)308 — de Sitter region/complement purity and CFT formula pitfalls; examples with entropy growth and no island, as well as an island that never dominates. **Existence != dominance** is known.
4. **Akhmedov & Bazarov**, *Backreaction issue for the black hole in de Sitter spacetime*, PRD 107, 105012 (2023), https://doi.org/10.1103/PhysRevD.107.105012 — generic thermal states of scalar fields on SdS have problematic horizon stress/backreaction. Their thermal state and matter theory are not the specified outgoing UdS free-Dirac state.
5. **Kames-King, Verheijden & Verlinde**, *No Page curves for the de Sitter horizon*, JHEP03(2022)040, https://doi.org/10.1007/JHEP03(2022)040 — catastrophic collection-induced backreaction at the Page time in partially reduced 3D dS; not a theorem against every neutral 4D-spherically reduced SdS black-hole-radiation algebra.
6. **Aalsma & Sybesma**, *The price of curiosity*, JHEP05(2021)291, https://doi.org/10.1007/JHEP05(2021)291 — a collector in dS JT can develop a singularity at island appearance, before Page time. Model-sensitive.
7. **Goswami & Narayan**, JHEP10(2022)031 and JHEP05(2024)016 — frozen-background SdS black-hole islands, different radiation regions and cosmological horizon treatments.
8. **Yadav & Joshi**, *Cosmological and black hole islands in multi-event horizon spacetimes*, PRD 107, 026009 (2023), https://doi.org/10.1103/PhysRevD.107.026009 — two-horizon SdS islands and Page curves using opaque membranes; opposite horizon treated frozen.
9. **Matsuo**, *Universal structure of islands in evaporating black holes*, JHEP03(2025)068, https://doi.org/10.1007/JHEP03(2025)068 — generic Unruh versus Hartle-Hawking QES horizon placement. Near-horizon position alone is not new.
10. **Shi, Turiaci & Wu**, *The fate of nucleated black holes in de Sitter quantum gravity*, JHEP08(2026)191, https://doi.org/10.1007/JHEP08(2026)191 — 2026 state regularity, nucleation and smooth semiclassical evaporation; does NOT calculate our fine-grained collector Page transition. Means "de Sitter backreaction inevitably prohibits physical evolution" would be too broad.
11. **Chowdhury, Paul & Gangopadhyay**, *New insights on mutual information in the island approach to the Page curve*, Eur Phys J C 86, 1031 (Sep 2026), https://doi.org/10.1140/epjc/s10052-026-16291-x — consistency/MI of island and radiation regions on Cauchy slices, including prior (A)dS discussion; not the time-dependent Polyakov SdS entropy with specified bath.
12. **Balasubramanian, Kar & Ugajin**, JHEP02(2021)072 — backreacted de Sitter JT coupled to a nongravitating system with replicas; not the complete SdS geometry.

## Corrected formulation: four-way classification, not just yes/no

Pick **one** complete microscopic model `M = (S_grav, free Dirac state, observer/collector coupling and its stress, purifier, renormalization, boundary/junction conditions)`. Define dimensionless SdS parameters `nu=3GM sqrt(Lambda)`, `gamma = cG/r_b^2`, plus physically fixed interaction/collector controls `tau` if an actual transparent/partially transmitting interface exists (not simply a free added knob).

For a fixed radiation von Neumann algebra `R_t`, evaluate:

- **E(t):** there exists an interior QES `q(t)` with `partial_u S_gen = partial_v S_gen=0` and endpoints separated from UV and horizon boundary-layer artifacts in a regular coordinate construction.
- **D(t):** the *on-shell* saddle entropies for the **same** `R_t` obey `S_island(t)<S_no-island(t)`; distinguish a mere QES from a dominant saddle.
- **C(t):** control is maintained by the chosen state, conservation and junction equations, a well-defined UV/gauge-invariant renormalized entropy difference, `cG/r_b^2 << 1`, `|dot kappa_b|/kappa_b^2 <<1`, and finite regular-frame curvature/stress corrections. Near Nariai the adiabatic expansion can be nonuniform and needs explicit special analysis.
- **P:** a crossing `t_P>0` has `S_island(t_P)=S_no(t_P)`, isolated simple zero of `D_{⁠diff}=S_no-S_island`, the no-island branch dominates before and island branch after, and C holds through `t_P`.

**Labels**: (i) controlled Page transition; (ii) QES exists but never dominates; (iii) possible crossing but prior semiclassical breakdown; (iv) no UV-controlled QES. These outcomes are different. A theorem about (iii) is not a theorem that no microscopic Page transition exists; it may only be outside semiclassical reach.

**Specific testable problem:** Determine whether an open subset `Omega_yes` of the `(nu,gamma,tau)` parameter domain satisfies E,D,C,P, and characterize its boundary against other labels. For a well-defined no-external-bath algebra omit tau, rather than use a fictitious collector. Compare `Lambda -> 0` against ordinary Schwarzschild and `nu -> 1` against near-Nariai results. Preserve fixed state, boundary and renormalization scheme across both branches.

## Main new dangers found by this re-check

- **Static sufficient curvature criterion already covers classical SdS**. Do not market the mere sign of `f''` or near-horizon existence as a novel theorem.
- **Existence != dominance** was studied directly by Ageev et al.
- **State/regularity and actual observer backreaction** are not optional; thermal state results differ from prescribed UdS results.
- **Euclidean/no-boundary vs UdS** are different physical initial state choices. 2026 Shi–Turiaci–Wu provides an important smooth-evaporation comparison but no collector-Page theorem.
- **Two EF outer expansions plus a smooth conformal bath collar do not create a physical gravity–bath junction or control the inner-horizon QES.** Existing repo still has no physical saddle, cutoff-independent crossing, or novelty proof.

## Research verdict

The **precisely defined existence+dominance+semiclassical-control phase boundary** for Easson's actual time-dependent neutral UdS SdS + a specified radiation algebra is **not located in this targeted literature search**. That is a **plausible narrowly scoped research opportunity**, contingent on a complete model and a real theorem/calculation. An uncontrolled "island yes/no" would not be original.

The search is not an exhaustive citation-graph audit; no proof of absence from all published/unpublished literature.
