# PI additions — controlled implementation plan

Source request: **Mechanistic Selectivity Roadmap — Computationally Efficient High-Impact Plan**

This branch starts from `mechanism-first`. The existing paper premise is not replaced. Every addition must strengthen the same physical question:

> When a MOF prefers CO2 over N2 at infinite dilution, what roles are played by the repulsive WCA reference, attractive Lennard-Jones interactions, and electrostatics?

The current central result remains:

`repulsive WCA baseline -> attractive creation of CO2 preference -> framework-specific electrostatic differentiation`

The PI roadmap is implemented only where it sharpens this story.

---

## Gate A — existing-data additions first

These require **no new molecular simulation** and are the first priority.

### A1. Continuous mechanism fingerprint
**Request:** define the continuous vector
`m_i = (Delta ln S_attr, Delta ln S_elec)`.

**Scientific question:** after attraction creates CO2 preference, how variable is the balance between attractive and electrostatic enhancement?

**Existing evidence:** the authoritative 238-framework table already contains both increments.

**Work:**
- add mechanism magnitude;
- add mechanism angle;
- add electrostatic-share descriptor only as a descriptive statistic;
- preserve continuous coordinates as primary evidence;
- do not introduce hard mechanism classes in the main paper.

**Placement:** source data + SI; selected summary may enter Figure 2 if useful.

**Changes frozen claim?** no.

**New simulation?** no.

**Decision:** IMPLEMENT.

### A2. Performance twins / mechanism strangers
**Request:** identify pairs with similar final selectivity but large distance in mechanism-vector space.

**Scientific question:** can similar endpoint performance arise from different balances of attraction and electrostatics?

**Selection rule:** freeze the rule before inspecting pair identities. Use the closest 5% of all pairwise differences in `ln S_full`, then rank those pairs by Euclidean distance in `(Delta ln S_attr, Delta ln S_elec)`.

**Important constraint:** this is secondary to the already-established sequential population mechanism. It must not replace the main paper premise unless the result is genuinely population-wide.

**Placement:** Figure 2 or SI depending effect strength.

**Changes frozen claim?** potentially adds a secondary claim, not the central claim.

**New simulation?** no.

**Decision:** IMPLEMENT AND TEST.

### A3. Mechanism degeneracy
**Request:** quantify how common endpoint-matched/mechanism-distant behavior is.

**Work:**
- distribution of mechanism distance among endpoint-matched pairs;
- median and upper quantiles;
- local spread of mechanism angle within narrow `ln S_full` bins;
- no categorical mechanism map unless boundaries are objective and stable.

**Placement:** SI by default; one concise main-text statistic only if the population effect is strong.

**New simulation?** no.

**Decision:** IMPLEMENT.

### A4. Upgrade ranking consequences
**Existing work already completed:**
- WCA -> LJ Spearman;
- LJ -> full Spearman;
- median absolute rank shifts;
- maximum promotions/demotions;
- counts moving >=50 positions.

**Additions still worth doing:**
- Kendall tau;
- fractions moving >10%, >25%, >50% of the 238-position range;
- top-10, top-25, and top-10% overlap across stages.

**Placement:** Figure 5 / SI.

**New simulation?** no.

**Decision:** IMPLEMENT SMALL EXTENSION ONLY.

---

## Gate B — statistical negative controls

### B1. Framework-level bootstrap
Bootstrap 5,000-10,000 resamples for:
- median attractive increment;
- median electrostatic increment;
- fraction positive;
- attraction/electrostatic Spearman correlation;
- ranking statistics;
- geometry-predictability metrics.

**Placement:** SI; confidence intervals may enter the main text.

**Decision:** IMPLEMENT.

### B2. Predictability gap
Define
`G = R2_attr - R2_elec`
using the same out-of-fold prediction framework and obtain an uncertainty interval.

This supports the physical question:
> geometry encodes attraction much more strongly than electrostatic differentiation.

It is a mechanistic control, not an ML benchmark.

**Decision:** IMPLEMENT.

### B3. Outlier/leverage controls
Current evidence already includes top-1% sensitivity and WCA precision exclusions.

Add only if needed:
- top-5% removal;
- removal of largest rank movers.

Do not proliferate variants.

**Decision:** IMPLEMENT MINIMAL EXTENSION.

### B4. Geometry-controlled permutation
The PI roadmap proposes within-geometry shuffling.

This is scientifically acceptable only if the permutation scheme is defined entirely from geometry and does not reuse target-derived grouping. It is **not** the previously rejected circular-permutation analysis.

**Decision:** IMPLEMENT AS A SECONDARY NULL TEST if the geometry clusters/bins can be defined cleanly.

### B5. Grouped holdouts
The current repository does not contain a validated grouping hierarchy sufficient to support chemistry-family or topology leave-one-group-out claims. Previous family-holdout work was explicitly not completed.

**Decision:** DO NOT CLAIM. Revisit only if trustworthy group metadata are recovered and audited.

---

## Gate C — new electrostatic design coordinate

The PI roadmap's strongest proposed new calculation is a zero-density accessible-pore electrostatic probe descriptor across the 238 frameworks.

**Scientific value:** high. It would convert the current result
> geometry does not adequately recover electrostatic variation
into a positive physical statement about what does recover it.

**But:** this is a genuinely new analysis/simulation stage and is not needed to defend the current paper.

**Decision gate:** run only after A and B are complete.

Proceed only if:
1. endpoint-mechanism non-uniqueness is population-wide enough to matter;
2. the attraction/electrostatic predictability asymmetry survives bootstrap/outlier controls;
3. no existing retained calculation can already provide the same descriptor.

If opened, the calculation must:
- use the production electrostatic convention;
- use pre-specified probe/orientation counts;
- freeze descriptor definitions before model comparison;
- compare geometry-only vs electrostatic-probe-only vs combined features on identical splits;
- avoid simplistic maximum-charge or mean-charge descriptors as the main result.

**Decision:** CLOSED — numerical preflight failed after the single permitted 10,000-position retry. Retained as a non-supporting SI diagnostic; no 238-framework production run or downstream probe-model comparison.

---

## Gate D — second-force-field bridge

Existing evidence:
- the 15-MOF DREIDING-linker/UFF-Co bridge is valid only for the WCA coordinate;
- the previous full/LJ bridge cannot support quantitative claims because the retained output layout lost seed-resolved observables.

The PI roadmap asks for a new 8-framework alternative-LJ mechanism bridge.

**Scientific value:** useful robustness if the high-impact story needs it.

**Cost:** modest but real new simulation.

**Rule:** if run, select the panel before seeing alternative-force-field results and preserve all non-LJ choices.

**Decision:** CLOSED by the Gate C hard-stop. Do not run unless explicitly requested by the PI or a reviewer.

---

## Gate E — CO2/CH4 transfer

The roadmap treats this as the final generalization experiment, not part of the core proof.

The current Story Charter excludes CO2/CH4 transfer from the frozen claim set.

**Decision:** CLOSED by the Gate C hard-stop. Reopen only if explicitly requested by the PI or a reviewer.

---

## Figure protection

The current five-figure logic remains the default:

1. **Figure 1:** population-wide nested interaction hierarchy.
2. **Figure 2:** attraction vs electrostatic mechanism map.
3. **Figure 3:** geometry as a mechanistic control.
4. **Figure 4:** M029/M085 matched mechanistic counterexample.
5. **Figure 5:** screening/ranking consequence.

PI additions may enrich these figures, but they may not create a new paper identity.

Potential changes after Gate A/B:
- Figure 2 may gain an objective endpoint-twin overlay if the result is strong.
- Figure 3 may gain a bootstrap/predictability-gap inset.
- Figure 5 may gain top-set-overlap or thresholded rank-movement summaries.

Potential changes after Gate C:
- Figure 3 could become the two-design-coordinate figure only if the electrostatic probe descriptor clearly recovers explanatory power.

Figure 4 remains M029/M085 unless an objective pre-specified selection rule produces a stronger mechanistic exemplar without cherry-picking.

---

## Hard anti-drift rule

The PI roadmap states a broader desired message: similar endpoint selectivity may arise from different mechanism vectors.

That is compatible with the current paper only as a **secondary result** unless the existing-data analysis proves it is widespread.

The central paper claim remains population-level and sequential:

> repulsive WCA baseline -> attractive creation of CO2 preference -> framework-specific electrostatic differentiation.

No new result may replace that premise without an explicit evidence-based decision recorded in `docs/DECISIONS.md`.
