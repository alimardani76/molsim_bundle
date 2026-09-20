# PI Gate B — balanced geometry-conditioned shuffle diagnostic

Status: REVIEW

## Why the original K-means version was rejected

The first implementation produced singleton geometry clusters. Rather than
merge those outliers after seeing the result, this version uses deterministic
equal-frequency strata in the first two principal components of the same
14 frozen geometry descriptors.

This keeps the diagnostic geometry-conditioned while ensuring every stratum
contains enough frameworks to shuffle.

## Geometry stratification

PCA explained-variance ratios:
- PC1: 0.6512
- PC2: 0.2774
- PC1+PC2: 0.9287

Sensitivity designs:
- K=8: 4 PC1 bins x 2 within-bin PC2 bins
- K=10: 5 x 2
- K=12: 4 x 3

Each target is shuffled 10,000 times within these geometry-only strata.

## Primary K=10

### Attractive increment

Observed Spearman(prediction, target):
**0.9239**

Conditional-shuffle 95% interval:
**[0.8698,
0.8937]**

Observed R2:
**0.8745**

Shuffle-null R2 95% interval:
**[0.7495,
0.8032]**

### Electrostatic increment

Observed Spearman(prediction, target):
**0.2929**

Conditional-shuffle 95% interval:
**[-0.0891,
0.1541]**

Observed R2:
**0.1880**

Shuffle-null R2 95% interval:
**[-0.1547,
0.0914]**

## Sensitivity contrast

- K=8: attraction excess-over-null=0.024; electrostatic excess-over-null=0.135; contrast=-0.111
- K=10: attraction excess-over-null=0.030; electrostatic excess-over-null=0.139; contrast=-0.109
- K=12: attraction excess-over-null=0.022; electrostatic excess-over-null=0.117; contrast=-0.095

## Interpretation boundary

This is a secondary shuffle diagnostic, not an exact conditional-randomization
p-value. The OOF geometry score was learned on the original dataset and is held
fixed during the shuffles.

Use it only if the same qualitative attraction-vs-electrostatics asymmetry
appears for K=8, K=10, and K=12.

The main claim remains supported independently by:
1. the observed OOF predictability gap;
2. the framework bootstrap of that gap;
3. outlier/leverage sensitivity;
4. geometry-matched mechanistic counterexamples.

Do not use this diagnostic to claim that electrostatics is independent of
geometry or that chemistry has been uniquely isolated.
