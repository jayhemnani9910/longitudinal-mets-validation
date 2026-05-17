# B9 Decision Tree Structure (from Shin et al. 2023 PLoS One Figure 8)

Reference document for `R/scores/b9_tree.R` implementation (Task 1.5 in the implementation plan).

## What we have (verbatim from Figure 8B in B9 paper)

The published Figure 8 has three panels:
- **(A)**: full graphical tree at depth 5, but rendered at low resolution. Internal split values not legible in the published PDF.
- **(B)**: text representation of the SUBTREE rooted at `BPWC_add > 0.66` only (the right subtree). Transcribed below.
- **(C)**: 2D MetS risk-map projection on the (WC, BP) plane.

### Right subtree (BPWC_add > 0.66) — extracted verbatim

```
BPWC_add > 0.66
├── BPWC_mul <= 0.31
│   ├── BPWC_dif <= 0.10
│   │   ├── BPWC_add <= 0.84
│   │   │   ├── BPWC_dif <= -0.25 → class 1
│   │   │   └── BPWC_dif >  -0.25 → class 0
│   │   └── BPWC_add > 0.84
│   │       ├── BPWC_mul <= 0.24 → class 1
│   │       └── BPWC_mul >  0.24 → class 1
│   └── BPWC_dif > 0.10
│       ├── BPWC_mul <= 0.17
│       │   ├── BPWC_dif <= 0.56 → class 0
│       │   └── BPWC_dif >  0.56 → class 0
│       └── BPWC_mul > 0.17
│           ├── BPWC_mul <= 0.22 → class 1
│           └── BPWC_mul >  0.22 → class 1
└── BPWC_mul > 0.31  → class 1 (terminal; subtree continues but consistently labels MetS)
```

### What's missing

- **Left subtree (BPWC_add <= 0.66)**: the "safety zone" the paper describes. Most leaves likely class 0. Splits at this depth are not visible in published Figure 8(A).
- **Leaf probabilities**: Figure 8(A) is colored by class (orange = non-MetS, blue = MetS, darker = higher probability) but the exact probability per leaf is not numerically labeled. The paper instead reports calibrated risk via Pozzolo's method with adjusted threshold = 0.137.

### Calibration

After tree training, B9 applies Pozzolo's calibration:
- Pre-calibration tree outputs are biased high because training used downsampling (1:1 MetS:non-MetS), but population prevalence is 13.6%.
- Pozzolo's correction: `p' = (β·p_s) / (β·p_s − p_s + 1)` where β = downsample selection probability and p_s = uncorrected prediction probability.
- Adjusted threshold for binary diagnosis: 0.137 (matches population prevalence).

### Worked example to verify implementation (Figure 7)

Input: 55-yo female, sBP=140, dBP=90, waist=89cm.
Intermediate: BP=0.84, WC=0.66.
Synthetic: BPWC_add=1.50, BPWC_mul=0.55, BPWC_dif=0.18.
Path: BPWC_add (1.50) > 0.66 → BPWC_mul (0.55) > 0.31 → class 1.
Expected calibrated probability: 0.31.
Risk ratio: 2.25× the population baseline.

This worked example confirms that subjects with high BPWC_mul (> 0.31) reach a terminal leaf early in the right subtree with calibrated probability ≈ 0.31.

## Implementation strategy for Task 1.5

### Option A (recommended): Faithful right-subtree + safety-zone marginal

For `BPWC_add > 0.66`: use the transcribed right-subtree splits verbatim. Assign calibrated leaf probabilities by stratifying NHANES test cohort by leaf and computing observed MetS rate (or use the paper's 0.31 worked-example value as a baseline for the matched leaf).

For `BPWC_add <= 0.66`: lump all such subjects into a single "safety zone" leaf with calibrated probability ≈ population baseline (use Pozzolo's β·prior / [β·prior − prior + 1] = 0.137 if NHANES population MetS prevalence matches).

Flag this in the manuscript as a faithfulness limitation. Address it via Option B as a sensitivity analysis.

### Option B (sensitivity): Refit a depth-5 CART on NHANES

This is plan task 5.2b (transportability test). Train a depth-5 CART with the same hyperparameters on NHANES baseline data using the same synthetic features. Compare tree structure to what's visible in Figure 8(A). Report this as a transportability analysis, not a faithful external validation.

### Option C (fallback): Email Shin/Oh for the trained tree object

If their published code includes the trained `rpart` object, we can use that directly. Add to the reachout email at week 7-8: "Would you be willing to share the trained decision tree model from B9 for an exact external validation?"

## Implementer's task

When implementing `R/scores/b9_tree.R`:
1. Encode the right-subtree splits exactly per the verbatim transcription above.
2. Encode the left subtree as a single safety-zone leaf with `p = 0.137`.
3. Pass the unit tests in `tests/scores/test-b9-tree.R`.
4. Add an integration test using the Figure 7 worked example (Input: sBP=140, dBP=90, waist=89, female. Expected: BPWC_add=1.50, BPWC_mul=0.55, BPWC_dif=0.18, leaf risk ≈ 0.31).
5. Document Option A's limitation in the manuscript's "Methods" section.
