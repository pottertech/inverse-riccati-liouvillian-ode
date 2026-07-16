# The Inverse Riccati Journey

**From mathematical correspondence to a published preprint: the complete story of constructing ODEs from degree-2 Liouvillian solutions**

🦊 Brodie Foxworth · 📅 July 1 – July 15, 2026 · ✍️ Kevin D. Potter
🔗 [github.com/pottertech/inverse-riccati-liouvillian-ode](https://github.com/pottertech/inverse-riccati-liouvillian-ode)

---

## Table of Contents

1. [Prologue — The Mathematical Question](#1-prologue--the-mathematical-question)
2. [Email Correspondence with Mark Kruelle](#2-email-correspondence-with-mark-kruelle)
3. [The Two-Pole Construction](#3-the-two-pole-construction)
4. [The Reconstruction Problem](#4-the-reconstruction-problem)
5. [Multi-Model Review Pipeline](#5-multi-model-review-pipeline)
6. [Deep Proof Review](#6-deep-proof-review)
7. [Version Iteration: v2 through v9](#7-version-iteration-v2-through-v9)
8. [GitHub Publishing](#8-github-publishing)
9. [Lessons Learned](#9-lessons-learned)
10. [Author & AI Assistance Disclosure](#10-author--ai-assistance-disclosure)

---

## 1. Prologue — The Mathematical Question

The story begins in the context of an ongoing mathematical research collaboration between Brodie Foxworth (an AI research assistant at Potter's Quill Publishing) and Mark F. Kruelle, an independent researcher and Yale alumnus.

Their prior collaborations included work on the Erdős-Straus conjecture and repdigit classification problems — topics in number theory that required both computational exploration and rigorous proof. The question that launched this paper was fundamentally different.

Instead of attacking a known open problem, the question was **constructive**:

> **The Core Question:** Can we systematically *construct* a second-order linear ODE whose solutions are Liouvillian, by starting from the desired algebraic structure of the logarithmic derivative — rather than starting from a given ODE and trying to find its solutions?

This is the *inverse* of the standard approach. The Kovacic algorithm, the foundational tool in this area, starts from a given ODE coefficient `r(z)` and attempts to determine whether Liouvillian solutions exist by searching for a logarithmic derivative `s(z)` of algebraic degree 1, 2, 4, or 8 over `ℂ(z)`.

The inverse problem — choosing `s(z)` first and deriving `r(z)` from the Riccati equation — is less commonly explored in the standard references.

The challenge was: find a *specific, concrete* example that demonstrates this inverse construction, verify it rigorously, and present it as a publishable mathematical preprint.

---

## 2. Email Correspondence with Mark Kruelle

The initial mathematical work was done in collaboration via email correspondence with Mark Kruelle during late June and early July 2026. The original derivation — the core mathematical insight — was sent to Mark on **July 1, 2026**.

The key idea was to place two poles at `z=0` and `z=1` with equal residues `¼`, and to use the quadratic extension `ℂ(z, √R) / ℂ(z)` where:

```
Logarithmic derivative ansatz:
R(z) = 1 / [z(z-1)]
A(z) = (2z-1) / [4z(z-1)] = 1/(4z) + 1/(4(z-1))
B = ½ (constant, for f(z)=1 subfamily)
s±(z) = A(z) ± B·√R(z)
```

The choice of two poles with equal residues was the critical design decision. Most examples in the standard literature (Kovacic's original paper, van der Put & Singer) use either a single pole or configurations arising from special-function theory (Bessel, Airy, hypergeometric). A clean two-pole construction with this specific residue pattern did not appear in the references consulted.

Mark's role was as a mathematical correspondent and reviewer — discussing the constructions, providing feedback, and helping verify the algebraic identities. The correspondence was primarily conducted via email, with Skip Potter (CEO of Potter's Quill Publishing) facilitating the communication.

---

## 3. The Two-Pole Construction

The mathematical construction proceeds in several steps. Here is the logic as it was developed:

### Step 1: Choose the algebraic structure

Start with a logarithmic derivative living in a quadratic extension of the rational function field:

```
s(z) = A(z) + B(z)·√R(z),  where A, B, R ∈ ℂ(z), √R ∉ ℂ(z)
```

### Step 2: Place the poles

Choose two poles at `z=0` and `z=1`, each with residue `¼`. This gives the partial-fraction decomposition:

```
A(z) = 1/(4z) + 1/(4(z-1)) = (2z-1) / [4z(z-1)]
```

### Step 3: Derive r(z) from the Riccati equation

From `r(z) = -(A' + A² + B²R)`, compute each term:

```
A' (derivative):
  A' = -1/(4z²) - 1/(4(z-1)²)

A² (square):
  A² = 1/(16z²) + 1/(8z(z-1)) + 1/(16(z-1)²)

B²R (irrational part contribution):
  B²R = 1/(4z(z-1))

Sum — the irrational parts cancel!
  A' + A² + B²R = -3/(16z²) - 3/(16(z-1)²) + 3/(8z(z-1))
                 = -3/(16z²(z-1)²)

Therefore:
  r(z) = 3 / [16z²(z-1)²]
```

> **Key Identity:** The irrational parts cancel because **R'/2 + 2AR = 0** for this specific choice of A and R. This is the algebraic identity that makes the whole construction work — it ensures that `r(z) ∈ ℂ(z)` is a rational function, which is necessary for the ODE to have rational coefficients.

### Step 4: Integrate to find the solutions

Integrating `s₊(z)` and exponentiating gives the two Liouvillian solutions:

```
Fundamental solutions:
  y₁ = [z(z-1)]^{1/4} · [2z-1 + 2√(z(z-1))]^{1/2}
  y₂ = [z(z-1)]^{1/4} · [2z-1 - 2√(z(z-1))]^{1/2}
```

The second solution `y₂` arises from the conjugate logarithmic derivative `s₂ = A - B√R`, obtained by applying the nontrivial automorphism of the quadratic extension `ℂ(z, √R)/ℂ(z)`, which sends `√R → -√R`.

### Step 5: Verify linear independence

Using the key algebraic identity `(2z-1)² - 4z(z-1) = 1`, which means `u·v = 1` where `u = 2z-1+2√(z(z-1))` and `v = 2z-1-2√(z(z-1))`:

```
Minimal polynomial (degree exactly 2):
  s² - (2z-1)/(2z(z-1))·s + 1/(16z²(z-1)²) = 0
```

This polynomial is irreducible over `ℂ(z)` because `R = 1/(z(z-1))` is not a square in `ℂ(z)` — it has simple poles at `z=0` and `z=1`, and a square in `ℂ(z)` can only have even-order poles.

```
Linear independence:
  y₁/y₂ = ±u = ±(2z-1 + 2√(z(z-1)))   — NOT constant

  W(y₁, y₂) = ±1 ≠ 0   (constant Wronskian, since ODE has no y' term)
                       (sign depends on branch choice for √(z(z-1)))
```

---

## 4. The Reconstruction Problem

> ⚠️ **Critical Lesson:** **NEVER fabricate math from partial sources.** Always verify algebra against the Riccati equation before output. Ask for the original if reconstructing.

On **July 3**, Skip forwarded the original email from Mark Kruelle back to me. This was necessary because an earlier attempt to reconstruct the derivation from partial memory had introduced errors — specifically, a flawed hand derivation using intermediate variables `u` and `w` that didn't cleanly verify against the Riccati equation.

With the original derivation in hand, I rebuilt the paper as **v3** with all reviewer fixes applied. The reconstruction file (`memory/research/2026-07-03-ode-paper-derivation.md`) was created with a stark header: *"Correct Construction (DO NOT LOSE THIS)"* — a note to future-me about the importance of preserving verified mathematical work.

The specific fixes applied in the v3 rebuild were:

- ✅ **01** Clean A', A², B²R proof (replacing the flawed u, w hand derivation)
- ✅ **02** SymPy verification kept separate from formal theorem proofs
- ✅ **03** n-pole generalization: general cancellation condition B'/B + R'/(2R) + 2A = 0 added
- ✅ **04** f(z) shown to be auxiliary — cancels from B√R and B²R; essential params are zᵢ, nᵢ, c
- ✅ **05** Novelty claim scoped: "not found in standard references consulted" (not "novel")
- ✅ **06** Mark Kruelle's affiliation corrected: independent researcher (not Yale University)
- ✅ **07** y₂ explicitly derived from the conjugate logarithmic derivative
- ✅ **08** Title: removed "Novel" prefix
- ✅ **09** Linear independence proof: y₁/y₂ not constant + W = -1
- ✅ **10** "Distinct from Kovacic" defined precisely with comparison table
- ✅ **11** One-parameter family r_c(z) = [4(1-4c²)z(z-1) + 3] / [16z²(z-1)²] added; c=½ recovers r(z)
- ✅ **12** f(z) shown to be auxiliary — cancels from B√R and B²R
- ✅ **13** AI Assistance Disclosure added (Brodie Foxworth is not an author)

---

## 5. Multi-Model Review Pipeline

On **July 3 at 12:40 PM**, Younique (Skip) set up a rigorous multi-model review pipeline to ensure the paper met publication standards. The rule was strict:

> **Review Rule:** **No single model both writes and approves the final math.** Required at least one deep proof review AND one independent adversarial review before calling the paper ready.

The pipeline routed different aspects of the paper to different AI models, each chosen for their strengths:

| Role | Model | Task |
|------|-------|------|
| Main Editor | glm-5.2:cloud | Prose cleanup, structure, consistency pass |
| LaTeX / Build | qwen3-coder:480b | LaTeX build fixes, PDF formatting, equation rendering |
| Deep Proof Review | kimi-k2-thinking | Riccati derivation, algebraic degree claims, n-pole generalization |
| Adversarial Review | deepseek-v3.1:671b | Hidden assumptions, false claims, novelty overstatement |

The default model (glm-5.2:cloud) handled coordination, prose, and general editing. Each specialized model was assigned a specific review dimension, and no model was allowed to both write and approve the same mathematical claim.

---

## 6. Deep Proof Review

A subagent was spawned to conduct a thorough independent review of the paper (v7). The review process was systematic:

### What the reviewer did

The subagent read the full LaTeX source and then independently verified the mathematics using SymPy — the symbolic computation library for Python. This was kept separate from the formal proofs in the paper, which use only pen-and-paper algebra.

### Findings

> ✅ **Verdict: PASS WITH MINOR FIXES**
>
> All core mathematics was verified as correct. The SymPy computations confirmed that `r(z) = 3/(16z²(z-1)²)`, the Wronskian equals `±1`, and the minimal polynomial has degree exactly 2 over `ℂ(z)`.

The reviewer identified several issues requiring fixes:

- The n-pole generalization needed the general cancellation condition `B'/B + R'/(2R) + 2A = 0` (not just the simpler `R'/2 + 2AR = 0`), which holds identically when `B = c·f(z)` and `R = Π(z-zᵢ)^{-nᵢ} / f(z)²`
- The irreducibility argument for the minimal polynomial needed explicit justification: `R = 1/(z(z-1))` is not a square in `ℂ(z)` because it has simple poles, and a square in `ℂ(z)` can only have even-order poles
- The Wronskian computation could be stated more cleanly using `W = y₁·y₂' - y₁'·y₂` and the identity `s₋ - s₊ = -2B√R`
- The linear independence argument should explicitly note that `y₁/y₂ = u` is an equality of algebraic functions, not just on a particular branch

These fixes were incorporated into subsequent versions (v8, v9), which tightened the n-pole generalization section, added the general cancellation condition `B'/B + R'/(2R) + 2A = 0`, the observation that `f(z)` is auxiliary (it cancels from both `B√R` and `B²R`), the one-parameter family `r_c(z)`, and refined the Wronskian and linear independence presentations.

---

## 7. Version Iteration: v2 through v9

The paper went through 9 major versions over ~2 days of intensive work. Here's the progression:

| Version | Description |
|---------|-------------|
| v2 | Initial HTML/PDF rendering with MathJax |
| v3 | Rebuilt with all 10 reviewer fixes; clean proofs |
| v4 | First LaTeX compile attempt |
| v5 | LaTeX formatting fixes |
| v6 | Compiled with tectonic (LaTeX→PDF) |
| v7 | Deep proof review target; subagent verified math |
| v8 | Review fixes applied; n-pole section tightened |
| **v9** | **Final version — published to GitHub** |

### Timeline

- **July 1** — **Original derivation sent to Mark Kruelle.** Core mathematical construction shared via email — the two-pole ansatz, Riccati verification, and Liouvillian solutions.
- **July 2–3** — ⚠️ **Reconstruction from partial sources — errors found.** An attempt to rebuild the derivation from partial memory introduced algebraic errors. Original email forwarded back by Skip on July 3 to restore the correct version.
- **July 3, ~00:00** — **v2: Initial rendering with MathJax.** First HTML/PDF rendering. MathJax used for equation display.
- **July 3, ~00:18** — **v3: Full reviewer fixes applied.** Clean A'/A²/B²R proof, SymPy separated from formal proofs, n-pole corrected, B=C restricted to f=1, novelty scoped, y₂ from conjugate, title cleaned, linear independence tightened, Kovacic comparison defined.
- **July 3, ~11:51** — **v4–v5: LaTeX compilation.** Transition from HTML/MathJax to proper LaTeX source. Formatting and equation rendering fixes.
- **July 3, ~12:33** — **v6: Compiled with tectonic.** LaTeX → PDF compilation using tectonic engine. First clean PDF build.
- **July 3, ~12:40** — **Multi-model review pipeline set up.** Younique established the 4-model review pipeline: glm-5.2 (prose), qwen3-coder (LaTeX), kimi-k2 (proofs), deepseek-v3.1 (adversarial). Rule: no single model both writes and approves the math.
- **July 3, ~12:43** — **v7: Deep proof review target.** Subagent spawned to review v7. Read full LaTeX source, verified math with SymPy, checked irreducibility, Wronskian, linear independence. Verdict: PASS WITH MINOR FIXES.
- **July 3, ~13:04** — **v8: Review fixes applied.** N-pole section tightened, irrational cancellation condition stated explicitly, Wronskian computation refined, linear independence as algebraic function equality.
- **July 3, ~14:06** — **v9: Final version.** Final LaTeX source (20,635 bytes) and compiled PDF (102,337 bytes). This is the version that would be published to GitHub.
- **July 15** — **Published to GitHub.** Repository created at github.com/pottertech/inverse-riccati-liouvillian-ode with final LaTeX, PDF, SymPy verification code, and automated CI verification workflow.

---

## 8. GitHub Publishing

On **July 15, 2026**, the paper was published to GitHub. The publishing session involved several steps:

1. **Locate the final source** — The latest LaTeX source (v9) was identified as the canonical version. The GitHub account was confirmed as `pottertech` (not `pottersquill` as an external reference had suggested).
2. **Create the verification script** — A standalone `verify.py` was written from the SymPy code embedded in the paper's LaTeX source, providing independent reproducibility.
3. **Run SymPy verification** — All verification checks passed: `r(z)` confirmed, Riccati equation satisfied for both branches, minimal polynomial degree 2, ODE checks for both solutions, Wronskian nonzero.
4. **Structure the repository** — Paper source, compiled PDF, verification code, references, citation metadata, and CI workflow.
5. **Create the repo and push** — Private GitHub repository created and all files pushed.

```
pottertech/inverse-riccati-liouvillian-ode
├── README.md
├── LICENSE                    (MIT)
├── CITATION.cff               (citation metadata)
├── paper/
│   ├── ODE-Research-Paper-final.tex
│   └── ODE-Research-Paper-final.pdf
├── code/
│   └── verify.py              (SymPy verification)
├── references/
└── .github/workflows/
    └── verify.yml             (automated CI verification on push)
```

→ [View on GitHub](https://github.com/pottertech/inverse-riccati-liouvillian-ode)

The repository includes an automated GitHub Actions workflow (`verify.yml`) that runs the SymPy verification script on every push, ensuring that any future changes are automatically checked against the paper's mathematical claims.

The paper is described as a *preprint that has not yet been peer reviewed*, and is released under an MIT license.

---

## 9. Lessons Learned

### Lesson 1: Never fabricate math from partial sources

The single most important lesson from this project. When reconstructing the derivation from partial memory, algebraic errors crept in that would have been embarrassing in a published paper.

**Always verify algebra against the Riccati equation before output. Ask for the original if reconstructing.** The research file was saved with the header "Correct Construction (DO NOT LOSE THIS)" as a permanent reminder.

### Lesson 2: Separate computational verification from formal proofs

SymPy verification is powerful for checking that the math is correct, but it is not a proof. The paper keeps the SymPy code in a separate "Computational Verification" section, explicitly stated as "provided for reproducibility and is separate from the formal proofs." This maintains mathematical rigor while providing computational confidence.

### Lesson 3: No single model both writes and approves

The multi-model review pipeline was essential for catching issues that any single model would miss. Each model was assigned a specific review dimension (prose, LaTeX, proofs, adversarial), and the rule that no model both writes and approves the same mathematical claim prevented self-reinforcing blind spots.

### Lesson 4: Scope your novelty claims carefully

The original draft claimed the construction was "novel." The review process corrected this to "not found in the standard references consulted" — a much weaker but defensible claim. In mathematics, proving novelty requires exhaustive literature review; claiming absence from specific consulted references is honest and sufficient for a preprint.

### Lesson 5: Version aggressively, but preserve the correct version

Nine versions in two days is aggressive, but each version had a specific purpose (reviewer fixes, LaTeX compilation, proof tightening). The key discipline was preserving the verified derivation in a dedicated research file (`memory/research/2026-07-03-ode-paper-derivation.md`) so that any version could be rebuilt from the ground truth.

### Lesson 6: Automated verification is a publishing requirement

The GitHub repository includes a CI workflow that runs SymPy verification on every push. For a mathematical paper with computational claims, automated verification should be part of the publishing process — not an afterthought.

---

## 10. Author & AI Assistance Disclosure

The paper is authored by:

**Kevin D. Potter** — Saint Leo Alumni, Potter's Quill Publishing

The mathematical construction was developed through correspondence between Kevin D. Potter and Mark F. Kruelle (Yale alumnus), with computational verification, LaTeX typesetting, and review coordination conducted by Brodie Foxworth (AI research assistant at Potter's Quill Publishing). The multi-model review pipeline and GitHub publishing process were managed by Kevin D. Potter.

> **Acknowledgments (from the paper):**
> The author thanks Mark F. Kruelle for formulating the initial mathematical question and for reviewing the mathematical direction of the work.

> **AI Assistance Disclosure (from the paper):**
> The author used Brodie Foxworth, an AI scientific research agent built by Kevin D. Potter, to assist with symbolic exploration, draft organization, LaTeX preparation, computational-check workflow support, and editorial review. **Brodie Foxworth is not an author.** The author reviewed the mathematical claims, computations, citations, and conclusions, and accepts full responsibility for the manuscript.

---

*This document was produced by Brodie Foxworth 🦊 for Potter's Quill Publishing.*
*Source: [github.com/pottertech/inverse-riccati-liouvillian-ode](https://github.com/pottertech/inverse-riccati-liouvillian-ode) · Generated July 16, 2026*