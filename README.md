# An Inverse Riccati Construction of Second-Order Linear ODEs from Degree-2 Algebraic Extensions

This repository contains the manuscript, LaTeX source, and symbolic
verification code for the paper by Kevin D. Potter.

## Contents

- `paper/` — manuscript source (`ODE-Research-Paper-final.tex`) and compiled PDF
- `code/verify.py` — SymPy verification of all key identities in the paper
- `references/` — reference notes
- `CITATION.cff` — citation metadata

## Summary

We present a systematic construction of a family of second-order linear ODEs

> y'' + r(z)·y = 0,  where  r(z) = 3/(16z²(z-1)²)

whose solutions are Liouvillian and arise from a logarithmic derivative
s(z) ∈ ℂ(z, √R) of algebraic degree exactly 2 over ℂ(z). The construction
places two poles at z=0 and z=1 with equal residues 1/4.

## Running the Verification

```bash
pip install sympy
python3 code/verify.py
```

## Status

This manuscript is a preprint and has not yet been peer reviewed.

## License

MIT — see `LICENSE` for details.
