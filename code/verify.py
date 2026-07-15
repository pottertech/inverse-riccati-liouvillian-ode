#!/usr/bin/env python3
"""
SymPy verification for:
  "An Inverse Riccati Construction of Second-Order Linear ODEs from Degree-2 Algebraic Extensions"

Verifies computationally:
  1. r(z) = 3/(16 z^2 (z-1)^2)  (from A' + A^2 + B^2 R)
  2. Riccati equation: s' + s^2 + r = 0  for both branches s_+, s_-
  3. Minimal polynomial candidate is quadratic in s
  4. Key identity: (2z-1)^2 - 4z(z-1) = 1
  5. ODE check: y'' + r*y = 0  for both solutions y_1, y_2
  6. Wronskian (via direct computation and logarithmic derivatives)
  7. Ratio y1/y2 is nonconstant (derivative is nonzero)

Note: Exact algebraic degree 2 (irreducibility of the minimal polynomial)
requires the proof that R is not a square in C(z), which is established
in the paper's Theorem 2. The script verifies the polynomial is quadratic
in s; irreducibility follows from the nonsquare argument in the proof.

Author: Kevin D. Potter
"""

import sympy as sp

z = sp.Symbol('z')

# --- Define A, B, R ---
A = (2*z - 1) / (4*z*(z - 1))
B = sp.Rational(1, 2)
R = 1 / (z*(z - 1))

# --- Define s_+ and s_- ---
s_plus  = A + B * sp.sqrt(R)
s_minus = A - B * sp.sqrt(R)

# --- 1. Verify r(z) = -(A' + A^2 + B^2*R) ---
r = -(sp.diff(A, z) + A**2 + B**2 * R)
r_simplified = sp.simplify(r)
expected_r = sp.Rational(3, 16) / (z**2 * (z - 1)**2)
print("1. r(z) =", r_simplified)
assert sp.simplify(r_simplified - expected_r) == 0, "r(z) mismatch!"
print("   ✓ r(z) = 3/(16 z^2 (z-1)^2)  CORRECT")
print()

# --- 2. Verify Riccati: s' + s^2 + r = 0 for both branches ---
ric_plus  = sp.simplify(sp.diff(s_plus, z)  + s_plus**2  + r)
ric_minus = sp.simplify(sp.diff(s_minus, z) + s_minus**2 + r)
print("2. Riccati s+ :", sp.simplify(ric_plus))
print("   Riccati s- :", sp.simplify(ric_minus))
assert ric_plus == 0, "Riccati s+ failed!"
assert ric_minus == 0, "Riccati s- failed!"
print("   ✓ Both branches satisfy Riccati equation")
print()

# --- 3. Verify minimal polynomial is degree 2 ---
s = sp.Symbol('s')
minpoly = s**2 - 2*A*s + (A**2 - B**2*R)
minpoly_simplified = sp.simplify(minpoly)
print("3. Minimal polynomial:", minpoly_simplified)
# Check it's degree 2 in s
poly_in_s = sp.Poly(minpoly_simplified, s)
print("   Degree in s:", poly_in_s.degree())
assert poly_in_s.degree() == 2, "Minimal polynomial not degree 2!"
print("   ✓ Candidate polynomial is monic and quadratic in s")
print("   Exact degree 2 additionally uses the proof that R is not a square in C(z)")
print()

# --- 4. Verify key identity: (2z-1)^2 - 4z(z-1) = 1 ---
identity = sp.expand((2*z - 1)**2 - 4*z*(z - 1))
print("4. (2z-1)^2 - 4z(z-1) =", identity)
assert identity == 1, "Key identity failed!"
print("   ✓ Identity holds: u·v = 1")
print()

# --- 5. Verify solutions satisfy ODE: y'' + r*y = 0 ---
y1 = (z*(z-1))**sp.Rational(1,4) * (2*z - 1 + 2*sp.sqrt(z*(z-1)))**sp.Rational(1,2)
y2 = (z*(z-1))**sp.Rational(1,4) * (2*z - 1 - 2*sp.sqrt(z*(z-1)))**sp.Rational(1,2)

ode_check_1 = sp.simplify(sp.diff(y1, z, 2) + r * y1)
ode_check_2 = sp.simplify(sp.diff(y2, z, 2) + r * y2)
print("5. ODE check y1:", sp.simplify(ode_check_1))
print("   ODE check y2:", sp.simplify(ode_check_2))
assert ode_check_1 == 0, f"ODE check for y1 failed: {ode_check_1}"
assert ode_check_2 == 0, f"ODE check for y2 failed: {ode_check_2}"
print("   ✓ Both solutions satisfy y'' + r(z)y = 0")
print()

# --- 6. Wronskian (W^2 = 1 via algebraic identity) ---
W = sp.simplify(y1 * sp.diff(y2, z) - sp.diff(y1, z) * y2)
print("6. Wronskian W(y1, y2):", W)
# W is branch-sensitive in explicit form. Use the algebraic identity:
#   y1*y2 = [z(z-1)]^(1/2) * (uv)^(1/2) = ±sqrt(z(z-1))
#   s_- - s_+ = -1/sqrt(z(z-1))
#   W = y1*y2*(s_- - s_+) = ±1
# So W^2 = z(z-1) * (s_- - s_+)^2 should be identically 1
W_squared = sp.simplify((z*(z-1)) * (s_minus - s_plus)**2)
print(f"   W^2 = z(z-1)*(s_- - s_+)^2 = {W_squared}")
assert W_squared == 1, f"W^2 should be 1, got {W_squared}"
print("   ✓ W^2 = 1 (branch-aware: W = ±1, constant and nonzero)")
print()

# --- 7. Ratio y1/y2 is nonconstant (via u' ≠ 0) ---
ratio = sp.simplify(y1 / y2)
print("7. y1/y2 =", ratio)
# Formal proof: y1/y2 = ±u where u = 2z-1+2√(z(z-1)), nonconstant
dratio = sp.simplify(sp.diff(ratio, z))
print("   d/dz(y1/y2) =", dratio)
# Stronger check: u is nonconstant
u = 2*z - 1 + 2*sp.sqrt(z*(z-1))
du = sp.simplify(sp.diff(u, z))
print(f"   u = {u}")
print(f"   u' = {du}")
assert not sp.simplify(du).equals(0), "u' = 0 — linear independence failed!"
print("   ✓ u' ≠ 0 → u nonconstant → y1/y2 nonconstant → linearly independent")
print()

print("=" * 50)
print("ALL CHECKS PASSED ✓")
print("=" * 50)
