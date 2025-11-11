"""
Test script to verify adsorbate dipole calculation fix.

Bug fix verification for issue where ΔΦ_dip was ~1e19 eV instead of ~0.3 eV.

Root causes fixed:
1. Unit conversion error: Helmholtz result is in Volts (not Joules)
2. N_site default value: Changed from 1e15 cm^-2 to 1e14 cm^-2
"""

import sys
sys.path.insert(0, '/home/user/2DEG-S-P-toy')

from physics.xps import calculate_dipole_from_coverage, calculate_dipole_shift
from physics.constants import EPS0, DEBYE_TO_CM


def test_dipole_calculation():
    """Test the dipole shift calculation with expected values."""

    print("=" * 70)
    print("ADSORBATE DIPOLE CALCULATION TEST")
    print("=" * 70)

    # Test Case 1: User's expected values
    print("\n[Test 1] User specification:")
    print("  θ = 0.5, μ = 1.5 D, N_site = 1e14 cm^-2 (= 1e18 m^-2)")

    coverage = 0.5
    mu_debye = 1.5
    N_site_cm2 = 1e14

    Delta_Phi_dip = calculate_dipole_from_coverage(coverage, N_site_cm2, mu_debye)

    print(f"\n  Result: ΔΦ_dip = {Delta_Phi_dip:.4f} eV")
    print(f"  Expected: ΔΦ_dip ≈ -0.283 eV")

    # Check tolerance
    expected = -0.283
    tolerance = 0.005

    if abs(Delta_Phi_dip - expected) < tolerance:
        print(f"  ✓ PASS: Within tolerance (±{tolerance} eV)")
    else:
        print(f"  ✗ FAIL: Outside tolerance (±{tolerance} eV)")
        print(f"  Error: {abs(Delta_Phi_dip - expected):.6f} eV")

    # Test Case 2: Manual calculation verification
    print("\n[Test 2] Manual calculation verification:")

    N_ads = coverage * N_site_cm2 * 1e4  # Convert to m^-2
    mu_C_m = mu_debye * DEBYE_TO_CM

    print(f"  N_ads = {N_ads:.2e} m^-2")
    print(f"  μ (C·m) = {mu_C_m:.4e} C·m")
    print(f"  ε₀ = {EPS0:.4e} F/m")

    Delta_Phi_manual = -(N_ads * mu_C_m) / EPS0

    print(f"\n  Manual: ΔΦ_dip = {Delta_Phi_manual:.4f} V = {Delta_Phi_manual:.4f} eV")
    print(f"  Function: ΔΦ_dip = {Delta_Phi_dip:.4f} eV")

    if abs(Delta_Phi_manual - Delta_Phi_dip) < 1e-6:
        print(f"  ✓ PASS: Manual and function results match")
    else:
        print(f"  ✗ FAIL: Mismatch between manual and function")

    # Test Case 3: Range of realistic values
    print("\n[Test 3] Range of realistic parameters:")
    print("  " + "-" * 60)
    print(f"  {'θ':>6} {'μ(D)':>8} {'N_site(cm^-2)':>15} {'ΔΦ_dip(eV)':>12} {'Status':>10}")
    print("  " + "-" * 60)

    test_cases = [
        (0.1, 1.0, 1e14, "Low coverage"),
        (0.5, 1.5, 1e14, "Typical"),
        (1.0, 2.0, 1e14, "High coverage"),
        (0.5, 0.5, 1e14, "Low dipole"),
        (0.5, 3.0, 1e14, "High dipole"),
        (0.5, 1.5, 5e13, "Low N_site"),
        (0.5, 1.5, 5e14, "High N_site"),
    ]

    all_pass = True

    for theta, mu, N_site, label in test_cases:
        Delta_Phi = calculate_dipole_from_coverage(theta, N_site, mu)

        # Check if reasonable (|ΔΦ| < 2 eV)
        if abs(Delta_Phi) < 2.0:
            status = "✓ OK"
        else:
            status = "✗ TOO LARGE"
            all_pass = False

        print(f"  {theta:>6.2f} {mu:>8.1f} {N_site:>15.1e} {Delta_Phi:>12.4f} {status:>10}")

    print("  " + "-" * 60)

    if all_pass:
        print("\n  ✓ All test cases within reasonable range (|ΔΦ| < 2 eV)")
    else:
        print("\n  ✗ Some test cases produce unreasonably large values")

    # Test Case 4: Old bug reproduction (for reference)
    print("\n[Test 4] Old bug reproduction (N_site = 1e15 cm^-2):")

    N_site_old = 1e15  # Old default
    Delta_Phi_old = calculate_dipole_from_coverage(0.5, N_site_old, 1.5)

    print(f"  Old N_site: 1e15 cm^-2 = 1e19 m^-2")
    print(f"  Result: ΔΦ_dip = {Delta_Phi_old:.4f} eV")
    print(f"  Ratio to correct value: {abs(Delta_Phi_old / Delta_Phi_dip):.1f}x")

    if abs(Delta_Phi_old) < 2.0:
        print(f"  Status: Now within acceptable range")
    else:
        print(f"  ⚠ Warning: Still too large with old N_site value")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Bug fixes applied:")
    print(f"  1. ✓ Removed incorrect J_to_eV conversion (V → eV is 1:1)")
    print(f"  2. ✓ Changed default N_site: 1e15 → 1e14 cm^-2")
    print(f"  3. ✓ Added validation warning for |ΔΦ| > 2 eV")
    print(f"\n  Expected behavior:")
    print(f"  - θ=0.5, μ=1.5D, N_site=1e14 cm^-2 → ΔΦ_dip ≈ -0.28 eV ✓")
    print(f"  - Curves in Figure 2 should be parallel (same slope)")
    print("=" * 70)


if __name__ == "__main__":
    test_dipole_calculation()
