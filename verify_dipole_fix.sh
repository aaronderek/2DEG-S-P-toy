#!/bin/bash
# Quick verification of dipole calculation fix

echo "========================================================================"
echo "ADSORBATE DIPOLE BUG FIX VERIFICATION"
echo "========================================================================"
echo ""
echo "Bug Description:"
echo "  - ΔΦ_dip was ~1e19 eV (should be ~0.3 eV)"
echo "  - Figure 2 curves became horizontal lines"
echo ""
echo "Root Causes:"
echo "  1. Unit conversion error: Used J_to_eV on Volts (wrong!)"
echo "     - Helmholtz result is in Volts, not Joules"
echo "     - For electric potential: 1 V = 1 eV (no conversion)"
echo ""
echo "  2. N_site too large: 1e15 cm^-2 = 1e19 m^-2 (too big!)"
echo "     - Should be 1e14 cm^-2 = 1e18 m^-2 for oxide surfaces"
echo ""
echo "========================================================================"
echo "MANUAL CALCULATION"
echo "========================================================================"
echo ""
echo "Test case: θ=0.5, μ=1.5 D, N_site=1e14 cm^-2"
echo ""

# Use Python for calculation
python3 << 'EOF'
# Constants
EPS0 = 8.854187817e-12  # F/m
DEBYE_TO_C_M = 3.33564e-30  # C·m
Q = 1.602176634e-19  # C

# Input parameters
theta = 0.5
mu_debye = 1.5
N_site_cm2 = 1e14

# Calculate
N_ads_cm2 = theta * N_site_cm2
N_ads_m2 = N_ads_cm2 * 1e4  # cm^-2 to m^-2

mu_C_m = mu_debye * DEBYE_TO_C_M

# Helmholtz equation
Delta_Phi_dip_V = -(N_ads_m2 * mu_C_m) / EPS0

# For electric potential: 1 V = 1 eV
Delta_Phi_dip_eV = Delta_Phi_dip_V

print(f"Step-by-step calculation:")
print(f"  N_ads = θ × N_site = {theta} × {N_site_cm2:.1e} cm^-2 = {N_ads_cm2:.1e} cm^-2")
print(f"  N_ads (m^-2) = {N_ads_m2:.2e} m^-2")
print(f"  μ (C·m) = {mu_debye} D × {DEBYE_TO_C_M:.4e} = {mu_C_m:.4e} C·m")
print(f"  ε₀ = {EPS0:.4e} F/m")
print(f"")
print(f"  ΔΦ_dip = -(N_ads × μ) / ε₀")
print(f"         = -({N_ads_m2:.2e} × {mu_C_m:.4e}) / {EPS0:.4e}")
print(f"         = {Delta_Phi_dip_V:.4f} V")
print(f"         = {Delta_Phi_dip_eV:.4f} eV  (1 V = 1 eV for potential)")
print(f"")
print(f"Expected result: ΔΦ_dip ≈ -0.283 eV")
print(f"Actual result:   ΔΦ_dip = {Delta_Phi_dip_eV:.3f} eV")
print(f"")

if abs(Delta_Phi_dip_eV - (-0.283)) < 0.005:
    print("✓ PASS: Result matches expected value!")
else:
    print("✗ FAIL: Result does not match expected value")

print("")
print("========================================================================")
print("BEFORE FIX (for reference)")
print("========================================================================")
print("")

# Old bug behavior
N_site_old = 1e15  # Old default
N_ads_old_cm2 = theta * N_site_old
N_ads_old_m2 = N_ads_old_cm2 * 1e4

# WRONG calculation (old code)
Delta_Phi_wrong_J = -(N_ads_old_m2 * mu_C_m) / EPS0  # This is actually in V
Delta_Phi_wrong_eV = Delta_Phi_wrong_J / Q  # Wrong! Dividing by Q

print(f"Old N_site: {N_site_old:.1e} cm^-2 = {N_site_old * 1e4:.1e} m^-2")
print(f"")
print(f"Old (WRONG) calculation:")
print(f"  1. Delta_Phi_'J' = -({N_ads_old_m2:.2e} × {mu_C_m:.4e}) / {EPS0:.4e}")
print(f"                   = {Delta_Phi_wrong_J:.4f} V  (but treated as J)")
print(f"  2. Delta_Phi_eV = {Delta_Phi_wrong_J:.4f} / {Q:.4e}  (WRONG!)")
print(f"                  = {Delta_Phi_wrong_eV:.2e} eV  (HUGE!)")
print(f"")
print(f"This caused:")
print(f"  - ΔΦ_dip ≈ {Delta_Phi_wrong_eV:.2e} eV (absurdly large)")
print(f"  - Figure 2 became horizontal line")
print(f"  - Sidebar showed garbage values")
print("")

# With corrected N_site but old conversion
Delta_Phi_half_wrong_V = -(N_ads_m2 * mu_C_m) / EPS0
Delta_Phi_half_wrong_eV = Delta_Phi_half_wrong_V / Q

print(f"With corrected N_site but old J_to_eV conversion:")
print(f"  ΔΦ_dip = {Delta_Phi_half_wrong_eV:.2e} eV (still wrong!)")
print(f"")
print("========================================================================")
print("SUMMARY")
print("========================================================================")
print("")
print("Fixes applied:")
print("  1. ✓ Removed J_to_eV conversion (V → eV is 1:1 for potential)")
print("  2. ✓ Changed N_site: 1e15 → 1e14 cm^-2")
print("  3. ✓ Added warning if |ΔΦ_dip| > 2 eV")
print("")
print("Expected behavior:")
print("  - Typical ΔΦ_dip: 0.1 - 1 eV (not 1e19 eV!)")
print("  - Figure 2: parallel lines (same slope, different intercept)")
print("  - 'With ads' curve: vertical shift of 'no ads' curve")
print("")
print("========================================================================"
print("")

EOF

echo "Verification complete!"
echo ""
echo "To test in the app:"
echo "  1. Run: streamlit run app.py"
echo "  2. Enable 'Show with adsorbates'"
echo "  3. Set: θ=0.5, μ=1.5D"
echo "  4. Check: ΔΦ_dip should show ≈ -0.28 eV (not 1e19 eV)"
echo "  5. Verify: Figure 2 shows two parallel lines"
echo ""
