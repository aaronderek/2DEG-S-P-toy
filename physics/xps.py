"""
XPS (X-ray Photoelectron Spectroscopy) modeling and adsorbate effects.
"""

import numpy as np
from .constants import Q, EPS0, DEBYE_TO_CM
from .units import nm_to_m, eV_to_J, J_to_eV, degrees_to_radians


def calculate_xps_weight(z_array, lambda_nm, theta_deg):
    """
    Calculate XPS sampling weight function.

    The weight function accounts for the finite escape depth of photoelectrons.

    Parameters
    ----------
    z_array : array
        Array of depth positions in meters
    lambda_nm : float
        Inelastic mean free path in nanometers
    theta_deg : float
        Detection angle in degrees (0° = normal emission)

    Returns
    -------
    w_z : array
        Normalized weight function at each depth
    """
    lambda_m = nm_to_m(lambda_nm)
    theta_rad = degrees_to_radians(theta_deg)

    # Effective escape depth
    lambda_eff = lambda_m * np.cos(theta_rad)

    # Weight: w(z) = (1/λeff) * exp(-z/λeff)
    w_z = (1 / lambda_eff) * np.exp(-z_array / lambda_eff)

    return w_z


def calculate_core_level_shift(V_z, z_array, lambda_nm, theta_deg):
    """
    Calculate XPS core level shift due to band bending.

    The shift is calculated as a weighted average of the potential profile.

    Parameters
    ----------
    V_z : array
        Potential energy profile in Joules at positions z_array
    z_array : array
        Array of depth positions in meters
    lambda_nm : float
        Inelastic mean free path in nanometers
    theta_deg : float
        Detection angle in degrees

    Returns
    -------
    Delta_E_CL : float
        Core level shift in eV
    eta : float
        Effective sampling factor (0 < η < 1)
    """
    # Get weight function
    w_z = calculate_xps_weight(z_array, lambda_nm, theta_deg)

    # Calculate weighted average: ΔE_CL = -∫ w(z)·V(z) dz
    # The negative sign is because we measure binding energy shift
    Delta_E_CL_J = -np.trapz(w_z * V_z, z_array)
    Delta_E_CL_eV = J_to_eV(Delta_E_CL_J)

    # Calculate effective sampling factor η
    # If V(z) is approximately linear near surface: V(z) ≈ -Φs for z<<W
    # Then η = ΔE_CL / Φs
    # We'll approximate Φs from the maximum of V(z)
    if np.any(V_z != 0):
        Phi_s_approx = J_to_eV(np.max(np.abs(V_z)))
        if Phi_s_approx > 1e-6:  # Avoid division by zero
            eta = abs(Delta_E_CL_eV) / Phi_s_approx
        else:
            eta = 0
    else:
        eta = 0

    return Delta_E_CL_eV, eta


def calculate_dipole_shift(N_ads_cm2, mu_debye):
    """
    Calculate work function shift due to surface adsorbate dipoles.

    Uses the Helmholtz equation for a dipole layer.

    Parameters
    ----------
    N_ads_cm2 : float
        Adsorbate density in cm⁻²
    mu_debye : float
        Perpendicular dipole moment in Debye

    Returns
    -------
    Delta_Phi_dip : float
        Dipole-induced work function shift in eV
    """
    # Convert units
    N_ads = N_ads_cm2 * 1e4  # cm⁻² to m⁻²
    mu_C_m = mu_debye * DEBYE_TO_CM  # Debye to C·m

    # Helmholtz equation: ΔΦdip = -(Nads * μ⊥) / ε₀
    # Result is in Volts (V), which equals eV for electron work function
    Delta_Phi_dip_V = -(N_ads * mu_C_m) / EPS0

    # For electric potential: 1 V = 1 eV (no conversion needed)
    Delta_Phi_dip_eV = Delta_Phi_dip_V

    return Delta_Phi_dip_eV


def calculate_dipole_from_coverage(coverage, N_site_cm2, mu_debye):
    """
    Calculate dipole shift from coverage fraction.

    Parameters
    ----------
    coverage : float
        Adsorbate coverage (0 to 1)
    N_site_cm2 : float
        Surface site density in cm⁻²
    mu_debye : float
        Dipole moment per adsorbate in Debye

    Returns
    -------
    Delta_Phi_dip : float
        Work function shift in eV
    """
    N_ads_cm2 = coverage * N_site_cm2
    return calculate_dipole_shift(N_ads_cm2, mu_debye)


class XPSModel:
    """
    Model for XPS measurements including band bending and adsorbate effects.
    """

    def __init__(self, lambda_nm=1.8, theta_deg=0, N_site_cm2=1e15):
        """
        Initialize XPS model.

        Parameters
        ----------
        lambda_nm : float
            Inelastic mean free path in nm
        theta_deg : float
            Detection angle in degrees
        N_site_cm2 : float
            Surface site density in cm⁻² (default ~1e15 for oxide surfaces)
        """
        self.lambda_nm = lambda_nm
        self.theta_deg = theta_deg
        self.N_site_cm2 = N_site_cm2

    def calculate_shift(self, potential_model, Phi_s_eV, z_max_factor=3):
        """
        Calculate XPS core level shift for a given potential model.

        Parameters
        ----------
        potential_model : object
            Model object with get_potential method
        Phi_s_eV : float
            Surface potential in eV
        z_max_factor : float
            Maximum depth as factor of lambda (default 3)

        Returns
        -------
        Delta_E_CL : float
            Core level shift in eV
        eta : float
            Effective sampling factor
        """
        # Create z array for integration
        lambda_m = nm_to_m(self.lambda_nm)
        z_max = z_max_factor * lambda_m
        z_array = np.linspace(0, z_max, 1000)

        # Get potential profile
        V_z = potential_model.get_potential(Phi_s_eV, z_array)

        # Calculate core level shift
        Delta_E_CL, eta = calculate_core_level_shift(
            V_z, z_array, self.lambda_nm, self.theta_deg
        )

        return Delta_E_CL, eta

    def add_adsorbate_shift(self, Delta_WF_bend, coverage, mu_debye):
        """
        Add adsorbate dipole shift to band bending contribution.

        Parameters
        ----------
        Delta_WF_bend : float or array
            Work function change from band bending in eV
        coverage : float
            Adsorbate coverage (0 to 1)
        mu_debye : float
            Dipole moment in Debye

        Returns
        -------
        Delta_WF_total : float or array
            Total work function change in eV
        """
        Delta_Phi_dip = calculate_dipole_from_coverage(
            coverage, self.N_site_cm2, mu_debye
        )

        return Delta_WF_bend + Delta_Phi_dip

    def set_parameters(self, lambda_nm=None, theta_deg=None):
        """Update XPS parameters."""
        if lambda_nm is not None:
            self.lambda_nm = lambda_nm
        if theta_deg is not None:
            self.theta_deg = theta_deg
