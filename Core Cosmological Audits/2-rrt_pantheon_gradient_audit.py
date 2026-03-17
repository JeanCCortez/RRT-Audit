import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import least_squares
import os

# ==============================================================================
# RRT CONFIGURATION: PANTHEON+ GRADIENT AUDIT
# Target: Type Ia Supernovae (Standard Candles)
# Goal: Detection of the Dipole Anisotropy Gradient (D0)
# ==============================================================================
DATA_FILE = "PantheonPlusSH0ES.csv"
REDSHIFT_CUTOFF = 0.02  # Minimizing local peculiar velocity bias
MC_ITERATIONS = 100     # Monte Carlo shuffles for significance testing

def run_pantheon_plus_gradient_audit(file_name):
    """
    Performs a high-precision audit on the Pantheon+ catalog to test
    the Referential Relativity Theory (RRT) phase-gradient model.
    """
    print("="*80)
    print(f"RRT SCIENTIFIC AUDIT: {file_name}")
    print("Protocol: Hubble Detrending + Monte Carlo Spatial Shuffle")
    print("="*80)

    if not os.path.exists(file_name):
        print(f"ERROR: {file_name} not found in the current directory.")
        return

    try:
        # 1. Data Ingestion & Column Mapping
        df = pd.read_csv(file_name, sep=None, engine='python')
        cols = {c.lower(): c for c in df.columns}
        
        # Standard Pantheon+ naming conventions
        z_col = cols.get('zcmb') or cols.get('z')
        mu_col = cols.get('mu_shoes') or cols.get('mu_pantheon') or cols.get('m_b_corr')
        ra_col = cols.get('ra')
        dec_col = cols.get('dec')
        err_col = cols.get('mu_err') or cols.get('mu_err_shoes')

        if mu_col is None:
            raise KeyError("Distance Modulus (MU) column not identified.")

        # 2. Pre-processing & Quality Cuts
        # Filtering for z > 0.02 to enter the Hubble flow regime (Phase 2/3 Transition)
        df = df.dropna(subset=[z_col, mu_col, ra_col, dec_col]).copy()
        df = df[df[z_col] > REDSHIFT_CUTOFF].copy()

        z = df[z_col].values.astype(float)
        mu_obs = df[mu_col].values.astype(float)
        mu_err = df[err_col].values.astype(float) if err_col else np.ones_like(z) * 0.15

        # 3. Coordinate Transformation (ICRS to Galactic)
        coords = SkyCoord(ra=df[ra_col].values*u.degree, dec=df[dec_col].values*u.degree, frame='icrs')
        l_gal = coords.galactic.l.radian
        b_gal = coords.galactic.b.radian
        
        # 4. Hubble Detrending
        # Isotropic baseline removal to isolate the anisotropic residual signal
        residuals = mu_obs - (5 * np.log10(z))
        residuals -= np.mean(residuals)
        weights = 1.0 / mu_err 

        # 5. RRT Dipole Model Function
        def rrt_cost_function(params, l, b, redshift, res, w):
            """
            Calculates the difference between the RRT prediction and observed residuals.
            Prediction: Delta_m = D0 * z * cos(theta)
            """
            d0, lp, bp = params
            # Spherical cosine law for angular separation from the dipole axis
            cos_theta = np.sin(b) * np.sin(bp) + np.cos(b) * np.cos(bp) * np.cos(l - lp)
            return ((d0 * redshift) * cos_theta - res) * w

        # 6. Primary Optimization (Real Data Fit)
        print("-> Optimizing Anisotropic Gradient (D0) on observed dataset...")
        # Initial guess based on RRT Vol I findings
        initial_guess = [0.1, np.radians(148), np.radians(-5)]
        fit_bounds = ([0, 0, -np.pi/2], [2.0, 2*np.pi, np.pi/2])
        
        results_real = least_squares(rrt_cost_function, initial_guess, 
                                    args=(l_gal, b_gal, z, residuals, weights), 
                                    bounds=fit_bounds)
        
        d0_final = results_real.x[0]
        l_final_deg = np.degrees(results_real.x[1])
        b_final_deg = np.degrees(results_real.x[2])

        # 7. Monte Carlo Robustness Test (Falsification Protocol)
        # Shuffling residuals to determine the probability of a chance alignment
        print(f"-> Starting Monte Carlo Null-Hypothesis Test ({MC_ITERATIONS} iterations)...")
        shuffled_d0_distribution = []
        residuals_shuffled = residuals.copy()
        
        for i in range(MC_ITERATIONS):
            np.random.shuffle(residuals_shuffled) 
            fit_mc = least_squares(rrt_cost_function, initial_guess, 
                                  args=(l_gal, b_gal, z, residuals_shuffled, weights), 
                                  bounds=fit_bounds)
            shuffled_d0_distribution.append(fit_mc.x[0])
            if (i+1) % 20 == 0: print(f"   Iteration {i+1}/{MC_ITERATIONS} complete.")

        # 8. Statistical Inference & Model Comparison
        # Calculating Significance (Sigma) and AIC (Akaike Information Criterion)
        z_score = (d0_final - np.mean(shuffled_d0_distribution)) / np.std(shuffled_d0_distribution)
        
        n_samples = len(z)
        rss_rrt = np.sum(results_real.fun**2)
        rss_iso = np.sum((residuals * weights)**2)
        
        # AIC: Lower values indicate better model efficiency
        aic_rrt = 2*3 + n_samples * np.log(rss_rrt/n_samples)
        aic_iso = 2*0 + n_samples * np.log(rss_iso/n_samples)
        delta_aic = aic_rrt - aic_iso

        print("\n" + "="*80)
        print("FINAL AUDIT VERDICT")
        print("="*80)
        print(f"STATISTICAL SIGNIFICANCE: {z_score:.2f} σ")
        print(f"DELTA AIC (RRT vs Isotropic): {delta_aic:.2f} (Negative favors RRT)")
        print(f"ANISOTROPY GRADIENT (D0): {d0_final:.6f}")
        print(f"DIPOLE DIRECTION: Galactic l={l_final_deg:.2f}°, b={b_final_deg:.2f}°")
        print(f"USEFUL SAMPLE SIZE: {n_samples} Supernovae")
        print("="*80)

    except Exception as e:
        print(f"CRITICAL AUDIT ERROR: {e}")

if __name__ == "__main__":
    run_pantheon_plus_gradient_audit(DATA_FILE)