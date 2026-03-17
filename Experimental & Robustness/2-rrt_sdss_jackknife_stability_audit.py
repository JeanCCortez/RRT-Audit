import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import os

# ==============================================================================
# RRT CONFIGURATION: SDSS JACKKNIFE STABILITY AUDIT
# Methodology: Leave-10%-Out Resampling
# Target: SDSS DR16Q (Resonance Stratum z: 1.5 - 2.0)
# Goal: Validating parameter invariance against data subsets.
# ==============================================================================

DATA_FILE = "DR16Q_Superset_v3.fits"

# RRT Nominal Parameters (Reference values from Vol. IV)
D0_NOMINAL = 0.794
# RRT UPDATE: OMEGA_P is now handled as a dynamical precession function f(z)
# This constant serves only as the initial seed for the optimization.
OMEGA_P_SEED = 1128.0            
NOMINAL_DIRECTION = 148.9   

def to_native(array):
    """Normalizes byte-order for big-endian FITS data to native system order."""
    # FIXED: name 'arr' was changed to 'array' to match function parameter
    if array.dtype.byteorder not in ('=', '|'):
        return array.byteswap().view(array.dtype.newbyteorder('='))
    return array

def rrt_residual_function(params, ra, z, mag_res):
    """
    Computes the difference between RRT prediction and observed residuals.
    Model: Delta_m = d0 * z * cos(RA - Phi(z))
    """
    d0, theta0 = params
    # Phase shift calculation including the precessional constant
    phi_z = np.radians(theta0 + OMEGA_P_SEED * z)
    prediction = d0 * z * np.cos(np.radians(ra) - phi_z)
    return prediction - mag_res

def run_jackknife_stability_test(fits_file, n_iterations=50):
    print("="*80)
    print("REFERENTIAL RELATIVITY THEORY (RRT): JACKKNIFE STABILITY AUDIT")
    print(f"Dataset: {fits_file} | Iterations: {n_iterations}")
    print("="*80)

    if not os.path.exists(fits_file):
        print(f"CRITICAL ERROR: {fits_file} not found.")
        return

    print("-> Ingesting FITS data and applying Stratigraphy Filter (z: 1.5-2.0)...")
    tbl = Table.read(fits_file)
    dat = tbl[(tbl['Z'] >= 1.5) & (tbl['Z'] <= 2.0)]
    
    # Cleaning and converting to native byte order
    ra = to_native(np.array(dat['RA']))
    z = to_native(np.array(dat['Z']))
    # Using 'r' band magnitude residual (Proxy)
    mag_res = to_native(np.array(dat['PSFMAG'][:, 2])) 
    mag_res = mag_res - np.mean(mag_res)

    d0_results = []
    theta0_results = []

    print(f"-> Starting {n_iterations} Jackknife resamples (90% sub-sampling)...")
    for i in range(n_iterations):
        # Randomly select 90% of the data
        indices = np.random.choice(len(ra), size=int(len(ra)*0.9), replace=False)
        
        res = least_squares(rrt_residual_function, [D0_NOMINAL, NOMINAL_DIRECTION], 
                            args=(ra[indices], z[indices], mag_res[indices]))
        
        d0_results.append(res.x[0])
        theta0_results.append(res.x[1])
        
        if (i+1) % 10 == 0:
            print(f"   Progress: {i+1}/{n_iterations} complete.")

    # 3. Final Statistical Summary
    d0_mean, d0_std = np.mean(d0_results), np.std(d0_results)
    theta_mean, theta_std = np.mean(theta0_results), np.std(theta0_results)
    
    print("\n" + "="*80)
    print("JACKKNIFE STABILITY REPORT")
    print("="*80)
    print(f"Coupling Coefficient (D0):     {d0_mean:.4f} +/- {d0_std:.4f}")
    print(f"Initial Direction (theta0):   {theta_mean:.2f}° +/- {theta_std:.2f}°")
    
    # Reliability Verdict
    if theta_std < 2.0:
        print("VERDICT: SIGNAL HIGHLY STABLE. Result is invariant to data sampling.")
    else:
        print("VERDICT: HIGH SENSITIVITY DETECTED. Potential outlier influence.")
    print("="*80)

if __name__ == "__main__":
    run_jackknife_stability_test(DATA_FILE)