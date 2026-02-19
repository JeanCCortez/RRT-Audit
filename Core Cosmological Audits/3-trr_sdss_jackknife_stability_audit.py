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

# Relative path for repository portability
DATA_FILE = "DR16Q_Superset_v3.fits"

# RRT Nominal Parameters (Reference values from Vol. IV)
D0_NOMINAL = 0.794
OMEGA_P = 1128.0            # Precession Constant (deg/z)
NOMINAL_DIRECTION = 148.9   # Primordial Axis (Degrees)

def to_native(array):
    """Normalizes byte-order for big-endian FITS data to native system order."""
    if array.dtype.byteorder not in ('=', '|'):
        return array.byteswap().view(arr.dtype.newbyteorder('='))
    return array

def rrt_residual_function(params, ra, z, mag_res):
    """
    Computes the difference between RRT prediction and observed residuals.
    Model: Delta_m = d0 * z * cos(RA - Phi(z))
    """
    d0, theta0 = params
    # Cortez Precession Model
    phase_z = (theta0 + (OMEGA_P / z)) % 360
    prediction = d0 * z * np.cos(np.radians(ra - phase_z))
    return prediction - mag_res

def run_jackknife_stability_test(file_path, n_iterations=50):
    """
    Executes the Jackknife audit by randomly removing 10% of the dataset
    in each iteration to check for parameter drift.
    """
    print("="*80)
    print("REFERENTIAL RELATIVITY THEORY (RRT): JACKKNIFE STABILITY AUDIT")
    print(f"Dataset: {file_path} | Iterations: {n_iterations}")
    print("="*80)

    if not os.path.exists(file_path):
        print(f"CRITICAL ERROR: {file_path} not found.")
        return

    # 1. Data Ingestion & Resonance Stratum Filtering
    print("-> Ingesting FITS data and applying Stratigraphy Filter (z: 1.5-2.0)...")
    dat = Table.read(file_path, format='fits')
    
    ra = to_native(np.array(dat['RA']))
    z = to_native(np.array(dat['Z']))
    mag_i = to_native(np.array(dat['PSFMAG'][:, 3]))
    
    df = pd.DataFrame({'ra': ra, 'z': z, 'mag': mag_i})
    # Entering the Phase 3 (Viscous) Resonance Layer
    df = df[(df['z'] >= 1.5) & (df['z'] <= 2.0) & (df['mag'] > 0)].copy()
    
    # Hubble Detrending
    df['mag_res'] = df['mag'] - (5 * np.log10(df['z']))
    
    print(f"   Total filtered sample: {len(df)} objects.")
    
    # Storage for Jackknife coefficients
    d0_results = []
    theta0_results = []
    
    # 2. Resampling Loop
    print(f"-> Starting resampling (Removing 10% data per cut)...")
    for i in range(n_iterations):
        # Sample 90% of the data without replacement
        df_jack = df.sample(frac=0.9)
        
        # Non-linear Least Squares optimization
        initial_guess = [D0_NOMINAL, NOMINAL_DIRECTION]
        res = least_squares(rrt_residual_function, initial_guess, 
                            args=(df_jack['ra'], df_jack['z'], df_jack['mag_res']))
        
        d0_results.append(res.x[0])
        theta0_results.append(res.x[1] % 360)
        
        if (i+1) % 10 == 0: 
            print(f"   Iteration {i+1}/{n_iterations} complete.")

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

    # 4. Visualization: Parameter Dispersion Map
    plt.figure(figsize=(9, 6))
    plt.scatter(theta0_results, d0_results, alpha=0.6, color='#3498db', edgecolor='#2980b9')
    plt.axvline(theta_mean, color='#e74c3c', linestyle='--', label='Mean Directional Phase')
    
    plt.xlabel('Axis Direction (theta0) [Degrees]', fontweight='bold')
    plt.ylabel('Signal Intensity (D0)', fontweight='bold')
    plt.title('RRT Jackknife Stability: Parameter Covariance (SDSS DR16Q)', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.2)
    
    plt.savefig("rrt_jackknife_stability_plot.png", dpi=300)
    print("-> Stability plot saved: rrt_jackknife_stability_plot.png")
    plt.show()

if __name__ == "__main__":
    run_jackknife_stability_test(DATA_FILE)