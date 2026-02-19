import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ==============================================================================
# RRT CONFIGURATION: LAGEOS-2 GRAVITATIONAL SHIELDING AUDIT
# Target: LAGEOS-2 Passive Satellite (NASA/ASI Ephemerides)
# Goal: Validating the Baryonic Neutrality Principle (BNP) in Phase 1 (Saturated).
# Logic: In high-density energy regimes, RRT must recover General Relativity.
# ==============================================================================

# RRT Theoretical Reference Axis
CORTEZ_AXIS_RA = 148.9
# Standard SP3 format file (NASA/ASI)
DATA_FILE = "asi.orb.lageos2.251220.v80.sp3"

def run_lageos_shielding_audit():
    """
    Executes a gravitational null-test on LAGEOS-2 orbit residuals.
    This audit verifies if sidereal anisotropy affects local high-density baryonic 
    matter, testing the 'shielding' efficiency of the RRT Phase 1 regime.
    """
    print("="*80)
    print("REFERENTIAL RELATIVITY THEORY (RRT): LAGEOS-2 GRAVITATIONAL AUDIT")
    print("Experimental Focus: Baryonic Neutrality Principle (BNP) Validation")
    print("="*80)
    
    if not os.path.exists(DATA_FILE):
        print(f"CRITICAL ERROR: File {DATA_FILE} not found.")
        print("Please ensure the NASA SP3 ephemerides file is in the working directory.")
        return

    positions = []
    
    # Precise extraction of satellite coordinates from SP3 format
    print(f"-> Extracting orbital vectors from {DATA_FILE}...")
    with open(DATA_FILE, 'r', encoding='latin-1') as f:
        for line in f:
            parts = line.split()
            # SP3 lines starting with 'P' contain satellite position data
            if len(parts) >= 5 and parts[0].startswith('P'):
                # Identification of satellite ID 52 (LAGEOS-2)
                if '52' in parts[0] or '52' in parts[1]:
                    try:
                        # Handle varied SP3 spacing (PL52 vs P L52)
                        idx_x = 1 if len(parts[0]) > 1 else 2
                        x = float(parts[idx_x])
                        y = float(parts[idx_x+1])
                        z = float(parts[idx_x+2])
                        positions.append([x, y, z])
                    except (ValueError, IndexError):
                        continue

    pos_vectors = np.array(positions)
    
    if len(pos_vectors) == 0:
        print("ERROR: Coordinate extraction failed. Verify Satellite ID formatting.")
        return

    print(f"   Success: {len(pos_vectors)} data points ingested.")

    # 1. Causal Geometry (Instantaneous Right Ascension)
    # Mapping the satellite's sidereal position relative to the Cortez Axis
    ra_instantaneous = np.degrees(np.arctan2(pos_vectors[:, 1], pos_vectors[:, 0])) % 360
    
    # 2. Gravitational Energy Residuals (Radial Variance)
    # Analyzing deviations from the mean orbital radius
    radii = np.linalg.norm(pos_vectors, axis=1)
    radial_residuals = radii - np.mean(radii)
    
    # 3. RRT Significance Calculation
    # Testing for correlation between radial residuals and the RRT dipole orientation
    alignment_factor = np.cos(np.radians(ra_instantaneous - CORTEZ_AXIS_RA))
    
    # Statistical dataframe for correlation analysis
    audit_df = pd.DataFrame({
        'ra': ra_instantaneous, 
        'residual': radial_residuals, 
        'alignment': alignment_factor
    })
    
    pearson_r = audit_df['alignment'].corr(audit_df['residual'])
    sigma_level = abs(pearson_r) * np.sqrt(len(audit_df))

    print("\n" + "="*80)
    print(f"FINAL AUDIT VERDICT: {sigma_level:.2f} SIGMA")
    print(f"Correlation with Cortez Axis (148.9°): {pearson_r:.4f}")
    
    # Logic: If sigma < 5, Einstein's Isotropy is preserved (Success for RRT BNP)
    if sigma_level < 5:
        print("STATUS: PHASE 1 SHIELDING CONFIRMED (Einsteinian Isotropy Prevails).")
        print("Result: Baryonic matter is insensitive to vacuum viscosity in this regime.")
    else:
        print("STATUS: LOCAL ISOTROPY BREACH DETECTED (BNP Failure).")
    print("="*80)

    # Visualization: Residual Analysis vs. Sidereal Direction
    plt.figure(figsize=(10, 6))
    plt.scatter(audit_df['ra'], audit_df['residual'], 
                c=audit_df['residual'], cmap='viridis', s=5, alpha=0.4)
    plt.axvline(CORTEZ_AXIS_RA, color='#00FFFF', linestyle='--', label='Cortez Axis Alignment')
    
    plt.title("LAGEOS-2 Orbital Residuals vs. Sidereal Direction\nBNP Validation - Local Shielding Test", fontsize=12)
    plt.xlabel("Right Ascension (Degrees)")
    plt.ylabel("Radial Displacement (km)")
    plt.colorbar(label='Residual Amplitude')
    plt.legend()
    plt.grid(True, alpha=0.2)
    
    plt.savefig("rrt_lageos_shielding_test.png", dpi=300)
    print("-> Plot saved: rrt_lageos_shielding_test.png")
    plt.show()

if __name__ == "__main__":
    run_lageos_shielding_audit()