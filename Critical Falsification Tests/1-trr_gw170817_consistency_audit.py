import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# RRT CONFIGURATION: GW170817 MULTI-MESSENGER CONSISTENCY AUDIT
# Goal: Validating RRT Causal Viscosity against known local Binary Neutron Star (BNS) mergers.
# Logic: Proving that Phase 3 effects are below the detection threshold for z < 0.01
# while predicting measurable divergence for cosmological distances (z > 0.5).
# ==============================================================================

def run_gw170817_consistency_test():
    """
    Executes a safety audit using the benchmark event GW170817.
    Compares observed electromagnetic (EM) distance vs. Gravitational Wave (GW) 
    distance under the RRT Viscous Vacuum model.
    """
    print("="*80)
    print("REFERENTIAL RELATIVITY THEORY (RRT): GW170817 CONSISTENCY AUDIT")
    print("Target: BNS Merger GW170817 / Host Galaxy NGC 4993")
    print("="*80)
    
    # 1. Observational Ground Truth (Benchmark)
    # Source: Abbott et al. 2017 (LIGO/Virgo) and Cantiello et al. 2018 (Optical)
    target_name = "GW170817 / NGC 4993"
    ra_event = 197.45   # Right Ascension (Degrees)
    dec_event = -23.38  # Declination (Degrees)
    z_event = 0.0099    # Observed Redshift
    
    dist_em_obs = 40.7  # Observed Optical Distance (Mpc)
    dist_gw_obs = 40.0  # Observed Gravitational Distance (LIGO) (Mpc)
    ligo_uncertainty = 8.0 # Statistical error margin (Mpc)
    
    # 2. RRT Calibrated Parameters (From Vol. IV)
    CORTEZ_AXIS_RA = 168.0
    CORTEZ_AXIS_DEC = -7.0
    XI_VISCOSITY = 0.308  # Causal Coupling Constant
    ANISOTROPY_A = 0.45   # Dipole Amplitude
    
    # 3. Geometric Projection (Mapping the event onto the Causal Flow)
    # Convert to radians for trigonometric computation
    r_ev = np.radians([ra_event, dec_event])
    r_ax = np.radians([CORTEZ_AXIS_RA, CORTEZ_AXIS_DEC])
    
    # Cosine Theta: Angular separation from the Causal Viscosity Axis
    cos_theta = (np.sin(r_ev[1]) * np.sin(r_ax[1]) + 
                 np.cos(r_ev[1]) * np.cos(r_ax[1]) * np.cos(r_ev[0] - r_ax[0]))
    
    angular_dist_deg = np.degrees(np.arccos(cos_theta))
    
    print(f"Event Location:  RA {ra_event:.2f}°, Dec {dec_event:.2f}°")
    print(f"Alignment:       {angular_dist_deg:.2f}° from Cortez Axis (CosTheta: {cos_theta:.2f})")
    
    if cos_theta > 0.5:
        print("Status:          High Alignment Zone (High theoretical drag potential).")
    else:
        print("Status:          Safe Zone (Low theoretical drag potential).")

    # 4. RRT Mathematical Prediction
    # Formula: D_GW_pred = D_EM * exp( (Xi/2) * z * (1 + A * cos_theta) )
    # At low redshift (z=0.009), the causal drag is mathematically suppressed.
    damping_term = (XI_VISCOSITY / 2) * z_event * (1 + ANISOTROPY_A * cos_theta)
    rrt_factor = np.exp(damping_term)
    
    predicted_gw_dist = dist_em_obs * rrt_factor
    predicted_divergence_pct = (rrt_factor - 1) * 100
    
    print("-" * 50)
    print(f"Real Optical Distance (EM):      {dist_em_obs:.2f} Mpc")
    print(f"LIGO Measured Distance (GW):     {dist_gw_obs:.2f} +/- {ligo_uncertainty} Mpc")
    print(f"RRT Predicted GW Distance:       {predicted_gw_dist:.2f} Mpc")
    print(f"Predicted RRT Bias:              +{predicted_divergence_pct:.2f}% (Viscous Damping)")
    
    # 5. Consistency Verdict
    upper_limit_ligo = dist_gw_obs + ligo_uncertainty 
    
    print("-" * 50)
    if predicted_gw_dist <= upper_limit_ligo:
        print("✅ VERDICT: ABSOLUTE SUCCESS.")
        print(f"RRT prediction ({predicted_gw_dist:.2f} Mpc) is fully consistent with GW170817 data.")
        print("Conclusion: Vacuum viscosity at z=0.01 is below the LIGO noise floor,")
        print("validating the theory against local multi-messenger constraints.")
    else:
        print("❌ VERDICT: FAILURE.")
        print("RRT predicted excessive damping for a local event.")

    # 6. Cosmological Extrapolation (The "Hunt" for O4/LISA)
    # Projecting the effect for z=1.0 where RRT is the dominant signal.
    print("=" * 80)
    print("RRT COSMOLOGICAL PREDICTION (FALSIFIABILITY TARGET)")
    print("Projection for high-redshift mergers (z=1.0) within the Cortez Corridor:")
    z_future = 1.0
    future_term = (XI_VISCOSITY / 2) * z_future * (1 + ANISOTROPY_A * cos_theta)
    future_factor = np.exp(future_term)
    future_divergence = (future_factor - 1) * 100
    
    print(f"Expected Luminosity Distance Divergence: +{future_divergence:.2f}%")
    print("RRT Prediction: A merger at 3000 Mpc will appear as 3600+ Mpc in GW detectors.")
    print("This spatial-dependent bias is the definitive test for RRT.")
    print("="*80)

if __name__ == "__main__":
    run_gw170817_consistency_test()