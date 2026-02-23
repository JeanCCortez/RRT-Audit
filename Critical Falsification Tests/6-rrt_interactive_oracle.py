"""
TRR Redshift Predictor (Interactive & Batch Mode)
=================================================
Calculates the anisotropic redshift prediction based on the Phase Drag 
equations of the Theory of Referential Relativity (TRR).

Usage:
    - Manual Mode: python trr_interactive_oracle.py
    - Batch Mode:  python trr_interactive_oracle.py catalog.csv
"""

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

# --- TRR CONSTANTS ---
CORTEZ_AXIS_L = 237.0  # Galactic Longitude of the Causal South Pole

def calculate_galactic_longitude(ra, dec):
    """Converts Equatorial coordinates (RA, DEC) to Galactic Longitude (l)."""
    coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    return coords.galactic.l.degree

def predict_trr_redshift(l_gal, z_base, amplitude):
    """Calculates the TRR Redshift incorporating vectorial phase drag."""
    delta_theta = np.radians(l_gal - CORTEZ_AXIS_L)
    z_trr = z_base - amplitude * np.cos(delta_theta)
    return z_trr

def interactive_manual_mode():
    """Mode 1: Manual Input Calculator for single targets."""
    print("\n" + "="*60)
    print(" TRR MANUAL CALCULATOR: SINGLE TARGET PREDICTION")
    print("="*60)
    print("Enter the parameters of the Supernova/Galaxy to calculate")
    print("the predicted redshift considering the causal viscosity.\n")
    
    try:
        ra = float(input("1. Right Ascension (RA) in degrees: "))
        dec = float(input("2. Declination (DEC) in degrees: "))
        z_base = float(input("3. Baseline Redshift (Isotropic Expectation): "))
        amp = float(input("4. Friction Amplitude D0 (e.g., 0.005): "))
        
        l_gal = calculate_galactic_longitude(ra, dec)
        z_trr = predict_trr_redshift(l_gal, z_base, amp)
        
        print("\n" + "-"*40)
        print(" RESULTS:")
        print("-"*40)
        print(f"Target Galactic Longitude (l): {l_gal:.3f} deg")
        print(f"Distance to Cortez Axis:       {abs(l_gal - CORTEZ_AXIS_L):.3f} deg")
        print(f"Standard Model Prediction:     z = {z_base:.5f}")
        print(f"TRR Anisotropic Prediction:    z = {z_trr:.5f}")
        print(f"Phase Drag Variation (Delta):  {z_trr - z_base:+.5f}")
        print("-"*40 + "\n")
        
    except ValueError:
        print("\n[ERROR] Invalid input. Please enter numerical values only.")

def interactive_batch_mode(file_path):
    """Mode 2: Safe CSV Processing with manual column mapping and visualization."""
    print("\n" + "="*70)
    print(f" TRR BATCH PROCESSING: {file_path}")
    print("="*70)
    
    try:
        # Robust reading for different CSV formats
        try:
            df = pd.read_csv(file_path, sep=None, engine='python')
        except:
            df = pd.read_csv(file_path, sep=r'\s+', engine='python')
            
        columns_list = list(df.columns)
        print("\nColumns detected in the catalog:")
        print(columns_list)
        print("\nPlease type the exact name of the columns to use (case-sensitive).")
        
        ra_col = input("Name of the RA column: ").strip()
        dec_col = input("Name of the DEC column: ").strip()
        z_col = input("Name of the REDSHIFT column: ").strip()
        
        # Validation
        for col in [ra_col, dec_col, z_col]:
            if col not in columns_list:
                print(f"\n[CRITICAL ERROR] Column '{col}' not found in the dataset.")
                sys.exit(1)
                
        try:
            amp = float(input("Enter Friction Amplitude D0 (e.g., 0.0048): "))
        except ValueError:
            print("[ERROR] Amplitude must be a number. Defaulting to 0.0048.")
            amp = 0.0048
            
        print("\n[PROCESSING] Calculating Galactic Coordinates...")
        df['TRR_L_GAL'] = calculate_galactic_longitude(df[ra_col].values, df[dec_col].values)
        
        print(f"[PROCESSING] Applying TRR Phase Drag Equation (D0 = {amp})...")
        df['TRR_PREDICTED_Z'] = predict_trr_redshift(df['TRR_L_GAL'], df[z_col], amp)
        df['TRR_DELTA_Z'] = df['TRR_PREDICTED_Z'] - df[z_col]
        
        output_file = "trr_enriched_catalog.csv"
        df.to_csv(output_file, index=False)
        
        print(f"\n[SUCCESS] Processing complete!")
        print(f"New catalog saved as: {output_file}")
        print("New columns added: TRR_L_GAL, TRR_PREDICTED_Z, TRR_DELTA_Z\n")

        # --- GRAPHICAL VISUALIZATION ---
        print("Generating diagnostic plot...")
        plt.figure(figsize=(10, 6))
        
        # Observed Data
        plt.scatter(df['TRR_L_GAL'], df[z_col], alpha=0.5, color='gray', edgecolor='black', label=f'Observed Data (N={len(df)})')
        
        # TRR Theoretical Curve
        z_mean = df[z_col].mean()
        l_range = np.linspace(0, 360, 100)
        z_trr_curve = predict_trr_redshift(l_range, z_mean, amp)
        
        plt.axhline(z_mean, color='blue', linestyle='--', linewidth=2, label=r'Isotropic Baseline ($\Lambda$CDM Expectation)')
        plt.plot(l_range, z_trr_curve, color='red', linewidth=3, label=f'TRR Phase Drag (D0 = {amp})')
        plt.axvline(CORTEZ_AXIS_L, color='darkred', linestyle=':', alpha=0.6, label='Cortez Axis (237 deg)')
        
        plt.title('TRR Predictive Anisotropy: Redshift vs Galactic Longitude', fontsize=14)
        plt.xlabel('Galactic Longitude $l$ (Degrees)', fontsize=12)
        plt.ylabel('Redshift ($z$)', fontsize=12)
        plt.legend(loc='best', framealpha=0.9)
        plt.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        
        plt.show()

    except FileNotFoundError:
        print(f"\n[ERROR] File '{file_path}' not found. Check the path and try again.")
    except Exception as e:
        print(f"\n[ERROR] An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        catalog_file = sys.argv[1]
        interactive_batch_mode(catalog_file)
    else:
        interactive_manual_mode()
