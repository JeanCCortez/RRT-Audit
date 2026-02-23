import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# Parâmetros Nominais da TRR
D0_NOMINAL = 0.794
OMEGA_P = 1128.0
DIRECAO_NOMINAL = 148.9

def to_native(array):
    if array.dtype.byteorder not in ('=', '|'):
        return array.byteswap().view(array.dtype.newbyteorder('='))
    return array

def residuo_trr(params, ra, z, mag_res):
    d0, theta0 = params
    # Modelo de Precessão de Cortez
    fase = (theta0 + (OMEGA_P / z)) % 360
    predicao = d0 * z * np.cos(np.radians(ra - fase))
    return predicao - mag_res

def executar_jackknife(caminho, n_cortes=50):
    print(f"Iniciando Teste Jackknife em {caminho}...")
    dat = Table.read(caminho, format='fits')
    
    # Extração e limpeza (Foco no estrato de ressonância z: 1.5 - 2.0)
    ra = to_native(np.array(dat['RA']))
    z = to_native(np.array(dat['Z']))
    mag_i = to_native(np.array(dat['PSFMAG'][:, 3]))
    
    df = pd.DataFrame({'ra': ra, 'z': z, 'mag': mag_i})
    df = df[(df['z'] >= 1.5) & (df['z'] <= 2.0) & (df['mag'] > 0)].copy()
    df['mag_res'] = df['mag'] - (5 * np.log10(df['z']))
    
    print(f"Amostra total: {len(df)} objetos.")
    
    # Listas para armazenar resultados das iterações
    d0_results = []
    theta0_results = []
    
    # Loop Jackknife: Remove 10% dos dados aleatoriamente em cada iteração
    for i in range(n_cortes):
        # Sampling: removemos uma fração dos dados
        df_jack = df.sample(frac=0.9)
        
        # Ajuste de Mínimos Quadrados
        x0 = [D0_NOMINAL, DIRECAO_NOMINAL]
        res = least_squares(residuo_trr, x0, args=(df_jack['ra'], df_jack['z'], df_jack['mag_res']))
        
        d0_results.append(res.x[0])
        theta0_results.append(res.x[1] % 360)
        
        if i % 10 == 0: print(f"Iteração {i}/{n_cortes} concluída...")

    # Estatística Final
    d0_mean, d0_std = np.mean(d0_results), np.std(d0_results)
    theta_mean, theta_std = np.mean(theta0_results), np.std(theta0_results)
    
    print(f"\n--- RELATÓRIO DE ESTABILIDADE JACKKNIFE ---")
    print(f"Coeficiente D0: {d0_mean:.4f} +/- {d0_std:.4f}")
    print(f"Direção Inicial theta0: {theta_mean:.2f}° +/- {theta_std:.2f}°")
    
    if theta_std < 2.0:
        print("VEREDITO: Sinal ALTAMENTE ESTÁVEL. Invariante a cortes de dados.")
    else:
        print("VEREDITO: Sinal sensível a outliers.")

    # Gráfico de Dispersão dos Parâmetros
    plt.figure(figsize=(8, 5))
    plt.scatter(theta0_results, d0_results, alpha=0.5, color='blue')
    plt.axvline(theta_mean, color='red', linestyle='--', label='Média Direcional')
    plt.xlabel('Direção do Eixo (theta0)')
    plt.ylabel('Intensidade (D0)')
    plt.title('Estabilidade Jackknife: TRR no SDSS DR16Q')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    executar_jackknife(r"C:\Users\JM\tese\novos_testes\DR16Q_Superset_v3.fits")
