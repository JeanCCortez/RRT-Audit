import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# CONFIGURAÇÃO TRR
EIXO_CORTEZ_RA = 148.9
ARQUIVO = "asi.orb.lageos2.251220.v80.sp3"

def auditoria_lageos_v3():
    print(f"--- TRR: AUDITORIA GRAVITACIONAL LAGEOS-2 (ID: L52) ---")
    
    if not os.path.exists(ARQUIVO):
        print(f"ERRO: Arquivo {ARQUIVO} não encontrado.")
        return

    posicoes = []
    
    # Lendo o arquivo com o ID correto identificado na sua amostra
    with open(ARQUIVO, 'r', encoding='latin-1') as f:
        for linha in f:
            partes = linha.split()
            # Procuramos linhas que começam com P e identificam o satélite 52 ou L52
            if len(partes) >= 5 and partes[0].startswith('P'):
                if '52' in partes[0] or '52' in partes[1]:
                    try:
                        # O formato SP3: P ID X Y Z
                        # Se o ID estiver grudado no P (ex: PL52), o X é partes[1]
                        # Se houver espaço (ex: P L52), o X é partes[2]
                        idx_x = 1 if len(partes[0]) > 1 else 2
                        x = float(partes[idx_x])
                        y = float(partes[idx_x+1])
                        z = float(partes[idx_x+2])
                        posicoes.append([x, y, z])
                    except (ValueError, IndexError):
                        continue

    pos = np.array(posicoes)
    
    if len(pos) == 0:
        print("ERRO: Ainda não foi possível extrair coordenadas. Verifique o ID nas linhas de posição.")
        return

    print(f"Sucesso: {len(pos)} pontos de dados extraídos.")

    # 1. Geometria Causal (Ângulo de Ascensão Reta)
    ra_inst = np.degrees(np.arctan2(pos[:, 1], pos[:, 0])) % 360
    
    # 2. Resíduos de Energia Gravitacional (Variação do Raio)
    raios = np.linalg.norm(pos, axis=1)
    residuos_r = raios - np.mean(raios)
    
    # 3. Cálculo da Significância TRR
    alinhamento = np.cos(np.radians(ra_inst - EIXO_CORTEZ_RA))
    df = pd.DataFrame({'ra': ra_inst, 'residuo': residuos_r, 'alinhamento': alinhamento})
    
    r_obs = df['alinhamento'].corr(df['residuo'])
    sigma = abs(r_obs) * np.sqrt(len(df))

    print("\n" + "="*60)
    print(f"VEREDITO LAGEOS-2 (UNIFICAÇÃO MACRO): {sigma:.2f} SIGMA")
    print(f"CORRELAÇÃO COM EIXO 148.9°: {r_obs:.4f}")
    print(f"RESULTADO: {'QUEBRA DE ISOTROPIA CONFIRMADA' if sigma > 5 else 'ISOTROPIA DE EINSTEIN PREVALECE'}")
    print("="*60)

    # Gráfico de Dispersão para sua tese
    plt.figure(figsize=(10,6))
    plt.scatter(df['ra'], df['residuo'], c=df['residuo'], cmap='magma', s=5, alpha=0.5)
    plt.axvline(EIXO_CORTEZ_RA, color='cyan', linestyle='--', label='Eixo de Cortez')
    plt.title("Resíduos de Órbita LAGEOS-2 (NASA/ASI) vs Direção Sideral")
    plt.xlabel("Ascensão Reta (Graus)")
    plt.ylabel("Desvio Radial (km)")
    plt.colorbar(label='Amplitude do Resíduo')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    auditoria_lageos_v3()
