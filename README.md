# Teoria da Relatividade Referencial (TRR) - Repositório de Auditoria Científica
# Referential Relativity Theory (RRT) - Scientific Audit Repository

---

## 🚀 Motor Cosmológico TRR / RRT Cosmological Engine (Interactive Audit)
Para facilitar a auditoria imediata sem necessidade de ambiente Python local, disponibilizamos o Motor TRR (Streamlit App).
*To facilitate immediate auditing without the need for a local Python environment, we provide the RRT Engine (Streamlit App).*

**🔗 Acesso / Access:** [https://rrt-motor.streamlit.app/](https://rrt-motor.streamlit.app/)

**Função / Function:** Validação de Dinâmica Galáctica (SPARC), Óptica Cosmológica Integrada, Predição Cega de Redshift (SLACS) e Predição Determinística de Ruptura em Correntes Estelares (Gaia). O motor emite relatórios técnicos de auditoria em PDF que quantificam as falhas matemáticas do modelo $\Lambda$CDM.
*Validation of Galactic Dynamics (SPARC), Integrated Cosmological Optics, Blind Redshift Prediction (SLACS), and Deterministic Stream Rupture Prediction (Gaia). The engine generates technical PDF audit reports quantifying the mathematical failures of the $\Lambda$CDM model.*

---

## Descrição da Obra / Work Description

### 🇧🇷 Português
Este repositório contém a infraestrutura computacional e os algoritmos de auditoria estatística utilizados para validar a **Teoria da Relatividade Referencial (TRR)**. A TRR propõe uma reformulação hidrodinâmica do espaço-tempo baseada em **Transições de Fase Termodinâmicas**. A teoria substitui entidades hipotéticas (Matéria e Energia Escuras) por um campo temporal viscoso ($\mathcal{T}_{\mu\nu}$) cuja interação com a matéria é governada pela densidade local de energia ($\rho$).

A tese está estruturada em **quatro volumes**, estabelecendo que o universo opera em regimes distintos de viscosidade causal:
1. **Fase 1 (Saturada):** Regime de alta densidade (Sistema Solar, CERN) onde a TRR é blindada, recuperando a Relatividade Geral e o Modelo Padrão.
2. **Fase 2 (Transição):** Regime de densidade crítica (Halos Galácticos e Lentes Gravitacionais) onde a viscosidade ($\beta$) gera o arrasto que sustenta órbitas e amplifica a deflexão óptica.
3. **Fase 3 (Viscosa):** Regime de vácuo profundo (Vazios Cósmicos) onde o fluxo temporal impulsiona a expansão acelerada e gera a Anisotropia Topológica.

### 🇺🇸 English
This repository hosts the computational infrastructure and statistical audit algorithms used to validate the **Referential Relativity Theory (RRT)**. RRT proposes a hydrodynamic reformulation of spacetime based on **Thermodynamic Phase Transitions**. The theory replaces hypothetical entities (Dark Matter and Dark Energy) with a dynamic, viscous temporal vector field ($\mathcal{T}_{\mu\nu}$), whose interaction with baryonic matter is strictly governed by the local energy density ($\rho$).

The thesis is structured across **four volumes**, establishing that the universe operates in distinct regimes of causal viscosity:
1. **Phase 1 (Saturated):** High-density regime (Solar System, CERN) where RRT is shielded, recovering General Relativity and the Standard Model.
2. **Phase 2 (Transition):** Critical density regime (Galactic Halos and Gravitational Lenses) where viscosity ($\beta$) generates the drag that sustains orbits and amplifies optical deflection.
3. **Phase 3 (Viscous):** Deep vacuum regime (Cosmic Voids) where temporal flow drives accelerated expansion and generates Topological Anisotropy.

---

## 📂 Organização dos Módulos / Module Organization

1. **Core Cosmological Audits (`/Core Cosmological Audits`):**
    * Algoritmos de processamento de grandes catálogos (SDSS DR16Q, Pantheon+, Planck) para extração de significância estatística e validação da Rotação de Cortez ($\omega_p$).

2. **Experimental & Robustness (`/Experimental & Robustness`):**
    * Testes de nulidade em ambientes de alta densidade (LAGEOS-2, CMS/CERN) e simulações de dinâmica galáctica (SPARC). Confirmação da **isotropia local** e da validade da Fase 1 (Saturação).

3. **Critical Falsification Tests (`/Critical Falsification Tests`):**
    * Algoritmos desenhados para testar os limites físicos do Modelo Padrão. Inclui testes de Causalidade de Eddington, Anisotropia de Ondas Gravitacionais e o **Oráculo Interativo de Anisotropia** (Eixo Cortez).

4. **Official Validation Reports (`/Official_Validation_Reports`):** *(NOVO)*
    * Repositório das **Predições Cegas ("Eclipse de 1919")** e calibrações geradas pelo Motor TRR. Contém os PDFs oficiais atestando a recuperação precisa de Redshifts (SLACS) e as coordenadas predatórias de ruptura (Gaps) em Correntes Estelares (Gaia).

---

## 💾 Declaração de Disponibilidade de Dados / Data Availability Statement

Para garantir a **reprodutibilidade independente**, este projeto utiliza exclusivamente dados públicos brutos de repositórios oficiais. Nenhum dado foi pré-processado manualmente para favorecer a teoria.
*To ensure **independent reproducibility**, this project exclusively uses raw public data from official repositories. No data was manually pre-processed to favor the theory.*

**Fontes Oficiais e Variáveis Extraídas / Official Sources and Extracted Variables:**
1. **Pantheon+SH0ES:** [GitHub Oficial](https://github.com/PantheonPlusSH0ES/Data_Release)
   * *Extracted:* Redshift ($z$), Distance Moduli ($\mu$), and Covariance Matrices.
2. **SDSS DR16Q:** [SDSS eBOSS Algorithms](https://www.sdss.org/dr16/algorithms/qso_catalog/)
   * *Extracted:* Right Ascension (RA), Declination (DEC), and Redshift ($z$) for mapping the $51.73\sigma$ $\mathcal{T}_{\mu\nu}$ phase gradient.
3. **SPARC Database:** [Case Western Reserve University](http://astroweb.cwru.edu/SPARC/)
   * *Extracted:* Radius ($Rad$), Observed Velocity ($V_{obs}$), Gas Velocity ($V_{gas}$), and Disk Velocity ($V_{disk}$) for hydrodynamic drag calculations.
4. **SLACS Lens Survey:** [Sloan Lens ACS Survey](https://www.slacs.org/)
   * *Extracted:* Salpeter IMF Total Mass, Source/Lens Redshifts.
5. **Gaia Stellar Streams:** [ESA/Gaia Consortium](https://www.cosmos.esa.int/web/gaia/data-releases)
   * *Extracted:* Astrometric Pericenters and Apocenters for viscous shear predictions.
6. **Planck Legacy Archive (ESA):** [ESA Planck Public Data](https://pla.esac.esa.int/)
   * *Extracted:* CMB temperature anisotropies and polarization data.
7. **LAGEOS-2 (ILRS):** [International Laser Ranging Service](https://ilrs.gsfc.nasa.gov/)
   * *Extracted:* Satellite Laser Ranging (SLR) orbital residuals for Phase 1 null-tests.
8. **CMS/CERN Open Data:** [CERN Open Data Portal](https://opendata.cern.ch/)
   * *Extracted:* High-density collision events for Phase 1 shielding validation.
9. **JWST Early Release Science:** [MAST Portal (STScI)](https://mast.stsci.edu/)
   * *Extracted:* Spectroscopic redshifts ($z > 5$) and structural maturity indicators.
10. **LIGO/Virgo/KAGRA (GWTC):** [GWOSC](https://gwosc.org/)
    * *Extracted:* Luminosity Distance ($D_L$), RA, and DEC from O4 alerts.

---

## 📋 Tabela de Scripts e Evidências Principais / Main Scripts & Evidence Table

| Script/Module | Alvo / Target | Fase (Regime) | Resultado / Result |
| :--- | :--- | :--- | :--- |
| `trr_sdss_dr16q_51sigma_audit.py` | SDSS DR16Q | **Fase 3** | **51.73σ (Anisotropy)** |
| `Motor TRR (Aba Redshift)` | Lentes SLACS | **Fase 2** | **Convergiu sem Matéria Escura ($\Delta z \le 0.13$)** |
| `Motor TRR (Aba Correntes)` | Correntes Gaia | **Fase 2** | **Predição de Coordenadas Exatas de Ruptura** |
| `trr_pantheon_plus_gradient.py` | Pantheon+ | **Fase 2/3** | **25.47σ (Gradient)** |
| `trr_sparc_rotation_curves.py` | SPARC | **Fase 2** | **1.33% Error (Residual)** |
| `trr_ruptura_cronologia.py` | Quasars $z > 5$ | **Fase 3** | **100% Causal Violation ($\Lambda$CDM)** |
| `trr_lageos_pnb_shielding.py` | LAGEOS-2 | **Fase 1** | **0.22σ (Shielded / Null)** |

> **Nota de Auditoria:** Os resultados preditivos para Lentes Gravitacionais (SLACS) dispensam o uso de halo escuro e exigem estritamente a "Massa Bariônica Total" (Salpeter) como matriz de arrasto fluido, provando a ruptura da Relatividade Geral Clássica na deflexão óptica.

---

### 🛠️ Requisitos Técnicos / Technical Requirements
Utilize **Python 3.11+** com as bibliotecas: `streamlit`, `numpy`, `scipy`, `pandas`, `astropy`, `matplotlib` e `fpdf`.

---
**Autor / Author:** Jean Coutinho Cortez  
**Local / Location:** Rio de Janeiro, Brasil 🇧🇷  
**Data / Date:** Fevereiro / February 2026
