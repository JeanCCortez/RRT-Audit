# Teoria da Relatividade Referencial (TRR) - Repositório de Auditoria Científica
# Referential Relativity Theory (RRT) - Scientific Audit Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19058853.svg)](https://doi.org/10.5281/zenodo.19058853)

---

## 🚀 Motor Cosmológico TRR / RRT Cosmological Engine (Interactive Audit)
Para facilitar a auditoria imediata sem necessidade de ambiente Python local, disponibilizamos o Motor TRR (Streamlit App).
*To facilitate immediate auditing without the need for a local Python environment, we provide the RRT Engine (Streamlit App).*

**🔗 Acesso / Access:** [https://rrt-motor.streamlit.app/](https://rrt-motor.streamlit.app/)

**Função / Function:** Validação de Dinâmica Galáctica (SPARC), Óptica Cosmológica Integrada, Predição Cega de Redshift (SLACS) e Predição Determinística de Ruptura em Correntes Estelares (Gaia). O motor emite relatórios técnicos de auditoria em PDF que quantificam as falhas matemáticas do modelo $\Lambda$CDM baseadas puramente na viscosidade do espaço-tempo.
*Validation of Galactic Dynamics (SPARC), Integrated Cosmological Optics, Blind Redshift Prediction (SLACS), and Deterministic Stream Rupture Prediction (Gaia). The engine generates technical PDF audit reports quantifying the mathematical failures of the $\Lambda$CDM model based purely on spacetime viscosity.*

---

## 📖 Descrição da Obra / Work Description

### 🇧🇷 Português
Este repositório contém a infraestrutura computacional e os algoritmos de auditoria empírica rigorosa utilizados para validar a **Teoria da Relatividade Referencial (TRR)**. A TRR propõe uma reformulação hidrodinâmica do espaço-tempo baseada em **Transições de Fase Termodinâmicas**. A teoria substitui entidades hipotéticas (Matéria e Energia Escuras) por um campo temporal viscoso ($\mathcal{T}_{\mu\nu}$) cuja interação com a matéria é estritamente governada pela densidade local de energia ($\rho$).

A tese está estruturada em **quatro volumes**, estabelecendo que o universo opera numa Escada de Transição de Fase (Phase Ladder) com diferentes regimes de viscosidade causal:
1. **Fase 1 (Saturada):** Regime de alta densidade (Sistema Solar, CERN) onde o fluxo causal sofre estagnação (0% de Ativação), blindando a TRR e recuperando perfeitamente a Relatividade Geral e o Modelo Padrão.
2. **Fase 2 (Transição):** Regime de densidade crítica (Halos Galácticos e Lentes Gravitacionais) onde a viscosidade bariônica ($\beta$) gera o arrasto que sustenta órbitas espirais e amplifica a deflexão óptica.
3. **Fase 3 (Viscosa):** Regime de vácuo profundo (Vazios Cósmicos) onde o fluxo temporal livre impulsiona a expansão e gera a Anisotropia Topológica do eixo de Cortez.

### 🇺🇸 English
This repository hosts the computational infrastructure and rigorous empirical audit algorithms used to validate the **Referential Relativity Theory (RRT)**. RRT proposes a hydrodynamic reformulation of spacetime based on **Thermodynamic Phase Transitions**. The theory replaces hypothetical entities (Dark Matter and Dark Energy) with a dynamic, viscous temporal vector field ($\mathcal{T}_{\mu\nu}$), whose interaction with baryonic matter is strictly governed by the local energy density ($\rho$).

The thesis is structured across **four volumes**, establishing that the universe operates on a Phase Transition Ladder with distinct regimes of causal viscosity:
1. **Phase 1 (Saturated):** High-density regime (Solar System, CERN) where causal flow stagnates (0% Activation), shielding RRT and perfectly recovering General Relativity and the Standard Model.
2. **Phase 2 (Transition):** Critical density regime (Galactic Halos and Gravitational Lenses) where baryonic viscosity ($\beta$) generates the drag that sustains spiral orbits and amplifies optical deflection.
3. **Phase 3 (Viscous):** Deep vacuum regime (Cosmic Voids) where free temporal flow drives expansion and generates the Topological Anisotropy of the Cortez axis.

---

## 📂 Organização dos Módulos / Module Organization

1. **Cosmology Supernovae (`/1_Cosmology_Supernovae`):**
    * Algoritmos de processamento do catálogo Pantheon+ para extração do gradiente causal via Hubble Detrending.

2. **Cosmology Quasars (`/2_Cosmology_Quasars`):**
    * Scripts de validação estratigráfica da Anisotropia Topológica no SDSS DR16Q (via Jackknife com Sigma Clipping orgânico) e cálculo da quebra de causalidade cronológica no crescimento de buracos negros.

3. **Orbital and Gravitational Mechanics (`/3_Orbital_and_Gravitational_Mechanics`):**
    * Testes de nulidade em ambientes de alta densidade (LAGEOS-2). Confirmação do Princípio da Neutralidade Bariônica (BNP) e validação do viés de fase na onda gravitacional GW170817.

4. **High Redshift Anomalies (`/4_High_Redshift_Anomalies`):**
    * Resolução do paradoxo das Galáxias Impossíveis do James Webb (JWST) através da equação de Dilatação Causal da TRR.

---

## ⚠️ Declaração de Disponibilidade de Dados / Data Availability Statement

Para garantir a **reprodutibilidade independente** e absoluta transparência científica, este projeto não hospeda os bancos de dados astronômicos brutos em seu repositório. Pesquisadores devem baixar os catálogos oficiais inalterados diretamente de suas fontes acadêmicas e alocá-los nas respectivas pastas dos scripts.
*To ensure **independent reproducibility** and absolute scientific transparency, this project does not host raw astronomical databases in its repository. Researchers must download the unaltered official catalogs directly from their academic sources and place them in the respective script folders.*

**Fontes Oficiais Utilizadas pelos Scripts / Official Sources Used by Scripts:**
1. **Pantheon+SH0ES:** [GitHub Oficial](https://github.com/PantheonPlusSH0ES/Data_Release)
2. **SDSS DR16Q:** [SDSS eBOSS Algorithms](https://www.sdss.org/dr16/algorithms/qso_catalog/)
3. **SPARC Database:** [Case Western Reserve University](http://astroweb.cwru.edu/SPARC/)
4. **SLACS Lens Survey:** [Sloan Lens ACS Survey](https://www.slacs.org/)
5. **Gaia Stellar Streams:** [ESA/Gaia Consortium](https://www.cosmos.esa.int/web/gaia/data-releases)
6. **Planck Legacy Archive (ESA):** [ESA Planck Public Data](https://pla.esac.esa.int/)
7. **LAGEOS-2 (ILRS):** [International Laser Ranging Service](https://ilrs.gsfc.nasa.gov/) (Arquivos `.sp3`)
8. **CMS/CERN Open Data:** [CERN Open Data Portal](https://opendata.cern.ch/)
9. **JWST Early Release Science:** [MAST Portal (STScI)](https://mast.stsci.edu/)
10. **LIGO/Virgo/KAGRA (GWTC):** [GWOSC](https://gwosc.org/)

---

## 📋 Tabela de Scripts e Evidências Principais / Main Scripts & Evidence Table

| Script/Module | Alvo / Target | Fase (Regime) | Resultado / Result |
| :--- | :--- | :--- | :--- |
| `trr_sdss_dr16q_topological_audit.py` | SDSS DR16Q | **Fase 3** | **30.36σ (Anisotropy)** |
| `trr_pantheon_plus_full_audit.py` | Pantheon+ | **Fase 3** | **23.24σ (Gradient)** |
| `trr_sparc_rotation_curves.py` | SPARC | **Fase 2** | **1.33% Error (Residual)** |
| `Motor TRR (Aba Redshift)` | Lentes SLACS | **Fase 2** | **Convergiu sem Matéria Escura** |
| `trr_lageos2_pnb_shielding_audit.py` | LAGEOS-2 | **Fase 1** | **0% Activation (Shielded)** |
| `5-trr_blackhole_growth_causality_audit.py` | Quasars $z > 5$ | **Fase 3** | **100% Failure Rate in $\Lambda$CDM (3,978 violations)** |
| `trr_jwst_cosmos_audit.py` | Galáxias JWST | **Fase 3** | **Anomaly Resolved via Causal Dilation** |
| `9_rrt_gwtc4_anisotropy_audit.py` | GWTC-4 (LIGO-Virgo) | **Fase 3** | **3.5% Directional Mass Inflation** |

**🔥 Latest Falsification: The GWTC-4 Algorithmic Illusion**
The RRT Cosmological Engine recently audited the `35 M_sun` black hole peak reported in the GWTC-4 catalog. The script `9_rrt_gwtc4_anisotropy_audit.py` proved that this "new family" of black holes is a directional artifact caused by spacetime friction along the Cortez Axis:
* **Anti-Axis (Alignment -0.9):** 0.1% Algorithmic Error (Mass is real).
* **Cortez Axis (Alignment +0.9):** 3.5% Algorithmic Error (Mass artificially inflated by LIGO pipeline to hide distance fatigue).

> **Nota de Auditoria:** Os resultados preditivos para Lentes Gravitacionais (SLACS) dispensam o uso de halo escuro e exigem estritamente a "Massa Bariônica Total" como matriz de arrasto fluido, provando a ruptura da Relatividade Geral Clássica na deflexão óptica. O Motor RRT consolida essas provas de forma interativa.

---

## 📚 Citação e Indexação Acadêmica / Citation and Indexing

Para citar a Teoria da Relatividade Referencial (TRR), os dados de auditoria associados ou o código-fonte em artigos acadêmicos, por favor utilize o DOI oficial registrado via Zenodo. Isso garante o rastreamento adequado no Google Scholar e outras bases de dados.
*To cite the Referential Relativity Theory (RRT), the associated audit data, or the source code in academic papers, please use the official DOI registered via Zenodo. This ensures proper tracking on Google Scholar and other databases.*

**BibTeX:**
```bibtex
@misc{cortez_2026_rrt,
  author       = {Cortez, Jean Coutinho},
  title        = {Referential Relativity Theory (RRT) - Scientific Audit Repository},
  month        = mar,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19058853},
  url          = {[https://doi.org/10.5281/zenodo.19058853](https://doi.org/10.5281/zenodo.19058853)}
}
