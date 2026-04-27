## 📘 Overview  
Companion to  
**“Multiomics gene regulatory network analysis reveals the distinct roles of chromatin architectural factors in gene expression regulation”**

This repository provides code for:

- Gene regulatory network (GRN) construction  
- Network centrality analysis  
- Network modularity analysis  

---

## 🌐 Resources  
For the web tool associated with the paper, see: http://tfcentral.wjklab.com/

---

## 📂 Code Description  

### 🔹 Network Analysis
  - `centrality.R`  : Calculate seven network centrality measures  

  - `module.R`  : Identify network modules  



### 🔹 Multi-omics-based GRNs

- `pySCENIC.ipynb`  : Construct GRNs using pySCENIC  

- `direct_net.R`  : Construct GRNs using DIRECT-NET  

- `celloracle_step1.R`  : scATAC-seq analysis using Cicero  

- `celloracle_step2.ipynb`  : Preparation of base GRN input data  

- `celloracle_step3.ipynb`  : scRNA-seq preprocessing  

- `celloracle_step4.ipynb`  : GRN construction using CellOracle  



### 🔹 In silico Perturbation-based GRNs  

- `GENKI.ipynb`  : Build GRNs via virtual knockout using GENKI  

- `scTenifoldKnk.ipynb`  : Build GRNs via virtual knockout using scTenifoldKnk  

---

## 💡 Summary  
This project integrates multiple computational frameworks to systematically construct and analyze gene regulatory networks, with a focus on revealing the functional roles of chromatin architectural factors.
