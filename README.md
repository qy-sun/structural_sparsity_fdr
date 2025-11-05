# **Code for "Data-Adaptive Inference with FDR Control for Structural Sparsity Problems"**

Under review and available on [arXiv](https://arxiv.org/abs/xxxx.xxxxx)

---

This repository contains code for reproducing the results in the paper *"Data-Adaptive Inference with FDR Control for Structural Sparsity Problems"*. 

It implements a unified framework for data-adaptive inference with false discovery rate (FDR) control under structural sparsity, including:

* **TLasso: Transformational Lasso**
* **StatKnock: Stability Transformational Knockoff**
* **SATLasso: Sequential Adjusted Transformational Lasso**
* **Two-Stage "SATLasso + StatKnock" Procedure**

---

## **Description**

The repository is organized into the following main directories:

| Folder | Description |
|:--|:--|
| **function/** | Core implementation of the proposed methods, including TLasso, SATLasso, StatKnock, and SATLasso+StatKnock. |
| **function/other methods/** | Comparative baseline methods such as SplitKnockoff, GKnockoff, DS, and BY for structural sparsity problems. |
| **simulation/** | Simulation scripts for reproducing all experiments in the paper and its supplement. |
| **simulation/helper/** | Graph difference operators construction and corresponding coefficient generation utilities. |
| **real/** | Scripts for ADNI sMRI real-data analysis with 5 level hierarchical brain regions. |
| **real/lookup/** | Lookup tables for hierarchical brain regions and the corresponding edge endpoints in the brain region graph. |

Each main file begins with a detailed commented header specifying purpose, inputs, and outputs.  

---

## **Simulation Studies**

- `simulation/TLasso_exp.R` — Reproduces **Table 1** and **Figure 4** in the main text.  
- `simulation/SATLasso+StatKnock_exp.R` — Reproduces **Figure 5** in the main text.
- `simulation/StatKnock_exp.R` — Reproduces **Figure S.4** and **Figure S.5** in the main text.
- `simulation/SATLasso_exp.R` — Reproduces **Table S.1** and **Table S.2** in the main text.
---

## **Real Data Analysis**


> **Note:** Due to data-sharing and privacy restrictions, the ADNI dataset are not included in this repository.


---

## **Additional Notes**

Some experiments were executed on the Duke Compute Cluster (DCC) using parallel jobs to accelerate replication. Minor numerical discrepancies may occur when running locally due to randomization and parallel execution order.
