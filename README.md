# **Code for "Data-Adaptive Inference with FDR Control for Structural Sparsity Problems"**

Under review and available on [arXiv](https://arxiv.org/abs/xxxx.xxxxx)

---

This repository contains code for reproducing the results in the paper *"Data-Adaptive Inference with FDR Control for Structural Sparsity Problems"*. 

It implements a unified framework for data-adaptive inference with false discovery rate (FDR) control under structural sparsity, including:

* **TLasso: Transformational Lasso**
* **StatKnock: Stability Transformational Knockoff**
* **SATLasso: Sequential Adjusted Transformational Lasso**
* **Two-Stage SATLasso + StatKnock Procedure**

---

## **Description**

The repository is organized into the following main directories:

| Folder | Description |
|:--|:--|
| **function/** | Core implementation of the proposed methods, including TLasso, SATLasso, StatKnock, and SATLasso+StatKnock. |
| **function/other methods/** | Comparative baseline methods such as SplitKnockoff, GKnockoff, DS, and BY for structural sparsity problems. |
| **simulation/** | Simulation scripts for reproducing all experiments in the paper and its supplement. |
| **simulation/helper/** | Graph difference operators construction and corresponding $\beta$ generation utilities. |
| **real_data/** | Scripts for ADNI sMRI real-data analysis with 5 level hierarchical brain-region graphs. |

Each main file begins with a detailed commented header specifying purpose, inputs, and outputs.  

---

## **Simulation Studies**


---

## **Real Data Analysis**

---

## **Additional Notes**

Some experiments were executed on the Duke Compute Cluster (DCC) using parallel jobs to accelerate replication. Minor numerical discrepancies may occur when running locally due to randomization and parallel execution order. All random seeds are set (set.seed(2025)) for baseline reproducibility.
