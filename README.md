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

The real data analysis in this work is based on the **Alzheimer's Disease Neuroimaging Initiative (ADNI)** structural MRI (sMRI) dataset, focusing on the Level 5 segmentation results derived from **MRICloud** processing. Due to data-sharing and privacy restrictions, the ADNI dataset are not included in this repository. In what follows, we provide a detailed description of the data acquisition, preprocessing, and integration workflow used to obtain the analysis-ready dataset.

- **ADNI data portal:** [https://adni.loni.usc.edu/data-samples/access-data/](https://adni.loni.usc.edu/data-samples/access-data/)  
- **MRICloud platform:** [https://mricloud.org/](https://mricloud.org/)  
- **MRIStudio platform:** [https://www.mristudio.org/installation.html](https://www.mristudio.org/installation.html)  
- **Reference tutorial:** [MRICloud T1 Tutorial](https://braingps.mricloud.org/docs/tutorials/mricloud.html)

### **Data Acquisition and Preprocessing Workflow**

The data preprocessing follows the official **MRICloud–T1 segmentation pipeline** with manual verification in **ROIEditor**. Below outlines the detailed workflow used to obtain the Level 5 regional volume data used in our analysis:

1. **Download from ADNI**  
   - Obtain sMRI data in **DICOM (dcm)** format from the ADNI portal.  
   - Select relevant metadata (e.g., subject ID, age, sex, diagnosis group, and cognitive assessments such as **MMSE**).  
   - Download and retain the associated `.csv` files that record acquisition information and assessment scores.  
   - Use 2D structural MRI scans for subsequent processing.

2. **Convert DICOM to Analyze Format**  
   - Use the conversion tool **`Dcm2Analyze_v3.exe`** (available via MRIStudio) to transform DICOM files into the **Image** and **.hdi** formats required by MRICloud.  
   - Confirm that file naming conventions remain consistent (RID or subject ID).

3. **Confirm Slice Type in ROIEditor**  
   - Open the converted images in **ROIEditor** (MRIStudio).  
   - Verify that the slice type is correct: *Sagittal*, *Axial*, or *Sagittal data converted to Axial*.  
   - In this study, the source files were **Sagittal** scans.

4. **Segmentation on MRICloud**  
   - Upload the preprocessed files to **MRICloud** for automatic segmentation.  
   - For batch uploads, combine up to **five subjects per .zip file** (no subfolders, and compressed using the system’s default “Compress to ZIP” function while avoiding any third-party compression software).  
   - Choose the appropriate **Slice Type** and **Atlas** version. Recommended: `Adult_286labels_10atlases_V5L` (latest version of M2_252).  
   - Submit and wait for processing to complete.  
   - After segmentation, download the result package and extract:  
     - `corrected_stats.txt` – contains regional quantitative measures (used as feature matrix `X`).  
     - `multilevel_lookup_table.txt` – hierarchical label mapping file (used to construct transformation matrices).

5. **Outcome Variable (Response `y`)**  
   - The cognitive outcome variable is the **ADNI_Mem** composite memory score from  
     **“UW - Neuropsych Summary Scores [ADNI1, GO, 2, 3]”**.  
   - This variable serves as a standardized cognitive performance indicator widely used in ADNI literature.  
   - Match each MRI sample to `ADNI_Mem` by aligning:
     - The last four digits of the file name with the **RID** field.  
     - The **EXAMDATE** with the imaging acquisition date.

6. **Integration**  
   - Combine `corrected_stats.txt` (regional measures) with the matched `ADNI_Mem` scores into a single dataset.  
   - Construct transformation matrices (`D1`, `D1`, `D3`) using the **lookup table**, which define the structural sparsity constraints used by the proposed SATLasso–StatKnock framework.

### **Remarks**

- The segmentation results can also be visualized directly on the **MRICloud** web interface or inspected locally using **ROIEditor**.  
- Due to privacy regulations, the processed ADNI data and lookup tables are **not publicly included** in this repository.  
- Reproduction of this section requires individual ADNI data access approval.

---

## **Additional Notes**

Some experiments were executed on the Duke Compute Cluster (DCC) using parallel jobs to accelerate replication. Minor numerical discrepancies may occur when running locally due to randomization and parallel execution order.
