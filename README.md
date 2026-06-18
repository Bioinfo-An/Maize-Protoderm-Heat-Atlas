# Single-cell and spatial multi-omics reveal a protoderm-driven survival trade-off during extreme heat stress in maize


## 🧬 Project Overview
By integrating high-resolution **snRNA-seq**, **hdWGCNA**, **DAP-seq**, and **In Situ Sequencing (ISS)**, we constructed a spatiotemporal expression atlas of maize seedling leaves under acute heat stress and sustained acclimation. We identified a *ZmHsftf11*-dictated "proteostatic checkpoint" specifically anchored in the protodermal layer, driving a systemic survival trade-off network.

## 📂 Repository Structure
The analysis workflow is modularized into the following R scripts:

* `01_snRNA_seq_QC_and_Clustering.R`: Quality control, DoubletFinder filtering, Harmony batch-effect correction, and unsupervised Louvain clustering (Seurat v4).
* `02_Pseudotime_and_CytoTRACE.R`: Cellular stemness inference and single-cell developmental trajectory reconstruction (Monocle 3).
* `03_hdWGCNA_and_GENIE3_Network.R`: High-dimensional weighted gene co-expression network analysis and directional gene regulatory network (GRN) inference.
* `04_Spatial_ISS_Quantification.R`: Spatial transcriptomic mapping and absolute quantification of core target genes.

## 💻 System Requirements
All computational analyses were performed utilizing **R version 4.2.0** or higher. 

**Core Dependencies:**
* `Seurat` (v4.3.0)
* `Harmony`
* `Monocle3`
* `hdWGCNA`
* `GENIE3`

## 📊 Data Availability
The raw sequencing data (snRNA-seq, bulk RNA-seq, and DAP-seq) generated in this study have been deposited in the Genome Sequence Archive (GSA) under accession number: **CRA043742**. 

## ✉️ Contact
For any questions regarding the code or data, please contact:
Ruilin An (Your Email Here) or Prof. Yufeng Jiang (Jiangyufeng@gxaas.net).
