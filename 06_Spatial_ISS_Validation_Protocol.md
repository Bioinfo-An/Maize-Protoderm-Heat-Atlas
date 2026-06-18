# Step 06: Spatial In Situ Sequencing (ISS) Validation Protocol

## 🔬 1. Overview of the ISS Pipeline
Unlike the snRNA-seq analysis which relied on custom R-based bioinformatic pipelines (Scripts `01` to `05`), the high-resolution spatial validation of critical survival trade-off genes (*ZmHsftf11*, *ZmPdk1*, etc.) within maize seedling leaf tissues was executed via a standardized, hardware-integrated imaging and cyclic decoding pipeline provided by **DynaSpatial (德运康瑞)**.

Because the raw multidimensional TIFF micro-images and automated combinatorial fluorescence decoding steps require specialized hardware workstation environments, we document the precise software parameters and workflow utilized via **inSituFocus** here to ensure absolute technical reproducibility.

## 💻 2. Image Acquisition and Signal Decoding
The automated acquisition of fluorescent signals and initial transcripts spot calling were performed on the DynaSpatial high-throughput spatial imaging platform.
* **Image Processing Software Platform:** inSituFocus (DynaSpatial)
* **Signal Decoding Algorithm:** Transcripts were decoded sequentially based on rolling circle amplification (RCA) and cyclic combinatorial barcode hybridization signals. Background tissue autofluorescence and morphological noise were computationally minimized using the built-in adaptive top-hat filtering and rolling-ball background subtraction modules within the **inSituFocus** suite.

## 🎯 3. ROI Selection and Absolute Quantification
To specifically validate the spatial coordinates and absolute transcript density of our target genes across the distinct cellular layers (specifically capturing the protoderm and underlying mesophyll boundary), the following strict graphical user interface (GUI) workflow was executed within **inSituFocus**:

1. **Cell Segmentation & Registration:** DAPI nucleic staining was utilized to define absolute nuclear boundaries, and morphological landmarks were used to clear the architectural orientation of the maize seedling leaf mid-sections (specifically at the three-leaf stage).
2. **Region of Interest (ROI) Micro-annotation:** Cellular layers were strictly partitioned into distinct ROIs representing the Protoderm (Pro) layer and the adjacent underlying Epidermis/Mesophyll boundary. These tissue boundaries were traced consistently across both control (T0) and extreme heat-treated (T2 and T24) biological replicates.
3. **Absolute Spot Counting:** The **inSituFocus** engine automatically computed the absolute number of decoded, single-molecule fluorescent spots (representing individual mRNA transcripts) falling strictly within the coordinate masks of each defined cell-layer ROI.
4. **Data Export:** The raw, matrix-aligned spot density counts per cellular layer were exported directly as standard tabular formats for downstream structural visualization and bar plot rendering.

## 📊 4. Data Availability
The raw high-resolution cyclic fluorescence image stacks (`.tiff`) and the aligned spatial coordinate data tables output by the **inSituFocus** system have been fully archived and are available from the corresponding author upon reasonable request or through the project data repository.
