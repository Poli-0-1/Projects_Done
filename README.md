# Projects Portfolio  💻📓 


__The following repository will act as a centralized archive for some of the computational projects developed through these last few academic years, and hopefully the ones which will follow. The overall aim tìis to have a collection of projects, solutions, and implementations of computational approaches related to multiple fields of studies, ranging from simple bioinformatic analysis to the development of real-world application of computational tools. The objective is to denote not only some of my previous work but also my expertise development oiver time__




## Repositories 📁

Each repository holds one or multiple files related to a specific project for a certain subject of these last academic years. Multiple types of fle may be present acting as supplementary information or data needed for that specific task/assignment.




## Coding Projects

1. **Thesis Project**
   - _**Aim**_ --> Develop a computational tool for the prediction of differentiation outcome in hiPSC-CMs cellular cultures, focusing on risk stratification for biological long-term outcome to reduce resources, efforts, and time spent on cellular cultures that will likely either not differentiate in the target lineage, or at all.
   - **_Files_**:
      - `Cellpose_code_ver4.ipynb` → Heavier version of the pipeline developed and executed on **Google Colab** (GPU T4 required). **CellPose** is implemented as the main segmentation approach for fluorescence images of hiPSC-CMs, handling nuclei and cell body detection across 2D cultures fairly well. Feeds extracted morphological features into an GRU-based classifier for differentiation assessment, although due to lack of data the architecture could not be tested, the eventual architecture can be seen within the last sections of the **Testing_Perfectionism** script.
      - `Testing_Perfectionism.ipynb` → Lightweight local version of the pipeline. Replaces CellPose with a classical **Kitler-Illingworth / Otsu fallback thresholding mechanism** for image segmentation, keeping hardware requirements minimal. Implements the full workflow: interactive CLI (InquirerPy), TIFF time-series import, Image Quality Analysis (IQA — blur, luminosity, contrast, edge density, sharpness via Laplacian/Sobel/Canny/Tenengrad), adaptive preprocessing (bilateral filtering, adaptive thresholding for low-quality frames), watershed segmentation, morphological feature extraction (regionprops), LSTM model training, and binary classification output (Differentiating / Non-Differentiating). As for the network, due to lack of data and project-related time constraints, the architecture could not be tested.
   - **_Languages_** --> Python, Linux(very few projects and instances in which t was implemented).
   - **_Resources/Packages_** --> Keras, TensorFlow, Seaborn, pyplot, tqdm, scikit, sklearn, tifff, cv2, mahoas, InquirerPy, rich, pandas, scipy, skimage, matplotlib, cellpose, pathlib, collections, gc, torch, and other supplementary.
  

<br>
<br>


2. **Bioinfo Analysis**
   - ***Aim*** --> Collection of scripts and tutorial notebooks used to grasp the essential operations in bioinformatic analysis using R/RStudio. The repository covers pipelines spanning from R language fundamentals, to microarray and RNA-seq expression analysis, up to gene set enrichment theory and interpretation.
   - ***Files***:
     - `code.Rmd` → R Markdown introductory notebook covering R basics: data types, vectors, matrices, data frames, I/O, plotting, tidyverse operations, and basic statistical analysis using the `GeneModel.txt` dataset as reference.
     - `Expression_Analysis_Lab.R` → Microarray-based differential expression analysis pipeline on prostate cancer data (GEO dataset **GSE55945**), from CEL file import and RMA normalization to DEG detection using both Wilcoxon test and linear modeling.
     - `Expression_Analysis_Lab_RNAseq_Source.R` → RNA-seq differential expression pipeline on a hypoxia-treated prostate cancer cell line (GEO dataset **GSE106305**), covering raw count filtering, TMM normalization, and DEG analysis via generalized linear models.
     - `GENE_set_enrichment.R` → Conceptual and theoretical notes on Gene Set Enrichment Analysis (GSEA), covering Gene Ontology (GO) structure, MSigDB, enrichment scoring methods, and the Fisher test approach.
     - `GeneModel.txt` → Supplementary tab-delimited gene annotation table used in the `code.Rmd` tidyverse exercises.
   - ***Outputs*** → Boxplots (pre/post normalization), PCA plots (sample clustering), heatmaps (DEG visualization), MA plots, Volcano plots, hierarchical clustering dendrograms, exported DEG tables (`.txt`).
   - ***Languages*** → R, R Markdown.
   - ***Resources/Packages*** → `GEOquery`, `oligo`, `pd.hg.u133.plus.2`, `hgu133plus2.db`, `genefilter`, `limma`, `pheatmap`, `RColorBrewer`, `edgeR`, `GenomicFeatures`, `ggplot2`, `stringr`, `tidyverse`, `AnnotationDbi`, `MASS`.

<br>
<br>


3. **Smith-Waterman**
   - ***Aim*** --> Develop from scratch the Smith-Waterman algorithm used for scored sequences alignment, returning the optimal aligned sequences score-wise, with custom settings dedicated to quality filtering and consecutive matches threshold.
   - ***Files***:
     - `final_modified_SW_2_ver.py` → Final and fully functional version of the modified Smith-Waterman algorithm. Implements the classic SW local alignment pipeline (scoring matrix + traceback via arrow symbols ←↑\) with two custom extensions: a **threshold filter** (60% of the maximum alignment score) to retain only high-quality alignments, and a **three-consecutive-matches filter** that counts and validates runs of at least 3 uninterrupted matches. Additionally computes and prints the **longest streak** of consecutive matches. The full process runs in a loop, allowing the user to test multiple sequence pairs without restarting.
   - ***Outputs*** → Alignment scoring and traceack matrices alongside the optimally aligned sequence, highest score, longest strak of consecutive matches, sequences with 60% (or more) score, sequences presenting only 3-consecutive matches.
   - ***Languages*** → Python.
   - ***Resources/Packages*** → `pandas`, `enum`, `numpy`.




## Future Improvements ⏭️

Some of the future prospects for this page include completing and revising every single file/project, and assign to each one of them a complementary description of their objective, development, features, and possible implementations, corrections, improvements, and personal comments related to their structure and notes that could be useful for the enhancing/fixing of said codes.


## Contacts 📱

- [Linkedin](www.linkedin.com/in/andrea-policano-455b68287)
- [email](a.policano2001@gmail.com)
