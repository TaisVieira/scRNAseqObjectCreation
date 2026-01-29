# scRNAseqObjectCreation
Scripts to create a Seurat object and analyze it.

### objectCreation.R
Uses the expression matrix outputed in the [scRNAseqProcessPipeline repository](https://github.com/TaisVieira/scRNAseqProcessPipeline) to create a Seurat object. 

### runScrublet.py
Runs Scrublet on the filtered matrix from `objectCreation`. Possibility of setting a manual threshold available.
