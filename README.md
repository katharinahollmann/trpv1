# TRPV1_public
🧬 TRPV1 RNA-seq Analysis Scripts

Bulk RNA-seq data analysis and visualization pipeline
These R scripts were used to analyze murine RNA-seq datasets — including normalization, differential expression (DESeq2), PCA, and functional enrichment (clusterProfiler).

All scripts are annotated and commented to support understanding and reproducibility, especially for medical researchers without extensive programming experience.

⸻

📂 Contents
	•	DESeq2 normalization & DEG identification
For processing raw count data and identifying differentially expressed genes.

	•	PCA analysis, heatmap, VENN, beeswarm scripts
For visualizing sample clustering and group separation.

	•	Functional enrichment (GO & KEGG)
Using clusterProfiler and Enrichr for pathway interpretation.

	•	Visualization scripts
For generating publication-ready plots (volcano plots, dotplots, etc.).

⸻

⚙️ Requirements

R ≥ 4.3
Bioconductor packages: DESeq2, clusterProfiler, org.Mm.eg.db, enrichR, ggplot2, dplyr, tidyverse, etc.
(Exact package usage is documented in each script header, if published they need to be cited !! use citation())

⸻

⚠️ Disclaimer

These scripts are shared for educational and reference purposes only.
They were developed during a medical research project and are not guaranteed to be error-free or optimized for all datasets.

I just am a medical student — these scripts reflect my learning process.

⸻

💡 Citation / Contact

If you use or adapt these scripts for your own research, please cite the respective R packages and tools 
(e.g., Love et al., DESeq2, Genome Biology 2014; Yu et al., clusterProfiler, OMICS 2012; Chen et al., enrichR, BMC Bioinformatics 2013).

Suggestions and feedback are very welcome — open an issue or pull request.
