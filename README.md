\# RNA-seq exploratory analysis (airway dataset)



Small project to practice RNA-seq QC and PCA in R using the airway dataset.



\## Contents



\- `data/`: airway\_scaledcounts.csv (RNA-seq counts) and airway\_metadata.csv (sample metadata)

\- `scripts/01\_qc\_airway.R`: R script for:

&nbsp; - loading counts and metadata

&nbsp; - library size QC

&nbsp; - log2 transformation

&nbsp; - PCA coloured by treatment and cell type

\- `results/`:

&nbsp; - `library\_sizes\_airway.png`

&nbsp; - `pca\_airway.png`

&nbsp; - `airway\_qc\_objects.rds`

