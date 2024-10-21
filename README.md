# Bioinformatics Project: Single Cell RNA-seq Analysis with Seurat

## Summary
This project focuses on analyzing single-cell RNA-seq data using the Seurat package. While Seurat provides insights into the overall percentage of cells expressing a specific gene within a chosen cluster, it lacks information about gene expression percentages for each individual cell across different clusters. To address this gap, we developed a tool to present the results from analyses conducted with Seurat in a user-friendly manner.

The project was implemented in R and utilizes Seurat for quality control, analysis, and exploration of single-cell RNA-seq data. We created a script that processes data to compute the gene expression percentages exceeding specified thresholds for each cell within a target cluster. We validated the tool using data from a study on enhancers related to leukemia, focusing on the BCL11B gene.

## Introduction
Seurat allows for two types of data analyses: spatially resolved data and integrative multimodal analysis. The aim of this project is to elucidate cellular characteristics within different clusters and identify differentially expressed genes. Our tool provides a clear and interpretable way to access detailed gene expression information that is often obscured by raw data. By normalizing and organizing the information, users can better understand the relationship between gene expression levels and other data points.

## Methods
- Implemented in R using Seurat for single-cell RNA-seq analysis.
- Followed Seurat's tutorial to create a Seurat object from the dataset.
- Filtered cells based on various criteria (e.g., unique gene counts, total molecule counts).
- Normalized the data and employed PCA for dimensionality reduction.
- Clustered cells using KNN graphs and identified features that define each cluster.
- Developed functions to assess gene expression levels and visualize results.

## Results and Discussion
We created a script that organizes data and allows users to retrieve information on gene expression levels for specific clusters. Key functions include:
- `Cells_in_clusters`: Returns cells in specified clusters.
- `FindMarkers_New`: Identifies differentially expressed genes.
- `ExpressionGenes`: Provides gene expression data for specific clusters.

The tool was tested on data from the article "Enhancer Hijacking Drives Oncogenic BCL11B Expression in Lineage-Ambiguous Stem Cell Leukemia," which discusses the role of enhancers in leukemia and focuses on the BCL11B gene.

## Future Directions
Future work could enhance the tool's capabilities to analyze scRNA-seq data for cells with significant expression levels within clusters. This could allow the identification of recurring patterns under certain conditions.

## References
- [Seurat Documentation](https://satijalab.org/seurat/)
- [Enhancer Hijacking Drives Oncogenic BCL11B Expression](https://aacrjournals.org/cancerdiscovery/article/11/11/2846/666485/Enhancer-Hijacking-Drives-Oncogenic-BCL11B)
- [GeneCards for BCL11B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=BCL11B)

## Data
- Download the dataset used in this project from [10X Genomics](https://www.10xgenomics.com/).
- GEO Data: [GSE173320](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173320)

