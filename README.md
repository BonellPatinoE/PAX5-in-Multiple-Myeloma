# PAX5-in-Multiple-Myeloma
PAX5-in-Myeloma RNAseq data analysis about PAX5 expression in primary samples from coMMpass (MMRF) IA19.

Prior bioinformatic analysis performed at Wiita Lab at the University of California, San Francisco (UCSF) on RNAseq data showed higher CD70 expression in relapsed myeloma patients. Based on a machine learning model run thereafter, we found that one of the negative regulators of CD70 is PAX5. Here we analyzed PAX5 expression in RNAseq data from coMMpass dataset (IA19) in CD138+ cells from primary samples.

This analysis was performed based on PAX5 expression based on well-characterized tumor genotypes and showed that several myeloma karyotypic subtypes, specifically those with t(4;14), gain of Chr1q, and loss of Chr1p, showed significantly downregulated CD70 expression compared to all others in CoMMpass dataset, inversely what it was found for CD70. However, PAX5, a TF closely linked to B-cell identity, is most highly expressed in the t(11;14) subtype of myeloma. This suggests that PAX5 is a negative regulator at least at the transcriptome level, and probably close to a B-cell lymphoma phenotype and cytogenetic profile.

This analysis was performed using DESeq2 (v1.34.0) package in R (v4.1.2) as well, using 'apeglm' for LFC shrinkage and “gprofiler” (v0.2.2) for gene list functional enrichment analysis and namespace conversion.
