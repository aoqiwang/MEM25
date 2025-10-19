Based on transcriptome expression data, we analyzed the highly expressed genes and KEGG pathway enrichment of MEM25 under different salt treatments.
For transcriptome expression data under different conditions, genes with FPKM greater than 60 and high expression were selected for KEGG analysis, while enriched pathways with corrected P-values (p.adjust) less than 0.05 were retained for enrichment analysis.

These highly expressed genes are stored in these files according to their salt condition groups:
Genes_of_Ctrl-3h.txt; Genes_of_Ctrl-24h.txt; Genes_of_SS-3h.txt; Genes_of_SS-24h.txt

The results of their KEGG enrichment analysis are saved as KEGG_desult. txt, which contains detailed KEGG annotation information, including pathways and their corresponding gene sets, as well as statistical indicators.

As a new model species, Chlorella sp. MEM25, we have constructed the necessary files for KEGG annotation of MEM25, namely term2gene.txt and term2name.txt, which record the correspondence between each gene of MEM25 and the KEGG ontology. Enrichment analysis of MEM25 can be performed using R.code.
