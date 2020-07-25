library(KEGGREST)
library(pathfindR)

# Obtain list of A.thaliana pathways
ath_kegg_descriptions <- keggList("pathway", "ath")

# Turn the identifiers of KEGGREST into KEGG-style pathway identifiers
kegg_ids <- sub("path:", "", names(ath_kegg_descriptions))

# Obtain and parse genes per each pathway
ath_kegg_genes <- sapply(kegg_ids, function(pwid){
  pw <- keggGet(pwid)
  pw <- pw[[1]]$GENE[c(TRUE, FALSE)] # get gene symbols, not descriptions
  #pw <- sub(";.+", "", pw) # discard any remaining description
  pw <- pw[grep("^[A-Za-z0-9_-]+(\\@)?$", pw)] # remove mistaken lines that cannot be gene symbols
  pw <- unique(pw) # keep unique symbols
  pw
})

## Filter list to exclude terms with 0 genes (metabolic pathways)
ath_kegg_genes <- ath_kegg_genes[sapply(ath_kegg_genes, length) != 0]

## Form the custom descriptions vector
names(ath_kegg_descriptions) <- sub("path:", "", names(ath_kegg_descriptions))
ath_kegg_descriptions <- sub(" - Arabidopsis thaliana \\(thale cress\\)", "", ath_kegg_descriptions)
ath_kegg_descriptions <- ath_kegg_descriptions[names(ath_kegg_descriptions) %in% names(ath_kegg_genes)]

## Save both as RDS files for later use
saveRDS(ath_kegg_genes, "ath_kegg_genes.RDS")
saveRDS(ath_kegg_descriptions, "ath_kegg_descriptions.RDS")

input_df <- tt[rownames(tt) %in% rscudo_probes$x, 2:3]
rscudo_agi <- anno$AGI[rownames(exp) %in% rownames(input_df)]
input_df <- cbind(rscudo_agi, input_df)
input_df[, 1] <- sub(" /// .+", "", input_df[, 1])


output_df <- run_pathfindR(input = input_df,
                           convert2alias = FALSE,
                           gene_sets = "Custom",
                           custom_genes = ath_kegg_genes,
                           custom_descriptions = ath_kegg_descriptions,
                           pin_name_path = "C:/Users/aless/Downloads/at.pin.filt.sif")
clustered <- cluster_enriched_terms(output_df, use_description = TRUE)
term_gene_graph(output_df, use_description = TRUE)
enrichment_chart(output_df)
