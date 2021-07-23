# Subset CCLF expression data 
library(tidyverse)
expression.cclf.20201005.genes <- load.from.taiga(data.name='rnaseq-expression-data-47e5', data.version=1, data.file='expression.cclf_20201005.genes')
rna_id_mapping <- readr::read_csv("/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_with_cclf/rna_smids_terraids_map_full_characterization_only.csv")

rna_samples_to_use <- rna_id_mapping %>% 
  .[["sampleID"]]

length(rna_samples_to_use)
cclf_expression_subset <- subset(expression.cclf.20201005.genes, rownames(expression.cclf.20201005.genes) %in% rna_samples_to_use)

# 123 58676
dim(cclf_expression_subset)

write.csv(cclf_expression_subset, '/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_with_cclf/cclf_expression_sample_subset.csv')
# cclf <- read.csv('/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_with_cclf/cclf_expression_sample_subset.csv')
cclf <- data.table::fread('/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_with_cclf/cclf_expression_sample_subset.csv') %>% 
  as.data.frame() %>%
  tibble::column_to_rownames('V1') %>%
  as.matrix()
