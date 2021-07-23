library(tidyverse)
library(taigr)

##### Create Juric metadata file for Celligner

# NOTE: we had Brian manually annotate the `participants` Taiga file 
# to add a "Celligner_subtype" column. This became the file he emailed Gwen and Javad, which I downloaded the second tab of and saved as `celligner-output_v27-classification-subtypes (Brian's Assignments).csv`
# participants <- load.from.taiga(data.name='juric-rapid-autopsy-data-91ac', data.version=4, data.file='participants')
brian_manual <- readr::read_csv("/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_tcga_depmap-dmc21q1_pdx_met500_juric/celligner-output_v27-classification-subtypes (Brian's Assignments).csv") %>% 
  dplyr::select(participant_id, Celligner_subtype)

juric_participant_ann <- brian_manual %>% 
  # readr::read_csv("/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_tcga_depmap-dmc21q1_pdx_met500_juric/juric_participant_metadata_fixed_disease.csv") %>% 
  tidyr::separate(col="Celligner_subtype", into = c('lineage','subtype'), sep=' - ', remove=T) %>% 
  # filter out the samples that Brian doesn't want us to use
  dplyr::filter(lineage != 'Do not use') %>%
  dplyr::mutate(subtype = case_when(subtype == 'all' ~ '', TRUE ~ subtype),
                type = 'juric-tumor') %>% 
  dplyr::select(participant_id,lineage, subtype, type)

# Load Juric data from Taiga
expression.tpm <- load.from.taiga(data.name='juric-rapid-autopsy-data-91ac', data.version=2, data.file='expression_tpm')
biospecimens <- load.from.taiga(data.name='juric-rapid-biopsy-data-87fd', data.version=2, data.file='biospecimens') %>% 
  dplyr::rename(sampleID = collaborator_sample_id)

# Generate an annotated participant file, like the `juric_celligner_metadata` file in Taiga
tmp_juric <- dplyr::left_join(juric_participant_ann, biospecimens, by='participant_id') %>% 
  # restrict to Juric samples that we have expression data for
  dplyr::filter(sampleID %in% (rownames(expression.tpm)))
tmp_juric <- tmp_juric %>% dplyr::select(participant_id, sampleID, lineage, subtype, everything())

# Save down the metadata. Before running in Celligner, we will need to double-check the subtype and lineage annotations. This is already built into my Celligner process.
readr::write_csv(tmp_juric, file = "/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_tcga_depmap-dmc21q1_pdx_met500_juric/juric_celligner_metadata.csv")
