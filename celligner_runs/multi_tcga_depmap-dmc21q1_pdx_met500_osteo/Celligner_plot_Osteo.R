library(taigr)
library(tidyverse)

# where to save the plots
save_output = '/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_tcga_depmap-dmc21q1_pdx_met500_osteo/plots'
dir.create(save_output, showWarnings = FALSE)


# load data
# df <- taigr::load.from.taiga(data.name='multi-dataset-celligner-tcga-depmap-dmc-21q1-pdx-met500-juric--5001', data.version=NULL, data.file='Celligner_alignment')
# biospecimens <- load.from.taiga(data.name='juric-rapid-autopsy-data-91ac', data.version=2, data.file='biospecimens') %>% 
#   dplyr::rename(sampleID = collaborator_sample_id)
df <- Celligner_res
df <- fread('/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_tcga_depmap-dmc21q1_pdx_met500_osteo/Celligner_multidataset_info.csv')
all_osteo_metadata <- readr::read_csv('/Users/gmiller/Documents/Work/GitHub/osteo/data/osteo_metadata_for_se_gene_exploration.csv')
  

# create a display name for plotting

depmap.sample.info <- load.from.taiga(data.name='dmc-21q1-0e11', data.version=29, data.file='sample_info')
df$display_name <- plyr::mapvalues(x = df$sampleID, 
                                   from = depmap.sample.info$DepMap_ID,
                                   to = depmap.sample.info$stripped_cell_line_name)

df$display_name <- plyr::mapvalues(x = df$display_name, 
                                   from = all_osteo_metadata$exp_sample_id,
                                   to = all_osteo_metadata$display_name)
min_tumor_samples_per_label <- 10

# group_by_lineage
lin_meds <- df %>%
  dplyr::filter(type == 'tcgaplus-tumor') %>%
  group_by(lineage) %>% 
  dplyr::summarise(UMAP_1 = median(UMAP_1),
                   UMAP_2 = median(UMAP_2),
                   n = n()) %>% 
  dplyr::filter(n >= min_tumor_samples_per_label)
subtype_meds <- df %>%
  dplyr::filter(type == 'tcgaplus-tumor') %>%
  group_by(subtype) %>% 
  dplyr::summarise(UMAP_1 = median(UMAP_1),
                   UMAP_2 = median(UMAP_2),
                   n = n()) %>% 
  dplyr::filter(n >= min_tumor_samples_per_label)


# relative sizes to use in plots
sizes <- c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `peddepCL-model` = 1.5, `peddepPDX-model` = 1.5, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )

#color by lineage and label by lineage
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(fill = lineage, size = type), color = 'white', alpha=0.6, pch=21)  +
  ggplot2::geom_point(data = df %>% filter(type == 'depmap-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=22, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=25, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=24, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepCL-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
  ggplot2::scale_size_manual(values=sizes) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Osteo_Celligner_linlabs.png'), width = 6, height = 5, dpi = 300)


#color by lineage and label by subtype
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(fill = lineage, size = type), color = 'white', alpha=0.6, pch=21)  +
  ggplot2::geom_point(data = df %>% filter(type == 'depmap-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=22, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=25, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=24, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepCL-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
  ggplot2::scale_size_manual(values=sizes) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggrepel::geom_text_repel(data = subtype_meds, aes(label = subtype), size = 2.0) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Osteo_Celligner_subtypelabs.png'), width = 6, height = 5, dpi = 300)


# TODO: this was originally "color by tissue_site and label by lineage" for Juric data. What is tissue_site? Metastasis, right?
# color by subtype and label by lineage
# ggplot2::ggplot(df, 
#                 ggplot2::aes(UMAP_1, UMAP_2)) +
#   ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(fill = subtype, size = type), color = 'white', alpha=0.6, pch=21)  +
#   ggplot2::geom_point(data = df %>% filter(type == 'depmap-model'), aes(fill = subtype, size = type), color = 'black', alpha=0.6, pch=21, stroke = 0.2)  +
#   ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model'), aes(fill = subtype, size = type), color = 'black', alpha=0.6, pch=22, stroke = 0.2)  +
#   ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model'), aes(fill = subtype, size = type), color = 'black', alpha=0.6, pch=25, stroke = 0.2)  +
#   ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor'), aes(fill = subtype, size = type), color = 'black', alpha=0.6, pch=24, stroke = 0.2)  +
#   ggplot2::geom_point(data = df %>% filter(type == 'peddepCL-model'), aes(fill = subtype, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
#   ggplot2::geom_point(data = df %>% filter(type == 'peddepPDX-model'), aes(fill = subtype, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
#   ggplot2::scale_size_manual(values=sizes) +
#   ggplot2::theme_classic() + 
#   ggplot2::theme(legend.position = 'bottom', 
#                  text=ggplot2::element_text(size=8),
#                  legend.margin =ggplot2::margin(0,0,0,0)) +
#   ggplot2::guides(size = FALSE) +
#   ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
#   ggplot2::xlab("UMAP 1") +
#   ggplot2::ylab("UMAP 2") 
# ggsave(file.path(save_output, 'Osteo_Celligner_subtype_labs.png'), width = 6, height = 5, dpi = 300)

avail_lins <- df %>% 
  dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model'), !is.na(lineage)) %>% 
  pull(lineage) %>% 
  unique()

for (cur_lin in avail_lins) {
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor', lineage == cur_lin), aes(size = type, fill = type), color = 'white', alpha=0.8, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'depmap-model', lineage == cur_lin), aes(fill = type, size = type), color = 'white', alpha=0.8, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model', lineage == cur_lin), aes(fill = type, size = type), color = 'white', alpha=0.8, pch=22, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model', lineage == cur_lin), aes(fill = type, size = type), color = 'white', alpha=0.8, pch=25, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor', lineage == cur_lin), aes(fill = type, size = type), color = 'white', alpha=0.8, pch=24, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepCL-model', lineage == cur_lin), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepPDX-model', lineage == cur_lin), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2)  +
  
  ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `peddepCL-model` = 2, `peddepPDX-model` = 2, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
    
    
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, sprintf('Osteo_Celligner_%s.png', cur_lin)), width = 6, height = 5, dpi = 300)

}


# Lineage specific plots with dset origin ("type") in legend
for (cur_lin in avail_lins) {
  ggplot2::ggplot(df, 
                  ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type, shape = type), fill = 'darkgray', color = 'white', alpha=0.6, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor', lineage == cur_lin), aes(size = type, fill = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'depmap-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'peddepCL-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'black', alpha=0.8, stroke = 0.2)  +
    ggplot2::geom_point(data = df %>% filter(type == 'peddepPDX-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'black', alpha=0.8, stroke = 0.2)  +
    ggplot2::scale_shape_manual(values=c(`depmap-model`=21, `tcgaplus-tumor`=21, `peddepCL-model` = 23, `peddepPDX-model` = 23, `novartisPDX-model` = 22, `pediatricPDX-model` = 25, `met500-tumor` = 24 )) +
    ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `peddepCL-model` = 2, `peddepPDX-model` = 2, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'bottom', 
                   text=ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0)) +
    ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
    # ggplot2::guides(fill=FALSE, size = FALSE) +
    ggplot2::guides(size = FALSE) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") 
  ggsave(file.path(save_output, sprintf('Osteo_Celligner_%s_with_legend.png', cur_lin)), width = 6, height = 5, dpi = 300)
  
}

tmp <- ggplot2::ggplot(df, 
                       ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type, shape = type), fill = 'darkgray', color = 'white', alpha=0.6, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor', lineage == cur_lin), aes(size = type, fill = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'depmap-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'white', alpha=0.8, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepCL-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'black', alpha=0.8, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'peddepPDX-model', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'black', alpha=0.8, stroke = 0.2)  +
  ggplot2::scale_shape_manual(values=c(`depmap-model`=21, `tcgaplus-tumor`=21, `peddepCL-model` = 23, `peddepPDX-model` = 23, `novartisPDX-model` = 22, `pediatricPDX-model` = 25, `met500-tumor` = 24 )) +
  ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `peddepCL-model` = 2, `peddepPDX-model` = 2, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  # ggplot2::guides(fill=FALSE, size = FALSE) +
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggplotly(tmp) 

# Just highlight the Osteo samples (background = TCGA only)
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'tcgaplus-tumor'), aes(size = type), fill = 'lightgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'depmap-model', subtype == 'osteosarcoma'), aes(fill = type, size = type), color = 'darkgray', alpha=0.8, pch=24, stroke = 0.2) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model')), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tcgaplus-tumor`=0.75, `depmap-model` = 1.25, `peddepCL-model` = 1.25, `peddepPDX-model` = 1.25)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Osteo_Celligner_highlight_osteo_samples.png'), width = 6, height = 5, dpi = 300)

# Highlight the Osteo samples and label by lineage (background = TCGA only)
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'depmap-model', subtype == 'osteosarcoma'), aes(fill = type, size = type), color = 'darkgray', alpha=0.8, pch=24, stroke = 0.2) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model')), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tcgaplus-tumor`=0.75, `depmap-model` = 1.25, `peddepCL-model` = 1.25, `peddepPDX-model` = 1.25)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Osteo_Celligner_highlight_osteo_samples_linlab.png'), width = 6, height = 5, dpi = 300)

# Highlight the Osteo samples and label the outlier sample (background = TCGA only)
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'depmap-model', subtype == 'osteosarcoma'), aes(fill = type, size = type), color = 'darkgray', alpha=0.8, pch=24, stroke = 0.2) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model')), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tcgaplus-tumor`=0.75, `depmap-model` = 1.25, `peddepCL-model` = 1.25, `peddepPDX-model` = 1.25)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = df %>% 
                             dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model', 'depmap-model'), 
                                           subtype == 'osteosarcoma'), aes(label = display_name), size = 2.5) +
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Osteo_Celligner_highlight_osteo_samples_label_outlier.png'), width = 6, height = 5, dpi = 300)

# Highlight the Osteo samples in zoomed-in plot and label the samples (background = TCGA only)
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'depmap-model', subtype == 'osteosarcoma'), aes(fill = type, size = type), color = 'darkgray', alpha=0.8, pch=24, stroke = 0.2) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model')), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tcgaplus-tumor`=0.75, `depmap-model` = 1.25, `peddepCL-model` = 1.25, `peddepPDX-model` = 1.25)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = df %>% 
                             dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model', 'depmap-model'), 
                                           subtype == 'osteosarcoma'), 
                           aes(label = display_name), size = 2.5, max.overlaps = Inf) +
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2")  +
  coord_cartesian(xlim = c(-11, -3), ylim = c(-7, -2))
ggsave(file.path(save_output, 'Osteo_Celligner_highlight_osteo_samples_label_individual_samples.png'), width = 6, height = 5, dpi = 300)


## Interactive version


interact <- ggplot2::ggplot(df, 
                            ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(label = display_name, size = type), fill = 'lightgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'depmap-model', subtype != 'osteosarcoma'), aes(label = display_name, size = type), fill = 'lightgray', color = 'darkgray', alpha=0.6, pch=24, stroke = 0.2) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'depmap-model', subtype == 'osteosarcoma'), aes(label = display_name, fill = type, size = type), color = 'darkgray', alpha=0.8, pch=24, stroke = 0.2) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type %in% c('peddepCL-model', 'peddepPDX-model')), aes(label = display_name, fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tcgaplus-tumor`=0.75, `depmap-model` = 1.25, `peddepCL-model` = 1.25, `peddepPDX-model` = 1.25)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::geom_text(data = lin_meds, aes(label = lineage), size = 2.5) +
  
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2")  +
  coord_cartesian(xlim = c(-11, -3), ylim = c(-7, -2))

ggplotly(interact)
