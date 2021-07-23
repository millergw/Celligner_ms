library(taigr)
library(tidyverse)

# where to save the plots
save_output = '/Users/gmiller/Documents/Work/GitHub/celligner_tasks/celligner_runs/multi_tcga_depmap-dmc21q1_pdx_met500_juric/plots'

# load data
df <- taigr::load.from.taiga(data.name='multi-dataset-celligner-tcga-depmap-dmc-21q1-pdx-met500-juric--5001', data.version=NULL, data.file='Celligner_alignment')
biospecimens <- load.from.taiga(data.name='juric-rapid-autopsy-data-91ac', data.version=2, data.file='biospecimens') %>% 
  dplyr::rename(sampleID = collaborator_sample_id)

df <- df %>% dplyr::left_join(biospecimens)


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

#color by lineage and label by lineage
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(fill = lineage, size = type), color = 'white', alpha=0.6, pch=21)  +
  ggplot2::geom_point(data = df %>% filter(type == 'depmap-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=22, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=25, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=24, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'juric-tumor'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
  ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `juric-tumor` = 1.5, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Juric_Celligner_linlabs.png'), width = 6, height = 5, dpi = 300)


#color by lineage and label by subtype
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(fill = lineage, size = type), color = 'white', alpha=0.6, pch=21)  +
  ggplot2::geom_point(data = df %>% filter(type == 'depmap-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=22, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=25, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor'), aes(fill = lineage, size = type), color = 'black', alpha=0.6, pch=24, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'juric-tumor'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
  ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `juric-tumor` = 1.5, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggrepel::geom_text_repel(data = subtype_meds, aes(label = subtype), size = 2.0) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Juric_Celligner_subtypelabs.png'), width = 6, height = 5, dpi = 300)


# color by tissue_site and label by lineage
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(fill = tissue_site, size = type), color = 'white', alpha=0.6, pch=21)  +
  ggplot2::geom_point(data = df %>% filter(type == 'depmap-model'), aes(fill = tissue_site, size = type), color = 'black', alpha=0.6, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'novartisPDX-model'), aes(fill = tissue_site, size = type), color = 'black', alpha=0.6, pch=22, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'pediatricPDX-model'), aes(fill = tissue_site, size = type), color = 'black', alpha=0.6, pch=25, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'met500-tumor'), aes(fill = tissue_site, size = type), color = 'black', alpha=0.6, pch=24, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'juric-tumor'), aes(fill = tissue_site, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.5)  +
  ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `juric-tumor` = 1.5, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(size = FALSE) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Juric_Celligner_tissue_site_labs.png'), width = 6, height = 5, dpi = 300)

avail_lins <- df %>% 
  filter(type == 'juric-tumor', !is.na(lineage)) %>% 
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
  ggplot2::geom_point(data = df %>% filter(type == 'juric-tumor', lineage == cur_lin), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2)  +

  ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `juric-tumor` = 2, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
    
    
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, sprintf('Juric_Celligner_%s.png', cur_lin)), width = 6, height = 5, dpi = 300)

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
    ggplot2::geom_point(data = df %>% filter(type == 'juric-tumor', lineage == cur_lin), aes(fill = type, size = type, shape = type), color = 'black', alpha=0.8, stroke = 0.2)  +
    ggplot2::scale_shape_manual(values=c(`depmap-model`=21, `tcgaplus-tumor`=21, `juric-tumor` = 23, `novartisPDX-model` = 22, `pediatricPDX-model` = 25, `met500-tumor` = 24 )) +
    ggplot2::scale_size_manual(values=c(`depmap-model`=1.25, `tcgaplus-tumor`=0.75, `juric-tumor` = 2, `novartisPDX-model` = 1.25, `pediatricPDX-model` = 1.25, `met500-tumor` = 1.25 )) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'bottom', 
                   text=ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0)) +
    ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
    # ggplot2::guides(fill=FALSE, size = FALSE) +
    ggplot2::guides(size = FALSE) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") 
  ggsave(file.path(save_output, sprintf('Juric_Celligner_%s_with_legend.png', cur_lin)), width = 6, height = 5, dpi = 300)
  
}



# Just highlight the Juric samples (backgroun = TCGA only)
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'juric-tumor'), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tcgaplus-tumor`=0.75, `juric-tumor` = 2)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Juric_Celligner_highlight_juric_samples.png'), width = 6, height = 5, dpi = 300)

# Highlight the Juric samples and color by lineage (background = TCGA only)
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tcgaplus-tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'juric-tumor'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tcgaplus-tumor`=0.75, `juric-tumor` = 2)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(file.path(save_output, 'Juric_Celligner_highlight_juric_samples_linlab.png'), width = 6, height = 5, dpi = 300)

