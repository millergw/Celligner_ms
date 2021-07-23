library(taigr)
library(tidyverse)

min_tumor_samples_per_label <- 10
df <- load.from.taiga(data.name='celligner-runs-b3e5', data.version=2, data.file='aligned_data_with_metadata')

df <- df %>% dplyr::mutate(type = case_when(type == 'CL' ~ 'depmap-model',
                                     type == 'tumor' ~ 'tcgaplus-tumor',
                                     type == 'CCLF' ~ 'cclf-model',
                                     TRUE ~ type))
# group_by_lineage
lin_meds <- df %>%
  dplyr::filter(type == 'tumor') %>%
  group_by(lineage) %>% 
  dplyr::summarise(UMAP_1 = median(UMAP_1),
                   UMAP_2 = median(UMAP_2),
                   n = n()) %>% 
  dplyr::filter(n >= min_tumor_samples_per_label)
subtype_meds <- df %>%
  dplyr::filter(type == 'tumor') %>%
  group_by(subtype) %>% 
  dplyr::summarise(UMAP_1 = median(UMAP_1),
                   UMAP_2 = median(UMAP_2),
                   n = n()) %>% 
  dplyr::filter(n >= min_tumor_samples_per_label)

#color by lineage and label by lineage
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tumor'), aes(fill = lineage, size = type), color = 'white', alpha=0.8, pch=21)  +
  ggplot2::geom_point(data = df %>% filter(type == 'CL'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'CCLF'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=24, stroke = 0.2)  +
  ggplot2::scale_size_manual(values=c(`CL`=1.25, `tumor`=0.75, `CCLF` = 1.5)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave('~/Desktop/CCLF_Celligner_linlabs.png', width = 6, height = 5, dpi = 300)


#color by lineage and label by subtype
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tumor'), aes(fill = lineage, size = type), color = 'white', alpha=0.8, pch=21)  +
  ggplot2::geom_point(data = df %>% filter(type == 'CL'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'CCLF'), aes(fill = lineage, size = type), color = 'black', alpha=0.8, pch=24, stroke = 0.2)  +
  ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=0.75, `CCLF` = 2)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggplot2::guides(fill=FALSE, size = FALSE) +
  ggrepel::geom_text_repel(data = subtype_meds, aes(label = subtype), size = 2.5) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave('~/Desktop/CCLF_Celligner_subtypelabs.png', width = 6, height = 5, dpi = 300)


avail_lins <- df %>% 
  filter(type == 'CCLF', !is.na(lineage)) %>% 
  pull(lineage) %>% 
  unique()

for (cur_lin in avail_lins) {
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df %>% filter(type == 'tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.8, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'tumor', lineage == cur_lin), aes(size = type, fill = type), color = 'white', alpha=0.8, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% filter(type == 'CL', lineage == cur_lin), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=21, stroke = 0.2)  +
  ggplot2::geom_point(data = df %>% filter(type == 'CCLF', lineage == cur_lin), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=24, stroke = 0.2)  +
  ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1, `CCLF` = 2)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
    # ggplot2::guides(fill=FALSE, size = FALSE) +
    ggplot2::guides( size = FALSE) +
    ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
ggsave(sprintf('~/Desktop/CCLF_Celligner_%s.png', cur_lin), width = 6, height = 5, dpi = 300)
}

# Create zoomed-in plots
for (cur_lin in avail_lins) {
  ggplot2::ggplot(df, 
                  ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = df %>% filter(type == 'tumor'), aes(size = type), fill = 'darkgray', color = 'white', alpha=0.8, pch=21, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'tumor', lineage == cur_lin), aes(size = type, fill = type), color = 'white', alpha=0.8, pch=21, stroke = 0.1)  +
    ggplot2::geom_point(data = df %>% filter(type == 'CL', lineage == cur_lin), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=21, stroke = 0.2)  +
    ggplot2::geom_point(data = df %>% filter(type == 'CCLF', lineage == cur_lin), aes(fill = type, size = type), color = 'black', alpha=0.8, pch=24, stroke = 0.2)  +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1, `CCLF` = 2)) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'bottom', 
                   text=ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0)) +
    ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
    # ggplot2::guides(fill=FALSE, size = FALSE) +
    ggplot2::guides( size = FALSE) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") 
  ggsave(sprintf('~/Desktop/CCLF_Celligner_%s.png', cur_lin), width = 6, height = 5, dpi = 300)
}

# Color plot by data source ("type" column)
ggplot2::ggplot(df, 
                ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(data = df, aes(color = type)) +
  # ggplot2::geom_point(data = df %>% dplyr::filter(type == 'tumor'), aes(size = type), fill = 'lightgray', color = 'white', alpha=0.6, pch=21, stroke = 0.1)  +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'CL'), aes(fill = type, size = type), alpha=0.8, pch=24, stroke = 0.2) +
  ggplot2::geom_point(data = df %>% dplyr::filter(type == 'CCLF'), aes(fill = type, size = type), alpha=0.8, pch=23, stroke = 0.2) +
  ggplot2::scale_size_manual(values=c(`tumor`=0.75, `CL` = 1.25, `CCLF` = 1.25)) +
  ggplot2::theme_classic() + 
  ggplot2::theme(legend.position = 'bottom', 
                 text=ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0)) +
  
  ggrepel::geom_text_repel(data = lin_meds, aes(label = lineage), size = 2.5) +
  ggplot2::guides(size = FALSE) +
  ggplot2::xlab("UMAP 1") +
  ggplot2::ylab("UMAP 2") 
# ggsave(file.path(save_output, 'Osteo_Celligner_highlight_osteo_samples.png'), width = 6, height = 5, dpi = 300)