
library(shiny)
library(plotly)
library(dplyr)
library(ggplot2)
library(ggrepel)
# library(pheatmap)
library(shinythemes)
library(taigr)
library(readr)
library(shinyWidgets)
library(ggpubr)
library(DT)


# Helper functions --------------------------------------------------------

tumor_CL_dist <- taigr::load.from.taiga(data.name='celligner-output-f163', data.version=22, data.file='tumor_CL_distances')

alignment <- taigr::load.from.taiga(data.name='celligner-output-f163', data.version=22, data.file='Celligner_alignment') %>%
    dplyr::mutate(subtype = ifelse(is.na(subtype), 'NA', subtype),
                  `Primary/Metastasis` = ifelse(is.na(`Primary/Metastasis`), 'NA', `Primary/Metastasis`),
                  cluster = factor(cluster)) %>%
    dplyr::mutate(point_lab = sprintf("sample id: %s<br>sample type: %s<br>tissue: %s<br>subtype: %s<br>origin: %s<br>cluster: %s",
                                      sampleID, type, tissue, subtype, `Primary/Metastasis`, cluster))

tumor_classification_subtypes <- load.from.taiga(data.name='celligner-output-f163', data.version=22, data.file='classification_subtypes') %>% set_colnames(c('tissue', 'subtype'))

CFEs <- load.from.taiga(data.name='cfes-f6b9', data.version=1, data.file='CFEs')

# App ---------------------------------------------------------------------
ui <- fluidPage(theme=shinytheme("flatly"),
                
                # Application title
                titlePanel(h1("Celligner", h4("Computational alignment of tumor and cell line transcriptional profiles"))),
                # Sidebar with different panels for the different tabs
                sidebarPanel(conditionalPanel(condition="input.conditionedPanels == 'Explore the data'",
                                              position="left",
                                         multiInput("tissue","Select tissue site", choices=c("All tumor types", 'adrenal', 'bile_duct',
                                                                                                    'bone', 'brain', 'breast',
                                                                                                    'cervix', 'colorectal', 'endocrine', 'esophagus',
                                                                                                    'eye', 'fibroblast', 'gastric', 'germ_cell', 'kidney',
                                                                                                    'leukemia', 'liver', 'lung', 'lymphoma', 'mesothelioma', 'neuroblastoma',
                                                                                                    'ovary', 'pancreas', 'prostate', 'rhabdoid', 'rhabdomyosarcoma', 'skin',
                                                                                                    'soft_tissue', 'thyroid', 'upper_aerodigestive', 'urinary_tract',
                                                                                                    'uterus'), selected="All tumor types"),
                                         selectInput("feature_type", label = h5("Select feature to color by"), choices = c('tissue', 'subtype', 'Primary/Metastasis', 'type', 'cluster'), selected='tissue'),
                                         prettyCheckbox(inputId = "show_tumors", label = "Show tumor samples", icon = icon("check"), value = TRUE),
                                         prettyCheckbox(inputId = "show_cell_lines", label = "Show cell line samples", icon = icon("check"), value = TRUE),
                                         prettyCheckbox(inputId = "show_psite_labels", label = "Show tissue site labels", icon = icon("check"), value = TRUE),
                                         prettyCheckbox(inputId = "no_webGL", label = "High-resolution (slower)", icon = icon("check"), value = FALSE),
                                         h4("Highlight a sample"), position="left", selectizeInput("sample_ID", "Enter sample ID", choices='', selected=NULL),
                                         prettyCheckbox(inputId = "highlight_sample", label = "Highlight sample", icon = icon("check"), value = FALSE)),
                             conditionalPanel(condition="input.conditionedPanels == 'Find most similar tumors for a given cell line'", h3("Choose cell line"), position="left",
                                              selectizeInput("cur_CL","Enter cell line ID", choices='', selected=NULL),  numericInput("k_val", "select k nearest neighbors:", 25, min = 5, max = 50)),
                             conditionalPanel(condition="input.conditionedPanels == 'Rank cell lines for selected tumors'", h5("Select a tumor type/subtype to get a ranked list of the cell lines that are most similar to that tumor type"), 
                                              position='left',
                                              selectInput("cur_tissue","Select tissue site", choices=c('')),
                                              selectInput("cur_subtype", label = h5("Select subtype"), choices = c('')), 
                                              selectInput("color_by", label = h5("color plot by"), choices = c('selected tumors', 'cell line distances'))),
                             conditionalPanel(condition="input.conditionedPanels == 'CFEs'", h3("Choose a cancer functional event"), position="left",
                                              selectizeInput("CFE","select a CFE", choices='', selected=NULL)),
                             conditionalPanel(condition="input.conditionedPanels == 'Add your own data'", h4("Upload files", h6('upload your own data to align to the tumor-cell line Celligner projection')), 
                                              position='left',
                                              fileInput("exp_file",h5("1. upload expression file", 
                                                                      h6("expression file must a csv file of log2(TPM + 1) expression values, rows must be sample IDs (these should be the the first unamed column of the file), 
                                                                         and columns must be Ensembl gene ids")), 
                                                                      multiple = FALSE, accept = c("text/csv",
                                                                                                                           "text/comma-separated-values,text/plain",
                                                                                                                           ".csv")),
                                              fileInput("ann_file", h5("2. upload annotation file",
                                                                       h6("annotation file must be a csv file containing the columns sampleID (must match rownames of expression matrix), tissue, subtype, and type")),
                                                        multiple = FALSE, accept = c("text/csv",
                                                                                                                          "text/comma-separated-values,text/plain",
                                                                                                                          ".csv")),
                                              numericInput("mnn_k_val", h5("3. select k value for mutual nearest neighbors:",
                                                                           h6("used for MNN alignment step, should be between 5 and 50 and 
                                                                              approximately equal to the size of the smallest subtype group in the data")),
                                                           20, min = 5, max = 50),
                                              actionButton("runScript", "Run"))),
                                                #downloadButton("downloadData", "Download")
                
                mainPanel(
                  tabsetPanel(
                    tabPanel("Explore the data", plotlyOutput('explore_alignment', height='100%', width = '100%')),
                    tabPanel("Find most similar tumors for a given cell line", plotOutput('tumor_cl_dist_plot'),
                             plotOutput('class_hist'), DTOutput("tumor_ranking")),
                    tabPanel("Rank cell lines for selected tumors", plotlyOutput('alignment_class'), DTOutput('Data')),
                    tabPanel("CFEs", plotOutput('CFE_plot')),
                    tabPanel("Add your own data", plotlyOutput('new_data_plot')),
                    id='conditionedPanels')
                )
)

# Load in data based on user input
# create plots
server <- function(input, output, session) {
  options(shiny.maxRequestSize=80*1024^2)
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # Plotly variables
  t <- list(
    family = "sans serif",
    size = 11,
    color = toRGB("grey50"))
  t2 <- list(
    family = "sans serif",
    size = 11,
    color = toRGB("red"))

  #color palette
  tissue_colors <- c(`adrenal`="coral2", `bile_duct`="goldenrod", `bone`="darkorchid4", 
                           `brain`="hotpink3", `central_nervous_system`="hotpink3",`blood`="steelblue", `breast`="plum3", `cervix`="olivedrab3", `colorectal`="darkolivegreen3", 
                           `endocrine`="palevioletred", `esophagus`="steelblue2", `eye`="darkolivegreen", 
                           `fibroblast`="hotpink3", `gastric`="green4", `kidney`='chocolate2', 
                           `leukemia`="steelblue", `liver`="darkslategray4", `lung`="coral", `lymphoma`="lightseagreen",`lymphocyte`="lightseagreen", `mesothelioma`="chartreuse3", 
                            `nasopharynx`="gold4", `neuroblastoma`="violetred", `peripheral_nervous_system`="violetred", `ovary`="lightcoral", 
                           `pancreas`="tomato3", `pineal`='tomato2',`plasma_cell`="dodgerblue4",`prostate`="mediumpurple3", `rhabdoid`="tomato2", 
                           `rhabdomyosarcoma`="darkviolet", `skin`="yellowgreen", `soft_tissue`="lightblue3", `teratoma`="deeppink3",
                           `testis`="violetred1", `embyro`="violetred1", `germ_cell`="violetred1", `thymus`="turquoise3", `thyroid`="aquamarine4", `upper_aerodigestive`="slateblue", 
                           `urinary_tract`="palegreen3", `uterus`="skyblue3")
  
  #all_tissues <-  reactive({sort(unique(alignment$tissue))})
  classification_categories <-  reactive({as.character(tumor_classification_subtypes$tissue)})
  observe({
    updateSelectInput(session, "cur_tissue",
                      choices = classification_categories()
    )})
  
  ## MAIN PLOT ##
  all_tissues <-  reactive({c('All tumor types', as.character(sort(filter(as.data.frame(table(alignment$tissue)), Freq>2)$Var1)))})
  
  feature_types <-  reactive({colnames(alignment)[4:ncol(alignment)]})
  
  all_samples <-  reactive({unique(alignment$sampleID)})
  observe({
    updateSelectInput(session, "sample_ID",
                      choices = all_samples()
    )})
  
  zoom_output <- reactive({
    zoom <- event_data("plotly_relayout", "source1")
    
    # if plot just rendered, event_data is NULL
    # if user double clicks for autozoom, then zoom$xaxis.autorange is TRUE
    # if user resizes page, then zoom$width is pixels of plot width
    num_points <- nrow(alignment)
    xmin <- -12
    xmax <- 23
    ymin <- -12
    ymax <- 18
    if (is.null(zoom) || names(zoom[1]) %in% c("xaxis.autorange", "width") || "xaxis.showspikes" %in% names(zoom)) {
      #stick with defaults
    } else if (length(zoom) == 4) {
      xmin <- zoom$`xaxis.range[0]`
      xmax <- zoom$`xaxis.range[1]`
      ymin <- zoom$`yaxis.range[0]`
      ymax <- zoom$`yaxis.range[1]`
      num_points <- nrow(filter(alignment, UMAP_1 < xmax & UMAP_1 > xmin & UMAP_2 > ymin & UMAP_2 < ymax))
    } else if (length(zoom) == 2 & names(zoom[1]) == 'yaxis.range[0]') {
      ymin <- zoom$`yaxis.range[0]`
      ymax <- zoom$`yaxis.range[1]`
      num_points <- nrow(filter(alignment, UMAP_2 > ymin,UMAP_2 < ymax))
    } else if (length(zoom) == 2 & names(zoom[1]) == 'xaxis.range[0]') {
      xmin <- zoom$`xaxis.range[0]`
      xmax <- zoom$`xaxis.range[1]`
      num_points <- nrow(filter(alignment, UMAP_1 > xmin, UMAP_1 < xmax))
    }
    
    c(num_points, xmin, xmax, ymin, ymax)
    #list(`num_points` = num_points, `xmin`=xmin, `xmax`=xmax, `ymin`=ymin, `ymax`=ymax)
  })
  
  output$relayout <- renderPrint({event_data("plotly_relayout", "source1")})
  
  output$explore_alignment <- renderPlotly({
    validate(need(input$feature_type, message = FALSE), 
             need(input$tissue, message = FALSE),
             need(!is.null(input$show_psite_labels), message = FALSE))
    tumor_size <- 0.45
    CL_size <- 1.8
    xmin <- zoom_output()[2]
    xmax <- zoom_output()[3]
    ymin <- zoom_output()[4]
    ymax <- zoom_output()[5]
    num_pts <- zoom_output()[1]
    if(num_pts < 1000) {
      if(num_pts <= 50) {
        tumor_size <- 3
        CL_size <- 5
      } else {
        tumor_size <- 3.14 - (.0027*num_pts)
        CL_size <- 5.17 - (.0034*num_pts)
      }
    }
    if (input$show_psite_labels) {
      avgs <- alignment %>% 
          dplyr::filter(!is.na(tissue)) %>% 
          dplyr::group_by(tissue) %>% 
          dplyr::filter(!tissue %in% c('adrenal_cortex', 'central_nervous_system', 'small_intestine', 'embryo',
                                       'peripheral_nervous_system', 'blood', 'lymphocyte', 'teratoma', 'unknown')) %>% 
          dplyr::filter(!subtype %in% c('osteosarcoma')) %>% 
          dplyr::filter(!sampleID %in% alignment$sampleID[grep('engineered', alignment$tissue)]) %>%
          dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm=T),
                           UMAP_2 = median(UMAP_2, na.rm=T))
      lab_aes <- ggplot2::geom_text(data = avgs, mapping = aes(x = UMAP_1, y = UMAP_2, label = tissue), size = 3)
    }
    
    if('All tumor types' %in% input$tissue) {
      g <- ggplot(alignment, aes(UMAP_1, UMAP_2))
      if (input$show_tumors) {
        g <- g + geom_point(data = alignment %>% filter(type == 'tumor'),
                            aes_string(color = input$feature_type, text='point_lab'),
                            size = tumor_size, alpha = 0.6) +
          # scale_size_manual(values = c(CL = 2, tumor = 0.4)) +
          # scale_alpha_manual(values = c(CL = 0.75, tumor = 0.6))
          ggtitle("Alignment of tumors and cell lines") +
          theme(legend.position = "none")
      }
      if (input$show_cell_lines) {
        g <- g + geom_point(data = alignment %>% filter(type == 'CL'),
                   aes_string(fill = input$feature_type, text='point_lab'),
                   pch = 21, color = 'black', size = CL_size, alpha = 0.7, stroke = 0.3)
      }
    } else {
      cur_shown_tissues <- input$tissue
      if('adrenal' %in% cur_shown_tissues) {
        cur_shown_tissues <- c(cur_shown_tissues, "adrenal_cortex")
      }
      if('brain' %in% cur_shown_tissues) {
        cur_shown_tissues <- c(cur_shown_tissues, "central_nervous_system")
      }
      if('neuroblastoma' %in% cur_shown_tissues) {
        cur_shown_tissues <- c(cur_shown_tissues, "peripheral_nervous_system")
      } 
      if('lymphoma' %in% cur_shown_tissues) {
        cur_shown_tissues <- c(cur_shown_tissues, "lymphocyte")
      } 
      if('leukemia' %in% cur_shown_tissues) {
        cur_shown_tissues <- c(cur_shown_tissues, "blood")
      }
      if('germ_cell' %in% cur_shown_tissues) {
        cur_shown_tissues <- c(cur_shown_tissues, "embryo")
      }
      g <- ggplot(alignment, aes(UMAP_1, UMAP_2)) +
        geom_point(data = filter(alignment, !tissue %in% cur_shown_tissues), color = "#BFBFBF", alpha = 0.25, size = 0.25) 
      if (input$show_tumors) {
       g <- g +  geom_point(data = alignment %>% filter(type == 'tumor', tissue %in% cur_shown_tissues),
                   aes_string(color = input$feature_type, text='point_lab'),
                   size = tumor_size, alpha = 0.6) +
        # geom_point(data = filter(alignment, tissue %in% input$tissue),
        #            aes_string(color = input$feature_type, size = 'type'), alpha = 0.7) +
        #scale_size_manual(values = c(CL = 2, tumor = 0.4)) +
        ggtitle("Alignment of tumors and cell lines") + theme(legend.position = "none")
      }
      if (input$show_cell_lines) {
        g <- g + geom_point(data = alignment %>% filter(type == 'CL', tissue %in% cur_shown_tissues),
                            aes_string(fill = input$feature_type, text='point_lab'),
                            pch = 21, color = 'black', size = CL_size, alpha = 0.7, stroke = 0.3)
      }
      }
      if(input$highlight_sample) {g <- g + geom_point(data=filter(alignment, sampleID==input$sample_ID), color="#333333", size=4, aes_string(text='point_lab'))}
      if (input$show_psite_labels) {g <- g + lab_aes}
      if(input$feature_type=='tissue') {g <- g + scale_color_manual(values=tissue_colors) + scale_fill_manual(values = tissue_colors)}
      g <- g + scale_x_continuous(limits=c(xmin, xmax)) + scale_y_continuous(limits=c(ymin, ymax))
      p <- ggplotly(g, tooltip = c('text'), source="source1")  
      if (!input$no_webGL) {
        p %<>% plotly::toWebGL()
      }
      layout(p, margin=list(t=40), height=650, showlegend = FALSE)
    })
  
 

  all_CFEs <-  reactive({unique(CFEs$`CFE nodes`)})
  
  observe({
    updateSelectInput(session, "CFE",
                      choices = all_CFEs()
    )})
  
  
 ## Cell line nearest tumor neighbors plots ##
  
  all_cell_lines <-  reactive({unique(filter(alignment, type=='CL')$sampleID)})

  observe({
    updateSelectInput(session, "cur_CL",
                      choices = all_cell_lines()
    )})
  
    cur_neighbors <-  reactive({validate(need(input$cur_CL, message = FALSE), need(input$k_val, message = FALSE))
      rownames(tumor_CL_dist)[FastKNN::k.nearest.neighbors(which(colnames(tumor_CL_dist)==input$cur_CL),
                                                                                       t(tumor_CL_dist), k=input$k_val)]})

    tumor_label <- reactive({names(sort(table(filter(alignment, sampleID %in% cur_neighbors())$tissue), decreasing=T)[1])})


    output$tumor_cl_dist_plot <- renderPlot({
      validate(need(input$cur_CL, message = FALSE), need(input$k_val, message=FALSE))
          ggplot(alignment, aes(UMAP_1, UMAP_2, size=type)) +
                     geom_point(data = alignment %>% filter(!(sampleID %in% cur_neighbors()) & sampleID != input$cur_CL), color = 'gray', alpha=0.05) +
                     geom_point(data = alignment %>% filter(sampleID %in% cur_neighbors()),  aes(color=tissue), alpha=0.7) +
                     geom_point(data = alignment %>% filter(sampleID == input$cur_CL), color='black', alpha=0.9) +
                     scale_size_manual(values = c(`CL` = 2, `tumor` = 1.5)) +
                      scale_color_manual(values=tissue_colors) +
                     theme(legend.position = "right")})


    tumor_distances <- reactive({cbind.data.frame(`dist`=as.numeric(tumor_CL_dist[,which(colnames(tumor_CL_dist)==input$cur_CL)]),
                                                  `tissue`=alignment[which(alignment$sampleID %in% rownames(tumor_CL_dist)),'tissue']) %>%
    mutate(class=ifelse(tissue==tumor_label(), tumor_label(), 'other'))})

    output$class_hist <- renderPlot({validate(need(input$cur_CL, message = FALSE), need(input$k_val, message=FALSE))
      ggdensity(tumor_distances(), x='dist', y='..count..', fill='class')})

    output$tumor_ranking = renderDT({
      validate(need(input$cur_CL, message = FALSE), need(input$k_val, message=FALSE))
      td <- tumor_distances()
      td$sampleID <- rownames(tumor_CL_dist)
      td <- td %>% arrange(dist) %>% dplyr::select(dist, tissue, sampleID)
      td <- left_join(td, alignment[,c('sampleID', 'subtype', 'Primary/Metastasis')], by='sampleID')
      td[1:100,]
    })
    
    
    ## cell line ranking by tumor type ##
    
    subtype_categories <-  reactive({as.character(filter(tumor_classification_subtypes, tissue == input$cur_tissue)$subtype)})
    
    observe({
      updateSelectInput(session, "cur_subtype",
                        choices = subtype_categories()
      )})
    
    output$alignment_class <- renderPlotly({
      validate(need(input$cur_tissue, message = FALSE))
      validate(need(input$cur_subtype, message = FALSE))
      
      cur_used_subtypes <- input$cur_subtype
      if(cur_used_subtypes == "all") {
        cur_used_subtypes <- unique(filter(alignment, tissue == input$cur_tissue)$subtype)
      }
      
     selected_tumors <- filter(alignment, type=='tumor' & tissue == input$cur_tissue & subtype %in% cur_used_subtypes)$sampleID
      
      CL_dist <- apply(tumor_CL_dist[selected_tumors,,drop=F], 2, FUN=median)

      if(input$color_by == 'selected tumors') {
        g <- ggplotly(ggplot(alignment, aes(UMAP_1, UMAP_2, size=type, color=tissue == input$cur_tissue & subtype %in% cur_used_subtypes)) + geom_point() +
                        scale_size_manual(values = c(`CL` = 2.5, `tumor` = 0.5)) + theme(legend.position = "none"))
      } else {
        g <- ggplotly(ggplot(alignment, aes(UMAP_1, UMAP_2, size=type)) +
                        geom_point(data = alignment %>% filter(!(sampleID %in% selected_tumors) & type=='tumor'), color = 'gray', alpha=0.2) +
                        geom_point(data = alignment %>% filter(sampleID %in% selected_tumors), color="#EE4B33", alpha=0.5) +
                        geom_point(data = alignment %>% filter(type=="CL"), aes(color=CL_dist), alpha=0.7) +
                        scale_size_manual(values = c(`CL` = 2, `tumor` = 0.5)) +
                        theme(legend.position = "right"))
      }
      #if(input$feature=='tissue' & length(tumor_samples)==0) {g <- g + scale_color_manual(values=tissue_colors)}
      g
    })
    
    
    output$Data = renderDT({
      validate(need(input$cur_tissue, message = FALSE))
      validate(need(input$cur_subtype, message = FALSE))
      
      cur_used_subtypes <- input$cur_subtype
      if(cur_used_subtypes == "all") {
        cur_used_subtypes <- unique(filter(alignment, tissue == input$cur_tissue)$subtype)
      }
      
      selected_tumors <- filter(alignment, type=='tumor' & tissue == input$cur_tissue & subtype %in% cur_used_subtypes)$sampleID
      
      
      if(length(selected_tumors)>0) {
        CL_dist <- apply(tumor_CL_dist[selected_tumors,,drop=F], 2, FUN=median)
        CL_dist <- cbind.data.frame(CL_dist, colnames(tumor_CL_dist))
        colnames(CL_dist) <- c('distance_to_tumor', 'sampleID')
        CL_dist <- CL_dist %>% arrange(distance_to_tumor)
        CL_dist <- left_join(CL_dist, alignment[,c('sampleID', 'tissue','subtype')], by='sampleID')
        CL_dist
        #filter(CL_dist, CL %in% filter(alignment, tissue==input$cur_tissue)$sampleID)
      } else {
        data.frame(
          distance_to_tumor = integer(),
          sampleID = character(),
          tissue = character(),
          subtype = character()
        )
      }
    })
    
    ## CFE data
    output$CFE_plot <- renderPlot({
      validate(need(input$CFE, message = FALSE))
      ggplot(alignment, aes(UMAP_1, UMAP_2, size=type, color=type, 
                            fill = Depmap_ID %in% filter(CFEs, `CFE nodes` %in% input$CFE)$sampleID)) +
        geom_point(alpha=0.7, pch=21) +
        scale_size_manual(values = c(`CL` = 1.5, `tumor` = 1)) +
        scale_color_manual(values = c(`CL` = 'black', `tumor` = 'white')) + theme_classic() +
        theme(legend.position = "none")})
    
    ## add data
    

    

    output$new_data_plot <- renderPlotly({
      exp_input_file <- input$exp_file
      ann_input_file <- input$ann_file
      
      if (is.null(exp_input_file) | is.null(ann_input_file))
        return(NULL)
      else{
        if(input$runScript) {
          #progress <- shiny::Progress$new()
          # Make sure it closes when we exit this reactive, even if there's an error
          #on.exit(progress$close()
          exp_df <- readr::read_csv(exp_input_file$datapath) %>% as.data.frame() %>% column_to_rownames('X1')
          ann_df <- readr::read_csv(ann_input_file$datapath) %>% as.data.frame()
          source("./add_data_script.R", local = TRUE)
          output_df <- add_data_to_Celligner(exp_df, ann_df, input$mnn_k_val)
          output_df$new_data <- output_df$type
          output_df[which(!ouput_df$new_data %in% c('tumor', 'CL')),]$new_data <- 'other'
          g <- ggplotly(ggplot(output_df, aes(UMAP_1, UMAP_2, fill = tissue, size = new_data, color=new_data, text=paste(sampleID, tissue, subtype, type))) +
                          geom_point(alpha=0.7, pch=21) +
                          scale_color_manual(values=c(`other`='black', `CL` = 'gray90', `tumor`='white')) +
                          scale_size_manual(values = c(`other`=1.75, `CL`=1.25, `tumor`=0.75)) +
                          theme(legend.position = 'none'))
          return(g)
        } else {
          return(NULL)
        }
      }
    })

    # Downloadable csv of selected dataset ----
    # output$downloadData <- downloadHandler(
    #   filename = function() {
    #     paste(input$dataset, ".csv", sep = "")
    #   },
    #   content = function(file) {
    #     write.csv(datasetInput(), file, row.names = FALSE)
    #   }
    # )
   
  
    
  }

# Run the application 
shinyApp(ui = ui, server = server)

