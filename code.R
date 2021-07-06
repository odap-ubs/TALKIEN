require(shiny, quietly = T, warn.conflicts = F)
require(shinycssloaders, quietly = T, warn.conflicts = F)
require(shinyjs, quietly = T, warn.conflicts = F)
require(shinyWidgets, quietly = T, warn.conflicts = F)
require(dplyr, quietly = T, warn.conflicts = F)
require(igraph, quietly = T, warn.conflicts = F)
require(visNetwork, quietly = T, warn.conflicts = F)
require(ggplot2, quietly = T, warn.conflicts = F)
require(DT, quietly = T, warn.conflicts = F)

# status v4.6 --- rdata prepared to load

# preparing files to proceed with the App...
folder <- "./data/"
talkien_version <- "v4"

# loading data
load(paste0(folder, "data_talkien.RData"))

# user interface defining
ui <- fluidPage(
  
  # for showing hidding objects
  useShinyjs(),
  # title of the App
  titlePanel(title=div(img(src="talkien.jpg", height = 80, width = 320), " -crossTALK bIpartitE Network-"), windowTitle="Talkien -crossTALK bIpartitE Network-"),
  
  #layout and input options
  sidebarLayout(
    
    sidebarPanel(width = 3,
                 
                 talkien_version,
                 
                 # upload bed file option
                 h5("Please input two files with header and one column length, otherwise column #1 will be used as reference."),
                 h5("Both files must have the same annotation."),
                 
                 # toy example
                 checkboxInput(inputId = "clx", label = "Load example data", 
                               value = FALSE),
                 
                 conditionalPanel(condition = "input.clx == 0",
                                  fileInput("file1", "Input list #1:",
                                            multiple = FALSE,
                                            accept = "tsv"),
                                  fileInput("file2", "Input list #2:",
                                            multiple = FALSE,
                                            accept = "tsv")),
                 
                 # buttons for annotation types and for network interactions. Only if example is not selected
                 selectInput(inputId = "netint",
                             label = "Choose Network DataBase",
                             choices = list("STRING", "HIPPIE", "both"),
                             selected = "both"),
                 
                 br(),
                 
                 selectInput(inputId = "source",
                             label = "Choose info source",
                             choices = list("Human Protein Atlas" = "hpa", "Uniprot" = "up"),
                             selected = "hpa"),
                 
                 conditionalPanel(condition = "input.clx == 0",
                                  radioButtons("annot", "Choose type annotation:",
                                               choices = c(GeneSymbol = "gene_symbol", 
                                                           entrez = "entrezgene", 
                                                           proteinName = "protein_name",
                                                           uniprot = "uniprot_id"), 
                                               selected = "gene_symbol")),
                 
                 
                 radioButtons("nettype", "Choose type of network:",
                              choices = c(Whole = "who",
                                          Crosstalk = "crt"),
                              selected = "who"),
                 
                 # score filter
                 sliderInput(inputId = "score_slider", 
                             label = "Score threshold", 
                             min = 400, 
                             max = 900, 
                             value = 700, 
                             step = 50),
                 
                 # show downstream nodes
                 checkboxInput(inputId = "graylinks", label = "show downstream interactions",
                               value = FALSE),
                 
                 # network layouts
                 conditionalPanel(condition = "input.nettype == 'crt'",
                                  radioButtons("grlay2", "Graphical Layout:",
                                               choices = c(Circular = "layout_in_circle",
                                                           ForceDirected = "layout_with_fr",
                                                           Sphere = "layout_on_sphere",
                                                           Bipartite = "layout_as_bipartite"),
                                               selected = "layout_as_bipartite")),
                 
                 conditionalPanel(condition = "input.nettype == 'who'",
                                  radioButtons("grlay", "Graphical Layout:",
                                               choices = c(Circular = "layout_in_circle",
                                                           ForceDirected = "layout_with_fr", 
                                                           Sphere = "layout_on_sphere",
                                                           largeGraph = "layout_with_drl"),
                                               selected = "layout_with_fr"))
    ),
    
    # main panel for output. 4 main tabs and help tab
    mainPanel(
      
      tabsetPanel(id = "tabs",
                  tabPanel("Plot", icon = icon("project-diagram", lib = "font-awesome"), value = 1,
                           span(htmlOutput("warning1"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           conditionalPanel(condition = "input.grlay2 == 'layout_as_bipartite' & input.nettype == 'crt'",
                                            dropdownButton(radioButtons("bipnet", "bipartite by:",
                                                                        choices = c(list = "list", 
                                                                                    location = "location")),
                                                           materialSwitch("direct", "Up-Down",
                                                                          status = "primary",
                                                                          right = F),
                                                           inputId = "plotOpt", circle = T, size = "xs", tooltip = T, label = "type of bipartite net", icon = icon("bar-chart-o"), status = "primary")),
                           visNetworkOutput(outputId = "network1", width = "100%", height = "900px") %>% withSpinner(type = 8),
                           h5(downloadButton('downloadNetwork', 'Download network as .html'), align = "center"),
                           br(),
                           h4(strong("Node's Pathways"), align = "center"),
                           h5("click on any node to display all pathways on which is involved", align = "center"),
                           span(htmlOutput("downstreamID"), style = "color:gray; font-size:20px"),
                           dataTableOutput(outputId = "downstream"),
                           h5(downloadButton('downloadPaths', "Download node's pathways list"), align = "center")),
                  
                  tabPanel("Network Interactions", icon = icon("list-alt", lib = "font-awesome"), value = 2,
                           span(htmlOutput("warning2"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           h4(strong("Network Interaction List"), align = "center"),
                           dataTableOutput(outputId = "netlist") %>% withSpinner(type = 8),
                           h5(downloadButton('downloadNetTable', 'Download data as .tsv'), align = "center")),#),
                  
                  tabPanel("Network Parameters", icon = icon("list-alt", lib = "font-awesome"), value = 3,
                           span(htmlOutput("warning3"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           h4(strong("Network Descriptive Parameters"), align = "center"),
                           span(htmlOutput("info2"), style = "color:gray; font-size:16px"),
                           br(),
                           tableOutput(outputId = "table1") %>% withSpinner(type = 8),
                           br(),
                           h4(strong("Network Topological Parameters"), align = "center"),
                           dataTableOutput(outputId = "topolist") %>% withSpinner(type = 8),
                           h5(downloadButton('downloadTopoTable', 'Download data as .tsv'), align = "center")),#),
                  
                  tabPanel("Enrichment", icon = icon("list", lib = "font-awesome"), value = 4,
                           span(htmlOutput("warning4"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           conditionalPanel(condition = "input.tabs == '4'",
                                            dropdownButton(materialSwitch("globEnr", "Add downstream nodes",
                                                                          status = "primary",
                                                                          right = F),
                                                           radioButtons("showCat", "Number of showed pathways:",
                                                                        choices = c("10" = 10, 
                                                                                    "15" = 15,
                                                                                    "20" = 20,
                                                                                    "ALL" = "all"),
                                                                        selected = "15"),
                                                           inputId = "enrPlotOpt", circle = T, size = "xs", tooltip = T, label = "Add downstream nodes", icon = icon("bar-chart-o"), status = "primary"),
                                            h4(strong("Enrichment plot"), align = "center"),
                                            plotOutput("oraplot") %>% withSpinner(type = 8),
                                            h5(downloadButton('downloadPlot', 'Download plot as .png'), align = "center"),
                                            tags$hr(),
                                            h4(strong("Over Representation Analysis"), align = "center"),
                                            dataTableOutput(outputId = "oea") %>% withSpinner(type = 8),
                                            h5(downloadButton('downloadNetTable2', 'Download data as .tsv'), align = "center"))),
                  
                  tabPanel("Help", icon = icon("book", lib = "font-awesome"), value = 4, 
                           tags$iframe(src = "talkien-help4.html", width = "100%", height = "800"))
                  
      )
    )
  )
)

# server function
server <- function(input, output) { 
  
  # beggining with hide objetcs
  shinyjs::hide("downstream")
  shinyjs::hide("downstreamID")
  shinyjs::hide("downloadPaths")
  
  output$info2 <- renderText({ "" })
  outora <- reactiveValues(reacplot = NULL)
  
  net <- reactive({
    net <- list()
    data <- list()
    lists <- NULL
    gene_anot <- NULL
    anot <- NULL
    out_genes <- NULL
    if (input$clx == F) {   
      tryCatch({
        data$df1 <- read.delim(input$file1$datapath,
                               header = T,
                               stringsAsFactors = F)
        data$df2 <- read.delim(input$file2$datapath,
                               header = T,
                               stringsAsFactors = F)
        lists$list1 <- gsub("\\..*", "", input$file1$name)
        lists$list2 <- gsub("\\..*", "", input$file2$name)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        validate(
          need(input$file1 != "", ""),
          need(input$file2 != "", "")
        )
        stop(safeError(e))
      })
      gene_annot <- input$annot
    } else {
      data$df1 <- clx_norm
      data$df2 <- clx_tum
      lists$list1 <- "Normal"
      lists$list2 <- "CRC"
      gene_annot <- "gene_symbol"
    }
    
    if (nrow(data$df1)==0 | nrow(data$df2)==0) {
      
      output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      net$data <- data
      
    } else {
      
      # beggining to analyze data. first, removing duplicate ids
      data$df1 <- gsub(" .*", "", data$df1[, 1])
      data$df2 <- gsub(" .*", "", data$df2[, 1])
      data$df1 <- data$df1[!duplicated(data$df1)]
      data$df2 <- data$df2[!duplicated(data$df2)]
      
      # getting common genes in both lists and removing them
      out_genes$comGenes <- data$df1[data$df1%in%data$df2]
      data$df1 <- data$df1[!data$df1%in%out_genes$comGenes]
      data$df2 <- data$df2[!data$df2%in%out_genes$comGenes]
      
      # first check, just a freak joke... 
      if (length(grep("quenya", data$df1))>0 | length(grep("quenya", data$df2))>0) {

        output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        net <- list(data)
        
      } else {
        
        # getting non common genes between lists... the ones we want
        out_genes$uncomGenes <- c(data$df1, data$df2)
        
        # second check, if everything is common, no analysis is possible
        if (length(out_genes$uncomGenes)==0) {
          
          output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          net <- list(data, out_genes)
          
        } else {
          
          # for the sake of clarity, base annotation is gene symbol
          data$df1 <- descr[match(data$df1, table = descr[, gene_annot]), "gene_symbol"]
          data$df2 <- descr[match(data$df2, table = descr[, gene_annot]), "gene_symbol"]
          data$df1 <- data$df1[!is.na(data$df1)]
          data$df2 <- data$df2[!is.na(data$df2)]
          
          # third check, format must be the same in both lists
          if(length(data$df1)==0 || length(data$df2)==0) {
            
            output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            net <- list(data, out_genes)
            
          } else {
            
            # getting all genes needed to produce the net: genes in lists with functions as Membrane receptor prot or Secreted prot
            ppi <- c(data$df1, data$df2)
            ppi <- ppi[!duplicated(ppi)]
            ppi <- as.data.frame(ppi)
            names(ppi) <- "gene_symbol"
            ppi_mapped <- merge(ppi, descr, all = F)
            if (input$source == "hpa") {
              anot$annot_genes <- c(as.character(ppi_mapped$gene_symbol[ppi_mapped$gene_symbol%in%membRecept$Gene]), as.character(ppi_mapped$gene_symbol[ppi_mapped$gene_symbol%in%secrProt$Gene]))
            } else {
              anot$annot_genes <- c(as.character(ppi_mapped$gene_symbol[ppi_mapped$gene_symbol%in%uni_memb$Gene]), as.character(ppi_mapped$gene_symbol[ppi_mapped$gene_symbol%in%uni_secr$Gene]))
            }
            
            # fourth check, needed receptor o secretor proteins... also in our list
            if(length(anot$annot_genes)==0) {
              
              output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              net <- list(data, out_genes, anot)
              
            } else {
              
              # getting only interactions between receptors and or secretors, the other ones to downstream purposes
              ppi_mapped$anot <- rep("downstream", nrow(ppi_mapped))
              ppi_mapped$anot[ppi_mapped$gene_symbol%in%anot$annot_genes] <- "main"
              
              # transforming to igraph subnetwork... and the same for downstream network 
              ppi_mapped <- ppi_mapped[ppi_mapped$string_id%in%names(V(ppiGen_graph)), ]
              net <- list(data, out_genes, anot, lists, ppi_mapped)
              
            }
          }
        }
      }    
    } 
    
  })
  
  net_clean <- reactive({
    
    net_clean <- list()
    net <- net()
    data <- net[[1]]


    if (length(data$df1)==0 | length(data$df2)==0) {
    
    } else {
      
      # first check, just a freak joke... 
      if (length(grep("quenya", data$df1))>0 | length(grep("quenya", data$df2))>0) {
      
      } else {
        
        out_genes <- net[[2]]
            
        # second check, if everything is common, no analysis is possible
        if (length(out_genes$uncomGenes)==0) {
        
        } else {
          
          # third check, format must be the same in both lists
          if(is.na(data$df1) || is.na(data$df2)) {
          
          } else {
            
            anot <- net[[3]]
                
            # fourth check, needed receptor o secretor proteins...
            if(length(anot$annot_genes)==0) {
            
            } else {
                  
              lists <- net[[4]]
              ppi_mapped <- net[[5]]
              
              if (input$netint == "STRING") {
                disnet <- names(V(ppiGen_graph))[names(V(ppiGen_graph))%in%ppi_mapped$string_id]
                ppi_subgraph <- induced_subgraph(ppiGen_graph, disnet)
              } else if (input$netint == "HIPPIE") {
                disnet <- names(V(ppiGen_graph2))[names(V(ppiGen_graph2))%in%ppi_mapped$string_id]
                ppi_subgraph <- induced_subgraph(ppiGen_graph2, disnet)
              } else {
                disnet <- names(V(ppiGen_graph))[names(V(ppiGen_graph))%in%ppi_mapped$string_id]
                ppi_subgraph <- induced_subgraph(ppiGen_graph, disnet)
                disnet2 <- names(V(ppiGen_graph2))[names(V(ppiGen_graph2))%in%ppi_mapped$string_id]
                ppi_subgraph2 <- induced_subgraph(ppiGen_graph2, disnet2)
                ppi_subgraph <- ppi_subgraph + ppi_subgraph2
              }
              
              # fifth check, needed receptor o secretor proteins... also in our list
              if(length(E(ppi_subgraph))==0) {
                
                output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                net_clean <- list(ppi_subgraph)
                
              } else {
                
                # setting gene names as labels
                V(ppi_subgraph)$labels <- as.character(ppi_mapped[match(names(V(ppi_subgraph)), table=ppi_mapped$string_id), "gene_symbol"])
                V(ppi_subgraph)$int <- as.character(ppi_mapped[match(names(V(ppi_subgraph)), table=ppi_mapped$string_id), "anot"])
                
                # creating a data frame with empty node attributes
                net_attrib <- data.frame(gene = V(ppi_subgraph)$labels, fun = rep("downstream", length(V(ppi_subgraph))), list = rep("", length(V(ppi_subgraph))), size = rep("", length(V(ppi_subgraph))), color = rep("", length(V(ppi_subgraph))), shape = rep("", length(V(ppi_subgraph))), stringsAsFactors = F)
                
                # ... populating it
                if (input$source == "hpa") {
                  net_attrib$fun[net_attrib$gene%in%membRecept$Gene] <- "Receptor"
                  net_attrib$fun[net_attrib$gene%in%secrProt$Gene] <- "Secreted"
                } else {
                  net_attrib$fun[net_attrib$gene%in%uni_memb$Gene] <- "Receptor"
                  net_attrib$fun[net_attrib$gene%in%uni_secr$Gene] <- "Secreted"
                }
                
                net_attrib$list[net_attrib$gene%in%data$df1] <- lists$list1
                net_attrib$list[net_attrib$gene%in%data$df2] <- lists$list2
                
                # transforming to visnetwork object
                visData <- toVisNetworkData(ppi_subgraph)
                net$nodes <- visData$nodes
                net$edges <- visData$edges
                
                if (input$netint == "both") {
                  net$edges$combined_score <- apply(net$edges[,3:4], 1, function(x) max(x, na.rm = T))
                  net$edges <- net$edges[, c(1,2,6)]
                }
                
                # reformating visnetworkdata object for plots
                net$nodes <- net$nodes[, -4]
                names(net$nodes)[2] <- "label"
                net$nodes <- merge(net$nodes, net_attrib, all = F, by.x = "label", by.y = "gene")
                
                # prune edges between secrprot--secrprot.
                net$edges$int <- rep("down", nrow(net$edges))
                net$edges$int[net$edges$from%in%net$nodes$id[net$nodes$fun=="Secreted"] & net$edges$to%in%net$nodes$id[net$nodes$fun=="Receptor"]] <- "main"
                net$edges$int[net$edges$to%in%net$nodes$id[net$nodes$fun=="Secreted"] & net$edges$from%in%net$nodes$id[net$nodes$fun=="Receptor"]] <- "main"
                
                # adding attributes
                net$edges$color <- rep("lightgrey", nrow(net$edges))
                
                # aditional filter by score depending on user filter
                pruned_network <- net$edges[net$edges$combined_score>=input$score_slider, ]
                
                # sixth check, selected interactions must be more than 0
                if(nrow(pruned_network)==0) {
                  
                  output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                  output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                  output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                  output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                  net_clean <- list(ppi_subgraph, pruned_network)
                  
                } else {
                  
                  vert <- net$nodes[net$nodes$id%in%unique(c(pruned_network$from, pruned_network$to)), ]
                  
                  # more information to the table
                  pruned_network$width <- pruned_network$combined_score/min(pruned_network$combined_score)
                  pruned_network[,7:12] <- rep("", NROW(pruned_network))
                  names(pruned_network)[7:12] <- c("interaction_typeA", "interaction_typeB", "interaction_locationA", "interaction_locationB", "ProteinA", "ProteinB")
                  
                  # computing interaction types
                  gene1 <- vert[match(pruned_network$from, vert$id), c(1,3,4,5)]
                  gene2 <- vert[match(pruned_network$to, vert$id), c(1,3,4,5)]
                  pruned_network$interaction_typeA <- gene1$list
                  pruned_network$interaction_typeB <- gene2$list
                  pruned_network$interaction_locationA <- gene1$fun
                  pruned_network$interaction_locationB <- gene2$fun
                  pruned_network$ProteinA <- gene1$label
                  pruned_network$ProteinB <- gene2$label
                  
                  # removing secreted - secreted, downstream - secreted, secreted - downstream (can be modified easily)
                  int <- paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")
                  pruned_network <- pruned_network[!int %in% c("Secreted_Secreted", "Secreted_downstream", "downstream_Secreted"), ]
                  # downstream proteins must be from the same list of their receptor activators, so prune edges between receptor-list1 and downstream-list2 and viceversa
                  int <- paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")
                  pruned_network[int%in%"downstream_Receptor", ] <- pruned_network[int%in%"downstream_Receptor", c(2,1,3:6,8,7,10,9,12,11)]
                  pruned_network[int%in%"Receptor_Secreted", ] <- pruned_network[int%in%"Receptor_Secreted", c(2,1,3:6,8,7,10,9,12,11)]
                  int <- paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")
                  lis <- paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = "_")
                  rem_edges <- int%in%"Receptor_downstream" & lis%in%c(paste(lists$list1, lists$list2, sep = "_"), paste(lists$list2, lists$list1, sep = "_"))
                  pruned_network <- pruned_network[!rem_edges, ]
                  # now, the same between downstream genes, so prune edges between downstream1-list2 -- downstream2-list1 and viceversa
                  int <- paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")
                  lis <- paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = "_")
                  rem_edges <- int%in%c("downstream_downstream") & lis%in%c(paste(lists$list1, lists$list2, sep = "_"), paste(lists$list2, lists$list1, sep = "_"))
                  pruned_network <- pruned_network[!rem_edges, ]
                  
                  # removing components with downstream genes only...
                  vert <- vert[vert$id%in%unique(c(pruned_network$from, pruned_network$to)), ]
                  downGenes <- vert$id[vert$int == "downstream"]
                  down_comps <- data.frame(membership = components(graph_from_data_frame(pruned_network[, 1:2]))$membership)
                  down_comps <- merge(vert[, c(2,4)], down_comps, by.y = "row.names", by.x = "id")
                  down_comps$rem_vert <- rep(T, nrow(down_comps))
                  # removing nodes where there are no receptors in the whole component
                  for (i in 1:length(table(down_comps$membership))) {
                    tmp <- down_comps$fun[down_comps$membership == i]
                    if("Receptor"%in%tmp == F) {
                      down_comps$rem_vert[down_comps$membership == i] <- F
                    } else {
                      if (!"Secreted"%in%names(table(tmp)) & unname(table(tmp)[names(table(tmp))%in%"Receptor"])==1) down_comps$rem_vert[down_comps$membership == i] <- F
                    }
                  }
                  keep_vert <- down_comps$id[down_comps$rem_vert == T]
                  pruned_network <- pruned_network[pruned_network$from%in%keep_vert | pruned_network$to%in%keep_vert, ]
                  int <- paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")
                  pruned_network$color[int%in%c("downstream_downstream", "Receptor_downstream")] <- 'rgba(0,0,0,0)'
                  vert <- vert[vert$id %in% pruned_network$from | vert$id %in% pruned_network$to, ]
                  # double links on downstream downstream and receptor receptor
                  bidir_downs <- pruned_network[paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")%in%c("downstream_downstream", "Receptor_Receptor"), c(2,1,3:6,8,7,10,9,12,11)]
                  names(bidir_downs) <- names(pruned_network)
                  
                  downstream_network <- rbind.data.frame(pruned_network, bidir_downs)
                  downcross_network <- pruned_network
                  
                  crtGenes <- c()
                  if (input$nettype == "crt") {
                    crtGenes <- vert$label[vert$int=="main"]
                    pruned_network <- pruned_network[!paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = "-")%in%paste(lists$list1, lists$list1, sep = "-"), ]
                    pruned_network <- pruned_network[!paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = "-")%in%paste(lists$list2, lists$list2, sep = "-"), ]
                    pruned_network <- pruned_network[paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "-") =="Receptor-Secreted" | paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "-") =="Secreted-Receptor", ]
                    vert_crt <- vert[vert$id%in%unique(c(pruned_network$from, pruned_network$to)), ]
                    crtGenes <- crtGenes[!crtGenes%in%unique(c(pruned_network$ProteinA, pruned_network$ProteinB))]
                    
                    # seventh check, if crosstalk selected, still interactions found...
                    if(nrow(pruned_network)==0) {
                      
                      output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      net_clean <- list(ppi_subgraph, pruned_network)
                      
                    } else {
                      
                      tmp_cross <- downstream_network[paste(downstream_network$interaction_locationA, downstream_network$interaction_locationB, sep = "_")%in%"downstream_downstream", ]
                      tmp_cross2 <- downstream_network[downstream_network$ProteinA%in%vert_crt$label, ]
                      tmp_cross2 <- tmp_cross2[tmp_cross2$interaction_locationB == "downstream", ]
                      down_comps <- data.frame(membership = components(graph_from_data_frame(rbind.data.frame(pruned_network, tmp_cross, tmp_cross2)[, 1:2]))$membership)
                      down_comps <- merge(vert[, c(2,4)], down_comps, by.y = "row.names", by.x = "id")
                      down_comps$rem_vert <- rep(F, nrow(down_comps))
                      # removing nodes where there are no receptors in the whole component
                      for (i in 1:length(table(down_comps$membership))) {
                        tmp <- down_comps$fun[down_comps$membership == i]
                        if("Receptor"%in%tmp == T  & "Secreted"%in%tmp == T) down_comps$rem_vert[down_comps$membership == i] <- T
                      }
                      keep_vert <- down_comps$id[down_comps$rem_vert == T]
                      downcross_network <- rbind.data.frame(pruned_network, tmp_cross, tmp_cross2)
                      downcross_network <- downcross_network[downcross_network$from%in%keep_vert | downcross_network$to%in%keep_vert, ]
                      int <- paste(downcross_network$interaction_locationA, downcross_network$interaction_locationB, sep = "_")
                      #downcross_network$color[int%in%c("downstream_downstream", "Receptor_downstream")] <- "lightgrey" #'rgba(0,0,0,0)'
                      vert <- vert[vert$id %in% downcross_network$from | vert$id %in% downcross_network$to, ]
                      # double links on downstream downstream and receptor receptor
                      bidir_downs <- downcross_network[paste(downcross_network$interaction_locationA, downcross_network$interaction_locationB, sep = "_")%in%c("downstream_downstream", "Receptor_Receptor"), c(2,1,3:6,8,7,10,9,12,11)]
                      names(bidir_downs) <- names(downcross_network)
                      downstream_network <- rbind.data.frame(downcross_network, bidir_downs)
                      
                    } # end of the seventh check
                    
                  } else {
                    pruned_network <- pruned_network[paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "-") =="Receptor-Receptor" | paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "-") =="Secreted-Receptor", ]
                  }
                  
                  if(nrow(pruned_network)==0) {
                    
                    output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                    output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                    output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                    output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                    net_clean <- list(ppi_subgraph, pruned_network)
                    
                  } else {
                    
                    output$warning1 <- renderText({ "" })
                    output$warning2 <- renderText({ "" })
                    output$warning3 <- renderText({ "" })
                    output$warning4 <- renderText({ "" })
                    
                    # formatting network plot
                    vert$shape <- "square"
                    vert$shape[vert$fun=="Receptor"] <- "triangleDown"
                    vert$shape[vert$fun=="Secreted"] <- "dot"
                    vert$color[vert$list==lists$list1] <- "coral"
                    vert$color[vert$list==lists$list2] <- "steelblue"
                    vert$color[vert$fun == "downstream"] <- 'rgba(0,0,0,0)'
                    vert$color[!vert$id%in%unique(c(pruned_network$from, pruned_network$to))] <- 'rgba(0,0,0,0)'
                    
                    # computing size of nodes (here from downstream, not from pruned, to be discussed)
                    downcross_graph <- graph_from_data_frame(downcross_network, directed = F)
                    size_vert <- degree(downcross_graph, v = V(downcross_graph), mode = "all")
                    size_vert <- data.frame(id = names(size_vert), size = unname(size_vert))
                    vert$size <- size_vert$size[match(vert$id, table = size_vert$id)]
                    
                    # reformating size for better visualization, and other vars
                    vert$size <- as.numeric(vert$size)
                    vert$size <- round(5*log2(vert$size+1))
                    vert$title <- vert$label
                    
                    # getting final network tables
                    nodes <- descr[match(vert$label, table = descr$gene_symbol), c(4,3,1)]
                    nodes[,4:12] <- vert
                    nodes <- nodes[, c(4:12,1:3)]
                    nodes <- nodes[order(nodes$color, nodes$shape), ]
                    nodes <- nodes[!duplicated(nodes$label), ]
                    
                    # getting final tables
                    pruned_network <- pruned_network[, c(11,12,9,10,7,8,3)]
                    downcross_network <- downcross_network[, c(11,12,9,10,7,8,3)]
                    
                    graph_stat <- graph_from_data_frame(pruned_network[, 1:2])
                    graph_stat2 <- graph_from_data_frame(downcross_network[, 1:2])
                    
                    # enrichment analysis
                    enr <- NULL
                    
                    if (input$graylinks) {
                      netstat <- toposummary(graph_stat2)
                    } else {
                      netstat <- toposummary(graph_stat) 
                    }
                    netstat <- merge(netstat, nodes, by = "label", all = F)
                    netstat <- netstat[order(netstat$degree, decreasing = T), c(1,2,6,7,8,4,9,14,15,22,20,21,12)]
                    netstat[, 4:7] <- round(netstat[, 4:7], digits = 3)
                    netstat[, 3] <- formatC(netstat[, 3], digits = 2, format = "E")
                    names(netstat) <- c("Gene_Symbol", "Degree", "Closeness", "Betweenness", "Eigen_Centrality", "Clustering_coef", "Eccentricity", "Location", "List", "Entrez_ID", "protein_name", "uniprot_ID", "String_ID")
                    net_clean <- list(ppi_subgraph, pruned_network, nodes, downcross_network, graph_stat, graph_stat2, lists, downstream_network, enr, netstat)
                  
                  } # end of the seventh check
                } # end of the sixth check
              }   # end of the fifth check
            }     # end of the fourth check
          }       # end of the third check
        }         # end of the second check
      }           # end of the first check
    }             # end of the sanity check
  })
  
  ora <- reactive({
    
    net_clean <- net_clean()
    ora <- list()
    
    if (length(net_clean)==10) {
      ora <- list()
      nodes <- net_clean[[3]]
      
      if (input$tabs == '4') {
        print("performing enrichment...")
        if (input$globEnr ==F) nodes <- nodes[nodes$fun != "downstream", ]
        enr <- nodes$entrez
        enr <- ReactomePA::enrichPathway(gene = enr, pvalueCutoff = 0.01, readable = T)
        ora <- list(enr)
      }
    }
  })
  
  output$table1 <- renderTable({
    
    net_clean <- net_clean()
    
    if (length(net_clean)==10) {
      
                  nodes <- net_clean[[3]]
                  pruned_network <- net_clean[[2]]
                  downcross_network <- net_clean[[4]]
                  graph_stat <- net_clean[[5]]
                  graph_stat2 <- net_clean[[6]]
                  
                  if (components(graph_stat)$no>1) {
                    output$info2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Disconected graph. Showing total nodes and edges, rest of the parameters computed only for giant component")) })
                  } else {
                    output$info2 <- renderText({ "" })
                  }
                  
                  netsum <- netsummary(graph_stat)
                  netsum2 <- netsummary(graph_stat2)
                  netsum[1:3] <- sprintf("%1.0f", netsum[1:3])
                  netsum2[1:3] <- sprintf("%1.0f", netsum2[1:3])
                  
                  if (input$graylinks) {
                    netsum <- rbind.data.frame(netsum, netsum2)
                    row.names(netsum) <- c("Secreted - Receptor", "Downstream")
                  }
                  netsum
    }
    
  })
  
  output$netlist <- renderDataTable({
    
    net_clean <- net_clean()
    
    if (length(net_clean)==10) {
      
                pruned_network <- net_clean[[2]]
                
                if(nrow(pruned_network)==0) {
                  
                } else {
                  
                  if (input$graylinks) {
                    edges <- net_clean[[4]]
                  } else {
                    edges <- net_clean[[2]]
                  }
                  datatable(edges, options = list(pageLength = 15))
                }
              }
  })
  
  output$topolist <- renderDataTable({
    
    net_clean <- net_clean()
    
    if (length(net_clean)==10) {
      
      netstat <- net_clean[[10]]
      datatable(netstat)
    }
    
  })
  
  output$network1 <- renderVisNetwork({
    
    net_clean <- net_clean()
    
    if (length(net_clean)==10) {
                  
      nodes <- net_clean[[3]]
      lists <- net_clean[[7]]
      downstream_network <- net_clean[[8]]
      pruned_network <- net_clean[[2]]
      # adding legend
      lnodes <- data.frame(label = c(paste(lists$list1, "Receptor", sep = "\n"), 
                                     paste(lists$list1, "Secreted", sep = "\n"), 
                                     paste(lists$list2, "Receptor", sep = "\n"), 
                                     paste(lists$list2, "Secreted", sep = "\n")), 
                           shape = c("triangleDown", "dot", "triangleDown", "dot"), 
                           color = c("coral", "coral", "steelblue", "steelblue"), stringsAsFactors = F)
      
      if (input$graylinks == T) {
        downstream_network$color <- "lightgrey"
        nodes$color[nodes$color == 'rgba(0,0,0,0)'] <- "lightgrey"
        lnodes <- rbind(lnodes, c("Downstream", "square", "lightgrey"), c("Downstream\nselected", "square", "gray25"))
      } else {
        nodes[nodes$color == 'rgba(0,0,0,0)', c(1,9)] <- ""
      }
      
      # if not bipartite... default layout
      if (input$nettype != 'crt') {
        lay <- input$grlay
      } else {
        lay <- input$grlay2
      }

      if (lay == "layout_in_circle" & input$graylinks == T | lay != "layout_as_bipartite" & lay != "layout_in_circle") {
        subnetwork1 <- visNetwork(nodes, downstream_network, main = paste(lists$list1, lists$list2, sep = " X ")) %>% #, width = "100%", height = "100%") %>%
          visIgraphLayout(layout = lay, randomSeed = 1234) %>%
          visNodes() %>%
          visPhysics(stabilization = F)
        subnetwork <- subnetwork1 %>%
          visLegend(addNodes = lnodes, useGroups = F, width = 0.2, main = "Legend", position = "right") %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = F, autoResize = T) %>%
          visInteraction(navigationButtons = TRUE) %>%
          visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}")
        
        # if bipartite, a bit different
      } else if (lay == "layout_in_circle" & input$graylinks == F) {
        subnetwork1 <- visNetwork(nodes[nodes$color != 'rgba(0,0,0,0)', ], downstream_network[downstream_network$color != 'rgba(0,0,0,0)', ], main = paste(lists$list1, lists$list2, sep = " X ")) %>% #, width = "100%", height = "100%") %>%
          visIgraphLayout(layout = lay, randomSeed = 1234) %>%
          visNodes() %>%
          visPhysics(stabilization = F)
        subnetwork <- subnetwork1 %>%
          visLegend(addNodes = lnodes, useGroups = F, width = 0.2, main = "Legend", position = "right") %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = F, autoResize = T) %>%
          visInteraction(navigationButtons = TRUE) %>%
          visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}")
      } else {
        # conditionalpanel reformating
        if(input$direct == T) {dir = "UD"} else {dir = "LR"}
        # reformating stuff to show the bipartite net
        if (input$bipnet == "list") {
          nodes$level <- as.numeric(as.factor(nodes$list))
          nodes$level[nodes$int=="downstream"] <- 3
          nodes$x <- as.numeric(as.factor(nodes$color))
          nodes <- nodes[order(nodes$level, nodes$x), ]
          downstream_left <- downstream_network$to[downstream_network$from%in%nodes$id[nodes$list==lists$list2 & nodes$fun == "Receptor"] & downstream_network$interaction_locationB=="downstream"]
          downstream_left <- unique(c(downstream_left, downstream_network$to[downstream_network$from%in%downstream_left]))
          if (length(downstream_left)>0) {
            nodes$level <- nodes$level+1
            nodes$level[nodes$id%in%downstream_left] <- 1
          }
        } else {
          nodes$level <- as.numeric(factor(nodes$fun, levels = c("Secreted", "Receptor", "downstream")))
          nodes$x <- as.numeric(factor(nodes$shape))
          nodes <- nodes[order(nodes$level, nodes$x), ]
        }
        
        subnetwork1 <- visNetwork(nodes, downstream_network, main = paste(lists$list1, lists$list2, sep = " X ")) %>% #, width = "100%", height = "100%") %>%
          visHierarchicalLayout(direction = dir, levelSeparation = 500, treeSpacing = 0.001, nodeSpacing = 25, blockShifting = T) %>%
          visEdges(scaling = list(min=0.2, max=0.5))
        subnetwork <- subnetwork1 %>%
          visOptions(highlightNearest = T, nodesIdSelection = F, autoResize = T) %>%
          visLegend(useGroups = F, addNodes = lnodes, position = "right", zoom = T, main = "Legend", width = 0.1) %>%
          visPhysics(stabilization = T) %>%
          visInteraction(navigationButtons = TRUE) %>%
          visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}")
        
        subnetwork
        
      }
    }
  })
  
  output$oraplot <- renderPlot({
    
    ora <- ora()
    
    if (length(ora)>0) {
      
      enr <- ora[[1]]
      if (input$showCat == "all") {
        print("im here")
        print(nrow(enr))
        numCat = nrow(enr)
      } else {
        numCat = as.numeric(input$showCat)
      }
      # enrichment analysis
      if(nrow(enr)==0) {
        output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "No enrichment terms found!")) })
      } else if (nrow(enr) < numCat) {
        output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Only found", nrow(enr), "enrichment terms")) })
        reacplot <- function(){
          enrichplot::dotplot(enr, showCategory = nrow(enr))
        }
        outora$reacplot <- reacplot()
        reacplot()
      } else {
        reacplot <- function(){
          enrichplot::dotplot(enr, showCategory = numCat)
        }
        outora$reacplot <- reacplot()
        reacplot()
      }
    }
  })
  
  output$oea <- renderDataTable({
    
    ora <- ora()
    
    if (length(ora)>0) {
      
      enr <- ora[[1]]
      oratable <- as.data.frame(enr@result)
      oratable[, 5:7] <- apply(oratable[,5:7], 2, function(x) formatC(x, format = "e", digits = 2))
      oratable$geneID <- gsub("/", " / ", oratable$geneID)
      IDs <- oratable[,1]
      oratable <- cbind(oratable, createLink(IDs))
      colnames(oratable)[ncol(oratable)] <- "link"
      outora$oratable <- oratable
      datatable(oratable, rownames = F,
                options = list(pageLength = 10,
                               autoWidth = T,
                               columnDefs = list(list(targets = 6, width = '150px')), scrollX = T),
                escape = F)
    }
  })
  
  output$downloadPlot = downloadHandler(
    filename = function() {'oraplot.png'},
    content = function(file) {
      reacplot <- outora$reacplot
      png(file, width = 1200, height = 700)
      print(reacplot +
              ggtitle("\nEnrichment Plot") +
              theme(plot.title = element_text(hjust = 0.5)))
      dev.off()
    })
  
  output$downloadNetTable2 <- downloadHandler(
    filename = function() {'oratable.txt'}, content = function(file) {
      oratable <- outora$oratable
      oratable[,10] <- as.character(oratable[,10])
      oratable[,10] <- substr(oratable[,10], 10, nchar(oratable[,10])-55)
      write.table(oratable, file, row.names = FALSE, quote = F, sep = "\t")
    })
  
  output$downloadNetTable <- downloadHandler(
    filename = function() {'netlist.txt'}, content = function(file) {
      write.table(net_clean()[[2]], file, row.names = FALSE, quote = F, sep = "\t")})
  
  output$downloadTopoTable <- downloadHandler(
    filename = function() {'nodelist.txt'}, content = function(file) {
      write.table(net_clean()[[10]], file, row.names = FALSE, quote = F, sep = "\t")})
  
  observe({
    
    net_clean <- net_clean()
    
    if (length(net_clean)==10) {
                  
      nodes <- net_clean[[3]]
      downstream_network <- net_clean[[8]]
      nodes_sel_entrez <- NULL
      nodes_in <- NULL
      if(length(input$click) == 1) {
        shinyjs::show("downstream")
        shinyjs::show("downstreamID")
        shinyjs::show("downloadPaths")
        visNetworkProxy("network1") %>%
          visGetSelectedNodes(input = "connected_nodes")
        filt <- grepl(paste("^", input$connected_nodes, "$", sep =""), downstream_network$from) & downstream_network$interaction_locationB!="downstream"
        downstream_network_tmp <- downstream_network[!filt, ]
        down_graph <- graph_from_data_frame(downstream_network_tmp, directed = T)
        if (!is.null(input$connected_nodes) & input$graylinks == T) {
          nodes$color[nodes$color == 'gray25' | nodes$color == "rgba(0,0,0,0)"] <- "lightgrey"
          if (input$connected_nodes%in%V(down_graph)$name) {
            down_sel <- unname(bfs(down_graph, root = input$connected_nodes, neimode = "out", unreachable = F)$order)
            nodes_sel <- V(down_graph)[down_sel[!is.na(down_sel)][-1]]$name
            nodes_sel_entrez <- descr$entrezgene[match(nodes_sel, descr$string_id)]
            nodes_sel_entrez <- nodes_sel_entrez[!is.na(nodes_sel_entrez)]
            sel_in <- length(nodes_sel_entrez)
            nodes <- changeColorOfNodes(nodes, nodes_sel)
            print(nodes$color)
            visNetworkProxy("network1") %>%
              visUpdateNodes(nodes) 
          }
        } else if (input$graylinks == F) {
          visNetworkProxy("network1") %>%
            visUpdateNodes(nodes[nodes$fun != "downstream",]) %>%
            visUpdateEdges(downstream_network[downstream_network$interaction_locationB != "downstream", ])
        }
        ent_sel <- descr$entrezgene[grep(paste("^", input$connected_nodes, "$", sep = ""), descr$string_id)]
        sym_sel <- descr$gene_symbol[grep(paste("^", input$connected_nodes, "$", sep = ""), descr$string_id)]
        hsa <- xx[grep(paste("^", as.character(ent_sel), "$", sep = ""), names(xx))]
        hsa <- names(xy[names(xy)%in%unname(unlist(hsa))])
        ids <- gsub("Homo sapiens: ", "", unlist(unname(xy[names(xy)%in%unname(unlist(hsa))])))
        reaclink <- createLink(hsa)
        if (!is.null(input$connected_nodes)) { 
          if(input$graylinks == T) {
            if (input$connected_nodes%in%V(down_graph)$name) {
              nodes_in <- do.call(rbind, lapply(hsa, FUN = function(x) nodes_sel_entrez%in%as.integer(reactome.db::reactomePATHID2EXTID[[x]])))
              if (length(nodes_in)!=0) {
                path_len <- do.call(rbind, lapply(hsa, function(x)length(reactome.db::reactomePATHID2EXTID[[x]])))
                univ_len <- nrow(descr)-path_len
                nodes_in <- apply(nodes_in, 1, function(x) nodes_sel_entrez[x])
                path_in <- do.call(rbind, lapply(nodes_in, length))
                nodes_in <- lapply(nodes_in, function(x) descr$gene_symbol[match(x, descr$entrezgene)])
                nodes_in <- unlist(lapply(nodes_in, function(x) paste(x, collapse = " / ")))
                contTable <- data.frame(univ_len, path_len, sel_in, path_in)
                pvals <- apply(contTable, 1, function(x) unlist(fisher.test(matrix(x, nrow = 2))[1]))
                output$downstreamID <- renderUI({ HTML(paste('<br/>', '<br/>', "<center>", "<b>", sym_sel, "</b>", "</center>"))})
                output$downstream <- renderDataTable({
                  downpaths <- cbind.data.frame(hsa, ids, path_len, pvals, nodes_in, reaclink)
                  downpaths <- downpaths[order(downpaths$pvals, decreasing = F), ]
                  downpaths$pvals <- formatC(downpaths$pvals, digits = 2, format = "E")
                  colnames(downpaths) <- c('ReactomeID', "Pathway_name", "pathway_length", "p-value", "nodes_in", "link")
                  downpaths <<- downpaths
                  return(downpaths)
                  
                }, escape = FALSE)
              }
            }
          }
          
          if (input$graylinks == F | length(nodes_in)==0 | !input$connected_nodes%in%V(down_graph)$name) {
            
            output$downstreamID <- renderUI({ HTML(paste('<br/>', '<br/>', "<center>", "<b>", sym_sel, "</b>", "</center>"))})
            output$downstream <- renderDataTable({
              
              downpaths <- cbind(hsa, ids, reaclink)
              colnames(downpaths) <- c('ReactomeID', "Pathway_name", "link")
              downpaths <<- downpaths
              return(downpaths)
              
            }, escape = FALSE)
          }
          
        } 
        
        # downloadable dynamic object
        output$downloadPaths <- downloadHandler(
          filename = function() {paste(sym_sel, 'pathways.tsv', sep="_")}, content = function(file) {
            downpaths[, ncol(downpaths)] <- substr(downpaths[, ncol(downpaths)], 10, nchar(downpaths[,ncol(downpaths)])-55)
            write.table(downpaths, file, row.names = FALSE, quote = F, sep = "\t")
          })
        
      } else {
        shinyjs::hide("downstream")
        shinyjs::hide("downstreamID")
        shinyjs::hide("downloadPaths")
      }
    }
  })
  
  output$downloadNetwork <- downloadHandler(
    filename = function() {
      paste('network-', Sys.Date(), '.html', sep='')
    },
    content = function(con) {
      visNetwork(nodes = net_clean()[[3]], edges = net_clean()[[8]], width = "100%") %>%
        visOptions(highlightNearest = TRUE) %>% visExport() %>%
        visPhysics(enabled = FALSE) %>% visEdges(smooth = FALSE) %>% visSave(con)
    }
  )
}

shinyApp(ui, server)
