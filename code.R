require(shiny, quietly = T, warn.conflicts = F)
require(shinycssloaders, quietly = T, warn.conflicts = F)
require(shinyjs, quietly = T, warn.conflicts = F)
require(shinyWidgets, quietly = T, warn.conflicts = F)
require(dplyr, quietly = T, warn.conflicts = F)
require(igraph, quietly = T, warn.conflicts = F)
require(visNetwork, quietly = T, warn.conflicts = F)
require(ggplot2, quietly = T, warn.conflicts = F)
require(DT, quietly = T, warn.conflicts = F)

# preparing files to proceed with the App...
folder <- ".data/"
talkien_version <- "v2.1"

# loading data
membRecept <- read.delim(paste0(folder, "membRecept_final.txt"), header = T, stringsAsFactors=F)
secrProt <- read.delim(paste0(folder, "secrProt_final.txt"), header = T, stringsAsFactors=F)
descr <- read.delim(paste0(folder, "genes_peptides.txt"), header = T, stringsAsFactors=F)
ppiGen_graph <- readRDS(paste0(folder, "ppiGen_graph_sub.rds"))
clx_tum <- read.table(paste0(folder, "tumor_sub.txt"), header = T, stringsAsFactors = F)
clx_norm <- read.table(paste0(folder, "norm_sub.txt"), header = T, stringsAsFactors = F)

# cleaning a bit the data
descr <- descr[!duplicated(descr$gene_symbol), ]
descr2 <- descr[descr$gene_symbol%in%secrProt$Gene | descr$gene_symbol%in%membRecept$Gene, ]

# getting entrez to pathwayID links and getting pathway names from pathwayIDs
xx <- as.list(reactome.db::reactomeEXTID2PATHID)
xy <- as.list(reactome.db::reactomePATHID2NAME)

# function to create links
createLink <- function(val) {
  sprintf('<a href="https://reactome.org/PathwayBrowser/#/%s" target="_blank" class="btn btn-primary">more info</a>',val)
}

# network descriptives
netsummary <- function(net) {
  nodes = vcount(net)
  edges = ecount(net)
  diameter = max(eccentricity(net, mode = "all"))
  comps = components(net)$no
  shortestpath = round(mean_distance(net, directed = F), 3)
  dens = round(graph.density(net), 3)
  av_neigh = round(mean(degree(net, v=V(net), mode = "all", loops = F)), 3)
  cc = round(transitivity(net, type = "average"), 3)
  centr = round(centr_degree(net)$centralization, 3)
  return(data.frame(Nodes = nodes, 
                    Edges = edges, 
                    Diameter = diameter, 
                    Shortest_Path = shortestpath,
                    Density = dens,
                    Average_Neighbors = av_neigh,
                    Clustering_Coefficient = cc,
                    Centralization = centr))
}

# network topologies
toposummary <- function(net) {
  cc <- transitivity(net, type = "local")
  cc[is.na(cc)] <- 0
  deg <- degree(net, v=V(net), mode = "all", loops = F, normalized = F)
  deg_norm <- degree(net, v=V(net), mode = "all", loops = F, normalized = T)
  clo <- closeness(net, vids = V(net), mode = "all", normalized = F)
  bet <- betweenness(net, v=V(net), directed = F, normalized = T)
  eig <- eigen_centrality(net, directed = F, scale = F)
  pr <- page_rank(net, vids = V(net), directed = F, weights = E(net)$weight)
  ecc <- eccentricity(net, vids = V(net), mode = c("all"))
  edgebet <- edge_betweenness(net, e = E(net), directed = FALSE)
  weights <- E(net)$weight
  node_at <- data.frame(degree = deg,
                        norm_degree = deg_norm,
                        clust_coef = cc,
                        cc = round(cc, 2),
                        closeness = clo,
                        betwenness = bet,
                        eigencent = eig$vector,
                        eccentricity = ecc,
                        pagerank = pr$vector,
                        pr = round(pr$vector, 3),
                        label = names(V(net)))
  return(node_at)
}


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
                 
                 # submit button to begin the analysis
                 actionButton("submit", 
                              label = "Load Data",
                              icon = icon("upload"),
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 br(),br(),br(),
                 
                 # buttons for annotation types and for network interactions. Only if example is not selected
                 conditionalPanel(condition = "input.clx == 0",
                                  radioButtons("annot", "Choose type annotation:",
                                               choices = c(GeneSymbol = "gene_symbol", 
                                                           entrez = "entrezgene", 
                                                           ensembl = "ensembl_id",
                                                           uniprot = "uniprot_id"), 
                                               selected = "gene_symbol")),
                 
                 
                 radioButtons("nettype", "Choose type of network:",
                              choices = c(Whole = "who",
                                          Crosstalk = "crt"),
                              selected = "who"),
                 
                 # score filter
                 sliderInput(inputId = "score_slider", 
                             label = "Score threshold", 
                             min = 150, 
                             max = 900, 
                             value = 400, 
                             step = 50),
                 
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
                                                           DH = "layout_with_dh"),
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
                 conditionalPanel(condition = "input.submit",  
                                  visNetworkOutput(outputId = "network1", width = "100%", height = "900px") %>% withSpinner(type = 8),
                                  h5(downloadButton('downloadNetwork', 'Download network as .html'), align = "center"),
                                  br(),
                                  h4(strong("Node's Pathways"), align = "center"),
                                  h5("click on any node to display all pathways on which is involved", align = "center"),
                                  span(htmlOutput("downstreamID"), style = "color:gray; font-size:20px"),
                                  dataTableOutput(outputId = "downstream"),
                                  h5(downloadButton('downloadPaths', "Download node's pathways list"), align = "center"),
                                  br(),
                                  h5(downloadButton('downloadOutgenes', 'Download out genes list'), align = "center"))),
        
        tabPanel("Network Interactions", icon = icon("list-alt", lib = "font-awesome"), value = 2,
                 span(htmlOutput("warning2"), style = "color:red; font-size:20px", align = "center"),
                 br(),
                 conditionalPanel(condition = "input.submit",
                                  h4(strong("Network Interaction List"), align = "center"),
                                  dataTableOutput(outputId = "netlist") %>% withSpinner(type = 8),
                                  h5(downloadButton('downloadNetTable', 'Download data as .tsv'), align = "center"))),
        
        tabPanel("Network Parameters", icon = icon("list-alt", lib = "font-awesome"), value = 3,
                 span(htmlOutput("warning3"), style = "color:red; font-size:20px", align = "center"),
                 br(),
                 conditionalPanel(condition = "input.submit",
                                  h4(strong("Network Descriptive Parameters"), align = "center"),
                                  span(htmlOutput("info2"), style = "color:gray; font-size:16px"),
                                  br(),
                                  tableOutput(outputId = "table1") %>% withSpinner(type = 8),
                                  br(),
                                  h4(strong("Network Topological Parameters"), align = "center"),
                                  dataTableOutput(outputId = "topolist") %>% withSpinner(type = 8),
                                  h5(downloadButton('downloadTopoTable', 'Download data as .tsv'), align = "center"))),
        
        tabPanel("Enrichment", icon = icon("list", lib = "font-awesome"), value = 4,
                 span(htmlOutput("warning4"), style = "color:red; font-size:20px", align = "center"),
                 br(),
                 conditionalPanel(condition = "input.tabs == '4' & input.submit",
                                  h4(strong("Enrichment plot"), align = "center"),
                                  plotOutput("oraplot"),
                                  h5(downloadButton('downloadPlot', 'Download plot as .png'), align = "center"),
                                  tags$hr(),
                                  h4(strong("Over Representation Analysis"), align = "center"),
                                  dataTableOutput(outputId = "oea") %>% withSpinner(type = 8),
                                  h5(downloadButton('downloadNetTable2', 'Download data as .tsv'), align = "center"))),
        
        tabPanel("Help", icon = icon("book", lib = "font-awesome"), value = 4, 
                 tags$iframe(src = "TALKIEN-help2.html", width = "100%", height = "800"))
        
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
  shinyjs::hide("oraplot")
  shinyjs::hide("oea")
  
  output$info2 <- renderText({ "" })
  
  
  # reactive data for input values
  data <- reactiveValues(df1 = NULL)
  lists <- reactiveValues(list1 = NULL)
  net <- reactiveValues(nodes = NULL)
  net2 <- reactiveValues(pruned_network = NULL)
  anot <- reactiveValues(annot_genes = NULL)
  out_genes <- reactiveValues(uncomGenes = NULL)
  visData <- reactiveValues(nodes = NULL)
  subnet <- reactiveValues(ppi_subgraph = NULL)
  
  observeEvent(input$submit, {
    
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
          need(input$file1 != "", "Please select a data set"),
          need(input$file2 != "", "Please select a data set")
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
    
    # sanity check. If header is not present or if empty list is uploaded, no analysis can be done
    if (nrow(data$df1)==0 | nrow(data$df2)==0) {
      
      output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, No header in lists or maybe no entries in lists")) })
      
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
      if (length(grep("quenya", data$df1))>0) {
        
        output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, Quenya annotation not supported!")) })
        
      } else {
      
        # getting non common genes between lists... the ones we want
        out_genes$uncomGenes <- c(data$df1, data$df2)
        
        # second check, if everything is common, no analysis is possible
        if (length(out_genes$uncomGenes)==0) {
          
          output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, all genes common between lists!")) })
          
        } else {
        
          # for the sake of clarity, base annotation is gene symbol
          data$df1 <- descr[match(data$df1, table = descr[, gene_annot]), "gene_symbol"]
          data$df2 <- descr[match(data$df2, table = descr[, gene_annot]), "gene_symbol"]
          
          # third check, format must be the same in both lists
          if(is.na(data$df1) || is.na(data$df2)) {
            
            output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, format not allowed")) })
            
          } else {
          
            # getting all genes needed to produce the net: genes in lists with functions as Membrane receptor prot or Secreted prot
            ppi <- c(data$df1, data$df2)
            ppi <- ppi[!duplicated(ppi)]
            ppi <- as.data.frame(ppi)
            names(ppi) <- "gene_symbol"
            ppi_mapped <- merge(ppi, descr, all = F)
            anot$annot_genes <- c(as.character(ppi_mapped$gene_symbol[ppi_mapped$gene_symbol%in%membRecept$Gene]),
                             as.character(ppi_mapped$gene_symbol[ppi_mapped$gene_symbol%in%secrProt$Gene]))
            
            # fourth check, needed receptor o secretor proteins... also in our list
            if(length(anot$annot_genes)==0) {
              
              output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, neither receptor nor secreted proteins found!")) })
              
            } else {
              
              # getting only interactions between receptors and or secretors
              ppi_mapped <- ppi_mapped[ppi_mapped$gene_symbol%in%anot$annot_genes, ]
              
              # subsetting interactions between our genes of interest
              ppi_mapped$gene_symbol <- as.character(ppi_mapped$gene_symbol)
              ppi_mapped$ensembl_id <- as.character(ppi_mapped$ensembl_id)
              ppi_mapped$uniprot_id <- as.character(ppi_mapped$uniprot_id)
              
              # transforming to igraph subnetwork
              ppi_mapped <- ppi_mapped[ppi_mapped$STRING_id%in%names(V(ppiGen_graph)), ]
              disnet <- V(ppiGen_graph)[ppi_mapped$STRING_id]
              ppi_subgraph <- induced_subgraph(ppiGen_graph, disnet)
              subnet$ppi_subgraph <- ppi_subgraph
              
              # fourth check, needed receptor o secretor proteins... also in our list
              if(length(E(ppi_subgraph))==0) {
                
                output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "Check your input data, no interactions found!")) })
                
              } else {
              
                # setting gene names as labels
                V(ppi_subgraph)$labels <- as.character(ppi_mapped[match(names(V(ppi_subgraph)), table=ppi_mapped$STRING_id), "gene_symbol"])
                
                # creating a data frame with empty node attributes
                net_attrib <- data.frame(gene = V(ppi_subgraph)$labels, fun = rep("", length(V(ppi_subgraph))), list = rep("", length(V(ppi_subgraph))), size = rep("", length(V(ppi_subgraph))), color = rep("", length(V(ppi_subgraph))), shape = rep("", length(V(ppi_subgraph))), stringsAsFactors = F)
                
                # ... populating it
                net_attrib$fun[net_attrib$gene%in%membRecept$Gene] <- "Receptor"
                net_attrib$fun[net_attrib$gene%in%secrProt$Gene] <- "Secreted"
                
                id1 <- net_attrib$gene[net_attrib$gene%in%membRecept$Gene]
                
                net_attrib$list[net_attrib$gene%in%data$df1] <- lists$list1
                net_attrib$list[net_attrib$gene%in%data$df2] <- lists$list2
                
                # transforming to visnetwork object
                visData <- toVisNetworkData(ppi_subgraph)
                net$nodes <- visData$nodes
                net$edges <- visData$edges
                
                # reformating visnetworkdata object for plots
                net$nodes$label <- net$nodes$labels
                net$nodes <- net$nodes[, -2]
                net$nodes <- merge(net$nodes, net_attrib, all = F, by.x = "label", by.y = "gene")
                
                # prune edges between secrprot--secrprot. To do this, for loop:
                pruned_network <- as.data.frame(matrix(nrow=nrow(net$edges), ncol=ncol(net$edges)))
                names(pruned_network) <- names(net$edges)
                count = 0
                
                for (i in 1:nrow(net$edges)) {
                  if (net$nodes[grep(net$edges$from[i], net$nodes$id), 3] == "Secreted" && net$nodes[grep(net$edges$to[i], net$nodes$id), 3] == "Secreted") {
                  } else {
                    count <- count + 1
                    pruned_network[count, ] <- net$edges[i, ]
                  }
                }
                
                # adding attributes
                pruned_network <- pruned_network[1:count, ]
                pruned_network$color <- rep("darkgrey", nrow(pruned_network))
                
                net2$pruned_network <- pruned_network
                
              } # end of fifth check
            } # end of fourth check
          } # end of third check
        } # end of second check
      } # end of first check
    } # end of sanity check
    
  })
  
  talkien <- observe({
    
    if (input$submit) {
      
      if (NROW(data$df1)==0 | NROW(data$df2)==0) {
      } else {
        # first check, just a freak joke... 
        if (length(grep("quenya", data$df1))>0) {
        } else {
          # second check, if everything is common, no analysis is possible
          if (length(out_genes$uncomGenes)==0) {
          } else {
            # third check, format must be the same in both lists
            if(is.na(data$df1) || is.na(data$df2)) {
            } else {
              # fourth check, needed receptor o secretor proteins... also in our list
              if(length(anot$annot_genes)==0) {
              } else {
                # fifth check, no interactions
                if(length(E(subnet$ppi_subgraph))==0) {
                } else {
                  
                  vert <- net$nodes
                  # default layout
                  lay <- input$grlay
                  
                  # aditional filter by score depending on user filter
                  pruned_network <- net2$pruned_network[net2$pruned_network$combined_score>=input$score_slider, ]
                  
                  # sixth check, selected interactions must be more than 0
                  if(nrow(pruned_network)==0) {
                    
                    output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                    output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                    output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                    output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found! - score filter too high -")) })
                    
                  } else {
                    
                    # more information to the table
                    pruned_network$width <- pruned_network$combined_score/min(pruned_network$combined_score)
                    pruned_network$interaction_typeA <- rep("", NROW(pruned_network))
                    pruned_network$interaction_typeB <- rep("", NROW(pruned_network))
                    pruned_network$interaction_locationA <- rep("", NROW(pruned_network))
                    pruned_network$interaction_locationB <- rep("", NROW(pruned_network))
                    pruned_network$ProteinA <- rep("", NROW(pruned_network))
                    pruned_network$ProteinB <- rep("", NROW(pruned_network))
                    
                    # computing interaction types
                    
                    gene1 <- vert[match(pruned_network$from, vert$id), c(1,3,4)]
                    gene2 <- vert[match(pruned_network$to, vert$id), c(1,3,4)]
                    pruned_network$interaction_typeA <- gene1$list
                    pruned_network$interaction_typeB <- gene2$list
                    pruned_network$interaction_locationA <- gene1$fun
                    pruned_network$interaction_locationB <- gene2$fun
                    pruned_network$ProteinA <- gene1$label
                    pruned_network$ProteinB <- gene2$label
                    
                    crtGenes <- c()
                    if (input$nettype == "crt") {
                      crtGenes <- unique(c(pruned_network$ProteinA, pruned_network$ProteinB))
                      lay <- input$grlay2
                      pruned_network <- pruned_network[!paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = "-")%in%paste(lists$list1, lists$list1, sep = "-"), ]
                      pruned_network <- pruned_network[!paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = "-")%in%paste(lists$list2, lists$list2, sep = "-"), ]
                      pruned_network <- pruned_network[!paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "-") =="Receptor-Receptor", ]
                      crtGenes <- crtGenes[!crtGenes%in%unique(c(pruned_network$ProteinA, pruned_network$ProteinB))]
                    }
                    
                    # seventh check, if crosstalk selected, still interactions found...
                    if(nrow(pruned_network)==0) {
                      
                      output$warning1 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      output$warning2 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      output$warning3 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "No interactions found!")) })
                      
                    } else {
                      
                      # formatting network plot
                      vert <- vert[vert$id %in% unique(c(pruned_network$from, pruned_network$to)), ]
                      vert$shape[vert$fun=="Receptor"] <- "triangleDown"
                      vert$shape[vert$fun=="Secreted"] <- "dot"
                      vert$color[vert$list==lists$list1] <- "coral"
                      vert$color[vert$list==lists$list2] <- "steelblue"
                      
                      # computing size of nodes
                      for (i in 1:NROW(vert)) {
                        vert$size[i] <- length(grep(vert[[2]][i], c(pruned_network$from, pruned_network$to)))
                      }
                      
                      # reformating size for better visualization, and other vars
                      vert$size <- as.numeric(vert$size)
                      vert$size <- round(5*log2(vert$size+1))
                      vert$title <- vert$label
                      
                      # getting final network tables
                      nodes <- descr[match(vert$label, table = descr$gene_symbol), c(4,3,5)]
                      nodes[,4:11] <- vert
                      nodes <- nodes[, c(4:11,1:3)]
                      nodes <- nodes[order(nodes$color, nodes$shape), ]
                      nodes <- nodes[!duplicated(nodes$label), ]
                      
                      # adding legend
                      lnodes <- data.frame(label = c(paste(lists$list1, "Receptor", sep = "\n"), 
                                                     paste(lists$list1, "Secreted", sep = "\n"), 
                                                     paste(lists$list2, "Receptor", sep = "\n"), 
                                                     paste(lists$list2, "Secreted", sep = "\n")), 
                                           shape = c("triangleDown", "dot", "triangleDown", "dot"), 
                                           color = c("coral", "coral", "steelblue", "steelblue"))
                      
                      # if not bipartite...
                      if (lay != "layout_as_bipartite") {
                        subnetwork1 <- visNetwork(nodes, pruned_network, main = paste(lists$list1, lists$list2, sep = " X "), width = "100%", height = "100%") %>%
                          visIgraphLayout(layout = lay, randomSeed = 1234) %>%
                          visNodes() %>%
                          visPhysics(stabilization = F)
                        subnetwork <- subnetwork1 %>%
                          visLegend(addNodes = lnodes, useGroups = F, width = 0.2, main = "Legend", position = "right") %>%
                          visOptions(highlightNearest = TRUE, nodesIdSelection = T, autoResize = T) %>%
                          visInteraction(navigationButtons = TRUE) %>%
                          visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}"
                          )
                        # if bipartite, a bit different
                      } else {
                        # conditionalpanel reformating
                        if(input$direct == T) {dir = "UD"} else {dir = "LR"}
                        # reformating stuff to show the bipartite net
                        if (input$bipnet == "list") {
                          nodes$level <- as.numeric(as.factor(nodes$list))
                          nodes$x <- as.numeric(as.factor(nodes$color))
                        } else {
                          nodes$level <- as.numeric(factor(nodes$fun, levels = c("Secreted", "Receptor")))
                          nodes$x <- as.numeric(factor(nodes$shape))
                        }
                        
                        nodes <- nodes[order(nodes$level, nodes$x), ]
                        
                        subnetwork1 <- visNetwork(nodes, pruned_network, width = "100%", height = "100%", main = paste(lists$list1, lists$list2, sep = " X ")) %>%
                          visHierarchicalLayout(direction = dir, levelSeparation = 1500, treeSpacing = 0.001, nodeSpacing = 25, blockShifting = T) %>%
                          visEdges(scaling = list(min=0.2, max=0.5))
                        subnetwork <- subnetwork1 %>%
                          visOptions(highlightNearest = T, nodesIdSelection = TRUE, autoResize = T) %>%
                          visLegend(useGroups = F, addNodes = lnodes, position = "right", zoom = T, main = "Legend", width = 0.1) %>%
                          visPhysics(stabilization = T) %>%
                          visInteraction(navigationButtons = TRUE) %>%
                          visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}"
                          )
                        
                        nodes <- nodes[, -c(10,11)]
                      }
                      
                      # getting final tables
                      pruned_network <- pruned_network[, c(23,24,21,22,19,20,16)]
                      
                      graph_stat <- graph_from_data_frame(pruned_network[, 1:2])
                      if (components(graph_stat)$no>1) {
                        output$info2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Disconected graph. Showing total nodes and edges, rest of the parameters computed only for giant component")) })
                      } else {
                        output$info2 <- renderText({ "" })
                      }
                      
                      netsum <- netsummary(graph_stat)
                      netsum[1:3] <- sprintf("%1.0f", netsum[1:3])
                      row.names(netsum) <- ""
                      netstat <- toposummary(graph_stat)
                      netstat <- merge(netstat, nodes, by = "label", all = F)
                      netstat <- netstat[order(netstat$degree, decreasing = T), c(1,2,6,7,8,4,9,16,17,12,13,14,15)]
                      netstat[, 4:7] <- round(netstat[, 4:7], digits = 3)
                      netstat[, 3] <- formatC(netstat[, 3], digits = 2, format = "E")
                      names(netstat) <- c("Gene_Symbol", "Degree", "Closeness", "Betweenness", "Eigen_Centrality", "Clustering_coef", "Eccentricity", "Location", "List", "Entrez_ID", "ensembl_ID", "uniprot_ID", "String_ID")
                      
                      
                      # enrichment analysis
                      enr = F
                      if (input$tabs == '4') {
                        de <- nodes$entrez
                        x <- ReactomePA::enrichPathway(gene = de, pvalueCutoff = 0.01, readable = T)
                        if(nrow(x)==0) {
                          output$warning4 <- renderUI({ HTML(paste('<br/>', '<br/>', "No enrichment terms found!")) })
                        } else {
                          enr = T
                          oratable <- as.data.frame(x@result)
                          oratable[, 5:7] <- apply(oratable[,5:7], 2, function(x) formatC(x, format = "e", digits = 2))
                          oratable$geneID <- gsub("/", " / ", oratable$geneID)
                          IDs <- oratable[,1]
                          reacplot <- function(){
                            enrichplot::dotplot(x, showCategory = 15)
                          }
                        }
                      }
                      
                      # genes not used
                      switch (input$annot,
                              "entrezgene" = {out_genes$uncomGenes <- out_genes$uncomGenes[!out_genes$uncomGenes%in%nodes$entrez]
                              crtGenes <- descr$entrezgene[match(crtGenes, descr$gene_symbol)]},
                              "gene_symbol" = {out_genes$uncomGenes <- out_genes$uncomGenes[!out_genes$uncomGenes%in%nodes$label]},
                              "ensembl_id" = {out_genes$uncomGenes <- out_genes$uncomGenes[!out_genes$uncomGenes%in%nodes$ensembl]
                              crtGenes <- descr$ensembl_id[match(crtGenes, descr$gene_symbol)]},
                              "uniprot_id" = {out_genes$uncomGenes <- out_genes$uncomGenes[!out_genes$uncomGenes%in%nodes$uniprot]
                              crtGenes <- descr$uniprot_id[match(crtGenes, descr$gene_symbol)]}
                      )
                      
                      com <- data.frame(gene=out_genes$comGenes, reason=rep("common in both lists", length(out_genes$comGenes)), stringsAsFactors = F)
                      func <- data.frame(gene=out_genes$uncomGenes, reason=rep("location (neither membRecept nor secreted)", length(out_genes$uncomGenes)), stringsAsFactors = F)
                      crt <- data.frame(gene=crtGenes, reason=rep("crosstalk (node is not involved in crosstalk)", length(crtGenes)), stringsAsFactors = F)
                      crt <- crt[!crt$gene%in%func$gene, ]
                      outGenes <- rbind.data.frame(com, func, crt)
                      
                      
                      #rendering outputs
                      output$network1 <- renderVisNetwork(subnetwork)
                      output$table1 <- renderTable(netsum, include.rownames = T, align = "r")
                      output$netlist <- renderDataTable(pruned_network, options = list(pageLength = 15))
                      output$topolist <- renderDataTable(netstat, options = list(pageLength = 10))
                      
                      output$warning1 <- renderText({ "" })
                      output$warning2 <- renderText({ "" })
                      output$warning3 <- renderText({ "" })
                      if (enr == T && nrow(x)!=0)  output$warning4 <- renderText({ "" })
                      
                      # downloadable objects
                      output$downloadNetwork <- downloadHandler(
                        filename = function() {paste("network", Sys.Date(), '.html', sep='')}, content = function(con) {
                          subnetwork %>% 
                            visLegend(useGroups = F, addNodes = lnodes, position = "right", zoom = T, main = "Legend", width = 0.2) %>%
                            visOptions(width = "100%", height = "900px", highlightNearest = T) %>%
                            visExport() %>%
                            visSave(con)})
                      
                      output$downloadOutgenes <- downloadHandler(
                        filename = function() {'outGenes.txt'}, content = function(file) {
                          write.table(outGenes, file, row.names = FALSE, quote = F, sep = "\t")})
                      
                      output$downloadNetTable <- downloadHandler(
                        filename = function() {'netlist.txt'}, content = function(file) {
                          write.table(pruned_network, file, row.names = FALSE, quote = F, sep = "\t")})
                      
                      output$downloadTopoTable <- downloadHandler(
                        filename = function() {'nodelist.txt'}, content = function(file) {
                          write.table(netstat, file, row.names = FALSE, quote = F, sep = "\t")})
                      
                      
                      # printing ORA results
                      if (enr == T && nrow(x)>0) {
                        shinyjs::show("oraplot")
                        shinyjs::show("oea")
                        output$oraplot <- renderPlot(reacplot())
                        oratable <- cbind(oratable, createLink(IDs))
                        colnames(oratable)[ncol(oratable)] <- "link"
                        output$oea <- DT::renderDataTable(
                          DT::datatable(oratable, rownames = F,
                                        options = list(pageLength = 10,
                                                                 autoWidth = T,
                                                                 columnDefs = list(list(targets = 6, width = '150px')), scrollX = T),
                                                      escape = F))
                        
                        output$downloadPlot = downloadHandler(
                          filename = function() {'oraplot.png'},
                          content = function(file) {
                            png(file, width = 1200, height = 700)
                            print(reacplot() +
                                    ggtitle("\nEnrichment Plot") +
                                    theme(plot.title = element_text(hjust = 0.5)))
                            dev.off()
                          })
                        
                        output$downloadNetTable2 <- downloadHandler(
                          filename = function() {'oratable.txt'}, content = function(file) {
                            oratable[,10] <- as.character(oratable[,10])
                            oratable[,10] <- substr(oratable[,10], 10, nchar(oratable[,10])-55)
                            write.table(oratable, file, row.names = FALSE, quote = F, sep = "\t")
                          })  
                      } else {
                        shinyjs::hide("oraplot")
                        shinyjs::hide("oea")
                      }
                      
                      
                      # the trick for downstream analysis...
                      observe({
                        if(length(input$click) == 1) {
                          shinyjs::show("downstream")
                          shinyjs::show("downstreamID")
                          shinyjs::show("downloadPaths")
                          visNetworkProxy("network1") %>%
                            visGetSelectedNodes(input = "connected_nodes")
                          
                          ent_sel <- descr$entrezgene[grep(paste("^", input$connected_nodes, "$", sep = ""), descr$STRING_id)]
                          sym_sel <- descr$gene_symbol[grep(paste("^", input$connected_nodes, "$", sep = ""), descr$STRING_id)]
                          hsa <- xx[grep(paste("^", as.character(ent_sel), "$", sep = ""), names(xx))]
                          hsa <- names(xy[names(xy)%in%unname(unlist(hsa))])
                          ids <- gsub("Homo sapiens: ", "", unlist(unname(xy[names(xy)%in%unname(unlist(hsa))])))
                          reaclink <- createLink(hsa)
                          
                          output$downstreamID <- renderUI({ HTML(paste('<br/>', '<br/>', "<center>", "<b>", sym_sel, "</b>", "</center>"))})
                          output$downstream <- renderDataTable({
                            
                            downpaths <<- cbind(hsa, ids, reaclink)
                            colnames(downpaths) <- c('ReactomeID', "Pathway_name", "link")
                            return(downpaths)
                            
                          }, escape = FALSE)
                          
                          # downloadable dynamic object
                          output$downloadPaths <- downloadHandler(
                            filename = function() {paste(sym_sel, 'pathways.tsv', sep="_")}, content = function(file) {
                              downpaths[,3] <- substr(downpaths[,3], 10, nchar(downpaths[,3])-55)
                              write.table(downpaths, file, row.names = FALSE, quote = F, sep = "\t")
                            })
                          
                        } else {
                          shinyjs::hide("downstream")
                          shinyjs::hide("downstreamID")
                          shinyjs::hide("downloadPaths")
                        }
                      })
        
                    } # end of the seventh check
                  } # end of the sixth check
                } # end of the fifth check
              } # end of fourth check
            } # end of third check
          } # end of second check
        } # end of first check
      } # end of sanity check
    }
    
  })
  
  
  }

shinyApp(ui, server)
