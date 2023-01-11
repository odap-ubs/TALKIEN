require(shiny, quietly = T, warn.conflicts = F)
require(shinycssloaders, quietly = T, warn.conflicts = F)
require(shinyjs, quietly = T, warn.conflicts = F)
require(shinyWidgets, quietly = T, warn.conflicts = F)
require(dplyr, quietly = T, warn.conflicts = F)
require(igraph, quietly = T, warn.conflicts = F)
require(visNetwork, quietly = T, warn.conflicts = F)
require(ggplot2, quietly = T, warn.conflicts = F)
require(DT, quietly = T, warn.conflicts = F)
require(ReactomePA, quietly = T, warn.conflicts = F)
require(RColorBrewer, quietly = T, warn.conflicts = F)
require(scales, quietly = T, warn.conflicts = F)

# status v1: working with few warnings. 
# static plot more or less fixed with visExport at higher resolution and download with rightclick
# DGIbd implemented, todo: links to chembl and PMIDs
# enrichments now not updating automatically... careful with swithces, at switching off, last input value keeps stored!
# displayed components done. Now only implemented on plot (and enrichment), but it might be implemented on tables
# inspect coloring behavior at downstream click: STATUS --- at 28/10/22 seems ok...


# preparing files to proceed with the App...
folder <- "./utils/"
talkien_version <- "v1.0"

# loading data
load(paste0(folder, "geneAnnot.Rdata"))
load(paste0(folder, "PPIs.Rdata"))
load(paste0(folder, "reactomeIDs.Rdata"))
load(paste0(folder, "LigRecInts.Rdata"))
load(paste0(folder, "ctnets.Rdata"))
load(paste0(folder, "DGBidb.Rdata"))
source(paste0(folder, "secondaryFunctions.R"))

# patch to put in reactomeFilters
ezLigRec <- idResource$entrez.id[match(V(crossTalkNetIgraph)$labels, idResource$gene.symbol)]
ezID2pathLR <- ezID2path[names(ezID2path)%in%ezLigRec]

# loading example files
samp <- read.delim("input_data/clx_input.txt", header = T, stringsAsFactors = F)

# user interface defining
ui <- fluidPage(
  
  # for showing hidding objects
  useShinyjs(),
  # title of the App
  titlePanel(title=div(img(src="talkien.jpg", height = 80, width = 320), " -crossTALK IntEraction Network-"), windowTitle="Talkien -crossTALK IntEraction Network-"),
  
  #layout and input options
  sidebarLayout(
    
    sidebarPanel(width = 3,
                 
                 talkien_version,
                 
                 # upload bed file option
                 h5("Please input one file with header and two columns length, first should contain genes / proteins, second cell types / lists to which it belongs to."),
                 
                 # toy example
                 checkboxInput(inputId = "clx", label = "Load example data", 
                               value = FALSE),
                 
                 # upload file
                 conditionalPanel(condition = "input.clx == 0",
                                  fileInput("file", "Input file:",
                                            multiple = FALSE,
                                            accept = "tsv")),
                 
                 # starting list array
                 uiOutput("arraySel"),
                 
                 # buttons for annotation types and for network interactions. Only if example is not selected
                 selectInput(inputId = "netint",
                             label = "Choose Ligand-Receptor DataBase:",
                             choices = list("CellChat", "CellPhoneDB", "iCellNet", "Ramilowsky", "Integrated"),
                             selected = "Integrated"),
                 
                 br(),

                 
                 conditionalPanel(condition = "input.clx == 0",
                                  radioButtons("annot", "Choose type annotation:",
                                               choices = c(GeneSymbol = "gene.symbol", 
                                                           EntrezGene = "entrez.id", 
                                                           ProteinName = "protein.name",
                                                           UniProt = "uniprot.ID"), 
                                               selected = "gene.symbol")),
                 
                 
                 radioButtons("nettype", "Choose type of network:",
                              choices = c(Full = "full",
                                          Crosstalk = "crt"),
                              selected = "full"),
                 
                 # show downstream nodes
                 h5(strong("Show downstream interactions:")), 
                 
                 switchInput(inputId = "graylinks", label = "", value = FALSE, size = "mini"),
                 
                 # select kind of downstream interactions to be shown
                 conditionalPanel(condition = "input.graylinks",# | input.graylinks && input.netint != 'STRING_exp'",
                                  selectInput(inputId = "curation", 
                                              label = "Choose interaction curation:",
                                              choices = list("HighScore", "Experimental"),
                                              selected = NULL)),
                 
                 # score filter
                 conditionalPanel(condition = "input.graylinks & input.curation == 'HighScore'",
                                  sliderInput(inputId = "score_sliderHi", 
                                              label = "Score threshold", 
                                              min = 700, 
                                              max = 950, 
                                              value = 950, 
                                              step = 50)),
                 
                 conditionalPanel(condition = "input.graylinks & input.curation == 'Experimental'",
                                  sliderInput(inputId = "score_sliderEx", 
                                              label = "Score threshold", 
                                              min = 300, 
                                              max = 950, 
                                              value = 300, 
                                              step = 50)),
                 
                 # network layouts
                 conditionalPanel(condition = "input.tabs == '1'",
                                  conditionalPanel(condition = "input.nettype == 'crt'",
                                                   radioButtons("grlay2", "Graphical Layout:",
                                                                choices = c(Circular = "layout_in_circle",
                                                                            ForceDirected = "layout_with_kk",
                                                                            Sphere = "layout_on_sphere",
                                                                            MDS = "layout_with_mds",
                                                                            Bipartite = "layout_as_bipartite"),
                                                                selected = "layout_as_bipartite")),
                                  
                                  conditionalPanel(condition = "input.nettype == 'full'",
                                                   radioButtons("grlay", "Graphical Layout:",
                                                                choices = c(Circular = "layout_in_circle",
                                                                            ForceDirected = "layout_with_kk", 
                                                                            Sphere = "layout_on_sphere",
                                                                            MDS = "layout_with_mds",
                                                                            largeGraph = "layout_with_lgl"),
                                                                selected = "layout_with_kk")))
    ),
    
    # main panel for output. 4 main tabs and help tab
    mainPanel(
      
      tabsetPanel(id = "tabs",
                  tabPanel("Plot", icon = icon("project-diagram", lib = "font-awesome"), value = 1,
                           span(htmlOutput("warning1"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           dropdownButton(radioButtons("centr", "node size by:",
                                                       choices = c(Degree = "deg",
                                                                   Closeness = "clo",
                                                                   Betweenness = "bet",
                                                                   EigenCentrallity = "eig",
                                                                   PageRank = "pr"),
                                                        selected = "deg"),
                                          uiOutput("showComp"),
                                          conditionalPanel(condition = "input.grlay2 == 'layout_as_bipartite' & input.nettype == 'crt'",
                                                           radioButtons("bipnet", "bipartite by:",
                                                                        choices = c(list = "list", 
                                                                                    location = "location")),
                                                           materialSwitch("direct", "Up-Down",
                                                                          status = "primary",
                                                                          right = F)),
                                          inputId = "filterWiew", circle = T, size = "s", tooltip = T, label = "Customize plot", icon = icon("cogs"), status = "primary"),
                           visNetworkOutput(outputId = "network1", width = "100%", height = "900px") %>% withSpinner(type = 8),
                           h5(downloadButton('downloadNetHTML', 'Download network as .html'), align = "center"),
                           h5(downloadButton('downloadNetCSV', 'Download network as .csv'), align = "center"),
                           h5(downloadButton('downloadNetGML', 'Download network as .graphML'), align = "center"),
                           br(),
                           h4(strong("Node's Pathways"), align = "center"),
                           h5("click on any node to display all pathways on which is involved", align = "center"),
                           span(htmlOutput("warning1_2"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           span(htmlOutput("downstreamID"), style = "color:gray; font-size:20px"),
                           dataTableOutput(outputId = "downstream"),
                           h5(downloadButton('downloadPaths', "Download node's pathways list"), align = "center")),
                  
                  tabPanel("Network Interactions", icon = icon("list-alt", lib = "font-awesome"), value = 2,
                           span(htmlOutput("warning2"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           h4(strong("Network Interaction List"), align = "center"),
                           dataTableOutput(outputId = "netlist") %>% withSpinner(type = 8),
                           h5(downloadButton('downloadNetTable', 'Download data as .tsv'), align = "center")),
                  
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
                           h5(downloadButton('downloadTopoTable', 'Download data as .tsv'), align = "center")),
                  
                  tabPanel("Enrichment", icon = icon("list", lib = "font-awesome"), value = 4,
                           span(htmlOutput("warning4"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           conditionalPanel(condition = "input.tabs == '4'",
                                                             dropdownButton(
                                                               conditionalPanel(condition = "input.graylinks",
                                                                                h5(strong("Add downstream nodes:")),
                                                                                switchInput("globEnr",
                                                                                             value = F,
                                                                                             size = "mini")),
                                                               radioButtons("showCat", "Number of showed pathways:",
                                                                            choices = c("10" = 10, 
                                                                                        "15" = 15,
                                                                                        "20" = 20,
                                                                                        "ALL" = "all"),
                                                                            selected = "15"),
                                                               h5(strong("component-wise enrichment:")),
                                                               switchInput("comp",
                                                                           value = F,
                                                                           size = "mini"),
                                                               uiOutput("enrComp"),
                                                               h5(strong("enrichent by list:")),
                                                               switchInput("tissue",
                                                                           value = F,
                                                                           size = "mini"),
                                                               uiOutput("enrTissue"),
                                                           inputId = "enrPlotOpt", circle = T, size = "s", tooltip = T, label = "Customize options", icon = icon("cogs"), status = "primary"),
                                            h4(strong("Enrichment plot"), align = "center"),
                                            plotOutput("oraplot") %>% withSpinner(type = 8),
                                            h5(downloadButton('downloadPlot', 'Download plot as .png'), align = "center"),
                                            tags$hr(),
                                            h4(strong("Over Representation Analysis"), align = "center"),
                                            dataTableOutput(outputId = "oea") %>% withSpinner(type = 8),
                                            h5(downloadButton('downloadNetTable2', 'Download data as .tsv'), align = "center"))),
                  
                  tabPanel("Drugs", icon = icon("list-alt", lib = "font-awesome"), value = 5,
                           span(htmlOutput("warning5"), style = "color:red; font-size:20px", align = "center"),
                           br(),
                           h4(strong("Receptor-Drug Interaction"), align = "center"),
                           dataTableOutput(outputId = "drugs") %>% withSpinner(type = 8),
                           h5(downloadButton('downloadDrugTable', 'Download data as .tsv'), align = "center")),
                  
                  tabPanel("Help", icon = icon("book", lib = "font-awesome"), value = 6, 
                           tags$iframe(src = "TALKIEN_readme.html", width = "100%", height = "800"))
                  
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
  drint <- reactiveValues(recs = NULL)
  cols <- brewer.pal(n=8, name="Set2")
  
  netIni <- reactive({

    netIni <- list()
    data <- lists <- gene_anot <- NULL
    
    if (input$clx == F) {   
      tryCatch({
        data <- read.delim(input$file$datapath,
                           header = T,
                           stringsAsFactors = F)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        validate(
          need(input$file != "", "")
        )
        stop(safeError(e))
      })
      gene_annot <- input$annot
    } else {
      data <- samp
      gene_annot <- "gene.symbol"
    }
    
    if (NCOL(data)<2) {
      
      output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnColumn)) })
      netIni <- list(data)
      
    } else {

      if (nrow(data)==0 | nrow(data[!is.na(data[,1]), ])==0) {
        
        output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnEntries)) })
        netIni <- list(data)
        
      } else {
        
        # beggining to analyze data. first, removing duplicate ids
        names(data)[1:2] <- c("entries", "list")
        data <- data[!is.na(data$entries), ]
        lists <- sort(unique(data$list))
        # zero check, number of lists (or column two is not a list variable)
        if (length(lists)>8) {
          
          output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnMaxlist)) })          
          netIni <- list(data)
          
        } else {

          netIni <- list(data, gene_annot, lists)
          
        }
      }
    }
  })
  
  output$arraySel <- renderUI({

    netIni <- netIni()
    req(length(netIni)==3)

        pairedLists <- combn(netIni[[3]], m = 2, simplify = T)
        pairedLists <- apply(pairedLists, 2, function(x) paste(x, collapse = " X "))
        pickerInput(inputId = "list", 
                    label = "Choose interacting lists", 
                    choices = pairedLists, 
                    multiple = T, 
                    selected = pairedLists, options = list(`none-selected-text` = "Please provide at least one pairwise comparison",
                                                           `actions-box` = TRUE))

  })
  
  net <- reactive({

    net <- list()
    netIni <- netIni()
    
    req(length(netIni)==3, input$list)

      data <- netIni[[1]]
      gene_annot <- netIni[[2]]
      selLists <- lists <- netIni[[3]]
      anot <- out_genes <- arrayList <- NULL

      if (!is.null(input$list)) {
        pairwiseList <- strsplit(input$list, split = " X ")
        selLists <- unique(unlist(pairwiseList))
        arrayList <- input$list
        if (any(!selLists%in%lists)) {
          selLists <- lists
          arrayList <- combn(selLists, m = 2, simplify = T)
          arrayList <- apply(arrayList, 2, function(x) paste(x, collapse = " X "))
        }

      }

      data[,"entries"] <- gsub(" .*", "", data[, "entries"])

      data <- filter(data, list %in% selLists)

      data <- suppressMessages(reshape2::dcast(data = data, entries~list, fill = 0, fun.aggregate = length))
      row.names(data) <- data$entries
      data$entries <- NULL
      data[data!=0] <- 1
      filtCom <- apply(data, 1, sum)

      # getting common genes in both lists and removing them
      out_genes$comGenes <- row.names(data)[filtCom>1]
      data <- data[filtCom==1, ]
      
      # first check, just a freak joke... 
      if (length(grep("quenya", row.names(data)))>0) {
        
        output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnQuenya)) })        
        net <- list(data)
        
      } else {
        
        # getting non common genes between lists... the ones we want
        out_genes$uncomGenes <- row.names(data)
        
        # second check, if everything is common, no analysis is possible
        if (length(out_genes$uncomGenes)==0) {
          
          output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnCommon)) })          
          net <- list(data, out_genes)
          
        } else {
          
          # for the sake of clarity, base annotation is gene symbol
          valEntries <- idResource[match(row.names(data), table = idResource[, gene_annot]), "gene.symbol"]
          data <- data[!is.na(valEntries), ]
          row.names(data) <- valEntries[!is.na(valEntries)]
          
          # third check, format must be the same in both lists
          if(NROW(data)==0) {
            
            output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnFormat)) })
            net <- list(data, out_genes)
            
          } else {
            
            for (i in 1:ncol(data)) { 
              tmp <- data[,i]==1
              data[, i] <- ifelse(tmp, names(data)[i], NA)
            }
            entryList <- apply(data, 2, function(x) row.names(data)[!is.na(x)])
            data$list <- gsub("NA", "", apply(data, 1, paste, collapse = ""))
            
            # getting all genes needed to produce the net: genes in lists with functions as Membrane receptor prot or Secreted prot
            ppi <- data.frame(gene.symbol=row.names(data))
            ppi_mapped <- left_join(ppi, idResource, by = "gene.symbol")
            anot$annot_genes <- ppi_mapped$gene.symbol[ppi_mapped$gene.symbol%in%row.names(ddbbs_merged)]
            
            # fourth check, needed receptor o secretor proteins... also in our list
            if(length(anot$annot_genes)==0) {
              
              output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnAnnot)) })              
              net <- list(data, out_genes, anot)
              
            } else {
              
              # transforming to igraph subnetwork... and the same for downstream network
              ppi_downstream <- as.undirected(ppi_hiscore, mode = "each") %>% simplify(remove.multiple = T, edge.attr.comb = max)
              
              if (input$graylinks && input$curation == "Experimental") {
                ppi_downstream <- as.undirected(ppi_exp, mode = "each") %>% simplify(remove.multiple = T, edge.attr.comb = max)
              }
              
              switch (input$netint,
                      "CellChat" = { 
                        ppiGen_graph <- crossTalkNetCellChat
                        nodeClass <- data.frame(nodeClass = ddbbs_annot$cellchat, gene.symbol = row.names(ddbbs_annot))
                      },
                      "CellPhoneDB" = { 
                        ppiGen_graph <- crossTalkNetCellPhone
                        nodeClass <- data.frame(nodeClass = ddbbs_annot$cellphone, gene.symbol = row.names(ddbbs_annot))
                      },
                      "iCellNet" = { 
                        ppiGen_graph <- crossTalkNetiCellNet
                        nodeClass <- data.frame(nodeClass = ddbbs_annot$icellnet, gene.symbol = row.names(ddbbs_annot))
                      },
                      "Ramilowsky" = { 
                        ppiGen_graph <- crossTalkNetRam
                        nodeClass <- data.frame(nodeClass = ddbbs_annot$ramilowsky, gene.symbol = row.names(ddbbs_annot))
                      },
                      "Integrated" = { 
                        ppiGen_graph <- crossTalkNetIgraph
                        nodeClass <- data.frame(nodeClass = ddbbs_annot$consensus, gene.symbol = row.names(ddbbs_annot))
                      }
              )
              
              # getting only interactions between receptors and or secretors, the other ones to downstream purposes
              nodeClass <- nodeClass[!is.na(nodeClass$nodeClass), ]
              
              ppi_mappedDS <- ppi_mapped[ppi_mapped$string.id%in%names(V(ppi_downstream)), ]
              
              ppi_mapped <- ppi_mapped[ppi_mapped$gene.symbol%in%nodeClass$gene.symbol, ]
              ppi_mapped <- ppi_mapped[ppi_mapped$string.id%in%names(V(ppiGen_graph)), ]
              ppi_mapped$anot <- nodeClass$nodeClass[match(ppi_mapped$gene.symbol, nodeClass$gene.symbol)]

              net <- list(data, arrayList, anot, lists, ppi_mapped, ppiGen_graph, nodeClass, ppi_downstream, ppi_mappedDS, out_genes)
              
            }
          }
        }
      }
  })
  
  net_clean <- reactive({
    
    net_clean <- list()
    net <- net()
    
    # checking previous reactive worked fine
    req(length(net)==10)
      
      data <- net[[1]]
      arrayList <- net[[2]]
      anot <- net[[3]]
      lists <- net[[4]]
      ppi_mapped <- net[[5]]
      ppiGen_graph <- net[[6]]
      nodeClass <- net[[7]]
      ppi_downstream <- net[[8]]
      ppi_mappedDS <- net[[9]]
      out_genes <- net[[10]]
    
      # initializing enrichment analysis
      enr <- NULL
      
      # get a network of CT interactions and removing isolated nodes
      # WARNING !!! check as.undirected is not affecting union graphs!!!
      disnet <- names(V(ppiGen_graph))[names(V(ppiGen_graph))%in%ppi_mapped$string.id]
      ppi_subgraph <- induced_subgraph(as.undirected(ppiGen_graph), disnet)
      ppi_subgraph <- delete.vertices(ppi_subgraph, degree(ppi_subgraph)==0)

      # fifth check, needed receptor or ligand proteins... also in our list
      if(length(E(ppi_subgraph))==0) {
        
        output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnInt)) })        
        net_clean <- list(ppi_subgraph)
        
      } else {
        
        # removing isolated nodes from nodes dataframe
        E(ppi_subgraph)$desc <- "ppi"
        V(ppi_subgraph)$fun <- ppi_mapped$anot[match(V(ppi_subgraph)$labels, ppi_mapped$gene.symbol)]
        ppi_mapped <- ppi_mapped[ppi_mapped$gene.symbol%in%V(ppi_subgraph)$labels, ]
        
        # get a network of PPI interactions with nodes included in input lists and removing isolated nodes
        ppi_mappedDS$anot <- ppi_mapped$anot[match(ppi_mappedDS$gene.symbol, ppi_mapped$gene.symbol)]
        ppi_mappedDS$anot[is.na(ppi_mappedDS$anot)] <- "downstream"
        ppi_mappedDS <- ppi_mappedDS[ppi_mappedDS$anot!="Ligand", ]
        
        disnetDS <- names(V(ppi_downstream))[names(V(ppi_downstream))%in%ppi_mappedDS$string.id]
        ppi_subgraphDS <- induced_subgraph(ppi_downstream, disnetDS)
        ppi_subgraphDS <- delete.vertices(ppi_subgraphDS, degree(ppi_downstream)==0)
        # getting node class
        V(ppi_subgraphDS)$labels <- ppi_mappedDS$gene.symbol[match(V(ppi_subgraphDS)$name, ppi_mappedDS$string.id)]
        V(ppi_subgraphDS)$fun <- ppi_mappedDS$anot[match(V(ppi_subgraphDS)$labels, ppi_mappedDS$gene.symbol)]
        E(ppi_subgraphDS)$desc <- "down"

        ppi_union <- igraph::union(ppi_subgraph, ppi_subgraphDS)
        V(ppi_union)$labels <- idResource$gene.symbol[match(names(V(ppi_union)), idResource$string.id)]
        V(ppi_union)$fun <- V(ppi_union)$fun_2
        V(ppi_union)$fun[is.na(V(ppi_union)$fun)] <- V(ppi_union)$fun_1[is.na(V(ppi_union)$fun)]
        E(ppi_union)$desc <- E(ppi_union)$desc_2
        E(ppi_union)$desc[is.na(E(ppi_union)$desc)] <- "main"
        
        ppi_union <- delete_vertex_attr(ppi_union, name= "fun_1") %>% 
          delete_vertex_attr(name="fun_2") %>% 
          delete_vertex_attr(name="labels_1") %>% 
          delete_vertex_attr(name="labels_2") %>%
          delete_edge_attr(name="desc_1") %>% 
          delete_edge_attr(name="desc_2")
        
        # creating a data frame with empty node attributes
        net_attrib <- data.frame(gene = V(ppi_union)$labels, 
                                 fun = V(ppi_union)$fun, 
                                 list = rep("", length(V(ppi_union))), 
                                 size = rep("", length(V(ppi_union))), 
                                 color = rep("", length(V(ppi_union))), 
                                 shape = rep("", length(V(ppi_union))),
                                 component = rep("", length(V(ppi_union))), stringsAsFactors = F)
        
        # ... populating it
        net_attrib$list <- data$list[match(net_attrib$gene, row.names(data))]
        
        # transforming to visnetwork object
        visData <- toVisNetworkData(ppi_union)
        net$nodes <- visData$nodes
        net$edges <- visData$edges

        # reformating visnetworkdata object for plots
        net$nodes <- net$nodes[, c("labels", "id")]
        net$nodes <- left_join(net$nodes, net_attrib, by = c("labels" = "gene"))
    
        # aditional filter by score depending on user filter
        pruned_network <- rbind.data.frame(net$edges[is.na(net$edges$weight), ], net$edges[net$edges$weight>=950 & !is.na(net$edges$weight), ])

        filtScore <- 950
        if (input$graylinks == T && input$curation == "HighScore") filtScore <- input$score_sliderHi
        if (input$graylinks == T && input$curation == "Experimental") filtScore <- input$score_sliderEx
        
        pruned_network <- rbind.data.frame(net$edges[is.na(net$edges$weight), ], net$edges[net$edges$weight>=filtScore & !is.na(net$edges$weight), ])
        
        # sixth check, selected interactions must be more than 0
        if(nrow(pruned_network)==0) {
          
          output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnScore)) })
          net_clean <- list(ppi_subgraph, pruned_network)
          
        } else {
          
          vert <- net$nodes[net$nodes$id%in%unique(c(pruned_network$from, pruned_network$to)), ]
          
          # computing interaction types
          gene1 <- vert[match(pruned_network$from, vert$id), c("labels", "fun", "list")]
          gene2 <- vert[match(pruned_network$to, vert$id), c("labels", "fun", "list")]
          
          # more information to the table
          pruned_network <- pruned_network %>% mutate(interaction_typeA = gene1$list) %>%
            mutate(interaction_typeB = gene2$list) %>%
            mutate(interaction_locationA = gene1$fun) %>%
            mutate(interaction_locationB = gene2$fun) %>%
            mutate(ProteinA = gene1$label) %>%
            mutate(ProteinB = gene2$label)
          
          # removing secreted - secreted, downstream - secreted, secreted - downstream (can be modified easily)
          int <- paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")
          pruned_network <- pruned_network[!int %in% c("Ligand_Ligand", "Ligand_downstream", "downstream_Ligand"), ]
          int <- int[!int %in% c("Ligand_Ligand", "Ligand_downstream", "downstream_Ligand")]
          # downstream proteins must be from the same list of their receptor activators, so prune edges between receptor-list1 and downstream-list2 and viceversa
          # we switch order from receptor-downstream interaction and receptor-ligand interactions
          pruned_network[int%in%"downstream_Receptor", ] <- pruned_network[int%in%"downstream_Receptor", c("to", "from", "weight", "desc", "interaction_typeB", "interaction_typeA", "interaction_locationB", "interaction_locationA", "ProteinB", "ProteinA")]
          pruned_network[int%in%"Receptor_Ligand", ] <- pruned_network[int%in%"Receptor_Ligand", c("to", "from", "weight", "desc", "interaction_typeB", "interaction_typeA", "interaction_locationB", "interaction_locationA", "ProteinB", "ProteinA")]
          # finally, removing pairwise interaction not selected

          if (!is.null(input$list)) {
              pairwiseList <- strsplit(arrayList, split = " X ")
              filtList <- unique(unlist(pairwiseList))
              pairwiseRevList <- lapply(pairwiseList, rev)
              pairwiseRevList <- unlist(lapply(pairwiseRevList, function(x) paste(x, collapse = " X ")))
              allFiltList <- c(arrayList, pairwiseRevList, paste(filtList, filtList, sep = " X "))
              filtPairwiseList <- paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = " X ")
              pruned_network <- pruned_network %>% filter(filtPairwiseList%in%allFiltList)
              lists <- lists[lists%in%filtList]
          }
          
          int <- paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")
          lis <- paste(pruned_network$interaction_typeA, pruned_network$interaction_typeB, sep = "_")
          lis <- unlist(lapply(strsplit(lis, split = "_"), function(x) length(x[!duplicated(x)])))
          recs <- unique(pruned_network %>% filter(interaction_locationB=="Receptor") %>% dplyr::select(ProteinB))
          recs2rem <- unique(pruned_network %>% filter(interaction_locationA=="Receptor") %>% dplyr::select(ProteinA))
          recs2rem <- recs2rem$ProteinA[!recs2rem$ProteinA%in%recs$ProteinB]
          
          rem_edges <- int=="downstream_downstream" & lis==2 | int=="Receptor_downstream" & lis==2 | int=="Receptor_Receptor" & lis==2
          
          pruned_network <- pruned_network %>% filter(!rem_edges) %>% filter(!ProteinA %in% recs2rem)
          # double links on downstream downstream and receptor receptor
          bidir_downs <- pruned_network[paste(pruned_network$interaction_locationA, pruned_network$interaction_locationB, sep = "_")%in%c("downstream_downstream", "Receptor_Receptor"), c("to", "from","weight","desc", "interaction_typeB", "interaction_typeA", "interaction_locationB", "interaction_locationA", "ProteinB", "ProteinA")]
          names(bidir_downs) <- names(pruned_network)
          
          pruned_network <- rbind.data.frame(pruned_network, bidir_downs)
          
          # removing components with downstream genes only... to remove if working fine
          removedComps1 <- remDownComp(edgeList = pruned_network, nodeList = vert)
          pruned_network <- removedComps1[[1]]
          pruned_network$color <- "darkgray"
          vert <- removedComps1[[2]]
          
          downcross_network <<- pruned_network
          
          if (input$nettype == "crt") {

            main_edges <- downcross_network[downcross_network$desc=="main" & downcross_network$interaction_typeA != downcross_network$interaction_typeB & downcross_network$interaction_locationA != downcross_network$interaction_locationB, ]
            down_edges <- downcross_network[downcross_network$desc=="down", ]
            downcross_network <- rbind.data.frame(main_edges, down_edges) %>% 
              filter(!(interaction_locationA=="Receptor" & interaction_locationB=="Receptor"))
            crtrecs <- unique(downcross_network %>% filter(interaction_locationB=="Receptor") %>% dplyr::select(ProteinB))
            crt2rem <- unique(downcross_network %>% filter(interaction_locationA=="Receptor") %>% dplyr::select(ProteinA))
            crt2rem <- crt2rem$ProteinA[!crt2rem$ProteinA%in%crtrecs$ProteinB]
            downcross_network <- downcross_network %>% filter(!ProteinA %in% crt2rem)
            vert <- vert[vert$id%in%unique(c(downcross_network$from, downcross_network$to)), ]

            # seventh check, if crosstalk selected, still interactions found...
            if(nrow(downcross_network[downcross_network$desc=="main", ])==0) {
              
              output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnCrosst)) })
              downcross_network <- downcross_network[downcross_network$desc=="main", ]
              net_clean <- list(ppi_subgraph, pruned_network)
              
            } else {
              
              removedComps2 <- remDownComp(edgeList = downcross_network, nodeList = vert)
              downcross_network <- removedComps2[[1]]
              vert <- removedComps2[[2]]
              # end of the seventh check
            }
          }
          
          if(nrow(downcross_network)==0) {
            
            output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderUI({ HTML(paste(warnCrosst)) })
            net_clean <- list(ppi_subgraph, pruned_network)
            
          } else {
            
            output$warning1 <- output$warning2 <- output$warning3 <- output$warning4 <- output$warning5 <- renderText({ "" })
            
            # formatting network plot
            cols <- data.frame(list=lists, col=cols[1:length(lists)])
            #cols <- data.frame(list=lists, col = hue_pal()(length(lists)))
            vert$shape <- "square"
            vert$shape[vert$fun=="Receptor"] <- "triangleDown"
            vert$shape[vert$fun=="Ligand"] <- "dot"
            vert$color <- cols$col[match(vert$list, cols$list)]
            vert$color[vert$fun == "downstream"] <- "lightgray"
            
            # getting final network tables
            nodes <- idResource[match(vert$labels, table = idResource$gene.symbol), c("uniprot.id", "protein.id", "entrez.id")]
            nodes <- cbind.data.frame(vert, nodes)
            nodes <- arrange(nodes, list,fun,labels,)
            
            # adding legend
            legendList <- rep(lists, each=2)
            legendType <- rep(c("Receptor", "Ligand"), length(lists))
            lnodes <- data.frame(label = paste(legendList, legendType, sep = "\n"), 
                                 shape = rep(c("triangleDown", "dot"), length(lists)), 
                                 color = rep(cols$col, each = 2), stringsAsFactors = F)
            
            # getting final tables
            pruned_network <- downcross_network[downcross_network$desc=="main", c("ProteinA", "ProteinB", "interaction_locationA", "interaction_locationB", "interaction_typeA", "interaction_typeB", "color")]
            downcross_network <- downcross_network[, c("ProteinA", "ProteinB", "interaction_locationA", "interaction_locationB", "interaction_typeA", "interaction_typeB", "color")]
            names(pruned_network)[1:2] <- names(downcross_network)[1:2] <- c("from", "to")
            
            graph_stat <- graph_from_data_frame(pruned_network, directed = F) %>% simplify(remove.multiple = T)
            graph_stat2 <- graph_from_data_frame(downcross_network, directed = F) %>% simplify(remove.multiple = T)
            
            if (input$graylinks) {
              plotGraph <- graph_stat2
              netstat <- suppressWarnings(toposummary(graph_stat2))
              netstat$components <- unname(components(graph_stat2)$membership)
              edges <- downcross_network
              edges <- anti_join(edges, bidir_downs, by = c("from" = "ProteinA", "to" = "ProteinB"))
              lnodes <- rbind(lnodes, c("Downstream", "square", "lightgrey"))
            } else {
              plotGraph <- graph_stat
              netstat <- suppressWarnings(toposummary(graph_stat))
              netstat$components <- unname(components(graph_stat)$membership)
              edges <- pruned_network
            }

            netstat <- sortComps(netstat)
            # computing size of nodes (here from downstream, not from pruned, to be discussed)
            size_vert <- size_vert <- degree(plotGraph, v = V(plotGraph), mode = "all")
            size_vert <- data.frame(id = names(size_vert), size = unname(size_vert))
            nodes$size <- size_vert$size[match(nodes$labels, table = size_vert$id)]
            
            # reformating size for better visualization, and other vars
            nodes$size <- round(5*log2(nodes$size+1))
            nodes$title <- nodes$label
            
            names(edges)[1:2] <- c("ProteinA", "ProteinB")
            edges <- dplyr::select(edges, -color)

            netstat <- left_join(netstat, nodes, by = "labels")
            netstat <- netstat[order(netstat$degree, decreasing = T), c("labels", "degree", "closeness", "betwenness", "eigencent", "clust_coef", "eccentricity", "components", "fun", "list", "entrez.id", "protein.id", "uniprot.id", "id")]#c(1,2,6,7,8,4,9,14,15,22,20,21,12)]
            netstat[, c("betwenness", "eigencent", "clust_coef")] <- round(netstat[, c("betwenness", "eigencent", "clust_coef")], digits = 3)
            netstat[, "closeness"] <- formatC(netstat[, "closeness"], digits = 2, format = "E")
            names(netstat) <- c("Gene_Symbol", "Degree", "Closeness", "Betweenness", "Eigen_Centrality", "Clustering_coef", "Eccentricity", "Component", "Location", "List", "Entrez_ID", "Protein_name", "UniProt_ID", "String_ID")
            names(nodes)[1:2] <- c("id", "string.id")
            
            # preparing prefix for downloaded objects
            tagFiles <- ifelse(input$clx, "example", substr(input$file$name, 1, nchar(input$file$name)-4))
            tagFiles <- paste(tagFiles, gsub(" ", "", paste(arrayList, collapse = "_")), sep = "_")
            tagFiles <- ifelse(input$netint != "Integrated", paste(tagFiles, input$netint, sep = "_"), tagFiles)
            tagFiles <- ifelse(input$nettype=="crt", paste(tagFiles, input$nettype, sep = "_"), tagFiles)
            tagFiles <- ifelse(input$graylinks, paste(tagFiles, input$curation, filtScore, sep = "_"), tagFiles)
            
            net_clean <- list(pruned_network, nodes, downcross_network, graph_stat, graph_stat2, lists, enr, netstat, edges, lnodes, tagFiles)

        }         # end of the sixth check
      }           # end of the fifth check
    }             # end of the sanity check
  })
  
  output$showComp <- renderUI({
    net_clean <- net_clean()
    req(length(net_clean)==11)
    netstat <- net_clean[[8]]
    comps2show <- "ALL"
    if (length(unique(netstat$Component))>1) { 
      comps2show <- c("ALL", sort(unique(netstat$Component)))
    }
    selectInput(inputId = "shcomp", label = "Show component:", choices = comps2show, selected = "ALL")
  })
  
  output$enrComp <- renderUI({
    net_clean <- net_clean()
    req(length(net_clean)==11, req(input$comp))
      netstat <- net_clean[[8]]
      comps <- sort(unique(netstat$Component))
      selectInput(inputId = "comps", label = "", choices = comps, selected = NULL)
  })
  
  output$enrTissue <- renderUI({
    net_clean <- net_clean()
    req(length(net_clean)==11, input$tissue)
      netstat <- net_clean[[8]]
      tissues <- as.list(sort(unique(netstat$List)))
      checkboxGroupInput(inputId = "tis", label = "", choices = tissues, selected = NULL)
  })
  
  ora <- reactive({

    net_clean <- net_clean()
    ora <- list()

    req(length(net_clean)==11)
    
      netstat <- net_clean[[8]]
      numComp <- "ALL"
      tis <- sort(unique(netstat$List))
      univ <- names(ezID2path)
      
      if (input$globEnr ==F) {
        netstat <- netstat[netstat$Location != "downstream", ]
        univ <- names(ezID2pathLR)
      }
      
      if (input$comp == T) {
        req(input$comps)
          numComp <- input$comps
          netstat <- netstat[netstat$Component%in%as.integer(numComp), ]
      }
      
      if (input$tissue == T) {
        req(input$tis)
          tis <- input$tis

        netstat <- netstat %>% filter(List %in% tis)
      }
   
      if (input$tabs == '4') {
        enr <- netstat$Entrez_ID
        enr <- ReactomePA::enrichPathway(gene = enr, pvalueCutoff = 0.05, universe = univ, readable = T)
        ora <- list(enr, numComp, tis)
    }
  })
  
  output$table1 <- renderTable({
    
    net_clean <- net_clean()

    req(length(net_clean)==11)
      
      nodes <- net_clean[[2]]
      pruned_network <- net_clean[[1]]
      downcross_network <- net_clean[[3]]
      graph_stat <- net_clean[[4]]
      graph_stat2 <- net_clean[[5]]
      
      if (components(graph_stat)$no>1) {
        output$info2 <- renderUI({ HTML(paste('<br/>', '<br/>', "Disconected graph. Showing total nodes and edges, rest of the parameters computed only for giant component")) })
      } else {
        output$info2 <- renderText({ "" })
      }
      
      netsum <- as.data.frame(netsummary(graph_stat), row.names = "")
      netsum2 <- as.data.frame(netsummary(graph_stat2), row.names = "")
      netsum[c("Nodes", "Edges", "Diameter")] <- sprintf("%1.0f", netsum[c("Nodes", "Edges", "Diameter")])
      netsum2[c("Nodes", "Edges", "Diameter")] <- sprintf("%1.0f", netsum2[c("Nodes", "Edges", "Diameter")])
      
      if (input$graylinks) {
        netsum <- rbind.data.frame(netsum, netsum2)
        row.names(netsum) <- c("Ligand-Receptor", "Downstream")
      }
      netsum
    
  }, rownames = T)
  
  output$netlist <- renderDataTable({

    net_clean <- net_clean()
    req(length(net_clean)==11)
    if (nrow(net_clean[[1]]) !=0 ) {
        datatable(net_clean[[9]], options = list(pageLength = 15))
    }
  })

  output$topolist <- renderDataTable({
    
    net_clean <- net_clean()
    req(length(net_clean)==11)
    
      datatable(net_clean[[8]])

  })
  
  output$drugs <- renderDataTable({
    
    net_clean <- net_clean()
    
    req(length(net_clean)==11)
    
      if (input$tabs == '5') {
        
        recs <- drint$recs
        recs <- net_clean[[2]] %>% filter(fun == "Receptor") %>% dplyr::select(id)
        recs <- inner_join(recs, dgbi, by = c("id" = "gene_name")) %>% 
          dplyr::arrange(dplyr::desc(interaction_group_score), drug_name)
        drint$recs <- recs
        
        if(nrow(recs)==0) { 
          output$warning5 <- renderUI({ HTML(paste(warnDrug)) })
        } else {
          datatable(recs)
        }
      }
  })

  output$network1 <- renderVisNetwork({
    
    net_clean <- net_clean()
    
    req(length(net_clean)==11)

    pruned_network <- net_clean[[1]]
    nodes <- net_clean[[2]]
    downstream_network <- net_clean[[3]]
    lists <- net_clean[[6]]
    netstat <- net_clean[[8]]
    lnodes <- net_clean[[10]]

    plotTitle <- paste(lists, collapse = " X ")
    plotSubTitle <- ""
    
    plotGraph <- graph_from_data_frame(downstream_network, directed = F)%>%simplify(remove.multiple=T)
    if (input$graylinks == F) {
      nodes <- dplyr::filter(nodes, fun!="downstream")
      plotGraph <- graph_from_data_frame(pruned_network, directed = F)%>%simplify(remove.multiple=T)
    }
    
    if (!is.null(input$shcomp)) {
      if(input$shcomp != "ALL") {
        plotSubTitle <- paste("#", input$shcomp, "Component")
        netstat <- dplyr::filter(netstat, Component == as.numeric(input$shcomp))
        nodes <- dplyr::filter(nodes, nodes$id %in% netstat$Gene_Symbol)
        downstream_network <- downstream_network %>% dplyr::filter(downstream_network$from %in% nodes$id)
        pruned_network <- pruned_network %>% dplyr::filter(pruned_network$from %in% nodes$id)
      }
    }
    
    switch(input$centr,
           "deg" = size_vert <- size_vert <- degree(plotGraph, v = V(plotGraph), mode = "all"),
           "clo" = size_vert <- closeness(plotGraph, v = V(plotGraph), mode = "all", normalized = F),
           "bet" = size_vert <- betweenness(plotGraph, v = V(plotGraph), directed = F, normalized = F),
           "eig" = size_vert <- eigen_centrality(plotGraph, directed = F, scale = T)$vector,
           "pr" = size_vert <- page_rank(plotGraph, vids = V(plotGraph), directed = F)$vector
    )
    size_vert <- data.frame(id = names(size_vert), size = unname(size_vert))
    nodes$size <- size_vert$size[match(nodes$id, table = size_vert$id)]
    
    # reformating size for better visualization, and other vars
    nodes$size <- round(5*log2(nodes$size+1))
    
    # if not bipartite... default layout
    if (input$nettype != 'crt') {
      lay <- input$grlay
    } else {
      lay <- input$grlay2
    }

    if (input$graylinks == T & lay != "layout_as_bipartite") {
      subnetwork1 <- visNetwork(nodes, downstream_network, main = plotTitle, submain = plotSubTitle, width = "100%") %>% #, width = "100%", height = "100%") %>%
        visIgraphLayout(layout = lay, randomSeed = 1234) %>%
        visNodes(shadow = F) %>%
        visPhysics(stabilization = F)
      subnetwork <- subnetwork1 %>%
        visLegend(addNodes = lnodes, useGroups = F, width = 0.2, main = "Legend", position = "right") %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = T, autoResize = T) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}")
      
      # if bipartite, a bit different
    } else if (input$graylinks == F & lay != "layout_as_bipartite") {
      subnetwork1 <- visNetwork(nodes, pruned_network, main = plotTitle, submain = plotSubTitle, width = "100%") %>% #, width = "100%", height = "100%") %>%
        visIgraphLayout(layout = lay, randomSeed = 1234) %>%
        visNodes(shadow = F) %>%
        visPhysics(stabilization = F)
      subnetwork <- subnetwork1 %>%
        visLegend(addNodes = lnodes, useGroups = F, width = 0.2, main = "Legend", position = "right") %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = T, autoResize = T) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}")
    } else if (lay == "layout_as_bipartite") {
      # conditionalpanel reformating
      if(input$direct == T) {dir = "UD"} else {dir = "LR"}
      # reformating stuff to show the bipartite net
      if (input$bipnet == "list") {
        nodes$level <- as.numeric(as.factor(nodes$list))
        nodes$level[nodes$fun=="downstream"] <- length(lists)+2
        downstream_nodes1 <- unique(downstream_network$to[downstream_network$to%in%nodes$id[nodes$fun=="downstream"] & downstream_network$interaction_locationA=="Receptor"])
        nodes$level[nodes$id%in%downstream_nodes1] <- length(lists)+1
        nodes$x <- as.numeric(as.factor(nodes$color))
      } else {
        nodes$level <- as.numeric(factor(nodes$fun, levels = c("Ligand", "Receptor", "downstream")))
        nodes$level[nodes$level==3] <- 4                          
        downstream_nodes <- unique(downstream_network$to[downstream_network$to%in%nodes$id[nodes$level==4] & downstream_network$interaction_locationA=="Receptor"])
        nodes$level[nodes$id%in%downstream_nodes] <- 3
        nodes$x <- as.numeric(factor(nodes$shape))
      }
      nodes <- nodes[order(nodes$level, nodes$x), ]
      subnetwork1 <- visNetwork(nodes, pruned_network, main = plotTitle, submain = plotSubTitle, width = "100%")
      
      if (input$graylinks == T) {
        subnetwork1 <- visNetwork(nodes, downstream_network, main = plotTitle, submain = plotSubTitle, width = "100%")
      }
      subnetwork1 <- subnetwork1 %>% #, width = "100%", height = "100%") %>%
        visHierarchicalLayout(direction = dir, levelSeparation = 200, treeSpacing = 1e-5, nodeSpacing = 50, blockShifting = T, edgeMinimization = T) %>%
        visEdges(scaling = list(min=0.2, max=0.5)) %>%
        visNodes(shadow = F)
      subnetwork <- subnetwork1 %>%
        visOptions(highlightNearest = T, nodesIdSelection = T, autoResize = T) %>%
        visLegend(useGroups = F, addNodes = lnodes, position = "right", zoom = T, main = "Legend", width = 0.1) %>%
        visPhysics(enabled = F) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes[0]);
                  ;}")
      subnetwork

    }
  })
  
  output$oraplot <- renderPlot({
    
    req(!is.null(ora))
      
      ora <- ora()
      enr <- ora[[1]]
      numComp <- ora[[2]]
      netstat <- net_clean()[[8]]
      tis <- unlist(ora[[3]])

      if (identical(tis, sort(unique(netstat$List)))) {
        tis <- ""
      } else {
        tis <- paste(unlist(ora[[3]]), collapse = " + ")
      }
      
      if(is.null(enr)) { 
        output$warning4 <- renderUI({ HTML(paste(warnEnr)) })
      } else {
        enr@result <- enr@result[enr@result$p.adjust <0.05, ]
        
        if (numComp != "" | numComp != "ALL") numComp <- paste("components: #", numComp)
        
        if (input$showCat == "all") {
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
            enrichplot::dotplot(enr, showCategory = nrow(enr)) +
              ggtitle(paste0(numComp, "\n", tis)) +
              theme(plot.title = element_text(hjust = 0.5))
          }
          outora$reacplot <- reacplot()
          reacplot()
        } else {
          output$warning4 <- renderText({ "" })
          reacplot <- function(){
            enrichplot::dotplot(enr, showCategory = numCat) +
              ggtitle(paste0(numComp, "\n", tis)) +
              theme(plot.title = element_text(hjust = 0.5))
          }
          outora$reacplot <- reacplot()
          reacplot()
        }
      }
  })
  
  output$oea <- renderDataTable({

    req(!is.null(ora))
      
      enr <- ora()[[1]]
      
      req(enr)
        oratable <- as.data.frame(enr@result) %>% dplyr::select(!qvalue) %>% filter(p.adjust<0.05)
        if (nrow(oratable)>0) {
          oratable[, c("pvalue", "p.adjust")] <- apply(oratable[,c("pvalue", "p.adjust")], 2, function(x) formatC(x, format = "e", digits = 2))
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
    filename = function() {
      tagFiles <- net_clean()[[11]]
      if (input$globEnr == F) tagFiles <- gsub("_Experimental_.*", "", gsub("_HighScore_.*", "", tagFiles))
      paste0("oraplot-", tagFiles, ".png")
      },
    content = function(file) {
      reacplot <- outora$reacplot
      png(file, width = 1200, height = 700)
      print(reacplot +
              ggtitle("\nEnrichment Plot") +
              theme(plot.title = element_text(hjust = 0.5)))
      dev.off()
    })
  
  output$downloadNetTable2 <- downloadHandler(
    filename = function() {
      tagFiles <- net_clean()[[11]]
      if (input$globEnr == F) tagFiles <- gsub("_Experimental_.*", "", gsub("_HighScore_.*", "", tagFiles))
      paste0("oratable-", tagFiles, ".txt")
    }, content = function(file) {
      oratable <- outora$oratable
      oratable[,9] <- as.character(oratable[,9])
      oratable[,9] <- substr(oratable[,9], 10, nchar(oratable[,9])-55)
      write.table(oratable, file, row.names = FALSE, quote = F, sep = "\t")
    })
  
  output$downloadNetTable <- downloadHandler(
    filename = function() {
      paste0('netlist-', net_clean()[[11]], ".txt")
      }, content = function(file) {
      write.table(net_clean()[[9]], file, row.names = FALSE, quote = F, sep = "\t")}
    )
  
  output$downloadTopoTable <- downloadHandler(
    filename = function() {
      paste0('nodelist-', net_clean()[[11]], ".txt")
      }, content = function(file) {
      write.table(net_clean()[[8]], file, row.names = FALSE, quote = F, sep = "\t")})
  
  output$downloadDrugTable <- downloadHandler(
    filename = function() {
      paste0('drugRec_ints-', net_clean()[[11]], ".txt")
      }, content = function(file) {
      recs <- drint$recs
      write.table(recs, file, row.names = FALSE, quote = F, sep = "\t")
    })
  
  vals <- reactiveValues(coords=NULL)
  observe( {
      invalidateLater(1000)
      net_clean <- net_clean()
      req(length(net_clean)==11)
      nodes <- net_clean[[2]]
      visNetworkProxy("network1") %>% visGetPositions()
      vals$coords <- if (!is.null(input$network1_positions))
        do.call(rbind.data.frame, input$network1_positions)
      return(vals)
  })
  
  output$downloadNetCSV <- downloadHandler(
    filename = function() {
      tagFiles <- net_clean()[[11]]
      if(!is.null(input$shcomp) && input$shcomp != "ALL") tagFiles <- paste0(tagFiles, "_comp", input$shcomp)
      paste0('network-', tagFiles, '.csv')
    }, content = function(file) { 
      nodes <- net_clean()[[2]]
      if (input$graylinks == F) {
        nodes <- nodes %>% filter(color != "lightgray")
      }
      write.table(nodes, file, row.names = F, quote = F, sep = ",")
  })
  
  output$downloadNetGML <- downloadHandler(
    filename = function() {
      tagFiles <- net_clean()[[11]]
      if(!is.null(input$shcomp) && input$shcomp != "ALL") tagFiles <- paste0(tagFiles, "_comp", input$shcomp)
      paste0('network-', tagFiles, '.GraphML')
    }, 
    content = function(file) {
      net_clean <- net_clean()
      write_graph(graph_from_data_frame(d = net_clean[[9]], vertices = net_clean[[2]]), file, format = "graphml")}
  )
  
  output$downloadNetHTML <- downloadHandler(
    filename = function() {
      tagFiles <- net_clean()[[11]]
      if(!is.null(input$shcomp) && input$shcomp != "ALL") tagFiles <- paste0(tagFiles, "_comp", input$shcomp)
      paste0('network-', tagFiles, ".html")
    }, 
    content = function(con) {
      edgesSmall <- net_clean()[[1]]
      nodes <- net_clean()[[2]]
      edges <- net_clean()[[3]]
      lists <- net_clean()[[6]]
      netstat <- net_clean()[[8]]
      lnodes <- net_clean()[[10]]
      
      if (input$graylinks == F) {
        nodes <- nodes[nodes$color != 'lightgray' , ]
        edges <- edgesSmall
      }
      
      if (!is.null(input$shcomp)) {
        if(input$shcomp != "ALL") {
          netstat <- dplyr::filter(netstat, Component == as.numeric(input$shcomp))
          nodes <- dplyr::filter(nodes, nodes$id %in% netstat$Gene_Symbol)
          edges <- edges %>% dplyr::filter(edges$from %in% nodes$id)
        }
      }
      
      plotGraph <- graph_from_data_frame(edges, directed = F)%>%simplify(remove.multiple=T)
      
      switch(input$centr,
             "deg" = size_vert <- size_vert <- degree(plotGraph, v = V(plotGraph), mode = "all"),
             "clo" = size_vert <- closeness(plotGraph, v = V(plotGraph), mode = "all", normalized = F),
             "bet" = size_vert <- betweenness(plotGraph, v = V(plotGraph), directed = F, normalized = F),
             "eig" = size_vert <- eigen_centrality(plotGraph, directed = F, scale = T)$vector,
             "pr" = size_vert <- page_rank(plotGraph, vids = V(plotGraph), directed = F)$vector
      )
      size_vert <- data.frame(id = names(size_vert), size = unname(size_vert))
      nodes$size <- size_vert$size[match(nodes$id, table = size_vert$id)]
      
      # reformating size for better visualization, and other vars
      nodes$size <- round(5*log2(nodes$size+1))
    
      vals$coords$id <- row.names(vals$coords)
      vals$coords <- left_join(nodes, vals$coords, by = "id")
      height <- max(10*nrow(nodes), 1000)

      visNetwork(vals$coords, edges, width = "100%", height = height, main = paste(lists, collapse = " X ")) %>% 
        visOptions(highlightNearest = T) %>% visLegend(useGroups = F, addNodes = lnodes, width = 0.2, main = "Legend", position = "right") %>%
        visExport() %>% 
        visPhysics(enabled = F) %>% 
        visSave(con)
    }
  )
  
  observe( {
    
    net_clean <- net_clean()
    
    req(length(net_clean)==11)
                  
      nodes <- net_clean[[2]]
      downstream_network <- net_clean[[3]]
      netstat <- net_clean[[8]]
      lnodes <- net_clean[[10]]
      nodes_sel_entrez <- NULL
      nodes_in <- NULL
      tagFiles <- net_clean[[11]]
      
      if(length(input$click) == 1) {
        shinyjs::show("downstream")
        shinyjs::show("downstreamID")
        shinyjs::show("downloadPaths")
        
        if (!is.null(input$shcomp)) {
          if (input$shcomp != "ALL") {
            netstat <- dplyr::filter(netstat, Component == as.numeric(input$shcomp))
            nodes <- dplyr::filter(nodes, nodes$id %in% netstat$Gene_Symbol)
            downstream_network <- downstream_network %>% dplyr::filter(downstream_network$from %in% nodes$id)
          }
        }
        
        visNetworkProxy("network1") %>%
          visGetSelectedNodes(input = "connected_nodes")
        filt <- grepl(paste("^", input$connected_nodes, "$", sep =""), downstream_network$from) & downstream_network$interaction_locationB!="downstream"
        downstream_network_tmp <- downstream_network[!filt, ]
        down_graph <- graph_from_data_frame(downstream_network_tmp, directed = T)
        if (!is.null(input$connected_nodes) & input$graylinks == T) {
          nodes$color[nodes$color == 'gray25'] <- "lightgrey"
          if (input$connected_nodes%in%V(down_graph)$name) {
            down_sel <- unname(bfs(down_graph, root = input$connected_nodes, neimode = "out", unreachable = F)$order)
            nodes_sel <- V(down_graph)[down_sel[!is.na(down_sel)][-1]]$name
            nodes_sel_entrez <- idResource$entrez.id[match(nodes_sel, idResource$gene.symbol)]
            nodes_sel_entrez <- nodes_sel_entrez[!is.na(nodes_sel_entrez)]
            sel_in <- length(nodes_sel_entrez)
            nodes <- changeColorOfNodes(nodes, nodes_sel)
            lnodes <- rbind.data.frame(lnodes, c("Downstream\nselected", "square", "gray25"))
          } else {
            nodes$color[nodes$color == 'gray25'] <- "lightgrey"
            lnodes <- lnodes[!grepl("gray25", lnodes$color), ]
          }
          visNetworkProxy("network1") %>%
            visUpdateNodes(nodes)
        } else if (input$graylinks == F) {
          visNetworkProxy("network1") %>%
            visUpdateNodes(nodes %>% filter(fun != "downstream")) %>%
            visUpdateEdges(downstream_network %>% filter(interaction_locationA == "Ligand") %>% filter(interaction_locationB == "Receptor"))
        } else if (is.null(input$connected_nodes)) {
          shinyjs::hide("downstream")
          shinyjs::hide("downstreamID")
          shinyjs::hide("downloadPaths")
        }
        ent_sel <- idResource$entrez.id[grep(paste("^", input$connected_nodes, "$", sep = ""), idResource$gene.symbol)]
        sym_sel <- idResource$gene.symbol[grep(paste("^", input$connected_nodes, "$", sep = ""), idResource$gene.symbol)]
        hsa <- ezID2path[grep(paste("^", as.character(ent_sel), "$", sep = ""), names(ezID2path))]
        hsa <- names(path2name[names(path2name)%in%unname(unlist(hsa))])
        ids <- gsub("Homo sapiens: ", "", unlist(unname(path2name[names(path2name)%in%unname(unlist(hsa))])))
        reaclink <- createLink(hsa)
        if (!is.null(input$connected_nodes)) { 
          if(input$graylinks == T) {
            if (input$connected_nodes%in%V(down_graph)$name) {
              nodes_in <- suppressWarnings(do.call(rbind, lapply(hsa, FUN = function(x) nodes_sel_entrez%in%as.integer(reactome.db::reactomePATHID2EXTID[[x]]))))
              if (length(nodes_in)!=0) {
                path_len <- do.call(rbind, lapply(hsa, function(x)length(reactome.db::reactomePATHID2EXTID[[x]])))
                univ_len <- nrow(idResource)-path_len
                nodes_in <- apply(nodes_in, 1, function(x) nodes_sel_entrez[x])
                path_in <- do.call(rbind, lapply(nodes_in, length))
                nodes_in <- lapply(nodes_in, function(x) idResource$gene.symbol[match(x, idResource$entrez.id)])
                nodes_in <- unlist(lapply(nodes_in, function(x) paste(x, collapse = " / ")))
                contTable <- data.frame(univ_len, path_len, sel_in, path_in)
                pvals <- apply(contTable, 1, function(x) unlist(fisher.test(matrix(x, nrow = 2))[1]))
                p.adjust <- p.adjust(pvals, method = "BH", n = length(pvals))
                output$downstreamID <- renderUI({ HTML(paste('<br/>', '<br/>', "<center>", "<b>", sym_sel, "</b>", "</center>"))})
                output$downstream <- renderDataTable({
                  downpaths <- cbind.data.frame(hsa, ids, path_len, pvals, p.adjust, nodes_in, reaclink)
                  downpaths <- downpaths %>% arrange(pvals) %>% filter(p.adjust < 0.05)
                  if (nrow(downpaths)==0) {
                    output$warning1_2 <- renderUI({ HTML(paste('<br/>', '<br/>', "No enriched pathways found")) })
                  } else {
                    output$warning1_2 <- renderText({ "" })
                    downpaths[, c("pvals", "p.adjust")] <- apply(downpaths[, c("pvals", "p.adjust")], 2, function(x) formatC(x, digits = 2, format = "E"))
                  }
                  colnames(downpaths) <- c('ReactomeID', "Pathway_name", "pathway_length", "p.value", "p.adjust", "nodes_in", "link")
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
          filename = function() {
            if (ncol(downpaths) == 3) {
              tagFiles <- paste0(tagFiles, "_involved")
            } else {
              tagFiles <- paste0(tagFiles, "_signif")
            }
            paste0(sym_sel, tagFiles, 'Pathways.txt')
            }, content = function(file) {
            downpaths[, ncol(downpaths)] <- substr(downpaths[, ncol(downpaths)], 10, nchar(downpaths[,ncol(downpaths)])-55)
            write.table(downpaths, file, row.names = FALSE, quote = F, sep = "\t")
          })

      } else {
      }
  })
}

shinyApp(ui, server)
