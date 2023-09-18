library(networkD3)
library(visNetwork)
library(pheatmap)
library(igraph,warn.conflicts = FALSE)
library(shiny)
library(DT)

###################Preprocessing#######################

#Here, we should load all the files
load("data/dataForNetwork2.RData")
load("data/functionalData.RData")
load("data/geneAnnot.RData")
ZTLookup=as.numeric(substring(colnames(diurnalExpression),1,2))
cellTypesInOrder=c(
  "epidermal", 
  "stressed epidermal",
  "mesophyll group 1", 
  "mesophyll group 2",
  "mesophyll group 3",
  "stressed mesophyll", 
  "bundle sheath",
  "phloem",
  "vascular",
  "companion",
  "guard",
  "sieve",
  "hydathode",
  "myrosinase",
  "unknown group 1",
  "unknown group 2",
  "unknown group 3")

currentDAPseq=dapSeq

####Set up table merge helper functions
mergeAnnot <-function(a) {
  a=cbind(a, rep("", length(a[,1])), rep("", length(a[,1])), rep("", length(a[,1])), rep("", length(a[,1])))
  colnames(a)=c("Source", "Target", "Name Source", "Name Target", "Description Source", "Description Target")
  idsS=which(a[,1] %in% rownames(geneAnnot))
  idsT=which(a[,2] %in% rownames(geneAnnot))
  a[idsS,c(3, 5)]=geneAnnot[a[idsS,1],]
  a[idsT,c(4, 6)]=geneAnnot[a[idsT,2],]
  a=rbind(colnames(a), a)
  a
}


####Set up pafway helper functions:
pafway <- function(GO, edges, GOtypes) {
  
  GOinNetwork = GO[unique(c(edges[, 1], edges[, 2]))]
  
  grepLen=sapply(GOtypes, function(i){
    length(grep(i, GOinNetwork))
  })
  names(grepLen)=GOtypes
 
  sapply(GOtypes, function(i) {
    if(grepLen[i]!=0){
      #find edges first:
      grepI=grep(i, GOinNetwork[edges[, 1]])
      
      sapply(GOtypes, function(j) {
        if(grepLen[j]!=0){
         
          a=length(grep(j, GOinNetwork[edges[grepI, 2]]))
          
          p_bot = grepLen[i]/length(GOinNetwork) * grepLen[j]/length(GOinNetwork)
         
          b = stats::binom.test(a, length(edges[, 1]), p = p_bot, alternative = c("greater"))
          b$p.value
        }else{1}
      })
    }else{rep(1, length(GOtypes))}
  })
}

transformToEdgeList<- function(mat, thresh=0.001){
  ind=which(mat<thresh, arr.ind = T)
  cbind(colnames(mat)[ind[,2]], rownames(mat)[ind[,1]], signif(diag(mat[ind[,1], ind[,2]]), digits=3))
}




###The UI

ui <- fluidPage(

  headerPanel( title=div(img(src="./images/logo.svg", width = 60), "AraLeTA")),
 # tags$img(src = "./images/logo.svg", width = "99px"),
 h2("Arabidopsis Leaf Time-dependent Atlas"),
  sidebarPanel(
    
    
    # Input: Slider for the number of bins ----
    sliderInput(inputId = "age",
                label = "Days post-germination",
                min = 4,
                max = 30,
                step = 2,
                value = c(4, 30)),
    p("\nNote that the G-to-M transition is around day 13 and the M-to-S transition is around day 21."),
    hr(style = "border-top: 1px solid #000000;"),
    sliderInput(inputId = "time",
                label = "Time of day (hours from dawn)",
                min = 0,
                max = 24,
                value = c(0, 24)),
    p("\nNote that the lights were turned off after the 10h sample."),
    hr(style = "border-top: 1px solid #000000;"),
    
    checkboxGroupInput(inputId = "cell", label="Cell type", selected = cellTypesInOrder,
      choices = cellTypesInOrder
    ),
    hr(style = "border-top: 1px solid #000000;"),
    sliderInput(inputId = "bulkThreshold",
                label = "Bulk RNA-seq filtering threshold",
                min = 0,
                max = 50,
                value = 5),
    sliderInput(inputId = "scThreshold",
                label = "scRNA-seq filtering threshold",
                min = 0.0,
                max = 2,
                step=0.01,
                value = 0.05),
    
    width=2
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Table View",
               DTOutput("outdata"),
               downloadButton("downloadPathcsv", "Download Network")
      ),
      tabPanel("Heatmap Network View",
               h3("Heatmap view of network"),
               p("Rows are sources and columns are targets.  Red implies an edge is present."),
               p("Please note that a random subsample of 2000 target genes are selected if network is too large."),
               plotOutput("simpleHeatmap")),
      tabPanel("Heatmap Cell Annotations",
               h3("Heatmap view of network"),
               p("Rows are sources and columns are targets.  Red implies an edge is present."),
               p("Please note that a random subsample of 2000 target genes are selected if network is too large."),
               plotOutput("heatmapAnnotCellTypes")),
       tabPanel("Plot IGraph",
                plotOutput("igraphplot")
       ),
       tabPanel("PAFway Table",
              h3("Associations between functional terms within network"),
              p("Please note that the PAFway network will only display if there are <100,000 edges.  Otherwise, the full network will be displayed.  Also note that PAFway tables can take up to 5 minutes to load."),
              DTOutput("pafwayNetworkData"),
              downloadButton("downloadPafwaycsv", "Download Network")
       ),
      tabPanel("PAFway Network",
              h3("Top 100 edges"),
              plotOutput("igraphplotpafway")
     ),
     tabPanel("Citations",
              h3("Papers to cite when using AraLeTA"),
              p("AraLeTA: Vong et al. (2023) BioRxiv.  DOI:"),
              p("PAFway: Mahjoub and Ezer (2020) Bioinformatics. DOI: 10.1093/bioinformatics/btaa639"),
              p("Developmental time series: Woo et al. (2016). DOI: 10.1104/pp.15.01929"),
              p("Diurnal time series: Hickman et al. (2017). DOI: 10.1105/tpc.16.00958"),
              p("Leaf scRNA-seq: Procko et al. (2022). DOI: 10.1093/plcell/koac167"),
              h3("More info"),
              tags$a(href="https://stressedplants.github.io/labwebsite/", "Learn more about our lab!"),
        
              
     ),

    ),
    width = 10
  ),
 # absolutePanel(
#    img(src='images/logo.svg', align = "right-top", width="20")
 #   , id = "input_date_control", class = "panel panel-default", fixed = TRUE
#    , draggable = FALSE, top = 'auto', left = 'auto', right = 0, bottom = 0
 #   , width = '5%', height = 'auto')
  
)

server <- function(input, output,session) {


 # cell<- reactive({
 #   input$cell
 # }) %>% debounce(1000)
  
  filter_network <- reactive({
    
    if(input$bulkThreshold==0){
      edgesToInclude1=c(1:dim(dapSeq)[1])
      edgesToInclude2=c(1:dim(dapSeq)[1])}else{
        #filter by developmental stage
       
        devRange=paste("rep1_X", seq(input$age[1], input$age[2]+1, 2), "D", sep="")
        
        if(length(devRange)>1){
        maxExpInRange=apply(developmentalExpression[,devRange], 1, max)
        
        genesInRange=names(which(maxExpInRange>input$bulkThreshold))
        
        edgesToInclude1=which( (dapSeq[,1] %in% genesInRange) & (dapSeq[,2] %in% genesInRange) )
        }else{edgesToInclude1=c()}
       
        #filter by ZT
        ztRange=which(ZTLookup>=input$time[1] & ZTLookup<=input$time[2])
        if(length(ztRange)>1){
        maxExpInRange=apply(diurnalExpression[,ztRange], 1, max)
        genesInRange=names(which(maxExpInRange>input$bulkThreshold))
        edgesToInclude2=which( (dapSeq[,1] %in% genesInRange) & (dapSeq[,2] %in% genesInRange) )
        }else{edgesToInclude2=c()}
      
      }
    
    #filter by cell type
    if(input$scThreshold!=0){
      
      
      
      if(length(input$cell)>1){
      
      #if(length(cell)>0){
      cellsToInclude=scRNA[,input$cell]
      maxExpInRange=apply(cellsToInclude, 1, max)
      genesInRange=names(which(maxExpInRange>input$scThreshold))
      edgesToInclude3=which( (dapSeq[,1] %in% genesInRange) & (dapSeq[,2] %in% genesInRange) )
     
      }else{
        
        if(length(input$cell)==1){
         
          maxExpInRange=scRNA[,input$cell]
       
          genesInRange=names(which(maxExpInRange>input$scThreshold))
         
          edgesToInclude3=which( (dapSeq[,1] %in% genesInRange) & (dapSeq[,2] %in% genesInRange) )
        }else{
        
        edgesToInclude3=c()}}
    }else{
      edgesToInclude3=edgesToInclude2
    }
    edgesToInclude=edgesToInclude1[which( (edgesToInclude1 %in% edgesToInclude2) & (edgesToInclude1 %in% edgesToInclude3))]
    currentDAPseq=dapSeq[edgesToInclude,]
    dapSeq[edgesToInclude,]
  }) %>% debounce(1000)
  
  
  output$outdata <- renderDT({
    a=filter_network()
    
    shiny::validate(
      need(dim(a)[1]>0, "No genes selected")
    )
    
    mergeAnnot(a)
    
    
    }, options = list(lengthChange = F))
  
  
  ###Now do the same with pafway
  
  pafway_network <- reactive({
    net=filter_network()
    if(length(net[,1])>100000){
      pafOuts=NA }else{
      pafOuts=pafway(GOconcat, net, goOfInterest)
      a=transformToEdgeList(pafOuts)
      a=a[order(as.numeric(a[,3])),]
      colnames(a)=c("Source", "Target", "P-value")
      a
    }
  })
  
  output$pafwayNetworkData <- renderDT({
    a=pafway_network()
    if(is.null(dim(a))){
      filter_network()
    }else{a}
  }, options = list(lengthChange = F))
  
  
  mainHeatmap <- reactive({
    net=filter_network()
    
    colsToUse=unique(net[,2])
    if(length(colsToUse)>2000){
      colsToUse=sample(colsToUse, 2000)
    }
    
    mat=matrix(0, nrow=length(unique(net[,1])), ncol=length(colsToUse), 
               dimnames=list(unique(net[,1]), colsToUse))
    
    netSub=net[which(net[,2] %in% colsToUse),]
    for(i in 1:length(netSub[,1])){
      mat[netSub[i,1], netSub[i,2]]=1
      }
    
    mat
    
   
  })
  
  output$simpleHeatmap <- renderPlot({
    pheatmap(mainHeatmap())
  })
  
  output$heatmapAnnotCellTypes <- renderPlot({
    temp=mainHeatmap()
    temp=temp[which(rownames(temp) %in% rownames(scRNA)),]
    
    if(length(input$cell)==1){
      annotRow=matrix(scRNA[rownames(temp), input$cell], byrow=T)
    }else{
    
    annotRow=scRNA[rownames(temp), input$cell]
    }
    
    #annotRow=scRNA[rownames(temp), cell]
    annotRow=apply(annotRow, 1, function(i){
      if((max(i)-mean(i))<0.001){0}else{
      (i-min(i))/(max(i)-min(i))}})
    
    pheatmap(temp, annotation_row = data.frame(t(annotRow)), legend=F)
  })
  
  output$downloadPathcsv <- downloadHandler(
    filename = function() {
      paste("dapseqfiltered-csv-age", input$age[1], "-", input$age[2], "_time", input$time[1], "-", input$time[2], ".csv", sep = "")
    },
    content = function(file) {
      a=filter_network()
     
      
      cat(paste("#Parameters:", 
                "AGE:", input$age[1], "to", input$age[2],
                "; TIME:", input$time[1], "to", input$time[2],
                "; CELLS:", paste(input$cell, collapse=","),
                "; SC_THRESH", input$scThreshold,
                "; BULK_THRESH", input$bulkThreshold,
                "; Columns:", collapse=","), file=file)
      write.table(mergeAnnot(a), file, row.names = FALSE, append=T, sep=",")
      
    })
  
  output$igraphplot <- renderPlot({
    
    #shiny::validate(
    #  need(input$file1, "Waiting for file!")
    #)
    
    graph=filter_network()
    graph=graph[which(graph[,2] %in% dapSeq[,1]),]
    ig <- graph_from_data_frame(d=data.frame(graph), vertices=NULL, directed=T)
    
    plot(ig, layout=layout_with_fr(ig), edge.arrow.size=.5, vertex.color="green", vertex.size=15, 
         
         vertex.frame.color="magenta", vertex.label.color="black", 
         
         vertex.label.cex=0.5, vertex.label.dist=2, edge.curved=0,zoom=TRUE) 
  })
  
  output$igraphplotpafway <- renderPlot({
    
    #shiny::validate(
    #  need(input$file1, "Waiting for file!")
    #)
    
    graph=pafway_network()
    if(is.null(dim(graph))){
      graph=filter_network()
      graph=graph[which(graph[,2] %in% dapSeq[,1]),]
      }
    
    ig <- graph_from_data_frame(d=data.frame(graph[1:min(100, length(graph[,1])),]), vertices=NULL, directed=T)
   
    plot(ig, layout=layout_with_kk(ig), edge.arrow.size=0.8, vertex.color="cornflowerblue", vertex.size=15, 
         
         vertex.frame.color="darkblue", vertex.label.color="black", 
         
         vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0,zoom=TRUE) 
  })
  
  output$pafwaySimpleNetwork <- renderPlot({
    
    graph=pafway_network()
    if(is.null(dim(graph))){
      graph=filter_network()
      graph=graph[which(graph[,2] %in% dapSeq[,1]),]
    }
    
    nodes <-  data.frame(
      id = unique(c(graph[,1], graph[,2])),
      label= unique(c(graph[,1], graph[,2]))
    )
    
    edges <- data.frame(
      from = graph[,1], 
      to = graph[,2],
      arrows = c("to"),
      label = graph[,3]
    )
    
    
    # arrows
    arrows = c("to", "from", "middle", "middle;to")
    visNetwork(nodes, edges, 
               height = "100%", width = "100%",
               main = "")
    
    })
  
  output$downloadPafwaycsv <- downloadHandler(
    filename = function() {
      paste("pafway-csv-age", input$age[1], "-", input$age[2], "_time", input$time[1], "-", input$time[2], ".csv", sep = "")
    },
    content = function(file) {
      
      cat(paste("#Parameters:", 
                "AGE:", input$age[1], "to", input$age[2],
                "; TIME:", input$time[1], "to", input$time[2],
                "; CELLS:", paste(input$cell, collapse=","),
                "; SC_THRESH", input$scThreshold,
                "; BULK_THRESH", input$bulkThreshold,
                "; Columns:", collapse=","), file=file)
      
      write.table(pafway_network(), file, row.names = FALSE, append=T, sep=",")
    })
  
}

shinyApp(ui, server)
