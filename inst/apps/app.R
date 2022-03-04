############################################################
# 
# author: Ludwig Geistlinger
# date: 2021-08-02 08:22:40
# 
# descr: scaffold for simple graph viewer
# 
############################################################

library(shiny)
library(BioPlex)
library(BioNet)
library(ggnetwork)

gr <- readRDS("ribo_graph.rds")
gr <- BioNet::largestComp(gr)
ig <- igraph::graph_from_graphnel(gr)

# quick fix
igraph::vertex_attr(ig, "ENTREZID") <- as.character(igraph::vertex_attr(ig, "ENTREZID"))

.getNumericVertexAttributes <- function(ig)
{
    vas <- igraph::vertex_attr(ig)
    ind <- vapply(vas, is.numeric, logical(1))
    return(names(vas)[ind])
}

.getNumericEdgeAttributes <- function(ig)
{
    vas <- igraph::edge_attr(ig)
    ind <- vapply(vas, is.numeric, logical(1))
    return(names(vas)[ind])
}

buildUI <- function(ig)
{
    nchoices <- .getNumericVertexAttributes(ig)
    echoices <- .getNumericEdgeAttributes(ig)
    
    ui <- fluidPage(
        
        # Application title
        titlePanel("BioPlex graph data viewer"),
        
        # Sidebar with a slider input for number of bins 
        sidebarLayout(
            sidebarPanel(
                selectInput("ndata",
                            label = "Node data:",
                            choices = nchoices,
                            selected = nchoices[1]),
                selectInput("edata",
                            label = "Edge data:",
                            choices = echoices,
                            selected = echoices[1])
            ),
            
            # TODO: support hover 
            mainPanel(
                div(
                    style = "position:relative",
                    plotOutput("ggnplot", hover = hoverOpts(id = "plot_hover",
                                                            delay = 100, 
                                                            delayType = "debounce")),
                    uiOutput("hover_info")
                ),
            )
        )
    )
    return(ui)
}

# Define server logic required to draw a histogram
server <- function(input, output) 
{
    df <- ggnetwork(ig)
    output$ggnplot <- renderPlot({
        ndata <- input$ndata
        edata <- input$edata
        ggplot(df, 
               aes(x = x, y = y, xend = xend, yend = yend)) +
            geom_nodes(aes_string(size = ndata)) +
            geom_nodelabel_repel(aes(label = SYMBOL), box.padding = unit(1, "lines")) +
            geom_edges(aes_string(color = edata), 
                       arrow = arrow(length = unit(6, "pt"), type = "closed")) +
            theme_blank()
    })
    
    output$hover_info <- renderUI({
        hover <- input$plot_hover
        point <- nearPoints(df, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
        if (nrow(point) == 0) return(NULL)
        
        # calculate point position INSIDE the image as percent of total dimensions
        # from left (horizontal) and from top (vertical)
        left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
        top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
        
        # calculate distance from left and bottom side of the picture in pixels
        left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
        top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
        
        # create style property fot tooltip
        # background color is set so tooltip is a bit transparent
        # z-index is set so we are sure are tooltip will be on top
        style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                        "left:", 10 , "px; top:", 300, "px;")
        
        # actual tooltip created as wellPanel
        ndata <- input$ndata
        nstr <- paste0(ndata, ":")
        nstr <- paste("<b>", nstr, "</b>")
        wellPanel(
            style = style,
            p(HTML(paste0("<b> Isoform: </b>", point$ISOFORM, "<br/>",
                          "<b> Symbol: </b>", point$SYMBOL, "<br/>",
                          nstr, signif(point[[ndata]], digits = 1))))
        )
    })
}


GraphViewer <- function(ig)
{
    shinyApp(ui = buildUI(ig), server = server)
}

# run the app
GraphViewer(ig)
