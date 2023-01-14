#' Plot a graph 
#' 
#' @description Plotting utility for graphs
#' @param gr graph. An object of class \code{graphNEL} from the graph package.
#' @param node.label character. A node data attribute that is used to label the 
#' nodes.
#' @param node.data character. A node data attribute that is used to overlay on
#' the nodes using different colors.
#' @param edge.data character. An edge data attribute that is used to overlay
#' on the edges using different line types.
#' @param node.color character. Default color for the nodes. Is overwritten by
#' \code{node.data} if not null.
#' @param edge.color character. Default color for the edges.
#' @param edge.mode character. Whether to display directed or undirected edges.
#' @param layout character. Circular layout or network layout.
#' @return A ggplot. 
#' @examples
#'      library(BioPlex)
#'      corum.df <- getCorum("core")
#'      corum.glist <- corum2graphlist(corum.df)
#'      gr <- corum.glist[["CORUM107_TFIIH_transcription_factor_complex"]]
#'
#'      plotGraph(gr, edge.mode = "undirected", 
#'                edge.color = "grey", layout = "circle")
#'
#' @export
plotGraph <- function(gr,
                      node.label = "SYMBOL",
                      node.data = NULL,
                      edge.data = NULL,   
                      node.color = "firebrick",
                      edge.color = "firebrick",
                      edge.mode = c("directed", "undirected"),
                      layout = c("fruchterman.reingold", "circle"))
{
    edge.mode <- match.arg(edge.mode)
    layout <- match.arg(layout)

    if(!is.null(node.data))
    {
        has.node.data <- node.data %in% names(graph::nodeDataDefaults(gr))
        stopifnot(is.character(node.data) && has.node.data)
    }

    if(!is.null(edge.data))
    {
        has.edge.data <- edge.data %in% names(graph::edgeDataDefaults(gr))
        stopifnot(is.character(edge.data) && has.edge.data)
    }

    ig <- igraph::graph_from_graphnel(gr)
    igraph::vertex_attr(ig, "ENTREZID") <- as.character(igraph::vertex_attr(ig, "ENTREZID"))

    # layout
    layout <- switch(layout,
                     fruchterman.reingold = igraph::nicely(),
                     circle = igraph::in_circle())        
    df <- ggnetwork::ggnetwork(ig, layout = layout)
    
    
    # plot
    gp <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend))

    # node layout
    if(is.null(node.data))
        gp <- gp + ggnetwork::geom_nodes(color = node.color, size = 3)
    else
        gp <- gp + ggnetwork::geom_nodes(ggplot2::aes_string(color = node.data),
                                         size = 3)

    # edge layout
    if(is.null(edge.data))
    {
        if(edge.mode  == "undirected")
            gp <- gp + ggnetwork::geom_edges(color = edge.color)
        else
            gp <- gp + ggnetwork::geom_edges(color = edge.color,
                                  arrow = grid::arrow(length = grid::unit(6, "pt"),
                                                type = "closed"))
    }
    else
    {
        if(edge.mode  == "undirected")
            gp <- gp + ggnetwork::geom_edges(ggplot2::aes_string(
                                                linetype = edge.data),
                                                color = edge.color)
        else
            gp <- gp + ggnetwork::geom_edges(
                                  ggplot2::aes_string(linetype = edge.data),
                                  color = edge.color,
                                  arrow = grid::arrow(
                                            length = grid::unit(6, "pt"),
                                            type = "closed"))
    }
    
    gp <- gp + ggnetwork::geom_nodetext_repel(
                    ggplot2::aes_string(label = node.label), size = 3)
    gp + ggnetwork::theme_blank()
}

