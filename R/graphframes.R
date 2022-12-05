#' @title Representation of BioPlex PPIs using a GraphFrames backend
#' @description Representation of BioPlex PPIs in a \code{GraphFrame} object
#' from the \code{graphframes} package.
#' @param gr an object of class \code{\linkS4class{graph}} storing the
#' BioPlex PPIs. Typically obtained via \code{bioplex2graph}.
#' @param spark.con Spark connection. Typically obtained via
#' \code{sparklyr::spark_connect}.
#' @return An object of class \code{GraphFrame}. 
#' @seealso \code{\link{graphframe2graph}}
#' @examples
#'
#' library(sparklyr)
#' library(graphframes)
#' library(BioPlex)
#'
#' # (1) Obtain the latest version of the 293T PPI network
#' bp.293t <- getBioPlex(cell.line = "293T", version = "3.0")
#' 
#' # (2) Turn the data into a graph 
#' bp.gr <- bioplex2graph(bp.293t)
#' 
#' # (3) Switch to a graphframes backend
#' sc <- spark_connect(master = "local", version = "3.0")
#' bp.gf <- graph2graphframe(bp.gr, sc) 
#'
#' @importFrom utils stack
#' @export
graph2graphframe <- function(gr, spark.con) 
{
    # node df
    nattrs <- names(graph::nodeDataDefaults(gr))
    .getNCol <- function(n) graph::nodeData(gr, graph::nodes(gr), n)
    cols <- lapply(nattrs, .getNCol)
    lens <- vapply(cols, function(x) all(lengths(x) == 1), logical(1))
    .collna <- function(y) ifelse(length(y) == 1 && is.na(y),
                                  NA_character_,
                                  paste(y, collapse = ","))
    .collapse <- function(x) vapply(x, .collna, character(1))
    cols[!lens] <- lapply(cols[!lens], .collapse)
    cols <- lapply(cols, unlist)
    names(cols) <- nattrs
    node.df <- data.frame(id = graph::nodes(gr), cols)
    colnames(node.df) <- tolower(colnames(node.df))
    node.df <- sparklyr::copy_to(spark.con, node.df)
    
    # edge df
    edge.df <- stack(graph::edges(gr))
    colnames(edge.df) <- c("dst", "src")
    edge.df <- edge.df[,2:1]
    for(i in 1:2) edge.df[,i] <- as.vector(edge.df[,i]) 
    
    eattrs <- names(graph::edgeDataDefaults(gr))
    eattrs <- setdiff(eattrs, "weight")
    .getECol <- function(n) graph::edgeData(gr, 
                                            from = edge.df$src, 
                                            to = edge.df$dst, n)
    cols <- lapply(eattrs, .getECol)
    cols <- lapply(cols, unlist)
    names(cols) <- eattrs
    edge.df <- data.frame(edge.df, cols)

    edge.df <- sparklyr::copy_to(spark.con, edge.df)
    bp.gf <- graphframes::gf_graphframe(node.df, edge.df)
    return(bp.gf)
}

#' @title Convert a GraphFrames object to a graphNEL object
#' @description Conversion of a GraphFrames object from the graphframes package
#' to a graphNEL object from the graph package.
#' @param gf An object of class \code{GraphFrame}. 
#' @return An object of class \code{\linkS4class{graph}}
#' @seealso \code{\link{graph2graphframe}}
#' @export
graphframe2graph <- function(gf)
{
    node.df <- data.frame(graphframes::gf_vertices(gf))
    edge.df <- data.frame(graphframes::gf_edges(gf))

    ft.cols <- c("src", "dst")
    ftm <- as.matrix(edge.df[,ft.cols])
    gr <- graph::ftM2graphNEL(ftm, edgemode = "directed")

    # node data
    ncols <- setdiff(colnames(node.df), "id")
    for(col in ncols)
    {
        graph::nodeDataDefaults(gr, col) <- NA
        graph::nodeData(gr, node.df$id, col) <- node.df[,col]
    }

    # edge data
    ecols <- setdiff(colnames(edge.df), ft.cols)
    for(col in ecols)
    {   
        graph::edgeDataDefaults(gr, col) <- numeric(0L)
        graph::edgeData(gr, edge.df[,"src"], edge.df[,"dst"], col) <- edge.df[,col]
    }   

    return(gr)
}
