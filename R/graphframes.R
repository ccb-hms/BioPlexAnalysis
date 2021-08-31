#' @title Representation of BioPlex PPIs using a GraphFrames backend
#' @description Representation of BioPlex PPIs in a \code{GraphFrame} object
#' from the \code{graphframes} package.
#' @param bp.gr an object of class \code{\linkS4class{graph}} storing the
#' BioPlex PPIs. Typically obtained via \code{\link{bioplex2graph}}.
#' @return An object of class \code{GraphFrame}. 
#' @references BioPlex: \url{https://bioplex.hms.harvard.edu/interactions.php}
#' @seealso \code{\link{bioplex2graph}}
#' @examples
#'
#' # (1) Obtain the latest version of the 293T PPI network
#' bp.293t <- getBioPlex(cell.line = "293T", version = "3.0")
#' 
#' # (2) Turn the data into a graph 
#' bp.gr <- bioplex2graph(bp.293t)
#' 
#' # (3) Switch to a graphframes backend
#' sc <- sparklyr::spark_connect(master = "local", version = "3.0")
#' bp.gf <- bioplex2graphframe(bb.gr, sc) 
#'
#' @export
bioplex2graphframe <- function(bp.gr, spark.con) 
{
    # node df
    nattrs <- names(graph::nodeDataDefaults(bp.gr))
    .getNCol <- function(n) graph::nodeData(bp.gr, graph::nodes(bp.gr), n)
    cols <- lapply(nattrs, .getNCol)
    lens <- vapply(cols, function(x) all(lengths(x) == 1), logical(1))
    .collna <- function(y) ifelse(length(y) == 1 && is.na(y),
                                  NA_character_,
                                  paste(y, collapse = ","))
    .collapse <- function(x) vapply(x, .collna, character(1))
    cols[!lens] <- lapply(cols[!lens], .collapse)
    cols <- lapply(cols, unlist)
    names(cols) <- nattrs
    node.df <- data.frame(id = graph::nodes(bp.gr), cols)
    colnames(node.df) <- tolower(colnames(node.df))
    node.df <- sparklyr::copy_to(spark.con, node.df)
    
    # edge df
    edge.df <- stack(graph::edges(bp.gr))
    colnames(edge.df) <- c("dst", "src")
    edge.df <- edge.df[,2:1]
    for(i in 1:2) edge.df[,i] <- as.vector(edge.df[,i]) 
    
    eattrs <- names(graph::edgeDataDefaults(bp.gr))
    .getECol <- function(n) graph::edgeData(bp.gr, 
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
