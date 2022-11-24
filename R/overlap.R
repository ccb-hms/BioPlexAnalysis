#' Assess coverage of a graph on a set of node IDs
#' 
#' @description This function performs a statistical assessment of the
#' overlap (number of edges) between a graph and a set of node IDs. This
#' is useful for eg assessing the overlap (number of edges) between a PPI
#' network and a set of node IDs representing protein sets of interests such
#' as eg. a CORUM complex, a KEGG pathway, or a set of proteins associated with
#' a disease.
#' @details The test is based on network randomization. In each replication,
#' the network is randomized/re-wired and the number of edges within the input
#' set of node IDs is compared against the observed number of edges. The observed
#' number of edges corresponds here to the number of edges within the input set of
#' node IDs based on the true / non-randomized version of the input graph.
#' @param ids character. A set of ids. Need at least three IDs present
#' as nodes in the graph.
#' @param gr graph. An object of class \code{graphNEL} from the graph package.
#' @param nr.reps integer. Number of replications. Defaults to 1000.
#' @return A p-value that is calculated as for a permutation test, ie as the relative
#' frequency of obtaining as many or more edges in a randomized setup when compared
#' to the observed number of edges within the input set of node IDs.
#' @examples
#'   library(graph)
#'   n <- LETTERS[1:10]
#'   edl <- vector("list", length = 10)
#'   names(edl) <- n
#'   for(i in 1:10) edl[[i]] <- list(edges = 11 - i)
#'   gr <- graphNEL(nodes = n, edgeL = edl)
#'   
#'   ids <- LETTERS[c(1:2, 9:10)]
#'   testCoverage(ids, gr, nr.reps = 100) 
#'
#' @export
testCoverage <- function(ids, gr, nr.reps = 1000)
{
    n <- intersect(graph::nodes(gr), ids)
    if(length(n) < 3) 
        stop("Need >= 3 IDs present as nodes in the provided graph")
    if(length(n) < length(ids))
        message(length(n), " of ", length(ids),
                " IDs present as nodes in the provided graph")
    sgr <- graph::subGraph(n, gr)
    obs.nr.edges <- graph::numEdges(sgr)
    igr <- igraph::graph_from_graphnel(gr)
    rand.nr.edges <- replicate(nr.reps, .countRandomizedEdges(n, igr))
    p <- (sum(rand.nr.edges >= obs.nr.edges) + 1) / (nr.reps + 1)
    return(p)
}

.countRandomizedEdges <- function(n, igr)
{
    rew.igr <- BiRewire::birewire.rewire.undirected(igr, verbose = FALSE)
    sgr <- igraph::subgraph(rew.igr, n)
    nr.edges <- igraph::ecount(sgr)
    return(nr.edges)
}
