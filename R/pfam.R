############################################################
# 
# author: Ludwig Geistlinger
# date: 2022-03-02 11:02:10
# 
# descr: PFAM domain-domain association analysis
# 
############################################################

#' @title PFAM domain-domain association analysis
#' @description Identifies statistically associated protein domain pairs
#' connected by a disproportionally high number of protein-protein interactions. 
#' @param gr An object of class \code{\linkS4class{graph}} storing the BioPlex PPIs.
#' Typically obtained via \code{\link{bioplex2graph}}.
#' @return A \code{data.frame} containing all tested PFAM domain pairs ordered by
#' association strength.
#' @details Given PFAM domain annotations for each node of the graph, the function 
#' assesses domain pairs connected by two or more interactions for significance
#' using Fisher’s exact test based on 1) the number of interactions connecting both
#' domains; 2) the numbers of interactions involving either domain individually;
#' and 3) the number of interactions not involving either domain. 
#'
#' It is important to note that although the domain association analysis identifies
#' pairs of PFAM domains whose parent proteins interact preferentially across the
#' network, this does not necessarily mean that these domains themselves are
#' responsible for the interaction. Many may simply be passengers that associate
#' as a consequence of interactions mediated by contacts elsewhere in the protein
#' sequence. 
#'
#' @references
#' Huttlin et al. Dual proteome-scale networks reveal cell-specific remodeling
#' of the human interactome. Cell, 184(11):3022-3040.e28, 2021.
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils combn
#' @export
testDomainAssociation <- function(gr)
{
    # create a map from PFAM to UNIPROT
    unip2pfam <- graph::nodeData(gr, graph::nodes(gr), "PFAM")
    pfam2unip <- stack(unip2pfam)
    pfam2unip <- split(as.character(pfam2unip$ind), pfam2unip$values)

    # create a map from UNIPROT to PFAMs of interaction partners
    .getPF <- function(n) graph::nodeData(gr, graph::edges(gr)[[n]], "PFAM")
    unip2iapfams <- lapply(graph::nodes(gr), .getPF)
    unip2iapfams <- lapply(unip2iapfams, unlist)
    names(unip2iapfams) <- graph::nodes(gr)

    # for each pfam, get the pfams of interacting proteins
    pfam2iapfams <- lapply(pfam2unip, function(ps) unlist(unip2iapfams[ps]))
    pfam2iapfams <- pfam2iapfams[lengths(pfam2iapfams) > 0]
    pfam2pfam <- stack(pfam2iapfams)
    pfam2pfam$ind <- as.character(pfam2pfam$ind)
    pfam2pfam <- pfam2pfam[!is.na(pfam2pfam$values),]
    pfam2pfam <- as.matrix(pfam2pfam)

    # compute an all-against-all matrix storing pfam2pfam interaction counts
    ## initialize + turn pfam2pfam mapping into an integer index mapping 
    pfams <- unique(as.vector(pfam2pfam))
    len <- length(pfams)
    ia.mat <- matrix(0, nrow = len, ncol = len)
    rownames(ia.mat) <- colnames(ia.mat) <- pfams
    ind1 <- match(pfam2pfam[,1], pfams)
    ind2 <- match(pfam2pfam[,2], pfams)
    pfam2pfam <- cbind(ind1, ind2)
    
    ## loop over all pfam-pfam pair and count the PPIs connecting them
    for(i in seq_len(len))
    {
        p2p <- pfam2pfam[pfam2pfam[,1] == i |  pfam2pfam[,2] == i,,drop = FALSE]
        if(nrow(p2p))
        {
            for(j in i:len)
            {
                if(i == j) ia.mat[i,j] <- sum(p2p[,1] == j & p2p[,2] == j)
                else ia.mat[i,j] <- sum(p2p[,1] == j | p2p[,2] == j)
            }
        }
    }

    # filter out pfams that have not at least one pfam with 2 interactions
    gr1 <- apply(ia.mat, 1, function(x) any(x > 1))
    ia.mat <- ia.mat[gr1 , gr1]

    # calculate 2x2 contigency tables for all domain pairs
    ## some precomputations of the margins for efficiency
    total <- sum(as.vector(ia.mat))
    m <- ia.mat
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    rs <- rowSums(m)

    ## get all possible domain pairs and restrict to those pairs connected
    ## by at least 2 PPIs
    combs <- combn(rownames(m), 2)
    ind <- apply(combs, 2, function(x) m[x[1], x[2]] > 1)
    combs <- combs[,ind]

    ## calculate contigency for all remaining domain pairs
    conts <- apply(combs, 2, function(x) .getContingency(x[1], x[2], m, rs, total))
    conts <- t(conts)

    # test each domain pair for over-representation of connecting PPIs
    # based on the hypergeometric distribution.
    ps <- apply(conts, 1, .fisherp)
    adjp <- p.adjust(ps, method = "BH")

    res <- data.frame(PFAM1 = combs[1,],
                  PFAM2 = combs[2,],
                  PVAL = ps,
                  ADJ.PVAL = adjp)
    res <- res[order(ps),]
    return(res)    
}

# helper function that computes for each domain pair:
# 1) the number of interactions connecting both domains;
# 2) the numbers of interactions involving either domain individually; and
# 3) the number of interactions not involving either domain
.getContingency <- function(p1, p2, m, rs, total)
{
    both <- m[p1,p2]
    one <- rs[p1] - both
    two <- rs[p2] - both
    all <- total - one - two - both
    unname(c(both, one, two, all))
}

.fisherp <- function(x) fisher.test(matrix(x, nrow = 2),
                                    alternative = "greater")$p.value


# Can we do some of this vectorized ?
.f <- function(x, rs, total)
{
  one <- rs[1] - x
  two <- rs[2] - x
  all <- total - one - two - both
  cbind(x, one, two, all)
}

