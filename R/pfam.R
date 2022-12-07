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
#' Typically obtained via \code{bioplex2graph}.
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
#'
#' @examples
#'
#'   # (1) obtain BioPlex PPI network for 293T cells
#'   library(BioPlex)
#'   df <- getBioPlex(cell.line = "HCT116", version = "1.0")
#'   hct.gr <- bioplex2graph(df)
#'  
#'   # (2) annotate PFAM domains
#'   library(AnnotationHub)
#'   ah <- AnnotationHub()
#'   orgdb <- query(ah, c("orgDb", "Homo sapiens"))
#'   orgdb <- orgdb[[1]]
#'   hct.gr <- annotatePFAM(hct.gr, orgdb)
#'
#'   # (3) domain-domain association analysis
#'   res <- testDomainAssociation(hct.gr)
#'
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils combn relist
#' @export
testDomainAssociation <- function(gr)
{
    # create a map from PFAM to UNIPROT
    unip2pfam <- graph::nodeData(gr, graph::nodes(gr), "PFAM")
    pfam2unip <- stack(unip2pfam)
    pfam2unip <- split(as.character(pfam2unip$ind), pfam2unip$values)

    # create a map from UNIPROT to PFAMs of interaction partners
    unip2iapfams <- graph::nodeData(gr, unlist(graph::edges(gr)), "PFAM")
    unip2iapfams <- relist(unip2iapfams, graph::edges(gr))
    unip2iapfams <- lapply(unip2iapfams, unlist)
    names(unip2iapfams) <- graph::nodes(gr)

    # for each pfam, get the pfams of interacting proteins
    pfam2iapfams <- lapply(pfam2unip, function(ps) unlist(unip2iapfams[ps]))
    pfam2iapfams <- pfam2iapfams[lengths(pfam2iapfams) > 0]
    pfam2pfam <- stack(pfam2iapfams)
    pfam2pfam$ind <- as.character(pfam2pfam$ind)
    pfam2pfam <- pfam2pfam[!is.na(pfam2pfam$values),]
    tab <- table(pfam2pfam)
    tab <- as.data.frame.table(tab)

    ig <- igraph::graph_from_data_frame(tab, directed = FALSE)
    m <- igraph::get.adjacency(ig, attr = "Freq", sparse = TRUE)
    m2 <- m <- as.matrix(m)
    m2[lower.tri(m2)] <- NA
    tab <- reshape2::melt(m2, na.rm = TRUE)

    # calculate 2x2 contigency tables for all domain pairs
    ## some precomputations of the margins for efficiency
    total <- sum(tab$value)
    rs <- rowSums(m)

    # restrict to domain pairs connected by at least 2 PPIs 
    tab <- subset(tab, value > 1)

    ## calculate contigency for all remaining domain pairs
    conts <- apply(tab, 1, function(x) .getContingency(x[1], x[2], m, rs, total))
    conts <- t(conts)

    # test each domain pair for over-representation of connecting PPIs
    # based on the hypergeometric distribution.
    ps <- apply(conts, 1, .fisherp)
    adjp <- p.adjust(ps, method = "BH")

    res <- data.frame(PFAM1 = tab[,1],
                  PFAM2 = tab[,2],
                  FREQ = tab[,3],  
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
  all <- total - one - two - x
  cbind(x, one, two, all)
}

