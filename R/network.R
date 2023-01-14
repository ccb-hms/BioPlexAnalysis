############################################################
# 
# author: Ludwig Geistlinger
# date: 2023-01-13 11:19:25
# 
# descr: functionality for working with network data 
# 
############################################################

#' Join multiple networks around a target node ID
#' 
#' @description This function joins multiple networks and resticts the joined
#' network to interactions involving a target node ID of interest 
#' @param id character. A gene or protein ID around which a joined network
#' should be constructed.
#' @param network list. Each element of the list is a \code{data.frame} representing
#' a network. That means each row of a \code{data.frame} represents an interaction.
#' @param intA.col character. Name of the column that stores the ID of the first
#' interactor in each \code{data.frame} of argument \code{network}.
#' @param intB.col character. Name of the column that stores the ID of the second
#' interactor in each \code{data.frame} of argument \code{network}.
#' @return A \code{data.frame} containing the joined network. 
#' @examples
#'  library(BioPlex)
#'  bp.293t <- BioPlex::getBioPlex(cell.line = "293T", version = "3.0")
#'  bp.hct <- BioPlex::getBioPlex(cell.line = "HCT116", version = "1.0")
#'  netl <- list(HEK = bp.293t, HCT = bp.hct)
#'
#'  jnet <- joinNetworks(id = "CDH2", network = netl)
#'
#' @export
joinNetworks <- function(id, network, intA.col = "SymbolA", intB.col = "SymbolB")
{
    stopifnot(is.list(network))
    stopifnot(all(vapply(network, is.data.frame, logical(1))))

    # get all interactions involving target ID
    .getIA <- function(id, net)
    { 
        ind <- net[,intA.col] == id | net[,intB.col] == id
        net[ind,] 
    }
    id.net <- lapply(network, .getIA, id = id)

    # expand by including interactions between interactors of target ID
    idi <- lapply(id.net, function(x) c(x[,intA.col], x[,intB.col]))
    idi <- Reduce(union, idi)
    idi <- setdiff(idi, id)
    .expandIA <- function(net)
    {
        ind <- net[,intA.col] %in% idi & net[,intB.col] %in% idi
        net[ind,]
    }
    exp.id.net <- lapply(network, .expandIA)

    # construct joined network of interactions involving target ID or one of its
    # interactors 
    jnet <- mapply(rbind, id.net, exp.id.net, SIMPLIFY = FALSE)
    for(i in seq_along(jnet)) jnet[[i]][["net"]] <- names(jnet)[i]
    jnet <- do.call(rbind, jnet)
    return(jnet)
}
