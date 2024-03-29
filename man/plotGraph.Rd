% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotGraph}
\alias{plotGraph}
\title{Plot a graph}
\usage{
plotGraph(
  gr,
  node.label = "SYMBOL",
  node.data = NULL,
  edge.data = NULL,
  node.color = "firebrick",
  edge.color = "firebrick",
  edge.mode = c("directed", "undirected"),
  layout = c("fruchterman.reingold", "circle")
)
}
\arguments{
\item{gr}{graph. An object of class \code{graphNEL} from the graph package.}

\item{node.label}{character. A node data attribute that is used to label the
nodes.}

\item{node.data}{character. A node data attribute that is used to overlay on
the nodes using different colors.}

\item{edge.data}{character. An edge data attribute that is used to overlay
on the edges using different line types.}

\item{node.color}{character. Default color for the nodes. Is overwritten by
\code{node.data} if not null.}

\item{edge.color}{character. Default color for the edges.}

\item{edge.mode}{character. Whether to display directed or undirected edges.}

\item{layout}{character. Circular layout or network layout.}
}
\value{
A ggplot.
}
\description{
Plotting utility for graphs
}
\examples{
     library(BioPlex)
     corum.df <- getCorum("core")
     corum.glist <- corum2graphlist(corum.df)
     gr <- corum.glist[["CORUM107_TFIIH_transcription_factor_complex"]]

     plotGraph(gr, edge.mode = "undirected", 
               edge.color = "grey", layout = "circle")

}
