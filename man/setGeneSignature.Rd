\name{setGeneSignature}
\alias{setGeneSignature}
\alias{setDirect}
\alias{setThresholds}
\alias{setGenes}
\alias{getDirect}
\alias{getThresholds}
\alias{getGenes}
\alias{getNGenes}

\title{geneSignature functions}

\description{
The geneSignature object contains the necessary elements defining the signaling environment on which a prognostic gene signature will be created.  This collection of functions are used to manipulate or retrieve the data slots of a given geneSignature object.}

\usage{
setGeneSignature(g, direct=NA, thresholds=c(0), genes=NA, mat=matrix())
setDirect(g, direct)
setThresholds(g, thresholds)
setGenes(g, genes)
getDirect(g)
getThresholds(g)
getGenes(g)
getNGenes(g)
}

\arguments{
	\item{g}{
		geneSignature object
		}
	\item{direct}{
		vector of -1s or 1s representing down- or up-regulation respectively
		}
	\item{thresholds}{
		vector of values containing thresholds for the geneSignature object
		}
	\item{genes}{
		character vector of gene names
		}
	\item{mat}{
		matrix of interactions between genes (unused)
		}
}

\value{All setting functions return objects of class geneSignature. getDirect yields a vector of -1s or 1s, getThesholds yields a vector of theshold values, getGenes yields a character vector of gene names, getNGenes yields the number of genes in the geneSignature}

\author{UnJin Lee}

\examples{
## Generate and read out values of a geneSignature object
gs <- setGeneSignature(new("geneSignature"), c(1, 1), c(0, 0), c("BACH1", "RKIP"), matrix())
getDirect(gs)
getThresholds(gs)
getGenes(gs)
getNGenes(gs)
}
