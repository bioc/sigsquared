\name{analysisPipeline}
\alias{analysisPipeline}

\title{Training of thresholds}

\description{
The analysisPipeline function is used to train a set of thresholds for predicting survival outcome within the context of a given signaling environment. This signaling environment is encoded in a geneSignature object.
}

\usage{
analysisPipeline(dataSet, geneSig, iterPerK=2500, k=3, rand=TRUE, newjpdf=FALSE, jpdf=FALSE, nJPDF=12500, disc=c(0.005, 0.01, 0.03, 0.05), MFS="MFS", met="met", optMeth="Nelder-Mead")
}

\arguments{
	\item{dataSet}{
		ExpressionSet object containing both expression data (exprs) and phenotypic survival data (pData)
		}
	\item{geneSig}{
		geneSignature object containing directions, thresholds, and gene symbols
		}
	\item{iterPerK}{
		integer number of optimization iterations for each k
		}
	\item{k}{
		integer k for k-fold cross-validation
		}
	\item{rand}{
		boolean determining whether the k subsets are randomly drawn (otherwise k subsets are selected ordinally)
		}
	\item{newjpdf}{
		boolean for generating a joint probability function for alternate smoothed cost function (not recommended)
		}
	\item{jpdf}{
		solnSpace object containing empirical joint probability function for alternate smoothed cost function (not recommended)
		}
	\item{nJPDF}{
		value determining the number of samples with which to estimate the empirical joint probability function for alternate smoothed cost function (not recommended)
		}
	\item{disc}{
		vector of discretation thresholds for discretized cost function
		}
	\item{MFS}{
		variable name for survival-time data in dataSet object
		}
	\item{met}{
		variable name for metastasis event data in dataSet object
		}
	\item{optMeth}{
		optimization method used by R function 'optim'
		}
}

\details{
The analysisPipeline function optimizes over a cost function designed to minize both type I and II error. There is a discretized and smoothed cost function available, however implementation of the smoothed cost function relies on sampling of the solution space. This sampling may be pre-computed and implemented through the 'jpdf' argument, however overall usage of the smoothed cost function is not recommended. 
}

\value{A geneSignature object containing newly trained thresholds}

\author{UnJin Lee}

\examples{
## Load in example data
data("BrCa443")

## Create initial geneSignature object
## Note it is not necessary to define thresholds at this point
gs <- setGeneSignature(g=new("geneSignature"), direct=c(-1,1,1,1,1,1,1), genes=c("RKIP", "HMGA2", "SPP1", "CXCR4", "MMP1", "MetaLET7", "MetaBACH1"))

## Generate thresholds
gs <- analysisPipeline(dataSet=BrCa443, geneSig=gs, iterPerK=50, k=2, rand=FALSE)
}

