\name{ensembleAdjustable}
\alias{ensembleAdjustable}

\title{Application of geneSignature object}

\description{
The ensembleAdjustable function applies a geneSignature object to a data matrix  containing expression values and gene symbols or an ExpressionSet object.
}

\usage{
ensembleAdjustable(dataSet, geneSig, index=F)
}


\arguments{
	\item{dataSet}{
		data set object, may be numeric matrix or an ExpressionSet
		}
	\item{geneSig}{
		geneSignature object containing directions, thresholds, and gene symbols
		}
	\item{index}{
		index to indicate which samples are to be subsetted, may be FALSE for no subsetting or a vector of column numbers
		}
}

\value{A logical vector with length equal to the number of samples (or samples subsetted), TRUE indicating a positive, FALSE indicating a negative}

\author{UnJin Lee}

\examples{
require(Biobase)
## Generate test geneSignature object with 0s for thresholds
gs <- setGeneSignature(g=new("geneSignature"), direct=c(1,1,1), genes=c("A", "B", "C"), thresholds=c(0, 0, 0))

## Generate randomly distributed matrix and ExpressionSet
mat <- matrix(rnorm(9, 0, 1), nrow=3)
rownames(mat) <- c("A", "B", "C")
posmat <- abs(mat)
expset <- new("ExpressionSet", exprs=mat)

## Apply geneSignature to matrices
ensembleAdjustable(mat, gs)
ensembleAdjustable(posmat, gs)

## Apply geneSignature to ExpressionSet
ensembleAdjustable(expset, gs)

## Apply geneSignature with subsetting
ensembleAdjustable(mat, gs, c(1, 3))
ensembleAdjustable(expset, gs, c(1, 3))
}
