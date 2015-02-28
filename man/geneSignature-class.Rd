\name{geneSignature-class}
\Rdversion{1.1}
\docType{class}
\alias{geneSignature-class}
\alias{analysisPipeline,ExpressionSet,geneSignature-method}
\alias{ensembleAdjustable,ExpressionSet,geneSignature-method}
\alias{ensembleAdjustable,matrix,geneSignature-method}
\alias{getDirect,geneSignature-method}
\alias{getGenes,geneSignature-method}
\alias{getNGenes,geneSignature-method}
\alias{getThresholds,geneSignature-method}
\alias{setDirect,geneSignature,numeric-method}
\alias{setGenes,geneSignature-method}
\alias{setGeneSignature,geneSignature-method}
\alias{setThresholds,geneSignature-method}
\alias{getGenes,solnSpace-method}
\alias{getThresholds,solnSpace-method}
\alias{setGenes,solnSpace-method}
\alias{setThresholds,solnSpace-method}

\title{Class \code{"geneSignature"}}
\description{
The geneSignature object contains the necessary elements defining the signaling environment on which a prognostic gene signature will be created.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("geneSignature", ...)}.
Objects all contain 4 slots - geneSet, geneDirect, thresholds, dirMat (unused).
}
\section{Slots}{
  \describe{
    \item{\code{geneSet}:}{Object of class \code{"character"} ~~ }
    \item{\code{geneDirect}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{thresholds}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dirMat}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{analysisPipeline}{\code{signature(dataSet = "ExpressionSet", geneSig = "geneSignature")}: ... }
    \item{ensembleAdjustable}{\code{signature(dataSet = "ExpressionSet", geneSig = "geneSignature")}: ... }
    \item{ensembleAdjustable}{\code{signature(dataSet = "matrix", geneSig = "geneSignature")}: ... }
    \item{getDirect}{\code{signature(g = "geneSignature")}: ... }
    \item{getGenes}{\code{signature(g = "geneSignature")}: ... }
    \item{getNGenes}{\code{signature(g = "geneSignature")}: ... }
    \item{getThresholds}{\code{signature(g = "geneSignature")}: ... }
    \item{setDirect}{\code{signature(g = "geneSignature", direct = "numeric")}: ... }
    \item{setGenes}{\code{signature(g = "geneSignature")}: ... }
    \item{setGeneSignature}{\code{signature(g = "geneSignature")}: ... }
    \item{setThresholds}{\code{signature(g = "geneSignature")}: ... }
	 }
}
\references{
Lee U, Frankenberger C, Yun J, Bevilacqua E, Caldas C, et al. (2013) A Prognostic Gene Signature for Metastasis-Free Survival of Triple Negative Breast Cancer Patients. PLoS ONE 8(12): e82125. doi:10.1371/journal.pone.0082125
}
\author{
UnJin lee
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("geneSignature")
}
\keyword{classes}
