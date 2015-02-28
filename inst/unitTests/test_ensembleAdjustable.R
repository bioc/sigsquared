test_ensembleAdjustable <- function(){
	data(BrCa443)
	
	gs <- new("geneSignature")
	gs@geneSet <- c("HMGA2", "SPP1", "CXCR4")
	gs@geneDirect <- c(-1, 1, 1)
	gs@thresholds <- c(0.3, 0.02, -0.17)
	
	nPos <- sum(ensembleAdjustable(BrCa443, gs))

	gs1 <- new("geneSignature")
	gs1@geneSet <- c("MMP1", "MetaBACH1", "RKIP")
	gs1@geneDirect <- c(-1, 1, -11)
	gs1@thresholds <- c(0.02, 0.314, -0.22)
	
	nPos1 <- sum(ensembleAdjustable(BrCa443, gs1))

	checkEqualsNumeric(nPos, 99)
	checkEqualsNumeric(nPos1, 37)
}
