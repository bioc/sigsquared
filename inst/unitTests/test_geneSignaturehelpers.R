test_getDirect<- function(){
	genedir <- c(-1, 1, -1, -1, 1)
	gs <- new("geneSignature")
	gs@geneDirect <- genedir

	checkIdentical(genedir, getDirect(gs))
}

test_getThresholds <- function(){
	genethresh <- c(0.25, 0.1, 0.23, -0.1)
	gs <- new("geneSignature")
	gs@thresholds <- genethresh

	checkIdentical(genethresh, getThresholds(gs))
}

test_getGenes <- function(){
	genenames <- c("A", "B", "C", "DEF")
	gs <- new("geneSignature")		
	gs@geneSet <- genenames

	checkIdentical(genenames, getGenes(gs))
}
	
test_getNGenes <- function(){
	genenames <- c("A", "B", "C", "DEF")
	gs <- new("geneSignature")
	gs@geneSet <- genenames		
	
	checkEqualsNumeric(4, getNGenes(gs))
}		

test_setDirect <- function(){
	genedir <- c(-1, 1, -1, -1, 1)
	gs <- new("geneSignature")
	gs@geneDirect <- genedir

	gs1 <- new("geneSignature")
	gs1 <- setDirect(gs1, genedir)

	checkIdentical(gs, gs1)
	}

test_setThresholds <- function(){
	genethresh <- c(0.25, 0.1, 0.23, -0.1)
	gs <- new("geneSignature")
	gs@thresholds <- genethresh

	gs1 <- new("geneSignature")
	gs1 <- setThresholds(gs1, genethresh)

	checkIdentical(gs, gs1)
}

test_setGenes <- function(){
	genenames <- c("A", "B", "C", "DEF")
	gs <- new("geneSignature")		
	gs@geneSet <- genenames

	gs1 <- new("geneSignature")
	gs1 <- setGenes(gs1, genenames)

	checkIdentical(gs, gs1)
}
	

