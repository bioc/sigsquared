require(Biobase)
require(survival)


#define a class called solnSpace which will hold information regarding our
#solution space, including input parameters, statistical properties, as well 
#as identifiers containing which samples were trained on and tested on.
#the indices of the objects in each slot should match
setClass("solnSpace", representation(thresholds="matrix", statsTrain="matrix",
	statsTest="matrix", identTrain="matrix", identTest="matrix"))
	
#function to add a new solution to a solnSpace object
setGeneric("addSoln", signature=c("s", "thresholds", "statsTrain", "statsTest", "identTrain", "identTest"),
	function(s, thresholds, statsTrain, statsTest, identTrain, identTest) standardGeneric("addSoln"))

setMethod("addSoln", c(s="solnSpace"),
	function(s, thresholds, statsTrain, statsTest, identTrain, identTest){
		s@thresholds <- rbind(s@thresholds, thresholds)
		s@statsTrain <- rbind(s@statsTrain, statsTrain)
		s@statsTest <- rbind(s@statsTest, statsTest)
		s@identTrain <- rbind(s@identTrain, identTrain)
		s@identTest <- rbind(s@identTrest, identTest)
		
		return(s)
	}) 

#function to combine two solnspace objects
setGeneric("combineSolnSpace", signature=c("x", "y"),
	function(x, y) standardGeneric("combineSolnSpace"))		

setMethod("combineSolnSpace", c(x="solnSpace", y="solnSpace"),
	function(x, y){
		x@thresholds <- rbind(x@thresholds, y@thresholds)
		x@statsTrain <- rbind(x@statsTrain, y@statsTrain)
		x@statsTest <- rbind(x@statsTest, y@statsTest)
		x@identTrain <- rbind(x@identTrain, y@identTrain)
		x@identTest <- rbind(x@identTest, y@identTest)
		
		return(x)
	})
	
#function to subset a solnspace object
setGeneric("subsetSolnSpace", signature=c("x", "index"),
	function(x, index) standardGeneric("subsetSolnSpace"))
	
setMethod("subsetSolnSpace", c(x="solnSpace", index="numeric"),
	function(x, index){
		x@thresholds <- x@thresholds[index,]
		x@statsTrain <- x@statsTrain[index,]
		x@statsTest <- x@statsTest[index,]
		x@identTrain <- x@identTrain[index,]
		x@identTest <- x@identTest[index,]
	
		return(x)
	})	
	
#helpers to return slots in a solnSpace object
setMethod("getThresholds", c(g="solnSpace"), 
	function(g){
		return(g@thresholds)
	})
	
setMethod("getGenes", c(g="solnSpace"),
	function(g){
		return(colnames(g@thresholds))
	})

setGeneric("getStatsTrain", signature = "s",
	function(s) standardGeneric("getStatsTrain"))
	
setMethod("getStatsTrain", c(s="solnSpace"), 
	function(s){
		return(s@statsTrain)
	})
	
setGeneric("getStatsTest", signature = "s",
	function(s) standardGeneric("getStatsTest"))
	
setMethod("getStatsTest", c(s="solnSpace"), 
	function(s){
		return(s@statsTest)
	})
	
setGeneric("getIdentTrain", signature = "s",
	function(s) standardGeneric("getIdentTrain"))
	
setMethod("getIdentTrain", c(s="solnSpace"), 
	function(s){
		return(s@identTrain)
	})	
	
setGeneric("getIdentTest", signature = "s",
	function(s) standardGeneric("getIdentTest"))
	
setMethod("getIdentTest", c(s="solnSpace"), 
	function(s){
		return(s@identTest)
	})
	
#helpers to set slots in a solnSpace object
setMethod("setGenes", c(g="solnSpace"),
	function(g, genes){
		colnames(g@thresholds) <- genes
		return(g)
	})	

setGeneric("setSolnSpace", signature = "s",
	function(s, thresholds=NA, statsTrain=NA, statsTest=NA, identTrain=NA, 
		identTest=NA) standardGeneric("setSolnSpace"))

setMethod("setSolnSpace", c(s="solnSpace"),
	function(s, thresholds, statsTrain, statsTest, identTrain, identTest){
		s@thresholds <- thresholds
		s@statsTest <- statsTest
		s@statsTrain <- statsTrain
		s@identTrain <- identTrain
		s@identTest <- identTest
		return(s)
	})

setMethod("setThresholds", c(g="solnSpace"), 
	function(g, thresholds){
		g@thresholds <- thresholds
		return(g)
	})

setGeneric("setStatsTrain", signature = "s",
	function(s, d) standardGeneric("setStatsTrain"))
	
setMethod("setStatsTrain", c(s="solnSpace"), 
	function(s, d){
		s@statsTrain <- d
		return(s)
	})
	
setGeneric("setStatsTest", signature = "s",
	function(s, d) standardGeneric("setStatsTest"))
	
setMethod("setStatsTest", c(s="solnSpace"), 
	function(s, d){
		s@statsTest <- d
		return(s)
	})
	
setGeneric("setIdentTrain", signature = "s",
	function(s, d) standardGeneric("setIdentTrain"))
	
setMethod("setIdentTrain", c(s="solnSpace"), 
	function(s, d){
		s@identTrain <- d
		return(s)
	})	
	
setGeneric("setIdentTest", signature = "s",
	function(s, d) standardGeneric("setIdentTest"))
	
setMethod("setIdentTest", c(s="solnSpace"), 
	function(s, d){
		s@identTest <- d
		return(s)
	})	

#this function will apply a solnSpace's set of thresholds onto a data set
#this function is parallelized, with maximum used cores as nCores
setGeneric("applySigSolnSpace", signature=c("dataSet", "geneSig", "sp"),
	function(dataSet, geneSig, sp, MFS="MFS", met="met", index=NA) standardGeneric("applySigSolnSpace"))
	
setMethod("applySigSolnSpace", c(dataSet ="ExpressionSet", geneSig="geneSignature", sp="solnSpace"),
	function(dataSet, geneSig, sp, MFS="MFS", met="met", index){
		nTh <- nrow(sp@thresholds)
		gs <- geneSig
		overlap <- geneSig@geneSet %in% colnames(sp@thresholds)
		gs@geneSet <- 	geneSig@geneSet[overlap]
		gs@geneDirect <- geneSig@geneDirect[overlap]
		
		result <- matrix(0, nrow=nTh, ncol=2)
		colnames(result) <- c("pval", "size")
		
		if(!is.integer(index)) index <- 1:ncol(exprs(dataSet))
		
		if(!mc){
			for(i in 1:nTh){
				gs@thresholds <- sp@thresholds[i,]
				result[i, "pval"] <- applySigLR(dataSet, gs, MFS, met, index)
				result[i, "size"] <- applySigSize(dataSet, gs, TRUE, index=index)
			}
		}
		
		else{
			processes <- list()
			iterPerCore <- floor(nTh/nCores)
			
			#x in applySig(x) is a vector of indices of sp@thresholds
			#over which we are computing
			applySig.local <- function(x) {
					result.local <- matrix(0, nrow=length(x), ncol=2)
					colnames(result.local) <- c("pval", "size")
					gs.local <- gs
					
					for(i in x){
						gs.local@thresholds <- sp@thresholds[i,]
						result.local[which(x==i), "pval"] <- applySigLR(dataSet, gs.local, MFS, met, index)
						result.local[which(x==i), "size"] <- applySigSize(dataSet, gs.local, TRUE, index=index)
					}
					
					return(result.local)
				}
			
			for(i in 1:(nCores-1)){
				processes[[i]] <- mcparallel(applySig.local(((i-1)*iterPerCore+1):(i*iterPerCore)))
			}
			
			processes[[nCores]] <- mcparallel(applySig.local(((nCores-1)*iterPerCore+1):nTh))
			
			resultCollect <- mccollect(processes)
	
			result <- do.call("rbind", resultCollect)	
		}
		
		return(result)
	})

#this function will take a solnSpace object and summarize the solnSpace data
setGeneric("summarizeSolnSpace", signature=c("sp"),
	function(sp, train=FALSE, test=TRUE) standardGeneric("summarizeSolnSpace"))
	
setMethod("summarizeSolnSpace", c(sp = "solnSpace"),
	function(sp, train, test){
		result <- vector()
		ind <- 1
		if(train){
			result[ind] <- mean(sp@statsTrain[,"pval"])
			names(result)[ind] <- "pvalTrainMean"
			ind <- ind+1
			result[ind] <- sd(sp@statsTrain[,"pval"])
			names(result)[ind] <- "pvalTrainSD"
			ind <- ind+1
			
			result[ind] <- mean(sp@statsTrain[,"size"])
			names(result)[ind] <- "sizeTrainMean"
			ind <- ind+1
			result[ind] <- sd(sp@statsTrain[,"size"])
			names(result)[ind] <- "sizeTrainSD"
			ind <- ind+1
		}
		
		if(test){
			result[ind] <- mean(sp@statsTest[,"pval"])
			names(result)[ind] <- "pvalTestMean"
			ind <- ind+1
			result[ind] <- sd(sp@statsTest[,"pval"])
			names(result)[ind] <- "pvalTestSD"
			ind <- ind+1
			
			result[ind] <- mean(sp@statsTest[,"size"])
			names(result)[ind] <- "sizeTestMean"
			ind <- ind+1
			result[ind] <- sd(sp@statsTest[,"size"])
			names(result)[ind] <- "sizeTestSD"
			ind <- ind+1
		}	
		
		return(result)
	})	
