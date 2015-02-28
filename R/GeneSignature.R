require(Biobase)
require(survival)

mc <<- F
##unused parallel functionality
##	if("parallel" %in% rownames(installed.packages())){
##		require(parallel)
##		mc <<- exists("mcparallel")
##	}

nCores <<- 1
##	if(mc){
##		require(survival)
##		nCores <<- detectCores()
##		nCores <<- max(nCores - 2, 1)
##	}


#define a class called geneSignature which will hold all of our relevant information
setClass("geneSignature", representation(geneSet="character", geneDirect="numeric",
	thresholds="numeric", dirMat="matrix"))

#helpers to return slots in a geneSignature object
setGeneric("getDirect", signature = "g",
	function(g) standardGeneric("getDirect"))
	
setMethod("getDirect", c(g="geneSignature"), 
	function(g){
		return(g@geneDirect)
	})

setGeneric("getThresholds", signature = "g",
	function(g) standardGeneric("getThresholds"))
	
setMethod("getThresholds", c(g="geneSignature"), 
	function(g){
		return(g@thresholds)
	})

setGeneric("getGenes", signature = "g",
	function(g) standardGeneric("getGenes"))
	
setMethod("getGenes", c(g="geneSignature"),
	function(g){
		return(g@geneSet)
	})
	
setGeneric("getNGenes", signature = "g",
	function(g) standardGeneric("getNGenes"))
	
setMethod("getNGenes", c(g="geneSignature"),
	function(g){
		return(length(g@geneSet))
	})		

setGeneric("getDirMat", signature = "g",
	function(g) standardGeneric("getDirMat"))
	
setMethod("getDirMat", c(g="geneSignature"),
	function(g){
		return(g@dirMat)
	})
	
#helpers to adjust slots in a geneSignature object
setGeneric("setGeneSignature", signature = "g",
	function(g, direct=NA, thresholds=c(0), genes=NA, mat=matrix()) standardGeneric("setGeneSignature"))
	
setMethod("setGeneSignature", c(g="geneSignature"), 
	function(g, direct, thresholds, genes, mat){
		g@geneDirect <- direct
		g@thresholds <- thresholds
		g@geneSet <- genes
		g@dirMat <- mat
		return(g)
	})
	
setGeneric("setDirect", signature = c("g", "direct"),
	function(g, direct) standardGeneric("setDirect"))
	
setMethod("setDirect", c(g="geneSignature", direct="numeric"), 
	function(g, direct){
		g@geneDirect <- direct	
		return(g)
	})

setGeneric("setThresholds", signature = c("g", "thresholds"),
	function(g, thresholds) standardGeneric("setThresholds"))
	
setMethod("setThresholds", c(g="geneSignature"), 
	function(g, thresholds){
		g@thresholds <- thresholds	
		return(g)
	})

setGeneric("setGenes", signature = c("g", "genes"),
	function(g, genes) standardGeneric("setGenes"))
	
setMethod("setGenes", c(g="geneSignature"),
	function(g, genes){
		g@geneSet <- genes	
		return(g)
	})	

setGeneric("setDirMat", signature = c("g", "d"),
	function(g, d) standardGeneric("setDirMat"))
	
setMethod("setDirMat", c(g="geneSignature"),
	function(g, d){
		g@dirMat <- d	
		return(g)
	})
	
#function implementing neural net model with or without thresholding 
setGeneric("neuralNetTh", signature = c("vminus", "T"),
	function(type = 1,vminus, T, R, a, decay) standardGeneric("neuralNetTh"))
	
setMethod("neuralNetTh", c(vminus="numeric", T="matrix"),
	function(type, vminus, T, R, a, decay){
		g <- function(x) x
		if(type == 2) g <- function(x) (1+tanh(x))/2
		
		vbar=R/2
		if(type == 1) {
			R <- 1
			vbar <- 0
		}

		mu <- (T %*% ((vminus-vbar)/vbar)) + a

		vplus <- R*g(mu) - decay*(vminus-vbar)
		if(type == 1) vplus <- vplus + vminus
		
		vplus <- as.numeric(vplus)
		return(vplus)
	})
	
#function to take in an edge matrix and apply the neural net model to model
#time evolution of a system given the function through time of gene products,
#predicting the dynamics of the rest of the gene products.
#as the model is discrete through time, each element/column in fcn corresponds 
#to a time point
#the 'type' variable determines what model to use - with or without thresholding
#type=1 is using threshold g(u)=u 
#type=2 is using g(u)=(1+tanh(u))/2 (more "realistic" & used by genGeneDirect)
#if vinit is non-null/is a vector, this function will generate the time evolution
#of each gene product
#if fcn is non-na/is a vector or matrix of gene products (with index fcnIndex),
#this function will generate the time evolution of the other gene products in this system
#if you forget both, the function will return(0), if you apply both, 
#the function will automatically select fcn and make vinit NULL
setGeneric("timeCourseNNetFcn", signature = c("T"),
	function(T, fcn=NA, fcnIndex=1, type=2, baseline=50, decay=0.1, vinit=NULL, time=20) standardGeneric("timeCourseNNetFcn"))

setMethod("timeCourseNNetFcn", c(T="matrix"),
	function(T, fcn, fcnIndex, type, baseline, decay, vinit, time){
		if(is.na(fcn) && is.null(vinit)){
			message("Incompatible 'fcn' and 'vinit' (Both Empty)")
			return(0)
		}
			
		if(is.numeric(fcn) && is.numeric(vinit)){
			message("Incompatible 'fcn' and 'vinit'. Proceeding using 'fcn'")
			vinit <- NULL
		}
		
		dim <- nrow(T)
		notIndex <- (1:dim)[-fcnIndex]
		
		if(is.matrix(fcn)) time <- ncol(fcn)
		else if(is.numeric(fcn)) time <- length(fcn)		
		
		vInit <- rep(baseline, dim)
		if(is.vector(vinit)) vInit <- vinit
		
		if(type == 1 && baseline == 50) baseline <- 0
		results <- matrix(baseline, nrow=dim, ncol=time)
		
		if(is.null(vinit)) results[fcnIndex,] <- fcn
		else results[,1] <- vinit

		for(i in 1:(time-1)){
			vplus <- neuralNetTh(type, results[,i], T, rep(2*baseline, dim), rep(0, dim), decay)
			
			if(is.null(vinit)) results[notIndex,i+1] <- vplus[notIndex]
			else results[,i+1] <- vplus
		}
		
		return(results)
	})
	
#function to take in a edge matrix (and optionally master regulator(s)) and 
#returns a vector with threshold directions
setGeneric("genGeneDirect", signature = c("dirMat"),
	function(dirMat, mReg=1, mRegDir=1, pump=1.5, baseline=50) standardGeneric("genGeneDirect"))
	
setMethod("genGeneDirect", c(dirMat="matrix"),
	function(dirMat, mReg, mRegDir, pump, baseline){
		nDim <- dim(dirMat)
		if(!nDim[1] == nDim[2]){
			message("Not a square matrix")
			return(0)
		}
		else nDim <- nDim[1]
	
		if(!length(mReg) == length(mRegDir)){
			message("Master regulator length mismatch")
			return(0)
		}
	
		for(i in 1:length(mReg)){
			if(length(pump) == 1) dirMat[i,i] <- pump
			else dirMat[i,i] <- pump[i]
		}

		vinit <- rep(baseline, nDim)
		vinit[mReg] <- (1+mRegDir*0.1)*baseline
				
		message("Simulating network")
		tCourse <- timeCourseNNetFcn(dirMat, type=2, decay=0,
			vinit=vinit, time=20*nDim, baseline=baseline)		

		time <- ncol(tCourse)
		th <- vector()
		
		for(i in 1:nDim){	
			dirLogicUp <- tCourse[i, -1] >= tCourse[i, -time]
			dirLogicDown <- tCourse[i, -1] <= tCourse[i, -time]

			if((sd(tCourse[i,]) == 0)) th[i] <- 0
			else if(sum(dirLogicUp) == (time-1)) th[i] <- 1
			else if(sum(dirLogicDown) == (time-1)) th[i] <- -1
		}
	
		return(th)
	})
	

#This function takes in arguments as numeric vectors and
#boolean list determining if a patient is positive or negative
setGeneric("ensembleAdjustable", signature=c("dataSet", "geneSig"),
	function(dataSet, geneSig, index=F) standardGeneric("ensembleAdjustable"))
	
setMethod("ensembleAdjustable", c(dataSet="matrix", geneSig="geneSignature"),
	function(dataSet, geneSig, index){
		if(class(index) %in% c("numeric", "integer")) dataSet <- dataSet[,index]
		th <- getThresholds(geneSig)
		dir <- getDirect(geneSig)
		g <- getGenes(geneSig)
		
		netReg <- dir*dataSet[g,] > dir*th
		
		return(colSums(netReg) == length(dir))
	})		

#This function takes an expressionset object, breaks the expression data into
#numerical vectors and applies the ensembleAdjustable function to it
setMethod("ensembleAdjustable", c(dataSet="ExpressionSet", geneSig="geneSignature"),
	function(dataSet, geneSig, index){
		if(class(index) %in% c("numeric", "integer")) z <- exprs(dataSet)[,index]
		else z <- exprs(dataSet)
		
		g <- getGenes(geneSig)
		
		result <- ensembleAdjustable(z[g,], geneSig)
		return(result)
	})
	
	
#Cost function to be minimized
#can be used to adjust relative weights of log-rank and size
#can also be used to estimate a joint pdf using random solutions
setGeneric("ensembleCostFcn", signature = c("dataSet", "geneSig", "jpdf"),
	function(dataSet, geneSig, jpdf=FALSE, disc=c(0.005, 0.01, 0.03, 0.05),
		index=NA, MFS="MFS", met="met") standardGeneric("ensembleCostFcn"))

#version called if not using jpdf
setMethod("ensembleCostFcn", c(dataSet="ExpressionSet", geneSig="geneSignature", jpdf="logical"),
	function(dataSet, geneSig, jpdf, disc, index, MFS, met){
		#perc is non-positive percent (to minimize)
		perc <- 1 - applySigSize(dataSet, geneSig, TRUE, index)

		pval <- applySigLR(dataSet, geneSig, MFS, met, index)

		#we are discretizing p-values using the list disc
		disc <- sort(disc)
		pvalCategory <- sum(pval < disc) + 1
		if(pvalCategory == 1) pval <- 2*pval
		else if(pvalCategory < length(disc)) pval <- disc[pvalCategory]
				
		#we estimate perc to be on [.9, 1] while desired pvals on (0, .05]
		#so we scale perc to be 0.1*(perc-0.9)
		#we then make adjustments to punish low percents
		percScaled <- 0.5*(perc-0.9)
		costfcn <- pval + percScaled

		return(costfcn)
	})

#version called if using jpdf
setMethod("ensembleCostFcn", c(dataSet="ExpressionSet", geneSig="geneSignature", jpdf="solnSpace"),
	function(dataSet, geneSig, jpdf, disc, index, MFS, met){
		ref <- getStatsTrain(jpdf)
		
		perc <- applySigSize(dataSet, geneSig, TRUE, index)
		pval <- applySigLR(dataSet, geneSig, MFS, met, index)
		
		jpval <- sum((ref[,"size"] > perc) & (ref[,"pval"] < pval))/nrow(ref)
		
		return(jpval)
	})

#this function returns the log-rank p-value of a signature in a data set
#values MFS are for survival data, and met are for censor data within the dataSet object
setGeneric("applySigLR", signature = c("dataSet", "geneSig", "MFS", "met"),
	function(dataSet, geneSig, MFS="MFS", met="met", index=NA) standardGeneric("applySigLR"))
	
setMethod("applySigLR", c(dataSet="ExpressionSet", geneSig="geneSignature"),
	function(dataSet, geneSig, MFS, met, index){
		positive <- ensembleAdjustable(dataSet, geneSig, index)
		if(class(index) == "integer") result <- summary(coxph(Surv(MFS, met) ~ positive, pData(dataSet)[index,]))$sctest["pvalue"]
		else result <- summary(coxph(Surv(MFS, met) ~ positive, pData(dataSet)))$sctest["pvalue"]
		return(result)
	})	
	
#this function returns the number of patients in a signature within a data set
#results can be either a percent or a raw number, controlled by parameter perc	
setGeneric("applySigSize", signature = c("dataSet", "geneSig", "perc"),
	function(dataSet, geneSig, perc=TRUE, index=NA) standardGeneric("applySigSize"))

setMethod("applySigSize", c(dataSet="ExpressionSet", geneSig="geneSignature", perc="logical"),
	function(dataSet, geneSig, perc, index){
		result <- sum(ensembleAdjustable(dataSet, geneSig, index))
		if(perc){
			if(class(index) == "integer") return(result/length(index))
			else return(result/ncol(exprs(dataSet)))
		}
		return(result)
	})
	
##this function gives a baseline to estimate the null distribution in 
##p-value/size space by using random sampling centered around 0, sd=1
setGeneric("eJPDF", signature = c("dataSet", "geneSig", "n"),
	function(dataSet, geneSig, n, index=NA, MFS="MFS", met="met") standardGeneric("eJPDF"))
	
setMethod("eJPDF", c(dataSet="ExpressionSet", geneSig="geneSignature", n="numeric"),
	function(dataSet, geneSig, n, index, MFS, met){
		result <- new("solnSpace")
		
		if(mc) n <- round(n/nCores)*nCores
		
		if(is.na(index)) index <- 1:ncol(exprs(dataSet))
		
		maxIndex <- ncol(exprs(dataSet))
		ident <- matrix(NA, nrow=n, ncol=maxIndex)
		statsTrain <- matrix(NA, nrow=n, ncol=2)

		genes <- getGenes(geneSig)
		dim <- length(genes)
		
		thresholds <- matrix(rnorm(n*dim, 0, 1), ncol=dim)
		colnames(thresholds) <- getGenes(geneSig)
		if(!mc){
			for(i in 1:n){
				geneSig <- setThresholds(geneSig, thresholds[i,])
			
				pval <- applySigLR(dataSet, geneSig, MFS, met, index)
				size <- applySigSize(dataSet, geneSig, TRUE, index)
			
				statsTrain[i,] <- c(pval, size)
				ident[i,] <- rep(TRUE, maxIndex)
			}
		}
		
		else{
			mcResult <- list()
			
			mcApplySigLRSize <- function(th){
				nr <- nrow(th)
				mcST <- matrix(0, nrow=nr, ncol=2)
				for(i in 1:nr){
					geneSig <- setThresholds(geneSig, th[i,])
				
					pval <- applySigLR(dataSet, geneSig, MFS, met, index)
					size <- applySigSize(dataSet, geneSig, TRUE, index)
			
					mcST[i,] <- c(pval, size)
				} 
				
				return(mcST)	
			}
			
			chunk <- n/nCores
			
			for(i in 1:nCores){
				mcResult[[i]] <- mcparallel(mcApplySigLRSize(thresholds[((i-1)*chunk+1):(i*chunk),]))
			}
			
			for(i in 1:n){
				ident[i,] <- rep(TRUE, maxIndex)
			}
			
			mcStatsTrain <- mccollect(mcResult)
			statsTrain <- do.call("rbind", mcStatsTrain)
			
		}

		colnames(statsTrain) <- c("pval", "size")
		
		result@thresholds <- thresholds
		result@statsTrain <- statsTrain
		result@statsTest <- statsTrain
		result@identTrain <- ident
		result@identTest <- ident
		
		return(result)
	})


##function that performs the optimizations in a given data set using a given 
##expressionSet object, training and testing over two given sets of indices
##output is a solnSpace object
setGeneric("optCF", signature = c("dataSet", "geneSig"),
	function(dataSet, geneSig, iter, indexTrain, indexTest, disc=c(0.005, 0.01, 0.03, 0.05), 
		MFS="MFS", met="met", optMeth="Nelder-Mead", usePDF=FALSE, jPDF=FALSE) standardGeneric("optCF"))
		
setMethod("optCF", c(dataSet="ExpressionSet", geneSig="geneSignature"),
	function(dataSet, geneSig, iter, indexTrain, indexTest, disc, MFS, met, optMeth, usePDF, jPDF){
		result <- new("solnSpace")
		nGenes <- getNGenes(geneSig)
		maxIndex <- ncol(exprs(dataSet))
		thresholds <- matrix(NA, nrow=iter, ncol=nGenes)
		colnames(thresholds) <- getGenes(geneSig)
		identTrain <- matrix(NA, nrow=iter, ncol=maxIndex)
		identTest <- matrix(NA, nrow=iter, ncol=maxIndex)
		statsTrain <- matrix(NA, nrow=iter, ncol=2)
		statsTest <- matrix(NA, nrow=iter, ncol=2)
		colnames(statsTrain) <- c("pval", "size")		
		colnames(statsTest) <- c("pval", "size")

		cf <- function(x) ensembleCostFcn(dataSet, 
			setThresholds(geneSig, x), jPDF, disc, indexTrain, MFS, met)
		
		for(j in 1:iter){
			initparms <- rnorm(nGenes, sd=1)
	 		fit <- optim(initparms, cf, method=optMeth)
			parms <- fit$par
			geneSig <- setThresholds(geneSig, parms)

			positiveTrain <- ensembleAdjustable(dataSet, geneSig, indexTrain)
			pvalTrain <- applySigLR(dataSet, geneSig, MFS, met, indexTrain)
			sizeTrain <- applySigSize(dataSet, geneSig, TRUE, indexTrain)
			
			positiveTest <- ensembleAdjustable(dataSet, geneSig, indexTest)
			pvalTest <- applySigLR(dataSet, geneSig, MFS, met, indexTest)
			sizeTest <- applySigSize(dataSet, geneSig, TRUE, indexTest)
			
			thresholds[j,] <- parms
			identTrain[j,] <- 1:maxIndex %in% indexTrain 
			identTest[j,] <- 1:maxIndex %in% indexTest
			statsTrain[j,] <- c(pvalTrain, sizeTrain)
			statsTest[j,] <- c(pvalTest, sizeTest)
		}
		
		result <- setSolnSpace(result, thresholds, statsTrain, statsTest, identTrain, identTest)
		
		return(result)
	})	
	
	
##function to generate the solution space objects using an optimizer options from R's optim fcn
##this function utilizes parallel processing, and by default uses k-fold cv
##parallelization is utilized by splitting up each k into iterPerK/nCores chunks
setGeneric("analysisPipeline", signature = c("dataSet", "geneSig"),
	function(dataSet, geneSig, iterPerK=2500, k=3, rand=TRUE, newjpdf=FALSE, jpdf=FALSE, nJPDF=12500, 
		disc=c(0.005, 0.01, 0.03, 0.05), MFS="MFS", met="met", optMeth="Nelder-Mead") standardGeneric("analysisPipeline"))
	
setMethod("analysisPipeline", c(dataSet="ExpressionSet", geneSig="geneSignature"),
	function(dataSet, geneSig, iterPerK, k, rand, newjpdf, jpdf, nJPDF, disc, MFS, met, optMeth){
		result <- list()
		mc <- FALSE
		
		nSamples <- ncol(exprs(dataSet))
		allIndices <- 1:nSamples
		indices <- allIndices
		if(rand) indices <- sample(1:ncol(dataSet), replace=F)
		pdf <- jpdf
		if(newjpdf) pdf <- eJPDF(dataSet, geneSig, nJPDF, allIndices, MFS, met)
		
		if(!mc){
			for(i in 0:(k-1)){
				kmeansIndices <- ((round(i*nSamples)/k)+1):(round((i+1)*nSamples/k))
				kmeansIndices <- indices[kmeansIndices]
				testIndices <- indices[!indices %in% kmeansIndices]
				
				resultPerK <- optCF(dataSet, geneSig, iterPerK, kmeansIndices, testIndices, 
					disc, MFS, met, optMeth, jpdf, pdf)
				
				result[[i+1]] <- resultPerK
			}
		}

		else{
			processes <- list()
			indexListTrain <- list()
			indexListTest <- list()
			
			iterPerKPerCore <- floor(iterPerK/nCores)
			iterLeftover <- iterPerK - nCores*iterPerKPerCore
			
			for(i in 0:(k-1)){
				kmeansIndices <- ((round(i*nSamples)/k)+1):(round((i+1)*nSamples/k))
				kmeansIndices <- indices[kmeansIndices]
				testIndices <- indices[!indices %in% kmeansIndices]
				
				indexListTrain[[i+1]] <- kmeansIndices
				indexListTest[[i+1]] <- testIndices
				
				for(j in 1:(nCores-1)){
					processes[[j]] <- mcparallel(optCF(dataSet, geneSig,
						iterPerKPerCore, indexListTrain[[i+1]], indexListTest[[i+1]],
						disc, MFS, met, optMeth, jpdf, pdf))
				}
				
				processes[[nCores]] <- mcparallel(optCF(dataSet, geneSig,
						iterPerKPerCore+iterLeftover, indexListTrain[[i+1]], indexListTest[[i+1]],
						disc, MFS, met, optMeth, jpdf, pdf))
				
				result <- c(result, mccollect(processes))
			}
						
			result <- unlist(result)
		}

		res <- result[[1]]

		for(i in 2:length(result)){
			res <- combineSolnSpace(res, result[[i]])
		}

		res <- summarizeCVCuts(res, geneSig)
		
		return(res)
	})
	
	
#This function takes in solnSpace results and selects solutions that produce 
#significance in both training and cv sets
setGeneric("getCVCuts", signature="cutoffResults",
	function(cutoffResults) standardGeneric("getCVCuts"))

setMethod("getCVCuts", c(cutoffResults = "solnSpace"),
	function(cutoffResults){
		statsTrain <- getStatsTrain(cutoffResults)
		statsTest <- getStatsTest(cutoffResults)
	
		trainPValues <- statsTrain[,1]
		testPValues <- statsTest[,1]

		p <- 0.05
		sigCV <- c()
		while(length(sigCV) <= 1){
			sigTrain <- which(trainPValues < p)
			sigTest <- which(testPValues < p)

			sigCV <- intersect(sigTrain, sigTest)
			if(length(sigCV) == 0) {
				message(paste("No solutions are significant for both training and testing sets, p=", p, sep=""))
				message("\nIncreasing threshold by 0.05")
			}
			p <- p + 0.05
		}

		return(subsetSolnSpace(cutoffResults, sigCV))
	})
	
#This function takes in solnSpace results and selects solutions that produce 
#significance in both training and cv sets, and summarizes those solns, returning
#a geneSignature object
setGeneric("summarizeCVCuts", signature="cutoffResults",
	function(cutoffResults, geneSig) standardGeneric("summarizeCVCuts"))

setMethod("summarizeCVCuts", c(cutoffResults = "solnSpace"),
	function(cutoffResults, geneSig){
		th <- getThresholds(getCVCuts(cutoffResults))
		th <- colMeans(th)
		
		result <- geneSig
		result <- setThresholds(result, th)
		return(result)	
	})

#this function will generate, for a gene signature of N genes, N solnSpace
#objects, each solnSpace object leaving out one gene product
#this function is pseudo-parallelized, as the important interactions happen 
#when calling analysisPipeline, and analysisPipeline IS parallelized
setGeneric("analysisLOOGen", signature=c("dataSet", "geneSig"),
	function(dataSet, geneSig, iterPerK=2500, k=3, rand=TRUE,  jpdf=NULL, newjpdf=TRUE, njpdf=12500) standardGeneric("analysisLOOGen"))
	
setMethod("analysisLOOGen", c(dataSet = "ExpressionSet", geneSig = "geneSignature"),
	function(dataSet, geneSig, iterPerK, k, rand, jpdf, newjpdf, njpdf){
		result <- list()
		
		if(is.null(jpdf) && newjpdf) jpdf <- eJPDF(dataSet, geneSig, njpdf)
		
		for(i in 1:length(geneSig@geneSet)){
			gs <- geneSig
			gs@geneSet <- gs@geneSet[-i]
			gs@geneDirect <- gs@geneDirect[-i]

			result[[i]] <- analysisPipeline(dataSet, gs, iterPerK, k, rand, jpdf=jpdf) 
		}
		
		return(result)
	})

	
#this function generates the bpms signature
setGeneric("genBPMSSig", signature="dataSet",
	function(dataSet) standardGeneric("genBPMSSig"))
	
setMethod("genBPMSSig", c(dataSet = "ExpressionSet"),
	function(dataSet){
	
	})
