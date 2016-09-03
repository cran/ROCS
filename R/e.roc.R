e.roc <- function(x, mu, method="RNA", bt.ci=TRUE, bt.nreps=100, do.plot=TRUE)
{
    ######### Larger scores tend to be positvie #########
    nLen0 <- length(x)  
    nLen1 <- length(mu)
    if(nLen0!=nLen1)	stop("Error in input arguments!")
      
    x0 <- x[mu<=0.5]  
    x1 <- x[mu>0.5]
    if(median(x0)>median(x1)){
		cat("Switch the Positive Class Label!")
		mu <- rep(1, nLen1) - mu
    }
    xObs <- x  
    yProb <- mu
    
	get.metrics <- function(xObs, yProb)
	{
		obs.num <- length(xObs)  ## sample size
		zProb <- 1 - yProb		 ## probabilities of being negative
		
		## Generating v vector
    	tVec <- rep(0, obs.num)
    	vec0 <- rep(0, obs.num )  
    	vec1 <- rep(0, obs.num )
    	for(iter in 0:(obs.num-2))	tVec[iter+1] <- 1/((iter+1)*(iter+2))
    	tVec[obs.num] <- 1/obs.num
    	for(i in 1:obs.num ){
      		val <- 0:(obs.num-1)
      		n.vec <- zProb[-i] 
      		p.vec <- yProb[-i]
      		prob.n <- ppoibin(val, pp=n.vec, method=method)
      		prob.p  <- ppoibin(val, pp=p.vec, method=method)
      		vec0[i] <- sum( prob.n*tVec )
      		vec1[i] <- sum( prob.p*tVec )
    	}
		## sort the samples
    	xVal <- unique(xObs)  
    	pLen <- length(xVal)
    	xSort <- sort(xVal)
    	## compute eFPR, eTPR, eTDR
    	get.fptptdr <- function(x,y,z)
    	{
      		l <- length(y)
      		new.y <- y*vec1
      		new.z <- z*vec0
      		l1<-sum(new.y)
      		l0<-sum(new.z)

      		bb<- NULL  #count of true positive being claimed positive
      		aa<- NULL  #count of true negative claimed positive
      		cc<- NULL 
      		for(i in 1:pLen){
	  			tpr <- sum( new.y[x>=xSort[i] ] )
	  			fpr <- sum( new.z[x>=xSort[i] ] )
	  			tdr <- sum( y[x>=xSort[i] ] )/sum(x>=xSort[i])
	  			if(tpr>1)	 tpr <- 1.0 
	  			if(fpr>1)	 fpr <- 1.0
        		bb <- c(bb, tpr )
        		aa <- c(aa, fpr )
	  			cc <- c(cc, tdr )
      		}
      		FP<-aa
      		TP<-bb
      		TDR<-cc
      		FP<-c(1, FP, 0)
      		TP<-c(1, TP, 0)
      		TDR<-c(min(TDR), TDR, 1)
      		return(list(FP=FP, TP=TP, TDR=TDR) )
    	}
		## compute eAUC
		get.eAuc <- function(x,y,z)
    	## x are values ordered from smallest to largest, y are group labels
    	{
      		FP.TP<-get.fptptdr(x,y,z)
      		FP<-FP.TP$FP
      		TP<-FP.TP$TP
      		l<-length(FP)

      		abs.d.TP<- -diff(TP)
      		mid.FP<- FP[1:(l-1)]+diff(FP)/2
      		auc<-sum(abs.d.TP*(1-mid.FP) )
      		auc
    	}
    	eAuc <- round(get.eAuc(xObs,yProb,zProb), digits=3)
    	fptptdr <- get.fptptdr(xObs,yProb,zProb)
		## return volume under surface
    	list(eAUC=eAuc, eTPR=fptptdr$TP, eFPR=fptptdr$FP, eTDR=fptptdr$TDR)
	}
	
	expected.acu.values <- get.metrics(xObs, yProb)
	expected.auc.ci <- rep(0, 2)
	if(bt.ci){
		totalObs <- length(xObs)
		alpha <- 0.05
		expected.auc <- NULL
		for(it in 1:bt.nreps){
  			indexVector <- sample(1:totalObs, size=totalObs, replace=TRUE)
  			dataVector <- xObs[indexVector]
  			probVector <- yProb[indexVector]
  			e.roc.bt <- get.metrics(xObs=dataVector, yProb=probVector)
  			expected.auc <- c(expected.auc, e.roc.bt$eAUC)
  		}
  		expected.auc <- sort(expected.auc)
  		expected.auc.ci <- c(expected.auc[floor(bt.nreps*alpha*0.5)], expected.auc[floor(bt.nreps*(1-alpha*0.5))]) 	 
	}
	## generate plots
    if(do.plot){
    	plot(expected.acu.values$eFPR, expected.acu.values$eTPR, type="l", lwd=2, col="red", xlab="eFPR",ylab="eTPR",main="eROC Curve")
    }

    list(eAUC=expected.acu.values$eAUC, eTPR=expected.acu.values$eTPR, eFPR=expected.acu.values$eFPR, eTDR=expected.acu.values$eTDR, CI=expected.auc.ci)
}