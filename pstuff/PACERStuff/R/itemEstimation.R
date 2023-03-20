

############################################################ Helper Functions ############################################################
 
irtProbs <- function(D, A, betas, C, machineValue){
	Da <- D * A
	pr <- C + (1 - C)/(1 + exp(Da * betas))
	ind1 <- which(pr > 1 - machineValue, arr.ind = TRUE)
	ind2 <- which(pr < machineValue, arr.ind = TRUE)
    if (nrow(ind1 > 0)) 
        pr[ind1] <- 1 - machineValue
    if (nrow(ind2 > 0)) 
        pr[ind2] <- machineValue
    pr
}

### Function to find parameters of beta distribution
### given an expected mean and standard deviation 
### findBeta(.2, sqrt(.0076)) is beta(4,16)

findBeta <- function(mu, sigma){
	a <- ((1-mu)/sigma^2 - 1/mu)*mu^2
	b <- a * (1/mu - 1)
	list(a = a, b = b)
}

fixItemConstraints <- function(model, position.a = NULL, position.b = NULL, position.c = NULL, fix.a.value = NULL, fix.b.value = NULL, fix.c.value = NULL, L){
	if(model == '1pl'){
		b.lower <- rep(-Inf, L)
		b.upper <- rep(Inf, L)		 
		if(!is.null(fix.b.value)) b.lower[position.b] <- b.upper[position.b] <- fix.b.value
		lower <- b.lower
		upper <- b.upper		
	}
	if(model == '2pl'){
		a.lower <- b.lower <- rep(-Inf, L)
		a.upper <- b.upper <- rep(Inf, L)		 
		if(!is.null(fix.a.value)) a.lower[position.a] <- a.upper[position.a] <- fix.a.value
		if(!is.null(fix.b.value)) b.lower[position.b] <- b.upper[position.b] <- fix.b.value
		lower <- c(a.lower, b.lower)
		upper <- c(a.upper, b.upper)		
	}
	if(model == '3pl'){
		a.lower <- b.lower <- c.lower <- rep(-Inf, L)
		a.upper <- b.upper <- c.upper <- rep(Inf, L)		 
		if(!is.null(fix.a.value)) a.lower[position.a] <- a.upper[position.a] <- fix.a.value
		if(!is.null(fix.b.value)) b.lower[position.b] <- b.upper[position.b] <- fix.b.value
		if(!is.null(fix.c.value)) c.lower[position.c] <- c.upper[position.c] <- fix.c.value
		lower <- c(a.lower, b.lower, c.lower)
		upper <- c(a.upper, b.upper, c.upper)		
	}
	list(lower = lower, upper = upper)
}

bigGPCM <- function(theta, d, score, a, D = 1.7){
	L <- length(score)
	Da <- D * a
	#numer <- score * theta - sapply(1:L, function(i) sum(d[1:score[i]]))
	numer <- score * theta - cumsum(d)[score]
	denom <- sum(exp(cumsum(Da * (theta - d))))
	exp(Da * (numer))/denom
}

bigGRM <- function(theta, d, score, a, D = 1.7, machineValue){
	maxD <- length(d)
	pos1 <- which(score == 0)
	pos2 <- which(score == maxD)
	pos3 <- which(score > 0 & score < maxD)

	pr.1 <- 1/(1 + exp(D*a*(theta - d[(score[pos1]+1)])))
	pr.2 <- 1/(1 + exp(-D*a*(theta - d[score[pos2]])))
	pr.3 <- 1/(1 + exp(D*a*(theta - d[(score[pos3]+1)]))) - 1/(1 + exp(D*a*(theta - d[score[pos3]])))

	### Now put them back together in order
	pr <- numeric(length(score))
	pr[pos1] <- pr.1
	pr[pos2] <- pr.2
	pr[pos3] <- pr.3
	ind1 <- which(pr > 1 - machineValue)
	ind2 <- which(pr < machineValue)
    if (any(ind1)) 
        pr[ind1] <- 1 - machineValue
    if (any(ind2))  
        pr[ind2] <- machineValue
    pr
}

############################################# Function #############################################
### Note: This version uses the Rcpp compiled version to compute the marginal likelihood
### To run the pure R version uncomment the lines below and comment out the irt_like function

irt <- function(dat, model = c('1pl', '2pl', '3pl'), ind.poly = NULL, p.model = c('gpcm', 'pcm'), a.prior = FALSE, b.prior = FALSE, c.prior = TRUE, control = list()){
	con <- list(Q = 21, D = 1, pop.mu = 0, pop.sigma = 1, startValues = NULL, tol = 1e-03, max.iter = 500, aPriorMean = 0, aPriorSD = 1, 
		bPriorMean = 0, bPriorSD = 1, c.mean = .2, c.sd = .087, missing = c('sparse', 'complete'), which.optim = c('nlminb','bobyqa'),
		position.a = NULL, position.b = NULL, position.c = NULL, fix.a.value = NULL, fix.b.value = NULL, fix.c.value = NULL, common.c = FALSE, get.sems = FALSE)
	con[names(control)] <- control
	qq <- gauss.quad.prob(con$Q, dist = 'normal', con$pop.mu, con$pop.sigma)
	nds <- qq$nodes
	wts <- qq$weights
	if(con$missing == 'complete') dat <- dat[complete.cases(dat), ]
	if(con$missing == 'sparse') dat[which(is.na(dat), arr.ind = TRUE)] <- 0
	if(length(ind.poly) > 0) {		## Split the data if necessary
		X <- dat[, -ind.poly]
		X2 <- dat[, ind.poly]
		X2 <- as.matrix(X2)
	} else {
		X <- dat
	}
	X <- as.matrix(X)
	mX <- 1 - X
	L <- ncol(X)
	
	mod.num <- ifelse(model == '1pl', 1, ifelse(model == '2pl', 2, 3))
	num.dich.params <- L * mod.num
	
	if(length(ind.poly) > 0) {
		X2 <- X2 + 1 #Account for way GPCM is parameterized
		S <- apply(X2,2,max) ## Max score by item
		M <- rep(1:length(S), S-1)
		T <- length(M) ## total number of steps to estimate
		I <- length(S) ### number of items
		G <- sum(S)
	}	
	
	if(c.prior) betaParms <- findBeta(con$c.mean, con$c.sd)
	if(is.null(con$startValues)) {		 
		if(model == '3pl')con$startValues <- c(rep(1,L*2), rep(con$c.mean,L))
		if(model == '2pl')con$startValues <- c(rep(1,L), rep(0,L))
		if(model == '1pl')con$startValues <- rep(0,L)
	}
	if(length(ind.poly) > 0) con$startValues <- c(con$startValues, c(rep(1,I), sort(rnorm(T)))) ### Always put poly part at end
	
	machineValue <- sqrt(.Machine$double.eps) ### recycle this value	
	fn <- function(params){
		if(model == '3pl'){
			a <- params[1:L]            
			b <- params[(L+1):(L*2)]        
			if(con$common.c) {
				c <- params[(L*2+1)]
			} else {
				c <- params[(L*2+1):(L*3)]
			}			
		}
		if(model == '2pl'){
			a <- params[1:L]            
			b <- params[(L+1):(L*2)]        
			c <- rep(0, L)
		}
		if(model == '1pl'){
			b <- params[1:L]            
			a <- rep(1,L)        
			c <- rep(0, L)
		}
		if(length(ind.poly) > 0){
			if(p.model == 'gpcm') a.params <- params[(num.dich.params+1):(I+num.dich.params)]            
			if(p.model == 'pcm') a.params <- rep(1, length(ind.poly))
			d.params <- params[(I+num.dich.params+1):(num.dich.params+G)]
			d.params <- split(d.params, M)
			d.params <- lapply(d.params, function(i) c(0,i)) ### prepend the 0	
		}		
		
		if(length(ind.poly) > 0){
			poly.pr <- exp(sapply(1:length(nds), function(j) rowSums(log(sapply(1:I, function(i) bigGPCM(nds[j],d.params[[i]],X2[,i],a.params[i]))))))
			poly.res <- poly.pr %*% wts			
		}
		
		poly.total <- if(length(ind.poly) > 0) {
			sum(log(poly.res))
		} else 0		
		
		A <- matrix(rep(a, con$Q), nrow = con$Q, byrow = TRUE)
		betas <- matrix(rep(b, con$Q), nrow = con$Q, byrow = TRUE) - nds
		if(con$common.c & model == '3pl'){
			C <- matrix(rep(rep(c, L), con$Q), nrow = con$Q, byrow = TRUE)
		} else {
			C <- matrix(rep(c, con$Q), nrow = con$Q, byrow = TRUE)
		}
		pr <- irtProbs(con$D, A, betas, C, machineValue)
		pr.t <- t(pr)
		a.prior.values <- if(a.prior & model %in% c('2pl', '3pl')) {
			dlnorm(a, meanlog = con$aPriorMean, sdlog = con$aPriorSD, log = TRUE)
		} else 0
		b.prior.values <- if(b.prior) {
			dnorm(b, con$bPriorMean, con$bPriorSD, log = TRUE)
		} else 0
		c.prior.values <- if(c.prior & model == '3pl') {
			pos <- which(!is.infinite(dbeta(c, betaParms$a, betaParms$b, log = TRUE))) ### This is needed in case user anchors a c-parameter at 0. It removes log(0) = Inf
			dbeta(c, betaParms$a, betaParms$b, log = TRUE)[pos]
		} else 0
		res <- exp(X %*% log(pr.t) + mX %*% log(1 - pr.t)) %*% wts
		#res <- irt_like(X, mX, pr, wts)
		total <- sum(log(res)) + sum(a.prior.values) + sum(b.prior.values) + sum(c.prior.values) + poly.total		
		-total		
	}  
	start <- Sys.time()		
	constraints <- fixItemConstraints(model, con$position.a, con$position.b, con$position.c, con$fix.a.value, con$fix.b.value, con$fix.c.value, L = L)
	if(length(ind.poly) > 0){
		constraints$lower <- c(constraints$lower, rep(-Inf, (T+I)))
		constraints$upper <- c(constraints$upper, rep(Inf, (T+I)))
	}
	if(con$which.optim == 'nlminb') res <- nlminb(con$startValues, fn, control = list(trace=0, rel.tol = con$tol, iter.max = con$max.iter), lower = constraints$lower, upper = constraints$upper)		
	if(con$which.optim == 'bobyqa') res <- bobyqa(con$startValues, fn, control = list(xtol_rel = con$tol, maxeval = con$max.iter))
	end <- Sys.time()
	timer <- end - start
	tmp1 <- data.frame(MODEL = rep('IRT3pl', L), ScorePoints = rep(1,L))
	tmp2 <- data.frame(matrix(NA, L,8)) #augmented columns for test scoring to work
	colnames(tmp2) <- paste('P', 3:10, sep='')
	if(model == '3pl' & con$common.c) coefs <- data.frame(tmp1,P0 = res$par[1:L], P1 = res$par[(L+1):(L*2)], P2 = rep(res$par[(L*2+1)], L),tmp2)
	if(model == '3pl' & con$common.c == FALSE) coefs <- data.frame(tmp1,P0 = res$par[1:L], P1 = res$par[(L+1):(L*2)], P2 = res$par[(L*2+1):(L*3)], tmp2)
	if(model == '2pl') coefs <- data.frame(tmp1, P0 = res$par[1:L], P1 = res$par[(L+1):(L*2)], P2 = 0, tmp2)
	if(model == '1pl') coefs <- data.frame(tmp1, P0 = 1, P1 = res$par[1:L], P2 = 0, tmp2)
	
	if(length(ind.poly) > 0){
		out <- split(res$par[(I+num.dich.params+1):(num.dich.params+G)], M)
		ScorePoints <- sapply(out, length)
		steps <- plyr::ldply(out, rbind)[,-1]
		colnames(steps) <- paste('P', 1:ncol(steps), sep='')
		a.params <- res$par[(num.dich.params+1):(I+num.dich.params)]
		tmp3 <- data.frame(matrix(NA, length(ind.poly), 10 - ncol(steps))) #augmented columns for test scoring to work
		colnames(tmp3) <- paste('P', (ncol(steps)+1):10, sep = '')
		poly.coefs <- data.frame(MODEL = 'IRTgpc', ScorePoints = ScorePoints, P0 = a.params,steps, tmp3)
		coefs <- rbind(coefs, poly.coefs)
	}
	
	if(con$get.sems){
		SE <- if(model == '3pl' & con$common.c){
		P <- L * 2 + 1
		A <- ncol(dat)*3-P
		r1 <- sqrt(diag(solve(hessian(fn, res$par)[1:P, 1:P])))
		c(r1, rep(r1[P], A))
	} else {
		theHessian <- hessian(fn, res$par)
		if(all(eigen(theHessian)$val > machineValue)){
			round(sqrt(diag(solve(theHessian))), 3)
		} else {
			rep(0, length(res$par))
		}
	}
	} else {
		SE.mat <- NULL
	}
	
	### Format the SE for printing in the UI
	if(con$get.sems){
		se.out <- data.frame(matrix(round(SE[1:num.dich.params], 3), L, mod.num))
		colnames(se.out) <- paste('P', 0:(mod.num-1), sep='') 
		tmp <- data.frame(matrix(NA, L, 10 - (mod.num-1)))
		colnames(tmp) <- paste('P', ncol(se.out):10, sep = '')
		SE.mat <- data.frame(se.out, tmp)

		if(length(ind.poly) > 0){
			out <- split(SE[(I+num.dich.params+1):(num.dich.params+G)], M)
			steps <- plyr::ldply(out, rbind)[,-1]
			colnames(steps) <- paste('P', 1:ncol(steps), sep='')
			a.params <- SE[(num.dich.params+1):(I+num.dich.params)]
			tmp3 <- data.frame(matrix(NA, length(ind.poly), 10 - ncol(steps))) #augmented columns for test scoring to work
			colnames(tmp3) <- paste('P', (ncol(steps)+1):10, sep = '')
			poly.se <- data.frame(P0 = a.params,steps, tmp3)
			SE.mat <- rbind(SE.mat, poly.se)
		}
	}	
	list('coefficients' = coefs, time = timer, logLike = res$objective * -1, iterations = res$iterations, convergence = res$convergence, SE = SE.mat)
}

export4TestScoring <- function(dat){
	dat$testID <- 1:nrow(dat)
	mdata <- melt(data.table(dat), id=c("testID")) 
	mdata <- mdata[order(mdata$testID),]
	names(mdata) <- c('testID', 'key', 'score')
	pos <- which(is.na(mdata$score))
	if(length(pos) > 0) mdata <- mdata[-pos,]
	out <- mdata[, c('testID', 'key', 'score')]
	return(out)
}




