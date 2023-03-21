# irt.plot() accepts a list of paramList objects and outputs xyplots for TCC, TIF, or SEM.
irt.plot <- function(params, type = c('TCC', 'TIF', 'SEM'), cuts = NULL, proportional_tcc = FALSE,
	show_extremum = FALSE, control = list(), ...) {

	type <- toupper(type)
	type <- match.arg(type)
	K <- length(params)
	LegendText <- paste('Form', 1:K, sep = ' ')
	con <- list(D = 1.7, A = 1, B = 0, lower = -8, upper = 8, mu = 0, sigma = 1)
	con[names(control)] <- control
	thetas <- seq(con$lower, con$upper, by = .1) # Domain of plot.

	groups <- rep(1:K, each = length(thetas)) # To distinguish multiple series being plotted.
	#trellis.par.set(fontsize = list(text = 17)) # Default text size is 12.

   if(type == 'TCC') {

      result <- lapply(1:K, function(i) tcc(params[[i]], thetas = thetas, control = con))

      # combine TCC object into one.
      tcc_merge <- structure(list('value' = do.call(c, lapply(result, function(x) x$value)),
         'thetas' = rep(thetas, K), 'group' = groups), class = 'TCC')

      return(xyplot(tcc_merge, proportional_tcc = proportional_tcc, cuts = cuts))

   } else if(type == 'TIF') {
		info_class <- 'TIF'
		result <- lapply(1:K, function(i) testinfo(params[[i]], thetas = thetas, return_value = "INFO", control = con))

	} else if(type == 'SEM') {
		info_class <- c('TIF', 'SEM')
		result <- lapply(1:K, function(i) testinfo(params[[i]], thetas = thetas, return_value = "SE", control = con))
	}
	# combine TIF objects into one.
   tif_merge <- structure(list('value' = do.call(c, lapply(result, function(x) x$value)),
      'thetas' = rep(thetas, K), 'extremum' = do.call(c, lapply(result, function(x) x$extremum)),
      'theta_extremum' = do.call(c, lapply(result, function(x) x$theta_extremum)),
      'group' = groups), class = info_class)

   return(xyplot(tif_merge, show_extremum = show_extremum, cuts = cuts))
}

scoreConversionTable <- function(params, method = c("MLE", "MAP", "EAP", "TCC"), model = c('GPCM', 'GRM'), control = list(), ...){

	method <- toupper(method)
    method <- match.arg(method)
	model <- match.arg(model)

	con <- list(D = 1.7, A = 1, B = 0, lower = -6, upper = 6, mu = 0, sigma = 1)
    con[names(control)] <- control

	if(method == 'TCC'){
		max <- sum(params$score.pts)
		raw <- seq(from = 0, to = max)
		theta <- sapply(1:length(raw), function(i) score.tcc(raw[i], params, model = model, D = con$D, lower = con$lower, upper = con$upper))
		se <- sapply(1:length(raw), function(i) testInfo(theta[i], params, D = con$D))
	} else {

		tot.possible <- params$score.pts
		max.points <- tot.possible
		theta <- numeric(sum(max.points))
		se <- numeric(sum(max.points))

		tmp <- matrix(0, nrow = length(tot.possible), ncol = length(tot.possible))
		### First just subtract 1
		for(i in 1:length(tot.possible)){
			tot.possible[i] <- 	tot.possible[i] - 1
			tmp[,i] <- tot.possible
		}

		### Now, iteratively subtract 1
		while(any(tot.possible > 0)){
			pos <- which(tot.possible > 0)[1]
			tot.possible[pos] <- tot.possible[pos] - 1
			tmp <- cbind(tmp, tot.possible)
		}

		### Finally bind the perfect score
		tmp <- cbind(tmp, max.points)
		tmp
		raw <- colSums(tmp)

		### Now do scoring
		for(i in 1:ncol(tmp)){
			score <- irt.ability(as.vector(tmp[,i]), params = params, ind.dichot = params$mcPos,
				std.err = 'Hessian', model = model, method, control = list(D=con$D, mu = con$mu, sigma = con$sigma), ...)
			theta[i] <- coef(score)
			se[i] <- sqrt(vcov(score))
		}
	}
	results <- data.frame(theta = theta, se = se, raw.score = raw, scale.score = round(theta * con$A + con$B))
	results <- results[order(results$raw.score) , ]
	results
}

testPlotly <- function(dat, myTitle, ylab){
	dat <- dat
	K <- ncol(dat)-1	
	#values <- do.call(c, c(dat[, 2:ncol(dat)]))	
	if(K==1){
		values <- dat[, 2]
	} else {		
		values <- do.call(c, c(dat[, 2:ncol(dat)]))
	}
	thetas <- as.vector(rep(dat$thetas,K))	
	ind <- gl(K, nrow(dat), label = paste('Form', 1:K, sep=' '))
	tmp <- data.frame(theta = as.numeric(thetas), values = as.numeric(values), ind = ind)	

	fig <- plot_ly(
			data = tmp,
			x = ~theta, 
			y = ~values,
			color = ~ind, 
			colors = "Set2",
			type = "scatter",
			mode = "lines"
		)	
	fig <- fig %>% layout(title = myTitle,
		 xaxis = list(title = "Theta"),
		 yaxis = list (title = ylab))
	fig
}

### Function for building test form
simpleFormMax <- function(dat, theta, Nitems, D = 1.7){
	infoValue <- numeric(nrow(dat))
	nms <- paste('p', 1:10, sep='')
	vars <- which(colnames(dat) %in% nms)
	for(i in 1:nrow(dat)){
		if(dat$scorepoints[i]==1) {
			infoValue[i] <- info.3pl(theta, dat[i, 'p0'], dat[i, 'p1'], dat[i, 'p2'], D = D) 
		} else {
			steps <- dat[i,][which(colnames(dat[i,]) %in% nms)]
			steps <- steps[which(!is.na(steps))]
			infoValue[i] <- info.gpcm(theta, dat[i, 'p0'], steps, D = D)
		}
	}
	df <- data.frame(dat$itemid, infoValue)
	df <- df[order(df[,2], decreasing = TRUE), ]
	df[1:Nitems, 1]
}


autoTest <- function(dat, blueprint, theta, constraintMat, D = 1.7){

	blueprint <- as.data.frame(blueprint)

	### First compute item information for all items in the bank
	infoValue <- numeric(nrow(dat))
	nms <- paste('p', 1:10, sep='')
	vars <- which(colnames(dat) %in% nms)
	for(i in 1:nrow(dat)){
		if(dat$scorepoints[i]==1) {
			infoValue[i] <- info.3pl(theta, dat[i, 'p0'], dat[i, 'p1'], dat[i, 'p2'], D = D) 
		} else {
			steps <- dat[i,][which(colnames(dat[i,]) %in% nms)]
			steps <- steps[which(!is.na(steps))]
			infoValue[i] <- info.gpcm(theta, dat[i, 'p0'], steps, D = D)
		}
	}
	
	f.obj <- infoValue
	
	### Define constraints
	K <- nrow(dat)
	tot.con <- matrix(rep(1,K), nrow = 1, byrow = TRUE)  # total number of items; 
	
	### One-hot encoding for categories
	X.list <- vector(mode = "list", length = ncol(blueprint))
	for(i in 1:ncol(blueprint)){
		X.list[[i]] <- t(model.matrix(~0+blueprint[,i]))
	}
	X.strand <- do.call(rbind,X.list) ### Constraints for the BP inputs
	
	### Put all constraints together into big matrix
	f.con <- rbind(tot.con , X.strand, diag(K)) ### diag(K) makes sure each item is used only once
	
	### How to treat constraints
	f.dir <- c(constraintMat$direction, rep("<=",K))
	
	f.rhs <- c(constraintMat$value,    ## number of items; 
           rep(1,K)    ## each item can be used at most once
           ) 
	
	result <- lp("max", f.obj, f.con, f.dir, f.rhs, all.bin=TRUE)
	dat[which(result$solution==1), 'itemid']
}


