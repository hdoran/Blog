# This R script defines various math functions used throughout the PADI environment.

## Returns expected score value for a gpcm item
poly.exp.gpcm <- function(theta, D, a, d) {
	gpcm <- function (theta, score) {
	  Da <- D * a
	  exp(Da * (score * theta - sum(d[1:floor(score)])))/
			sum(exp(cumsum(Da * (theta - d))))
	}
	sum(sapply(2:length(d), function(i) gpcm(theta, i) * (i-1)))
}

## Returns a vector of probabilities for each possible scorepoint for gpcm items
poly.pr.gpcm <- function(theta, D, a, d) {
	gpcm <- function (theta, score) {
	  Da <- D * a
	  exp(Da * (score * theta - sum(d[1:floor(score)])))/
			sum(exp(cumsum(Da * (theta - d))))
	}
	return( sapply(2:length(d), function(i) gpcm(theta, i)) )
}

## Returns expected score value for grm item
poly.exp.grm <- function(theta, D, a, d) {
   grm <- function (theta, score) {
	   maxD <- length(d)
	   if(score == 0) {
			pr <- 1/(1 + exp(D*a*(theta - d[(score+1)])))
	   } else if(score == maxD) {
			pr <- 1/(1 + exp(-D*a*(theta - d[score])))
	   } else {
			pr <- 1/(1 + exp(D*a*(theta - d[(score+1)]))) - 1/(1 + exp(D*a*(theta - d[score])))
	   }
	   pr
   }
   sum(sapply(1:length(d), function(i) grm(theta, i) * i))
}

# Returns a vector of probabilities for each possible scorepoint for grm items
poly.pr.grm <- function(theta, D, a, d) {
	grm <- function (theta, score) {
		maxD <- length(d)
		if(score == 0) {
			pr <- 1/(1 + exp(D*a*(theta - d[(score+1)])))
		} else if(score == maxD) {
			pr <- 1/(1 + exp(-D*a*(theta - d[score])))
		} else {
			pr <- 1/(1 + exp(D*a*(theta - d[(score+1)]))) - 1/(1 + exp(D*a*(theta - d[score])))
		}
		pr
	}
	return( sapply(1:length(d), function(i) grm(theta, i)) )
}

## Returns probability (which is also expected score) that a 3pl item is answered correctory (scorepoint 1)
pr.3pl <- function(theta, D, a_uni, b_uni, c_uni) {
	return( c_uni + (1 - c_uni) / (1 + exp(-D * a_uni * (theta - b_uni))) )
}

# 3pl probability give some theta and nuisance parameter u (for cluster assertions)
pr.latentu <- function(a_assert, b_assert, c_assert, u_nuisance, D, theta) c_assert + (1 - c_assert) / (1 + exp(-D * a_assert * (theta + u_nuisance - b_assert)))

## Pr for each assertion within cluster item, at each theta
margprob4cluster <- function(params, thetas, control = list()) {

   con = list(D = 1.0, Q = 50, sigma = 1)
   con[names(control)] <- control

   gq <- gauss.quad.prob(con$Q, dist = 'normal', sigma = con$sigma)
   nodes <- gq$nodes * sqrt(params$variance)
   whts <- gq$weights
   a_par <- params$`3pl`$a
   b_par <- params$`3pl`$b
   c_par <- params$`3pl`$c

   # Loop over 1) Each assertion, 2) Each theta (ability estimate) of interest, and then 3) each quadrature node. Matrix returned, thetas BY assertions
   probs <- sapply(1:length(params$items), function(k)
      sapply(thetas, function(x) sum(
         sapply(1:con$Q, function(m) pr.latentu(a_assert = a_par[k], b_assert = b_par[k], c_assert = c_par[k], u_nuisance = nodes[m], D = con$D, theta = x) * whts[m]))))

	probs <- rbind(probs)
	rownames(probs) <- formatC(thetas) # To prevent, e.g., a # like -.4 to be printed as -.39999999999
	colnames(probs) <- paste0("assertion_", 1:length(params$items))
   class(probs) <- 'margprob'
	ES <- rowSums(probs)

   return( list('probs' = probs, 'ES' = ES) )
}

## Returns the derivative value for a 3pl item, requires the probability of correct response from pr.3pl()
deriv.3pl <- function(a_uni, c_uni, ps_3pl, D) {
	return ( D*a_uni*(1 - ps_3pl)*(ps_3pl - c_uni) / (1 - c_uni) )
}

## Returns the derivative value for a gpcm item, requires the probabilities for each scorepoint from poly.pr.gpcm()
deriv.gpcm <- function(a_gpcm, ps_gpcm, ES_gpcm, D) {
	return( sapply(1:length(ps_gpcm), function(x) (D*a_gpcm[x]) * ( sum(ps_gpcm[[x]] * (1:length(ps_gpcm[[x]]))^2) - ES_gpcm[x]^2)) )
}

# information for single 3pl item.
info.3pl <- function(thetas, a_uni, b_uni, c_uni, D) {
   pr <- pr.3pl(theta = thetas, D = D, a_uni = a_uni, b_uni = b_uni, c_uni = c_uni)
   qr <- 1 - pr
   info <- (D^2) * (a_uni^2) * ( (qr / pr) * ((pr - c_uni) / (1 - c_uni))^2)
   return(info)
}

# information for single gpcm item at a single theta.
info.gpcm <- function(theta, a_gpcm, d_gpcm, D) {
   ### Important Note
   ### Under the parameterization of the GPCM in this code
   ### I need to get rid of the base category which has a
   ### pseudo step of 0

   j <- 1:(length(d_gpcm)-1)
   m <- length(d_gpcm)
   d_gpcm <- d_gpcm[2:m]
   Da <- D * a_gpcm
   kernal <- exp(cumsum(Da * (theta - d_gpcm)))
   return( Da^2 * (sum(j^2 * kernal) / (1 + sum(kernal)) - (sum(j * kernal) / (1 + sum(kernal)))^2) )
}

# information for single cluster item. params is a cluster element in a paramList.
info.cluster <- function(thetas, params, control = list()) {

   con = list(D = 1.0, Q = 50, sigma = 1)
   con[names(control)] <- control

   gq <- gauss.quad.prob(con$Q, dist = 'normal', sigma = 1)
   nodes <- gq$nodes * sqrt(params$variance)
   whts <- gq$weights

   a_par <- params$`3pl`$a
   b_par <- params$`3pl`$b
   c_par <- params$`3pl`$c

   iif <- numeric(length(thetas))
   # Consider one theta at a time.
   for(i in 1:length(iif)) {

      ES_cluster <- margprob4cluster(params = params, thetas = thetas[i], control = con)$ES

      # num_vec has one element for each node
      num_vec <- sapply(1:length(nodes), function(x) {
         p <- sapply(1:length(params$items), function(k)
            pr.latentu(a_assert = a_par[k], b_assert = b_par[k], c_assert = c_par[k], D = con$D, u_nuisance = nodes[x], theta = thetas[i]))
         return( (sum(p * (1-p)) * whts[x]) )
      })

      num <- sum(num_vec)^2

      # denom_vec also has one element for each node.
      denom_vec <- sapply(1:length(nodes), function(x) {
         p <- sapply(1:length(params$items), function(k)
            pr.latentu(a_assert = a_par[k], b_assert = b_par[k], c_assert = c_par[k], D = con$D, u_nuisance = nodes[x], theta = thetas[i]))
         return( (sum(p)^2 + sum(p * (1-p))) * whts[x])
      })

      denom <- sum(denom_vec) - (ES_cluster)^2

      iif[i] <- (num / denom)
   }
   return(iif)
}

### A slightly modified version of base R function qqline
myQQline <- function (y, datax = FALSE, distribution = qnorm, probs = c(0.25, 
    0.75), qtype = 7, ...) 
{
    stopifnot(length(probs) == 2, is.function(distribution))
    y <- as.vector(quantile(y, probs, names = FALSE, type = qtype, 
        na.rm = TRUE))
    x <- distribution(probs)
    if (datax) {
        slope <- diff(x)/diff(y)
        int <- x[[1L]] - slope * y[[1L]]
    }
    else {
        slope <- diff(y)/diff(x)
        int <- y[[1L]] - slope * x[[1L]]
    }
    list(int=int, slope=slope, ...)
}