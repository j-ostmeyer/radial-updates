update_radial = function(x, prob, d, sigma, mode){
## Use one of the implemented substitutions to perform the radial update

	gamma = rnorm(1, sd = sigma)

	if(mode == "subst_exp"){
		# x = exp(sinh(z))
		z = asinh(log(x))
		zz = z + gamma
		y = exp(sinh(zz))
		p_acc = exp(prob(y) - prob(x) + (d-1)*log(y/x) + sinh(zz)-sinh(z) + log(cosh(zz)/cosh(z)))
	}else if(mode == "subst_sinh"){
		# x = sinh(sinh(z))
		z = asinh(asinh(x))
		zz = z + gamma
		y = sinh(sinh(zz))
		p_acc = exp(prob(y) - prob(x) + (d-1)*log(abs(y/x)) + log(cosh(sinh(zz)))-log(cosh(sinh(z))) + log(cosh(zz)/cosh(z)))
	}else if(mode == "power"){
		gamma = exp(gamma)
		y = x^gamma
		p_acc = exp(prob(y) - prob(x) + (d-1)*log(y/x) + log(gamma) + (gamma-1)*log(x))
	}else if(mode == "exp"){
		y = x * exp(gamma)
		p_acc = exp(prob(y) - prob(x) + d*gamma)
	}else{
		y = x + gamma
		p_acc = exp(prob(y) - prob(x) + (d-1)*log(abs(y/x)))
	}

	acc = isTRUE(runif(1) < p_acc) & !is.na(y)
	#cat(acc+0, "\n")

	return(ifelse(acc, y, x))
}

sample_radial = function(prob, n, x0, d = 1, sigma = 1/sqrt(2), mode = "exp"){
## Produce a Markov chain using the radial update

	x = c(x0)
	for(i in 2:n) x[[i]] = update_radial(x[[i-1]], prob, d, sigma, mode)
	return(x)
}

plot_radial = function(x, prob, d, bins = 10000, trafo=function(z){z}, d.trafo=function(z){1}, i.trafo=function(z){z}, scale.range=1, ...){
## Plot the time series in a histogram

	y = trafo(x)
	data = hist(y, freq=FALSE, ...)
	range = seq(min(y), max(y)*scale.range, length.out = bins*scale.range)
	target = exp(prob(i.trafo(range)) + (d-1)*log(abs(i.trafo(range))) - log(d.trafo(i.trafo(range))))
	target = ifelse(is.na(target), 1, target)
	#length(target)
	norm = sum(target)*(max(y)-min(y))/bins
	#print(norm)
	lines(range, target/norm, lw=3, col="red")
	legend("topright", legend=c("target", "samples"), col=c("red", "black"), lw=c(3, 1), bty="n")
	invisible(data)
}

## Different negative potentials

prob_poly = function(alpha){
	#return(function(x) return(-log(1 + abs(x)^alpha)))
	return(function(x) return(-log(abs(x)^alpha)))
	#return(function(x) return(ifelse(x>=0, -log(1 + x^alpha), -Inf)))
}

prob_exp = function(alpha, z){
	return(function(x) return(-alpha*abs(x)^z))
}

prob_exp_exp = function(alpha, z){
	return(function(x) return(-alpha*cosh(z*x)))
}

prob_cos = function(alpha, z){
	return(function(x) return(-alpha*x^z + log(cos(x)^2)))
}

prob_log = function(alpha){
	return(function(x) return(-alpha*ifelse(x > 1, log(x)*log(log(x)), 0)))
}

prob_sinc = function(alpha){
	return(function(x) return(alpha*log(abs(sin(x)/x))))
}
