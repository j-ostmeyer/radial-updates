source("radial_updates.R")

library(tikzDevice)
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ), "\\usepackage{amsmath}", "\\DeclareMathOperator{\\asinh}{asinh}"))

# V(r) = r, 100D, get histogram

res = sample_radial(prob_exp(1, 1), 100000, 100, d=100, sigma=sqrt(2/100), mode = "exp")

tikz('hist_V-r_100D.tex', width=2.7, height=2.7)
par(mar = c(4,4,0,0))
plot_radial(res, prob_exp(1, 1), 100, breaks=100, xlab="$r$", ylab="$\\rho(r)$", main=NULL)
dev.off()

# V(r) = ln(1+|r|^1.01), 1D, get histogram

res = sample_radial(prob_poly(1.01), 100000, 1, sigma=sqrt(2), mode = "subst_exp")

tikz('hist_V-ln1001_1D.tex', width=2.7, height=2.7)
par(mar = c(4,4,0,0))
plot_radial(res, prob_poly(1.01), 1,
			#trafo=asinh, d.trafo=function(z) 1/sqrt(1+z^2), i.trafo=sinh,
			trafo=log10, d.trafo=function(z) 1/z/log(10), i.trafo=function(z) 10^z,
			breaks=100, xlab="$\\log_{10}(r)$", ylab="$\\rho(\\log_{10}(r))$", main=NULL)
dev.off()

# V(r) = ln(1+|r|^1.1), 1D, get LINEAR histogram and time series

res = sapply(1:3, function(i) sample_radial(prob_poly(1.1), 100000, 1, sigma=sqrt(2), mode = "lin"))
write.table(res, "ts_V-ln11_1D-lin.txt", row.names=FALSE, col.names=FALSE)

tikz('hist_V-ln11_1D-lin.tex', width=2.7, height=2.7)
par(mar = c(4,4,0,0))
plot_radial(c(res), prob_poly(1.1), 1,
			#trafo=asinh, d.trafo=function(z) 1/sqrt(1+z^2), i.trafo=sinh,
			trafo=log10, d.trafo=function(z) 1/z/log(10), i.trafo=function(z) 10^z, bins=3000, scale.range=5,
			breaks=60, xlab="$\\log_{10}(r)$", ylab="$\\rho(\\log_{10}(r))$", main=NULL, xlim=c(-3,5))
dev.off()

# V(r) = 1/2 r^2, 100D, get time series

dim0 = 100
sigma = 10^(c(-20:20)/10)/sqrt(dim0)
dims = round(10^(c(-10:10)/10)*dim0)
#dims = c(10,100,1000)

for(dim in dims){
	res = sapply(sigma, function(s) sample_radial(prob_exp(.5, 2), 10000, sqrt(dim), d=dim, sigma=s, mode = "exp"))

	write.table(sigma, sprintf("sigma_V-r2_%dD.txt", dim), row.names=FALSE, col.names=FALSE)
	write.table(res, sprintf("ts_V-r2_%dD.txt", dim), row.names=FALSE, col.names=FALSE)
}

# V(r) = exp(r), 1D, get drift
res = sapply(1:3, function(i) sample_radial(prob_exp_exp(1, 1), 10000, 200, d=1, sigma=1, mode = "lin"))
write.table(cosh(res), "ts_V-exp_1D.txt", row.names=FALSE, col.names=FALSE)

# V(r) = r, 1D, get drift
res = sapply(1:3, function(i) sample_radial(prob_exp(1, 1), 10000, 200, d=1, sigma=1, mode = "lin"))
write.table(abs(res), "ts_V-r_1D.txt", row.names=FALSE, col.names=FALSE)

# V(r) = sqrt(r), 1D, get drift
res = sapply(1:3, function(i) sample_radial(prob_exp(1, 0.5), 100000, 200, d=1, sigma=1, mode = "lin"))
write.table(sqrt(abs(res)), "ts_V-sqrt_1D.txt", row.names=FALSE, col.names=FALSE)

# V(r) = ln(1+r^2), 1D, get drift
res = sapply(1:3, function(i) sample_radial(prob_poly(2), 100000, 200, d=1, sigma=1, mode = "lin"))
write.table(log(1+res^2), "ts_V-log_1D.txt", row.names=FALSE, col.names=FALSE)
