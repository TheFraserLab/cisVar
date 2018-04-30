#!/usr/bin/R

#        AUTHOR: Ashley Tehranchi, tehranchi576@gmail.com
#  ORGANIZATION: Stanford University

ldCalculation <- function(sig, outfile) {

  # Max distance between two SNPs to use for FDR calculation
  fdr.window = 200000

  r2.sig = vector("numeric", length=nrow(sig))
  holdnum = nrow(sig)-1
  maxcol = ncol(sig)

  # If two SNPs are within FDR window (default 200 kb) do correlation between
  # their genotypes Otherwise just set the r-squared to 0
  for (i in 1:holdnum) {
      if (sig[i,2] - sig[i+1,2] < fdr.window) {
            r2.sig[i] <- cor(t(sig[i,20:maxcol]), t(sig[i+1,20:maxcol]), method="pearson")
    }
    else{ r2.sig[i] <- 0}
  }

  # Add r-squared
  cat("Adding r-squared\n"); flush.console()
  sig$r2 <- (r2.sig)^2
  ldfileOUTsig = sig[, c("Chr", "position", "r2")]
  ldfileOUTsigFinal = ldfileOUTsig[complete.cases(ldfileOUTsig),]

  cat("Writing out LD data\n"); flush.console()
  write.table(ldfileOUTsigFinal, file=outfile, sep="\t", quote=F, row.names=F, col.names=T)
}

linearRegression <- function(pd, gt, tmppfx, indiv) {

  ###############
  #  Constants  #
  ###############

  # MAF
  MAF.low = 0.02
  MAF.high = 1-MAF.low

  # P-Value, used to split sig/nonsig for LD calculation
  P = 0.05

  # Force postfreq to be between eps and 1-eps (e.g 0.0001 and 0.9999)
  eps = 0.0001

  # Figure sizing/format defautls
  fig.width  = 1200
  fig.height = 1200
  fig.font   = 40

  ##################
  #  Read in Data  #
  ##################

  cat("..cisVAR Regression..\n"); flush.console()
  options("scipen"=100)
  n.indi = as.numeric(indiv)
  print(pd)
  cat("Loading POST..\n"); flush.console()
  postTotal = read.table(pd, sep="\t", header=T)
  cat("Loading Genotypes\n"); flush.console()
  genofile = matrix(scan(file = gt, sep = "\t", what = double(), nlines = n.indi), byrow = TRUE, nrow = n.indi)
  cat("Genotype Dimensions: "); flush.console()
  print(dim(genofile))

  genoTotal = 1 - (0.5 * genofile)

  ############################
  #  Restrict to MAF Cutoff  #
  ############################


  # We bind the genotype and transposed post matrices together and
  # then split them again so that MAF filtering applies equally to
  # both.

  # Boundaries of the two matrices
  post.start = 1
  post.end=15
  geno.start=16
  geno.end=15 + n.indi

  # Bind the matrices
  cat("Combining\n"); flush.console()
  a = t(genoTotal)
  b = cbind(postTotal,a)

  # MAF Filter
  cat(sprintf("Filtering by MAF %f\n", MAF.low)); flush.console()
  mafsub = b[apply(b[,geno.start:geno.end], MARGIN = 1, function(x) mean(as.numeric(x)) >= MAF.low & mean(as.numeric(x)) <= MAF.high), ]

  # Split the matrices again
  cat("Splitting\n"); flush.console()
  post = mafsub[,post.start:post.end]
  genotypes = as.matrix(mafsub[,geno.start:geno.end])
  # postTemp = mafsub

  ################
  #  Regression  #
  ################

  n.snps          = nrow(post)
  # genotypes       = as.matrix(t(genos))
  weights         = rep(1/n.indi,n.indi)
  real.props      = genotypes %*% as.matrix(weights)
  depths          = post$Depth
  postProp        = post$POSTfreq
  estimated.props = postProp

  G = genotypes[,-n.indi] - genotypes[,n.indi]
  Y = postProp - genotypes[,n.indi]

  # Force postfreq to be between eps and 1-eps (e.g 0.0001 and 0.9999)
  cat(sprintf("Filtering POST frequencies by eps of %f\n", eps)); flush.console()
  props.for.weights = pmin(1-eps,pmax(estimated.props,eps))

  # weight = depth / (adjusted.post * (1-adjusted.post))
  cat("Calculating Weights\n"); flush.console()
  regression.weights = depths / (props.for.weights * (1-props.for.weights) )
  good = which( (postProp>0.1) & (postProp < 0.9))

  # Regression is here:
  m = lm(Y[good] ~ G[good,]-1,weights=regression.weights[good])  ## run without intercept
  coefs = m$coef
  s = summary(m)
  cov.mat = s$cov.unscaled * s$sigma^2
  big.cov.mat = matrix(NA,n.indi,n.indi)
  big.cov.mat[-n.indi,-n.indi] = cov.mat
  big.cov.mat[n.indi,n.indi] = sum(cov.mat)
  big.cov.mat[n.indi,-n.indi] = big.cov.mat[-n.indi,n.indi] = -rowSums(cov.mat)
  cat("Regression...\n"); flush.console()
  vars = sapply(1:n.snps, function(i) genotypes[i,] %*% big.cov.mat %*% genotypes[i,])
  vars[vars < 0] <- 0
  all.coeffs = c(coefs, 1-sum(coefs))
  preProps = genotypes %*% as.matrix(all.coeffs)
  preProps[preProps > 1] <- 1
  preVars = vars
  postProps = estimated.props
  postProps[postProps > 1] <- 1
  postVars = 1/regression.weights
  postVars[postVars > 1] <- 1
  cat("Calculating Z values\n"); flush.console()
  Zs = (postProps - preProps)/sqrt(preVars + postVars)
  cat("Calculating P values\n"); flush.console()
  p.values = 2*(1-pnorm(abs(Zs)))

  analytic.pv = function(preProp,preVar,postProp,depth){
    possible.props = 0:depth / depth
    extreme.props = (0:depth)[abs(possible.props - preProp) >= abs(postProp-preProp)]
    if(preVar < 0.0001){
      return(sum(dbinom(extreme.props,depth,preProp)))
    }
    f = function(x){
      denom = pnorm(1,mean=preProp,sd=sqrt(preVar)) - pnorm(0,mean=preProp,sd=sqrt(preVar))
      dnorm(x,mean=preProp,sd=sqrt(preVar)) * sum(dbinom(extreme.props,depth,x)) / denom
    }
    integrate(f,0,1)$value
  }

  an.pvs<-sapply(1:n.snps, function(i) analytic.pv(preProps[i],preVars[i],postProps[i],depths[i]))

  cat("Adding Statistics to Matrix\n"); flush.console()
  post$prechipfreq <- preProps
  post$pvalue <- an.pvs
  post$zvalue <- Zs
  post$prevar <- preVars
  post$postvar <- postVars
  post$SNPpostfreq <- 1 - post$POSTfreq
  post$SNPprefreq <- 1 - post$prechipfreq

  ##############################################
  #  LD Calculation for Significance Estimate  #
  ##############################################

  cat("LD Calculation\n"); flush.console()
  combined = cbind(post, genotypes)

  sig = subset(combined, pvalue <= P)
  totaloutsig = paste(tmppfx, ".sigLD", sep="")
  ldCalculation(sig, totaloutsig)
  rm(sig)

  ########################################
  #  LD Calculation for Non-Significant  #
  ########################################

  nonsig = subset(combined, pvalue > P)
  totaloutnonsig = paste(tmppfx, ".nonsigLD", sep="")
  ldCalculation(nonsig, totaloutnonsig)
  rm(nonsig)
  rm(combined)

  ###################
  #  Write Outputs  #
  ###################

  cat("Filtering pre-freq between 0 and 1 only\n"); flush.console()
  ## remove SNPs with pre = 0 or 1
  # print(dim(post))
  temp = subset(post, prechipfreq>0 & prechipfreq<1)
  ## must be sorted by pvalue
  x = temp[with(temp, order(pvalue)), ]
  # print(dim(x))
  xs = x
  options("scipen"=100)
  xs$start <- xs$position - 1

  cat("Writing outputs\n"); flush.console()
  bed = xs[,c("Chr", "start", "position")]  ## reorder by column numbers
  sortedBed = bed[with(bed,order(Chr, start, position)),]
  finaloutFile = paste(tmppfx, ".total.txt", sep="")
  finalcoeffs = paste(tmppfx, ".coefficients.txt", sep="")
  finalbed = paste(tmppfx, ".final.bed", sep="")
  write.table(x, file=finaloutFile, sep="\t", quote=F, row.names=F, col.names=T)
  write.table(sortedBed, file=finalbed, sep="\t", quote=F, row.names=F, col.names=F)
  write.table(all.coeffs, file=finalcoeffs, sep="\t", quote=F, row.names=F, col.names=T)

  ###################
  #  Summary Plots  #
  ###################

  cat("Making summary plots\n"); flush.console()
  png(paste(tmppfx, ".summaryPlots.png", sep=""), width=fig.width, height=fig.height, pointsize=fig.font)
  par(mfrow=c(2,2))
  plot(table(round(p.values,d=2)))
  qqplot(p.values,runif(n.snps), )#,xlim=c(0,0.1),ylim=c(0,0.1))
  abline(0,1)
  plot(table(round(an.pvs,d=2)))
  qqplot(an.pvs,runif(n.snps))#,xlim=c(0,0.1),ylim=c(0,0.1))
  abline(0,1)
  dev.off()

  png(paste(tmppfx, ".qqplot.png", sep=""), width=fig.width, height=fig.height, pointsize=fig.font)
  qqplot(runif(n.snps), an.pvs, xlab="Expected", ylab="Observed")#,xlim=c(0,0.1),ylim=c(0,0.1))
  abline(0,1)
  dev.off()
  png(paste(tmppfx, ".qqplot.log.png", sep=""), width=fig.width, height=fig.height, pointsize=fig.font)
  qqplot(runif(n.snps), an.pvs, xlab="Expected", ylab="Observed")#,xlim=c(0,0.1),ylim=c(0,0.1))
  abline(0,1)
  dev.off()

  cat("Done!\n"); flush.console()


}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Rscript regression_qtls.R <post> <genos> <prefix> <indiv>.", call.=FALSE)
}
post  = args[1]
genos = args[2]
prefix = args[3]
indiv = as.numeric(args[4])

linearRegression(post, genos, prefix, indiv)
