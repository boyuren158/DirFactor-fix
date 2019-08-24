library(pspline)
library(splines)

args = commandArgs(trailingOnly=TRUE)

ListtoArray = function( ls.input ){
  # convert list of matrices to an array
  array( unlist( ls.input ), dim = c(dim(ls.input[[1]]),length(ls.input)) )
}

res.file = sprintf("res_mimix/rep%s/%s", args[1], args[2])
out.path = sprintf("trend_mimix/rep%s/%s_trend.rds", args[1], args[2])
y.raw = read.csv(sprintf("data/rep%s/covariates.csv", args[1]), header = F)
y.fix.inter = bs(y.raw$V1)
# population trends for DirFactor restults
y.seq = seq(-2,2,length.out = 20)
n.sub = 50
print(res.file)

process.MIMIX = function(path, L){
  b = as.matrix(read.csv(sprintf("%s/b.tsv", path), sep="\t"))
  g = as.matrix(read.csv(sprintf("%s/g.tsv", path), sep="\t"))
  Lambda = as.matrix(read.csv(sprintf("%s/Lambda-postmean.tsv", path), sep="\t", header = F))
  mu = as.matrix(read.csv(sprintf("%s/mu.tsv", path), sep="\t"))
  theta_var = as.matrix(read.csv(sprintf("%s/theta_var.tsv", path), sep="\t"))
  
  n.iter = nrow(theta_var)
  
  ListtoArray(lapply(1:n.iter, function(i){
    b.i = matrix(c(b[i,]),nrow=L)
    g.i = matrix(g[i,],nrow=L)
    mu.i = mu[i,]
    theta_var.i = theta_var[i,]
    
    tmp = sapply(y.seq, function(y.ind){
      y.ind.bs = predict(y.fix.inter, y.ind)
      f = c(b.i%*%t(y.ind.bs)) + g.i + matrix(rnorm(L*n.sub),nrow = L)
      theta = t(Lambda)%*%f + mu.i + matrix(rnorm(ncol(Lambda)*ncol(f), sd = sqrt(c(theta_var.i))),nrow=ncol(theta_var))
      
      Q = exp(theta)
      rowMeans(t(t(Q)/colSums(Q)))
    })
  }))
}

pop.mean.trend = process.MIMIX(res.file, L=10)
saveRDS(pop.mean.trend, out.path)