args = commandArgs(trailingOnly=TRUE)
library( pspline )
library( splines )

#args[1] replicates
#args[2] sparsity
#args[3] pure-error
#args[4] read-depth overdispersion

all.spar = as.numeric(strsplit(args[2], split = ",", fixed = T)[[1]])
all.er = as.numeric(strsplit(args[3], split = ",", fixed = T)[[1]])
all.nover = as.numeric(strsplit(args[4], split = ",", fixed = T)[[1]])

gen.data = function(rep, spar, er, n.over){
  n = 300
  n.sub = 50
  p = 100
  if(n.over==0){
    tot.count = rpois(n, 1e5)
  }else{
    tot.count = rnbinom(n, mu = 1e5, size = 1/n.over)
  }
  
  #subject identifier
  subject.id = data.frame( id=as.factor(rep(1:n.sub, each = n/n.sub)) )
  sub.design = t(model.matrix( ~.-1, subject.id ))
  sub.Z = data.frame( rep(1:n.sub, each = n/n.sub) )
  
  sim.res = readRDS( sprintf("data/rep%s/sim_data", rep) )
  
  sigma = sim.res$sigma
  y.raw = sim.res$y.raw
  v.fix.inter = sim.res$v.fix.inter
  X = sim.res$X
  Y.sub = sim.res$Y.sub
  y.fix.inter = t(model.matrix( ~ysmall-1, data = y.raw))
  y.spline = bs( y.fix.inter[1,] )
  
  Q = t(X)%*%Y.sub%*%sub.design + t(v.fix.inter)%*%y.fix.inter + matrix( rnorm( n*p, sd = sqrt(er) ), nrow = p )
  
  # weights
  w.DF = Q^2*(Q>0)*sigma
  w.MI = exp(Q)
  # prob
  w.DF.n = t(t(w.DF)/colSums(w.DF))
  w.MI.n = t(t(w.MI)/colSums(w.MI))
  # truncation
  w.DF.n.t = w.DF.n*(w.DF.n>=spar)
  w.MI.n.t = w.MI.n*(w.MI.n>=spar)
  # reads
  data.DF = sapply(1:ncol(w.DF.n.t), function(i) rmultinom(1, tot.count[i], w.DF.n.t[,i]))
  data.MI = sapply(1:ncol(w.MI.n.t), function(i) rmultinom(1, tot.count[i], w.MI.n.t[,i]))
  
  # save the data
  # first for DirFactor
  write.table(data.DF, sprintf("data/rep%s/count_DF_spar_%s_er_%s_ns_%s.csv", rep, spar, er, n.over), sep=",", row.names = F, col.names = F)
  write.table(data.MI, sprintf("data/rep%s/count_MI_spar_%s_er_%s_ns_%s.csv", rep, spar, er, n.over), sep=",", row.names = F, col.names = F)
  write.table(y.raw, sprintf("data/rep%s/covariates.csv", rep), sep = ",", row.names = F, col.names = F)
  
  # then for MIMIX
  dir1 = sprintf("count_DF_spar_%s_er_%s_ns_%s", spar, er, n.over)
  dir2 = sprintf("count_MI_spar_%s_er_%s_ns_%s", spar, er, n.over)
  if(!dir.exists(file.path(sprintf("data/rep%s", rep), dir1))){
    dir.create(file.path(sprintf("data/rep%s", rep), dir1))
  }
  if(!dir.exists(file.path(sprintf("data/rep%s", rep), dir2))){
    dir.create(file.path(sprintf("data/rep%s", rep), dir2))
  }
  write.table(t(data.DF), paste(sprintf("data/rep%s/", rep), dir1, "/Y.csv", sep=""), sep=",", row.names = F, col.names = F)
  write.table(t(data.MI), paste(sprintf("data/rep%s/", rep), dir2, "/Y.csv", sep=""), sep=",", row.names = F, col.names = F)
  write.table(y.spline, paste(sprintf("data/rep%s/", rep), dir1, "/X.csv", sep=""), sep=",", row.names = F, col.names = F)
  write.table(y.spline, paste(sprintf("data/rep%s/", rep), dir2, "/X.csv", sep=""), sep=",", row.names = F, col.names = F)
  write.table(sub.Z, paste(sprintf("data/rep%s/", rep), dir1, "/Z.csv", sep=""), sep=",", row.names = F, col.names = F)
  write.table(sub.Z, paste(sprintf("data/rep%s/", rep), dir2, "/Z.csv", sep=""), sep=",", row.names = F, col.names = F)
}

for( spar in all.spar ){
  for( er in all.er ){
    for( nover in all.nover ){
      gen.data(args[1], spar, er, nover)
    }
  }
}