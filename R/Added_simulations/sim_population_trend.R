library(pspline)
library(splines)

args = commandArgs(trailingOnly=TRUE)

ListtoArray = function( ls.input ){
  # convert list of matrices to an array
  array( unlist( ls.input ), dim = c(dim(ls.input[[1]]),length(ls.input)) )
}

res.files = sprintf("res/rep%s/%s/res_%s", args[1], args[2], seq(20099,99999,100))
out.path = sprintf("trend/rep%s/%s_trend.rds", args[1], args[2])
y.raw = read.csv(sprintf("data/rep%s/covariates.csv", args[1]), header = F)
y.fix.inter = bs(y.raw$V1)
# population trends for DirFactor restults
y.seq = seq(-2,2,length.out = 20)
n.sub = 50

get.pop.mean = function( y.fix, sim.ind ){
  Q = t(sim.ind$X)%*%sim.ind$Y.sub + t(sim.ind$x)%*%y.fix  + 
    matrix( rnorm( ncol(sim.ind$X)*ncol(sim.ind$Y.sub), sd = sqrt(sim.ind$er) ), nrow = ncol(sim.ind$X) )
  weight = sim.ind$sigma*Q^2*(Q>0)
  apply( t(weight)/colSums(weight), 2, mean )
}

pop.mean.trend = ListtoArray( lapply( res.files, function(res.file){
  x = readRDS(res.file)
  sapply( y.seq, function(xx){
    y.tmp.1 = predict(y.fix.inter, xx)
    y.tmp = matrix( rep(y.tmp.1,n.sub), nrow=3 )
    get.pop.mean(y.tmp, x)} )} ) )

saveRDS(pop.mean.trend, out.path)