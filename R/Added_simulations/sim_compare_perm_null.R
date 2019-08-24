args = commandArgs(trailingOnly=TRUE)

#args[1] input path
#args[2] repeat index: rep*permutation times

idx = as.numeric(args[2])
rep = floor((idx-1)/100)+1
perm = (idx-1)%%100 + 1

source("/n/hutlab12_nobackup/users/bren/DirFactor_fix_compare/src/mcmc_sampler.R")
n = 300
n.sub = 50
p = 100

subject.id = data.frame( id=as.factor(rep(1:n.sub, each = n/n.sub)) )
sub.design = t(model.matrix( ~.-1, subject.id ))
hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, 
              sub.design = sub.design, alpha = 10, beta = 0, a.x.sigma=5, b.x.sigma=5 )
sigma.value = seq(0.001,0.999,0.001)
tmp = c(0,pbeta( sigma.value, 10/68, 1/2-10/68 ))
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)

y.fix.inter = read.csv(sprintf("%s/rep%s/covariates.csv", args[1], rep), header = F)
if(perm<100){
  y.fix.inter = matrix(sample(y.fix.inter$V1),nrow = 1)
}else{
  y.fix.inter = matrix(y.fix.inter$V1,nrow = 1)
}

data = as.matrix(read.csv(sprintf("%s/rep%s/count.csv", args[1], rep), header = F))

start.sub = list( er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ), 
                  sigma = sample( sigma.value, size = p, replace = T, prob = sigma.prior ), 
                  T.aug = rgamma( n, 10, 1 ),
                  Q = matrix( 0.5, nrow = p, ncol = n ),
                  X = matrix( rnorm( hyper$m*p ), nrow = hyper$m ),
                  Y.sub = matrix( 0, nrow = hyper$m, ncol = n.sub ),
                  delta = c( rgamma( 1, shape = hyper$a1, rate = 1 ), rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) ),
                  phi = matrix( rgamma( hyper$m*n.sub, shape = 3/2, rate = 3/2 ), nrow = n.sub ),
                  y = y.fix.inter, 
                  x = matrix( rnorm( p*nrow(y.fix.inter) ), nrow=nrow(y.fix.inter) ),
                  x.sigma = rgamma( nrow(y.fix.inter), 1, 1 ) )

if(!dir.exists(sprintf("%s/rep%s/perm%s", args[1], rep, perm))){
  dir.create(sprintf("%s/rep%s/perm%s", args[1], rep, perm), recursive = T)
}

nod.free.mcmc.new( data, start.sub, hyper, sigma.value, sigma.prior, sprintf("%s/rep%s/perm%s/res", args[1], rep, perm), burnin = 0.2, step = 1e5, thin = 500 )
