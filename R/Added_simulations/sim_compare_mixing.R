args = commandArgs(trailingOnly=TRUE)

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


y.fix.inter = read.csv(sprintf("data/rep%s/covariates.csv", args[1]), header = F)
data = as.matrix(read.csv(sprintf("data/rep%s/%s.csv", args[1], args[2]), header = F))
y.fix.inter = t(bs( y.fix.inter[,1] ))

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

nod.free.mcmc.new( data, start.sub, hyper, sigma.value, sigma.prior, sprintf("mixing/rep%s/%s/chain%s/sim",args[1],args[2],args[3]), burnin = 0.2, step = 2e5, thin = 100 )
