tot.count = 1e5

#try to fit the data to our model
source("MCMC.R")
source("utilities.R")
library( MASS )
library( truncnorm )
library( mvtnorm )
library( MCMCpack )
library( mvtnorm )
library( tmvtnorm )
library( pspline )
library( splines )

nod.free.mcmc.new = function( data, start, hyper, 
                          sigma.value, sigma.prior,
                          save_path, 
                          save.obj = c("sigma", "Q", "T.aug", "X", "Y.sub", "delta", "phi", "x", "er", "x.sigma" ), 
                          burnin = 0.2, thin = 5, step = 1000 ){
  nv = hyper$nv
  a.er = hyper$a.er
  b.er = hyper$b.er
  #er is the variance of the pure error
  a1 = hyper$a1
  a2 = hyper$a2
  m = hyper$m
  #x.sigma variance of the regression coef
  a.x.sigma = hyper$a.x.sigma
  b.x.sigma = hyper$b.x.sigma
  #design matrix for subject
  #dim = n.sub*n.sample
  sub.design = hyper$sub.design
  sub.rep = rowSums(sub.design)
  
  n = ncol( data )
  n.sub = nrow( sub.design )
  p = nrow( data )
  sum.species = rowSums( data )
  sum.sample = colSums( data )
  sv.log = log( sigma.value )
  sp.log = log( sigma.prior )
  
  #fixed effect
  x.old = start$x
  x.sigma.old = start$x.sigma
  y.fix = start$y
  
  #random effect
  sigma.old = as.vector( start$sigma )
  Q.old = start$Q
  T.aug.old = as.vector( start$T.aug )
  er.old = start$er
  X.old = start$X
  Y.sub.old = start$Y.sub
  Y.old = Y.sub.old%*%sub.design
  
  delta.old = as.vector(start$delta)
  phi.old = start$phi
  
  all.cache = c()
  t.t = proc.time()
  for( iter in 2:step ){
        # browser()
    if( iter %% 1000 == 0 )
      print( iter )
    #sample sigma
    #this is a major time consuming step
    all.cache$sigma = sigma.old
    #     set.seed(1)
    
    sigma.old = sigma.gibbs(
      sigma.old, T.aug.old,
      Q.old, sum.species,
      sigma.value, sigma.prior,
      sv.log, sp.log
    )
    
    #     Q.pos = Q.old*(Q.old>0);
    #     tmp.weights = Q.pos^2%*%T.aug.old;
    #     sigma.old.t = sigma_gibbs( length( sigma.value ), 
    #                              p, 
    #                              sum.species, 
    #                              tmp.weights, 
    #                              sigma.value, 
    #                              sigma.prior )
    
    
    #sample T
    all.cache$T.aug = T.aug.old
    T.aug.old = T.aug.gibbs( 
      T.aug.old, sigma.old, 
      Q.old, sum.sample
    )
    
    #sample Q
    #matrix of all parameters
    #notice the mean of the prior conditional on fixed effect and random effect is changed
    Q.para.all = cbind( c(data), c(t(X.old)%*%Y.old+t(x.old)%*%y.fix), 
                        rep( sigma.old, n ), c(Q.old),
                        rep( T.aug.old, each = p ), rep( sqrt(er.old), n*p ) )
    all.cache$Q = Q.old
    #get samples
    
    labels = (Q.para.all[,1]==0)
    Q.0 = Q.gibbs.vec.0( Q.para.all[labels, c(2,3,5,6)] )
    Q.n0 = Q.gibbs.vec.n0( Q.para.all[!labels, c(1,2,3,5,6)], 
                           Q.para.all[!labels, 4] )
    #save it into proper locations
    Q.tmp = rep( 0, nrow( Q.para.all ) )
    Q.tmp[labels] = Q.0
    Q.tmp[!labels] = Q.n0
    #convert it back to matrix
    Q.old = matrix( Q.tmp, nrow = p )
    # browser()
    
    #sample Y.sub
    #this is a major time consuming step
    #get tau
    tau.tmp = cumprod( delta.old )
    all.cache$Y.sub = Y.sub.old
    Y.mean.pre = X.old%*%(Q.old-t(x.old)%*%y.fix)%*%t(sub.design)/er.old
    
    Y.sub.old = vapply( 1:n.sub, function(x){
      Sigma.Y = solve( chol( sub.rep[x]*X.old%*%t(X.old)/er.old + diag( phi.old[x,]*tau.tmp ) ) )
      mu.Y = t(Sigma.Y)%*%Y.mean.pre[,x,drop=F]
      Sigma.Y%*%( rnorm( length(mu.Y) ) + mu.Y )
    }, rep( 0, nrow( Y.sub.old ) ) )
    Y.old = Y.sub.old%*%sub.design
    #     Y.old = matrix( 0, nrow = nrow(Y.old), ncol = ncol( Y.old ) )
    
    #sample er
    #prior is ga(1,10)
    all.cache$er = er.old
    er.old = 1/rgamma( 1, shape = n*p/2+a.er, 
                       rate = sum((Q.old-t(X.old)%*%Y.old-t(x.old)%*%y.fix)^2)/2+b.er )
    
    #sample X
    Sigma.X = solve( Y.old%*%t(Y.old)/er.old + diag( m ) )
    all.cache$X = X.old
    #     X.mean = vapply( 1:p, function(x){
    #       Sigma.X%*%Y.old%*%Q.old[x,]/er.old
    #     }, rep( 0, nrow( X.old ) ) )
    X.mean = t(Sigma.X)%*%Y.old%*%t(Q.old-t(x.old)%*%y.fix)/er.old
    X.old = X.mean + t( rmvnorm( p, sigma = Sigma.X ) )
    
    #sample phi
    #prior is ga(nv/2,nv/2)
    #nv = 3
    all.cache$phi = phi.old
    phi.old = matrix( rgamma( n.sub*m, shape = (nv+1)/2, rate = c( Y.sub.old^2*tau.tmp + nv )/2 ), 
                      nrow = n.sub, byrow = T )
    
    #sample delta
    #we choose a1 = 3, a2 = 4
    delta.tmp = delta.old
    all.cache$delta = delta.old
    for( k in 1:m ){
      if( k == 1 ){
        delta.prime = cumprod( delta.tmp )/delta.tmp[1]
        delta.tmp[k] = rgamma( 1, shape = n.sub*m/2 + a1, 
                               rate = 1 + sum(colSums( t(Y.sub.old^2)*phi.old )*delta.prime)/2
        )
      }
      else{
        delta.prime = cumprod( delta.tmp )[-(1:(k-1))]/delta.tmp[k]
        delta.tmp[k] = rgamma( 1, shape = n.sub*(m-k+1)/2 + a2, 
                               rate = 1 + sum( colSums( t(Y.sub.old^2)*phi.old )[-(1:(k-1))]*delta.prime)/2
        )
      }
    }
    delta.old = delta.tmp
    
    #sample the fixed effect reg coef variances
    all.cache$x.sigma = x.sigma.old
    a.x.sigma.use = a.x.sigma + p/2
    b.x.sigma.use = b.x.sigma + rowSums(x.old^2)/2
    x.sigma.old = 1/rgamma(length(b.x.sigma.use), shape = a.x.sigma.use, rate = b.x.sigma.use)
    
    #sample the species specific effect x
    #x is a matrix, dimension = K*I, where K is the number of covariates
    all.cache$x = x.old
    x.cov = solve( (y.fix%*%t(y.fix))/er.old + diag(1/x.sigma.old, nrow=length(x.sigma.old), ncol=length(x.sigma.old)) )
    x.mean.all = x.cov%*%y.fix%*%t(Q.old-t(X.old)%*%Y.old)/er.old
    x.old = x.mean.all + t( rmvnorm( ncol(x.mean.all), sigma = x.cov ) )
    
    if( iter > burnin*step & iter%%thin == 0 )
      saveRDS( all.cache[save.obj], file = paste( save_path, iter-1, sep = "_") )
  }
  
  #save last run
  all.cache = list( sigma = sigma.old, T.aug = T.aug.old, Q = Q.old, Y.sub = Y.sub.old,
                    X = X.old, phi = phi.old, delta = delta.old, x = x.old, 
                    x.sigma = x.sigma.old, er = er.old )
  saveRDS( all.cache[save.obj], file = paste( save_path, iter, sep = "_") )
  
  print( "Total Time:")
  print( proc.time()-t.t )
  #   list( Q = Q, T.aug = T.aug, sigma = sigma, Y = Y, er = er, X = X, phi = phi, delta = delta )
}

n = 300
n.sub = 50
p = 100

Y.sub = matrix( rnorm( n.sub*4 ), nrow = 4 )
Y.sub[1:2, (n.sub/2+1):n.sub] = 0
Y.sub[3:4, 1:(n.sub/2)] = 0
X = matrix( rnorm( p*4 ), nrow = 4 )

subject.id = data.frame( id=as.factor(rep(1:n.sub, each = n/n.sub)) )#data.frame( id = as.factor(sample( rep(1:50, each = 6), size = 300, replace = F )) )
sub.design = t(model.matrix( ~.-1, subject.id ))

y.raw = data.frame( ysmall = rnorm(n) )
y.fix.inter = t(model.matrix( ~ysmall-1, data = y.raw))

v.fix.inter = matrix(c( rep(5,1), rep(-5,1), rep(0, p-2) ), nrow = 1)

hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, 
              sub.design = sub.design, alpha = 10, beta = 0, a.x.sigma=5, b.x.sigma=5 )
sigma.value = seq(0.001,0.999,0.001)
tmp = c(0,pbeta( sigma.value, 10/68, 1/2-10/68 ))
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)

sigma = sample( sigma.value, p, replace = T, prob = sigma.prior )
sigma.normal = rbeta( p, 1/2, 1 )

Q = t(X)%*%Y.sub%*%sub.design + t(v.fix.inter)%*%y.fix.inter + matrix( rnorm( n*p ), nrow = p )
weight = exp(Q)
weight.norm = t(t(weight)/colSums(weight))
weight.norm.cut = weight.norm*(weight.norm>=0)
data = apply( weight.norm.cut, 2, function(x) rmultinom( 1, 1e5, x ) )

y.fix.inter = t(bs( y.fix.inter[1,] ))

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

nod.free.mcmc.new( data, start.sub, hyper, sigma.value, sigma.prior, "sim_exp_0/sim", burnin = 0.2, step = 1e5, thin = 50 )

mcmc.res.simple = lapply( paste("sim_exp_0/sim", seq(20049,99999,50), sep = "_" ), readRDS )

get.pop.mean = function( y.fix, sim.ind ){
  Q = t(sim.ind$X)%*%sim.ind$Y.sub + t(sim.ind$x)%*%y.fix  + 
    matrix( rnorm( ncol(sim.ind$X)*ncol(sim.ind$Y.sub), sd = sqrt(sim.ind$er) ), nrow = ncol(sim.ind$X) )
  weight = sim.ind$sigma*Q^2*(Q>0)
  apply( t(weight)/colSums(weight), 2, mean )
}

get.pop.mean.tru = function( y.fix, sim.ind ){
  Q = t(sim.ind$X)%*%sim.ind$Y.sub + t(sim.ind$x)%*%diag(sim.ind$d,nrow=length(sim.ind$d),ncol=length(sim.ind$d))%*%y.fix  + 
    matrix( rnorm( ncol(sim.ind$X)*ncol(sim.ind$Y.sub), sd = sqrt(sim.ind$er) ), nrow = ncol(sim.ind$X) )
  weight = exp(Q)
  apply( t(weight)/colSums(weight), 2, mean )
}

truth.curve = apply( ListtoArray( lapply( 1:1000, function(iii) {sapply( seq(-2,2,length.out = 20), function(xx){
  y.tmp = matrix( rep(xx,n.sub), nrow = 1 )
  get.pop.mean.tru(y.tmp, list(X=X, Y.sub=Y.sub, sigma = sigma.normal, d = 1, x = v.fix.inter, er = 1))} )} ) ), 1:2, mean )

tmp.1 = ListtoArray( lapply( mcmc.res.simple, function(x){
  sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp.1 = predict(y.fix.inter, xx)
    y.tmp = matrix( rep(y.tmp.1,n.sub), nrow=3 )
    get.pop.mean(y.tmp, x)} )} ) )
tmp.2 = ListtoArray( lapply( mcmc.res.simple.1, function(x){
  sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp.1 = predict(y.fix.inter, xx)
    y.tmp = matrix( rep(y.tmp.1,n.sub), nrow=3 )
    get.pop.mean(y.tmp, x)} )} ) )
tmp.3 = ListtoArray( lapply( mcmc.res.simple.2, function(x){
  sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp.1 = predict(y.fix.inter, xx)
    y.tmp = matrix( rep(y.tmp.1,n.sub), nrow=3 )
    get.pop.mean(y.tmp, x)} )} ) )
tmp.1.mean = apply(tmp.1, 1:2, mean)
tmp.2.mean = apply(tmp.2, 1:2, mean)
tmp.3.mean = apply(tmp.3, 1:2, mean)
tmp.1.sd = apply(tmp.1, 1:2, quantile, probs = c(0.025,0.975))
plot(seq(-2,2,length.out = 20),tmp.1.mean[2,])
points(seq(-2,2,length.out = 20),tmp.2.mean[2,], col="red")
points(seq(-2,2,length.out = 20),tmp.3.mean[2,], col="blue")
lines(seq(-2,2,length.out = 20), truth.curve[2,])

rowMeans((tmp.1.mean-truth.curve)^2)
rowMeans((tmp.2.mean-truth.curve)^2)
rowMeans((tmp.3.mean-truth.curve)^2)