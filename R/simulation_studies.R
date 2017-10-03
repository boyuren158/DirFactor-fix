source("MCMC.R")
source("utilities.R")

# fixed quantities
n = 300
n.sub = 50
p = 100

subject.id = data.frame( id=as.factor(rep(1:n.sub, each = n/n.sub)) )#data.frame( id = as.factor(sample( rep(1:50, each = 6), size = 300, replace = F )) )
sub.design = t(model.matrix( ~.-1, subject.id ))
v.fix.inter = rbind( c(rep( 5, 8 ), rep( -5, 8 ), rep(0,p-16)), c(rep(5,4), rep(-5,4), rep(5,4), rep(-5,4), rep(0,p-16)),
                     c( rep(c(10,-5,-5,-10),2), rep(c(-10,5,5,10),2), rep(0,p-16) ) )
sigma = rbeta( p, 0.1, 0.4 )

##############################################################################
################## Correctly specified model simulation ######################
##############################################################################
# 50 replicates
for( sim.rep in 1:50 ){
  # simulate data
  Y.sub = matrix( rnorm( n.sub*4 ), nrow = 4 )
  Y.sub[1:2, (n.sub/2+1):n.sub] = 0
  Y.sub[3:4, 1:(n.sub/2)] = 0
  X = matrix( rnorm( p*4 ), nrow = 4 )
  
  y.raw = data.frame( ysmall = rnorm(n), y2 = rep( sample(0:1,n.sub,replace = T), each = n/n.sub ) )
  y.fix.inter = t(model.matrix( ~ysmall+y2+ysmall*y2-1, data = y.raw) )
  hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, 
                sub.design = sub.design, y.fix = y.fix.inter, alpha = 10, beta = 0 )
  Q = t(X)%*%Y.sub%*%sub.design + t(v.fix.inter)%*%y.fix.inter + matrix( rnorm( n*p ), nrow = p )
  weight = sigma*Q^2*(Q>0)
  data = apply( weight, 2, function(x) rmultinom( 1, 1e5, x ) )
  # save the truth
  data.save = list( Q=Q, X=X, Y.sub=Y.sub, sub.design=sub.design, v.fix.inter = v.fix.inter, y.fix.inter=y.fix.inter, sigma = sigma )
  saveRDS( data.save, paste("rep_", sim.rep, "/truth",sep="") )
  
  DirFactor.fix( data, hyper, save.path = paste("rep_", sim.rep, "/sim", sep=""), thin = 50, step = 1e5, step.disp = 1e3 )
}

# retrieve simulations
res.all = lapply( 1:50, function(x){
  folder = paste("rep", x, sep="_")
  mcmc.res.simple = lapply( paste(folder, "/sim_", seq(20050,1e5,500), sep = "" ), readRDS )
  mcmc.truth = readRDS( paste("rep_", x, "/truth", sep = "" ) )
  list( res = mcmc.res.simple, truth = mcmc.truth )
})

# rescaled estimates of v
mt.v.all = lapply( res.all, function(x){
  v.truth = x$truth$v.fix.inter
  Sigma.truth = t(x$truth$Y.sub%*%sub.design)%*%(x$truth$Y.sub%*%sub.design) + diag(rep(1,n))
  sigma.prop = x$truth$sigma/sum(x$truth$sigma)
  v.est = ListtoArray( lapply( x$res, function(y){
    Sigma.est = t(y$Y.sub%*%sub.design)%*%(y$Y.sub%*%sub.design) + diag(rep(y$er,n))
    sigma.prop.est = y$sigma/sum(y$sigma)
    t(t(diag(y$d)%*%y$x/sqrt(sum(diag(Sigma.est))))/sqrt(sigma.prop/sigma.prop.est))
  }) )
  list(truth = v.truth/sqrt(sum(diag(Sigma.truth))), est = v.est )
})
# mse of estimates
mse.v = sapply( mt.v.all, function(x){
  colMeans(x$truth - apply(x$est,1:2,mean))^2
})
# posterior distribution of v in a single replicate
mt.v.ind = mt.v.all[[1]]
plt.data = data.frame( species = as.factor(rep(rep(1:16, each=3),dim(mt.v.ind$est)[3])), v = c(mt.v.ind$est[,1:16,]),
                       v.tru = c(c(mt.v.ind$truth[,1:16]), rep(rep(NA,48),dim(mt.v.ind$est)[3]-1) ),
                       terms = rep(rep(factor(c("Continuous", "Binary", "Interaction"), levels = c("Continuous", "Binary", "Interaction")), 16), dim(mt.v.ind$est)[3]),
                       group = as.factor(rep(1:48,dim(mt.v.ind$est)[3])) )
pdf("v_est_single.pdf", width = 10, height = 4)
ggplot(data = plt.data ) + geom_boxplot(aes(x=species,y=v,group=group), outlier.shape = NA) + geom_point(aes(x=species,y=v.tru)) + facet_grid(.~terms) + theme_bw() + 
  ylab(expression(v[l][i]/sqrt(trace(Sigma))))
dev.off()

##################################################################
# estimates of S between subjects
mt.accu = sapply( res.all, function(x){
  mt.truth = cov2cor( t(x$truth$Y.sub)%*%x$truth$Y.sub + diag(rep(1,ncol(x$truth$Y.sub))) )
  mt.mean = apply( ListtoArray(lapply(x$res, function(x) cov2cor(t(x$Y.sub)%*%x$Y.sub+diag(rep(x$er,ncol(x$Y.sub)))))), 1:2, mean )
  list( est = mt.mean, truth = mt.truth, rv = rv.coef( mt.truth, mt.mean ) )
})
mean(sapply(mt.accu,function(x) x$rv))
sd(sapply(mt.accu,function(x) x$rv))

# plot of S in a single replicate
corr.est.single = mt.accu[[1]]$est
corr.est[lower.tri(corr.est)] = mt.accu[[1]]$truth[lower.tri(corr.trut)]
pdf("S_est_single.pdf", width = 5, height = 5)
gg.heatmap( corr.est ) + geom_vline(xintercept = 25.5,size=1.5) + geom_hline(yintercept = 25.5, size=1.5)
dev.off()

##################################################################
# estimates of derivatives
all.deriv = lapply( res.all, function(x){
  n = ncol(x$truth$Q)
  p = nrow(x$truth$Q)
  local.der.tru = apply( ListtoArray(lapply( 1:10, function(iii){
    Q = t(x$truth$X)%*%x$truth$Y.sub%*%sub.design + t(x$truth$v.fix.inter)%*%x$truth$y.fix.inter + 
      matrix( rnorm(n*p), ncol = n )
    CalcDerivativeAnalytic( Q, x$truth$sigma, x$truth$v.fix.inter, 
                            derivative.mt = rbind( 1, 0, x$truth$y.fix.inter[2,] ) )
  })), 1:2, mean )
  local.der.est.all = ListtoArray( lapply( x$res, function(y){
    apply( ListtoArray(lapply( 1:10, function(iii){
      Q =  t(y$X)%*%y$Y.sub%*%sub.design + t(y$x)%*%diag(y$d)%*%x$truth$y.fix.inter + 
        matrix( rnorm(n*p,sd=sqrt(y$er)), ncol = n )
      CalcDerivativeAnalytic( Q, y$sigma, diag(y$d)%*%y$x, rbind( 1, 0, x$truth$y.fix.inter[2,] ) ) 
    })), 1:2, mean )
  } ) )
  
  list( truth = local.der.tru, est = local.der.est.all )
} )

# MSE of estimated derivatives
mse.model = sapply( all.deriv, function(x){
  est.mean = apply( x$est, 1:2, mean )
  rowMeans( (x$truth - est.mean)^2 )
} )

# visualization of derivatives in a single replicate
plot.deriv = all.deriv[[1]]
plot.param = res.all[[1]]$truth
est.mean = apply( plot.deriv$est, 1:2, mean )
est.ci = apply( plot.deriv$est, 1:2, quantile, prob = c(0.025,0.975) )
weight = plot.param$Q^2*(plot.param$Q>0)*plot.param$sigma
obs.abd = t(t(weight)/colSums(weight))
plt.idx = sample(1:300, 150)

p1 = plot.deriv.simple( w = res.all[[1]]$truth$y.fix.inter[1,], est = apply( all.deriv[[1]]$est, 1:2, mean)[1,], 
                        est.ci[,1,], 
                        all.deriv[[1]]$truth[1,], obs.abd = obs.abd[1,], 
                        w2 = as.factor(plot.param$y.fix.inter[2,]), plt.idx = plt.idx ) + theme_bw()
p2 = plot.deriv.simple( w = res.all[[1]]$truth$y.fix.inter[1,], est = apply( all.deriv[[1]]$est, 1:2, mean)[2,], 
                        est.ci[,2,], 
                        all.deriv[[1]]$truth[2,], obs.abd = obs.abd[2,], 
                        w2 = as.factor(plot.param$y.fix.inter[2,]), plt.idx = plt.idx ) + theme_bw()
p3 = plot.deriv.simple( w = res.all[[1]]$truth$y.fix.inter[1,], est = apply( all.deriv[[1]]$est, 1:2, mean)[9,], 
                        est.ci[,9,], 
                        all.deriv[[1]]$truth[9,], obs.abd = obs.abd[9,], 
                        w2 = as.factor(plot.param$y.fix.inter[2,]), plt.idx = plt.idx ) + theme_bw()

pdf("derivatives_est_single.pdf", width = 16, height = 4)
multiplot( p1, p2, p3, cols = 3)
dev.off()

##################################################################
# estimates of population trends

all.pop.trend = lapply( res.all, function(x){
  # for binary = 0
  tmp = ListtoArray( lapply( x$res, function(y){
    sapply( seq(-2,2,length.out = 20), function(xx){
      y.tmp = rbind( rep(xx,n.sub), 0, 0 )
      get.pop.mean(y.tmp, y)} )} ) )
  # for binary = 1
  tmp.1 = ListtoArray( lapply( x$res, function(y){
    sapply( seq(-2,2,length.out = 20), function(xx){
      y.tmp = rbind( rep(xx,n.sub), 1, rep(xx,n.sub) )
      get.pop.mean(y.tmp, y)} )} ) )
  tmp.mean = apply( tmp, 1:2, mean )
  tmp.quan = apply( tmp, 1:2, quantile, probs = c(0.025, 0.975) )
  
  tmp.mean.1 = apply( tmp.1, 1:2, mean )
  tmp.quan.1 = apply( tmp.1, 1:2, quantile, probs = c(0.025, 0.975) )
  
  # get the truth
  truth.curve = apply( ListtoArray( lapply( 1:1000, function(iii) {sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp = rbind( rep(xx,n.sub), 0, 0 )
    get.pop.mean(y.tmp, list(X=x$truth$X, Y.sub=x$truth$Y.sub, sigma = x$truth$sigma, d = c(1,1,1), x = x$truth$v.fix.inter, er = 1))} )} ) ), 1:2, mean )
  
  truth.curve.1 = apply( ListtoArray( lapply( 1:1000, function(iii) {sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp = rbind( rep(xx,n.sub), 1, rep(xx,n.sub) )
    get.pop.mean(y.tmp, list(X=x$truth$X, Y.sub=x$truth$Y.sub, sigma = x$truth$sigma, d = c(1,1,1), x = x$truth$v.fix.inter, er = 1))} )} ) ), 1:2, mean )
  
  list( est = list(b0 = tmp.mean, b1=tmp.mean.1), truth = list(b0=truth.curve,b1=truth.curve.1), est.ci = list(b0=tmp.quan,b1=tmp.quan.1) )
} )

est.0.inter = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$est$b0
}) ), 1:2, quantile, prob = c(0.025,0.975) )
est.1.inter = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$est$b1
}) ), 1:2, quantile, prob = c(0.025,0.975) )
truth.0 = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$truth$b0
}) ), 1:2, mean )
truth.1 = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$truth$b1
}) ), 1:2, mean )

trend.est0 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$est$b0 )), 1:2, quantile, prob = c(0.025,0.975) )
trend.est1 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$est$b1 )), 1:2, quantile, prob = c(0.025,0.975) )
trend.tru0 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$truth$b0 )), 1:2, mean )
trend.tru1 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$truth$b1 )), 1:2, mean )

# population trends in one replicates and over all replicates
i = 1
p1 = plot.trend.simple( seq(-2,2,length.out = 20), all.pop.trend[[1]]$est.ci$b0[,i,], 
                        all.pop.trend[[1]]$est.ci$b1[,i,],
                        all.pop.trend[[1]]$truth$b0[i,], all.pop.trend[[1]]$truth$b1[i,] )
p1.all = plot.trend.simple( seq(-2,2,length.out = 20), trend.est0[,i,], trend.est1[,i,],
                            trend.tru0[i,], trend.tru1[i,] )

i = 2
p2 = plot.trend.simple( seq(-2,2,length.out = 20), all.pop.trend[[1]]$est.ci$b0[,i,], 
                        all.pop.trend[[1]]$est.ci$b1[,i,],
                        all.pop.trend[[1]]$truth$b0[i,], all.pop.trend[[1]]$truth$b1[i,] )
p2.all = plot.trend.simple( seq(-2,2,length.out = 20), trend.est0[,i,], trend.est1[,i,],
                            trend.tru0[i,], trend.tru1[i,] )

i = 11
p3 = plot.trend.simple( seq(-2,2,length.out = 20), all.pop.trend[[1]]$est.ci$b0[,i,], 
                        all.pop.trend[[1]]$est.ci$b1[,i,],
                        all.pop.trend[[1]]$truth$b0[i,], all.pop.trend[[1]]$truth$b1[i,] )
p3.all = plot.trend.simple( seq(-2,2,length.out = 20), trend.est0[,i,], trend.est1[,i,],
                            trend.tru0[i,], trend.tru1[i,] )

# plots for one replicate
pdf( "population_trends_single.pdf", width = 16, height = 4 )
multiplot( p1, p2, p3, cols = 3 )
dev.off()

# plots for all replicates
pdf( "population_trends_all.pdf", width = 16, height = 4 )
multiplot( p1.all, p2.all, p3.all, cols = 3 )
dev.off()

##############################################################################
################## Misspecified model simulation #############################
##############################################################################
for( sim.rep in 1:50 ){
  # simulate data
  Y.sub = matrix( rnorm( n.sub*4 ), nrow = 4 )
  Y.sub[1:2, (n.sub/2+1):n.sub] = 0
  Y.sub[3:4, 1:(n.sub/2)] = 0
  X = matrix( rnorm( p*4 ), nrow = 4 )
  
  y.raw = data.frame( ysmall = rnorm(n), y2 = rep( sample(0:1,n.sub,replace = T), each = n/n.sub ) )
  y.fix.inter = t(model.matrix( ~ysmall+y2+ysmall*y2-1, data = y.raw) )
  hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, 
                sub.design = sub.design, y.fix = y.fix.inter, alpha = 10, beta = 0 )
  Q = t(X)%*%Y.sub%*%sub.design + t(v.fix.inter)%*%y.fix.inter + matrix( rnorm( n*p ), nrow = p )
  weight = sigma*Q*(Q>0)
  data = apply( weight, 2, function(x) rmultinom( 1, 1e5, x ) )
  # save the truth
  data.save = list( Q=Q, X=X, Y.sub=Y.sub, sub.design=sub.design, v.fix.inter = v.fix.inter, y.fix.inter=y.fix.inter, sigma = sigma )
  saveRDS( data.save, paste("rep_mis_", sim.rep, "/truth",sep="") )
  
  DirFactor.fix( data, hyper, save.path = paste("rep_mis_", sim.rep, "/sim", sep=""), thin = 50, step = 1e5, step.disp = 1e3 )
}

# retrieve simulations
res.all = lapply( 1:50, function(x){
  folder = paste("rep_mis", x, sep="_")
  mcmc.res.simple = lapply( paste(folder, "/sim_", seq(20050,1e5,500), sep = "" ), readRDS )
  mcmc.truth = readRDS( paste(folder, "/truth", sep = "" ) )
  list( res = mcmc.res.simple, truth = mcmc.truth )
})

##################################################################
# estimates of S between subjects
mt.accu = sapply( res.all, function(x){
  mt.truth = cov2cor( t(x$truth$Y.sub)%*%x$truth$Y.sub + diag(rep(1,ncol(x$truth$Y.sub))) )
  mt.mean = apply( ListtoArray(lapply(x$res, function(x) cov2cor(t(x$Y.sub)%*%x$Y.sub+diag(rep(x$er,ncol(x$Y.sub)))))), 1:2, mean )
  list( est = mt.mean, truth = mt.truth, rv = rv.coef( mt.truth, mt.mean ) )
})
mean(sapply(mt.accu,function(x) x$rv))
sd(sapply(mt.accu,function(x) x$rv))

# plot of S in a single replicate
corr.est.single = mt.accu[[1]]$est
corr.est[lower.tri(corr.est)] = mt.accu[[1]]$truth[lower.tri(corr.trut)]
pdf("S_est_single.pdf", width = 5, height = 5)
gg.heatmap( corr.est ) + geom_vline(xintercept = 25.5,size=1.5) + geom_hline(yintercept = 25.5, size=1.5)
dev.off()

##################################################################
# estimates of derivatives
all.deriv = lapply( res.all, function(x){
  n = ncol(x$truth$Q)
  p = nrow(x$truth$Q)
  local.der.tru = apply( ListtoArray(lapply( 1:10, function(iii){
    Q = t(x$truth$X)%*%x$truth$Y.sub%*%sub.design + t(x$truth$v.fix.inter)%*%x$truth$y.fix.inter + 
      matrix( rnorm(n*p), ncol = n )
    CalcDerivativeAnalytic.mis( Q, x$truth$sigma, x$truth$v.fix.inter, 
                                derivative.mt = rbind( 1, 0, x$truth$y.fix.inter[2,] ) )
  })), 1:2, mean )
  local.der.est.all = ListtoArray( lapply( x$res, function(y){
    apply( ListtoArray(lapply( 1:10, function(iii){
      Q =  t(y$X)%*%y$Y.sub%*%sub.design + t(y$x)%*%diag(y$d)%*%x$truth$y.fix.inter + 
        matrix( rnorm(n*p,sd=sqrt(y$er)), ncol = n )
      CalcDerivativeAnalytic( Q, y$sigma, diag(y$d)%*%y$x, rbind( 1, 0, x$truth$y.fix.inter[2,] ) ) 
    })), 1:2, mean )
  } ) )
  
  list( truth = local.der.tru, est = local.der.est.all )
} )

# MSE of estimated derivatives
mse.model = sapply( all.deriv, function(x){
  est.mean = apply( x$est, 1:2, mean )
  rowMeans( (x$truth - est.mean)^2 )
} )

# visualization of derivatives in a single replicate
plot.deriv = all.deriv[[1]]
plot.param = res.all[[1]]$truth
est.mean = apply( plot.deriv$est, 1:2, mean )
est.ci = apply( plot.deriv$est, 1:2, quantile, prob = c(0.025,0.975) )
weight = plot.param$Q^2*(plot.param$Q>0)*plot.param$sigma
obs.abd = t(t(weight)/colSums(weight))
plt.idx = sample(1:300, 150)

p1 = plot.deriv.simple( w = res.all[[1]]$truth$y.fix.inter[1,], est = apply( all.deriv[[1]]$est, 1:2, mean)[1,], 
                        est.ci[,1,], 
                        all.deriv[[1]]$truth[1,], obs.abd = obs.abd[1,], 
                        w2 = as.factor(plot.param$y.fix.inter[2,]), plt.idx = plt.idx ) + theme_bw()
p2 = plot.deriv.simple( w = res.all[[1]]$truth$y.fix.inter[1,], est = apply( all.deriv[[1]]$est, 1:2, mean)[2,], 
                        est.ci[,2,], 
                        all.deriv[[1]]$truth[2,], obs.abd = obs.abd[2,], 
                        w2 = as.factor(plot.param$y.fix.inter[2,]), plt.idx = plt.idx ) + theme_bw()
p3 = plot.deriv.simple( w = res.all[[1]]$truth$y.fix.inter[1,], est = apply( all.deriv[[1]]$est, 1:2, mean)[9,], 
                        est.ci[,9,], 
                        all.deriv[[1]]$truth[9,], obs.abd = obs.abd[9,], 
                        w2 = as.factor(plot.param$y.fix.inter[2,]), plt.idx = plt.idx ) + theme_bw()

pdf("derivatives_est_single_mis.pdf", width = 16, height = 4)
multiplot( p1, p2, p3, cols = 3)
dev.off()

##################################################################
# estimates of population trends

all.pop.trend = lapply( res.all, function(x){
  # for binary = 0
  tmp = ListtoArray( lapply( x$res, function(y){
    sapply( seq(-2,2,length.out = 20), function(xx){
      y.tmp = rbind( rep(xx,n.sub), 0, 0 )
      get.pop.mean(y.tmp, y)} )} ) )
  # for binary = 1
  tmp.1 = ListtoArray( lapply( x$res, function(y){
    sapply( seq(-2,2,length.out = 20), function(xx){
      y.tmp = rbind( rep(xx,n.sub), 1, rep(xx,n.sub) )
      get.pop.mean(y.tmp, y)} )} ) )
  tmp.mean = apply( tmp, 1:2, mean )
  tmp.quan = apply( tmp, 1:2, quantile, probs = c(0.025, 0.975) )
  
  tmp.mean.1 = apply( tmp.1, 1:2, mean )
  tmp.quan.1 = apply( tmp.1, 1:2, quantile, probs = c(0.025, 0.975) )
  
  # get the truth
  truth.curve = apply( ListtoArray( lapply( 1:1000, function(iii) {sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp = rbind( rep(xx,n.sub), 0, 0 )
    get.pop.mean.mis(y.tmp, list(X=x$truth$X, Y.sub=x$truth$Y.sub, sigma = x$truth$sigma, d = c(1,1,1), x = x$truth$v.fix.inter, er = 1))} )} ) ), 1:2, mean )
  
  truth.curve.1 = apply( ListtoArray( lapply( 1:1000, function(iii) {sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp = rbind( rep(xx,n.sub), 1, rep(xx,n.sub) )
    get.pop.mean.mis(y.tmp, list(X=x$truth$X, Y.sub=x$truth$Y.sub, sigma = x$truth$sigma, d = c(1,1,1), x = x$truth$v.fix.inter, er = 1))} )} ) ), 1:2, mean )
  
  list( est = list(b0 = tmp.mean, b1=tmp.mean.1), truth = list(b0=truth.curve,b1=truth.curve.1), est.ci = list(b0=tmp.quan,b1=tmp.quan.1) )
})

est.0.inter = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$est$b0
}) ), 1:2, quantile, prob = c(0.025,0.975) )
est.1.inter = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$est$b1
}) ), 1:2, quantile, prob = c(0.025,0.975) )
truth.0 = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$truth$b0
}) ), 1:2, mean )
truth.1 = apply( ListtoArray( lapply( all.pop.trend, function(x){
  x$truth$b1
}) ), 1:2, mean )

trend.est0 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$est$b0 )), 1:2, quantile, prob = c(0.025,0.975) )
trend.est1 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$est$b1 )), 1:2, quantile, prob = c(0.025,0.975) )
trend.tru0 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$truth$b0 )), 1:2, mean )
trend.tru1 = apply( ListtoArray(lapply( all.pop.trend, function(x) x$truth$b1 )), 1:2, mean )

# population trends in one replicates and over all replicates
i = 1
p1 = plot.trend.simple( seq(-2,2,length.out = 20), all.pop.trend[[1]]$est.ci$b0[,i,], 
                        all.pop.trend[[1]]$est.ci$b1[,i,],
                        all.pop.trend[[1]]$truth$b0[i,], all.pop.trend[[1]]$truth$b1[i,] )
p1.all = plot.trend.simple( seq(-2,2,length.out = 20), trend.est0[,i,], trend.est1[,i,],
                            trend.tru0[i,], trend.tru1[i,] )

i = 2
p2 = plot.trend.simple( seq(-2,2,length.out = 20), all.pop.trend[[1]]$est.ci$b0[,i,], 
                        all.pop.trend[[1]]$est.ci$b1[,i,],
                        all.pop.trend[[1]]$truth$b0[i,], all.pop.trend[[1]]$truth$b1[i,] )
p2.all = plot.trend.simple( seq(-2,2,length.out = 20), trend.est0[,i,], trend.est1[,i,],
                            trend.tru0[i,], trend.tru1[i,] )

i = 11
p3 = plot.trend.simple( seq(-2,2,length.out = 20), all.pop.trend[[1]]$est.ci$b0[,i,], 
                        all.pop.trend[[1]]$est.ci$b1[,i,],
                        all.pop.trend[[1]]$truth$b0[i,], all.pop.trend[[1]]$truth$b1[i,] )
p3.all = plot.trend.simple( seq(-2,2,length.out = 20), trend.est0[,i,], trend.est1[,i,],
                            trend.tru0[i,], trend.tru1[i,] )

# plots for one replicate
pdf( "population_trends_single_mis.pdf", width = 16, height = 4 )
multiplot( p1, p2, p3, cols = 3 )
dev.off()

# plots for all replicates
pdf( "population_trends_all_mis.pdf", width = 16, height = 4 )
multiplot( p1.all, p2.all, p3.all, cols = 3 )
dev.off()