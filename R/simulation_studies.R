source("MCMC.R")

n = 300
n.sub = 50
p = 100

Y.sub = matrix( rnorm( n.sub*4 ), nrow = 4 )
Y.sub[1:2, (n.sub/2+1):n.sub] = 0
Y.sub[3:4, 1:(n.sub/2)] = 0
X = matrix( rnorm( p*4 ), nrow = 4 )

subject.id = data.frame( id=as.factor(rep(1:n.sub, each = n/n.sub)) )#data.frame( id = as.factor(sample( rep(1:50, each = 6), size = 300, replace = F )) )
sub.design = t(model.matrix( ~.-1, subject.id ))

y.raw = data.frame( ysmall = rnorm(n), y2 = rep( sample(0:1,n.sub,replace = T), each = n/n.sub ) )
y.fix.inter = t(model.matrix( ~ysmall+y2+ysmall*y2-1, data = y.raw) )
v.fix.inter = rbind( c(rep( 5, 8 ), rep( -5, 8 ), rep(0,p-16)), c(rep(5,4), rep(-5,4), rep(5,4), rep(-5,4), rep(0,p-16)),
                     c( rep(c(10,-5,-5,-10),2), rep(c(-10,5,5,10),2), rep(0,p-16) ) )

hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, 
              sub.design = sub.design, y.fix = y.fix.inter, alpha = 10, beta = 0 )
sigma.normal = rbeta( p, 1/2-10/p, 1-10/p )

Q = t(X)%*%Y.sub%*%sub.design + t(v.fix.inter)%*%y.fix.inter + matrix( rnorm( n*p ), nrow = p )
weight = sigma.normal*Q^2*(Q>0)
data = apply( weight, 2, function(x) rmultinom( 1, 1e5, x ) )

start.sub = list( er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ), 
                  sigma = sample( sigma.value, size = p, replace = T, prob = sigma.prior ), 
                  T.aug = rgamma( n, 10, 1 ),
                  Q = matrix( 0.5, nrow = p, ncol = n ),
                  X = matrix( rnorm( hyper$m*p ), nrow = hyper$m ),
                  Y.sub = matrix( 0, nrow = hyper$m, ncol = n.sub ),
                  delta = c( rgamma( 1, shape = hyper$a1, rate = 1 ), rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) ),
                  phi = matrix( rgamma( hyper$m*n.sub, shape = 3/2, rate = 3/2 ), nrow = n.sub ),
                  y = y.fix.inter.large, 
                  x = matrix( rnorm( p*nrow(y.fix.inter) ), nrow=nrow(y.fix.inter) ),
                  d = rgamma( nrow(y.fix.inter), 1, 1 ) )

source("MCMC.R")
test = DirFactor.fix( data, hyper, step = 1e2 )
tmp = readRDS("/private/var/folders/n1/q055hxqj61s_54c4g9lgvsg80000gn/T/RtmpcHZyeL/sim_100")
nod.free.mcmc( data, start.sub, hyper, sigma.value, sigma.prior, "C:/Users/Boyu/Desktop/simple_test_inter_nointercept_normal_sig_large/sim", burnin = 0.2, step = 1e5, thin = 50 )

mcmc.res.simple = lapply( paste("C:/Users/Boyu/Desktop/simple_test_inter_nointercept_normal_sig/sim", seq(20049,99999,50), sep = "_" ), readRDS )
corr.simple = apply( ListtoArray(lapply(mcmc.res.simple, function(x) cov2cor(t(x$Y.sub)%*%x$Y.sub+diag(rep(x$er,ncol(x$Y.sub)))))), 1:2, mean )

gg.heatmap( corr.simple )
gg.heatmap( cov2cor(t(Y.sub)%*%Y.sub+diag(ncol(Y.sub))) )

plt.corr = cov2cor(t(Y.sub)%*%Y.sub+diag(ncol(Y.sub)))
plt.corr[upper.tri(plt.corr)] = corr.simple[upper.tri(corr.simple)]
p.corr.comp = gg.heatmap( plt.corr, y.lab = "Estimated", x.lab = "Truth" ) + 
  geom_vline(xintercept = n.sub/2+0.5, size = 1.5) + 
  geom_hline(yintercept = n.sub/2+0.5, size = 1.5)

ggsave("subject/sim_sub_corr.pdf", p.corr.comp, width = 6, height = 5)

# now consider the estimates of x
plt.data = c()
for( iii in 1:16 ){
  ttt = sapply( mcmc.res.simple, function(x){
    # ratio = sqrt(mean( diag(t(Y.sub)%*%Y.sub+diag(n.sub))/diag(t(x$Y.sub)%*%x$Y.sub+diag(rep(x$er,ncol(x$Y.sub)))) ))
    # x$x[,iii]*x$d*ratio
    x$d*x$x[,iii]/sqrt(sum(x$d^2*x$x[,iii]^2))*sqrt(sum(v.fix.inter[,iii]^2))
  } )
  plt.data = rbind( plt.data, cbind( c(t(ttt)), c(do.call(rbind, c(list(v.fix.inter[,iii]),rep(list(c(NA,NA,NA)),ncol(ttt)-1)))), iii, 
                                     rep(1:3,each=ncol(ttt)) ) )
}
colnames(plt.data) = c("x", "x.tru", "species", "index")
plt.data = as.data.frame(plt.data)
plt.data$index = factor(plt.data$index, labels = c("Continuous","Binary","Interaction"))
plt.data$species = as.factor(plt.data$species)
p.coef = ggplot() + geom_boxplot( data = plt.data, aes(y=x, x=species), outlier.colour = NA ) + geom_point( data = plt.data, aes(y=x.tru,x=species)) +
  facet_grid(.~index) + ylab("v.fix") + theme_bw()
ggsave( "subject/sim_reg_coef.pdf", p.coef, width = 10, height = 4 )

# let's check the population average trajectory
get.pop.mean = function( y.fix, sim.ind ){
  Q = t(sim.ind$X)%*%sim.ind$Y.sub + t(sim.ind$x)%*%diag(sim.ind$d)%*%y.fix  + 
    matrix( rnorm( ncol(sim.ind$X)*ncol(sim.ind$Y.sub), sd = sqrt(sim.ind$er) ), nrow = ncol(sim.ind$X) )
  weight = sim.ind$sigma*Q^2*(Q>0)
  apply( t(weight)/colSums(weight), 2, mean )
}

# for binary = 0
tmp = ListtoArray( lapply( mcmc.res.simple, function(x){
  sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp = rbind( rep(xx,n.sub), 0, 0 )
    get.pop.mean(y.tmp, x)} )} ) )
# for binary = 1
tmp.1 = ListtoArray( lapply( mcmc.res.simple, function(x){
  sapply( seq(-2,2,length.out = 20), function(xx){
    y.tmp = rbind( rep(xx,n.sub), 1, rep(xx,n.sub) )
    get.pop.mean(y.tmp, x)} )} ) )
tmp.mean = apply( tmp, 1:2, mean )
tmp.quan = apply( tmp, 1:2, quantile, probs = c(0.025, 0.975) )

tmp.mean.1 = apply( tmp.1, 1:2, mean )
tmp.quan.1 = apply( tmp.1, 1:2, quantile, probs = c(0.025, 0.975) )

# get the truth
truth.curve = apply( ListtoArray( lapply( 1:1000, function(iii) {sapply( seq(-2,2,length.out = 20), function(xx){
  y.tmp = rbind( rep(xx,n.sub), 0, 0 )
  get.pop.mean(y.tmp, list(X=X, Y.sub=Y.sub, sigma = sigma.normal, d = c(1,1,1), x = v.fix.inter, er = 1))} )} ) ), 1:2, mean )

truth.curve.1 = apply( ListtoArray( lapply( 1:1000, function(iii) {sapply( seq(-2,2,length.out = 20), function(xx){
  y.tmp = rbind( rep(xx,n.sub), 1, rep(xx,n.sub) )
  get.pop.mean(y.tmp, list(X=X, Y.sub=Y.sub, sigma = sigma.normal, d = c(1,1,1), x = v.fix.inter, er = 1))} )} ) ), 1:2, mean )

pdf("subject/population_mean_inter.pdf", width = 6, height = 5)
for( idx in 1:p ){
  plt.data = data.frame( x = rep(seq(-2,2,length.out = 20),2),
                         y = c( tmp.mean[idx,], tmp.mean.1[idx,]),
                         y.tru = c(truth.curve[idx,],truth.curve.1[idx,]),
                         ymin = c( tmp.quan[1,idx,], tmp.quan.1[1,idx,]),
                         ymax = c( tmp.quan[2,idx,], tmp.quan.1[2,idx,]),
                         group = as.factor(rep(c(0,1), each = 20)) )
  print(ggplot() + geom_point( data = plt.data, aes(x=x,y=y.tru,color=group,shape = "Truth") ) + 
          geom_line( data = plt.data, aes(x=x,y=y,color=group) ) +
          geom_ribbon( data = plt.data, aes(x=x,ymin=ymin,ymax=ymax,fill=group), alpha = 0.2 ) +
          xlab("covariate") + ylab("Relative abundance") + ggtitle(sprintf("Species %d", idx)) + theme_bw())
}
dev.off()

# local derivative plot
local.der.tru = apply( ListtoArray(lapply( 1:10, function(iii){
  Q = t(X)%*%Y.sub%*%sub.design + t(v.fix.inter)%*%y.fix.inter + 
    matrix( rnorm(n*p), ncol = n )
  CalcDerivativeAnalytic( Q, sigma.normal, v.fix.inter, 
                          derivative.mt = rbind( 1, 0, y.fix.inter[2,] ) )
})), 1:2, mean )

local.der.est = lapply( mcmc.res.simple, function(x){
  #apply( ListtoArray(lapply( 1:1, function(iii){
  Q =  t(x$X)%*%x$Y.sub%*%sub.design + t(x$x)%*%diag(x$d)%*%y.fix.inter + 
    matrix( rnorm(n*p,sd=sqrt(x$er)), ncol = n )
  CalcDerivativeAnalytic( Q, x$sigma, diag(x$d)%*%x$x, rbind( 1, 0, y.fix.inter[2,] ) ) 
  #  })), 1:2, mean )
} )

local.der.est.mean = apply( ListtoArray(local.der.est), 1:2, mean )
local.der.est.ci = apply( ListtoArray(local.der.est), 1:2, quantile, probs = c(0.025,0.975))

plt.idx = sample(1:n, 150)
pdf("subject/sim_population_mean_derivative.pdf", width = 6, height = 5)
for( i in 1:p ){
  plt.data = data.frame(x=rep(y.fix.inter[1,plt.idx],2), y = c(local.der.est.mean[i,plt.idx], local.der.tru[i,plt.idx]),
                        obs.abd = rep( data[i,plt.idx]/1e5, 2 ),
                        group = rep(as.factor(y.fix.inter[2,plt.idx]),2), type = rep(c("Est.", "Truth"), each = length(plt.idx)) )
  plt.band = data.frame( x = y.fix.inter[1,plt.idx], ymin = local.der.est.ci[1,i,plt.idx], ymax = local.der.est.ci[2,i,plt.idx] )
  
  print( ggplot() + geom_point( data=plt.data, aes(x=x,y=y,color=group,shape=type,size=obs.abd) ) + 
           geom_errorbar( data = plt.band, aes(x=x,ymin=ymin,ymax=ymax), size = 0.25 ) +
           scale_shape_manual(values=c(19,4)) + theme_bw() + xlab("covariate") + ylab("Relative abundance") + 
           ggtitle( sprintf("Species %d", i)))
}