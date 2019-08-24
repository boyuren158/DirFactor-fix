args = commandArgs(trailingOnly=TRUE)

n = 300
n.sub = 50
p = 100

Y.sub = matrix( rnorm( n.sub*4 ), nrow = 4 )
Y.sub[1:2, (n.sub/2+1):n.sub] = 0
Y.sub[3:4, 1:(n.sub/2)] = 0
X = matrix( rnorm( p*4 ), nrow = 4 )

subject.id = data.frame( id=as.factor(rep(1:n.sub, each = n/n.sub)) )
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

data.out = list( sigma=sigma, y.raw=y.raw, v.fix.inter = v.fix.inter, X = X, Y.sub = Y.sub)
saveRDS(data.out, sprintf("%s/rep%s/sim_data", args[1],args[2]))
