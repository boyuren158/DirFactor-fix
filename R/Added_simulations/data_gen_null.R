args = commandArgs(trailingOnly=TRUE)

n = 300
n.sub = 50
p = 100
tot.count = 1e5

#subject identifier
subject.id = data.frame( id=as.factor(rep(1:n.sub, each = n/n.sub)) )
sub.design = t(model.matrix( ~.-1, subject.id ))
sub.Z = data.frame( rep(1:n.sub, each = n/n.sub) )

sim.res = readRDS( sprintf("%s/rep%s/sim_data", args[1], args[2]) )

sigma = sim.res$sigma
y.raw = sim.res$y.raw
v.fix.inter = sim.res$v.fix.inter
X = sim.res$X
Y.sub = sim.res$Y.sub
y.fix.inter = t(model.matrix( ~ysmall-1, data = y.raw))

Q = t(X)%*%Y.sub%*%sub.design + matrix( rnorm( n*p, sd = sqrt(1) ), nrow = p )

w.DF = Q^2*(Q>0)*sigma
w.DF.n = t(t(w.DF)/colSums(w.DF))

# reads
data.DF = sapply(1:ncol(w.DF.n), function(i) rmultinom(1, tot.count, w.DF.n[,i]))

# save the data
# first for DirFactor
write.table(data.DF, sprintf("%s/rep%s/count.csv", args[1], args[2]), sep=",", row.names = F, col.names = F)
write.table(y.raw, sprintf("%s/rep%s/covariates.csv", args[1], args[2]), sep = ",", row.names = F, col.names = F)
