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

y.raw = data.frame( ysmall = rnorm(n), ylarge = rnorm(n,sd = 10), y2 = rep( sample(0:1,n.sub,replace = T), each = n/n.sub ) )
y.fix = t(model.matrix( ~ysmall+y2, data = y.raw) )
y.fix.large = t(model.matrix( ~ylarge+y2-1, data = y.raw) )
y.fix.inter = t(model.matrix( ~ysmall+y2+ysmall*y2-1, data = y.raw) )
y.fix.inter.large = t(model.matrix( ~ylarge+y2+ylarge*y2-1, data = y.raw))
v.fix = rbind( rnorm(p), c(rep( 5, 8 ), rep( -5, 8 ), rep(0,p-16)), c(rep(5,4), rep(-5,4), rep(5,4), rep(-5,4), rep(0,p-16)) )
v.fix.inter = rbind( c(rep( 5, 8 ), rep( -5, 8 ), rep(0,p-16)), c(rep(5,4), rep(-5,4), rep(5,4), rep(-5,4), rep(0,p-16)),
                     c( rep(c(10,-5,-5,-10),2), rep(c(-10,5,5,10),2), rep(0,p-16) ) )
v.fix.inter.small = rbind( c(rep( 1, 8 ), rep( -1, 8 ), rep(0,p-16)), c(rep(1,4), rep(-1,4), rep(1,4), rep(-1,4), rep(0,p-16)),
                           c( rep(c(2,0,-1,-2),2), rep(c(-2,0,1,2),2), rep(0,p-16) ) )

hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, 
              sub.design = sub.design, 
              y.fix = y.fix.inter,
              alpha = 10, beta = 0 )
sigma.value = seq(0.001,0.999,0.001)
tmp = c(0,pbeta( sigma.value, 10/68, 1/2-10/68 ))
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)

sigma = sample( sigma.value, p, replace = T, prob = sigma.prior )
sigma.normal = rbeta( p, 1/2, 1 )

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
dev.off()

# now try to see the difference between two groups

covariate.cut = cut( y.fix.inter[1,], 
                     breaks = c(0, quantile( y.fix.inter[1,], probs = seq(0.2,0.8,0.2) ), 
                                max(y.fix.inter[1,])) )
y.linear.0 = rbind( y.fix.inter[1,], 0, 0 )
y.linear.1 = rbind( y.fix.inter[1,], 1, y.fix.inter[1,] )

linear.diff.truth = calc.fix.covariate( list(Q=Q,x=v.fix.inter,sigma=sigma), y.fix.inter, y.linear.1 )-
  calc.fix.covariate( list(Q=Q,x=v.fix.inter,sigma=sigma), y.fix.inter, y.linear.0 )
linear.diff.seg = aggregate(t(linear.diff.truth), by = list(covariate.cut), FUN = mean )

linear.diff.all = lapply( mcmc.res.simple, function(x){
  list( Fin = calc.fix.covariate( x, y.fix.inter, y.linear.0),
        EST = calc.fix.covariate( x, y.fix.inter, y.linear.1),
        RUS = calc.fix.covariate( x, y.fix.inter, y.linear.0) )
} )

plot.fix.covariate( y.fix.inter[1,], linear.diff.all, 1:p, "subject/sim_binary.pdf", linear.diff.seg[,-1] )

########################## Real data analysis ##############################

raw.data = read.csv("../../data/diabimmune/diabimmune_genus.txt", sep = "\t", row.names = 1 )
load( "../../data/diabimmune/DIABIMMUNE_Karelia_metadata_full.RData" )
metadata.wgs = metadata[(!is.na(metadata$mgx_reads_filtered))&(!is.na(metadata$gid_wgs)),]
row.names(metadata.wgs) = metadata.wgs$gid_wgs
use.data = raw.data[,names(raw.data)%in%row.names(metadata.wgs)]
#mgx_reads_filtered is the number of reads/1e6
metadata.use = metadata.wgs[names(use.data),]
use.data.count = t(round(t(use.data/100)*metadata.use$mgx_reads_filtered*1e6))

#get rid of NA in seroconverted
use.data.count = use.data.count[,!is.na(metadata.use$seroconverted)]
metadata.use = metadata.use[!is.na(metadata.use$seroconverted),]

get.dbdata.ready = function( use.data.count, metadata.use, idx, use.formula ){
  metadata.spline = metadata.use[idx,]
  db.sub.mt = t(model.matrix(~subjectID-1,data=metadata.spline))
  age.raw = metadata.spline$age_at_collection
  age.stand = (metadata.spline$age_at_collection-mean(metadata.spline$age_at_collection))/sd(metadata.spline$age_at_collection)
  metadata.spline$age.stand = age.stand
  db.design.linear.mt = t(model.matrix(use.formula, data=metadata.spline))
  
  use.data.spline = use.data.count[,idx]
  use.data.spline= use.data.spline[rowMeans(use.data.spline>0)>0.1,]
  use.data.spline.norm = t(t(use.data.spline)/colSums(use.data.spline))
  
  all.names = sapply( row.names(use.data.spline), function(x) gsub("[kpcofg]__", "", strsplit(x,split="|",fixed=T)[[1]] ) )
  
  return( list( sub.mt = db.sub.mt, design.mt = db.design.linear.mt, 
                data.mt = use.data.spline, data.mt.norm = use.data.spline.norm,
                raw.metadata = metadata.spline, all.taxa = all.names ) )
}

data.full = get.dbdata.ready( use.data.count, metadata.use, 1:ncol(use.data.count),
                              ~country+age.stand+seroconverted+age.stand*country )
get.mcmc.start = function( data.all ){
  hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, sub.design = data.all$sub.mt )
  start = list( er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ), 
                sigma = sample( sigma.value, size = nrow(data.all$data.mt), replace = T, prob = sigma.prior ), 
                T.aug = rgamma( ncol(data.all$data.mt), 10, 1 ),
                Q = matrix( 0.5, nrow = nrow(data.all$data.mt), ncol = ncol(data.all$data.mt) ),
                X = matrix( rnorm( hyper$m*nrow(data.all$data.mt) ), nrow = hyper$m ),
                Y.sub = matrix( 0, nrow = hyper$m, ncol = nrow(hyper$sub.design) ),
                delta = c( rgamma( 1, shape = hyper$a1, rate = 1 ), rgamma( hyper$m-1, shape = hyper$a2, rate = 1 ) ),
                phi = matrix( rgamma( hyper$m*nrow(hyper$sub.design), shape = 3/2, rate = 3/2 ), nrow = nrow(hyper$sub.design) ),
                y = data.all$design.mt, 
                x = matrix( rnorm( nrow(data.all$data.mt)*nrow(data.all$design.mt) ), nrow=nrow(data.all$design.mt) ),
                d = rgamma( nrow(data.all$design.mt), 1, 1 ) )
  return( list( hyper = hyper, start = start ) )
}

# try the full model
start.full = get.mcmc.start( data.full )
nod.free.mcmc( data.full$data.mt, start.full$start, start.full$hyper, sigma.value, sigma.prior, 
               "~/Research/diabimmune/db_full_med_full/sim", burnin = 0.2, 
               step = 1e5, thin = 50 )

#try to see the stratified results
#adding noise on
get.pop.mean = function( y.fix, sim.ind ){
  Q.mean = t(sim.ind$X)%*%sim.ind$Y.sub + t(sim.ind$x)%*%diag(sim.ind$d)%*%y.fix
  Q = Q.mean + matrix( rnorm(length(Q.mean), sd = sqrt(sim.ind$er)), nrow = nrow(Q.mean) )
  weight = sim.ind$sigma*Q^2*(Q>0)
  apply( t(weight)/colSums(weight), 2, mean )
}

plot.pop.curve = function( raw.x, mcmc.res, data.all, fig.path, ... ){
  sub.traj.FIN = ListtoArray(lapply(mcmc.res, function(x){
    sapply( seq( min(raw.x), max(raw.x), length.out = 20 ), 
            function(xx) get.pop.mean(rbind(1, 0, 0, rep(xx,nrow(data.all$sub.mt)), 0, 0, 0 ),x) )
  }))
  sub.traj.EST = ListtoArray(lapply(mcmc.res, function(x){
    sapply( seq( min(raw.x), max(raw.x), length.out = 20 ), 
            function(xx) get.pop.mean(rbind(1, 1, 0, rep(xx,nrow(data.all$sub.mt)), 0, rep(xx,nrow(data.all$sub.mt)), 0 ),x) )
  }))
  sub.traj.RUS = ListtoArray(lapply(mcmc.res, function(x){
    sapply( seq( min(raw.x), max(raw.x), length.out = 20 ), 
            function(xx) get.pop.mean(rbind(1, 0, 1, rep(xx,nrow(data.all$sub.mt)), 0, 0, rep(xx,nrow(data.all$sub.mt)) ),x) )
  }))
  
  # let's try to make the population average plot
  sub.traj.FIN.mean = apply( sub.traj.FIN, 1:2, mean )
  sub.traj.FIN.interv = apply( sub.traj.FIN, 1:2, quantile, probs = c(0.025,0.975) )
  sub.traj.EST.mean = apply( sub.traj.EST, 1:2, mean )
  sub.traj.EST.interv = apply( sub.traj.EST, 1:2, quantile, probs = c(0.025,0.975) )
  sub.traj.RUS.mean = apply( sub.traj.RUS, 1:2, mean )
  sub.traj.RUS.interv = apply( sub.traj.RUS, 1:2, quantile, probs = c(0.025,0.975) )
  
  pdf( fig.path, ... )
  for( i in 1:nrow(data.all$data.mt) ){
    plt.data = data.frame( x = rep( seq( min(raw.x), max(raw.x), length.out = 20 )*sd(data.all$raw.metadata$age_at_collection)+mean(data.all$raw.metadata$age_at_collection), 3),
                           y = c(sub.traj.FIN.mean[i,],sub.traj.EST.mean[i,],sub.traj.RUS.mean[i,]),
                           ymin = c(sub.traj.FIN.interv[1,i,],sub.traj.EST.interv[1,i,],sub.traj.RUS.interv[1,i,]),
                           ymax = c(sub.traj.FIN.interv[2,i,],sub.traj.EST.interv[2,i,],sub.traj.RUS.interv[2,i,]),
                           country = rep( c("FIN","EST","RUS"), each = 20 ) )
    plt.raw = data.frame( x = data.all$raw.metadata$age_at_collection, y = data.all$data.mt.norm[i,], country =  data.all$raw.metadata$country )
    
    print(ggplot() + geom_line( data = plt.data, aes( x = x, y = y, color = country ) ) + 
            geom_ribbon( data = plt.data, aes(x=x,ymin=ymin,ymax=ymax,fill=country),alpha = 0.2 ) + 
      theme_bw() + ggtitle( data.all$all.taxa[6,i] ) + 
        geom_point( data = plt.raw, aes(x=x,y=y, color=country)) + xlab("Age (days)") + ylab("Relative abundance") +
      theme( axis.title = element_text(size=14), axis.text = element_text(size=12) ))
  }
  dev.off()
}

db.res.full = lapply( paste("~/../Desktop/db_full_med_full/sim", seq(20049,99999,50), sep = "_"),
                      function(x){
                        res = readRDS(x)
                        res} )

plot.pop.curve( data.full$design.mt[4,], db.res.full, data.full, "~/full_strat_trend.pdf", width = 6, height = 5 )

# dendrogram

genus.names = data.full$all.taxa[6,]

library( phyloseq )
genus.tree = read_tree("../../data/genus_tree.tre.txt")
genus.tree.label = genus.tree$tip.label
genus.name.rematched = sapply( genus.names, function(x){
  idx = grep( paste("_",x,sep=""), genus.tree.label, ignore.case = T )
  genus.tree.label[idx[1]] } )
genus.name.nomatch = genus.name.rematched[is.na( genus.name.rematched )]

genus.otu.table = diag( length(genus.tree$tip.label) )
row.names(genus.otu.table) = genus.tree.label
colnames(genus.otu.table) = genus.tree.label
phyloseq.otu = otu_table( genus.otu.table, taxa_are_rows = T )
phyloseq.genus = phyloseq( phyloseq.otu )
phyloseq.genus = merge_phyloseq( phyloseq.genus, genus.tree )

phyloseq.final = prune_taxa( genus.name.rematched[!is.na(genus.name.rematched)],
                             phyloseq.genus )
phyloseq.final = prune_samples( sample_sums(phyloseq.final)>0, phyloseq.final )

phylo.dist = UniFrac( phyloseq.final, weighted = F )
phylo.dendro = hclust( phylo.dist )
phylo.dendro$labels = sapply( phylo.dendro$labels, function(x) strsplit(x,"_")[[1]][2] )
phylo.ord = cmdscale(phylo.dist,k=43,eig = T)
phylo.sim = cov2cor((phylo.ord$points%*%t(phylo.ord$points)))

colnames(phylo.sim) = sapply( strsplit(colnames(phylo.sim), split = "_", fixed=T), function(x) x[2] )
row.names(phylo.sim) = colnames(phylo.sim)

dendro.labels.ordered = phylo.dendro$labels[phylo.dendro$order]
phylo.sim = phylo.sim[dendro.labels.ordered,dendro.labels.ordered]

spe.corr.all = lapply( db.res.full, function(x) 
  cov2cor(t(x$X)%*%(x$X)) )
  #cor(x$X) )
  #cor(t(x$Q-t(x$x)%*%diag(x$d)%*%data.full$design.mt)) )
  #cov2cor((x$Q-t(x$x)%*%diag(x$d)%*%data.full$design.mt)%*%t(x$Q-t(x$x)%*%diag(x$d)%*%data.full$design.mt)) )
#spe.corr.inter = apply( apply( ListtoArray(spe.corr.all), 1:2, quantile, probs = c(0.25,0.75) ), 2:3, function(x) prod(x)>0 )
spe.corr.mean = apply( ListtoArray(spe.corr.all), 1:2, mean )
#spe.corr.mean = spe.corr.mean*spe.corr.inter

spe.corr.sub = spe.corr.all[sample(1:length(spe.corr.all),800,replace = T)]
spe.dist.sub = lapply( spe.corr.all, function(x) 1-x )
spe.corr.raw = cor( t( data.full$data.mt.norm ) )
spe.corr.fix = apply( ListtoArray(lapply(db.res.full, function(x)
  #cov2cor(t(x$x)%*%diag(x$d)%*%diag(x$d)%*%x$x))), 1:2, mean )
  #cor(diag(x$d)%*%x$x))), 1:2, mean )
  #cor(t(x$Q)))), 1:2, mean )
  #cor(t(x$Q-t(x$X)%*%x$Y.sub%*%data.full$sub.mt)))), 1:2, mean )
#cov2cor((x$Q-t(x$X)%*%x$Y.sub%*%data.full$sub.mt)%*%t(x$Q-t(x$X)%*%x$Y.sub%*%data.full$sub.mt)))), 1:2, mean )
#cor(t(t(x$x)%*%diag(x$d)%*%data.full$design.mt+t(x$X)%*%x$Y.sub%*%data.full$sub.mt))) ), 1:2, mean )

tmp = hclust(as.dist(1-plot.genus.corr.fix))
par(mfrow=c(2,1))
plot(tmp)
plot(phylo.dendro)


row.names(spe.corr.mean) = genus.names
colnames(spe.corr.mean) = genus.names
plot.genus.corr.mean = spe.corr.mean[dendro.labels.ordered,dendro.labels.ordered]

row.names(spe.corr.raw) = genus.names
colnames(spe.corr.raw) = genus.names
plot.genus.corr.raw = spe.corr.raw[dendro.labels.ordered,dendro.labels.ordered]

row.names(spe.corr.fix) = genus.names
colnames(spe.corr.fix) = genus.names
plot.genus.corr.fix = spe.corr.fix[dendro.labels.ordered,dendro.labels.ordered]

plot.genus.corr.mean[lower.tri(plot.genus.corr.mean)] = phylo.sim[lower.tri(plot.genus.corr.mean)]
plot.genus.corr.raw[lower.tri(plot.genus.corr.mean)] = phylo.sim[lower.tri(plot.genus.corr.mean)]
plot.genus.corr.fix[lower.tri(plot.genus.corr.mean)] = phylo.sim[lower.tri(plot.genus.corr.mean)]

p1 = gg.heatmap(plot.genus.corr.mean)
p2 = gg.heatmap(plot.genus.corr.raw)
p3 = gg.heatmap(plot.genus.corr.fix)
multiplot(p1,p2,p3,cols=3)

rv.coef( plot.genus.corr.mean, phylo.sim )
rv.coef( plot.genus.corr.raw, phylo.sim )
rv.coef( plot.genus.corr.fix, phylo.sim )


row.names( plot.genus.corr.mean ) = NULL
colnames( plot.genus.corr.mean ) = NULL
dendro.labels.num = sapply( dendro.labels.ordered, function(x) which(genus.names==x) )

library(ggplot2)
library(reshape2)
library(ggdendro)

phylo.dendro = as.dendrogram( phylo.dendro )
ddata_xy = dendro_data( phylo.dendro )

theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)

#try to get rid of noname and unclassified
p1 = PlotStatis( lapply(spe.corr.sub,function(x) x[dendro.labels.num,dendro.labels.num]), types = data.full$all.taxa[2,dendro.labels.num], 
                 labels = data.full$all.taxa[6,dendro.labels.num], levels.min = T, size = 1.5 )
ggsave("tmp.pdf", p1[[1]], width = 10, height = 8 )

mt.rnd = spe.corr.mean[dendro.labels.num,dendro.labels.num]
mt.raw = spe.corr.raw[dendro.labels.num,dendro.labels.num]
mt.fix = spe.corr.fix[dendro.labels.num,dendro.labels.num]


plt.rvr = mt.rnd
plt.rvr[upper.tri(plt.rvr)] = mt.raw[upper.tri(mt.raw)]
plt.rvf = mt.rnd
plt.rvf[upper.tri(plt.rvf)] = mt.fix[upper.tri(mt.fix)]

row.names( plt.rvr ) = dendro.labels.ordered
colnames( plt.rvr ) = dendro.labels.ordered
row.names( plt.rvf ) = dendro.labels.ordered
colnames( plt.rvf ) = dendro.labels.ordered

p1.1 = gg.heatmap( t(plt.rvr), y.ticks = dendro.labels.ordered, x.ticks = dendro.labels.ordered, y.lab = "Correlation of random effect", 
                   x.lab = "Correlation of raw data" )
p2.1 = gg.heatmap( t(plt.rvf), y.ticks = dendro.labels.ordered, x.ticks = dendro.labels.ordered, y.lab = "Correlation of random effect", 
                   x.lab = "Correlation of fixed effect" )
pdf("subject/spe_corr_all.pdf", width = 16, height = 7 )
multiplot(p1.1,p2.1,cols =2)
dev.off()

p1 <- gg.heatmap( mt.rnd, y.ticks = dendro.labels.ordered ) + 
  theme(legend.position = "left")
p2 <- ggplot(segment(ddata_xy)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  theme_none + theme( plot.margin = unit(c(0,0,-1,0), 'cm') )

p3 <- p2 + coord_flip() + theme( plot.margin = unit(c(-1,-1,-1,-1), 'cm'))

pdf("subject/dendro_full.pdf", width = 8, height = 7)
grid.newpage()
print(p1, vp=viewport(0.8, 0.8, x=0.35, y=0.4))
print(p2, vp=viewport(0.515, 0.2, x=0.45, y=0.87))
print(p3, vp=viewport(0.15, 0.575, x=0.789, y=0.44))
dev.off()


tmp = PlotClustering( spe.dist.sub )
tmp = tmp + scale_x_continuous(breaks=1:length(genus.names),labels=genus.names, expand=c(0,0))+
  scale_y_continuous(breaks=1:length(genus.names),labels=genus.names, expand=c(0,0))+
  theme( axis.title = element_text(face="bold"), 
         plot.title = element_text(face="bold"),
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=6),
         axis.text.y = element_text(size=6)) + coord_fixed() + xlab("") + ylab("")
dev.off()
ggsave("subject/clustering_together.pdf",tmp, width = 12, height = 10)


tmp = db.res.full[[100]]

pdf("testes.pdf")
for( sub.id in 1:10 ){
  tmp.all = apply( ListtoArray(lapply( 1:100, function(iii){
    Q = t(tmp$X)%*%tmp$Y.sub%*%data.full$sub.mt + t(tmp$x)%*%diag(tmp$d)%*%data.full$design.mt + 
      matrix( rnorm(length(tmp$Q),sd = sqrt(tmp$er)), ncol = ncol(tmp$Q) )
    CalcDerivativeAnalytic( Q, tmp$sigma, diag(tmp$d)%*%tmp$x, 
                            derivative.mt = rbind( 0, 0, 0, rep(1,ncol(data.full$design.mt)), 0, data.full$design.mt[2:3,] ) )
    })), 1:2, mean )
  
  fff = lapply( db.res.full, function(tmp){
    Q = t(tmp$X)%*%tmp$Y.sub%*%data.full$sub.mt + t(tmp$x)%*%diag(tmp$d)%*%data.full$design.mt + 
      matrix( rnorm(length(tmp$Q),sd = sqrt(tmp$er)), ncol = ncol(tmp$Q) )
    CalcDerivativeAnalytic( Q, tmp$sigma, diag(tmp$d)%*%tmp$x, 
                            derivative.mt = rbind( 0, 0, 0, rep(1,ncol(data.full$design.mt)), 0, data.full$design.mt[2:3,] ) )
  })
  fff.mean = apply( ListtoArray(fff), 1:2, mean )
  
  plot(  data.full$raw.metadata$age_at_collection, tmp.all[5,], col = data.full$raw.metadata$country )
  
  plt.data = data.frame( x = data.full$raw.metadata$age_at_collection,
                         y = fff.mean[5,],
                         country = data.full$raw.metadata$country )
  ggplot( data = plt.data, aes(x=x,y=y,color=country)) + geom_point() + stat_smooth(method = "loess", formula = y ~ x, size = 1)
}
dev.off()

# Derivative plot
db.local.der.est = lapply( db.res.full, function(x){
  apply( ListtoArray(lapply( 1:10, function(iii){
  Q =  t(x$X)%*%x$Y.sub%*%data.full$sub.mt + t(x$x)%*%diag(x$d)%*%data.full$design.mt +
    matrix( rnorm(length(x$Q),sd=sqrt(x$er)), ncol = ncol(x$Q) )
  CalcDerivativeAnalytic( Q, x$sigma, diag(x$d)%*%x$x, rbind( 0, 0, 0, 1, 0, data.full$design.mt[2:3,] ) ) 
   })), 1:2, mean )
} )

db.local.der.est.1 = lapply( db.res.full, function(x){
  # apply( ListtoArray(lapply( 1:10, function(iii){
    # Q =  t(x$X)%*%x$Y.sub%*%data.full$sub.mt + t(x$x)%*%diag(x$d)%*%data.full$design.mt +
      # matrix( rnorm(length(x$Q),sd=sqrt(x$er)), ncol = ncol(x$Q) )
    CalcDerivativeAnalytic( x$Q, x$sigma, diag(x$d)%*%x$x, rbind( 0, 0, 0, 1, 0, data.full$design.mt[2:3,] ) ) 
  # })), 1:2, mean )
} )

db.local.der.est.mean = apply( ListtoArray(db.local.der.est), 1:2, mean )
db.local.der.est.ci = apply( ListtoArray(db.local.der.est), 1:2, quantile, probs = c(0.025,0.975))

db.plt.idx = sample( 1:ncol(data.full$data.mt), 150 )
pdf("subject/db_d_erprior_linear_sub_derivative.pdf", width = 6, height = 5)
for( i in 1:nrow(data.full$data.mt) ){
  plt.data = data.frame( x=data.full$raw.metadata$age_at_collection[db.plt.idx], y = db.local.der.est.mean[i,db.plt.idx],
                        obs.abd = data.full$data.mt.norm[i,db.plt.idx],
                        country = data.full$raw.metadata$country[db.plt.idx],
                        ymin = db.local.der.est.ci[1,i,db.plt.idx], 
                        ymax = db.local.der.est.ci[2,i,db.plt.idx] )
  
  print( ggplot() + geom_point( data=plt.data, aes(x=x,y=y,color=country,size=obs.abd) ) + 
           geom_errorbar( data = plt.data, aes(x=x,ymin=ymin,ymax=ymax), size = 0.25 ) +
           theme_bw() + xlab("age") + ylab("Relative abundance") + 
           ggtitle( sprintf("Genus %s", data.full$all.taxa[6,i])) )
}
dev.off()


# linear design matrix
db.linear.Fin = rbind( 1, 0, 0, data.full$design.mt[4,], data.full$design.mt[5,], 0, 0 )
db.linear.EST = rbind( 1, 1, 0, data.full$design.mt[4,], data.full$design.mt[5,], data.full$design.mt[4,], 0 )
db.linear.RUS = rbind( 1, 0, 1, data.full$design.mt[4,], data.full$design.mt[5,], 0, data.full$design.mt[4,] )

# sero linear
db.linear.sero = rbind( data.full$design.mt[1:4,], 1, data.full$design.mt[6:7,] )
db.linear.nosero = rbind( data.full$design.mt[1:4,], 0, data.full$design.mt[6:7,] )

# linear model
linear.nation = lapply( db.res.full, function(x){
  list( Fin = calc.fix.covariate( x, data.full$design.mt, db.linear.Fin),
        EST = calc.fix.covariate( x, data.full$design.mt, db.linear.EST),
        RUS = calc.fix.covariate( x, data.full$design.mt, db.linear.RUS))
} )

linear.sero = lapply( db.res.full, function(x){
  list( sero = calc.fix.covariate( x, data.full$design.mt, db.linear.sero),
        nosero = calc.fix.covariate( x, data.full$design.mt, db.linear.nosero))
})

pp.nation = plot.fix.covariate( data.full$raw.metadata$age_at_collection, linear.nation, data.full$all.taxa[6,], "subject/full_strat_nation.pdf" )
pp.sero = plot.fix.covariate.new( data.full$raw.metadata$age_at_collection, linear.sero, data.full$all.taxa[6,], "subject/full_strat_sero.pdf" )

# figures for slides
fig.mt = t(v.fix.inter[,1:16])
row.names(fig.mt) = as.factor(1:16)
colnames(fig.mt) = c("Continuous", "Binary", "Interaction")
fig.mt.long = melt(fig.mt)
fig.mt.long$Var1 = as.factor(fig.mt.long$Var1)
fig.mt.long$Var2 = factor(as.character(fig.mt.long$Var2), levels = c("Interaction","Binary","Continuous"))
fig.mt.long$value = as.factor(fig.mt.long$value)

p.v = ggplot(data=fig.mt.long, aes(x=Var1,y=Var2,fill=value)) + geom_tile(color="black") +
  geom_text( aes(x=Var1,y=Var2,label=value) )+ coord_fixed(ratio=1) +
  scale_fill_manual(values=c("blue","skyblue","pink","red"),guide=F) + scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) + xlab("Species") + ylab("") +
  theme(axis.text = element_text(size=12,face="bold"),axis.title = element_text(size=14,face="bold"))
  
ggsave("~/Dropbox (Huttenhower Lab)/Lorenzo_slides/dissertation/figures/v_spec.pdf", width = 8, height = 2)
