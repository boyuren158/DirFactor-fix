source("MCMC.R")
source("utlities.R")
library(phyloseq)
library(ggplot2)
library(reshape2)
library(ggdendro)

# pre-processing
# OTU table
raw.data = read.csv("../data/diabimmune_genus.txt", sep = "\t", row.names = 1 )
# Metadata table
load( "../data/DIABIMMUNE_metadata_full.RData" )

metadata.wgs = metadata[(!is.na(metadata$mgx_reads_filtered))&(!is.na(metadata$gid_wgs)),]
row.names(metadata.wgs) = metadata.wgs$gid_wgs
use.data = raw.data[,names(raw.data)%in%row.names(metadata.wgs)]
#mgx_reads_filtered is the number of reads/1e6
metadata.use = metadata.wgs[names(use.data),]
# convert to count data
use.data.count = t(round(t(use.data/100)*metadata.use$mgx_reads_filtered*1e6))

#get rid of NA in seroconverted
use.data.count = use.data.count[,!is.na(metadata.use$seroconverted)]
metadata.use = metadata.use[!is.na(metadata.use$seroconverted),]

data.full = get.dbdata.ready( use.data.count, metadata.use, 1:ncol(use.data.count),
                              ~country+age.stand+seroconverted+age.stand*country )
# MCMC simulations
hyper = list( nv = 30, a.er = 1, b.er = 0.3, a1 = 10, a2 = 20, m = 10, 
              sub.design = data.full$sub.mt, y.fix = data.full$design.mt,alpha = 10, beta = 0)
mcmc.res = DirFactor.fix( data.full$data.mt, hyper, step = 1e5, thin = 50, step.disp = 1e3 )
# retrieve the MCMC results
db.res.full = lapply( paste(mcmc.res$save.path, seq(20050,1e5,50), sep = "_"), readRDS )

#############################################################
# Derivative plots for all species
db.local.der.est = lapply( db.res.full, function(x){
  apply( ListtoArray(lapply( 1:10, function(iii){
    Q =  t(x$X)%*%x$Y.sub%*%data.full$sub.mt + t(x$x)%*%diag(x$d)%*%data.full$design.mt +
      matrix( rnorm(length(x$Q),sd=sqrt(x$er)), ncol = ncol(x$Q) )
    CalcDerivativeAnalytic( Q, x$sigma, diag(x$d)%*%x$x, rbind( 0, 0, 0, 1, 0, data.full$design.mt[2:3,] ) ) 
  })), 1:2, mean )
} )

db.local.der.est.mean = apply( ListtoArray(db.local.der.est), 1:2, mean )
db.local.der.est.ci = apply( ListtoArray(db.local.der.est), 1:2, quantile, probs = c(0.025,0.975))

db.plt.idx = sample( 1:ncol(data.full$data.mt), 150 )
pdf("db_derivative.pdf", width = 6, height = 5)
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

#############################################################
# population trends for all species
plot.pop.curve( data.full$design.mt[4,], db.res.full, data.full, "db_trend_age_by_nation.pdf", width = 6, height = 5 )

#############################################################
# species relations
# phylogenetic tree of genera in the dataset
genus.names = data.full$all.taxa[6,]
genus.tree = read_tree("../data/genus_tree.tre")
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

# Unifrac distances between genera
phylo.dist = UniFrac( phyloseq.final, weighted = F )
phylo.dendro = hclust( phylo.dist )
phylo.dendro$labels = sapply( phylo.dendro$labels, function(x) strsplit(x,"_")[[1]][2] )
phylo.ord = cmdscale(phylo.dist,k=43,eig = T)
# Gram matrix between genera based on Unifrac distances
phylo.sim = cov2cor((phylo.ord$points%*%t(phylo.ord$points)))
colnames(phylo.sim) = sapply( strsplit(colnames(phylo.sim), split = "_", fixed=T), function(x) x[2] )
row.names(phylo.sim) = colnames(phylo.sim)

dendro.labels.ordered = phylo.dendro$labels[phylo.dendro$order]
phylo.sim = phylo.sim[dendro.labels.ordered,dendro.labels.ordered]

spe.corr.all = lapply( db.res.full, function(x) cov2cor(t(x$X)%*%(x$X)) )
spe.corr.mean = apply( ListtoArray(spe.corr.all), 1:2, mean )
spe.corr.sub = spe.corr.all[sample(1:length(spe.corr.all),800,replace = T)]
spe.dist.sub = lapply( spe.corr.all, function(x) 1-x )
spe.corr.raw = cor( t( data.full$data.mt.norm ) )
spe.corr.fix = apply( ListtoArray(lapply(db.res.full, function(x) cov2cor(t(x$x)%*%diag(x$d)%*%diag(x$d)%*%x$x))), 1:2, mean )

row.names(spe.corr.mean) = genus.names
colnames(spe.corr.mean) = genus.names
plot.genus.corr.mean = spe.corr.mean[dendro.labels.ordered,dendro.labels.ordered]

row.names(spe.corr.raw) = genus.names
colnames(spe.corr.raw) = genus.names
plot.genus.corr.raw = spe.corr.raw[dendro.labels.ordered,dendro.labels.ordered]

row.names(spe.corr.fix) = genus.names
colnames(spe.corr.fix) = genus.names
plot.genus.corr.fix = spe.corr.fix[dendro.labels.ordered,dendro.labels.ordered]

plot.genus.corr.raw[lower.tri(plot.genus.corr.mean)] = plot.genus.corr.mean[lower.tri(plot.genus.corr.mean)]
plot.genus.corr.fix[lower.tri(plot.genus.corr.mean)] = plot.genus.corr.mean[lower.tri(plot.genus.corr.mean)]


dendro.labels.num = sapply( dendro.labels.ordered, function(x) which(genus.names==x) )
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
ggsave("species_ordination.pdf", p1[[1]], width = 10, height = 8 )

p1.1 = gg.heatmap( t(plot.genus.corr.raw), y.ticks = dendro.labels.ordered, x.ticks = dendro.labels.ordered, y.lab = "Correlation of random effect", 
                   x.lab = "Correlation of raw data" )
p2.1 = gg.heatmap( t(plot.genus.corr.fix), y.ticks = dendro.labels.ordered, x.ticks = dendro.labels.ordered, y.lab = "Correlation of random effect", 
                   x.lab = "Correlation of fixed effect" )
pdf("species_similarities.pdf", width = 16, height = 7 )
multiplot(p1.1,p2.1,cols =2)
dev.off()

p2 <- ggplot(segment(ddata_xy)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  theme_none + theme( plot.margin = unit(c(0,0,-1,0), 'cm') )

pdf("phylogenetic_tree.pdf", width = 8, height = 7)
grid.newpage()
print(p2, vp=viewport(0.515, 0.2, x=0.45, y=0.87))
dev.off()

rv.coef( plot.genus.corr.mean, phylo.sim )
rv.coef( plot.genus.corr.raw, phylo.sim )
rv.coef( plot.genus.corr.fix, phylo.sim )