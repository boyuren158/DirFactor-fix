library(grid)
library(ggplot2)
library(reshape2)

mt.trace = function( A ){
  sum(diag(A))
}

rv.coef = function( A, B ){
  mt.trace(A%*%B)/sqrt(mt.trace(A%*%A)*mt.trace(B%*%B))
}

ListtoArray = function( ls.input ){
  # convert list of matrices to an array
  array( unlist( ls.input ), dim = c(dim(ls.input[[1]]),length(ls.input)) )
}

multiplot <- function(..., plotlist=NULL, file, cols=1, byrow = F, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),byrow = byrow,
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

gg.heatmap <- function( cormt, x.ticks = NULL, y.ticks = NULL, x.lab="", y.lab="", title="", border = "black" ){
  plot.dat = reshape2::melt( cormt )
  plot.dat$Var1 = as.factor( plot.dat$Var1 )
  plot.dat$Var2 = as.factor( plot.dat$Var2 )
  p = ggplot( data = plot.dat, aes( x=Var1, y=Var2, fill=value ) ) + geom_tile(color = border) + 
    scale_fill_gradient2( low="steelblue", mid = "white", high = "red", limits=c(-1, 1) ) + 
    theme_bw() + coord_fixed() + 
    scale_y_discrete(labels = y.ticks, expand = c(0,0)) + 
    scale_x_discrete(labels = x.ticks, expand = c(0,0)) +
    theme( axis.title = element_text(face="bold"), 
           plot.title = element_text(face="bold"),
           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=6),
           axis.text.y = element_text(size=6)) + 
    xlab(x.lab) + ylab(y.lab) + ggtitle( title )
  if( all(is.null(x.ticks)) ){
    p = p + theme( axis.ticks.x = element_blank() )
  }
  if( all(is.null(y.ticks)) ){
    p = p + theme( axis.ticks.y = element_blank() )
  }
  p
}

summarize.mcmc.fixed = function( filepath, index, y.fix, P.proj ){
  res.raw = lapply( paste( filepath, index, sep = "_" ), readRDS )
  
  #Sigma given by YY'
  Sigma.mean = tryCatch({
    Sigma.raw = lapply( res.raw, function(x) P.proj%*%t(x$Y)%*%x$Y%*%P.proj + diag( n ) )
    apply( array( unlist(Sigma.raw), dim = c(dim(Sigma.raw[[1]]), length(Sigma.raw)) ), 1:2, mean )},
    error = function(err){
      return(NA)
    })
  
  #global scale factor d, matrix, each row for a component in d
  d.all = tryCatch({
    sapply( res.raw, function(x) x$d )},
    error = function(err){
      return(NA)
    })
  
  #individual effects v, 3d dimension, 3rd dim is posterior samples
  v.all = tryCatch({
    array( unlist( lapply( res.raw, function(x) x$x ) ), dim = c( dim( res.raw[[1]]$x ), length( res.raw ) ) )},
    error = function(err){
      return(NA)
    })
  
  #individual effects beta, 3d dimension, 3rd dim is posterior samples
  beta.all = tryCatch({
    beta.raw = lapply( res.raw, function(x) diag(x$d, nrow=length(x$d), ncol=length(x$d))%*%x$x )
    array( unlist( beta.raw ), dim = c( dim(beta.raw[[1]]), length( beta.raw ) ) )},
    error = function(err){
      return(NA)
    })
  
  return( list( Sigma = Sigma.mean, d = d.all, beta = beta.all, v = v.all ) )
}

boxplot.coef = function( v.all, beta.all, truth, legend, glm.res = NA ){
  plt.data = data.frame( est = c(c( v.all ), c( beta.all )), coef = as.factor(rep( 1:(prod(dim(v.all)[1:2])), dim(v.all)[3] )), 
                         group = rep( legend, each = length(v.all) ) )
  plt.truth = data.frame( x = as.factor(1:(prod(dim(v.all)[1:2]))), y = truth)
  
  p1 = ggplot() + geom_boxplot(data = plt.data, aes( y=est, x = coef, fill = group ),outlier.size = 0) +
    scale_y_continuous(limits = c(-2.5,2.5)) + geom_point( data = plt.truth, aes( x = x, y = y ), size = 3 )
  if( !is.na(glm.res) ){
    plt.glm = data.frame( x = as.factor(1:(prod(dim(v.all)[1:2]))), y = c(glm.res) )
    p1 = p1 + geom_point( data = plt.glm, aes( x = x, y = y ), size = 3, shape = 2 )
  }
  p1
}

boxplot.coef.simple = function( v.all, truth ){
	plt.data = data.frame( est = c(v.all), coef = as.factor( rep( 1:(prod(dim(v.all)[1:2])), dim(v.all)[3]) ) )
	plt.truth = data.frame( x = as.factor(1:(prod(dim(v.all)[1:2]))), y = truth)
	ggplot() + geom_boxplot(data = plt.data, aes( y=est, x = coef ),outlier.size = 0) +
    scale_y_continuous(limits = c(-2.5,2.5)) + geom_point( data = plt.truth, aes( x = x, y = y ), size = 3 )
}

# visualize confidence band
plot.conf.band = function( x, pred.mean, pred.ci, tru.mean ){
  require(ggplot2)
  plt.data = data.frame( x = x, y = pred.mean, 
                         y.true = tru.mean,
                         ymin = pred.ci[1,],
                         ymax = pred.ci[2,])
  ggplot(data = plt.data) + geom_line( aes(x=x,y=y,color="Predicted") ) +
    geom_line( aes(x=x,y=y.true,color="Truth") ) +
    geom_ribbon( aes(x=x,ymin=ymin,ymax=ymax,fill = "Credible band"), alpha = 0.2 ) +
    ylab("Relative abundance") + theme_bw()
}

# localized derivative
# calculate the truth
CalcDerivativeTruth = function( sim.truth, yinc, i.var ){
  # for sim.truth, you need the true Q, true sigma, true v, true abundance and true y.fix
  # i.var is the index of the covariate you want to consider
  # this fnc can only deal linear latent trend
  y.new = sim.truth$y
  y.new[i.var,] = y.new[i.var,] + yinc
  Q.new = sim.truth$Q + t(sim.truth$v)%*%(y.new-sim.truth$y)
  
  w.new = Q.new^2*(Q.new>0)*sim.truth$sigma
  abd.new = t(t(w.new)/colSums(w.new))
  
  (abd.new-sim.truth$abd)/yinc
}

CalcDerivativeSpline = function( mcmc.ind, spline.basis, y.obs, y.inc ){
  y.new = y.obs + y.inc*(1-2*(y.obs==max(y.obs)))
  spline.new = predict( spline.basis, y.new )
  
  Q.new = mcmc.ind$Q + t(mcmc.ind$x)%*%diag(mcmc.ind$d,nrow=length(mcmc.ind$d))%*%(rbind(1,t(spline.new)) - rbind(1,t(spline.basis)))
  w.new = Q.new^2*(Q.new>0)*mcmc.ind$sigma
  w.basis = mcmc.ind$Q^2*(mcmc.ind$Q>0)*mcmc.ind$sigma
  
  t((t(w.new)/colSums(w.new) - t(w.basis)/colSums(w.basis))/(y.inc*(1-2*(y.obs==max(y.obs)))))
}

# don't forget the scaling factor d for the derivative.mt
CalcDerivativeAnalytic = function( Q, sigma, v, derivative.mt ){
  Q.pos = Q*(Q>0)
  final.weights = sigma*Q.pos^2
  t((2*t(Q.pos*(t(v)%*%derivative.mt)*sigma)*colSums(final.weights) - 
       2*t(final.weights)*colSums(sigma*Q.pos*(t(v)%*%derivative.mt)))/
      colSums(final.weights)^2)
}

CalcDerivativeAnalytic.mis = function( Q, sigma, v, derivative.mt ){
  Q.pos = Q*(Q>0)
  final.weights = sigma*Q.pos
  t((t((t(v)%*%derivative.mt)*sigma*(Q>0))*colSums(final.weights) - 
       t(final.weights)*colSums(sigma*(Q>0)*(t(v)%*%derivative.mt)))/
      colSums(final.weights)^2)
}

PlotDerivative = function( plt.idx, all.covariate, obs.abd, tru.derivative, mean.derivative, ci.derivative, 
                           time, var.exp, plot.name, metadata.raw, min.0 = 1e-6 ){
  if( is.na( tru.derivative ) ){
    plt.data = data.frame( x = time,
                           y = mean.derivative[plt.idx,],
                           ymin = ci.derivative[1,plt.idx,],
                           ymax = ci.derivative[2,plt.idx,],
                           obs.abd = obs.abd[plt.idx,],
                           country = metadata.raw$country )
    p = ggplot() + geom_point( data = plt.data, aes( x=x, y=y, color = country, size = obs.abd ) ) + 
      geom_errorbar(data = plt.data, aes(x=x, ymin=ymin, ymax=ymax ) )
    p = p + geom_hline( yintercept = 0, linetype = 2, color = "black", size = 1.5)
  }else{
    plt.data = data.frame( x = rep(time,2),
                           y = c(mean.derivative[plt.idx,], tru.derivative[plt.idx,]),
                           obs.abd = rep(obs.abd[plt.idx,], 2),
                           group = rep(as.factor(all.covariate[2,]),2),
                           type = rep( c("Est.","Obs."), each = length(time) ) )
    plt.ribbon = data.frame( x = time, ymin = ci.derivative[1,plt.idx,],
                             ymax = ci.derivative[2,plt.idx,] )
    p = ggplot() + geom_point( data = plt.data, aes(x=x,y=y,color=group,size=obs.abd,shape=type) ) +
      geom_errorbar( data = plt.ribbon, aes(x=x,ymin=ymin,ymax=ymax) ) + 
      geom_hline( yintercept = 0, linetype =2 , color = "black", size = 1.5) + scale_shape_manual( values = c(16,4) )
  }
  # cov.prop = mean( (tru.derivative[plt.idx,]*(abs(tru.derivative[plt.idx,])>min.0)>=plt.data$ymin)&(tru.derivative[plt.idx,]*(abs(tru.derivative[plt.idx,])>min.0)<=plt.data$ymax) )
  # if(!smooth){
  
  # }else{
  #   g1 = ggplot(data = plt.data) + stat_smooth(aes(x=x,y=ymin), size=0, method = "loess", se = FALSE) +
  #     stat_smooth(aes(x=x,y=ymax), size=0, alpha = 0.2, method = "loess", se = FALSE)
  #   gg1 = ggplot_build(g1)
  #   plt.data2 = data.frame(x = gg1$data[[1]]$x,
  #                          ymin = gg1$data[[1]]$y,
  #                          ymax = gg1$data[[2]]$y) 
  #   p = g1 + geom_ribbon(data = plt.data2, aes(x = x, ymin = ymin, ymax = ymax, fill = "Credible band"), alpha = 0.2) +
  #     stat_smooth(aes(x=x,y=y,color="Estimated"), method = "loess", se = FALSE)
  # }
  # adding the baseline
  
  p1 = p + ylab("Derivative of probability of species") + xlab("covariate") + 
    ggtitle( sprintf("Species %s, diff.var=%f, overall diff.var=%f", 
                     plot.name[plt.idx], var.exp$ind[plt.idx], var.exp$overall) ) + theme_bw()
  list(p1)
  # if(!is.na(tru.derivative)){
  #   list(p1,p2)
  # }
  # else{
  #   list(p1)
  # }
  # }
  # else{
  #   p + stat_smooth(data=plt.tru, aes(x=x, y=y, color = "Truth"), method = "loess", se = FALSE) + ylab("Derivative of probability of species") +
  #     xlab("covariate") + ggtitle( sprintf("Species %d", plt.idx) ) + theme_bw()
  # }
}

DerivativeSimSummary = function( obs.abd, mcmc.res, covariate, covariate.raw, 
                                 covariate.var, covariate.const, tru.derivative, 
                                 fprime.mt, file.name, sample.plt,
                                 plot.name, plot.ctrl = rep(T,7), metadata.raw = NA, 
                                 prior = F){
  const.abd = apply( ListtoArray( 
    lapply( mcmc.res, function(x) calc.median.p.post( covariate.var, covariate.const, x, prior ) ) ), 
    1:2, mean )
  var.stat = change.variability( obs.abd, mcmc.res, covariate.var, covariate.const, prior )
  var.exp = list( overall = var.stat$diff.all/var.stat$raw.all,
                  ind = var.stat$diff.ind/var.stat$raw.ind )
  
  latent.var.exp = sapply( mcmc.res, function(x){
    Q.var = apply(x$Q,1,var)
    if(any(is.null(x[["d"]]))){
      Q.nofixed = apply( x$Q - t(x$x)%*%covariate.var, 1, var )
    }else{
      Q.nofixed = apply( x$Q - t(x$x)%*%diag(x$d)%*%covariate.var, 1, var )
    }
    1-Q.nofixed/Q.var
  })
  
  if(plot.ctrl[1]){
    pdf(paste(file.name,"_constvsRaw.pdf",sep=""), width = 10, height = 5)
    par(mfrow=c(1,2))
    for( idx in 1:nrow(obs.abd) ){
      hist( obs.abd[idx,],freq = F, main = sprintf("Species %s", plot.name[idx]) )
      hist( const.abd[idx,], freq = F, add = T, col = "gray" )
      boxplot( list(raw = obs.abd[idx,], const = const.abd[idx,]) )
    }
    par(mfrow=c(1,1))
    dev.off()
  }
  
  #check if the latent variables are correlated with the covariate
  #first we consider the overall correlation
  if(plot.ctrl[2]){
    corr.overall.covariate.mt = sapply( mcmc.res, function(x){
      Q.latent = t(x$X)%*%x$Y
      apply( Q.latent, 1, function(x) cor(x,covariate) )
    })
    pdf(paste(file.name,"_allvsCovariate.pdf",sep=""), width = 5, height = 5)
    for( idx in seq_len(nrow(corr.overall.covariate.mt)) ){
      hist( corr.overall.covariate.mt[idx,], freq = F,
            main = sprintf("Correlation between covariate and latent effects for species %s", 
                           plot.name[idx]), xlab = "Correlation" )
    }
    dev.off()
  }
  if(plot.ctrl[3]){
    #then check individual latent factors
    pdf(paste(file.name,"_YvsCovariate.pdf",sep=""), width = 5, height = 5)
    corr.Y.covariate.mt = sapply( mcmc.res, function(x) 
      apply( x$Y, 1, function(x.row) cor(x.row,covariate)) )
    for( Y.idx in seq_len( nrow(corr.Y.covariate.mt) ) ){
      hist( corr.Y.covariate.mt[Y.idx,], freq = F, 
            main = sprintf("Correlation between covariate and %dth sample factor", Y.idx), xlab = "Correlation" )
    }
    dev.off()
  }
  
  #species correlation
  if(plot.ctrl[4]){
    corr.species.mean = apply( ListtoArray( lapply( mcmc.res, function(x) cor(x$X) ) ), 1:2, mean )
    corr.species.fixed.mean = apply( ListtoArray( lapply( mcmc.res, function(x){
      cor(rbind(diag(x$d)%*%x$x,x$X))
    }) ), 1:2, mean )
    corr.species.raw = cor(t(obs.abd))
    colnames(corr.species.raw) = NULL
    row.names(corr.species.raw) = NULL
    pdf(paste(file.name,"_speciesCor.pdf",sep=""), width = 7, height = 5)
    p1 = gg.heatmap( corr.species.raw, title = "Raw correlation" )
    p2 = gg.heatmap( corr.species.mean, title = "Correlation of factors" )
    p3 = gg.heatmap( corr.species.fixed.mean, title = "Correlation of factors+reg coefs")
    multiplot( p1, p2, p3, cols = 3 )
    dev.off()
  }
  
  #derivative
  if( plot.ctrl[5] ){
    local.der.all = ListtoArray( lapply( mcmc.res, function(x){
      if(!is.null(x[["d"]])){
        fprime.mt = diag(x$d,nrow=length(x$d))%*%fprime.mt
      }
      CalcDerivativeAnalytic( x$Q, x$sigma, x$x, fprime.mt )
    }) )
    local.der.mean = apply( local.der.all, 1:2, mean )
    local.der.ci = apply( local.der.all, 1:2, quantile, probs = c(0.025,0.975) )
    pdf(paste(file.name,"_derivative.pdf",sep=""), width = 6, height = 5)
    for( i in 1:nrow(local.der.mean) ){
      if(!is.na(tru.derivative)){
        p.deriv = PlotDerivative( plt.idx = i, all.covariate = covariate.var, obs.abd = obs.abd[,sample.plt], tru.derivative = tru.derivative[,sample.plt], 
                              mean.derivative = local.der.mean[,sample.plt], 
                              ci.derivative = local.der.ci[,,sample.plt], 
                              time = covariate.raw[sample.plt], var.exp = var.exp, plot.name = plot.name )
        # plt.raw.data = data.frame( time = covariate.raw[sample.plt], relative.abd = obs.abd[i,sample.plt] )
        # p.raw = ggplot(data = plt.raw.data,aes(x=time, y = relative.abd )) + geom_point(size=1.5) + 
        #   ggtitle(sprintf("Species %s, %f variance explained by fix", plot.name[i], round(latent.var.exp[i],3))) + theme_bw()
      }
      else{
        p.deriv = PlotDerivative( plt.idx = i, all.covariate = covariate.var, obs.abd = obs.abd[,sample.plt], tru.derivative = NA, 
                                  mean.derivative = local.der.mean[,sample.plt], 
                                  ci.derivative = local.der.ci[,,sample.plt], 
                                  time = covariate.raw[sample.plt], var.exp = var.exp, plot.name = plot.name,
                                  metadata.raw = data.full$raw.metadata )
        country = rep( "Finland", ncol(covariate.var) )
        country[as.logical(covariate.var[2,])] = "Estonia"
        country[as.logical(covariate.var[3,])] = "Russia"
        plt.raw.data = data.frame( time = covariate.raw[sample.plt], relative.abd = obs.abd[i,sample.plt],
                                   country = country[sample.plt], sero = as.factor(covariate.var[5,sample.plt]) )
        # p.raw = ggplot(data = plt.raw.data,aes(x=time, y = relative.abd, color = country, shape = sero )) + geom_point(size=1.5) + 
        #   ggtitle(sprintf("Species %s, %f variance explained by fix", plot.name[i], round(latent.var.exp[i],3))) + theme_bw()
      }
      # p.deriv[[length(p.deriv)+1]] = p.raw
      multiplot(plotlist = p.deriv, cols = length(p.deriv))
    }
    dev.off()
  }
  
  #overall changes proportion
  if( plot.ctrl[6] ){
    pdf(paste(file.name,"_constPropall.pdf", sep = ""), width = 7, height = 5)
    boxplot( var.stat$diff.all.post/var.stat$raw.all, ylab = "Overall proportion of changes" )
    dev.off()
  }
  if( plot.ctrl[7] ){
    pdf(paste(file.name,"_constPropind.pdf", sep = ""), width = 7, height = 5)
    for( i in 1:nrow(obs.abd) ){
      boxplot( var.stat$diff.ind.post[i,]/var.stat$raw.ind[i], ylab = "Individual proportion of changes", main = plot.name[i] )
    }
    dev.off()
  }
}

#checked the variance explained
calc.var.exp = function( y.fix, v, Q ){
  all.var = apply( Q, 1, var )
  fix.var = apply( t(v)%*%y.fix, 1, var )
  fix.var/all.var
}

calc.median.p = function( y.fix, v, X, Y, sigma, er ){
  error = matrix( rnorm( ncol(X)*ncol(Y), sd = sqrt(er) ), nrow = ncol(X) )
  
  fix.median = t(apply( y.fix, 1, median )%*%v)[,rep(1,ncol(y.fix))]
  Q.median = fix.median + t(X)%*%Y + error
  weight.median = sigma*Q.median^2*(Q.median>0)
  t(t(weight.median)/colSums(weight.median))
}

calc.median.p.post = function( y.fix, y.median, sim, prior = F ){
  library(truncnorm)
  if(prior){
    sim$x = matrix( rnorm(length(sim$x)), nrow = nrow(sim$x) )
    if(!is.null(sim[["d"]])){
      sim$d = rtruncnorm( length(sim$d), a = 0 )
    }
  }
  if(is.null(sim[["d"]])){
    Q.nofix = sim$Q - t(sim$x)%*%y.fix
    Q.median = Q.nofix + t(y.median%*%sim$x)[,rep(1,ncol(y.fix))]
  }else{
    Q.nofix = sim$Q - t(sim$x)%*%diag(sim$d)%*%y.fix
    Q.median = Q.nofix + t(y.median%*%diag(sim$d)%*%sim$x)[,rep(1,ncol(y.fix))]
  }
  weight.median = sim$sigma*Q.median^2*(Q.median>0)
  t(t(weight.median)/colSums(weight.median))
}

distance.disperse = function( mt ){
  mt.mean = rowMeans(mt)
  apply( mt, 2, function(x) vegdist(rbind(x,mt.mean)) )
}

# barplot visualization
plot.bar = function( abd, covariate, n.plot ){
  require( reshape2 )
  #sort by covariate
  #sort by avg species abd
  
  plt.data = melt( abd[1:n.plot,] )
  #plt.data$Var1 = as.factor( plt.data$Var1 )
  plt.data$Var1 = as.factor( paste("Species",plt.data$Var1,sep="") )
  ggplot( data = plt.data, aes( x=Var2, y=value, fill = Var1 ) ) + geom_bar(stat = "identity") + ylim(c(0,1)) + theme_bw() +
    xlab("Biological sample") + ylab("Relative abundance")
}

# Plot distance dispersion
dist.disperse.const = function( mcmc.res, y.fix, y.const ){
  constant.abd = apply( ListtoArray( 
    lapply( mcmc.res, function(x) calc.median.p.post( y.fix, y.const,x ) ) ), 
    1:2, mean )
  distance.disperse( constant.abd )
}

# Plot the variability of changes
change.variability = function( raw.data, mcmc.res, y.fix, y.const, prior = F ){
  raw.diff = raw.data - rowMeans(raw.data)
  raw.var.all = sd( apply( raw.diff, 2, function(x) sqrt(sum(x^2))) )
  raw.var.ind = apply( raw.data, 1, sd )
  diff.all = ListtoArray( lapply( mcmc.res, function(x) calc.median.p.post( y.fix, y.const, x, prior ) - raw.data ) )
  diff.mt = apply( diff.all, 1:2, mean )
  diff.var.all = mean( apply( diff.mt, 2, function(x) sqrt(sum(x^2)) ) )
  diff.var.ind = apply( diff.mt, 1, function(x) mean(abs(x)) )
  diff.var.all.post = apply( diff.all, 3, function(x) mean( apply( x, 2, function(xx) sqrt(sum(xx^2)) ) ) )
  diff.var.ind.post = apply( diff.all, 3, function(x) apply( x, 1, function(xx) mean(abs(xx)) ) )
  
  list( raw.all = raw.var.all, raw.ind = raw.var.ind,
        diff.all = diff.var.all, diff.ind = diff.var.ind,
        diff.all.post = diff.var.all.post, diff.ind.post = diff.var.ind.post )
}

calc.fix.covariate = function( res, y.raw, y.fix ){
  if(!any(is.null(res[["d"]]))){
    Q.new = res$Q + t(res$x)%*%diag(res$d)%*%(y.fix-y.raw)
  }
  else{
    Q.new = res$Q + t(res$x)%*%(y.fix-y.raw)
  }
  weight.new = res$sigma*Q.new^2*(Q.new>0)
  t(t(weight.new)/colSums(weight.new))
}

plot.fix.covariate = function( covariate, res.all, plot.name, file.name, truth = NA ){
  res.Fin = ListtoArray( lapply( res.all, function(x) x$Fin ) )
  res.EST = ListtoArray( lapply( res.all, function(x) x$EST ) )
  res.RUS = ListtoArray( lapply( res.all, function(x) x$RUS ) )
  
  covariate.cut = cut( covariate, breaks = c(0, quantile( covariate, probs = seq(0.2,0.8,0.2) ), max(covariate)) )
  res.EST.Fin.time = ListtoArray( apply( res.EST - res.Fin, 3, function(x){
    aggregate(x=t(x),by=list(covariate.cut),FUN=mean)[,-1]
  } ) )
  res.RUS.Fin.time = ListtoArray( apply( res.RUS - res.Fin, 3, function(x){
    aggregate(x=t(x),by=list(covariate.cut),FUN=mean)[,-1]
  } ) )
  
  # try to get the posterior probability
  res.EST.Fin.pp = apply( res.EST.Fin.time, 1:2, function(x) min(mean(x>0),mean(x<=0)) )
  res.RUS.Fin.pp = apply( res.RUS.Fin.time, 1:2, function(x) min(mean(x>0),mean(x<=0)) )
  
  pdf( file.name, width = 6, height = 5 )
  for( i in 1:(dim(res.EST.Fin.time)[2]) ){
    plot.data.time = data.frame( compare = rep(c("EST vs. FIN","RUS vs. FIN"), each = dim(res.EST.Fin.time)[1]*dim(res.EST.Fin.time)[3]),
                                 time = factor( rep(levels(covariate.cut), 2*dim(res.EST.Fin.time)[3]), levels = levels(covariate.cut)),
                                 change = c(c(res.EST.Fin.time[,i,]),
                                            c(res.RUS.Fin.time[,i,])))
    plt.annote = data.frame( compare = rep(c("EST vs. FIN", "RUS vs. FIN"), each = dim(res.EST.Fin.time)[1] ),
                             time = factor( rep(levels(covariate.cut), 2 ) ),
                             posterior.p = as.character(round(c(c(res.EST.Fin.pp[,i]),c(res.RUS.Fin.pp[,i])),digits = 3)),
                             y = c(apply(res.EST.Fin.time,1:2,max)[,i], apply(res.RUS.Fin.time,1:2,max)[,i]) )
    p2 = ggplot() +
      geom_boxplot(data = plot.data.time, aes(x=time, y=change, fill=compare), outlier.shape=NA) + ylab("Change in relative abundance") +
      ggtitle(plot.name[i]) + 
      geom_text(data=plt.annote, aes(x=time, y=y, label = posterior.p, group=compare), position = position_dodge(width = 0.75) ) + theme_bw()
    
    if(!any(is.na(truth))){
      plt.truth = data.frame(time=1:nrow(truth), y = truth[,i])
      p2 = p2 + geom_point( data = plt.truth, aes(x=time,y=y), shape = 17 )
    }
    
    # multiplot( p1, p2, cols = 2 )
    print(p2)
  }
  dev.off()
}

plot.fix.covariate.new = function( covariate, res.all, plot.name, file.name ){
  res.sero = ListtoArray( lapply( res.all, function(x) x$sero ) )
  res.nosero = ListtoArray( lapply( res.all, function(x) x$nosero ) )
  
  covariate.cut = cut( covariate, breaks = c(-Inf, quantile( covariate, probs = seq(0.2,0.8,0.2) ), Inf) )
  res.sero.time = ListtoArray( apply( res.sero - res.nosero, 3, function(x){
    aggregate(x=t(x),by=list(covariate.cut),FUN=mean)[,-1]
  } ) )
  
  # try to get the posterior probability
  res.sero.pp = apply( res.sero.time, 1:2, function(x) min(mean(x>0),mean(x<=0)) )
  
  pdf( file.name, width = 6, height = 5 )
  for( i in 1:dim(res.sero.time)[2] ){
    plot.data.time = data.frame( compare = rep(c("Sero vs. Nonsero"), each = dim(res.sero.time)[1]*dim(res.sero.time)[3]),
                                 time = factor( rep(levels(covariate.cut), 2*dim(res.sero.time)[3]), levels = levels(covariate.cut)),
                                 change = c(res.sero.time[,i,]) )
    plt.annote = data.frame( time = factor( levels(covariate.cut) ),
                             posterior.p = as.character(round(c(res.sero.pp[,i]),digits = 3)),
                             y = apply(res.sero.time,1:2,max)[,i] )
    p2 = ggplot() +
      geom_boxplot(data = plot.data.time, aes(x=time, y=change, fill=compare), outlier.shape=NA) + ylab("Change in relative abundance") +
      ggtitle(plot.name[i]) + 
      geom_text(data=plt.annote, aes(x=time, y=y, label = posterior.p ) ) + theme_bw()
    
    # multiplot( p1, p2, cols = 2 )
    print(p2)
  }
  dev.off()
}

plot.deriv.simple = function( w, est, est.ci, truth, obs.abd, w2, plt.idx ){
  plt.data = data.frame( w = w[plt.idx], truth = truth[plt.idx], est = est[plt.idx], 
                         est.l = est.ci[1,plt.idx], 
                         est.u = est.ci[2,plt.idx], 
                         obs.abd = obs.abd[plt.idx], 
                         group = w2[plt.idx] )
  ggplot() + geom_point( data = plt.data, aes(x=w, y = est, color = group, size = obs.abd) ) + 
    geom_point( data = plt.data, aes(x=w, y = truth, color = group, size = obs.abd), shape = 4 ) +
    geom_errorbar( data = plt.data, aes(x=w, ymin = est.l, ymax=est.u) )
}

get.pop.mean = function( y.fix, sim.ind ){
  Q = t(sim.ind$X)%*%sim.ind$Y.sub + t(sim.ind$x)%*%diag(sim.ind$d)%*%y.fix  + 
    matrix( rnorm( ncol(sim.ind$X)*ncol(sim.ind$Y.sub), sd = sqrt(sim.ind$er) ), nrow = ncol(sim.ind$X) )
  weight = sim.ind$sigma*Q^2*(Q>0)
  apply( t(weight)/colSums(weight), 2, mean )
}

get.pop.mean.mis = function( y.fix, sim.ind ){
  Q = t(sim.ind$X)%*%sim.ind$Y.sub + t(sim.ind$x)%*%diag(sim.ind$d)%*%y.fix  + 
    matrix( rnorm( ncol(sim.ind$X)*ncol(sim.ind$Y.sub), sd = sqrt(sim.ind$er) ), nrow = ncol(sim.ind$X) )
  weight = sim.ind$sigma*Q*(Q>0)
  apply( t(weight)/colSums(weight), 2, mean )
}

plot.trend.simple = function( w, est0, est1, tru0, tru1 ){
  plt.data = data.frame( w = rep(w, 2), est.lower = c(est0[1,],est1[1,]), 
                         est.upper = c(est0[2,], est1[2,]),
                         truth = c(tru0,tru1),
                         group = as.factor(rep(c(0,1), each = length(w))) )
  ggplot(data = plt.data) + geom_ribbon(aes(x = w, ymin = est.lower, ymax = est.upper, fill = group ), alpha = 0.2) +
    geom_line( aes(x=w, y = truth, color = group ) ) + theme_bw()
}

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

gg.color.hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

PlotStatis = function( all.cov, n.dim = 2, levels.min = FALSE, labels = NA, types = NA, dist = F, ... ){
  all.cov.array = array( unlist( all.cov ), dim = c(dim(all.cov[[1]]), length(all.cov)) )
  all.statis.res = DistatisR::distatis( all.cov.array, Distance = dist, nfact2keep = n.dim )
  compromise.ev = eigen(all.statis.res$res4Splus$Splus)$values
  compromise.prop = compromise.ev[1:n.dim]/sum(compromise.ev)
  compromise.coord = apply( all.statis.res$res4Splus$PartialF, 2, rbind )
  n.rep = nrow( compromise.coord )/nrow(all.cov[[1]])
  
  apply( combn(1:n.dim, 2), 2, function(axis.idxs){
    axis.i = axis.idxs[1]
    axis.j = axis.idxs[2]
    plot.data = data.frame( x = compromise.coord[,axis.i], y = compromise.coord[,axis.j], 
                            group = rep( 1:nrow(all.cov[[1]]), n.rep) )
    
    if(!any(is.na(types))){
      types = as.factor(types)
      types.color = gg.color.hue( length(levels(types) ) )
      plot.data$types = rep(types, n.rep)
      contr = ggplot() + geom_density2d(data = plot.data, aes( x=x, y=y, group = group, color = types)) + 
        scale_color_manual( values = types.color ) + theme_bw()
    }
    else{
      contr = ggplot() + geom_density2d (data = plot.data, aes( x=x, y=y, group = group) ) +
        theme_bw()
    }
    if(levels.min){
      contr.data = ggplot_build(contr)$data[[1]]
      contr.min = contr.data[grep("-001", contr.data$group, fixed = T),]
      contr = ggplot(NULL) + geom_path( data = contr.min, aes(x=x,y=y,group=group,color=colour), ... ) + 
        scale_color_manual( values = sort(unique(contr.min$colour)), labels = unique(types)[order(unique(contr.min$colour))] ) + 
        theme_bw()
    }
    
    if(!any(is.na(labels))){
      x.annote = tapply( plot.data$x, plot.data$group, mean )
      y.annote = tapply( plot.data$y, plot.data$group, mean )
      annote.data = data.frame( x=x.annote, y=y.annote, 
                                labels = as.character( labels ) )
      contr = contr + geom_point( data = annote.data, aes( x=x, y=y )  ) +
        with(annote.data, annotate(geom="text", x = x+0.01 , y = y, label = labels, size = 6) )
    }
    contr + xlab(sprintf("Compromise axis %d (%.2f%%)", axis.i, compromise.prop[axis.i]*100)) +
      ylab(sprintf("Compromise axis %d (%.2f%%)", axis.j, compromise.prop[axis.j]*100))
  } )
}
