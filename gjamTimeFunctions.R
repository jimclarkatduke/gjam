
colF <- colorRampPalette( c('#8c510a','#d8b365','#c7eae5','#5ab4ac','#01665e','#2166ac') )

mergeList <- function( list1, list2 ){
  
  # update elements of list1 if contained list2
  
  stopifnot(is.list(list1), is.list(list2))
  
  n2 <- names(list2)
  
  for(k in n2){
    if(!k %in% names(list1)){
      list1 <- append( list1, list2[k] )
    }else{
      list1[k] <- list2[k]
    }
  }
  list1
}

gjamSimTime <- function(S, Q = 0, nsite, ntime = 50, termB, termR, termA, obsEffort = 1,
                        predPrey = NULL, zeroAlpha = NULL, PLOT = FALSE){
  
  # observations are counts ('DA' in gjam)
  # S           - no. species in Y
  # Q           - no. predictors in X
  # nsite       - no. groups/plots/sites
  # ntime       - mean length of time series
  # termB, termR, termA are logical for inclusion of immigration/emigration X%*%B, 
  #               DI growth VL, DD growth UA
  # obsEffort   - effort used in gjam
  # predPrey    - interactions assumed to be competition, unless predPrey = c(pred, prey), where
  #               prey and pred are the rows/columns in alpha 
  # zeroAlpha   - second column does not affect first column
  # PLOT        - logical to plot series
  
  cols <- colF( S )
  
  wdata <- NULL
  if(Q < 1)Q <- 1
  
  if(Q <= 1 | !termB){
    form <- '~ 1'
  }else{
    form <- paste0( 'x', c(2:Q), collapse = '+')
    form <- as.formula( paste(' ~', form ) )
  }
  
  if(!is.null(predPrey) & length(predPrey) == 2)   predPrey <- matrix(predPrey, nrow=1 )
  if(!is.null(zeroAlpha) & length(zeroAlpha) == 2)zeroAlpha <- matrix(zeroAlpha, nrow=1 )
  
  beta <- rhoTrue <- alphaTrue <- wstar <- NULL
  
  nt     <- 2 + rpois(nsite, ntime - 2)
  ntot   <- sum(nt)
  groups <- rep(1:nsite, nt)
  
  times  <- 1:nt[1]
  if(nsite > 1){
    for(k in 2:nsite)times <- c(times, c(1:nt[k]))
  }
  
  w <- matrix(1000,ntot,S)
  colnames(w) <- paste('s', 1:S, sep='')
  
  if(termB){ # environmental immigration/emigration
    
    bsd <- 5
    if(termR)bsd <- 5
    if(termA)bsd <- 20
    
    x     <- matrix( 0, ntot, Q)
    beta  <- matrix( rnorm(Q*S, 0, bsd), Q, S )
    beta[1,] <- rnorm(S, 0, 1)
    colnames(x) <- rownames(beta) <- paste('x',1:Q,sep='')
    rownames(beta)[1] <- 'intercept'
    colnames(beta) <- paste('s', 1:S, sep='')
    
  }else{
    
    x <- matrix(1, ntot, 1)
    colnames(x) <- 'intercept'
    
  }
  
  if(termR){ # growth rate rho
    
    gam <- runif(S, -.03, .03)       
    if(termB){
      gam <- runif(S, -.03, .05)
    }
    if(termA){
      gam <- runif(S, .04, .2)
    }
    rhoTrue <- gam
  }
  
  if(termA){ # alpha matrix
    
    intrWt <- 1.2
    if(S > 10)intrWt <- 10
    
    allPos <- F
    np     <- 0
    
    while(!allPos){
      daa <- runif(S,-1/100000,-1/200000)  # competition
      ir  <- runif(S^2, min(daa), 0)
      aa  <- matrix(ir, S, S)
      diag(aa) <- intrWt*daa
      
      if(!is.null(predPrey))aa[ predPrey ] <- -aa[ predPrey ] # prey has positive effect on pred
      if(!is.null(zeroAlpha))aa[ zeroAlpha ] <- 0             # no effect
      alphaTrue <- aa
      
      wstar <- -solve(crossprod(aa))%*%t(aa)%*%gam # carrying capacity
      if( all(wstar > 0) )allPos <- T
      np <- np + 1
    }
    print( paste(np, 'iterations for termA') )
  }
  
  # residual covariance small
  sigma <- diag(1, S)
  XB    <- 0
  
  for(k in 1:nsite){
    
    ww     <- matrix(10, nt[k], S)
    ww[1,] <- runif(S, 200, 400)
    if(!termR & !termA)ww[1,] <- beta[1,]
    if(termR & termA)  ww[1,] <- runif(S, wstar*.1, wstar*2)
    
    xx <- matrix( rnorm(nt[k]*Q), nt[k], Q)
    xx[,1] <- 1
    if(termB)XB <- xx%*%beta
    
    for(t in 2:nt[k]){
      
      ep    <- .rMVN(1, 0, sigma)
      ww[t,] <- ww[t-1,] + t( ep )
      if(termB) ww[t,] <- ww[t,] + XB[t,] 
      if(termR) ww[t,] <- ww[t,] + ww[t-1,]*gam
      if(termA) ww[t,] <- ww[t,] + diag(ww[t-1,])%*%aa%*%t(ww[drop=F,t-1,]) 
      
      if(sum(ww[t,]) > 1000000)stop('try again')
    }
    wk <- which(groups == k)
    w[wk,] <- ww
    if(termB)x[wk,] <- xx
  }
  
  wkeep <- which( is.finite(rowSums(w)) & apply(w, 1, min) > 0 )
  w <- w[wkeep,]
  if(termB){
    x <- x[wkeep,]
    x <- x[,!colnames(x) == 'x']
  }
  times <- times[wkeep]
  groups <- groups[wkeep]
  
  y <- round( w*obsEffort )
  colnames(y) <- colnames(w)
  
  if(PLOT){
    ylim <- c(0, 1.5*max(y))
    xlim <- c(0, max(nt))
    par(bty='n', cex=1.5)
    xlab <- expression( paste("Time ", italic(t), sep=''))
    ylab <- expression( paste("Count ", italic(y[s][t]), sep=''))
    
    wk <- which(groups == 1)
    plot(times[wk], y[wk,1],type='l', xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    
    for(k in 1:nsite){
      wk <- which(groups == k)
      for(s in 1:S){
        lines(times[wk], y[wk,s],col = cols[s],lwd=2)
        if(termA)abline(h=wstar[s]*obsEffort, lty=2, col= cols[s], lwd=2)
      }
    }
  }
  rownames(w) <- paste(groups, times, sep='-')
  
  if(termA){
    print( 'eigenvalues and carrying capacities' )
    wdata <- data.frame( eigenvalues = eigen(aa)$values, carryingCapacity = wstar )
    print( wdata )
  }
  
  rho <- matrix(0, S, S)
  diag(rho) <- rhoTrue
  
  trueValues <- list(beta = beta, rho = rho, alpha = alphaTrue, 
                     sigma = sigma, w = w)
  
  wt <- which( sapply(trueValues, length) > 0 )
  trueValues <- trueValues[ wt ]
  
  xdata <- data.frame( groups = groups, times = times, x, stringsAsFactors = F)
  
  list(xdata = xdata, ydata = y, edata = y*0 + obsEffort, formula = form, 
       groups = groups, times = times, trueValues = trueValues,
       wdata = wdata)
}

gjamTimePrior <- function( xdata, ydata, edata, priorList ){
  
  bp <- lp <- ap <- NULL
  formulaBeta <- formulaRho <- alphaSign <- NULL
  termB <- termR <- termA <- FALSE
  
  betaPrior <- rhoPrior <- alphaPrior <- NULL
  
  for(k in 1:length(priorList))assign( names(priorList)[k], priorList[[k]] )
  
  S <- ncol(ydata)
  w <- ydata/edata
  
  if(!is.null(betaPrior))termB <- TRUE
  if(!is.null(rhoPrior)) termR <- TRUE
  if(!is.null(alphaSign))termA <- TRUE
  
  if( 'formulaBeta' %in% names(priorList) )termB <- TRUE
  
  if(termB){  # only XB
    
    timeZero <- grep('-0', rownames(ydata)) # rownames from gjamFillMissing
    timeLast <- c( timeZero - 1, nrow(ydata))[-1]
    
    tmp <- model.frame(formulaBeta, data=xdata, na.action=NULL) # standardized design
    
    x   <- model.matrix(formulaBeta, data=tmp)
    
    # do not center intercept or factors
    wf <- names( which( attr(attributes(tmp)$terms, 'dataClasses') == 'factor' ) )
    wf <- which( startsWith(colnames(x), wf ) )
    wf <- c(1, wf)
    wp <- c(1:ncol(x))[-wf]
    
    if(length(wp) > 0){
      xm  <- colMeans(x, na.rm=T)
      xs  <- apply(x, 2, sd, na.rm=T)
      x[,wp]   <- t( (t(x[,wp]) - xm[wp])/xs[wp]  )
    }
    
    kname <- as.vector( outer( colnames(x), colnames(w), FUN = paste, sep='-') )
    
    xnames <- attributes(x)$dimnames[[2]]
    blo <- matrix(Inf, ncol(x), ncol(ydata))
    rownames(blo) <- colnames(x)
    colnames(blo) <- colnames(ydata)
    bhi <- -blo
    
    # variance in dy
    gr  <- unique(xdata$groups)
    nr  <- length(gr)
    
    sumw <- sumw2 <- rep(0, ncol(w))
    sumx <- sumx2 <- rep(0, ncol(x))
    
    for(j in gr){
      wj <- which(xdata[,'groups'] == j & is.finite(rowSums(x)) &
                    is.finite(rowSums(w)))
      wj <- wj[ !wj %in% timeZero ]
      
      if(length(wj) <= ncol(x)) next
      dw <- w[ wj[ -1 ], ] - w[ wj[ -length(wj) ], ]
      xj <- x[wj,][-1,]
      
      wv <- which( apply(xj, 2, var) > 0 )
      if(length(wv) == 0)next
      bb <- solve( crossprod(xj[,wv]) )%*%crossprod(xj[,wv], dw)
      
      bx <- blo
      bx[rownames(bb),] <- bb
      blo[bx < blo] <- bx[bx < blo]
      
      bx <- bhi
      bx[rownames(bb),] <- bb
      bhi[bx > bhi] <- bx[bx > bhi]
      
      sumw  <- sumw + colSums(dw)
      sumw2 <- sumw2 + colSums(dw^2)
      sumx  <- sumx + colSums(xj)
      sumx2 <- sumx2 + colSums(xj^2)
    }
    
    # variables that are fixed within a group
    wsd <- sqrt( sumw2/nrow(x) - (sumw/nrow(x))^2 )
    xsd <- sqrt( sumx2/nrow(x) - (sumx/nrow(x))^2 )
    xsd[1] <- 1
    
    wx  <- matrix(wsd, length(xsd), length(wsd), byrow=T)/matrix(xsd, length(xsd), length(wsd) )
    
    blo[ blo == 0 | !is.finite(blo) ] <- -wx[ blo == 0 | !is.finite(blo) ]
    bhi[ bhi == 0 | !is.finite(bhi) ] <-  wx[ bhi == 0 | !is.finite(bhi) ]
    
    bl <- blo - 2*abs(blo)
    bh <- bhi + 2*abs(bhi)
    rownames(bl)[1] <- rownames(bh)[1] <- 'intercept'
    
    if( is.list(betaPrior) ){
      gg  <- gjamPriorTemplate(formulaBeta, xdata, ydata = ydata, 
                               lo = betaPrior$lo, hi = betaPrior$hi)
      blo <- gg[[1]]
      bhi <- gg[[2]]
      
      blo[ !rownames(blo)%in% names(betaPrior$lo), colnames(blo) ] <- 
        bl[ !rownames(blo)%in% names(betaPrior$lo),  ]
      
      bhi[ !rownames(bhi)%in% names(betaPrior$hi), colnames(bhi) ] <- 
        bh[ !rownames(bhi)%in% names(betaPrior$hi),  ]
      
    }else{
      blo <- bl
      attr(blo,'formula') <- formulaBeta
      bhi <- bh
    }
    bp  <- list(lo = blo, hi = bhi )
  }
  
  if(termR){  # rho
    if( !is.list(rhoPrior) ){
      lp <- NULL
    }else{
      gg <- gjamPriorTemplate(formulaRho, xdata, ydata = ydata, 
                              lo = rhoPrior$lo, hi = rhoPrior$hi)
      lp <- list(lo = gg[[1]], hi = gg[[2]])
    }
  }
  if(termA){  # alpha
    
    if( !is.list(rhoPrior) ){
      if(is.null(lp))stop(' must have rhoPrior if there is alphaSign' )
      ap <- NULL
    }else{
      
      rho <- (lp$lo['intercept',] + lp$hi['intercept',1])/2 + .01
      
      wstar <- colMeans( ydata/edata, na.rm=T )
      aa <- (1 - rho)/wstar                  # crude carrying capacity
      
      aa <- rho/wstar
      
      a1 <- matrix(-aa, S, S)
      a2 <- matrix(-aa, S, S, byrow=T)
      a1[ a2 > a1 ] <- a2[ a2 > a1 ]
      
      a1 <- 2*a1
      
      #alphaSign
      
      alo <- a1
      ahi <- alo*0
      ww  <- which(alphaSign < 0)
      if(length(ww) > 0){
        alo[ww] <- a1[ww]
        ahi[ww] <- 0
      }
      ww <- which(alphaSign > 0)
      if(length(ww) > 0){
        alo[ww] <- 0
        ahi[ww] <- -a1[ww]
      }
      ww <- which(alphaSign == 0)
      if(length(ww) > 0)alo[ww] <- ahi[ww] <- NA
      
      ap <- list(lo = alo, hi = ahi)
    }
  }
  
  list(betaPrior = bp, rhoPrior = lp, alphaPrior = ap)
}

foodWebDiagram <- function(S, guildList = NULL, predPrey = NULL, zeroAlpha = NULL,
                           intraComp = 1:S, label = NULL, PLOT = TRUE, layout = 'rr'){
  
  # S - no. species
  # default interaction is negative, arrows only for negative interactions
  # guildList - overrides default, only members of the same guild compete
  # predPrey  - matrix with 2 columns, second column is prey of first column
  # zeroAlpha - matrix with 2 columns, second column does not affect first column
  # intraComp - needed for intraspecific comp if guildList is specified
  # layout can be 'tree', 'rr', ...
  
  require( DiagrammeR )
  
  pp <- numeric(0)
  
  fromTo <- as.matrix( expand.grid( c(1:S), c(1:S) ) )
  ft <- .pasteCols( fromTo )
  
  qq <- numeric(0)
  if( !is.null(guildList) ){
    
    fromTo <- numeric(0)
    for(k in 1:length(guildList)){
      ft <- as.matrix(expand.grid(guildList[[k]], guildList[[k]]))
      fromTo <- rbind(fromTo, ft)
    }
    if(!is.null(intraComp)){
      ft <- cbind(1:S, 1:S)
      fromTo <- rbind(fromTo, ft)
    }
    ft <- .pasteCols( fromTo )
    ww <- which(!duplicated(ft))
    ft <- ft[ww]
    fromTo <- fromTo[ww,]
    zeroAlpha <- NULL
  }
  
  if(!is.null(predPrey)){
    pp <- .pasteCols( predPrey[drop=F,,c(2,1)] )
    qq <- .pasteCols( predPrey)
    #   fromTo <- fromTo[ !ft %in% pp, ]
    fromTo <- rbind(fromTo, predPrey, predPrey[drop=F,,c(2,1)])
    ft <- .pasteCols( fromTo )
    ww <- which(!duplicated(ft))
    fromTo <- fromTo[ww,]
    ft <- ft[ww]
  }
  
  if(!is.null(zeroAlpha)){
    za <- .pasteCols( zeroAlpha[drop=F,,c(1,2)])
    
    fromTo <- fromTo[ !ft %in% za, ]
    ft <- ft[ !ft %in% za ]
  }
  
  if( length(qq) > 0 ){
    ft <- .pasteCols( fromTo )
    fromTo <- rbind( fromTo[ft %in% qq,], fromTo[!ft %in% qq, ] )
  }
  ft <- .pasteCols( fromTo )
  pp <- pp[ pp %in% ft ]
  
  if(!PLOT){
    return( fromTo )
  }else{
    ecol  <- rep( 'tan', length(ft) )
    #  ncol  <- rep( 'tan', S )
    ncol  <- colF(S)
    shape <- rep( 'rectangle', S)
    
    if( length(qq) > 0 ){
      
      wt <- which( ft %in% qq )
      ecol[ wt ] <- 'brown'
      
    }
    if( length(pp) > 0 ){
      
      wt <- which( ft %in% pp )
      ecol[ wt ] <- 'blue'
    }
    
    if( is.null(layout) )layout <- "tree"
    if(is.null(label))label <- paste('s', 1:S, sep='')
    nodes <- create_node_df(n = S, label = label, style = "filled", color = ncol, shape = shape)
    edges <- create_edge_df(from = fromTo[,1], to = fromTo[,2], color = ecol )
    graph <- create_graph(nodes_df = nodes, edges_df = edges)
    render_graph( graph,  layout = layout )
  }
}


