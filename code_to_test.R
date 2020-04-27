#########code for testing the gjamplot, will return error. 

outputBeetles   <- readRDS('min_obs_100.rds')
plotPars        <- list(PLOTALLY=T, SAVEPLOTS = F, outFolder = 'gjamOutputLakes1')
.gjamPlot(outputBeetles,plotPars)


#########code for testing buildguildstrait, yname not defined

traits          <- read.csv("../makeMastOnJimClark/traitsByGroup/carabidTraits.csv", 
                            stringsAsFactors = F)
#tmp             <- readRDS('data.rds')

tmp             <- readRDS('data_min_obs_70.rds')
xdata           <- tmp$xdata
ydata           <- tmp$ydata
edata           <- tmp$edata
timeList        <- tmp$timeList


######################
min_obs         <- 80
ydata           <- gjamTrimY(ydata, minObs = min_obs, OTHER = FALSE)$y
ynames          <- colnames(ydata)

S               <- ncol(ydata)
n               <- nrow(ydata)

xdata$u1        <- sin(xdata$slope)
xdata$u2        <- sin(xdata$slope)*sin(xdata$aspect)
xdata$u3        <- sin(xdata$slope)*cos(xdata$aspect)

formulaBeta     <- as.formula(~ deficit.JJA + tmmn.DJF + u1 + u2 + u3)
formulaRho      <- as.formula(~ deficit.JJA + tmmn.DJF + u1 + u2 + u3)

#############################
formulaBeta     <- as.formula(~ deficit.JJA + tmmn.DJF )
formulaRho      <- as.formula(~ deficit.JJA + tmmn.DJF)


traitNames      <- c('trophic','activity','burrow','climb')

gtmp            <- buildGuildTraits(ydata, traitTable = traits, traitNames, 
                                    traitCode = 'code6', minCooccur = 5,
                                    diagTrait = -1)
guildList       <- gtmp$guildList
zeroAlpha       <- gtmp$zeroAlpha

foodWebDiagram(S, guildList = guildList, label = ynames, zeroAlpha = zeroAlpha,
               intraComp = 1:S,  layout='rr')
			   
			   
			   
species    <- colnames(ydata)
edata      <- edata[,species] %>% as.data.frame()


lo         <- list(intercept = -.01)
hi         <- list(intercept = .3)
rhoPrior   <- list(lo = lo, hi = hi)

lo         <- list(intercept = -2)
hi         <- list(intercept = 2)
betaPrior  <- list(lo = lo, hi = hi)



alphaSign              <- matrix(-1, S, S)
colnames(alphaSign)    <- rownames(alphaSign) <- ynames
alphaSign[ zeroAlpha ] <- 0

priorList              <- list( formulaBeta = formulaBeta, formulaRho = formulaRho,
                                betaPrior = betaPrior, rhoPrior = rhoPrior, alphaSign = alphaSign)
tmp                    <- gjamTimePrior( xdata, ydata, edata, priorList)

timeList               <- mergeList(timeList, tmp)

effort     <- list(columns = 1:S, values = edata)
modelList  <- list( typeNames = 'DA', ng = 500, burnin=200,
                    timeList=timeList, effort = effort)
outputBeetles  <- .gjam(formulaBeta, xdata=xdata, ydata=ydata, modelList=modelList)
saveRDS(outputBeetles, sprintf('min_obs_%d.rds', min_obs))
plotPars       <- list(PLOTALLY=T, SAVEPLOTS = F, outFolder = sprintf('gjamOutput_min_obs_%d.rds', min_obs))
.gjamPlot(outputBeetles,plotPars)			  





