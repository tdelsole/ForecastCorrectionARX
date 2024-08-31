rm(list=ls())

lplotfile     = FALSE
region        = '60S60N'
model.forcing = 'gfdl'
dir.figs      = './figures/'

forcing.type  = 'annual.pooled'
forcing.type  = 'monthly.pooled'
forcing.type  = 'annual.lagged'
nlag.forcing  = 4  ### NLAG.FORCING = 1 MEANS 0 LAG; IGNORED FOR MONTHLY FORCING

p.order       = 2
nharm         = 5
mperiod       = 1

train.str     = as.Date('15-01-1991',format='%d-%m-%Y')
train.end     = as.Date('15-12-2023',format='%d-%m-%Y')
train.end     = as.Date('15-12-2020',format='%d-%m-%Y')

dir.Rlib      = '/Users/delsole/documents/R/delsole_tools/'
dir.Rlib      = NULL


source(paste(dir.Rlib,'pdf.eps.R',sep=''))
source(paste(dir.Rlib,'timeseries2arx.cyclo.R',sep=''))
source(paste(dir.Rlib,'n.to.monthly.R',sep=''))
source(paste(dir.Rlib,'traditional.correction.R',sep=''))
source(paste(dir.Rlib,'simulate.arx.cyclo.R',sep=''))
source(paste(dir.Rlib,'arx.correction.R',sep=''))
source(paste(dir.Rlib,'lm.TrendPlusCycle.R',sep=''))

library(lubridate)
library(xtable)
library(rworldmap)


#####################################
######### META DATA
#####################################
if (forcing.type == 'annual.pooled') {
	forcing.suffix = '.12mrunning'; forcing.part = 'pooled'
} else if (forcing.type == 'annual.lagged') {
	forcing.suffix = '.12mrunning'; forcing.part = 'lagged'
} else if (forcing.type == 'monthly.pooled') {
	forcing.suffix = NULL; forcing.part = 'skip'; nlag.forcing = 1
} else stop("do not recognize forcing.type")

################################################################
############ READ DATA
################################################################
# fread = paste('../RCodes/figures/data.SPEAR.',forcing.type,'.RData',sep='')
fread = paste('data.SPEAR.',forcing.type,'.RData',sep='')
load(fread)

################################################################
################################################################
### KEY DATA: 
### 	pc.nmme     [nlead.nmme,nstart.nmme,emax.nmme]				NMME SPEAR TIME SERIES
###		targ.nmme   [nlead.nmme,nstart.nmme,emax.nmme]				NMME SPEAR TARGET DATES
### 	pc.hist     [ntot.hist,emax.hist]							LENS SPEAR TIME SERIES
###		targ.hist   [ntot.hist,emax.hist]							LENS SPEAR TARGET DATES
###		pc.forcing  [ntot.forcing,nfor]								SPEAR FORCING TIME SERIES (NFOR = 3, GHG,AER,NAT)
###     targ.forcing[ntot.forcing]									SPEAR FORCING TARGET DATES
###     pc.hadcrut  [ntot.hadcrut]									HADCRUT TIME SERIES
###     targ.hadcrut[ntot.hadcrut]									HADCRUT TARGET DATES
###     pc.noaaobs  [ntot.noaaobs]									NOAAOBS TIME SERIES
###     targ.noaaobs[ntot.noaaobs]									NOAAOBS TARGET DATES
###     pc.best     [ntot.best]										BERKELEY TIME SERIES
###     targ.best   [ntot.best]										BERKELEY TARGE DATES
###     pc.best.clim[12]											BERKELEY CLIMATOLOGY
###     region[1]													NAME OF REGION
################################################################
################################################################


#################################
########## DEFINE METADATA
#################################
emax.nmme     = spear.list$emax.nmme
nlead.nmme    = spear.list$nlead.nmme
nstart.nmme   = length(spear.list$pc.nmme)/nlead.nmme/emax.nmme
emax.hist     = spear.list$emax.hist
ntot.hist     = length(spear.list$targ.hist)/emax.hist
forcing.say   = spear.list$forcing.say
nforcing      = length(forcing.say)

earliest.nmme = spear.list$targ.nmme[which.min(spear.list$targ.nmme)]
latest.nmme   = spear.list$targ.nmme[which.max(spear.list$targ.nmme)]

#####################################################################################
### COMPUTE ENSEMBLE MEAN AND ENSEMBLE STANDARD DEVIATION
#####################################################################################
dim(spear.list$pc.hist)      = c(ntot.hist,emax.hist)
dim(spear.list$pc.nmme)      = c(nlead.nmme*nstart.nmme,emax.nmme)
hist.emean                   = rowMeans(spear.list$pc.hist)
nmme.emean                   = rowMeans(spear.list$pc.nmme)
hist.stdv                    = sqrt(rowMeans((spear.list$pc.hist - hist.emean)^2))
nmme.stdv                    = sqrt(rowMeans((spear.list$pc.nmme - nmme.emean)^2))
dim(nmme.stdv)               = c(nlead.nmme,nstart.nmme)

#################################
########## RESCALE FORCINGS (ERF, W/m2)
#################################
scale.forcing          = 0.1
spear.list$pc.forcing  = scale.forcing * spear.list$pc.forcing

##################################################################
########## ADD CLIMATOLOGY TO BERKELEY DATA SET; CLIP TO INTEGER MULTIPLE OF 12
##################################################################
dum = spear.list$targ.best
if (length(dum) > 0) {
	spear.list$pc.best   = spear.list$pc.best + spear.list$pc.best.clim[month(spear.list$targ.best)] + 273.15
	nstop                = floor(length(spear.list$pc.best)/12)*12
	spear.list$pc.best   = spear.list$pc.best[1:nstop]
	spear.list$targ.best = spear.list$targ.best[1:nstop]
	ntot.best            = length(spear.list$targ.best)
}


###################################################################
### PLOT PC.FORCING TIME SERIES
###################################################################
ncol.plot = 3
nrow.plot = ceiling(nforcing/ncol.plot)
icnt      = 0
fout      = paste(dir.figs,'SPEAR.',model.forcing,'.forcing.ts',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
par(mfrow=c(nrow.plot,ncol.plot),mar=c(3,3,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
for (nf in 1:nforcing) {
	icnt = icnt + 1
	plot(spear.list$targ.forcing,spear.list$pc.forcing[,nf],type='l',col=icnt)
	title(main=paste(forcing.say[nf]))
}
if (lplotfile) dev.off()	

##################################################################
##################################################################
########## PLOT RAW TIME SERIES: VERIFY ANNUAL CYCLES LINE UP
##################################################################
##################################################################
dim(spear.list$targ.nmme) = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(spear.list$pc.nmme)   = c(nlead.nmme,nstart.nmme,emax.nmme)
yrange                    = range(spear.list$pc.nmme,spear.list$pc.best)
earliest.plot             = as.Date('15-01-1990',format='%d-%m-%Y')
latest.plot               = as.Date('15-01-1994',format='%d-%m-%Y')
pos.mon                   = c(1,2,3,4)

fout = paste(dir.figs,'SPEAR.RawTS',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
plot(c(earliest.plot,latest.plot),c(1,1),type='n',ylim=yrange,xlab='',ylab='Average 2m-Temperature')

for (ne in 1:1) for (ny in 1:5) for (nm in c(1,4,7,10)) {
	year.pic = year(earliest.nmme) + ny - 1
	npic     = which(spear.list$targ.nmme[1,,1] == as.Date(paste('15',nm,year.pic,sep='-'),format='%d-%m-%Y'))
	if (length(npic) == 1) {
		lines(spear.list$targ.nmme[,npic,ne],spear.list$pc.nmme[,npic,ne],col=(nm-1)/3+2,lwd=2)
		if (ne == 1 & ny == 1) text(spear.list$targ.nmme[1,npic,ne],spear.list$pc.nmme[1,npic,ne],month.abb[nm],col=(nm-1)/3+2,pos=pos.mon[(nm-1)/3+1])
	}		
}

for (ne in 1:1) lines (spear.list$targ.hist[,ne],spear.list$pc.hist[,ne],col='grey50',lwd=3)
for (ne in 1:1) points(spear.list$targ.hist[,ne],spear.list$pc.hist[,ne],col='grey50',pch=19)
lines(spear.list$targ.best,spear.list$pc.best,lwd=2,col='black')
legend.say = c('LENS','Berkeley','NMME Jan Start','NMME Apr Start','NMME Jul Start','NMME Oct Start')
legend.col = c('grey50','black',2,3,4,5)
legend.pch = c(19,rep(NA,6))
legend('topleft',legend=legend.say,lwd=2,col=legend.col,pch=legend.pch,bg='white',ncol=3)
# legend('topleft',legend=c('LENS','NMME (Jan,Apr,Jul,Oct)','Berkeley'),lwd=2,col=c('grey50','red','black'),pch=c(19,NA,NA),bg='white')

ftitle.top = paste('SPEAR Large-Ensemble and NMME Forcasts and Hindcasts')
ftitle.bot = paste('2mT',region,'(1 ensemble member)')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()


####################################################################################
####################################################################################
################## COMPUTE TRADITIONAL CLIMATOLOGIES 
####################################################################################
####################################################################################
source(paste(dir.Rlib,'traditional.correction.R',sep=''))
old.nmme.best    = traditional.correction(spear.list$pc.nmme,spear.list$targ.nmme,spear.list$pc.best,spear.list$targ.best,
                     nlead.nmme,nstart.nmme,emax.nmme,train.str=train.str,train.end=train.end)
old.nmme.hadcrut = traditional.correction(spear.list$pc.nmme,spear.list$targ.nmme,spear.list$pc.hadcrut,spear.list$targ.hadcrut,
                     nlead.nmme,nstart.nmme,emax.nmme,train.str=train.str,train.end=train.end)

best.old.cor = old.nmme.best$obs.anom
nmme.old.cor = old.nmme.best$forecast.anom
anom.hadcrut = old.nmme.hadcrut$obs.anom
ntot.hadcrut = length(spear.list$targ.hadcrut)
err.nmme.old = old.nmme.best$err.corrected

########## REMOVE CLIMATOLOGY OF HISTORICAL
if (ntot.hist %% 12 != 0) stop('historical is not an integral multiple of 12')
dim(spear.list$targ.hist) = c(12,ntot.hist/12*emax.hist)
dim(spear.list$pc.hist)   = c(12,ntot.hist/12*emax.hist)
lgood                     = spear.list$targ.hist >= train.str & spear.list$targ.hist <= train.end
dum                       = spear.list$pc.hist
dum[!lgood]               = NA
anom.hist                 = spear.list$pc.hist - rowMeans(dum,na.rm=TRUE)

####################################################################################
####################################################################################
################## COMPUTE ARX CLIMATOLOGIES 
####################################################################################
####################################################################################
source(paste(dir.Rlib,'arx.correction.R',sep=''))
arx.nmme.best = arx.correction(spear.list$pc.nmme,spear.list$targ.nmme,spear.list$pc.best,spear.list$targ.best,
                     nforcing,spear.list$pc.forcing,spear.list$targ.forcing,spear.list$forcing.say,
                     nlead.nmme,nstart.nmme,emax.nmme,p.order,nharm,train.str,train.end)
                     
n.best.nmmeformat     = arx.nmme.best$n.obs.forecast.format
pc.forcing.bestformat = arx.nmme.best$pc.forcing.obs.format
beta.nmme.stat        = arx.nmme.best$beta.forecast
beta.best.stat        = arx.nmme.best$beta.obs
nvar.best.stat        = arx.nmme.best$nvar.obs
nvar.nmme.stat        = arx.nmme.best$nvar.forecast
lm.best.stat          = arx.nmme.best$lm.obs
lm.nmme.stat          = arx.nmme.best$lm.forecast
nmme.arx.cor          = arx.nmme.best$forecast.arx.anom
best.arx.cor          = arx.nmme.best$o.forecast.format.anom
err.uncorrected       = arx.nmme.best$err.uncorrected
err.nmme.arx          = arx.nmme.best$err.corrected 
simlist.best.nmmeformat.stat = arx.nmme.best$simlist.obs
simlist.nmme.nmmeformat.stat = arx.nmme.best$simlist.forecast





#################################
########## PLOT ANOMALY TIME SERIES
#################################
dim(nmme.old.cor)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(anom.hist)            = c(ntot.hist,emax.hist)
dim(best.old.cor)         = c(ntot.best)
dim(anom.hadcrut)         = c(ntot.hadcrut)
yrange                    = range(nmme.old.cor,anom.hist)
yrange                    = c(-1,1)

dim(spear.list$targ.hist) = c(ntot.hist,emax.hist)
dim(spear.list$targ.best) = c(ntot.best)
dim(spear.list$pc.best)   = c(ntot.best)

earliest.plot             = as.Date('15-01-1960',format='%d-%m-%Y')
latest.plot               = as.Date('15-01-2024',format='%d-%m-%Y')

year.plot.str             = 1991
year.plot.end             = 2020

fout = paste(dir.figs,'SPEAR.AnomalyTS',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.2)
plot(c(earliest.plot,latest.plot),c(1,1),type='n',ylim=yrange,xlab='',ylab='2mT Anomaly')

for (ne in 1:1) for (y.base in c(year.plot.str,year.plot.end)) for (ny in c(1:5)) for (nm in c(1,4,7,10)) {
	year.pic = y.base + ny - 1
	npic     = which(spear.list$targ.nmme[1,,1] == as.Date(paste('15',nm,year.pic,sep='-'),format='%d-%m-%Y'))
	if (length(npic) == 1) lines(spear.list$targ.nmme[,npic,ne],nmme.old.cor[,npic,ne],col='red',lwd=2)
}

# for (ne in 1:1) lines(spear.list$targ.hist[,ne],anom.hist[,ne],col='blue',lty='solid',lwd=2)
lines(spear.list$targ.best,best.old.cor,lwd=2)
# points(spear.list$targ.hadcrut,anom.hadcrut,pch=19,col='darkgreen',cex=0.4)

# legend('topleft',legend=c('NMME','LENS','Berkeley','HadCRUT'),lwd=c(2,2,2,NA),pch=c(NA,NA,NA,19),col=c('blue','red','black','darkgreen'),cex=1.3)
# legend('topleft',legend=c('NMME','LENS','Berkeley'),lwd=c(2,2,2,NA),pch=c(NA,NA,NA,19),col=c('red','blue','black'),cex=1.3)
legend('topleft',legend=c('SPEAR','Berkeley'),lwd=c(2,2,2,NA),pch=c(NA,NA,19),col=c('red','black'),cex=1.3)
# ftitle.top = paste("Anomalies of NMME, LENS, Observations (Traditional Bias Correction)")
ftitle.top = paste("Anomalies of SPEAR and Observations (Traditional Bias Correction)")
ftitle.bot = paste('2mT',region,'(Jan, Apr, Jul, Oct Starts; 1 ensemble member)')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()

##################################################################
##################################################################
########## INITIAL ERRORS
##################################################################
##################################################################
err.init              = err.uncorrected
dim(err.init)         = c(nlead.nmme,nstart.nmme,emax.nmme)
err.init              = err.init[1,,]
err.init.emean        = rowMeans(err.init)

err.init.emean.targ   = spear.list$targ.nmme[1,,1]
trend.init.error.nmme = lm.TrendPlusCycle(err.init.emean,err.init.emean.targ,train.str,train.end)
trend.init.err        = round(c(trend.init.error.nmme$trend, 2*trend.init.error.nmme$sterr),3)
sigma.init.err        = signif(trend.init.error.nmme$sigma,2)
y.pred                = trend.init.error.nmme$y.pred


##################################################################
########## PLOT LEAD-1 ERRORS: TIME SERIES
##################################################################
fout = paste(dir.figs,'SPEAR.',model.forcing,'.lead1error',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
iyst.plot  = year(err.init.emean.targ[1])
iynd.plot  = year(err.init.emean.targ[length(err.init.emean.targ)])
yrange     = range(err.init.emean,na.rm=TRUE)
xrange     = c(iyst.plot,iynd.plot+1)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='',ylab='error')
for (m in 1:12) {
	nget = which( month(err.init.emean.targ) == m)
	lines(year(err.init.emean.targ)[nget],err.init.emean[nget],col=m,lwd=2)
	lines(year(err.init.emean.targ)[nget],y.pred[nget],col=m)
	# text (iynd.plot-1,err.init.emean[nget[length(nget)-1]],month.abb[m],pos=4,col=m)
	text (iynd.plot-1,y.pred[nget[length(nget)-1]],month.abb[m],pos=4,col=m)
}
model.stuff = bquote(.(trend.init.err[1])%+-%.(trend.init.err[2])*' C/decade; StErr='~.(sigma.init.err)*'C')
text(iynd.plot,yrange[1],model.stuff,pos=2,cex=1.1)
ftitle.top = paste("Initial SPEAR minus Berkeley (lead = 0.5; ensemble mean)")
ftitle.bot = paste('ARX parameters: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()


#####################################################################################
########## PLOTS OF ARX SIMULATIONS TRAINED ON BERKELEY AND NMME: NO ANNUAL CYCLE
#####################################################################################
beta.nmme.nocycle = beta.nmme.stat
beta.nmme.nocycle[p.order+1:(2*nharm)] = 0
simlist.nmme.nocycle = as.numeric(simulate.arx.cyclo(
   spear.list$pc.best[1:p.order],
   spear.list$targ.best,
   pc.forcing.bestformat,forcing.say,
   beta.nmme.nocycle,ntot.best,p.order,nharm,mperiod))

beta.best.nocycle = beta.best.stat
beta.best.nocycle[p.order+1:(2*nharm)] = 0
simlist.best.nocycle = as.numeric(simulate.arx.cyclo(
   spear.list$pc.best[1:p.order],
   spear.list$targ.best,
   pc.forcing.bestformat,forcing.say,
   beta.best.nocycle,ntot.best,p.order,nharm,mperiod))
  


earliest.plot  = as.Date('15-01-1990',format='%d-%m-%Y')
latest.plot    = as.Date('15-12-2023',format='%d-%m-%Y')
nst.best = which(spear.list$targ.best      == earliest.plot)
nnd.best = which(spear.list$targ.best      == latest.plot)
fout = paste(dir.figs,'SPEAR.',model.forcing,'.simBEST&NMME.stat',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
nmme.offset = 1.5
yrange = range(simlist.nmme.nocycle[nst.best:nnd.best]+ nmme.offset,simlist.best.nocycle[nst.best:nnd.best])
yrange[2] = yrange[2]
plot (spear.list$targ.best[nst.best:nnd.best  ], simlist.nmme.nocycle[nst.best:nnd.best]+ nmme.offset,type='l',xlab='',ylab='2mT (K)',col='red',lwd=2,ylim=yrange)
lines(spear.list$targ.best[nst.best:nnd.best  ], simlist.best.nocycle[nst.best:nnd.best],type='l',lwd=2,col='black')
legend('topleft',legend=c('Noise-Free ARX Trained on Berkeley',paste("Noise-Free ARX Trained on SPEAR +",nmme.offset)),col=c('black','red'),lwd=2,cex=1.3,bg='white')
ftitle.top = paste("Noise-Free ARX Model trained on Berkeley and SPEAR-NMME")
ftitle.bot = paste('w/o Annual cycle; ARX parameters: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()

##################################################################
##################################################################
########## PLOT REGRESSION COEFFICIENTS AND THEIR STANDARD ERRORS
##################################################################
##################################################################
ncases = 2
case.say = character(ncases)
FIRST  = TRUE
for (nc in 1:ncases) {
	if (nc == 1) {
		lm.pic = lm.best.stat; nvar.pic = nvar.best.stat; case.say[nc] = 'Berkeley'
	} else if (nc == 2) {
		lm.pic = lm.nmme.stat; nvar.pic = nvar.nmme.stat; case.say[nc] = 'SPEAR-NMME'
	} else stop('do not understand ncase')
	
	if (FIRST) {
		FIRST = FALSE
		beta       = lm.pic$coef
		beta.names = names(beta)
		nbeta      = length(beta)
		beta.all   = array(NA,dim=c(nbeta,ncases))
		stdv.all   = array(NA,dim=c(nbeta,ncases))
		nvar.all   = as.numeric(rep(NA,ncases))
		ipic       = which(endsWith(beta.names,'intercept'))
		if (length(ipic) != 1) stop('cannot find intercept')
	}
	beta.all[,nc] = coef(lm.pic)
	stdv.all[,nc] = summary(lm.pic)[['coefficients']][,'Std. Error']
	nvar.all[ nc] = nvar.pic
}

yrange = range(beta.all[-ipic,] + 2 * stdv.all[-ipic,], beta.all[-ipic,] - 2 * stdv.all[-ipic,])

colnames(beta.all) = case.say
rownames(beta.all) = names(beta.best.stat)
names(nvar.all)    = case.say

print(round(rbind(beta.all,nvar.all),4))

dev.type = 'coef'
caption.say   = paste(dev.type,'; H= ',nharm,'; P=',p.order,sep='')
ftable = paste(dir.figs,'table.',model.forcing,'.',dev.type,'.tex',sep='')
table.say   = xtable(rbind(beta.all,nvar.all),caption=caption.say,label=paste('tab:',dev.type,sep=''))
if (lplotfile) capture.output(table.say,file=ftable)


### PLOT
intercept.scale = 0.01
beta.plot = beta.all
stdv.plot = stdv.all
intercept.pic = which(rownames(beta.plot) == 'intercept')
beta.plot[intercept.pic,] = beta.plot[intercept.pic,]*intercept.scale
stdv.plot[intercept.pic,] = stdv.plot[intercept.pic,]*intercept.scale
beta.names[1:p.order] = paste('lag',1:p.order,sep='')

fout = paste(dir.figs,'SPEAR.',model.forcing,'.beta.stat',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
xoff   = 0.15
xshift = seq(from=-xoff,to=xoff,length.out=ncases) 
par(mfcol=c(1,1),mar=c(8,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
plot(1,1,type='n',xlim=c(1,length(beta.names)),ylim=yrange,xaxt='n',xlab='',ylab='regression coefficient')
axis(1,at=1:nbeta,beta.names,las=2)
for (nc in 1:ncases) {
	points(1:nbeta+xshift[nc],beta.plot[,nc],pch=19,cex=0.8,col=nc)
	arrows(1:nbeta+xshift[nc],beta.plot[,nc]+2 * stdv.plot[,nc],1:nbeta+xshift[nc],beta.plot[,nc]-2 * stdv.plot[,nc],code=3,angle=90,length=0.02,col=nc)
}
legend('right',legend=case.say,pch=19,col=1:ncases,cex=1.5)
abline(h=0,lwd=2,col='grey50')
abline(v=p.order+0.5,col='grey70',lty='dashed')
abline(v=p.order + 2*nharm + 0.5,col='grey70',lty='dashed')
abline(v=p.order + 2*nharm + 1 + 0.5, col='grey70',lty='dashed')
text((p.order+1)/2,yrange[1],'AR',col='blue')
text(p.order + (2*nharm+1)/2,yrange[1],'Annual Cycle',col='blue')
text(p.order + 2*nharm + 1 + (nforcing+1)/2,yrange[1],'GHG/AER/NAT',col='blue')
ftitle.top = paste("Coefficients of ARX")
ftitle.bot = paste('ARX Fit: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()

###########################################
############# PLOT ANOMALIES
###########################################
dim(nmme.old.cor)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(nmme.arx.cor)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(best.arx.cor)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(spear.list$targ.nmme) = c(nlead.nmme,nstart.nmme,emax.nmme)

yrange                    = c(-1,1)


for (lead.pic in c(1,p.order + 1)) for (mpic in c(1,4,7,10)) {
	fout = paste(dir.figs,'SPEAR.',model.forcing,'.lead',lead.pic,'old.anom.',month.abb[mpic],sep='')
	if (lplotfile) pdf.eps(fout,'pdf')
	par(mfcol=c(1,1),mar=c(8,5,3,1))
	par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)

	nget       = month(spear.list$targ.nmme[1,,1]) == mpic 
	nget.clim  = nget & spear.list$targ.nmme[1,,1] >= train.str & spear.list$targ.nmme[1,,1] <= train.end
	nbest      = month(spear.list$targ.best) == (mpic -1 + lead.pic - 1) %% 12 + 1 
	nbest.clim = nbest & spear.list$targ.best >= train.str & spear.list$targ.best <= train.end
	# if (abs(mean(best.old.cor[nbest.clim])) > 1.e-10) stop('Berkeley climatology has nonzero mean')
	# if (abs(mean(nmme.old.cor[lead.pic,nget.clim,])) > 1.e-10) stop('NMME climatology has nonzero mean')

	yrange = range(nmme.old.cor[lead.pic,nget,])
	plot  (spear.list$targ.nmme[lead.pic,nget,],nmme.old.cor[lead.pic,nget,],cex=0.5,xlab='',ylab='anomaly (degrees C)',col='red',ylim=yrange,pch=19)
	points(spear.list$targ.best[nbest],best.old.cor[nbest],col='black',cex=1,pch=17)

	y           = best.old.cor[nbest.clim]
	x           = year(spear.list$targ.best[nbest.clim])/10
	best.sum    = summary(lm(y~x))
	y           = as.numeric(nmme.old.cor[lead.pic,nget.clim,])
	x           = as.numeric(year(spear.list$targ.nmme[lead.pic,nget.clim,]))/10
	nmme.sum    = summary(lm(y~x))
	best.trend  = round(best.sum$coef[2,1:2],2)
	nmme.trend  = round(nmme.sum$coef[2,1:2],2)
	legend.best = bquote('Berkeley ('*.(best.trend[1])%+-%.(best.trend[2])*' C/decade)')
	legend.nmme = bquote('SPEAR   ('*.(nmme.trend[1])%+-%.(nmme.trend[2])*' C/decade)')
	
	x = as.numeric(spear.list$targ.nmme[lead.pic,nget.clim,])
	y = as.numeric(nmme.old.cor[lead.pic,nget.clim,])
	abline(lm(y~x),col='red',lwd=2)
	x = as.numeric(spear.list$targ.best[nbest.clim])
	y = as.numeric(best.old.cor[nbest.clim])
	abline(lm(y~x),col='black',lwd=2)

	ftitle.top = paste("Lead-",lead.pic-0.5," Anomalies Based on Traditional Bias Correction",sep='')
	ftitle.bot = paste('Target is ',month.abb[mpic],'; ',region,sep='')
	title(ftitle.top,line=2.0)
	title(ftitle.bot,line=0.5)
	abline(h=0,col='grey',lty='dashed',lwd=2)
	legend('topleft',legend=c(legend.nmme,legend.best),pch=c(19,17),col=c('red','black'),cex=1.2)
	if (lplotfile) dev.off()
}


lead.pic = p.order + 1
for (mpic in c(1,4,7,10)) {
	fout = paste(dir.figs,'SPEAR.',model.forcing,'.lead',lead.pic,'arx.anom.',month.abb[mpic],sep='')
	if (lplotfile) pdf.eps(fout,'pdf')
	par(mfcol=c(1,1),mar=c(8,5,3,1))
	par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)

	nget       = month(spear.list$targ.nmme[1,,1]) == mpic 
	nget.clim  = nget & spear.list$targ.nmme[1,,1] >= train.str & spear.list$targ.nmme[1,,1] <= train.end
	nbest      = month(spear.list$targ.best) == (mpic -1 + lead.pic - 1) %% 12 + 1 
	nbest.clim = nbest & spear.list$targ.best >= train.str & spear.list$targ.best <= train.end
	# if (abs(mean(best.old.cor[nbest.clim])) > 1.e-10) stop('Berkeley climatology has nonzero mean')
	# if (abs(mean(nmme.old.cor[lead.pic,nget.clim,])) > 1.e-10) stop('NMME climatology has nonzero mean')

	yrange = range(nmme.old.cor[lead.pic,nget,])

	plot  (spear.list$targ.nmme[lead.pic,nget,],nmme.arx.cor[lead.pic,nget,],cex=0.5,xlab='',ylab='anomaly (degrees C)',col='red',ylim=yrange,pch=19)
	abline(h=0,col='grey',lty='dashed',lwd=2)
	points(spear.list$targ.nmme[lead.pic,nget,],best.arx.cor[lead.pic,nget,],col='black',cex=1,pch=17)
	# points(spear.list$targ.best[nbest],best.old.cor[nbest],col='black',cex=1,pch=17)

	y           = as.numeric(best.arx.cor[lead.pic,nget.clim,]) 
	x           = as.numeric(year(spear.list$targ.nmme[lead.pic,nget.clim,]))/10
	best.sum    = summary(lm(y~x))
	y           = as.numeric(nmme.arx.cor[lead.pic,nget.clim,])
	x           = as.numeric(year(spear.list$targ.nmme[lead.pic,nget.clim,]))/10
	nmme.sum    = summary(lm(y~x))
	best.trend  = round(best.sum$coef[2,1:2],3)
	nmme.trend  = round(nmme.sum$coef[2,1:2],3)
	legend.best = bquote('Berkeley ('*.(best.trend[1])%+-%.(best.trend[2])*' C/decade)')
	legend.nmme = bquote('SPEAR   ('*.(nmme.trend[1])%+-%.(nmme.trend[2])*' C/decade)')

	ftitle.top = paste("Lead-",lead.pic-0.5," Anomalies Based on ARX Bias Correction",sep='')
	ftitle.bot = paste('Target is ',month.abb[mpic],'; ',region,sep='')
	title(ftitle.top,line=2.0)
	title(ftitle.bot,line=0.5)
	legend('topleft',legend=c(legend.nmme,legend.best),pch=c(19,17),col=c('red','black'),cex=1.2)
	if (lplotfile) dev.off()
}

###########################################
############# PLOT ERRORS
###########################################
dim(nmme.old.cor)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(nmme.arx.cor)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(best.arx.cor)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(n.best.nmmeformat)    = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(spear.list$targ.nmme) = c(nlead.nmme,nstart.nmme,emax.nmme)

dim(err.nmme.old)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(err.nmme.arx)         = c(nlead.nmme,nstart.nmme,emax.nmme)

yrange = c(-0.5,0.5)
for (mpic in c(1,4,7,10)) {
	fout = paste(dir.figs,'SPEAR.',model.forcing,'.lead',lead.pic,'error.',month.abb[mpic],sep='')
	if (lplotfile) pdf.eps(fout,'pdf',height=5)
	matrix.layout = matrix(c(1,1,2,3),byrow=TRUE,ncol=2)
	layout(matrix.layout,heights=c(1,3))
	par(mar=c(0.1,0.1,0.1,0.1))
	par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
	plot(1,1,type='n',axes=FALSE,xlab='',ylab='')
	ftitle.top = paste("Lead-",lead.pic-0.5," Errors of SPEAR",sep='')
	ftitle.bot = paste('Target is ',month.abb[mpic],'; ARX Fit: H= ',nharm,'; P=',p.order,'; ',region,sep='')
	text(1,1,paste(ftitle.top,'\n',ftitle.bot),cex=2)
	nget = which( month(spear.list$targ.nmme[1,,1]) == mpic )
	# yrange = range(err.nmme.old[lead.pic,nget,],err.nmme.arx[lead.pic,nget,],na.rm=TRUE)
	par(mar=c(8,5,3,1))
	plot  (spear.list$targ.nmme[lead.pic,nget,],err.nmme.old[lead.pic,nget,],cex=0.5,xlab='',ylab='error',col='blue',ylim=yrange,pch=19)
	abline(h=0,col='grey',lty='dashed',lwd=2)
	title('traditional bias correction',line=0.5)
	plot  (spear.list$targ.nmme[lead.pic,nget,],err.nmme.arx[lead.pic,nget,],cex=0.5,xlab='',ylab='error',col='red' ,ylim=yrange,pch=19)
	abline(h=0,col='grey',lty='dashed',lwd=2)
	title('ARX bias correction',line=0.5)
	# title(ftitle.top,line=2.0)
	# title(ftitle.bot,line=0.5)
	#legend('topleft',legend=c('traditional bias correction','ARX bias correction'),pch=c(1,4),col=c('blue','red'),cex=1.2)
	if (lplotfile) dev.off()
}

###########################################
############# PLOT MSE
###########################################
############
dim(err.nmme.old)              = c(nlead.nmme*nstart.nmme,emax.nmme)
dim(err.nmme.arx)              = c(nlead.nmme*nstart.nmme,emax.nmme)

err.nmme.old.emean             = rowMeans(err.nmme.old)
err.nmme.arx.emean             = rowMeans(err.nmme.arx)
dim(err.nmme.old.emean)        = c(nlead.nmme,nstart.nmme)
dim(err.nmme.arx.emean)        = c(nlead.nmme,nstart.nmme)
dim(best.arx.cor)              = c(nlead.nmme,nstart.nmme,emax.nmme)

mse.best.arx                   = rowMeans( best.arx.cor[,,1]^2,na.rm=TRUE)
mse.nmme.old                   = rowMeans(err.nmme.old.emean^2,na.rm=TRUE)
mse.nmme.arx                   = rowMeans(err.nmme.arx.emean^2,na.rm=TRUE)


fout = paste(dir.figs,'SPEAR.',model.forcing,'.mse',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
yrange = sqrt(range(mse.best.arx[-(1:p.order)],mse.nmme.old ,mse.nmme.arx[-(1:p.order)]))
xrange = c(1,nlead.nmme)-0.5
par(mfcol=c(1,1),mar=c(8,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='lead (months)',ylab='Square Root of MSE (degrees C)',xaxt='n')
axis(1,1:nlead.nmme-0.5,labels=1:nlead.nmme-0.5)
lead.clip = (1+p.order):nlead.nmme
lines(1:nlead.nmme-0.5,sqrt(mse.nmme.old),col='black',lwd=2)
lines(lead.clip-0.5,sqrt(mse.best.arx[lead.clip]),col='blue',lwd=2)
lines(lead.clip-0.5,sqrt(mse.nmme.arx[lead.clip]),col='red',lwd=2)
points(1:nlead.nmme-0.5,sqrt(mse.nmme.old),col='black',pch=19,cex=0.6)
points(lead.clip-0.5,sqrt(mse.best.arx[lead.clip]),col='blue',pch=19,cex=0.6)
points(lead.clip-0.5,sqrt(mse.nmme.arx[lead.clip]),col='red',pch=19,cex=0.6)
# lines(mse.ideal,col='blueviolet',lwd=2)
legend('topleft',legend=c(
  'Traditional Bias Correction',
  'ARX Trained on Berkeley',
  'ARX SPEAR Correction: Stationary'),
  # 'Theoretical Lower Bound: Stationary'),
  col=c('black','blue','red'),lwd=2,bg='white',cex=1.1,
  lty=c('solid','solid','solid'))
ftitle.top = paste("Square Root of Mean Square Errors")
ftitle.bot = paste('ARX Fit: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()

##################################################
############# PLOT MSE: TRAINING AND VALIDATION
##################################################
mse.nmme.old.train = as.numeric(rep(NA,nlead.nmme))
mse.nmme.old.valid = as.numeric(rep(NA,nlead.nmme))
mse.nmme.arx.train = as.numeric(rep(NA,nlead.nmme))
mse.nmme.arx.valid = as.numeric(rep(NA,nlead.nmme))
mse.best.arx.train = as.numeric(rep(NA,nlead.nmme))
mse.best.arx.valid = as.numeric(rep(NA,nlead.nmme))
for (lead in 1:nlead.nmme) {
	npic = spear.list$targ.nmme[lead,,1] >= train.str & spear.list$targ.nmme[lead,,1] <= train.end
	mse.nmme.old.train[lead] = mean(err.nmme.old.emean[lead,npic]^2,na.rm=TRUE)
	mse.nmme.arx.train[lead] = mean(err.nmme.arx.emean[lead,npic]^2,na.rm=TRUE)
	mse.best.arx.train[lead] = mean(best.arx.cor[lead,npic,1]^2,na.rm=TRUE)
	npic = spear.list$targ.nmme[lead,,1] >= train.end
	mse.nmme.old.valid[lead] = mean(err.nmme.old.emean[lead,npic]^2,na.rm=TRUE)
	mse.nmme.arx.valid[lead] = mean(err.nmme.arx.emean[lead,npic]^2,na.rm=TRUE)
	mse.best.arx.valid[lead] = mean(best.arx.cor[lead,npic,1]^2,na.rm=TRUE)
}

fout = paste(dir.figs,'SPEAR.',model.forcing,'.mse.trainvalid',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
yrange = sqrt(range(mse.nmme.arx.train[-(1:p.order)],mse.nmme.old.train,mse.nmme.arx.valid[-(1:p.order)],mse.nmme.old.valid))
xrange = c(1,nlead.nmme)-0.5
par(mfcol=c(1,1),mar=c(8,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='lead (months)',ylab='Square Root of MSE (degrees C)',xaxt='n')
axis(1,1:nlead.nmme-0.5,labels=1:nlead.nmme-0.5)
lead.clip = (1+p.order):nlead.nmme
lines(1:nlead.nmme-0.5,sqrt(mse.nmme.old.train),col='black',lwd=2)
lines(lead.clip-0.5,sqrt(mse.best.arx.train[lead.clip]),col='blue',lwd=2)
lines(lead.clip-0.5,sqrt(mse.nmme.arx.train[lead.clip]),col='red',lwd=2)
points(1:nlead.nmme-0.5,sqrt(mse.nmme.old.train),col='black',pch=19,cex=0.6)
points(lead.clip-0.5,sqrt(mse.best.arx.train[lead.clip]),col='blue',pch=19,cex=0.6)
points(lead.clip-0.5,sqrt(mse.nmme.arx.train[lead.clip]),col='red',pch=19,cex=0.6)

lines(1:nlead.nmme-0.5,sqrt(mse.nmme.old.valid),col='black',lwd=2,lty='dashed')
lines(lead.clip-0.5,sqrt(mse.best.arx.valid[lead.clip]),col='blue',lwd=2,lty='dashed')
lines(lead.clip-0.5,sqrt(mse.nmme.arx.valid[lead.clip]),col='red',lwd=2,lty='dashed')
points(1:nlead.nmme-0.5,sqrt(mse.nmme.old.valid),col='black',pch=15,cex=0.6)
points(lead.clip-0.5,sqrt(mse.best.arx.valid[lead.clip]),col='blue',pch=15,cex=0.6)
points(lead.clip-0.5,sqrt(mse.nmme.arx.valid[lead.clip]),col='red',pch=15,cex=0.6)

legend('topleft',legend=c(
  'Traditional Bias Correction',
  'ARX Trained on Berkeley',
  'ARX SPEAR Correction: Stationary'),
  # 'Theoretical Lower Bound: Stationary'),
  col=c('black','blue','red'),lwd=2,bg='white',cex=1.1,
  lty=c('solid','solid','solid'))
ftitle.top = paste("Square Root of Mean Square Errors")
ftitle.bot = paste('ARX Fit: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()


###############################################################################
############# CHECKERBOARD PLOT OF MSE AVERAGED OVER YEARS
###############################################################################

### COMPUTE MSE OF ENSEMBLE MEAN  
dim(err.nmme.old.emean)        = c(nlead.nmme,nstart.nmme)
dim(err.nmme.arx.emean)        = c(nlead.nmme,nstart.nmme)
mse.nmme.old.emean.checker     = array(NA,dim=c(nlead.nmme,12))
mse.nmme.arx.emean.checker     = array(NA,dim=c(nlead.nmme,12))
for (m in 1:12) for (lead in 1:nlead.nmme) {
	mpic = month(spear.list$targ.nmme[lead,,1]) == m
	mse.nmme.old.emean.checker[lead,m] = mean(err.nmme.old.emean[lead,mpic]^2,na.rm=TRUE)
	mse.nmme.arx.emean.checker[lead,m] = mean(err.nmme.arx.emean[lead,mpic]^2,na.rm=TRUE)
}
mse.nmme.arx.emean.checker[1:p.order,] = NA


fout = paste(dir.figs,'SPEAR.mse.checker',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=5,width=8.5)
par(mfrow=c(1,3),mar=c(12,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
breaks      = seq(from=0,to=0.06,by=0.01)
ncol        = length(breaks) - 1
col.checker = rev(heat.colors(ncol))
col.checker = colorRampPalette(c('orange','red','darkred'),space='rgb')(ncol)
col.checker = rev(gray.colors(ncol))
col.checker[1] = 'white'
image(1:nlead.nmme-0.5,1:12,mse.nmme.old.emean.checker,breaks=breaks,col=col.checker,xlab='lead (months)',ylab='target month')
addMapLegend(cutVector=breaks,colourVector=col.checker,legendLabels='all',legendMar=2)
title(main='MSE for Traditional Climatology')

image(1:nlead.nmme-0.5,1:12,mse.nmme.arx.emean.checker,breaks=breaks,col=col.checker,xlab='lead (months)',ylab='target month')
addMapLegend(cutVector=breaks,colourVector=col.checker,legendLabels='all',legendMar=2)
title(main='MSE for ARX Climatology')
x.grid = rep(1:p.order,each=12)
y.grid = rep(1:12,p.order)
points(x.grid,y.grid,pch=19,cex=0.5)

arx.minus.old =  mse.nmme.arx.emean.checker - mse.nmme.old.emean.checker
breaks        = seq(from=-0.03,to=0.02,by=0.01)
n.neg         = sum(breaks < 0)
n.pos         = sum(breaks > 0)
n.tot         = n.neg + n.pos + 1
if (length(breaks) != n.tot) print('breaks in checker plot do not make sense')
col.neg       = colorRampPalette(c('blue','lightblue'),space='rgb')(n.neg)
col.pos       = colorRampPalette(c('orange','red','darkred'),space='rgb')(n.pos)
col.checker   = c(col.neg,col.pos)

image(1:nlead.nmme-0.5,1:12,arx.minus.old,breaks=breaks,col=col.checker,xlab='lead (months)',ylab='target month')
addMapLegend(cutVector=breaks,colourVector=col.checker,legendLabels='all',legendMar=2)
title(main='ARX-Traditional MSEs')
points(x.grid,y.grid,pch=19,cex=0.5)
if (lplotfile) dev.off()

###############################################################################
############# PLOT RAW, ARX PREDICTIONS, AND BIAS CORRECTION FOR A FEW CASES
###############################################################################         
dim(spear.list$targ.nmme)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(spear.list$pc.nmme)           = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(simlist.best.nmmeformat.stat) = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(simlist.nmme.nmmeformat.stat) = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(nmme.emean)                   = c(nlead.nmme,nstart.nmme)
dim(nmme.stdv)                    = c(nlead.nmme,nstart.nmme)

earliest.plot                     = as.Date('15-01-1991',format='%d-%m-%Y')
latest.plot                       = as.Date('15-01-1998',format='%d-%m-%Y')
nst.best                          = which(spear.list$targ.best == earliest.plot)
nnd.best                          = which(spear.list$targ.best == latest.plot)
if (length(nst.best) != 1 | length(nnd.best) != 1) stop('cannot find unique start/end times in BEST')
nst.best.old = nst.best

for (nm in c(1,4,7,10)) {

	fout = paste(dir.figs,'SPEAR.RawTS.ARXpred.',month.abb[nm],sep='')
	
	mpic.nmme = which( month(spear.list$targ.nmme[1,,1]) == nm & 
	                    year(spear.list$targ.nmme[1,,1]) >= year(earliest.plot) &
	                    year(spear.list$targ.nmme[1,,1]) <= year(latest.plot) )
	
	yrange = range(simlist.nmme.nmmeformat.stat[,mpic.nmme,],
	               simlist.best.nmmeformat.stat[,mpic.nmme,],
	               nmme.emean[,mpic.nmme]+nmme.stdv[,mpic.nmme],
	               nmme.emean[,mpic.nmme]-nmme.stdv[,mpic.nmme],
	               spear.list$pc.best[nst.best:nnd.best])                    
	yrange[2] = yrange[2] + 1

	if (lplotfile) pdf.eps(fout,'pdf')
	par(mfcol=c(1,1),mar=c(5,5,3,1))
	par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
	plot(c(earliest.plot,latest.plot),c(1,1),type='n',ylim=yrange,xlab='',ylab='Average 2m-Temperature')
	
	     
	for (m in mpic.nmme) {
		arx.nmme.emean = rowMeans(simlist.nmme.nmmeformat.stat[,m,])
		arx.best.emean = rowMeans(simlist.best.nmmeformat.stat[,m,])
		x.poly = c(spear.list$targ.nmme[,m,1],rev(spear.list$targ.nmme[,m,1]))
		y.poly = c(nmme.emean[,m]+nmme.stdv[,m],rev(nmme.emean[,m]-nmme.stdv[,m]))
		
		polygon(x.poly,y.poly,col='green1',border='green1')
		lines(spear.list$targ.nmme[,m,1],arx.nmme.emean,col='darkgreen',lwd=2)
		lines(spear.list$targ.nmme[,m,1],arx.best.emean,col='grey',lwd=2)
		# lines(spear.list$targ.nmme[,m,1],nmme.emean[,m] - arx.nmme.emean + arx.best.emean,col='red',lwd=2)
	}

	nst.best = nst.best.old + nm - 1
	# lines(spear.list$targ.best[nst.best:nnd.best],spear.list$pc.best[nst.best:nnd.best],col='black',lwd=2)
	points(spear.list$targ.best[nst.best:nnd.best],spear.list$pc.best[nst.best:nnd.best],col='black',pch=19,cex=0.7)
	
	legend('topleft',legend=
	   # c('raw NMME','ARX-NMME','raw Berkeley','ARX-Berkeley','SPEAR ARX Correction'),
	   # col=c('green1','darkgreen','black','grey','red'),lwd=2,bg='white',cex=1.0,ncol=3)
	   c('SPEAR','ARX-SPEAR','Berkeley','ARX-Berkeley'),
	   col=c('green1','darkgreen','black','grey'),pch=c(NA,NA,19,NA),lwd=c(2,2,NA,2),bg='white',cex=1.0,ncol=3)
	
	ftitle.top = paste("SPEAR-NMME and Berkeley Data Sets with ARX Predictions")
	ftitle.bot = paste(month.abb[nm],' start; ARX Fit: H= ',nharm,'; P=',p.order,'; ',region,sep='')
	title(ftitle.top,line=2.0)
	title(ftitle.bot,line=0.5)
	if (lplotfile) dev.off()

}

###########################################
############# COMPARE LEAD-DEPENDENT TRENDS
###########################################
best.nmmeformat                   = spear.list$pc.best[n.best.nmmeformat]
best.old.clim.nmme                = spear.list$pc.best[n.best.nmmeformat] - best.old.cor[n.best.nmmeformat]
best.arx.clim.nmme                = simlist.best.nmmeformat.stat

dim(spear.list$targ.nmme)         = c(nlead.nmme,nstart.nmme,emax.nmme)
dim(nmme.old.cor)                 = c(nlead.nmme*nstart.nmme,emax.nmme)
dim(nmme.arx.cor)                 = c(nlead.nmme*nstart.nmme,emax.nmme)
dim(best.old.clim.nmme)           = c(nlead.nmme*nstart.nmme,emax.nmme)
dim(best.arx.clim.nmme)           = c(nlead.nmme*nstart.nmme,emax.nmme)
dim(best.nmmeformat)              = c(nlead.nmme*nstart.nmme,emax.nmme)

nmme.old.cor.emean                = rowMeans(nmme.old.cor)
nmme.arx.cor.emean                = rowMeans(nmme.arx.cor)
best.old.clim.nmme                = rowMeans(best.old.clim.nmme)
best.arx.clim.nmme                = rowMeans(best.arx.clim.nmme)
best.nmmeformat.emean             = rowMeans(best.nmmeformat)

dim(nmme.old.cor.emean)           = c(nlead.nmme,nstart.nmme)
dim(nmme.arx.cor.emean)           = c(nlead.nmme,nstart.nmme)
dim(best.old.clim.nmme)           = c(nlead.nmme,nstart.nmme)
dim(best.nmmeformat.emean)        = c(nlead.nmme,nstart.nmme)

nmme.total.arx.cor                = nmme.arx.cor.emean + best.arx.clim.nmme  
nmme.total.old.cor                = nmme.old.cor.emean + best.old.clim.nmme   ### best.old.clim.nmme is 12 x 401, repeating a 12 x 12 matrix

### INCLUDE ONLY DATA IN TRAINING PERIOD
lgood = 


ctype.all                 = c('old','arx','best')
trend.coef                = array(NA,dim=c(nlead.nmme,length(ctype.all)),dimnames = list(paste('lead',1:nlead.nmme,sep=''),ctype.all))
trend.stdv                = array(NA,dim=c(nlead.nmme,length(ctype.all)),dimnames = list(paste('lead',1:nlead.nmme,sep=''),ctype.all))
for (ctype in ctype.all) for (lead in 1:nlead.nmme) {
	x             = as.numeric(spear.list$targ.nmme[lead,,1])
	x             = x/365.64/10
	if (nharm > 0) for (nh in 1:nharm) x = cbind(x,
	     cos(2*pi*month(spear.list$targ.nmme[lead,,1])*nh/12),
	     sin(2*pi*month(spear.list$targ.nmme[lead,,1])*nh/12))
	if (ctype == 'old') {
		y  = nmme.total.old.cor[lead,]	
	} else if (ctype == 'arx') {
		y  = nmme.total.arx.cor[lead,]			
	} else if (ctype == 'best') {
		y  = best.nmmeformat.emean[lead,]
	}
	
	lgood = spear.list$targ.nmme[lead,,1] >= train.str & spear.list$targ.nmme[lead,,1] <= train.end
	y[!lgood]              = NA
	trend.summary          = summary(lm(y~x))
	trend.coef[lead,ctype] = trend.summary$coef[2,'Estimate']
	trend.stdv[lead,ctype] = trend.summary$coef[2,'Std. Error'] 
}


fout = paste(dir.figs,'SPEAR.',model.forcing,'.trend',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=5,width=8.5)
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
yrange = range(trend.coef + trend.stdv, trend.coef - trend.stdv)
# yrange[2] = yrange[2] + 0.05
xrange = c(1,nlead.nmme)-0.5
xoff   = 0.1
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='lead (months)',ylab='trend (C/decade)',xaxt='n')
axis(1,1:nlead.nmme     -0.5,labels=1:nlead.nmme     -0.5)
arrows(1:nlead.nmme-xoff-0.5,(trend.coef + trend.stdv)[,'old' ],1:nlead.nmme-xoff-0.5,(trend.coef - trend.stdv)[,'old' ],code=3,angle=90,length=0.02,lwd=2,col='black')
arrows(1:nlead.nmme     -0.5,(trend.coef + trend.stdv)[,'arx' ],1:nlead.nmme     -0.5,(trend.coef - trend.stdv)[,'arx' ],code=3,angle=90,length=0.02,lwd=2,col='red')
arrows(1:nlead.nmme+xoff-0.5,(trend.coef + trend.stdv)[,'best'],1:nlead.nmme+xoff-0.5,(trend.coef - trend.stdv)[,'best'],code=3,angle=90,length=0.02,lwd=2,col='blue')
legend('right',legend=c('SPEAR: traditional correction','SPEAR: ARX correction','Observation: Berkeley'),col=c('black','red','blue'),lwd=2,bg='white',cex=1.3)
ftitle.top = paste("Trends in the Corrected SPEAR Forecasts and Observations",sep='')
ftitle.bot = paste('ARX Fit: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()


###############################################################################
###############################################################################
###############################################################################
############# KALMAN FILTER
###############################################################################
###############################################################################
###############################################################################
beta.kf       = beta.nmme.stat
suffix        = paste('Lag',1:nlag.forcing-1,sep='')

### BERKELEY COEFFICIENTS TRAINED ON ALL THE DATA
# beta.kf[paste('GHG',suffix,sep='')] = beta.best.stat[paste('GHG',suffix,sep='')]
# beta.kf[paste('AER',suffix,sep='')] = beta.best.stat[paste('AER',suffix,sep='')]
# beta.kf[paste('NAT',suffix,sep='')] = beta.best.stat[paste('NAT',suffix,sep='')]

### DO NOT USE 'SHORT'-- CONFIDENCE INTERVALS CROSS ZERO
# beta.kf[paste('GHG',suffix,sep='')] = beta.best.short[paste('GHG',suffix,sep='')]
# beta.kf[paste('AER',suffix,sep='')] = beta.best.short[paste('AER',suffix,sep='')]
# beta.kf[paste('NAT',suffix,sep='')] = beta.best.short[paste('NAT',suffix,sep='')]

#### LENS RESPONSE COEFFICIENTS
# beta.kf[paste('GHG',suffix,sep='')] = beta.hist.stat[paste('GHG',suffix,sep='')]
# beta.kf[paste('AER',suffix,sep='')] = beta.hist.stat[paste('AER',suffix,sep='')]
# beta.kf[paste('NAT',suffix,sep='')] = beta.hist.stat[paste('NAT',suffix,sep='')]



dim(spear.list$pc.hist)      = c(ntot.hist,emax.hist)

nvar.kf       = nvar.nmme.stat
rvar.kf       = nvar.kf*15
# rvar.kf       = nvar.kf*0.5
# rvar.kf       = nvar.kf*10
init.date.kf  = as.Date('15-01-1990',format='%d-%m-%Y')
targ.date.kf  = as.Date('15-01-2024',format='%d-%m-%Y')
ntot.kf       = month(targ.date.kf) - month(init.date.kf) + 12 * (year(targ.date.kf) - year(init.date.kf))

nst.hist      = which(spear.list$targ.hist[,1] == init.date.kf)
if (length(nst.hist) != 1) stop("cannot find unique LENS date")

nst.best      = which(spear.list$targ.best == init.date.kf)
if (length(nst.best) != 1) stop("cannot find unique BEST date")

nst.forc      = which(spear.list$targ.forcing == init.date.kf)
if (length(nst.forc) != 1) stop("cannot find unique Forcing date")


init.kf       = spear.list$pc.hist  [nst.hist+(1:p.order)-1,1]
obs.kf        = spear.list$pc.best  [nst.best+1:ntot.kf-1]
date.kf       = spear.list$targ.best[nst.best+1:ntot.kf-1]
forcing.kf    = spear.list$pc.forcing[nst.forc+1:ntot.kf-1,]
phase         = (month(date.kf) - 1) %% 12 + 1


mean.kf.a    = as.numeric(rep(NA,ntot.kf))
covm.kf.a    = as.numeric(rep(NA,ntot.kf))
mean.kf.a[1:p.order] = obs.kf[1:p.order]

#### THIS IS A NOISE-FREE SOLUTION.  THE INTEGRATION SHOULD REPRODUCE THIS
init.nmme.pc           = spear.list$pc.nmme[1,,1]
init.nmme.targ         = spear.list$targ.nmme[1,,1]
# obs.kf               = simlist.nmme.stat[nst.hist + 1:ntot.kf-1]
# mean.kf.a[1:p.order] = simlist.nmme.stat[nst.hist + 1:p.order-1]
# covm.kf.a[1:p.order] = nvar.kf * 2

for (n in (1+p.order):ntot.kf) {
	mean.kf.b = beta.kf['intercept']
	if (p.order  > 0) for (lag in 1:p.order ) mean.kf.b = mean.kf.b + mean.kf.a[n-lag] * beta.kf[lag]
	if (nharm    > 0) for (nh  in 1:nharm   ) mean.kf.b = mean.kf.b + cos(2*pi*phase[n]*nh/12) * beta.kf[p.order + 2*nh - 1]
	if (nharm    > 0) for (nh  in 1:nharm   ) mean.kf.b = mean.kf.b + sin(2*pi*phase[n]*nh/12) * beta.kf[p.order + 2*nh    ]
	if (nforcing > 0) for (nf  in 1:nforcing) mean.kf.b = mean.kf.b + forcing.kf[n,nf] * beta.kf[p.order + 2*nharm + 1 + nf]
	
	covm.kf.b    = nvar.kf
	mean.kf.a[n] = mean.kf.b + covm.kf.b/(covm.kf.b + rvar.kf) * ( obs.kf[n] - mean.kf.b)
	covm.kf.a[n] = covm.kf.b * rvar.kf / (covm.kf.b + rvar.kf) 
}

fout = paste(dir.figs,'SPEAR.',model.forcing,'.kf',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=5,width=8.5)
rvar.say = round(rvar.kf/nvar.kf)
kf.say   = paste('Mean Analysis (r= ',rvar.say,')',sep='')
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
yrange = range(mean.kf.a,obs.kf,na.rm=TRUE)
yrange = range(mean.kf.a,na.rm=TRUE)
xrange = c(init.date.kf,targ.date.kf)
plot(xrange,yrange,type='n',xlim=xrange,ylim=yrange,xlab='',ylab='')

lines(date.kf,mean.kf.a,lwd=2,col='red')
lines(init.nmme.targ,init.nmme.pc,col='black')
legend('topleft',legend=c('SPEAR lead 0.5',kf.say),col=c('black','red'),lwd=2,bg='white',cex=1.1)

# lines(date.kf,obs.kf,lwd=2,col='darkgreen')
# legend('topleft',legend=c('Berkeley','ARX (noise free)','NMME lead 0.5',kf.say),col=c('darkgreen','blue','black','red'),lwd=2,bg='white',cex=1.1)
ftitle.top = paste("Data Assimilation with Berkeley Data and ARX trained on SPEAR",sep='')
ftitle.bot = paste('ARX Fit: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()

nst.kf   = which(date.kf == init.nmme.targ[1])
nst.nmme = 1
nnd.kf   = length(date.kf)
nnd.nmme = which(init.nmme.targ == date.kf[length(date.kf)])
print(all(date.kf[nst.kf:nnd.kf] == init.nmme.targ[nst.nmme:nnd.nmme]))
print(paste('rvar/nvar=',rvar.kf/nvar.kf,'MSE KF=',mean((init.nmme.pc[nst.nmme:nnd.nmme]-mean.kf.a[nst.kf:nnd.kf])^2)))

### COMPARE TRENDS
x.kf = as.numeric(date.kf)/364.24/10
if (nharm > 0) for (nh in 1:nharm) x.kf = cbind(x.kf,cos(2*pi*phase*nh/12),sin(2*pi*phase*nh/12))
trend.obs = summary(lm(obs.kf[nst.kf:nnd.kf] ~ x.kf[nst.kf:nnd.kf,]))
trend.kf  = summary(lm(mean.kf.a[nst.kf:nnd.kf] ~ x.kf[nst.kf:nnd.kf,]))
print('obs')
print(trend.obs[['coefficients']][2,1:2])
print('KF analysis')
print(trend.kf[['coefficients']][2,1:2])


##################################################################
########## PLOT LEAD-1 ERRORS OF KALMAN FILTER: TIME SERIES
##################################################################
err.kf              = mean.kf.a - obs.kf
trend.init.error.kf = lm.TrendPlusCycle(err.kf,date.kf,train.str,train.end)
trend.init.err.kf   = round(c(trend.init.error.kf$trend, 2* trend.init.error.kf$sterr),3)
sigma.init.err.kf   = signif(trend.init.error.kf$sigma,2)

##################################################################
########## PLOT LEAD-1 ERRORS: TIME SERIES
##################################################################
fout = paste(dir.figs,'SPEAR.',model.forcing,'.lead1error',sep='')
if (lplotfile) pdf.eps(fout,'pdf')
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
iyst.plot  = year(date.kf[1])
iynd.plot  = year(date.kf[length(err.init.emean.targ)])
yrange     = range(err.kf,na.rm=TRUE)
xrange     = c(iyst.plot,iynd.plot+1)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='',ylab='error')
for (m in 1:12) {
	nget = which( month(date.kf) == m)
	lines(year(date.kf)[nget],err.kf[nget],col=m,lwd=2)
	lines(year(date.kf)[nget],trend.init.error.kf$y.pred[nget],col=m)
	# text (iynd.plot-1,err.init.emean[nget[length(nget)-1]],month.abb[m],pos=4,col=m)
	text (iynd.plot,trend.init.error.kf$y.pred[nget[length(nget)]],month.abb[m],pos=4,col=m)
}
model.stuff = bquote(.(trend.init.err.kf[1])%+-%.(trend.init.err.kf[2])*' C/decade; StErr='~.(sigma.init.err.kf)*'C')
text(iynd.plot,yrange[1],model.stuff,pos=2,cex=1.1)
ftitle.top = paste("Initial Error from Data Assimilation (lead = 0.5; ensemble mean)")
ftitle.bot = paste('ARX parameters: H= ',nharm,'; P=',p.order,'; ',region,sep='')
title(ftitle.top,line=2.0)
title(ftitle.bot,line=0.5)
if (lplotfile) dev.off()

