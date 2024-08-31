arx.correction = function(f.ts,f.date,o.ts,o.date,nforcing,forcing.ts,forcing.date,forcing.names,
                          nlead,nstart,nens,order.pic,nharm,
                          train.str=NA, train.end=NA, mperiod=1) {
## COMPUTE ARX-CORRECTION BASED ON FORECAST AND OBSERVATION TIME SERIES
## INPUT:
##	F.TS   [NLEAD,NSTART,NENS]: FORECAST TIME SERIES
## 	F.DATE [NLEAD,NSTART,NENS]: VERIFICATION DATES OF THE FORECAST
##	O.TS   [NOBS]:				OBSERVATION TIME SERIES
##	O.DATE [NOBS]:				VERIFIATION DATES OF OBSERVATIONS
##	NFORCING:					NUMBER OF FORCINGS; SET TO ZERO TO IGNORE ALL FORCINGS
##	FORCING.TS[FTIME,NFORCING]:	FORCING TIME SERIES, EACH OF LENGTH FTIME; IGNORED IF FORCING=NA
##	FORCING.DATE[FTIME]:		VERIFICATION DATES OF FORCING; IGNORED IF FORCING=NA
##  FORCING.NAMES[NFORCING]:	NAME OF FORCINGS
##	NLEAD:						NUMBER OF LEAD TIMES IN FORECAST
##	NSTART:						NUMBER OF START TIMES IN FORECAST
##	NENS:						NUMBER OF ENSEMBLE MEMBERS PER LEAD,START
##	ORDER.PIC:					ORDER P OF THE ARX MODEL
##	NHARM:						NUMBER OF ANNUAL HARMONIC FORCINGS
##  TRAIN.STR					BEGINNING START DATE OF TRAINING PERIOD (IGNORED IF SET TO NA)
##  TRAIN.END 					ENDING START DATE OF TRAINING PERIOD (INGORED IF SET TO NA)
##	MPERIOD:					PERIOD OF CYCLOSTATIONARITY.  1 = STATIONARY, 12=CYCLOSTIONARY
##
## DEPENDENCIES:
##	LUBRIDATE 
##	N.TO.MONTHLY
##  TIMESERIES2ARX.CYCLO
##  SIMULATE.ARX.CYCLO

if (length(f.ts) != nlead*nstart*nens) stop('f.ts incorrectly dimensioned')
if (length(o.ts) != length(o.date))    stop('o.ts and o.date have inconsistent dimensions')

### FORCING IN FORECAST FORMAT: 
if (nforcing != 0) {
	ftime        = length(forcing.date)
	n.forcing.forecast.format = n.to.monthly(f.date,forcing.date)
	if (any(f.date != forcing.date[n.forcing.forecast.format])) stop('forcing and forecast dates misaligned')
	pc.forcing.forecast.format = forcing.ts[n.forcing.forecast.format,]
	
	n.forcing.obs.format      = n.to.monthly(o.date,forcing.date)
	if (any(o.date != forcing.date[n.forcing.obs.format])) stop('forcing and obs dates misaligned')
	pc.forcing.obs.format     = forcing.ts[n.forcing.obs.format,]
}

### OBSERVATION IN FORECAST FORMAT
n.obs.forecast.format   = n.to.monthly(f.date,o.date)
if (any(as.numeric(f.date) != as.numeric(o.date[n.obs.forecast.format]),na.rm=TRUE)) stop('obs and forecast dates misalined')

################################################
############ IDENTIFY TRAINING PERIOD
### INCLUDE ALL LEAD TIMES IF THE *LAST* LEAD TIME FALLS IN VERIFICATION PERIOD
################################################
dim(f.date)             = c(nlead,nstart*nens)
if (is.na(train.str)) train.str = -Inf
if (is.na(train.end)) train.end =  Inf
l.f.train.str = f.date[nlead,] >= train.str
l.f.train.end = f.date[nlead,] <= train.end
l.o.train.str = o.date         >= train.str
l.o.train.end = o.date         <= train.end

l.f.train = l.f.train.str & l.f.train.end
l.o.train = l.o.train.str & l.o.train.end

################################################
############ TRAIN ON FORECAST TIME SERIES
################################################
#### FIT FORECAST TIME SERIES TO ARX MODEL
dim(f.ts)               = c(nlead,nstart*nens)
dim(f.date)             = c(nlead,nstart*nens)
dim(pc.forcing.forecast.format) = c(nlead,nstart*nens,nforcing)
arx.forecast.mat        = timeseries2arx.cyclo(f.ts[,l.f.train],f.date[,l.f.train],pc.forcing.forecast.format[,l.f.train,],
                                               forcing.names,nlead,order.pic,nharm,mperiod)
y.forecast              = arx.forecast.mat$y.lhs
x.forecast              = cbind(arx.forecast.mat$x.ar, arx.forecast.mat$x.cyc, arx.forecast.mat$jvec, arx.forecast.mat$x.for)
lm.forecast             = lm(y.forecast ~ x.forecast - 1)
beta.forecast           = coef(lm.forecast)
nvar.forecast           = sum(residuals(lm.forecast)^2)/lm.forecast$df.residual
names(beta.forecast)    = colnames(x.forecast)
names(lm.forecast$coef) = colnames(x.forecast)

### NOISE-FREE SOLUTIONS OF ARX MODEL
simlist.forecast = simulate.arx.cyclo(f.ts[1:order.pic,],f.date,pc.forcing.forecast.format,forcing.names,
  beta.forecast,nlead,order.pic,nharm,mperiod)
  
## COMPUTE ARX CORRECTION FOR FORECAST
forecast.arx.anom = f.ts - simlist.forecast

################################################
############ TRAIN ON OBSERVATION TIME SERIES
################################################
nobs                     = length(o.date)
nobs.train               = length(o.date[l.o.train])
arx.obs.mat              = timeseries2arx.cyclo(o.ts[l.o.train],o.date[l.o.train],pc.forcing.obs.format[l.o.train],
                                                forcing.names,nobs.train,order.pic,nharm,mperiod)
y.obs                    = arx.obs.mat$y.lhs
x.obs                    = cbind(arx.obs.mat $x.ar, arx.obs.mat $x.cyc, arx.obs.mat $jvec, arx.obs.mat $x.for)
lm.obs                   = lm(y.obs ~ x.obs - 1)
beta.obs                 = coef(lm.obs)
nvar.obs                 = sum(residuals(lm.obs)^2)/lm.obs$df.residual
names(beta.obs)          = colnames(x.obs)
names(lm.obs$coef)       = colnames(x.obs)

#### NOISE-FREE SOLUTIONS OF ARX MODEL
o.forecast.format        = o.ts[n.obs.forecast.format]
dim(o.forecast.format)   = c(nlead,nstart*nens)
simlist.obs              = simulate.arx.cyclo(o.forecast.format[1:p.order,],f.date,pc.forcing.forecast.format,forcing.names,
                           beta.obs,nlead,order.pic,nharm,mperiod)
                           
## COMPUTE ARX CORRECTION FOR OBSERVATIONS
o.forecast.format.anom    = o.forecast.format - simlist.obs

## COMPUTE ERROR OF THE ARX CORRECTED FORECAST 
err.corrected             = forecast.arx.anom - o.forecast.format.anom

## COMPUTE ARX CORRECTED FORECAST
f.corrected               = forecast.arx.anom + simlist.obs

## COMPUTE ERROR OF UNCORRECTED FORECAST
err.uncorrected           = f.ts - o.forecast.format

list(pc.forcing.obs.format      = pc.forcing.obs.format,
     beta.forecast              = beta.forecast, 
     simlist.forecast           = simlist.forecast, 
     simlist.obs                = simlist.obs,
     forecast.arx.anom          = forecast.arx.anom,
     beta.obs                   = beta.obs,
     train.str                  = train.str,
     train.end                  = train.end,
     o.forecast.format.anom     = o.forecast.format.anom,
     n.obs.forecast.format      = n.obs.forecast.format,
     lm.obs                     = lm.obs,
     lm.forecast                = lm.forecast,
     nvar.obs                   = nvar.obs,
     nvar.forecast              = nvar.forecast,
     err.corrected              = err.corrected,
     f.corrected                = f.corrected,
     err.uncorrected            = err.uncorrected)

}