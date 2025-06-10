traditional.correction = function(f.ts,f.date,o.ts,o.date,nlead,nstart,nens,train.str=NA,train.end=NA,ldetrend=TRUE) {
## COMPUTE ARX-CORRECTION BASED ON FORECAST AND OBSERVATION TIME SERIES
## INPUT:
##	F.TS   [NLEAD,NSTART,NENS]: FORECAST TIME SERIES
## 	F.DATE [NLEAD,NSTART,NENS]: VERIFICATION DATES OF THE FORECAST
##	O.TS   [NOBS]:				OBSERVATION TIME SERIES
##	O.DATE [NOBS]:				VERIFIATION DATES OF OBSERVATIONS
##	NLEAD:						NUMBER OF LEAD TIMES IN FORECAST
##	NSTART:						NUMBER OF START TIMES IN FORECAST
##	NENS:						NUMBER OF ENSEMBLE MEMBERS PER LEAD,START
##  TRAIN.STR					BEGINNING START DATE OF TRAINING PERIOD (IGNORED IF SET TO NA)
##  TRAIN.END 					ENDING START DATE OF TRAINING PERIOD (INGORED IF SET TO NA)
##
## DEPENDENCIES:
##	LUBRIDATE 
##	N.TO.MONTHLY


if (length(f.ts) != nlead*nstart*nens) stop('f.ts incorrectly dimensioned')
if (length(o.ts) != length(o.date))    stop('o.ts and o.date have inconsistent dimensions')

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
############ CORRECTION WITH LINEAR TREND
################################################
nbeta         = ifelse(ldetrend,2,1)
beta.for.save = array(NA,dim=c(nbeta,nlead,12))
dim(f.ts)     = c(nlead,nstart*nens)
dim(f.date)   = c(nlead,nstart*nens)
forecast.anom = array(NA,dim=c(nlead,nstart*nens))
forecast.clim = array(NA,dim=c(nlead,nstart*nens))
for (mon in 1:12) for (lead in 1:nlead) {
  npic       = which(month(f.date[lead,]) == mon)
  npic.train = which(month(f.date[lead,]) == mon & l.f.train)
  if (ldetrend) {
    y          =            f.ts  [lead,npic.train]
    x1         = as.numeric(f.date[lead,npic.train])
    beta.for   = lm(y~x1)$coefficients
    forecast.clim[lead,npic] = beta.for[1] + beta.for[2]*as.numeric(f.date[lead,npic])
    forecast.anom[lead,npic] = f.ts[lead,npic] - forecast.clim[lead,npic]
  } else {
    beta.for                 = mean(f.ts[lead,npic.train])
    forecast.clim[lead,npic] = beta.for
    forecast.anom[lead,npic] = f.ts[lead,npic] - forecast.clim[lead,npic]
  }
  beta.for.save[,lead,mon] = beta.for
}

beta.obs.save        = array(NA,dim=c(nbeta,12))
obs.anom             = as.numeric(rep(NA,length(o.date)))
obs.clim             = as.numeric(rep(NA,length(o.date)))
for (mon in 1:12) {
  npic       = which(month(o.date) == mon)
  npic.train = which(month(o.date) == mon & l.o.train)
  if (ldetrend) {
    y          =            o.ts  [npic.train]
    x1         = as.numeric(o.date[npic.train])
    beta.obs   = lm(y~x1)$coefficients
    obs.clim[npic] = beta.obs[1] + beta.obs[2] * as.numeric(o.date[npic])   
  } else {
    beta.obs   = mean(o.ts[npic.train])
    obs.clim[npic] = beta.obs
  }
  obs.anom[npic] = o.ts[npic] - obs.clim[npic]
  beta.obs.save[,mon] = beta.obs
}


################################################
################################################
################################################
### PUT OBSERVATIONS IN FORECAST FORMAT
n.obs.forecast.format   = n.to.monthly(f.date,o.date)
if (any(as.numeric(f.date) != as.numeric(o.date[n.obs.forecast.format]),na.rm=TRUE)) stop('obs and forecast dates misalined')
obs.anom.forecast.format = obs.anom[n.obs.forecast.format]

### COMPUTE ERROR OF TRADITIONAL CORRECTION
err.corrected           = forecast.anom - obs.anom.forecast.format

### COMPUTE CORRECTED FORECAST
f.corrected             = forecast.anom + obs.clim[n.obs.forecast.format]

### COMPUTE UNCORRECTED ERROR
err.uncorrected         = as.numeric(f.ts) - o.ts[n.obs.forecast.format]


list(forecast.anom            = forecast.anom, 
     obs.anom                 = obs.anom, 
     train.str                = train.str, 
     train.end                = train.end,
     obs.anom.forecast.format = obs.anom.forecast.format,
     err.corrected            = err.corrected,
     f.corrected              = f.corrected,
     err.uncorrected          = err.uncorrected,
     beta.for                 = beta.for.save,
     beta.obs                 = beta.obs.save)

}