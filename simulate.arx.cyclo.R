simulate.arx.cyclo = function(initial,y.date,f.ts,fnames,beta.hat,nlead,order.pic,nharm,mperiod,noise.var=NA) {
### SIMULATE TIME SERIES USING FITTED CYCLOSTATIONARY MODEL
### INPUT:
###		INITIAL [ORDER.PIC,NENS]: INITIAL VALUES OF THE AR PROCESS
###		Y.DATE  [NLEAD,NENS]: DATES ASSOCIATED WITH PROCESS (THE FIRST ORDER.PIC LEADS ARE ICS)
###		F.TS    [NLEAD,NENS,NFOR]: FORCING; IGNORED IF F.TS=NA; NFOR = NUMBER OF FORCINGS
###     FNAMES  [NFOR]: FORCING NAMES
###		BETA.HAT[1+ORDER.PIC*MPERIOD+2*NHARM+NFOR]: COEFFICIENTS OF ARX-CYCLOSTATIONARY MODEL
###     NLEAD:		NUMBER OF LEADS
###     ORDER.PIC:	ORDER OF THE AUTOREGRESSIVE MODEL
###     NHARM:		NUMBER OF ANNUAL HARMONICS
###  	MPERIOD:	PERIOD OF THE CYCLOSTATIONARY CYCLE 
### 	NOISE.VAR:  VARIANCE OF THE NOISE (DEFAULT = NA GIVES NOISE-FREE SOLUTION)

if (length(y.date) %% nlead != 0) stop('y.ts dimensions not integer multiple of nlead')

nens        = length(y.date)/nlead
dim(y.date) = c(nlead,nens)

if (nens != length(initial)/order.pic) stop('ensemble size is inconsistent in simulate.arx.cyclo')

if (all(is.na(f.ts))) nfor = 0 else {
	nfor      = length(fnames)
	dim(f.ts) = c(nlead,nens,nfor)
}

if (length(beta.hat) != 1 + order.pic * mperiod + 2 * nharm + nfor) stop('beta.hat inconsistent with other dimensions')

if (order.pic > 0) {
	beta.ar.get  = array(NA,dim=c(mperiod,order.pic))
	for (lag in 1:order.pic) for (m in 1:mperiod) beta.ar.get[m,lag] = beta.hat[names(beta.hat) == paste('lag',lag,'M',m,sep='')]
}
if (nharm     > 0) {
	beta.cyc.get = array(NA,dim=c(nharm,2))
	for (nh in 1:nharm) beta.cyc.get[nh,1] = beta.hat[names(beta.hat) == paste('cos',nh,sep='')]
	for (nh in 1:nharm) beta.cyc.get[nh,2] = beta.hat[names(beta.hat) == paste('sin',nh,sep='')]
}
if (nfor      > 0) {
	beta.for.get = array(NA,dim=c(nfor))
	for (nf in 1:nfor) {
		npic             = which(names(beta.hat) == fnames[nf])
		if (length(npic) != 1) stop('cannot find unique forcing in simulate.arx.cyclo')
		beta.for.get[nf] = beta.hat[npic]
	}
}

y.sim = array(NA,dim=c(nlead,nens))
y.sim[1:order.pic,] = initial


for (n in (1+order.pic):nlead) {
	y.sim[n,] = beta.hat['intercept']
	if (order.pic > 0) for (lag in 1:order.pic) {
		phase = (month(y.date[n,]) - 1) %% mperiod + 1
		y.sim[n,] = y.sim[n,] + y.sim[n-lag,] * beta.ar.get[phase,lag]
	}
 	if (nharm > 0) for (nh in 1:nharm) {
		omega.t   = 2 * pi * month(y.date[n,]) * nh / 12
		y.sim[n,] = y.sim[n,] + cos(omega.t) * beta.cyc.get[nh,1]
		y.sim[n,] = y.sim[n,] + sin(omega.t) * beta.cyc.get[nh,2]
	}
	if (nfor > 0) for (nf in 1:nfor) {
		y.sim[n,] = y.sim[n,] + f.ts[n,,nf] * beta.for.get[nf]
	}
	if (!is.na(noise.var)) y.sim[n,] = y.sim[n,] + rnorm(nens,sd=sqrt(noise.var))
	# if (any(is.na(y.sim[n,]))) 	stop('undefined values detected in simulation')
}

y.sim

}