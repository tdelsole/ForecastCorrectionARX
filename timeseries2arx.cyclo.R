timeseries2arx.cyclo = function(y.ts,y.date,f.ts,fnames,nlead,order.pic,nharm,mperiod=12) {
#### GIVEN TIMESERIES Y, CREATE MATRICES FOR ARX(P) WITH ANNUAL CYCLE AND FORCING
#### INPUT:
####	Y.TS   [NLEAD,NENS]: TIME SERIES OF LENGTH NLEAD; NENS NUMBER OF ENSEMBLE MEMBERS
####    Y.DATE [NLEAD,NENS]: DATES ASSOCIATED WITH Y
####	F.TS   [NLEAD,NENS,NFOR]: FORCING; IGNORED IF F.TS = NA; NFOR = NUMBER OF FORCINGS
####    FNAMES[NFOR]: NAME OF FORCINGS; IGNORED IF F.TS = NA
####    NLEAD:		NUMBER OF LEADS
#####   ORDER.PIC:	ORDER OF THE AUTOREGRESSIVE MODEL
#####   NHARM:		NUMBER OF ANNUAL HARMONICS
#####	MPERIOD:	PERIOD OF THE CYCLOSTATIONARY CYCLE (DEFAULT = 12)
##### OUTPUT:LIST$
#####	Y.LHS[NTOT]: TIME SERIES AFTER OMITTING FIRST ORDER.PIC STEPS
#####	Y.LAG[NTOT,ORDER.PIC*PERIOD]: LAGGED PREDICTOR MATRIX FOR Y; EQUALS NULL IF ORDER.PIC = 0
#####	Y.CYC[NTOT,2*NHARM]:    COSINE/SINE PREDICTOR MATRIX; EQUALS NULL IF NHARM = 0
#####	Y.FOR[NTOT,NFOR]: FORCING TIME SERIES
#####	JVEC [NTOT]: VECTOR OF ONES (FOR THE INTERCEPT)
##### COMMENTS
##### 1) THE ANNUAL HARMONIC IS ALWAYS PERIOD 12.  IT DOES NOT NECESSARILY EQUAL MPERIOD

if (length(y.ts) %% nlead != 0) stop('y.ts dimensions not integer multiple of nlead')
if (length(y.ts) != length(y.date)) stop('y.ts and y.date must have same dimensions')

x.names     = NULL
nens        = length(y.ts)/nlead
dim(y.ts)   = c(nlead,nens)
dim(y.date) = c(nlead,nens)


if (all(is.na(f.ts))) nfor = 0 else {
	nfor      = length(fnames)
	dim(f.ts) = c(nlead*nens,nfor)
}

phase       = (month(y.date) - 1) %% mperiod + 1
lead.all    = array(1:nlead,dim=c(nlead,nens))

y.lhs         = NULL
num.per.phase = as.numeric(rep(NA,mperiod))
for (m in 1:mperiod) {
	npic    = which( phase == m & lead.all > order.pic)
	y.lhs   = c(y.lhs,y.ts[npic]) 
	num.per.phase[m] = length(npic)
}

x.ar        = NULL
names.ar    = NULL
if (order.pic > 0) for (m in 1:mperiod) {
	npic    = which( phase == m & lead.all > order.pic)
	x.all.lag = NULL
	for (lag in 1:order.pic) x.all.lag = cbind(x.all.lag,y.ts[npic-lag])
	
	if (m == 1      ) zero.pre = NULL else zero.pre = array(0,dim=c(length(npic),order.pic*(m-1)))
	if (m == mperiod) zero.suf = NULL else zero.suf = array(0,dim=c(length(npic),order.pic*(mperiod-m)))
	x.ar = rbind(x.ar,cbind(zero.pre,x.all.lag,zero.suf))

	names.ar = c(names.ar,paste('lag',1:order.pic,'M',m,sep=''))
}

x.cyc     = NULL
names.cyc = NULL
if (nharm > 0) for (m in 1:mperiod) {
	npic    = which(phase == m & lead.all > order.pic)
	x.temp  = NULL
	omega.t = 2 * pi * (month(y.date[npic])) / 12
	for (nh in 1:nharm) x.temp = cbind(x.temp,cos(omega.t * nh),sin(omega.t * nh))
	x.cyc = rbind(x.cyc,x.temp)
	if (m == 1) for (nh in 1:nharm) names.cyc = c(names.cyc,paste('cos',nh,sep=''),paste('sin',nh,sep=''))		
}

x.for     = NULL
if (nfor > 0) for (m in 1:mperiod) {
	npic  = which(phase == m & lead.all > order.pic)
	x.for = rbind(x.for,f.ts[npic,,drop=FALSE])
}

jvec = as.matrix(rep(1,length(y.lhs)))

if (order.pic > 0) colnames(x.ar)  = names.ar
if (nharm     > 0) colnames(x.cyc) = names.cyc
if (nfor      > 0) colnames(x.for) = fnames
colnames(jvec)  = 'intercept'


list(y.lhs=y.lhs,x.ar=x.ar,x.cyc=x.cyc,x.for=x.for,jvec=jvec,
   num.per.phase = num.per.phase)

}