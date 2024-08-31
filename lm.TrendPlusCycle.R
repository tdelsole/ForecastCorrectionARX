lm.TrendPlusCycle = function(y.ts,y.date,str.date,end.date,nharm=5) {
### FUNCTION TO FIT TIME SERIES TO TREND PLUS ANNUAL CYCLE
### AND THEN PLOT THE TIME SERIES WITH TREND FOR EACH MONTH
### INPUT:
###  Y.TS  [NTIME]:	TIME SERIES OF LENGTH NTIME
###  Y.DATE[NTIME]:	DATE FOR EACH STEP
###  STR.DATE:		STARTING DATE FOR FITTING REGRESSION MODEL
###  END.DATE:		ENDING DATE FOR FITTING REGRESSION MODEL
###  NHARM:			NUMBER OF ANNUAL HARMONICS (DEFAULT = 5)
### COMMENT: 
###   1) USES "DAY" TO MEASURE TREND. ANSWER DIFFERS SLIGHTLY IF "MONTH" IS USED TO MEASURE TREND. 

if (length(y.ts) != length(y.date)) stop('y.date and y.ts have different lengths')
ntot = length(y.ts)

### ALERT IF STR/END DATES EXCEED TIME SERIES PERIOD
# if (y.date[1] < str.date | y.date[ntot] > end.date) print('START/END DATES EXCEED TIME SERIES PERIOD')

### GENERATE LINEAR FUNCTION OF TIME, RE-SCALE TO DECADES, CENTER
y.trend = as.numeric(y.date) / 365.24 / 10
y.trend = y.trend - mean(y.trend)

### SPECIFY ANNUAL HARMONICS
y.harm = NULL
if (nharm >= 1) {
	phase = month(y.date)
	for (nh in 1:nharm) y.harm = cbind(y.harm,cos(2*pi*nh*phase/12),sin(2*pi*nh*phase/12))
}

y      = as.numeric(y.ts)
lgood  = !is.na(y) & y.date >= str.date & y.date <= end.date
ngood  = sum(lgood)

y      = y[lgood]
x      = cbind(y.trend[lgood],y.harm[lgood,])
lm.xy  = lm(y~x)
y.pred = as.numeric( cbind(1,y.trend,y.harm) %*% coef(lm.xy) )
lm.sum = summary(lm.xy)
r.sqr  = lm.sum$r.square
sigma  = sqrt(sum(residuals(lm.xy)^2)/lm.xy$df.residual)
trend  = lm.sum$coef[2,1]
sterr  = lm.sum$coef[2,2]

list(y.pred = y.pred, trend = trend, sterr = sterr, sigma = sigma, r.sqr = r.sqr, ngood = ngood, ntot = ntot)

}