n.to.monthly = function(date.random,date.ordered,lcheck=FALSE) {
### DETERMINE INDEX N SUCH THAT DATE.ORDERED[N] = DATE.RANDOM
### DATE.ORDERED IS ASSUMED TO BE MONTHLY, SEQUENTIAL, AND CONTAINS ALL TIMES IN DATE.RANDOM
### DATE.RANDOM HAS NO ASSUMED ORDER
### INPUT:
###   DATE.RANDOM[NFROM]: DATES OF THE 'RANDOM' DATA
###   DATE.ORDERED  [NTO]: DATES OF THE 'ORDERED' DATA

## FIND EARLIEST DATE IN DATE.RANDOM
n.first  = which.min(as.numeric(date.random))

## COUNT NUMBER OF MONTHS FROM EARLIEST DATE.RANDOM TO ANY DATE.RANDOM
n.ordered     = 12 * ( year(date.random) - year(date.random[n.first]) ) + (month(date.random) - month(date.random[n.first]))

## FIND DATE.ORDERED THAT MATCHES THE EARLIEST DATE IN DATE.RANDOM
ist      = date.ordered == date.random[n.first]
if (sum(ist) != 1) stop('cannot identify unique date.ordered that matches earliest date.random')

## FIND THE INDICES SUCH THAT DATE.ORDERED[N.ORDERED] = DATE.RANDOM
n.ordered     = which(ist) + n.ordered

if (lcheck) print(all.equal(date.ordered[n.ordered],date.random))

n.ordered
	
}