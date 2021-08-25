
##FOLLOWING 'A system of model fertility schedules with graphically intuitive parameters' BY Carl Schmertmann (2003)
##	(QUADRATIC SPLINE FERTILITY SCHEDULE)
##	VIA https://www.demographic-research.org/volumes/vol9/5/
##	INCLUDING ITS Appendix B: Non-matrix calculation of spline coefficients
##	AND ADDITIONAL EXCEL FILE VIA https://www.demographic-research.org/volumes/vol9/5/files/FILE%20I%20QS%20Model.xls
##	AUGUST 2021

##MAKE A FUNCTION FOR THE MODEL
QSFert<-function(R,alpha,P,H)	{
	##KNOTS (t) 
	W<-ifelse(.75 < .25+.025*(P-alpha), 
		.75,
		.25+.025*(P-alpha))
	
	beta<-ifelse(H+(1/3)*(H-P) > 50, 
		H+(1/3)*(H-P),
		ifelse(H+(3)*(H-P) < 50,
			H+(3)*(H-P),
			50))
	
	t0<-alpha
	t1<-(1-W)*alpha+W*P
	t2<-P
	t3<-(P+H)/2
	t4<-(H+beta)/2
	t5<-beta
	
	##COEFFICIENTS (FOR KNOTS 0 TO 4) BY NON-MATRIX FORMULAS (theta)
	theta0<-1/(W*(P-alpha)^2)
	theta1<--theta0/(1-W)
	
		##INTERMEDIATE VARIABLES
		Za<-(H-t2)^2
		Zb<-(H-t3)^2
		Zc<-(beta-t2)^2
		Zd<-(beta-t3)^2
		Ze<-(beta-t4)^2
		Zf<-2*(beta-t2)
		Zg<-2*(beta-t3)
		Zh<-2*(beta-t4)
		Z1<-0.5-(theta0*(H-alpha)^2 + theta1*(H-t1)^2)
		Z2<-0-(theta0*(beta-alpha)^2 + theta1*(beta-t1)^2)
		Z3<-0-(2*theta0*(beta-alpha) + 2*theta1*(beta-t1))
		DENOM<-Za*Zd*Zh - Za*Ze*Zg - Zc*Zb*Zh + Zf*Zb*Ze
	
	theta2<-(Z1*(Zd*Zh - Ze*Zg) - Z2*(Zb*Zh) + Z3*(Zb*Ze)) / DENOM
	theta3<-(Z1*(Ze*Zf - Zc*Zh) + Z2*(Za*Zh) - Z3*(Za*Ze)) / DENOM
	theta4<-(Z1*(Zc*Zg - Zd*Zf) + Z2*(Zb*Zf - Za*Zg) + Z3*(Za*Zd - Zb*Zc)) / DENOM
	
	t<-c(t0,t1,t2,t3,t4,t5)
	theta<-c(theta0,theta1,theta2,theta3,theta4)
	
	##APPLY THE MODEL
	x<-array(seq(0,65,1))
	BasisFunctions<-data.frame(x)
	BasisFunctions$t0<-BasisFunctions$t1<-BasisFunctions$t2<-BasisFunctions$t3<-BasisFunctions$t4<-0
	fx<-array(0,length(x))
	
	for (i in 1:nrow(BasisFunctions)) {
		for (j in 2:6) {
			ifelse(BasisFunctions$x[i]>t[j-1],
				BasisFunctions[i,j]<-(BasisFunctions$x[i]-t[j-1])^2,0)}}
	BasisMatrix<-as.matrix(BasisFunctions[,2:6])
	
	for (i in 1:nrow(BasisFunctions)) {
		ifelse(BasisFunctions$x[i]<t5,
			fx[i]<-R*(BasisMatrix[i,]%*%theta),0)}
	
	return(c(data.frame(t),data.frame(theta),data.frame(x),data.frame(fx)))}
	
##MODEL (AND QSFert()) INPUTS 
R<-.3		#LEVEL (FOR TFR)
alpha<-15	#LOWEST AGE
P<-23		#PEAK AGE
H<-31		#HALF PEAK AGE

##APPLY THE MODEL WITH THE QSFert FUNCTION
QSApplied<-QSFert(R,alpha,P,H)

##PRINT AND PLOT THE OUTPUT
QSApplied$t
QSApplied$theta
sum(QSApplied$fx)
plot(QSApplied$fx,main="Quadratic Spline Fertility Schedule",xlab="Age",ylab="Fertility Rate",axes=FALSE)
  axis(side=2,cex.axis=0.75)
  axis(side=1,at=1:length(QSApplied$x),labels=QSApplied$x,cex.axis=0.75,las=1)

