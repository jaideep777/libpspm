### initialisation

# initialisation of age-specific vectors
L<-Linf*(1-exp(-k*(1:ages-a0)))		# length-at-age
W<-theta*L^g						# weight-at-age
# natural mortality - fill in with the last value
M<-c(M,rep(M[length(M)],times=ages-length(M))); names(M)<-1:ages
# maturity-at-age - fill in 1's
mat<-c(mat,rep(1,times=ages-length(mat))); names(mat)<-1:ages

# initialisation of other data structures
N<-matrix(0,nrow=years+1,ncol=ages)	# numbers at age
N[1,1]<-N.ini
for(i in 2:ages) N[1,i]<-N.ini*exp(-sum(M[1:(i-1)])) 	# numbers at age, year 1
Res<-data.frame(year=1:years,N=NA,TSB=NA,SSB=NA,Y=NA,R=NA,E=NA,Ereal=NA,Fratio=NA,F.hypo=NA,
	D.sea=NA,D.shr=NA,D.tot=NA,P.sea=NA,P.shr=NA,P.tot=NA) # results
  
# find out carrying capacity (Z=M)
# independent of fishing parameters, so this needs to be done only once
N.<-matrix(0,nrow=1051,ncol=ages)	# this calculation should not be influenced by the length of the main simulation
N.[1,]<-N[1,]
s<-exp(-M)		# survival
for(year in 1:(nrow(N.)-1)){				
	N.[year+1,2:ages]<-N.[year,1:(ages-1)]*s[1:(ages-1)]
	SSB<-sum(N.[year,]*W*mat)
	N.[year+1,1]<-r0*SSB/(1+SSB/Bhalf)
}
K<-apply(N.[51:1050,],2,mean) # skip first 50 years as the transient
rm(N.)

### functions
# image that uses gray scale as the default
Image<-function(...){image(...,col=gray((32:0)/32)); box()}

# function needed to find fishing mortality corresponding to a certain effort
simul<-function(F.max,parms){ # parms=N.rel,E.real
	E<-parms[1]^(1-b)*F.max*(exp(-(F.max+M[ages])*(1-b))-1)/(q*(F.max+M[ages])*(b-1))
	(E-parms[2])^2
}

# Baranov catch equation - solve F when M, Y and B0 are known
Bce<-function(F,parms){
	Y<-parms[1,1]
	N0<-parms[,2]
	Y.<-sum(F/(F+M)*(1-exp(-F-M))*N0*W)
	sum((Y-Y.)^2)
}

# main simulation function
simulate<-function(h.prop,LF50,display=FALSE){
	F.max<- -log(1-h.prop)
	age.fish<-min((1:ages)[1/(1+exp(-SF*(L-LF50)))>0.5])
	F2E.factor<-ifelse(b==1,F.max/q,F.max*(exp(-(F.max+M[ages])*(1-b))-1)/(q*(F.max+M[ages])*(b-1))) # Liu & Heino 2014
	for(year in 1:years){ # 
		if(h.prop>0){
			# determine the effective fishing mortality
			N.rel<-mean(N[year,age.fish:ages]/K[age.fish:ages])
			# cat("N.rel = ", N.rel, "\n")
			E<-N.rel^(1-b)*F2E.factor			# effort corresponding to F (F2E.factor contains only terms that are year-independent)
			D.sea.demand<-d.sea*E				# corresponding employment demand
			D.eff<-d.max*D.sea.demand/(d.max+D.sea.demand)	# effective man hours
			E.real<-D.eff/d.sea					# realized effort (~vessel days)
			F.real<-ifelse(b==1,q*E.real,optimize(f=simul,interval=c(0,F.max),parms=c(N.rel,E.real),tol=0.0001)$minimum)
			F.ratio = F.real/F.max
			#F.real=F.max  # JAI: Remove effort dynamics for testing
			# population dynamics
			F<-F.real/(1+exp(-SF*(L-LF50)))		# fishing mortality at age
		}
		else{
			F<-0; E<-0; E.real<-0; F.real<-0
			D.sea.demand<-0
			F.ratio = 0
		}	
		Z<-M+F								# total mortality
		s<-exp(-Z)							# survival
		Y<-sum(F/Z*(1-s)*N[year,]*W)
		if(Y>0){ # calculate the hypothetical F corresponding to unselective harvesting of all age groups
			F.hypo<-optimize(f=Bce,interval=c(0,F.max),parms=data.frame(Y,N[year,]))$minimum
		}else F.hypo<-0
		TSB<-sum(N[year,]*W)
		SSB<-sum(N[year,]*W*mat)
		N[year+1,1]<-r0*SSB/(1+SSB/Bhalf)#*exp(rnorm(n=1,mean=-sf^2/2,sd=sf))  # JAI:: Remove environmental stoch for testing
		N[year+1,2:ages]<-N[year,1:(ages-1)]*s[1:(ages-1)]
		# economic calculations
		profit.sea<-Y*price.sea - scale.catch*(D.sea.demand*s.sea + E*varC.sea + fixedC.sea)
		profit.shr<-Y*(price.shr-price.sea) - Y*d.shr*s.shr - scale.catch*fixedC.shr
		Res[year,]<-c(year,sum(N[year,]),TSB,SSB,Y,N[year+1,1],E,E.real,F.ratio,F.hypo,
			D.sea.demand,Y*d.shr,D.sea.demand+Y*d.shr,profit.sea,profit.shr,profit.sea+profit.shr)
	}
	if(display==TRUE){
		op<-par(mfrow=c(2,2),pty="s",mar=c(4,4,.5,.5))
		plot(apply(N,1,sum)/sum(K),type="l",xlab="",ylab="Total numbers (relative to K)")
		plot(Res$SSB,type="l",xlab="",ylab="SSB (kg)")
		plot(Res$Y,type="l",xlab="",ylab="Total catch weight")
		Image(x=1:years,y=1:ages,z=log(N[1:years,]),xlab="",ylab=""); legend("top","Age structure",bty="n") # non-normalized
		par(op)
	}
	list(N=N[1:years,],summaries=Res)
}
