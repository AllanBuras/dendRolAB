###a function to conveniently generate master chronos per tree



###function chron by stem

chron.by.stem<-function(rwl,str.id1,str.id2)
	{
	IDS<-substr(colnames(rwl),str.id1,str.id2)
	AV.RWL<-matrix(nrow=nrow(rwl),ncol=length(unique(IDS)))
	for(i in 1:ncol(AV.RWL))
		{
		AV.RWL[,i]<-chron(as.matrix(rwl[,which(IDS==unique(IDS)[i])]),biweight=F)[,1]
		}
	rownames(AV.RWL)<-rownames(rwl)
	colnames(AV.RWL)<-unique(IDS)
	return(as.data.frame(AV.RWL))
	}

###meta-function needed for PCGA

.polar.trans<-function(matrix)
	{
	POL_MAT<-matrix(nrow=nrow(matrix),ncol=2)
	for(i in 1:nrow(matrix))
		{
		POL_MAT[i,2]<-dist(rbind(c(0,0),matrix[i,1:2]),method="euclidean")
		POL_MAT[i,1]<-atan2(matrix[i,2],matrix[i,1])
		}
	return(POL_MAT)
	}

###code for PCGA, provide this function with an rwl-data.frame

pcga<-function(rwl,plot.CI=T)
	{
	CI<-common.interval(rwl,make.plot=plot.CI)
	PCA<-prcomp(CI,scale=T)
	IMP<-summary(PCA)$imp[,1:4]
	POL_COORD<-.polar.trans(PCA$rot[,1:2])
	POL_COORD2<-.polar.trans(-PCA$rot[,1:2])
	if(max(dist(POL_COORD[,1]))>max(dist(POL_COORD2[,1])))
		{
		POL_COORD<-POL_COORD2
		PCA$rot<--PCA$rot
		}
	list(pca=PCA,imp=IMP,rank=order(POL_COORD[,1]),pol.coord=POL_COORD,pop=rwl[,colnames(CI)],period=rownames(CI))
	}

plot_pcga<-function(x,col.vec=NULL,...)
	{
	plot(0,0,type="n",xlim=range(-c(x$pca$rot[,1],x$pca$rot[,2]),c(x$pca$rot[,1],x$pca$rot[,2])),
	ylim=range(-c(x$pca$rot[,1],x$pca$rot[,2]),c(x$pca$rot[,1],x$pca$rot[,2])),xlab=paste("PC1 (rsq):",round(x$imp[2,1],2),sep=" "),ylab=paste("PC2 (rsq):",round(x$imp[2,2],2),sep=" "),...)
	if(length(col.vec)==0)
		{
		for(i in 1:nrow(x$pca$rot))
			{
			arrows(0,0,x$pca$rot[x$rank[i],1],x$pca$rot[x$rank[i],2],col=hsv(i/(nrow(x$pca$rot)*1.2)),lwd=2,length=0.1)
			}
		}
	if(length(col.vec)>0)
		{
		for(i in 1:nrow(x$pca$rot))
			{
			arrows(0,0,x$pca$rot[x$rank[i],1],x$pca$rot[x$rank[i],2],col=col.vec[x$rank[i]],lwd=2,length=0.1)
			}
		}
	}


identRank<-function(pcga.obj)
	{
	plot_pcga(pcga.obj)
	print("please select tree by clicking on the head of the loading-arrows")
	LOC<-locator(1)
	DIST<-vector(mode="numeric")
	for(i in 1:nrow(pcga.obj$pca$rot[,1:2]))
		{
		DIST[i]<-dist(rbind(LOC,pcga.obj$pca$rot[pcga.obj$rank[i],1:2]))
		}
	text(x=LOC[1],y=LOC[2],labels=which.min(DIST),pos=4)
	text(x=pcga.obj$pca$rot[pcga.obj$rank[c(1,length(pcga.obj$rank))],1],
	y=pcga.obj$pca$rot[pcga.obj$rank[c(1,length(pcga.obj$rank))],2],labels=c(1,length(pcga.obj$rank)),pos=4)
	return(which.min(DIST))
	}


###BTFS

.conf.test<-function(x,y)
	{
	SPLIT<-length(x)/2
	SAMPLE<-rbind(sample(1:SPLIT,size=length(x),replace=T),sample((SPLIT+1):length(x),size=length(x),replace=T))
	Y1<-y[SAMPLE[1,]]
	Y2<-y[SAMPLE[2,]]
	X1<-x[SAMPLE[1,]]
	X2<-x[SAMPLE[2,]]
	M1<-lm(Y1~X1)
	M2<-lm(Y2~X2)
	STATS_MAT<-matrix(nrow=2,ncol=6)
	if(length(summary(M1)$fstat)>0&&length(summary(M2)$fstat)>0)
		{
		STATS_MAT[1,1]<-1-pf(summary(M1)$fstat[1],summary(M1)$fstat[2],summary(M1)$fstat[3])
		STATS_MAT[2,1]<-1-pf(summary(M2)$fstat[1],summary(M2)$fstat[2],summary(M2)$fstat[3])
		}
	else
		{
		STATS_MAT[,1]<-rep(1,2)
		}
	STATS_MAT[,2]<-c(M1$coef[1],M2$coef[1])		
	STATS_MAT[,3]<-c(M1$coef[2],M2$coef[2])		
	STATS_MAT[,4]<-c(summary(M1)$r.sq,summary(M2)$r.sq)
	Pred_M1_P2<-M1$coef[1]+M1$coef[2]*X2
	Pred_M2_P1<-M2$coef[1]+M2$coef[2]*X1
	SSY_1<-sum((Y1-mean(Y1,na.rm=T))^2,na.rm=T)
	SSY_2<-sum((Y2-mean(Y2,na.rm=T))^2,na.rm=T)
	SSYRE_1<-sum((Y1-mean(Y2,na.rm=T))^2,na.rm=T)
	SSYRE_2<-sum((Y2-mean(Y1,na.rm=T))^2,na.rm=T)
	SSR1<-sum((Y1-Pred_M2_P1)^2,na.rm=T)
	SSR2<-sum((Y2-Pred_M1_P2)^2,na.rm=T)
	RE1<-1-SSR1/SSYRE_1
	RE2<-1-SSR2/SSYRE_2
	CE1<-1-SSR1/SSY_1
	CE2<-1-SSR2/SSY_2
	STATS_MAT[,5]<-c(RE1,RE2)
	STATS_MAT[,6]<-c(CE1,CE2)
	return(STATS_MAT)
	}

btfs<-function(x,y,boot.n=1000)
	{
	BOOT_MAT<-matrix(ncol=6,nrow=boot.n)
	COMP_MAT<-matrix(ncol=6,nrow=boot.n)
	for(i in 1:boot.n)
		{
		ITER<-.conf.test(x,y)
		COMP_MAT[i,1]<-ITER[1,2]
		COMP_MAT[i,2]<-ITER[2,2]
		COMP_MAT[i,3]<-ITER[1,3]
		COMP_MAT[i,4]<-ITER[2,3]
		COMP_MAT[i,5]<-ITER[1,4]
		COMP_MAT[i,6]<-ITER[2,4]
		BOOT_MAT[i,1]<-max(ITER[,1])
		BOOT_MAT[i,2]<-(ITER[1,2]/ITER[2,2])
		BOOT_MAT[i,3]<-(ITER[1,3]/ITER[2,3])
		BOOT_MAT[i,4]<-(ITER[1,4]/ITER[2,4])
		BOOT_MAT[i,5]<-min(ITER[,5])
		BOOT_MAT[i,6]<-min(ITER[,6])
		}
	colnames(COMP_MAT)<-c("intercept P1","intercept P2","slope P1","slope P2","rsq P1","rsq P2")
	MEANS<-apply(BOOT_MAT,2,median)
	STATS<-length(which(BOOT_MAT[,1]<=0.05))/boot.n
	ECDF2<-ecdf(BOOT_MAT[,2])
	ECDF3<-ecdf(BOOT_MAT[,3])
	ECDF4<-ecdf(BOOT_MAT[,4])
	ECDF5<-ecdf(BOOT_MAT[,5])
	ECDF6<-ecdf(BOOT_MAT[,6])
	if(MEANS[2]>=1){STATS[2]<-2*ECDF2(1)}
	if(MEANS[2]<1){STATS[2]<-2*(1-ECDF2(1))}
	if(MEANS[3]>=1){STATS[3]<-2*ECDF3(1)}
	if(MEANS[3]<1){STATS[3]<-2*(1-ECDF3(1))}
	if(MEANS[4]>=1){STATS[4]<-2*ECDF4(1)}
	if(MEANS[4]<1){STATS[4]<-2*(1-ECDF4(1))}
	STATS[5]<-1-ECDF5(0)
	STATS[6]<-1-ECDF6(0)
	MEANS<-round(MEANS,digits=2)
	STATS<-round(STATS,digits=3)
	OUT<-rbind(MEANS[2:6],STATS[2:6])
	colnames(OUT)<-c("intercept-ratio","slope-ratio","rsq-ratio","RE","CE")
	rownames(OUT)<-c("bootstrapped estimate","p-value")
	if(MEANS[1]>0.05){print(paste("Warning: period-specific regressions are in ",(1-STATS[1])*100,"% of cases insignificant",sep=""))}
	list(stats=OUT,boot.mat=COMP_MAT)
	}


plot_btfs<-function(x,...)
	{
	par(mfrow=c(1,3))
	for(i in 1:3)
		{
		boxplot(x$boot.mat[,(i*2)-(1:0)],main=c("intercept","slope","rsq")[i],cex.lab=1.5,cex.axis=1.5,cex.main=2)
		VEC<-as.vector(x$boot.mat[,(i*2)-(1:0)])
		if(x$stats[2,i]>0.05){SIGN<-""}
		if(x$stats[2,i]<0.05){SIGN<-"*"}
		if(x$stats[2,i]<0.01){SIGN<-"**"}
		if(x$stats[2,i]<0.001){SIGN<-"***"}
		text(x=1.5,y=seq(min(VEC),max(VEC),diff(range(VEC))/100)[95],SIGN,cex=3)
		}
	}



###BSGC as implemented in Buras et al. 2022:

####BSGC-code 20220314


.conflate<-function(p1,p2,w1,w2)
{
  p12<-(p1^w1*p2^w2)/(p1^w1*p2^w2+(1-p1)^w1*(1-p2)^w2)
  return(p12)
}

.bsgcP<-function(rwl,p.thresh=0.05,make.plot=T,maxlag=20,w1=0.93,w2=0.53,...)
{
  CHRON<-apply(rwl,1,mean)
  GC<-apply(rwl,2,diff)
  M.GC<-apply(GC,1,mean)
  PGSGC.MAT<-matrix(nrow=nrow(rwl)-1,ncol=maxlag)
  rownames(PGSGC.MAT)<-rownames(rwl)[2:nrow(rwl)]
  PSGC.MAT<-matrix(nrow=nrow(rwl)-1,ncol=maxlag)
  rownames(PSGC.MAT)<-rownames(rwl)[2:nrow(rwl)]
  for(i in 1:maxlag)
  {
    GC.I<-apply(rwl,2,diff,lag=i)
    M.GC.I<-apply(GC.I,1,mean)
    PSGC.I<-pnorm(scale(M.GC.I))
    PSGC.MAT[1:(nrow(rwl)-i),i]<-(PSGC.I)
    PGSGC.ALL<-vector(mode="numeric")
    for(j in 1:(nrow(rwl)-i))
    {
      DIFF.ALL<-CHRON[-j]-CHRON[j]
      PGSGC.ALL[j]<-pnorm((M.GC.I[j]-mean(DIFF.ALL))/sd(DIFF.ALL))
    }
    PGSGC.MAT[1:(nrow(rwl)-i),i]<-PGSGC.ALL
  }
  COMB.P<-.conflate(PGSGC.MAT,PSGC.MAT,w1,w2)
  rownames(COMB.P)<-rownames(rwl)[2:nrow(rwl)]
  list(pvals=COMB.P)
}



bsgc<-function(rwl,p.thresh=0.05,make.plot=T,maxlag=20,rm.succ=TRUE,w1=0.93,w2=0.53,...)###bias-corrected standardized growth change
{
  if(length(which(is.na(rwl)))!=0)
    {
    stop("rwl must not contain NAs. Only apply BSGC to common overlap period!")
    }
  if(maxlag>nrow(rwl)/2){maxlag<-nrow(rwl)/2;warning("maximum lag was set to half the number of observations!")}
  CHRON<-apply(rwl,1,mean)
  GC<-apply(rwl,2,diff)
  M.GC<-apply(GC,1,mean)
  SHAP<-shapiro.test(M.GC)
  if(SHAP$p.val<0.05){warning("shapiro-test p.value of growth changes < 0.05. Derived parameters are not reliable!")}
  PGSGC.MAT<-matrix(nrow=nrow(rwl)-1,ncol=maxlag)
  rownames(PGSGC.MAT)<-rownames(rwl)[2:nrow(rwl)]
  colnames(PGSGC.MAT)<-paste("lag",0:(ncol(COMB.P)-1))
  PSGC.MAT<-matrix(nrow=nrow(rwl)-1,ncol=maxlag)
  rownames(PSGC.MAT)<-rownames(rwl)[2:nrow(rwl)]
  colnames(PSGC.MAT)<-paste("lag",0:(ncol(COMB.P)-1))
  for(i in 1:maxlag)
  {
    GC.I<-apply(rwl,2,diff,lag=i)
    M.GC.I<-apply(GC.I,1,mean)
    PSGC.I<-pnorm(scale(M.GC.I))
    PSGC.MAT[1:(nrow(rwl)-i),i]<-(PSGC.I)
    PGSGC.ALL<-vector(mode="numeric")
    for(j in 1:(nrow(rwl)-i))
    {
      DIFF.ALL<-CHRON[-j]-CHRON[j]
      PGSGC.ALL[j]<-pnorm((M.GC.I[j]-mean(DIFF.ALL))/sd(DIFF.ALL))
    }
    PGSGC.MAT[1:(nrow(rwl)-i),i]<-PGSGC.ALL
  }
  COMB.P<-.conflate(PGSGC.MAT,PSGC.MAT,w1,w2)
  rownames(COMB.P)<-rownames(rwl)[2:nrow(rwl)]
  colnames(COMB.P)<-paste("lag",0:(ncol(COMB.P)-1))
  if(length(which(COMB.P[,1]<p.thresh/2))>0)
  {
    NEG<-COMB.P[which(COMB.P[,1]<p.thresh/2),1]	
    names(NEG)<-rownames(COMB.P)[which(COMB.P[,1]<p.thresh/2)]
  }
  else(NEG<-NULL)
  if(length(which(COMB.P[,1]>(1-p.thresh/2)))>0)
  {
    POS<-COMB.P[which(COMB.P[,1]>(1-p.thresh/2)),1]
    names(POS)<-rownames(COMB.P)[which(COMB.P[,1]>(1-p.thresh/2))]
  }
  else(POS<-NULL)
  ###now correct computations for successive pointer years to lower the likelihood of artefacts
  if(rm.succ==T)
  {
    if(length(NEG)>0&&length(POS)>0)
    {
      DIFF.POS<-vector(mode="numeric")
      DIFF.NEG<-vector(mode="numeric")
      for(i in 1:length(NEG))
      {
        DIFF.NEG[i]<-length(intersect(as.numeric(names(NEG)[i])-as.numeric(names(POS)),1))
      }
      names(DIFF.NEG)<-names(NEG)
      for(i in 1:length(POS))
      {
        DIFF.POS[i]<-length(intersect(as.numeric(names(POS)[i])-as.numeric(names(NEG)),1))
      }	
      names(DIFF.POS)<-names(POS)
      if(length(which(DIFF.NEG==1))>0)
      {
        EXCLUDE<-vector(mode="numeric")
        for(i in which(DIFF.NEG==1))
        {
          RWL<-rwl
          RWL[as.character(as.numeric(names(NEG)[i])-1),]<-1
          COMB.P.TEST<-.bsgcP(RWL,maxlag=maxlag,...)
          if(COMB.P.TEST$pvals[names(NEG)[i],1]>p.thresh/2){COMB.P[names(NEG)[i],]<-COMB.P.TEST$pvals[names(NEG)[i],];EXCLUDE<-c(EXCLUDE,i)}
        }
        if(length(EXCLUDE)>0){NEG<-NEG[-EXCLUDE]}
      }
      if(length(which(DIFF.POS==1))>0)
      {
        EXCLUDE<-vector(mode="numeric")
        for(i in which(DIFF.POS==1))
        {
          RWL<-rwl
          RWL[as.character(as.numeric(names(POS))-1)[i],]<-1
          COMB.P.TEST<-.bsgcP(RWL,maxlag=maxlag,...)
          if(COMB.P.TEST$pvals[names(POS)[i],1]<(1-p.thresh/2)){COMB.P[names(POS)[i],]<-COMB.P.TEST$pvals[names(POS)[i],];EXCLUDE<-c(EXCLUDE,i)}
        }
        if(length(EXCLUDE)>0){POS<-POS[-EXCLUDE]}
      }
    }
  }
  if(length(NEG)>0)
  {
    LEG.NEG<-rep(0,length(NEG))
    names(LEG.NEG)<-names(NEG)
    for(i in 1:(length(LEG.NEG)))
    {
      for(j in 1:max(which(is.finite(COMB.P[names(LEG.NEG)[i],]))))
      {
        if(COMB.P[names(LEG.NEG)[i],j]>p.thresh/2)
        {
          LEG.NEG[i]<-j-2
          break
        }
        if(j==max(which(is.finite(COMB.P[names(LEG.NEG)[i],])))){LEG.NEG[i]<-j-1}
      }
    }
  }
  else(LEG.NEG<-NULL)
  if(length(POS)>0)
  {
    LEG.POS<-rep(0,length(POS))
    names(LEG.POS)<-names(POS)
    for(i in 1:(length(LEG.POS)))
    {
      for(j in 1:max(which(is.finite(COMB.P[names(LEG.POS)[i],]))))
      {
        if(COMB.P[names(LEG.POS)[i],j]<(1-p.thresh/2))
        {
          LEG.POS[i]<-j-2
          break
        }
        if(j==max(which(is.finite(COMB.P[names(LEG.POS)[i],])))){LEG.POS[i]<-j-1}
      }
    }
  }
  else(LEG.POS<-NULL)
  if(length(which(LEG.NEG>0))>0)
  {
    LONG.LEGS<-which(LEG.NEG>0)
    YEARS<-vector(mode="numeric")
    for(i in LONG.LEGS)
    {
      YEARS<-c(YEARS,intersect(as.character(as.numeric(names(NEG)[i])+1:LEG.NEG[i]),names(NEG)[-i]))
    }
    YEARS<-unique(YEARS)
    if(length(YEARS)>0)
    {
      for(i in 1:length(YEARS))
      {
        NEG<-NEG[-which(names(NEG)==YEARS[i])]
        LEG.NEG<-LEG.NEG[-which(names(LEG.NEG)==YEARS[i])]
      }
    }
  }
  if(length(which(LEG.POS>0))>0)
  {
    LONG.LEGS<-which(LEG.POS>0)
    YEARS<-vector(mode="numeric")
    for(i in LONG.LEGS)
    {
      YEARS<-c(YEARS,intersect(as.character(as.numeric(names(POS)[i])+1:LEG.POS[i]),names(POS)[-i]))
    }
    YEARS<-unique(YEARS)
    if(length(YEARS)>0)
    {
      for(i in 1:length(YEARS))
      {
        POS<-POS[-which(names(POS)==YEARS[i])]
        LEG.POS<-LEG.POS[-which(names(LEG.POS)==YEARS[i])]
      }
    }
  }
  if(make.plot==T)
  {
    plot(CHRON~rownames(rwl),type="l",lwd=2,...)
    points(as.numeric(names(NEG)),CHRON[names(NEG)],pch=16,col="red",cex=2)
    points(as.numeric(names(POS)),CHRON[names(POS)],pch=16,col="blue",cex=2)
    for(i in which(LEG.NEG>0))
    {
      points((as.numeric(names(NEG)[i])+1):(as.numeric(names(NEG)[i])+LEG.NEG[i]),
             CHRON[as.character((as.numeric(names(NEG)[i])+1):(as.numeric(names(NEG)[i])+LEG.NEG[i]))],pch=16,col="orange",cex=2)
    }
    for(i in which(LEG.POS>0))
    {
      points((as.numeric(names(POS)[i])+1):(as.numeric(names(POS)[i])+LEG.POS[i]),
             CHRON[as.character((as.numeric(names(POS)[i])+1):(as.numeric(names(POS)[i])+LEG.POS[i]))],pch=16,col="dodgerblue",cex=2)
    }
  }
  list(neg=NEG,pos=POS,pvals.sgc=PSGC.MAT,pvals.gsgc=PGSGC.MAT,pvals.conf=COMB.P,def.neg=LEG.NEG,def.pos=LEG.POS)
}















