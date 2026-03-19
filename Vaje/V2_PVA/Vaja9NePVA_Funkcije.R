#################FUNKCIJE, POTREBNE ZA VAJE GES (INBREEDING, DRIFT, Ne...)
#################PRITISNI "SOURCE" (ZGORAJ DESNO), DA NALOŽIŠ FUNKCIJE V R
#################   1. 12. 2014 ######################################
#Časovno okno
get.window=function(data,width){
  #data = genind object
  #width = širina časovnega okna
  #ne obvlada NA v letu... odstrani.
  require(adegenet)
  Hobs=NULL
  Hexp=NULL
  minyear = NULL
  maxyear = NULL
  Fi = NULL
  A = NULL
  medA=NULL
  Aall=NULL
  medyear = NULL
  meanyear=NULL
  firstsamp=NULL
  lastsamp=NULL
  hwtest=NULL
  arlequin.data=NULL
  
  pb <- txtProgressBar(min = 1, max = length(indNames(data))-width+1, style = 3)
  
  for (i in 1:(length(indNames(data))-width+1)){
    data.sub=data[i:(i+width-1)]
    #print(cat("Okno:",as.character(i)))
    setTxtProgressBar(pb, i)
    
    other=data.sub$other
    data.sub=genind2df(data.sub,sep=";")
    
    sub.temp=data.sub
    sub.temp$pop = rep(i,length(sub.temp$pop))
    sub.temp$pop=as.character(sub.temp$pop)
    
    poprow=sub.temp[1,]
    row.names(poprow)="pop"
    poprow[1,1:ncol(poprow)]=rep("",ncol(poprow))
    
    
    arlequin.data=rbind(arlequin.data,poprow)
    arlequin.data=rbind(arlequin.data,sub.temp)
    rm(poprow)
    
    data.sub=df2genind(data.sub[,2:ncol(data.sub)],sep=";")
    data.sub$other=other
    
    params = summary(data.sub,verbose=F)
    #hwtest=rbind(hwtest, HWE.test.genind(data.sub,res="matrix",permut=T,nsim=10000))
    
    minyear=c(minyear,data.sub$other$Year[1])
    maxyear=c(maxyear,data.sub$other$Year[width])
    medyear=c(medyear,median(data.sub$other$Year))#median year of the subsample
    meanyear=c(meanyear,mean(data.sub$other$Year))#mean year of the subsample
    firstsamp=c(firstsamp,indNames(data.sub)[1])
    lastsamp=c(lastsamp,indNames(data.sub)[nInd(data.sub)])
    Hobs = c(Hobs,mean(params$Hobs))
    Hexp = c(Hexp,mean(params$Hexp))
    Fi = c(Fi, mean(sapply(inbreeding(data.sub,N=100),mean)))
    #A = c(A,mean(params$loc.nall))
  #   medA = c(medA,median(params$loc.nall))
    #Aall=rbind(Aall,params$loc.nall)
    
    
  }#for
  close(pb)
  #return(list(data.frame(minyear,maxyear,medyear,meanyear,A,medA,Hobs,Hexp,Fi,firstsamp,lastsamp),Aall,hwtest,arlequin.data))
  return(list(data.frame(minyear,maxyear,medyear,meanyear,Hobs,Hexp,Fi,firstsamp,lastsamp),arlequin.data))
  
}#function

########################################33

DriftGraf = function(t,R,Ne){
  #Simulira genetski zdrs (drift) po Wright-Fisherjevem modelu
  #t = čas
  #R = število simulacij
  #Ne = efektivna velikost populacije
  
  p=0.5 #frekvenca alela
  freq=as.numeric();
  for( i in 1:t ){
    A1=rbinom(1,2*Ne,p)
    p=A1/(Ne*2);
    freq[length(freq)+1]=p;
  }
  plot(freq,type="l",ylim=c(0,1),col=3,xlab="t",ylab=expression(p(A[1])))
  for(u in 1:R){
    freq1=as.numeric();
    p=0.5
    for( j in 1:t ){
      A1=rbinom(1,2*Ne,p)
      p=A1/(Ne*2);
      freq1[length(freq1)+1]=p;
    }
    random=sample(1:1000,1,replace=F)
    randomcolor=colors()[random] 
    lines(freq1,type="l",col=(randomcolor))
  }
}

#Drift_graph(100,50,250)
##############################

HFGraf = function(Ne,t,R=100){
  #Simulira genetski zdrs (drift) po Wright-Fisherjevem modelu
  #t = čas
  #R = število lokusov
  #Ne = efektivna velikost populacije
  
  ps=0.5 #frekvenca alela
  
  Het=NULL
  
  for( i in 1:R ){
    inHet = NULL
    p=ps
    for(j in 1:t) {
      A1=rbinom(1,2*Ne,p)
      p=A1/(Ne*2)
      He=2*p-2*p^2
      inHet=c(inHet,He)
    }
    Het=rbind(Het,inHet)
  }
  plHet=apply(Het,2,mean)
  
  Hs = 2*ps-2*ps^2
  
  plInbr=1-(plHet/Hs)
  par(mfrow=c(1,2))
  
  plot(plHet/(2*ps*(1-ps)),type="l",ylim=c(0,1),col="blue",xlab="t",ylab= expression(H [t]/H [0]), main = paste("delež preostale heterozigotnosti\nv času t (Ht/H0), Ne=",Ne),cex.main=0.8)
  plot(plInbr,type="l",ylim=c(0,1),col="blue",xlab="t",ylab= "F", main = paste("rast inbreedinga, Ne=",Ne),cex.main=0.8)
  abline(h=0.2,col="red",lty="dashed")
  par(mfrow=c(1,1))
  #DriftGraf(t,R,Ne)
  
}


####################################################



HFDataSim = function(Ne,R=1000,StartYear, EndYear, G, Hs=0.5, Htot=0.5){
  #Simulira genetski zdrs (drift) po Wright-Fisherjevem modelu
  #glede na začetno stanje heterozigotnosti.
  #t = čas
  #R = število lokusov
  #Ne = efektivna velikost populacije
  #StartYear = začetno leto simulacije
  #EndYear = končno leto simulacije
  #Hs = začetna heterozigotnost
  #G = generacijski čas.
  #Htot = skupna heterozigotnost (Za Wright-ov Fit)
  
  ps = (2-sqrt(4-8*Hs))/4 #frekvenca alela
  t=floor((EndYear-StartYear)/G)
  
  Het=NULL
  
  for( i in 1:R ){
    inHet = NULL
    p=ps
    for(j in 1:t) {
      A1=rbinom(1,2*Ne,p)
      p=A1/(Ne*2)
      He=2*p-2*p^2
      inHet=c(inHet,He)
    }
    Het=rbind(Het,inHet)
  }
  plHet=c(Hs,apply(Het,2,mean))
  
  #Hs = 2*ps-2*ps^2
  
  plInbr=1-(plHet/Htot)
  year = seq(from=StartYear,by=G,length.out=t+1)
  return(data.frame(plHet,plInbr,year))
  #par(mfrow=c(1,2))
  
  #plot(plHet,type="l",ylim=c(0,1),col="blue",xlab="t",ylab= expression(H [t]), main = paste("heterozigotnost\nv času t (Ht), Ne=",Ne),cex.main=0.8)
  #plot(plInbr,type="l",ylim=c(0,1),col="blue",xlab="t",ylab= "F", main = paste("rast inbreedinga, Ne=",Ne),cex.main=0.8)
  
  #abline(h=0.2,col="red",lty="dashed")
  par(mfrow=c(1,1))
  #DriftGraf(t,R,Ne)
  
}

##########################
HFDataSimMulti= function(Ne,R=1000,StartYear, EndYear, G, Hs=0.5, Htot=0.5,Nsim=100){
#Wrapper za HFDataSim za več iteracij.
#Nsim = število simulacij
  out = NULL  
  for (i in 1:Nsim){
      oneSim = HFDataSim(Ne=Ne, R=R,StartYear=StartYear, EndYear=EndYear, G=G, Hs = Hs, Htot = Htot)
      oneSim = cbind(oneSim, iter=rep(i,nrow(oneSim)))
      out=rbind(out,oneSim)
  }
  return(out)
}

    
########################


HoPadecSim = HFDataSim(Ne=20, R=40,StartYear = 1990, EndYear = 2010, G=3, Hs = 0.495, Htot = 0.594)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}