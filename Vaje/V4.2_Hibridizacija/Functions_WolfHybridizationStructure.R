WriteGenePop = function (genotypes1, genotypes2, path, note){
  require(adegenet)
#writes geneind as genepop file
  write.table(note,path,quote=F,col.names=F,row.names=F,sep="")
  write.table(paste(c(locNames(genotypes1),"POP"),sep="\n"),path,col.names=F,row.names=F,sep="",quote=F,append=T)
  outgenotypes=genind2df(genotypes1,sep="")
  outgenotypes[is.na(outgenotypes)]="000000"
  names=indNames(genotypes1)
  names=paste(names,",",sep="")
  outgenotypes = cbind(names,outgenotypes[2:ncol(outgenotypes)])
  
  write.table(outgenotypes,path,sep=" ",row.names=F,col.names=F,append=T,quote=F)

  #the other object
  cat("POP\n",file=path,append=T)
  outgenotypes=genind2df(genotypes2,sep="")
  outgenotypes[is.na(outgenotypes)]="000000"
  names=indNames(genotypes2)
  names=paste(names,",",sep="")
  outgenotypes = cbind(names,outgenotypes[2:ncol(outgenotypes)])
  
  write.table(outgenotypes,path,sep=" ",row.names=F,col.names=F,append=T,quote=F)
    
}

#export data for STRUCTURE
writeStructure = function (genotipi,pop,path){
  require(adegenet)
#STRUCTURE export
#genotipi = genotipi genind object
#pop = population
#path = output file
  gen=genind2df(genotipi,sep=" ")
  gen[is.na(gen)]="-9 -9"
  output = cbind(id=row.names(gen),pop=as.numeric(pop),gen[2:ncol(gen)])
  cat(locNames(genotipi),"\n",file=path)
  write.table(output,path,sep=" ",row.names=F,col.names=F,append=T,quote=F)#,eol="\r")
  outLegend=data.frame(id=1:length(levels(pop)),Pop=levels(pop))
  print(outLegend)
  return(outLegend)
}


##Plot background map
#path=file.choose()
library(sf)
borders.dinarics=st_read("GISData/SouthernEurope.shp") %>% 
  st_geometry()



#background map
plotmap=function(maintitle, xlim =c(377000,774000),ylim=c(-280000,158000),scaleloc=c(370000,-280000),scalelen=150000){
  require(sf)
  plot(borders.dinarics,xlim=xlim, ylim=ylim, col="grey90")
  scalebar(loc=scaleloc,length=scalelen)
  title(main=maintitle, cex.main=0.8)
  box()
}


#drawing a scalebar
scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])/1000
  labels <- append(labels,paste((x[5]-loc[1])/1000,unit))
  text(x[c(1,3,5)],y[4]+100,labels=labels,adj=.5,cex=division.cex)
}
