

##read otu table
otutab <- read.delim('otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

otutab[is.na(otutab)]<-0
##keep 12 otus
counts<-rowSums(otutab>0)

comm<-t(otutab)
sp.ra<-colMeans(comm)/17188   #relative abundance of each species

###### there are two choices to get the correlation matrix #######
###### choice 2: directely read in correlation matrix downloaded from MENAP. MENAP downloaded correlation matrix is an upper triangle. Need to make it into a symetric matrix.
cormatrix2<-read.table("net.txt",header = T,row.names=1,sep="\t")
###### end of the two choices of correlation matrix #######
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.
network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched

rand.remov2.once<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    
    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0) 
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]
      
      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }
      
      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }
      
    }
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}
rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}
node.attri<-read.csv("keystone.csv",header = T, sep=",")
currentdat = NULL
for (i in 1:100) {
  set.seed(5) # 设置随机种子
  node.attri<-read.csv("keystone.csv",header = T, sep=",")
  node.attri <- node.attri[sample(nrow(node.attri),17), ] 
  module.hub<-as.character(node.attri$v_name)

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=length(module.hub)),
                 time=rep(i,,2*length(module.hub)))
currentdat = rbind(currentdat,dat1)
write.csv(currentdat, "removal_keystone.csv")
}
