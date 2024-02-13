# functions_cpp.R

# Functions for detection calling from HG-U219 array
# The detection calling function needs to redesigned for readability and maintainability.
# Ideally create an R package with the following code.
# Hint: 1) Package the algib wicoxon hypothesis testing library (vectorize it).
#       2) Define a dbag2call function around it.


# Defining the R interfaces (dabg2call() & mas5gccall()) for the respective
# - HG-U219 detection calling algorithm implemented in C++.
# Further a higher level function, getu219dcal(), that depends on the above mentioned
# - R interfaces is also defined.
# See code/cpp/*


dyn.load("code/cpp/dabg2call.so") # dabg2call() depends on this shared object
dyn.load("code/cpp/mas5gccall.so") # mas5gccall() depends on this shared object


dabg2call<-function(aBatch,probeGcList,refProbeSetType="medGC",fixedRefGC=NULL){

  # R function to call the relevant dabg2 detection calling algorithm for
  # - the PM-only (Perfect match only) HG-U219 array implemented in C++.
  # Based on wilcoxon rank-sum test
  # Depends on shared object dabg2call.so
  # See Jose et.al. 2018, https://doi.org/10.1371/journal.pone.0203346


  # aBatch > AffyBatch object
  # probeGcList   > probeset named list with each entry corresponds to a vector of respective probe GC count.
  # refProbeSetType   > minGC/medGC/maxGC/relGC/allGC/fixGC
  # fixedRefGC > custom specified reference

  if(refProbeSetType=="fixGC"&&is.null(fixedRefGC)){
    stop("getDABG2: fixedRefGC unspecified.")
  }
  data=intensity(aBatch)
  rm(aBatch)
  gc()

  refIndex=grepl(x=names(probeGcList),pattern="GC",fixed=T)
  gcList=lapply(probeGcList[!refIndex],function(x){x$gc})
  indexList=lapply(probeGcList[!refIndex],function(x){x$indexIntensity})
  indexList_ref=lapply(probeGcList[refIndex],function(x){x$indexIntensity})
  indexList_ref=indexList_ref[sort(names(indexList_ref))]
  # sorted indexList_ref: index1->gc3; 2->gc4; ... 23->gc25.
  # index=gc-2

  # setting indexList_ref2 according to refProbeSetType
  if(refProbeSetType=="medGC"){
    indexList_ref2=lapply(gcList,
                          function(x,indexList_ref){
                            gc=as.integer( floor(median(x)) )
                            return(indexList_ref[[gc-2]])
                          },indexList_ref)
    names(indexList_ref2)=names(gcList)

  }else if(refProbeSetType=="maxGC"){
    indexList_ref2=lapply(gcList,
                          function(x,indexList_ref){
                            gc=as.integer(max(x))
                            return(indexList_ref[[gc-2]])
                          },indexList_ref)
    names(indexList_ref2)=names(gcList)

  }else if(refProbeSetType=="relGC"){
    indexList_ref2=lapply(gcList,
                          function(x,indexList_ref){
                            gc=unique(x)
                            index=lapply(gc,
                                         function(xx,indexList_ref){
                                           indexList_ref[[xx-2]]
                                         },indexList_ref)
                            return(unlist(index))
                          },indexList_ref)
    names(indexList_ref2)=names(gcList)

  }else if(refProbeSetType=="allGC"){
    indexList_ref2=list(allGC=unlist(indexList_ref))

  }else if(refProbeSetType=="fixGC"){
    indexList_ref2=lapply(fixedRefGC,
                          function(x,indexList_ref){
                            indexList_ref[[x-2]]
                          },indexList_ref)
    indexList_ref2=list(fixGc=unlist(indexList_ref2))
  }

  dCal=.Call("dabg2call",data,indexList,indexList_ref2,refProbeSetType)
  colnames(dCal)=colnames(data)
  rownames(dCal)=names(indexList)
  gc()
  return(dCal)
}



# Test/Back-up code for dabg2call()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# system.time(outMedGc2<-dabg2call(aBatch=bgCor.rma.u219.ff.sen,
#                                  probeGcList=hgu219.probeInfo.gcIndexList,
#                                  refProbeSetType="medGC",
#                                  fixedRefGC=NULL))
# system.time(outRelGc2<-dabg2call(aBatch=bgCor.rma.u219.ff.sen,
#                                  probeGcList=hgu219.probeInfo.gcIndexList,
#                                  refProbeSetType="relGC",
#                                  fixedRefGC=NULL))
# system.time(outMaxGc2<-dabg2call(aBatch=bgCor.rma.u219.ff.sen,
#                                  probeGcList=hgu219.probeInfo.gcIndexList,
#                                  refProbeSetType="maxGC",
#                                  fixedRefGC=NULL))
# system.time(outFixGc2<-dabg2call(aBatch=bgCor.rma.u219.ff.sen,
#                                  probeGcList=hgu219.probeInfo.gcIndexList,
#                                  refProbeSetType="fixGC",
#                                  fixedRefGC=c(12,13)))
# system.time(outAllGc2<-dabg2call(aBatch=bgCor.rma.u219.ff.sen,
#                                  probeGcList=hgu219.probeInfo.gcIndexList,
#                                  refProbeSetType="allGC",
#                                  fixedRefGC=NULL))

# # validation
# data=intensity(bgCor.rma.u219.ff.sen)
# samp1=data[,1]
# samp2=data[,2]
# #medGc samp1=16; samp2=11; samp3=11;
# # "AFFX-Nonspecific-GC16_at"; "AFFX-Nonspecific-GC11_at"
# wilcox.test(x=samp1[hgu219.probeInfo.gcIndexList[[1]]$indexIntensity],
#             y=samp1[hgu219.probeInfo.gcIndexList[["AFFX-Nonspecific-GC16_at"]]$indexIntensity],
#             alternative="greater")#p=0.1866; 0.8137
# wilcox.test(x=samp1[hgu219.probeInfo.gcIndexList[[2]]$indexIntensity],
#             y=samp1[hgu219.probeInfo.gcIndexList[["AFFX-Nonspecific-GC11_at"]]$indexIntensity],
#             alternative="greater")#p=0.9708; 0.02924
# wilcox.test(x=samp1[hgu219.probeInfo.gcIndexList[[3]]$indexIntensity],
#             y=samp1[hgu219.probeInfo.gcIndexList[["AFFX-Nonspecific-GC11_at"]]$indexIntensity],
#             alternative="greater")#p=0.9995; 0.0005379
#
# log2(samp1[hgu219.probeInfo.gcIndexList[[3]]$indexIntensity[1:10]])
# log2(samp1[hgu219.probeInfo.gcIndexList[["AFFX-Nonspecific-GC11_at"]]$indexIntensity[1:10]])
#
# #allGC
# index=sapply(antiGenomic,function(x,id){which(x==id)},
#                 id=names(hgu219.probeInfo.gcIndexList))
# #samp1
# refIndex=lapply(index,function(x,id){id[[x]]$indexIntensity},
#            id=hgu219.probeInfo.gcIndexList)
# refIndex=unlist(refIndex)
# curSamp=data[,1,drop=T]
#
# wilcox.test(x=curSamp[hgu219.probeInfo.gcIndexList[[1]]$indexIntensity],
#             y=curSamp[refIndex],alternative="less")#p=0.2566
# wilcox.test(x=curSamp[hgu219.probeInfo.gcIndexList[[2]]$indexIntensity],
#             y=curSamp[refIndex],alternative="less")#p=0.768
# wilcox.test(x=curSamp[hgu219.probeInfo.gcIndexList[[3]]$indexIntensity],
#             y=curSamp[refIndex],alternative="less")#p=0.2987





mas5gccall<-function(aBatch,probeGcList,tau=0.015){

  # R function to call the mas5gc detection calling algorithm for
  # - the PM-only (Perfect match only) HG-U219 array implemented in C++
  # Based on wilcoxon signed-rank test
  # Depends on shared object mas5gccall.so
  # See Jose et.al. 2018, https://doi.org/10.1371/journal.pone.0203346


  # aBatch > AffyBatch object
  # probeGcList   > probeset named list with each entry corresponds to a vector of respective probe GC count.
  # tau > same as MAS5 tau

  data=intensity(aBatch)
  rm(aBatch)
  gc()

  refIndex=grepl(x=names(probeGcList),pattern="GC",fixed=T)
  gcList=lapply(probeGcList[!refIndex],function(x){x$gc})
  indexList=lapply(probeGcList[!refIndex],function(x){x$indexIntensity})
  indexList_ref=lapply(probeGcList[refIndex],function(x){x$indexIntensity})
  indexList_ref=indexList_ref[sort(names(indexList_ref))]
  # sorted indexList_ref: index1->gc3; 2->gc4; ... 23->gc25.
  # index=gc-2

  # Finding antigenomic-probeset median per sample
  dataNcMed=matrix(data=NA,ncol=ncol(data),nrow=length(indexList_ref))
  for(i in 1:ncol(dataNcMed)){
    for(j in 1:nrow(dataNcMed)){
      dataNcMed[j,i]=median(data[indexList_ref[[j]],i])
    }
  }
  colnames(dataNcMed)=colnames(data)
  rownames(dataNcMed)=names(indexList_ref)

  # #validation
  # dCal=matrix(data=NA,ncol=ncol(data),nrow=100)
  # colnames(dCal)=colnames(data)
  # rownames(dCal)=names(indexList)[1:100]
  #
  # for(i in 1:ncol(data)){
  #     curSamp=data[,i]
  #     curNcMed=dataNcMed[,i]
  #     for(j in 1:100){
  #         curIndex=indexList[[j]]
  #         curGcIndex=gcList[[j]]-2
  #         x=curSamp[curIndex]-curNcMed[curGcIndex]
  #         dCal[j,i]=wilcox.test(x=x,alternative="greater",mu=0.015,correct = TRUE)$p.value
  #     }
  #
  # }

  dCal=.Call("mas5gccall",data,dataNcMed,indexList,gcList,tau=tau)
  colnames(dCal)=colnames(data)
  rownames(dCal)=names(indexList)
  gc()
  return(dCal)
}



# Test/Back-up code for mas5gccall()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# #example
# tmp=mas5gccall(aBatch=bgCor.rma.u219.ff.sen,
#                probeGcList=hgu219.probeInfo.gcIndexList,tau=0.015)


# #testing
# tmpTest=mas5gccall(aBatch=bgCor.rma.u219.ff.sen,probeGcList=hgu219.probeInfo.gcIndexList)




getu219dcal <- function(aBatch,gcIndexList,type="medGC"){

  # A generalized higher level function for HG-U219 detection calling which
  # depends on dabg2call() and mas5gccall().

  # if(!grepl(cdfName(aBatch),pattern="hgu219cdf",fixed=T))
  #   stop("Unknown array CDF. Only hgu219cdf is supported.")
  if(!grepl(cdfName(aBatch), pattern = "u219", ignore.case = TRUE))
    stop("Unknown array CDF. Only hgu219cdf is supported.")

  print(type)

  if(type=="medGC"){
    dCal=dabg2call(aBatch=aBatch,probeGcList=gcIndexList,
                   refProbeSetType="medGC",
                   fixedRefGC=NULL)
  }else if(type=="relGC"){
    dCal=dabg2call(aBatch=aBatch,probeGcList=gcIndexList,
                   refProbeSetType="relGC",
                   fixedRefGC=NULL)
  }else if(type=="maxGC"){
    dCal=dabg2call(aBatch=aBatch,probeGcList=gcIndexList,
                   refProbeSetType="maxGC",
                   fixedRefGC=NULL)
  }else if(type=="allGC"){
    dCal=dabg2call(aBatch=aBatch,probeGcList=gcIndexList,
                   refProbeSetType="allGC",
                   fixedRefGC=NULL)
  }else if(type=="fixGC"){
    dCal=dabg2call(aBatch=aBatch,probeGcList=gcIndexList,
                   refProbeSetType="fixGC",
                   fixedRefGC=c(11,12,13))
  }else if(type=="mas5GC"){
    dCal=mas5gccall(aBatch=aBatch,probeGcList=gcIndexList,tau=0.015)
  }

  return(dCal)
}


