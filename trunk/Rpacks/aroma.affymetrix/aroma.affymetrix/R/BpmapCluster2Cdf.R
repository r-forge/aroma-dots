
BpmapCluster2Cdf<-function(filename,cdfName,nProbes=30,gapDist=3000,rows=NULL,cols=NULL,field="fullname",verbose=10,stringRemove="Hs:Sun.Nov.19.17:25:02.2006;") {
  require(affxparser)
  require(R.utils)

  BpmapUnit2df<-function(u) {
    o<-order(u[["startpos"]])
    mmx<-u[["mmx"]]
    mmy<-u[["mmy"]]
    mmx<-if( all(mmx==0) | is.null(mmx)) 0 else mmx
    mmy<-if( all(mmy==0) | is.null(mmy)) 0 else mmy
    data.frame(seqname=u$seqInfo[[field]],u[c("pmx","pmy")],mmx=mmx,mmy=mmy,u[c("probeseq","strand","startpos","matchscore")],stringsAsFactors=FALSE)[o,]
  }
  if (verbose) cat("Reading BPMAP file",filename,"... ")
  bpmaplist<-readBpmap(filename,readMatchScore=TRUE)
  if (verbose) cat("extracting X/Y locations ... ")
  bpmapdflist<-lapply(bpmaplist,BpmapUnit2df)
  if (verbose) cat("Done.\n")

  rm(bpmaplist); gc()  # save some space

  if(is.null(rows)) {
    rows<-max(sapply(bpmapdflist,FUN=function(u) max(c(u$mmx,u$pmx))))
    cat("NB: 'rows' of CDF are being set as ",rows,". If this is not correct, stop now and specify 'rows' argument.",sep="")
  }
  if(is.null(cols)) {
    cols<-max(sapply(bpmapdflist,FUN=function(u) max(c(u$mmy,u$pmy))))
    cat("NB: 'cols' of CDF are being set as ",cols,". If this is not correct, stop now and specify 'col' argument.",sep="")
  }
    
  e<-vector("list",1)
  startps<-l<-ll<-vector("list",50000)
  if (verbose) cat("Creating structure for",length(bpmapdflist),"units (dot=250):\n")
  count<-0
  nm<-rep(NA,length(l))
  for(i in  1:length(bpmapdflist)) {
    ch<-gsub(stringRemove,"",bpmapdflist[[i]]$seqname[1])
    if(all(bpmapdflist[[i]]$mmx==0) & all(bpmapdflist[[i]]$startpos>0)) {
	  sp<-bpmapdflist[[i]]$startpos
      d<-diff(sp)
	  w<-which(d>gapDist)
	  ends<-c(w,nrow(bpmapdflist[[i]]))
	  starts<-c(1,w+1)
	  k<-which( (ends-starts)>nProbes )
	  if(verbose) cat(length(k)," ROIs for ",names(bpmapdflist)[i]," (dot=250):",sep="")
      for(j in 1:length(k)) {
	    w<-starts[k[j]]:ends[k[j]]
        np<-length(w)
		# this assumes only PM probes -- is this true for all bpmaps?
        e[[1]]<-list(x=bpmapdflist[[i]]$pmx[w],y=bpmapdflist[[i]]$pmy[w],pbase=rep("A",np),tbase=rep("T",np),
		             atom=0:(np-1),indexpos=0:(np-1),groupdirection="sense",natoms=np,ncellsperatom=1)
        names(e)<-paste(ch,"FROM",sp[starts[k[j]]],"TO",sp[ends[k[j]]],sep="")
        na<-sum(unlist(sapply(e,FUN=function(u) u$natoms)))
        nc<-sum(unlist(sapply(e,FUN=function(u) u$natoms*u$ncellsperatom)))
		count<-count+1
        l[[count]]<-list(unittype=1,unitdirection=1,groups=e,natoms=na,ncells=nc,ncellsperatom=nc/na,unitnumber=i)
		startps[[count]]<-sp[w]
		nm[count]<-names(e)
		if (verbose) { if(count %% 250==0) cat(".") }
      }
	  if(verbose) cat("\n")
	} else { # keep all probes
        np<-nrow(bpmapdflist[[i]])
	    if(verbose) cat("Skipping all ",np," probes for ",names(bpmapdflist)[i],".\n",sep="")
		next
	}
  }
  if (verbose) cat("\n")
  l<-l[1:count]
  startps<-startps[1:count]
  names(startps)<-names(l)<-nm[1:count]
  saveObject(startps,file=paste(cdfName,"pps",sep="."))
  hdr<-list(probesets=length(l),qcprobesets=0,reference="",chiptype=cdfName,filename=paste(cdfName,"cdf",sep="."),
           nqcunits=0,nunits=length(l),rows=rows,cols=cols,refseq="",nrows=rows,ncols=cols)
  writeCdf(hdr$filename, cdfheader=hdr, cdf=l, cdfqc=NULL, overwrite=TRUE, verbose=verbose)
  invisible(list(cdfList=l,cdfHeader=hdr))
}







