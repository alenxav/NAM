# This function allows you to export files into the fastPHASE format, perform imputations, and import back to R in the proper format

fast=function(gen,chr=NULL,act="EXP"){
  
  if((act=="EXP"||act=="IMP")==FALSE){stop(paste("'act' must to be set as either EXP or IMP"))}
  if(is.null(chr)){chr=1}
  
  if(act=="EXP"){
    #Generating the fastPHASE files
    outfast=function(gen,chr){
      if(is.null(chr)){chr=1}
      input=function(gen){
        id=c(1:nrow(gen))
        hash=data.matrix(paste("#",id,sep=""))
        gen1=data.frame(gen)
        gen1[gen1==2]=1
        gen1[is.na(gen1)]="?"
        gen1=data.matrix(apply(gen1,1,paste,collapse=""))
        gen2=data.frame(gen)
        gen2[gen2==1]=0
        gen2[gen2==2]=1
        gen2[is.na(gen2)]="?"
        gen2=data.matrix(apply(gen2,1,paste,collapse=""))
        gen3=c()
        for(i in id){
          gen3=rbind(gen3,
                     rbind(data.matrix(hash[i,1]),data.matrix(gen1[i,1]),data.matrix(gen2[i,1])))}
        gen3=data.frame(gen3)
        P=paste(c("P",paste(id,collapse=" ")), collapse=" ")
        S=paste(rep("S",ncol(gen)),collapse="")
        intro=data.frame(rbind(nrow(gen),ncol(gen),P,S));colnames(intro)="gen3"
        final=rbind(intro,gen3)
        return(final)}
      CHR=c();for(i in 1:length(chr)){a=rep(i,chr[i]);CHR=c(CHR,a)};rm(a)
      for(i in 1:length(chr)){
        a=gen[,which(CHR==i)]
        b=paste("CHR",if(i<10){as.raw(i)}else{i},sep="")
        c=input(a)
        write.table(c,b,quote=F,row.names=F,col.names=F)}}
    outfast(gen,chr)
    # Script to run fastPHASE: ./fastPHASE -T10 -oCHR01 -H-4 -Z -K8 CHR##
    b=c();for(i in 1:length(chr)){
      a=data.frame(paste("./fastPHASE -T10 -oCHR",if(i<10){as.raw(i)}else{i}," -H-4 -Z -K8 CHR",if(i<10){as.raw(i)}else{i},sep=""));colnames(a)="col"
      b=rbind(b,a)}
    write.table(b,"script",quote=F,row.names=F,col.names=F)
  }
  
  if(act=="IMP"){
    # Importing files back
    infast=function(chr){
      Dir=getwd();gen=c();
      for(i in 1:length(chr)){
        File=paste(Dir,"/CHR",if(i<10){as.raw(i)}else{i},"_genotypes.out",sep="")
        out <- read.table(File, quote="\"")
        splitting=function(XXX){a1=as.character(XXX);a2=strsplit(a1,"");a3=unlist(a2);a4=as.vector(as.numeric(a3));return(a4)}
        out1=t(apply(out,1,splitting));
        #out2=c()
        #for(j in 1:nrow(out)){out2=rbind(out2,out1[[j]])}
        out2=out1 
        even=which((1:nrow(out2)%%2)==0)
        odd=which((1:nrow(out2)%%2)!=0)
        out3=data.matrix(out2[even,]+out2[odd,])
        colnames(out3)=rep(i,ncol(out3))
        gen=cbind(gen,out3)}
    return(gen)}
    gen=infast(chr);return(gen)}
  
}
  
