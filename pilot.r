#data_positives = read.table("Self CTL positive hydro.csv",header=T,sep=",")
#data_positives = read.table("All_negative_unique.csv",header=T,sep=",")
#data_positives = read.table("All CTL.csv",header=T,sep=",")
#data_negatives = read.table("All_negative_unique.csv",header=T,sep=",")
#data_negatives = read.table("Self-eluted-part1&2-JP.csv",header=T,sep=",")
#data_negatives = read.table("Human_Self_Eluted.csv",header=T,sep=",") 
#data_positives = read.table("All CTL minus A2.csv",header=T,sep=",")
#data_negatives = read.table("All Self Eluted minus A2.csv",header=T,sep=",")
#data_positives = read.table("SelfAntigens positive complete data_unique.csv",header=T,sep=",")
#data_negatives = read.table("All Self Eluted minus A2.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A2 CTL.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A2 CTL positive minus linear invitrostim.csv",header=T,sep=",")
#data_positives = read.table("All pathogens H2Kb CTL.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted H2-Kb.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A2 CTL.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A2 CTL negative.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A3 CTL.csv",header=T,sep=",")
#data_negatives = read.table("All self eluted-A3-IEDB-syfpethi hydro.csv",header=T,sep=",")
#data_positives = read.table("All CTL A2.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted.csv",header=T,sep=",")
#data_positives = read.table("All CTL A3.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted A3.csv",header=T,sep=",")
#data_positives = read.table("All pathogens H2Db CTL.csv",header=T,sep=",")
#data_positives = read.table("All pathogens H2Db CTL.csv",header=T,sep=",")
#data_positives = read.table("All CTL H2-Db 9mers.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted H2-Db.csv",header=T,sep=",")
#data_positives = read.table("All CTL H2-Kb.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted H2-Kb.csv",header=T,sep=",")
#data_positives = read.table("H2Kb eluted-test hydro 8mers.csv",header=T,sep=",")
#data_negatives = read.table("H2Kb CTL-test hydro 8mers.csv",header=T,sep=",")
#data_negatives = read.table("All Self Eluted minus A2.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted A2.csv",header=T,sep=",")
#data_positives = read.table("All pathogens CTL.csv",header=T,sep=",")
#data_positives = read.table("All pathogens CTL Ag-type organism.csv",header=T,sep=",")
#data_positives = read.table("All pathogens CTL positive-noA2-only HLAs.csv",header=T,sep=",")
#data_negatives = read.table("All Self Eluted-noA2-only HLAs.csv",header=T,sep=",")
#data_positives = read.table("All CTL H2-Db 9mers.csv",header=T,sep=",y1")
#data_negatives = read.table("All Self-Eluted H2-Db 9mers.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A3 CTL.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted A3.csv",header=T,sep=",")
#data_positives = read.table("human-mouse self CTL 9mers.csv",header=T,sep=",")
#data_positives = read.table("human-mouse self CTL A2.csv",header=T,sep=",")
#data_positives = read.table("A2 full-1.csv",header=T,sep=",")
#data_positives = read.table("All pathogens H2Db CTL pos.csv",header=T,sep=",")
#data_negatives = read.table("All pathogens H2Db CTL neg.csv",header=T,sep=",")
#data_negatives = read.table("All Self-Eluted H2-Db.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A2 CTL pos.csv",header=T,sep=",")
#data_positives = read.table("All pathogens A2 CTL neg.csv",header=T,sep=",")
#data_positives = read.table("All pathogens H2Kb CTL pos.csv",header=T,sep=",")
data_positives = read.table("All pathogens H2Kb CTL neg.csv",header=T,sep=",")
data_negatives = read.table("All Self-Eluted H2-Kb.csv",header=T,sep=",")


HydrophoTable = read.table("HydrophobicityTable.csv",header=T,sep=",")   


lengths_peptides_positives=data_positives$Epitope.End-data_positives$Epitope.Start
lengths_peptides_positives=lengths_peptides_positives+1

lengths_peptides_negatives= data_negatives$Epitope.End-data_negatives$Epitope.Start
lengths_peptides_negatives=lengths_peptides_negatives+1


mer=8

index_peptides_positives=which(lengths_peptides_positives==mer)
peptides_positives=data_positives$Epitope.Linear.Sequence[index_peptides_positives]

index_peptides_negatives=which(lengths_peptides_negatives==mer)
peptides_negatives=data_negatives$Epitope.Linear.Sequence[index_peptides_negatives]

cc_positives=vector()                                       

for (i in 1:length(peptides_positives)){

chars_positives=strsplit(gsub("([[:alnum:]]{1})", "\\1 ", peptides_positives[i]), " ")[[1]]

for (j in 1:length(chars_positives)){

index_amino_positives=which(chars_positives[j]==HydrophoTable$Amino.acid)
chars_positives[j]=HydrophoTable$Hydrophobicity[index_amino_positives]

}

cc_positives=c(cc_positives,chars_positives)
cc_positives=as.numeric(cc_positives)

}

matrix_positives=matrix(cc_positives,nrow=mer,ncol=length(peptides_positives))
matrix_positives=t(matrix_positives)

#write.table(matrix_positives,file ="matrix_positives_file.txt", sep = "\t", col.names = NA, qmethod = "double")


cc_negatives=vector()

for (i in 1:length(peptides_negatives)){
  
chars_negatives=strsplit(gsub("([[:alnum:]]{1})", "\\1 ", peptides_negatives[i]), " ")[[1]]
  
for (j in 1:length(chars_negatives)){

index_amino_negatives=which(chars_negatives[j]==HydrophoTable$Amino.acid)
chars_negatives[j]=HydrophoTable$Hydrophobicity[index_amino_negatives]

}

cc_negatives=c(cc_negatives,chars_negatives)
cc_negatives=as.numeric(cc_negatives)

}

matrix_negatives=matrix(cc_negatives,nrow=mer,ncol=length(peptides_negatives))
matrix_negatives=t(matrix_negatives)


#ErrorBarExample <- function(){

#add.error.bars <- function(X,Y,SE,w,col=1){                                                    
#X0 = X; Y0 = (Y-(SE*1.96)); X1 =X; Y1 = (Y+(SE*1.96));
#arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col);
#}

y1=apply(matrix_positives,2,mean)

y2=apply(matrix_negatives,2,mean)

#y1.sd=apply(matrix_positives,2,sd)
#y2.sd=apply(matrix_negatives,2,sd)

plot(1:mer,y1,type="o",xlim=c(1,8),ylim=c(-3.5,4),xlab=paste("Residue position"),ylab="Mean hydrophobicity",main=paste("H2Kb 8mers"),font.main=1,cex.lab=1.5,cex.main=1.5,lwd=1.3,cex=1.3,pch=15, bty='l',las=1,xaxt="n",cex.axis=1.3)
#plot(1:mer,y1,type="o",xlim=c(1,8),ylim=c(-3,4),xlab=paste("Residue position"),ylab="Mean hydrophobicity",main=bquote(MHC~H2K^b ~ 8*mers),font.main=1,cex.lab=1.5,cex.main=1.5,lwd=1.3,cex=1.3,pch=15, bty='l',las=1,xaxt="n",cex.axis=1.3)
#add.error.bars(1:mer,y2,y2.sd/sqrt(1000),0.1,col="black");
par(new=TRUE)
plot(1:mer,y2,col="black",type="o",xlim=c(1,8),ylim=c(-3.5,4),ann=FALSE,cex.lab=1.5,lwd=1.5,cex=1.3,pch=22,bty='l',las=1, xaxt="n",cex.axis=1.3)
axis(1, at=1:10, labels=c(1,2,3,4,5,6,7,8,9,10),cex.axis=1.3)
#axis(1, at=1:8, labels=c(1,2,3,4,5,6,7,8),cex.axis=1.3)
#add.error.bars(1:mer,y1,y1.sd/sqrt(1000),0.1,col="black");
legend("topleft",legend=c("T cell negative","Self-eluted"),pch=c(15,22),bty="n",lwd=1.5,cex=1.1)
#legend("topleft",legend=c("Immunogenic (n=288)","Non-immunogenic (n=77)"),pch=c(15,22),bty="n",lwd=1.2,cex=1.1)

#}                                                                            

#ErrorBarExample();




#x1=apply(matrix_negatives,2,median)
#plot(apply(matrix_positives,2,median),type="b",ylim=c(min(x1)-0.5,4.5),xlab="Residue position",ylab="Mean hydrophobicity value",lwd=1.5)
#par(new=TRUE)
#plot(apply(matrix_negatives,2,median),col="red",type="b",ylim=c(min(x1)-0.5,4.5),xlab="Residue position",ylab="Mean hydrophobicity value",lwd=1.5)



#legend("topleft",c("CTL positive peptides","Self-eluted peptides"),col=c("black","red"),bty="n",lwd=5,cex=1)


#p.value.wilcox=wilcox.test(matrix_positives[,4],matrix_negatives[,4],correct=FALSE)
#p.value.t=t.test(matrix_positives[,i],matrix_negatives[,i])
