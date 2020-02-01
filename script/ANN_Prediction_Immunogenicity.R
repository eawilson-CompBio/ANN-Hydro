a<-commandArgs(trailingOnly = T)
##A Computational modelling framework for predicting immunogenic T-cell epitopes

##check that required packages are installed
req.packages<-c("neuralnet","stringr")

new.packages<-req.packages[-which(req.packages%in%installed.packages()[,"Package"])]
if(length(new.packages)>0){
    sapply(new.packages,function(package){
        install.packages(package)
    })
}

library(neuralnet)
library(stringr)
cat("reading in data\n\n")
data_training = read.csv(a[1],)
data_testing = read.csv(a[2]) 
load("./data/hydro_index.R")




convert_to_numeric<-function(x){
    data.frame(epitope=x,apply(str_split(x,pattern = "",simplify = T),2,function(w){
        hydro[w]
    }),stringsAsFactors = F)
}

training_set<-data.frame(convert_to_numeric(data_training[,1]),Immunogenicity=data_training[,2],stringsAsFactors = F)

testing_set<-convert_to_numeric(data_testing[,1])
                                                            

cat("training neural network\n\n")
predicted_immunogenicity<-sapply(1:a[3],function(count){
    if(count%%10==0){
        per<-round((count/as.numeric(a[3]))*100,2)
        cat(paste(per,"% done!\n",sep=""))
    }
nn <- neuralnet(Immunogenicity ~. , data=training_set[,-1], hidden=c(3), linear.output=FALSE, rep=1)
compute(nn, testing_set)$net.result
})

pred_out<-data.frame(peptide=testing_set[,1],predicted_immunogenicity=apply(predicted_immunogenicity,1,mean))

if(ncol(data_testing)==2){
    cat("\nuser provided immunogenicity data for testing peptides. calculating accuracy metrics\n")
    pred_out$immunogenicity<-data_testing[,2]
    if(is.na(a[5])){
        cat("No threshold provided for predictions. setting to 0.5\n")
        a[5]<-.05
    }
    tab<-table(predicted_immuno=as.numeric(pred_out$predicted_immunogenicity>a[5]),measured_immuno=pred_out$immunogenicity)
    ACC<-(tab[4]+tab[1])/sum(tab)
    PER<-(tab[4])/sum(tab[4]+tab[2])
    REC<-(tab[4])/sum(tab[4]+tab[3])
    FPR<-(tab[2])/sum(tab[1:2])
    cat(paste(paste("Accuracy:",ACC),paste("Precision:",PER),paste("Recall",REC),paste("False Positive Rate: ",FPR,"\n"),sep="\n"))
    write.csv(pred_out,file=paste("./Predictions/",a[4],sep = ""),row.names=F)
    
}else{
    write.csv(pred_out,file=paste("./Predictions/",a[4],sep = ""),row.names=F)
}    



