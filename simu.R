library(fda);library(dplyr)
library(Rmpfr)

source("funs_mixFPCA.R")
library(dplyr)
library(reshape2)
library(data.table)

################## Simulation study

load("first3.RData")



library(fda)
library(fdapace)
mixFPCA<-function(observed,timepoints){

	sigma = 1

	fmeanfit = findmean(observed, timepoints, minit=6,gamma=0,threshold=5e-3)

	observedmu1 = lapply(1:length(observed),function(i){
		 (eval.fd(timepoints[[i]], fmeanfit$pc_fit1))[,1]
	})
	observedmu2 = lapply(1:length(observed),function(i){
		 (eval.fd(timepoints[[i]], fmeanfit$pc_fit2))[,1]
	})

	
     

	if (mean(sapply(observedmu1,mean))<mean(sapply(observedmu2,mean))){
	tempp  = observedmu1
	observedmu1 = observedmu2
	observedmu2 = tempp
	
	mu1 = fmeanfit$pc_fit2
	mu2 = fmeanfit$pc_fit1
	fmeanfit_sigmak = rev(fmeanfit$sigmak)
	}else{
		mu1 = fmeanfit$pc_fit1
	mu2 = fmeanfit$pc_fit2
	fmeanfit_sigmak=fmeanfit$sigmak
	}
	
	res_pace<- FPCA(observed, timepoints,list(dataType='Sparse',error=TRUE, kernel='epan', verbose=TRUE,methodBwCov="GCV",methodBwMu="GCV"))


	
	timepts=1:10;
	norder=4 
	nbasis=norder+length(timepts)-2;
	spline_basis=create.bspline.basis(rangeval=c(1,10),nbasis=6,norder=4)

	coef_mat0 = coef(Data2fd(argvals = res_pace$workGrid,y=res_pace$phi,spline_basis))
	betainit=coef_mat0[,1]

	fpc1s=first_FPC(coef_mat0[,1],coef_mat0[,1],observed,observedmu1,observedmu2,timepoints,fmeanfit_sigmak,threshold=1e-3,maxit=2,gamma=0)

	previous_beta = list()
	previous_beta[[1]] = fpc1s$betak

	 
	fpc2s=third_FPC_conditional(coef_mat0[,2],coef_mat0[,2], 2, observed,observedmu1,observedmu2,timepoints,fmeanfit_sigmak, betalist=previous_beta, threshold=1e-3,maxit=50,gamma=0)

	return(list(
	mu1 = mu1,
	mu2 = mu2,

	pc11 = fpc1s$pc_fit1,
	pc12 = fpc1s$pc_fit2,

	pc21 = fpc2s$pc_fit1,
	pc22 = fpc2s$pc_fit2,

	sigma = fpc2s$sigmak,
	pik = fpc2s$pik,
	wikmat=fpc2s$wikmat
	))

}


ssize = 300


timegrid = seq(1,10,by=1)




for(it in 1:100){


clusters = lapply(1:ssize,function(i){
x= rbinom(1,1,0.4)
ifelse(x==1,1,2)
})


observed = lapply(1:ssize,function(i){

	mu1 = eval.fd(timegrid, meanfit$pc_fit1) 
	mu2 = eval.fd(timegrid, meanfit$pc_fit2)
    mu1 = mu1 - mean(mu1) - 5
	mu2 = mu2 - mean(mu2) + 5
	
	pc1add1 = eval.fd(timegrid, pc1s$pc_fit1)*rnorm(1,0,40)
	pc1add2 = eval.fd(timegrid, pc1s$pc_fit2)*rnorm(1,0,40)

	pc2add1 = eval.fd(timegrid, pc2s$pc_fit1)*rnorm(1,0,10)
	pc2add2 = eval.fd(timegrid, pc2s$pc_fit2)*rnorm(1,0,10)


	err1=rnorm(1,mean=0,sd=1)
	err2=rnorm(1,mean=0,sd=1)

	if(clusters[i]==1){
	return(as.numeric(mu1+pc1add1+pc2add1+err1))
	}else{
	return(as.numeric(mu2+pc1add2+pc2add2+err2))
	}

})

timepoints = lapply(1:ssize,function(i){
timegrid
})

if(TRUE){

res_pace<- FPCA(observed, timepoints,list(dataType='Sparse',error=TRUE, kernel='epan', verbose=TRUE,methodBwCov="GCV",methodBwMu="GCV",methodSelectK = 2))

cluster = kmeans(res_pace$xiEst,2)

estcluster = as.numeric(unlist(cluster[1]))
clustertrue = clusters%>%do.call(c,.)



pacecorrectrate  = max(length(which(3-estcluster==clustertrue))/length(estcluster),length(which(estcluster==clustertrue))/length(estcluster))

#mean
mean_correctrate=0
if(FALSE){
fmeanfit = findmean(observed, timepoints, minit=6,gamma=0,threshold=1e-5)
owik=fmeanfit$wikmat
classest = ifelse(owik[,2]>0.5,1,2)
classtrue = clusters%>%do.call(c,.)
mean_correctrate = length(which(classest==classtrue))/ssize
}

}

#inprod(pacemu-meanfit$pc_fit,pacemu-meanfit$pc_fit)

res_cfpca = mixFPCA(observed,timepoints)
errorcfpca_mu1 = inprod(res_cfpca$mu1-meanfit$pc_fit1,res_cfpca$mu1-meanfit$pc_fit1)
errorcfpca_mu2 = inprod(res_cfpca$mu2-meanfit$pc_fit2,res_cfpca$mu2-meanfit$pc_fit2)

errorcfpca_pc11= min(inprod(res_cfpca$pc11-pc1s$pc_fit1,res_cfpca$pc11-pc1s$pc_fit1),inprod(res_cfpca$pc11+pc1s$pc_fit1,res_cfpca$pc11+pc1s$pc_fit1))
errorcfpca_pc12= min(inprod(res_cfpca$pc12-pc1s$pc_fit2,res_cfpca$pc12-pc1s$pc_fit2),inprod(res_cfpca$pc12+pc1s$pc_fit2,res_cfpca$pc12+pc1s$pc_fit2))

errorcfpca_pc21 = min(inprod(res_cfpca$pc21-pc2s$pc_fit1,res_cfpca$pc21-pc2s$pc_fit1),inprod(res_cfpca$pc21+pc2s$pc_fit1,res_cfpca$pc21+pc2s$pc_fit1))
errorcfpca_pc22 = min(inprod(res_cfpca$pc22-pc2s$pc_fit2,res_cfpca$pc22-pc2s$pc_fit2),inprod(res_cfpca$pc22+pc2s$pc_fit2,res_cfpca$pc22+pc2s$pc_fit2))

o_sigma = res_cfpca$sigma
o_pik = res_cfpca$pik
owik = res_cfpca$wikmat

classest = ifelse(owik[,2]>0.5,1,2)

classtrue = clusters%>%do.call(c,.)

correctrate = max(length(which(classest==classtrue))/ssize,1-length(which(classest==classtrue))/ssize)


output = c(	
errorcfpca_mu1,errorcfpca_mu2,
errorcfpca_pc11,errorcfpca_pc12,
errorcfpca_pc21,errorcfpca_pc22,
o_sigma,o_pik,correctrate,mean_correctrate,pacecorrectrate
)


if (it==1){
outputmat = output
}else{
outputmat = rbind(outputmat,output)
}

}














