

findmean = function(observed,timepoints,threshold=1e-3,minit=3,gamma=0){


yall = as.matrix(do.call(c,observed))



xalpha = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati

})%>%do.call(rbind,.)

xall=xalpha
mu1all = rep(0,length(yall))
mu2all = rep(0,length(yall))
muall = cbind(mu1all,mu2all)

groupindi  = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])	
	as.matrix(rep(subj_index,length(timei)))

})%>%do.call(rbind,.)


sigma0 = sqrt( mean(lm.fit(xall,yall)$residual^2))
coef0=lm.fit(xall,yall)$coef

pik = c(0.5,0.5)
sigmak= c(sigma0,sigma0*2)
betak=list(as.matrix(coef0),as.matrix(coef0))

numofgroup = max(groupindi)

for(iter in 1:100){

##E step
wik = list()
for(i in 1:numofgroup){
tempwi=list()
for(k in 1:2){

ytemp = yall[which(groupindi==i),]
xtemp  = xall[which(groupindi==i),] %*% betak[[k]] + muall[which(groupindi==i),k]

output=Reduce('*',lapply(1:length(ytemp),function(ind){
yy=mpfr(ytemp[ind],500)
xx=xtemp[ind]
(dnorm(yy,xx,sigmak[k]))
}))

tempwi[[k]] =  output* pik[k]
}
tempwi2 = tempwi
for(k in 1:2){
tempwi2[[k]]=tempwi[[k]]/Reduce("+",tempwi)
}
wik[[i]]=sapply(tempwi2,as.numeric)
}


#M step
##Update pik

wikmat = Reduce(rbind,wik)
wikmat= cbind(seq(1,length(wik)),wikmat)
colnames(wikmat)=c('index',paste0("w",seq(1,2)))
colnames(groupindi)='index'

nitbl = as.data.frame(table(groupindi))
colnames(nitbl)=c('index','ni')
ni = merge(groupindi,nitbl,by='index')[,2]
weightvec = merge(groupindi,wikmat,by='index')

wiinfo = merge(nitbl,wikmat,by='index')%>%arrange(index)

for(k in 1:2){
pik[k] = sum(wiinfo[,k+2] /wiinfo$ni)/sum(1/wiinfo$ni)
}


for(k in 1:2){

model = lm.wfit(xall,yall-muall[,k],weightvec[,k+1]/ni)

betak[[k]] = as.matrix(unlist(model$coef),ncol=1)

sigmak[k] = sqrt(sum(weightvec[,k+1]/ni*model$residual^2)/sum(weightvec[,k+1]/ni))

W = diag(weightvec[,k+1]/ni)
Y=yall-muall[,k]
#solve(t(xall) %*% W %*% xall) %*% t(xall) %*% W %*% Y
#residual = yall-muall[,k] - xall %*% betak[[k]]
#sqrt(sum(weightvec[,k+1]/ni*residual^2)/sum(weightvec[,k+1]/ni))

}



	if(iter==1){
	betabefore = Reduce('rbind',betak)
	}else{
	thresh =  mean(abs(Reduce('rbind',betak)-betabefore))
	print(thresh)
	betabefore = Reduce('rbind',betak)
	if (thresh<threshold){
	break
	}
	}

	
}#for iter in 1:100 end

	
pc_fit1 = fd(betak[[1]], spline_basis)
pc_fit2 = fd(betak[[2]], spline_basis)

#betak = coef(pc_fit)%>%as.numeric

return(list(betak=betak,pik=pik,sigmak=sigmak,wikmat=wikmat,pc_fit1 = pc_fit1,pc_fit2 = pc_fit2))

}


#sigmainit = meanfit$sigmak
#beta1 = betainit
#beta2 = betainit
first_FPC  = function(beta1,beta2,observed,observedmu1,observedmu2,timepoints,sigmainit,threshold=1e-3,maxit=50,gamma=0){


pc_fit1 = fd(beta1, spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
beta1 = coef(pc_fit1)%>%as.numeric

pc_fit2 = fd(beta2, spline_basis)
pc_fit2 = (1/sqrt(inprod(pc_fit2, pc_fit2))*pc_fit2)
beta2 = coef(pc_fit2)%>%as.numeric

yall = as.matrix(do.call(c,observed))
mu1all = as.matrix(do.call(c,observedmu1))
mu2all = as.matrix(do.call(c,observedmu2))
muall = cbind(mu1all,mu2all)


groupindi  = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])	
	as.matrix(rep(subj_index,length(timei)))

})%>%do.call(rbind,.)

sigmak=sigmainit
pik = c(0.5,0.5)
betak=list(as.matrix(beta1),as.matrix(beta2))
numofgroup = max(groupindi)

#initialize
alpha_fit1 = function(subj_index){
	(timei =timepoints[[subj_index]])	
	xmati = eval.fd(timei, pc_fit1)	
	model = lm.fit(x=as.matrix(xmati),y=as.matrix(observed[[subj_index]]-observedmu1[[subj_index]]))						
	model$coefficient%>%as.numeric			
	#lm(observed[[subj_index]]~0+xmati)%>%coef%>%as.numeric
}
#initialize
sfit1= sapply(1:length(observed),alpha_fit1)
#initialize
alpha_fit2 = function(subj_index){
	(timei =timepoints[[subj_index]])	
	xmati = eval.fd(timei, pc_fit2)	
	model = lm.fit(x=as.matrix(xmati),y=as.matrix(observed[[subj_index]]-observedmu2[[subj_index]]))						
	model$coefficient%>%as.numeric			
	#lm(observed[[subj_index]]~0+xmati)%>%coef%>%as.numeric
}
#initialize
sfit2= sapply(1:length(observed),alpha_fit2)

#initialize
xalpha1 = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati*sfit1[subj_index]

})%>%do.call(rbind,.)
#initialize
xalpha2 = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati*sfit2[subj_index]

})%>%do.call(rbind,.)

xall = list(xalpha1,xalpha2)


for(iter in 1:maxit){
	
	
##E step
wik = list()
for(i in 1:numofgroup){
tempwi=list()
for(k in 1:2){

ytemp = yall[which(groupindi==i),]
xtemp  = xall[[k]][which(groupindi==i),] %*% betak[[k]] + muall[which(groupindi==i),k]

output=Reduce('*',lapply(1:length(ytemp),function(ind){
yy=mpfr(ytemp[ind],120)
xx=xtemp[ind]
(dnorm(yy,xx,sigmak[k]))
}))

tempwi[[k]] =  output* pik[k]
}
tempwi2 = tempwi
for(k in 1:2){
tempwi2[[k]]=tempwi[[k]]/Reduce("+",tempwi)
}
wik[[i]]=sapply(tempwi2,as.numeric)
}



wikmat = Reduce(rbind,wik)
wikmat= cbind(seq(1,length(wik)),wikmat)
colnames(wikmat)=c('index',paste0("w",seq(1,2)))
colnames(groupindi)='index'

nitbl = as.data.frame(table(groupindi))
colnames(nitbl)=c('index','ni')
ni = merge(groupindi,nitbl,by='index')[,2]
weightvec = merge(groupindi,wikmat,by='index')

wiinfo = merge(nitbl,wikmat,by='index')%>%arrange(index)

for(k in 1:2){
pik[k] = sum(wiinfo[,k+2] /wiinfo$ni)/sum(1/wiinfo$ni)
}


	
for (subiter in 1:10){
	
alpha_fit1 = function(subj_index){
	(timei =timepoints[[subj_index]])	
	xmati = eval.fd(timei, pc_fit1)	
	model = lm.fit(x=as.matrix(xmati),y=as.matrix(observed[[subj_index]]-observedmu1[[subj_index]]))						
	model$coefficient%>%as.numeric			
	#lm(observed[[subj_index]]~0+xmati)%>%coef%>%as.numeric
}

sfit1= sapply(1:length(observed),alpha_fit1)

alpha_fit2 = function(subj_index){
	(timei =timepoints[[subj_index]])	
	xmati = eval.fd(timei, pc_fit2)	
	model = lm.fit(x=as.matrix(xmati),y=as.matrix(observed[[subj_index]]-observedmu2[[subj_index]]))						
	model$coefficient%>%as.numeric			
	#lm(observed[[subj_index]]~0+xmati)%>%coef%>%as.numeric
}

sfit2= sapply(1:length(observed),alpha_fit2)


yem1 = do.call(c,observed)
yem2 = do.call(c,observed)


xalpha1 = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati*sfit1[subj_index]

})%>%do.call(rbind,.)

xalpha2 = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati*sfit2[subj_index]

})%>%do.call(rbind,.)


xall = list(xalpha1,xalpha2)

for(k in 1:2){

model = lm.wfit(xall[[k]],yall-muall[,k],weightvec[,k+1]/ni)
betak[[k]] = as.matrix(unlist(model$coef),ncol=1)
sigmak[k] = sqrt(sum(weightvec[,k+1]/ni*model$residual^2)/sum(weightvec[,k+1]/ni))

#W = diag(weightvec[,k+1]/ni)
#Y=yall-muall[,k]
#solve(t(xall) %*% W %*% xall) %*% t(xall) %*% W %*% Y
#residual = yall-muall[,k] - xall %*% betak[[k]]
#sqrt(sum(weightvec[,k+1]/ni*residual^2)/sum(weightvec[,k+1]/ni))
}


pc_fit1 = fd(betak[[1]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
betak[[1]] = coef(pc_fit1)%>%as.numeric

pc_fit2 = fd(betak[[2]], spline_basis)
pc_fit2 = (1/sqrt(inprod(pc_fit2, pc_fit2))*pc_fit2)
betak[[2]] = coef(pc_fit2)%>%as.numeric

}#for subiter end

xall = list(xalpha1,xalpha2)


	if(iter==1){
	betabefore = Reduce('rbind',betak)
	}else{
	thresh =  mean(abs(Reduce('rbind',betak)-betabefore))
	print(thresh)
	betabefore = Reduce('rbind',betak)
	if (thresh<threshold){
	break
	}
	}
	
}# for iter end




print("Done!")

return(list(betak=betak,pc_fit1 = pc_fit1,pc_fit2 = pc_fit2,sfit1 = sfit1,sfit2 = sfit2,sigmak=sigmak,pik=pik,wikmat=wikmat))
}


#beta1=coef_mat0[,3]
#beta2=coef_mat0[,3]
#betalist=previous_beta
third_FPC_conditional  = function(beta1,beta2, pc_index, observed,observedmu1,observedmu2,timepoints,sigmainit, betalist, threshold=1e-4,maxit=50,gamma=0){

observedx=observed
observedxmu1=observedmu1
observedxmu2=observedmu2
observedx[which(sapply(observed,length)<=pc_index)]=NULL
observedxmu1[which(sapply(observed,length)<=pc_index)]=NULL
observedxmu2[which(sapply(observed,length)<=pc_index)]=NULL
timepointsx = timepoints
timepointsx[which(sapply(observed,length)<=pc_index)]=NULL


pc_fit1 = fd(beta1, spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
beta1 = coef(pc_fit1)%>%as.numeric

pc_fit2 = fd(beta2, spline_basis)
pc_fit2 = (1/sqrt(inprod(pc_fit2, pc_fit2))*pc_fit2)
beta2 = coef(pc_fit2)%>%as.numeric


mu1all = as.matrix(do.call(c,observedxmu1))
mu2all = as.matrix(do.call(c,observedxmu2))
muall = cbind(mu1all,mu2all)


groupindi  = lapply(1:length(observedx),function(subj_index){

	(timei =timepointsx[[subj_index]])	
	as.matrix(rep(subj_index,length(timei)))

})%>%do.call(rbind,.)

sigmak=sigmainit
pik = c(0.5,0.5)
betak=list(as.matrix(beta1),as.matrix(beta2))
numofgroup = max(groupindi)


pc_fits_previous1 = lapply(1:length(betalist), function(x){
pc_fit1 = fd(betalist[[x]][[1]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
pc_fit1
})

pc_fits_previous2 = lapply(1:length(betalist), function(x){
pc_fit1 = fd(betalist[[x]][[2]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
pc_fit1
})


alpha_fit1 = function(subj_index){
	(timei =timepointsx[[subj_index]])
	
	xmat_previous = lapply(pc_fits_previous1, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit1)
	
	model = lm.fit(x=cbind(xmat_previous,as.matrix(xmati,ncol=1)),y=as.matrix(observedx[[subj_index]]-observedxmu1[[subj_index]]))						
	model$coefficient%>%as.numeric				
	
}

sfit1= lapply(1:length(observedx),alpha_fit1)%>%do.call(rbind,.)

alpha_fit2 = function(subj_index){
	(timei =timepointsx[[subj_index]])
	
	xmat_previous = lapply(pc_fits_previous2, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit2)
	
	model = lm.fit(x=cbind(xmat_previous,as.matrix(xmati,ncol=1)),y=as.matrix(observedx[[subj_index]]-observedxmu2[[subj_index]]))						
	model$coefficient%>%as.numeric				
	
}

sfit2= lapply(1:length(observedx),alpha_fit2)%>%do.call(rbind,.)



xalpha1 = lapply(1:length(observedx),function(subj_index){

	(timei =timepointsx[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit1[subj_index,ncol(sfit1)])#/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)

xalpha2 = lapply(1:length(observedx),function(subj_index){

	(timei =timepointsx[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit2[subj_index,ncol(sfit2)])#/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)


yem1 = lapply(1:length(observedx), function(i){
(timei =timepointsx[[i]])
yfits_previous = lapply(1:length(pc_fits_previous1), function(x){
	sfit1[i,x]*eval.fd(timei,pc_fits_previous1[[x]])%>%as.numeric
	})%>%do.call(rbind,.)%>%colSums

(observedx[[i]] - yfits_previous)
})%>%do.call(c,.)

yem2 = lapply(1:length(observedx), function(i){
(timei =timepointsx[[i]])
yfits_previous = lapply(1:length(pc_fits_previous2), function(x){
	sfit2[i,x]*eval.fd(timei,pc_fits_previous2[[x]])%>%as.numeric
	})%>%do.call(rbind,.)%>%colSums

(observedx[[i]] - yfits_previous)
})%>%do.call(c,.)

xall = list(xalpha1,xalpha2)

yall=list(as.matrix(yem1),as.matrix(yem2))


for(iter in 1:maxit){
	
	
##E step
wik = list()
for(i in 1:numofgroup){
tempwi=list()
for(k in 1:2){

ytemp = yall[[k]][which(groupindi==i),]
xtemp  = xall[[k]][which(groupindi==i),] %*% betak[[k]] + muall[which(groupindi==i),k]

output=Reduce('*',lapply(1:length(ytemp),function(ind){
yy=mpfr(ytemp[ind],120)
xx=xtemp[ind]
(dnorm(yy,xx,sigmak[k]))
}))

tempwi[[k]] =  output* pik[k]
}
tempwi2 = tempwi
for(k in 1:2){
tempwi2[[k]]=tempwi[[k]]/Reduce("+",tempwi)
}
wik[[i]]=sapply(tempwi2,as.numeric)
}



wikmat = Reduce(rbind,wik)
wikmat= cbind(seq(1,length(wik)),wikmat)
colnames(wikmat)=c('index',paste0("w",seq(1,2)))
colnames(groupindi)='index'

nitbl = as.data.frame(table(groupindi))
colnames(nitbl)=c('index','ni')
ni = merge(groupindi,nitbl,by='index')[,2]
weightvec = merge(groupindi,wikmat,by='index')

wiinfo = merge(nitbl,wikmat,by='index')%>%arrange(index)

for(k in 1:2){
pik[k] = sum(wiinfo[,k+2] /wiinfo$ni)/sum(1/wiinfo$ni)
}


	
for (subiter in 1:10){
	
	alpha_fit1 = function(subj_index){
	(timei =timepointsx[[subj_index]])
	
	xmat_previous = lapply(pc_fits_previous1, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit1)
	
	model = lm.fit(x=cbind(xmat_previous,as.matrix(xmati,ncol=1)),y=as.matrix(observedx[[subj_index]]-observedxmu1[[subj_index]]))						
	model$coefficient%>%as.numeric				
	
	}

sfit1= lapply(1:length(observedx),alpha_fit1)%>%do.call(rbind,.)

alpha_fit2 = function(subj_index){
	(timei =timepointsx[[subj_index]])
	
	xmat_previous = lapply(pc_fits_previous2, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit2)
	
	model = lm.fit(x=cbind(xmat_previous,as.matrix(xmati,ncol=1)),y=as.matrix(observedx[[subj_index]]-observedxmu2[[subj_index]]))						
	model$coefficient%>%as.numeric				
	
}

sfit2= lapply(1:length(observedx),alpha_fit2)%>%do.call(rbind,.)



xalpha1 = lapply(1:length(observedx),function(subj_index){

	(timei =timepointsx[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit1[subj_index,ncol(sfit1)])#/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)

xalpha2 = lapply(1:length(observedx),function(subj_index){

	(timei =timepointsx[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit2[subj_index,ncol(sfit2)])#/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)


yem1 = lapply(1:length(observedx), function(i){
(timei =timepointsx[[i]])
yfits_previous = lapply(1:length(pc_fits_previous1), function(x){
	sfit1[i,x]*eval.fd(timei,pc_fits_previous1[[x]])%>%as.numeric
	})%>%do.call(rbind,.)%>%colSums

(observedx[[i]] - yfits_previous)
})%>%do.call(c,.)

yem2 = lapply(1:length(observedx), function(i){
(timei =timepointsx[[i]])
yfits_previous = lapply(1:length(pc_fits_previous2), function(x){
	sfit2[i,x]*eval.fd(timei,pc_fits_previous2[[x]])%>%as.numeric
	})%>%do.call(rbind,.)%>%colSums

(observedx[[i]] - yfits_previous)
})%>%do.call(c,.)

xall = list(xalpha1,xalpha2)
yall=list(as.matrix(yem1),as.matrix(yem2))

for(k in 1:2){

#model = lm.wfit(xall[[k]],yall[[k]]-muall[,k],weightvec[,k+1]/ni)
#betak[[k]] = as.matrix(unlist(model$coef),ncol=1)
#sigmak[k] = sqrt(sum(weightvec[,k+1]/ni*model$residual^2)/sum(weightvec[,k+1]/ni))

W = diag(weightvec[,k+1]/ni)
Y=yall[[k]]-muall[,k]

qmat = 2*t(xall[[k]]) %*% W %*% xall[[k]]  #inverse 
pmat =as.numeric(2*( t(xall[[k]]) %*% W %*% Y * (-1)))

betamat = lapply(1:length(betalist), function(x){
betalist[[x]][[k]]
})%>%do.call(rbind,.)

E = inprod(spline_basis, spline_basis,0,0)
cmat = betamat%*%E
betatemp  = lsei::qp(qmat, pmat, c=cmat, d=rep(0,length(betalist)))

#betatemp= lsei::qp(qmat, pmat, e=cmat, f=c(rep(-1e-3,length(betalist)),rep(-1e-3,length(betalist))))

temppc_fit = fd(betatemp, spline_basis)
temppc_fit = (1/sqrt(inprod(temppc_fit, temppc_fit))*temppc_fit)
betak[[k]]   = coef(temppc_fit)%>%as.numeric

residual = yall[[k]]-muall[,k] - xall[[k]] %*% betak[[k]]

sigmak[k]=sqrt(sum(weightvec[,k+1]/ni*residual^2)/sum(weightvec[,k+1]/ni))
}


pc_fit1 = fd(betak[[1]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
betak[[1]] = coef(pc_fit1)%>%as.numeric

pc_fit2 = fd(betak[[2]], spline_basis)
pc_fit2 = (1/sqrt(inprod(pc_fit2, pc_fit2))*pc_fit2)
betak[[2]] = coef(pc_fit2)%>%as.numeric

}#for subiter end

xall = list(xalpha1,xalpha2)
yall=list(as.matrix(yem1),as.matrix(yem2))

	if(iter==1){
	betabefore = Reduce('rbind',betak)
	}else{
	thresh =  mean(abs(Reduce('rbind',betak)-betabefore))
	print(thresh)
	betabefore = Reduce('rbind',betak)
	if (thresh<threshold){
	break
	}
	}
	
}# for iter end




print("Done!")

return(list(betak=betak,pc_fit1 = pc_fit1,pc_fit2 = pc_fit2,sfit1 = sfit1,sfit2 = sfit2,sigmak=sigmak,pik=pik,wikmat=wikmat))
}


third_FPC_conditional2_notused  = function(beta3, observed, timepoints, betalist =previous_beta , threshold=1e-4,maxit=50,gamma=0,range_score=1e5){

thresh = 1
it = 1
R = inprod(spline_basis, spline_basis,2,2)
E = inprod(spline_basis, spline_basis,0,0)

pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

pc_fits_previous = lapply(1:length(betalist), function(x){
pc_fit1 = fd(betalist[[x]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
pc_fit1
})
value = -1e14

observed2 =observed
# observed2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints2 = timepoints
# timepoints2[which(sapply(observed,length)<=pc_index)]=NULL

beta3_before = rep(0, length(betalist[[1]]))
while (thresh >  threshold|it<minit){
	
	beta3_before = beta3
	value_before = value

alpha_fit = function(subj_index){
	(timei =timepoints2[[subj_index]])
	xmati = eval.fd(timei, pc_fit)
	sf = lm(observed2[[subj_index]]~xmati+0)%>%coef%>%as.numeric
	if(sf>range_score) sf =  range_score
	if(sf < -range_score) sf= - range_score
	return(sf)
}

sfit= lapply(1:length(observed2),alpha_fit)%>%do.call(rbind,.)

residual_fit = function(subj_index){
	(timei =timepoints2[[subj_index]])
	xmati = eval.fd(timei, pc_fit)
	lm(observed2[[subj_index]]~xmati+0)%>%residuals

}


rfit0= lapply(1:length(observed2),residual_fit)
rfit  = lapply(rfit0, function(x) {x^2%>%mean})%>%do.call(c,.)

value = as.numeric(mean(rfit)+gamma*inprod(pc_fit,pc_fit,2,2))
# if(abs((value_before - value)/value_before) < threshold) break
N = sapply(observed2,length)%>%sum

yem = lapply(1:length(observed2), function(i){
(timei =timepoints2[[i]])
(observed2[[i]])/sqrt(length(timei))/sqrt(N)
})%>%do.call(c,.)


xalpha = lapply(1:length(observed2),function(subj_index){

	(timei =timepoints2[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit[subj_index,ncol(sfit)])/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)


A = xalpha
qmat = 2*(t(A)%*%A+gamma*R)
pmat = as.numeric(-2*t(yem)%*%A)
betamat=  do.call(rbind,betalist)
cmat = rbind(betamat%*%E,-betamat%*%E)
# beta3  = lsei::qp(qmat, pmat, c=cmat, d=rep(0,length(betalist)))
beta3  = lsei::qp(qmat, pmat, e=cmat, f=c(rep(-1e-3,length(betalist)),rep(-1e-3,length(betalist))))
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
# plot(pc_fit)

if(inprod(pc_fit)<0) {
	pc_fit = - pc_fit
	
}
beta3 = coef(pc_fit)%>%as.numeric
thresh =  max(abs(beta3_before-beta3))

# loss_value= mean((yem - A%*%beta3)^2)+as.numeric(gamma*t(beta3)%*%R%*%beta3)
#print(value)
it = it+1
if(it%%10==0) {print(it);print(thresh)}

if(it>maxit){
break
}

}
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric
print("Done!")
cat(sprintf('The threshold is %s . \n',thresh))
return(list(beta=beta3,pc_fit = pc_fit, sfit = sfit, previous_beta =betalist, residual= rfit0, thresh = thresh,it=it,value=value,gamma=gamma))
}


pred_SOAP_notused = function(betalist,ylist, tlist,spline_basis,nminus=1){

cfits = lapply(1:length(ylist), function(x){
	print(x)
	timei = tlist[[x]]
	xmat = lapply(1:length(betalist), function(i){
		pc_fit  = fd(betalist[[i]], spline_basis)
		eval.fd(timei, pc_fit)%>%as.numeric
	})%>%do.call(cbind,.)
	index_pc = min(length(ylist[[x]])-nminus, length(betalist))
	if(index_pc<=1) {
	    index_pc=1
		cfit = rep(0,length(betalist))
		# cfit[which.min(abs(ylist[[x]]/xmat))]=(ylist[[x]]/xmat)[which.min(abs(ylist[[x]]/xmat))]

		cfit_temp  = sapply(1:length(betalist),function(j) {
			lm(ylist[[x]]~0+xmat[,j])%>%coef%>%as.numeric
		})
		cfit[which.min(abs(cfit_temp))] = cfit_temp[which.min(abs(cfit_temp))]
		
	} else {
			
			cfit = lm(ylist[[x]]~0+xmat[,1:index_pc])%>%coef%>%as.numeric
	}
	if(length(cfit)<length(betalist)) cfit = c(cfit, rep(0,length(betalist)-length(cfit)))
	cfit
})%>%do.call(rbind,.)

yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	rowSums(mapply("*",cfit, betalist,SIMPLIFY=TRUE))
}
)%>%do.call(cbind,.)

yfitsfd = fd(yfits, spline_basis)


residuals = sapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	resid = ylist[[x]] - as.numeric(eval.fd(timei,yfitsfd[x]))
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
	
})
sigmahat = mean(residuals[!is.na(residuals)])
list(sigmahat = sigmahat, yfd_fit=  yfitsfd,sfit = cfits)

}