fit.y = function(x,y,p,q,bL,bR,a,Z=NULL,bZ=NULL) #Z is an N by qZ matrix
{
	n = length(y)
	X.L = matrix(1, n, p+1)
	X.R = matrix(1, n, q+1)
	if (p>0)
	{
		for (i in 1:p) X.L[,i+1] = x^i
	}
	if (q>0)
	{	
		for (i in 1:q) X.R[,i+1] = x^i
	}
	X.L[x>a,] = 0
	X.R[x<=a,] = 0
	H = cbind(X.L,X.R,Z)
	b.est = c(bL,bR,bZ)
	y.hat = H%*%b.est
	return(y.hat)
}

# to calculate mean-squared-error
mse.loss = function(x,y,p,q,bL,bR,a,Z=NULL,bZ=NULL) 
{
	n = length(y)
	y.hat = fit.y(x,y,p,q,bL,bR,a,Z,bZ)
	mse = sum((y-y.hat)^2)/n
	return(mse)
}

# to estimate the regression coefficient vectors
b.est = function(x,y,p,q,aL,aR,Z=NULL)
{
	solve.disturb = 0.00001
	n = length(x)
	X.L = matrix(1, n, p+1)
	X.R = matrix(1, n, q+1)
	if (p>0)
	{
		for (i in 1:p) X.L[,i+1] = x^i
	}
	if (q>0)
	{	
		for (i in 1:q) X.R[,i+1] = x^i
	}
	X.L[x>aL,] = 0
	X.R[x<=aR,] = 0
	if (is.null(Z)) 
	{
		Z.L=Z.R= NULL
	} else {
		Z.L = Z.R = Z
		Z.L[x>aL,]=0
		Z.R[x<=aR,]=0
	}
	H = cbind(X.L,X.R,Z.L+Z.R)
	temp = t(H)%*%(H) 
	diag(temp) = diag(temp) +solve.disturb
	b.hat = solve(temp)%*%t(H)%*%y	
	return (list(bL = b.hat[1:(p+1),1],bR=b.hat[(p+2):(p+q+2)],bZ=b.hat[-c(1:(p+q+2)),1]))
}

# to estimate change-point under M2 (discontinuous change) 
change.est = function(x,y,aL,aR,p,q,Z=NULL)
{
	n = length(x)
	temp.b = b.est(x,y,p,q,aL,aR,Z)
	bL = temp.b$bL
	bR = temp.b$bR
	if (is.null(Z)) bZ=NULL else bZ = temp.b$bZ
	test = 1
	previous.mse = Inf
	while(test<=3)
	{			
		#optim method
		change=optimize(mse.loss,c(min(x)-0.0001,max(x)+0.0001),x=x,y=y,p=p,q=q,bL=bL,bR=bR,Z=Z,bZ=bZ)
		temp.a = change$minimum
		temp.b = b.est(x,y,p,q,temp.a,temp.a,Z)
		bL = temp.b$bL
		bR = temp.b$bR
		if (is.null(Z)) bZ=NULL else bZ = temp.b$bZ
		temp.mse = mse.loss(x,y,p,q,bL,bR,temp.a,Z,bZ)
		if (previous.mse > temp.mse)
		{
			previous.mse = temp.mse
			a = temp.a
			test = 1
		} else test = test + 1			
	}
	final.b = b.est(x,y,p,q,a,a,Z)
	bL = final.b$bL
	bR = final.b$bR
	if (is.null(Z)) bZ=NULL else bZ = final.b$bZ
	return (list(a=a,bL=bL,bR=bR,bZ=bZ))
}

mse.only.a = function(x,y,p,q,a,Z=NULL)
{
	temp.b = b.est(x,y,p,q,a,a,Z)
	bL = temp.b$bL
	bR = temp.b$bR
	if (is.null(Z)) bZ=NULL else bZ = temp.b$bZ
	rL = rR =0
	for (j in 1:(p+1)) rL = rL + bL[j]*(a^{j-1})
	for (j in 1:(q+1)) rR = rR + bR[j]*(a^{j-1})
	bR[1] = bR[1] + (rL-rR)
	mse = mse.loss(x,y,p,q,bL,bR,a,Z,bZ)
	return (mse)
}
# to estimate the break-point under M1 (continuous break)
break.est = function(x,y,p,q,a.est,Z=NULL)
{
	temp.a = a.est
	a = temp.a
	n = length(x)
	test = 1
	previous.mse = mse.only.a(x,y,p,q,a,Z)
	temp.b = b.est(x,y,p,q,a,a,Z)
	bL = temp.b$bL
	bR = temp.b$bR
	if (is.null(Z)) bZ=NULL else bZ = temp.b$bZ
	roots = Re(polyroot(bL-bR))
	cand.roots = roots[roots >x[1] & roots < x[n]]
	r = length(cand.roots)
	if (r!=0) 
	{
		compare.mse = array(0,r)
		for (i in 1:r)
		{
			temp.a = cand.roots[i]
			temp.b = b.est(x,y,p,q,temp.a,temp.a,Z)
			bL = temp.b$bL
			bR = temp.b$bR
			if (is.null(Z)) bZ=NULL else bZ = temp.b$bZ
			rL = rR =0
			for (j in 1:(p+1)) rL = rL + bL[j]*(temp.a)^(j-1)
			for (j in 1:(q+1)) rR = rR + bR[j]*(temp.a)^(j-1)
			bR[1] = bR[1] + (rL-rR)
			compare.mse[i] = mse.loss(x,y,p,q,bL,bR,temp.a,Z,bZ)
		}
		temp.mse = min(compare.mse)
		if (previous.mse > temp.mse)
		{
			temp.a = cand.roots[which(compare.mse==temp.mse)]
			a = temp.a[1]
		}
	}
	temp.b = b.est(x,y,p,q,a,a,Z)
	bL = temp.b$bL
	bR = temp.b$bR
	if (is.null(Z)) bZ=NULL else bZ = temp.b$bZ
	rL = rR =0
	for (j in 1:(p+1)) rL = rL + bL[j]*(a^{j-1})
	for (j in 1:(q+1)) rR = rR + bR[j]*(a^{j-1})
	bR[1] = bR[1] + (rL-rR)
	return (list(a=a,bL=bL,bR=bR,bZ=bZ))							
}

mean.model = function(x,y,p,q,lambda,aL,aR,Z=NULL) 
{
	n = length(x)
	M2.est = change.est(x,y,aL,aR,p,q,Z)
	M1.est = break.est(x,y,p,q,M2.est$a,Z)
	ma = max(p,q)
	X = matrix(1, n, ma+1)
	if (ma>0)
	{
		for (i in 1:ma) X[,i+1]=x^i
	}
	solve.disturb = 0.00001
	H = cbind(X,Z)
	temp = t(H)%*%H
	diag(temp) = diag(temp)+ solve.disturb
	M0.est = solve(temp)%*%t(H)%*%y
	M0.y = (H%*%M0.est)[,1]
	M0.sse = sum((y-M0.y)^2)
	M1.mse = mse.loss(x,y,p,q,M1.est$bL,M1.est$bR,M1.est$a,Z,M1.est$bZ)
	M2.mse = mse.loss(x,y,p,q,M2.est$bL,M2.est$bR,M2.est$a,Z,M2.est$bZ)
	crit = array(0,3)
	crit[1] = M0.sse/M2.mse + (p+1)*lambda
	crit[2] = n*M1.mse/M2.mse + (p+q+2)*lambda 
	crit[3] = n + (p+q+3)*lambda
	select = which(crit == min(crit))[1]-1
	return ( list(M.hat = select, M0.b = M0.est,M0.mse = M0.sse/n, M1.a = M1.est$a, M1.bL = M1.est$bL, M1.bR = M1.est$bR,M1.bZ=M1.est$bZ, M1.mse=M1.mse, M2.a = M2.est$a, M2.bL = M2.est$bL, M2.bR = M2.est$bR,M2.bZ=M2.est$bZ, M2.mse = M2.mse, crits = crit) )
}

var.model = function(x,y,p,q,lambda,aL,aR,aL.var,aR.var,Z.mean=NULL,Z.var=NULL)
{
	n = length(x)
	M2.est = change.est(x,y,aL,aR,p,q,Z.mean)
	M1.est = break.est(x,y,p,q,M2.est$a,Z.mean)
	M1.mse = mse.loss(x,y,p,q,M1.est$bL,M1.est$bR,M1.est$a,Z.mean,M1.est$bZ)
	M2.mse = mse.loss(x,y,p,q,M2.est$bL,M2.est$bR,M2.est$a,Z.mean,M2.est$bZ)
	ma = max(p,q)
	X = matrix(1, n, ma+1)
	if (ma>0)
	{
		for (i in 1:ma) X[,i+1]=x^i
	}
	solve.disturb = 0.00001
	H = cbind(X,Z.mean)
	temp = t(H)%*%H 
	diag(temp) = diag(temp) +solve.disturb
	M0.est = solve(temp)%*%t(H)%*%y
	M0.y = (H%*%M0.est)[,1]
	M0.sse = sum((y-M0.y)^2)
	a.est.all = matrix(NA,9,2)
	p.var = p
	q.var = q
	lambda.var = lambda
	crits = matrix(4,3,3)
	M0.sq.res = (y - M0.y)^2
	temp.M0 = mean.model(x,M0.sq.res,p.var,q.var,lambda.var,aL.var,aR.var,Z.var)
	a.est.all[2,2] = temp.M0$M1.a
	a.est.all[3,2] = temp.M0$M2.a
	crits[1,] = temp.M0$crits + M0.sse/M2.mse + (p+1)*lambda
	M1.y = fit.y(x,y,p,q,M1.est$bL,M1.est$bR,M1.est$a,Z.mean,M1.est$bZ)
	M1.sq.res = (y - M1.y)^2
	temp.M1 = mean.model(x,M1.sq.res,p.var,q.var,lambda.var,aL.var,aR.var,Z.var)
	a.est.all[4:6,1] = M1.est$a
	a.est.all[5,2] = temp.M1$M1.a
	a.est.all[6,2] = temp.M1$M2.a
	crits[2,] = temp.M1$crits + n*M1.mse/M2.mse + (p+q+2)*lambda 
	M2.y = fit.y(x,y,p,q,M2.est$bL,M2.est$bR,M2.est$a,Z.mean,M2.est$bZ)
	M2.sq.res = (y - M2.y)^2
	temp.M2 = mean.model(x,M2.sq.res,p.var,q.var,lambda.var,aL.var,aR.var,Z.var)
	a.est.all[7:9,1] = M2.est$a
	a.est.all[8,2] = temp.M2$M1.a
	a.est.all[9,2] = temp.M2$M2.a
	crits[3,] = temp.M2$crits + n + (p+q+3)*lambda
	select = which(crits == min(crits),arr.ind=TRUE)-1
	return( list(M.hat = select, a.est.all = a.est.all) )	
}

forward.detection = function(x,y,p,cut.set,penalty=1,Z=NULL)
{
	M = length(cut.set)
	cut.temp = c(min(x)-0.0001,cut.set,max(x)+0.0001)
	cut.final = NULL
	total.mse = 0
	for (i in 1:(M+1))
	{
		index.set = which(x>cut.temp[i] & x<=cut.temp[i+1])	
		Z.temp = Z[index.set,]
		x.temp = x[index.set]
		y.temp = y[index.set]
		n = length(x.temp)
		temp.select = mean.model(x.temp,y.temp,p,p,lambda=penalty*log(n),aL=x.temp[0.3*n],aR=x.temp[0.7*n],Z=Z.temp) 
		M.hat = temp.select$M.hat[1]
		if (M.hat == 1)
		{
			cut.final = c(cut.final,temp.select$M1.a,cut.temp[i+1])
			total.mse = total.mse + n*temp.select$M1.mse
		} else if (M.hat == 2) {
			cut.final = c(cut.final,temp.select$M2.a,cut.temp[i+1])
			total.mse = total.mse + n*temp.select$M2.mse
		} else {
			cut.final = c(cut.final,cut.temp[i+1])
			total.mse = total.mse + n*temp.select$M0.mse
		}
	}
 	cut.final = cut.final[-length(cut.final)]
	return( list(cut.final = cut.final, total.mse = total.mse) )
}

backward.deletion = function(x,y,p,cut.set,penalty=1,Z=NULL)
{
	cond = 1
	cut.temp = cut.set
	while (cond != 0)
	{
		M = length(cut.temp)
		cut.append = c(min(x)-0.0001,cut.temp,max(x)+0.0001)
		total.mse = 0
		for (i in 2:(M+1))
		{
			index.set = which(x>cut.append[i-1] & x<=cut.append[i+1])	
			Z.temp = Z[index.set,]
			x.temp = x[index.set]
			y.temp = y[index.set]
			n = length(x.temp) 
			temp.select = mean.model(x.temp,y.temp,p,p,lambda=penalty*log(n),aL=x.temp[0.3*n],aR=x.temp[0.7*n],Z=Z.temp) 
			M.hat = temp.select$M.hat[1]
			if (M.hat == 1)
			{
				cut.append[i] = temp.select$M1.a
				total.mse = total.mse + n*temp.select$M1.mse
			} else if (M.hat == 2) {
				cut.append[i] = temp.select$M2.a 
				total.mse = total.mse + n*temp.select$M2.mse

			} else {
				cut.append[i] = cut.append[i-1]
				total.mse = total.mse + n*temp.select$M0.mse
			}
		}
		cut.append
		cut.temp = unique(cut.append)[-1]
		cut.temp = cut.temp[-length(cut.temp)]
		if (M == length(cut.temp) || length(cut.temp)==0) cond = 0
	}
	return ( list(cut.final = cut.temp, total.mse = total.mse) )
}

for.back = function(x,y,p,penalty=1,Z=NULL)
{
	n = length(x)
	cut.seed = change.est(x,y,aL=x[0.3*n],aR=x[0.7*n],p=p,q=p,Z=Z)$a
	count = 1
	cut.set = forward.detection(x,y,p,cut.seed,penalty,Z)
	cut.min.mse = cut.set$total.mse
	cut.min.cut = cut.set$cut.final
	while (count <=3)
	{
		cut.temp = backward.deletion(x,y,p,cut.set$cut.final,penalty,Z)
		if (length(cut.temp$cut.final)==0) cut.temp$cut.final = runif(1,min(x),max(x))
		cut.set = forward.detection(x,y,p,cut.temp$cut.final,penalty,Z)
		if (cut.set$total.mse < cut.min.mse) 
		{
			count = 1
			cut.min.cut = cut.set$cut.final
			cut.min.mse = cut.set$total.mse
		} else count = count +1
	}
	cut.temp = backward.deletion(x,y,p,cut.set$cut.final,penalty,Z)
	cut.min.cut = cut.temp$cut.final
	return (list(cut.final = cut.min.cut))
}

multi.model.detecion = function(x,y,p,cut.set,penalty=1,Z=NULL)
{
	M = length(cut.set)
	cut.append = c(min(x)-0.0001,cut.set,max(x)+0.0001)
	model.final = array(0,M)
	for (i in 2:(M+1))
	{
		index.set = which(x>cut.append[i-1] & x<=cut.append[i+1])	
		Z.temp = Z[index.set,]
		x.temp = x[index.set]
		y.temp = y[index.set]
		n = length(x.temp) 
		temp.select = mean.model(x.temp,y.temp,p,p,lambda=penalty*log(n),aL=x.temp[0.3*n],aR=x.temp[0.7*n],Z=Z.temp) 
		model.final[i-1] = temp.select$M.hat
	}
	return (list(model.final = model.final))
}


