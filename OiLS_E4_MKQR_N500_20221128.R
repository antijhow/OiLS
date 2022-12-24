library(MultiKink)
start = proc.time()
loop = 100
n = 500
p = 1
qz=1
cut.set = c(-1,2)
M = length(cut.set)
coef.matrix = matrix(0,p+1,M+1)
coef.matrix[1,] = c(1,-2,-10)
coef.matrix[2,] = c(1,-2,2)
bZ0 = c(1)
jump = NULL
for (i in 1:M)
{
	temp = coef.matrix[1,i]+coef.matrix[2,i]*cut.set[i]-(coef.matrix[1,i+1]+coef.matrix[2,i+1]*cut.set[i])
	jump = c(jump,temp)
}
jump
#coef.matrix[1,1] = coef.matrix[1,1]-jump[1]
#coef.matrix[1,M+1]= coef.matrix[1,M+1]+jump[M]
#coef.matrix
est.cut.m = est.cut.5 = NULL
model.cut.5 = NULL
est.cut= NULL
count.m = count = 0
seed=c(1:13,15:21,28,31,33:38,41:43,45:48,50:53,55:57,59:61, 63:71,73:77,79:83,85:102,104:116,118,120:1000)
for (i in 1:loop)
{
	set.seed(seed[i])
	#x = sort(rnorm(n,0,1))
	x = sort(runif(n,-5,5))
	mu = NULL
	cut.temp = c(min(x)-0.0001,cut.set,max(x)+0.0001)
	for ( k in 1:(M+1) )
	{
		index.set = which(x>cut.temp[k] & x<=cut.temp[k+1])
		x.temp = x[index.set] 
		X.temp = matrix(1, length(x.temp), p+1)
		if (p>0)
		{
			for (j in 1:p) X.temp[,j+1] = (x.temp)^j
		}
		mu.seg = X.temp%*%coef.matrix[,k]
		mu = c(mu,mu.seg)
	}
	error = rnorm(n,0,1)
	Z = matrix(rnorm(n*qz,1,1),n,qz)
	y = mu + Z%*%bZ0+ error
	fit = mkqr.bea(y,x,Z,0.5,Cn=log(n))
	if (fit$n.psi == 2)
	{
		count.m = count.m+1
		est.cut.m = rbind(est.cut.m,fit$psi.est)
	}

	#final.det = for.back(x,y,p,penalty=2,Z=Z)
	#cut.number = length(final.det$cut.final)
	#if (cut.number == 2) 
	#{
	#	count = count + 1
	#	est.cut.5 = rbind(est.cut.5,final.det$cut.final)
	#} 
	#est.cut = c(est.cut,final.det$cut.final)
	print(i)
}

count.m/loop
#count/loop
proc.time()-start



apply(est.cut.m,2,mean)
apply(est.cut.m,2,sd)

#apply(est.cut.5,2,mean)
#apply(est.cut.5,2,sd)


