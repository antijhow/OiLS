start = proc.time()
loop = 1000
n = 500
p = 1
qz=3
cut.set = c(-3,0,3,6)
M = length(cut.set)
coef.matrix = matrix(0,p+1,M+1)
coef.matrix[1,] = c(1,-8,-8,4,-16)
coef.matrix[2,] = c(1,-2,2,-2,2)

bZ0 = c(1,-2,3)

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
est.cut.5 = NULL
model.cut.5 = NULL
est.cut= NULL
count = 0
for (i in 1:loop)
{
	set.seed(2000+i)
	x = sort(runif(n,-5,8))
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
	error = (1+0.2*x)*rnorm(n,0,1)
	#error =  (1+0.2*x)*rt(n,3)
	
	Z = matrix(rnorm(n*qz,0,1),n,qz)

	y = mu + Z%*%bZ0+ error
	
	final.det = for.back(x,y,p,penalty=2,Z=Z)
	cut.number = length(final.det$cut.final)
	if (cut.number == 4) 
	{
		count = count + 1
		est.cut.5 = rbind(est.cut.5,final.det$cut.final)
		model.cut = multi.model.detecion(x,y,p,final.det$cut.final,penalty=2,Z)
		model.cut.5 = rbind(model.cut.5,model.cut$model.final)
	} else {m = i}
	est.cut = c(est.cut,final.det$cut.final)
	print(i)
}

count/loop
proc.time()-start




plot(x, y, lty="solid", col = "grey",type="l",main="N=500",cex.main=3,cex.axis=1.5,ylab="",xlab="")              # Create first plot
abline(v=est.cut,col="blue",lty="dashed")
abline(v=cut.set,col="red",lty="dotted",lwd=4)
legend("bottomright",c(expression(y),expression("dist"*(hat(alpha)[m])),expression(alpha[m]^{(0)})),lty=c("solid","dashed","dotted"),col=c("grey","blue","red"),cex=1.1 ,lwd=4,bg="white")


apply(est.cut.5,2,mean)
apply(est.cut.5,2,sd)

length(which(model.cut.5[,1]==1))/loop
length(which(model.cut.5[,2]==1))/loop
length(which(model.cut.5[,3]==1))/loop
length(which(model.cut.5[,4]==2))/loop

