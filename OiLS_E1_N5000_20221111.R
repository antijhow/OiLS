set.seed(1000)
start = proc.time()
loop = 100
n = 5000
p = 1
q = 1

# Experiment I
bL0 = c(7,-7)
bR0 = c(-1,5)

a = 0.6 
b = 0.4
model.sel = array(0,loop)
sel.final = matrix(4,loop,2)
M1.est = M2.est = array(0,loop)
lambda= log(n)
mean.break2 = mean.break = var.break = array(0,loop)
for (i in 1:loop)
{
	x = sort(runif(n,0,1))
	x.L = x[x<=a]
	x.R = x[x>a]
	X.L = matrix(1, length(x.L), p+1)
	X.R = matrix(1, length(x.R), q+1)
	if (p>0)
	{
		for (j in 1:p) X.L[,j+1] = (x.L)^j
	}
	if (q>0)
	{
		for (j in 1:q) X.R[,j+1] = (x.R)^j
	}
	m.L = X.L%*%bL0
	m.R = X.R%*%bR0
	mu = c(m.L,m.R)
	
	error = rnorm(n,0,1)
   	error=  c(sqrt(2-1*x[x<=b]),sqrt(5-2*x[x>b]))*error
		
	y = mu + error
	aL.var = aL = x[0.3*n]
	aR.var = aR = x[0.7*n]
	temp = var.model(x,y,p,q,lambda,aL,aR,aL.var,aR.var)
	sel.final[i,] = temp$M.hat
	mean.break[i] = temp$a.est.all[5,1]
	mean.break2[i] = temp$a.est.all[9,1]
	var.break[i] = temp$a.est.all[9,2]
	print(i)
}

length(which(sel.final[,1]==2 & sel.final[,2]==2)) # depends on the underlying model
length(which(sel.final[,1]==2))



mean(mean.break2)
sd(mean.break2)

mean(mean.break)
sd(mean.break)

mean(var.break)
sd(var.break)	

proc.time()-start



par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(x, y, lty=1, col = "grey",type="l",main="N=5000",xlim=c(0,1),ylim=c(-2,12),xlab="",ylab="",cex.main = 3,cex.axis=1.5)              # Create first plot

par(new = TRUE)                             # Add new plot
plot(density(var.break),xlim=c(0,1), lty=2,lwd=4,axes=FALSE,xlab="",ylab="",col="red",main="",cex.main =3,cex.axis=1.5)
lines(density(mean.break2),lty=5,col="blue",lwd=4)
abline(v=a,lty=4,col="orange",lwd=4)
abline(v=b,lty=3,col="green",lwd=4)
axis(side = 4,  at = c(0,75,150,225,300),cex.axis=1.5)      # Add second axis
legend("topleft", c(expression(y),expression("dist"*(hat(a)["2,2"])),expression("dist"*(hat(alpha)[2])),expression(alpha^{(0)}),expression(a^{(0)})),lty=c(1,2,5,4,3),col=c("grey","red","blue","orange","green"),cex=1.1,lwd=4 )







