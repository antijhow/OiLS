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
	solve.disturb = 0.0000001
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
	y.hat = as.vector(H %*% b.hat)
	mse = mean((y - y.hat)^2)
	return(list(bL = b.hat[1:(p+1),1],bR = b.hat[(p+2):(p+q+2),1],bZ = b.hat[-c(1:(p+q+2)),1],mse = mse))
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

# Fast MSE for a given threshold and fixed regression coefficients
mse.loss.fixed = function(x, y, p, q, bL, bR, a, Z=NULL, bZ=NULL)
{
  n = length(y)
  
  # Evaluate the left and right polynomial functions once
  y.L = rep(bL[1], n)
  y.R = rep(bR[1], n)
  
  if (p > 0) {
    for (j in 1:p)
      y.L = y.L + bL[j + 1] * x^j
  }
  
  if (q > 0) {
    for (j in 1:q)
      y.R = y.R + bR[j + 1] * x^j
  }
  
  # Select the appropriate regime
  y.hat = ifelse(x <= a, y.L, y.R)
  
  # Add the common Z component
  if (!is.null(Z))
    y.hat = y.hat + as.vector(Z %*% bZ)
  
  mse = mean((y - y.hat)^2)
  
  return(mse)
}

# to estimate change-point under M2 (discontinuous change) 
# OiLS with grid search:
# coefficients are fixed during each threshold-update step
change.est.grid = function(x, y, aL, aR, p, q, Z=NULL)
{
  n = length(x)
  temp.b = b.est(x, y, p, q, aL, aR, Z)
  bL = temp.b$bL
  bR = temp.b$bR
  if (is.null(Z))
    bZ = NULL
  else
    bZ = temp.b$bZ
  test = 1
  previous.mse = Inf
  ord = order(x)
  x.sort = x[ord]
  y.sort = y[ord]
  if (is.null(Z))
    Z.sort = NULL
  else
    Z.sort = Z[ord, , drop=FALSE]
  grid = sort(unique(x))
  while (test <= 3)
  {
    # ========================================================
    # Threshold update:
    # regression coefficients are fixed
    # ========================================================
    
    # Left fitted values
    fit.L = rep(bL[1], n)
    if (p > 0)
    {
      for (j in 1:p) fit.L = fit.L + bL[j+1] * x.sort^j
    }
    # Right fitted values
    fit.R = rep(bR[1], n)
    if (q > 0)
    {
      for (j in 1:q) fit.R = fit.R + bR[j+1] * x.sort^j
    }
    
    # Common Z component
    if (!is.null(Z.sort))
    {
      fit.Z = as.vector(Z.sort %*% bZ)
      fit.L = fit.L + fit.Z
      fit.R = fit.R + fit.Z
    }
    
    # Squared residuals under left/right assignments
    loss.L = (y.sort - fit.L)^2
    loss.R = (y.sort - fit.R)^2
    
    # --------------------------------------------------------
    # Cumulative RSS
    # --------------------------------------------------------
    
    cum.L = cumsum(loss.L)
    
    # Total right loss minus observations already assigned left
    cum.R = rev(cumsum(rev(loss.R)))
    
    loss.all = numeric(n)
    
    # threshold = x.sort[i]:
    # observations 1,...,i belong to the left regime
    # observations i+1,...,n belong to the right regime
    if (n > 1)
    {
      loss.all[1:(n-1)] = cum.L[1:(n-1)] + cum.R[2:n]
    }
    
    # threshold = largest x:
    # all observations belong to the left regime
    loss.all[n] = cum.L[n]
    
    # Convert RSS to MSE to match mse.loss.fixed()
    loss.all = loss.all / n
    
    temp.a = x.sort[which.min(loss.all)]
    
    # ========================================================
    # Coefficient update:
    # exactly the same as the original version
    # ========================================================
    
    temp.b = b.est(x, y, p, q,
                   temp.a, temp.a, Z)
    
    bL = temp.b$bL
    bR = temp.b$bR
    
    if (is.null(Z))
      bZ = NULL
    else
      bZ = temp.b$bZ
    
    temp.mse =
      mse.loss.fixed(
        x, y, p, q,
        bL, bR,
        temp.a,
        Z, bZ
      )
    
    # ========================================================
    # Original stopping rule
    # ========================================================
    
    if (previous.mse > temp.mse)
    {
      previous.mse = temp.mse
      a = temp.a
      test = 1
    }
    else
    {
      test = test + 1
    }
  }
  
  # ------------------------------------------------------------
  # Final coefficient update
  # ------------------------------------------------------------
  
  final.b = b.est(x, y, p, q,
                  a, a, Z)
  
  bL = final.b$bL
  bR = final.b$bR
  
  if (is.null(Z))
    bZ = NULL
  else
    bZ = final.b$bZ
  
  return(
    list(
      a = a,
      bL = bL,
      bR = bR,
      bZ = bZ
    )
  )
}
change.est.standard.grid.fast =
  function(x, y, p, q, Z=NULL,
           lower.quantile=0.10,
           upper.quantile=0.90)
  {
    n = length(y)
    
    # ------------------------------------------------------------
    # Sort observations according to x
    # ------------------------------------------------------------
    ord = order(x)
    
    x = x[ord]
    y = y[ord]
    
    if (!is.null(Z))
      Z = Z[ord, , drop=FALSE]
    
    # ------------------------------------------------------------
    # Construct polynomial bases
    # ------------------------------------------------------------
    XL = matrix(1, n, p+1)
    XR = matrix(1, n, q+1)
    
    if (p > 0)
    {
      for (j in 1:p)
        XL[,j+1] = x^j
    }
    
    if (q > 0)
    {
      for (j in 1:q)
        XR[,j+1] = x^j
    }
    
    # ------------------------------------------------------------
    # Full design vectors:
    #
    # left observation:
    # (XL, 0, Z)
    #
    # right observation:
    # (0, XR, Z)
    # ------------------------------------------------------------
    if (is.null(Z))
    {
      d = (p+1) + (q+1)
    }
    else
    {
      qz = ncol(Z)
      d = (p+1) + (q+1) + qz
    }
    
    # ------------------------------------------------------------
    # Candidate threshold range
    # ------------------------------------------------------------
    lower = as.numeric(quantile(x, lower.quantile))
    upper = as.numeric(quantile(x, upper.quantile))
    
    valid =
      which(x >= lower &
              x <= upper)
    
    # We need observations on both sides
    valid = valid[valid < n]
    
    # ------------------------------------------------------------
    # Construct row vectors corresponding to assigning each
    # observation to the LEFT or RIGHT regime
    # ------------------------------------------------------------
    if (is.null(Z))
    {
      H.L = cbind(
        XL,
        matrix(0, n, q+1)
      )
      
      H.R = cbind(
        matrix(0, n, p+1),
        XR
      )
    }
    else
    {
      H.L = cbind(
        XL,
        matrix(0, n, q+1),
        Z
      )
      
      H.R = cbind(
        matrix(0, n, p+1),
        XR,
        Z
      )
    }
    
    # ------------------------------------------------------------
    # Cumulative sufficient statistics
    #
    # For threshold x[i]:
    #
    # H'H =
    # sum_{j<=i} hL_j hL_j'
    # +
    # sum_{j>i} hR_j hR_j'
    #
    # H'y =
    # sum_{j<=i} hL_j y_j
    # +
    # sum_{j>i} hR_j y_j
    # ------------------------------------------------------------
    
    # Arrays for cumulative H'H
    cum.HH.L = array(0, dim=c(d,d,n))
    cum.HH.R = array(0, dim=c(d,d,n))
    
    # Matrices for cumulative H'y
    cum.Hy.L = matrix(0, n, d)
    cum.Hy.R = matrix(0, n, d)
    
    # ----- left cumulative statistics -----
    
    HH = matrix(0, d, d)
    Hy = rep(0, d)
    
    for (i in 1:n)
    {
      hi = H.L[i,]
      
      HH = HH + tcrossprod(hi)
      Hy = Hy + hi * y[i]
      
      cum.HH.L[,,i] = HH
      cum.Hy.L[i,] = Hy
    }
    
    # ----- right cumulative statistics -----
    
    HH = matrix(0, d, d)
    Hy = rep(0, d)
    
    for (i in n:1)
    {
      hi = H.R[i,]
      
      HH = HH + tcrossprod(hi)
      Hy = Hy + hi * y[i]
      
      cum.HH.R[,,i] = HH
      cum.Hy.R[i,] = Hy
    }
    
    # ------------------------------------------------------------
    # Total y'y does not depend on threshold
    # ------------------------------------------------------------
    yy = sum(y^2)
    
    # ------------------------------------------------------------
    # Profile RSS for each candidate threshold
    # ------------------------------------------------------------
    rss.all = rep(NA, length(valid))
    
    beta.all = matrix(NA,
                      nrow=length(valid),
                      ncol=d)
    
    solve.disturb = 1e-7
    
    for (k in seq_along(valid))
    {
      i = valid[k]
      
      # Left: observations 1,...,i
      # Right: observations i+1,...,n
      
      HH =
        cum.HH.L[,,i] +
        cum.HH.R[,,i+1]
      
      Hy =
        cum.Hy.L[i,] +
        cum.Hy.R[i+1,]
      
      # Same numerical stabilization as b.est()
      HH.temp =
        HH +
        solve.disturb * diag(d)
      
      beta =
        solve(HH.temp, Hy)
      
      beta.all[k,] = beta
      
      # --------------------------------------------------------
      # RSS = y'y - 2 beta'H'y + beta'H'H beta
      #
      # Use the unperturbed HH in the objective.
      # --------------------------------------------------------
      rss.all[k] =
        yy -
        2 * sum(beta * Hy) +
        as.numeric(
          t(beta) %*% HH %*% beta
        )
    }
    
    # ------------------------------------------------------------
    # Select the threshold minimizing profile RSS
    # ------------------------------------------------------------
    k.min = which.min(rss.all)
    
    i.min = valid[k.min]
    
    a = x[i.min]
    
    beta = beta.all[k.min,]
    
    # ------------------------------------------------------------
    # Extract coefficients
    # ------------------------------------------------------------
    bL =
      beta[1:(p+1)]
    
    bR =
      beta[(p+2):(p+q+2)]
    
    if (is.null(Z))
    {
      bZ = NULL
    }
    else
    {
      bZ =
        beta[(p+q+3):d]
    }
    
    return(
      list(
        a = a,
        bL = bL,
        bR = bR,
        bZ = bZ,
        mse = rss.all[k.min] / n,
        rss = rss.all[k.min]
      )
    )
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
# using the continuity-based estimator constructed from the M2 fit
break.est = function(x, y, p, q, M2.est, Z=NULL)
{
  # ------------------------------------------------------------
  # Step 1. Fix the unrestricted M2 coefficient estimates
  # ------------------------------------------------------------
  bL.2 = M2.est$bL
  bR.2 = M2.est$bR
  if (is.null(Z)) bZ.2 = NULL else bZ.2 = M2.est$bZ
  
  # Polynomial coefficients of
  # D(a) = h(a; bL.2) - g(a; bR.2).
  # polyroot() uses coefficients in increasing powers of a.
  deg = max(p, q)
  coef.D = rep(0, deg + 1)
  coef.D[1:(p + 1)] = coef.D[1:(p + 1)] + bL.2
  coef.D[1:(q + 1)] = coef.D[1:(q + 1)] - bR.2
  
  # Admissible threshold interval.
  # This can later be replaced by the exact A_0 used in the paper
  # if a more restrictive interval is desired.
  a.lower = min(x)
  a.upper = max(x)
  
  # ------------------------------------------------------------
  # Step 2. Construct stationary points of D(a)^2
  #
  # d D(a)^2 / da = 2 D(a) D'(a),
  #
  # so candidates may arise from either
  # D(a) = 0 or D'(a) = 0.
  # ------------------------------------------------------------
  tol.imag = 1e-8
  tol.root = 1e-8
  
  # Roots of D(a)
  roots.D = polyroot(coef.D)
  roots.D = Re(roots.D[abs(Im(roots.D)) < tol.imag])
  
  # Roots of D'(a)
  if (deg >= 1) {
    coef.D1 = (1:deg) * coef.D[2:(deg + 1)]
    
    # polyroot() should not be called on an identically zero
    # derivative polynomial.
    if (any(abs(coef.D1) > tol.root)) {
      roots.D1 = polyroot(coef.D1)
      roots.D1 = Re(roots.D1[abs(Im(roots.D1)) < tol.imag])
    } else {
      roots.D1 = numeric(0)
    }
  } else {
    roots.D1 = numeric(0)
  }
  
  cand.stat = c(roots.D, roots.D1)
  # Keep only finite stationary points in the admissible interval
  cand.stat = cand.stat[is.finite(cand.stat) & cand.stat > a.lower & cand.stat < a.upper]
  # Remove numerical duplicates
  if (length(cand.stat) > 0) {
    cand.stat = sort(cand.stat)
    cand.stat = cand.stat[c(TRUE, diff(cand.stat) > tol.root)]
  }
  # ------------------------------------------------------------
  # Step 3. Retain only local minima of D(a)^2
  # ------------------------------------------------------------
  D.fun = function(a)
  {
    powers = a^(0:deg)
    sum(coef.D * powers)
  }
  D2.fun = function(a)
  {
    D.fun(a)^2
  }
  cand.min = numeric(0)
  if (length(cand.stat) > 0) 
  {
    for (aa in cand.stat) 
    {
      # Use a small neighborhood determined by the distances
      # to the adjacent stationary points and boundaries.
      left.points = cand.stat[cand.stat < aa]
      right.points = cand.stat[cand.stat > aa]
      left.bound =
        if (length(left.points) == 0)
          a.lower
      else
        max(left.points)
      
      right.bound =
        if (length(right.points) == 0)
          a.upper
      else
        min(right.points)
      
      h = min(aa - left.bound, right.bound - aa) / 2
      if (is.finite(h) && h > 0) 
      {
        D2.left = D2.fun(aa - h)
        D2.mid = D2.fun(aa)
        D2.right = D2.fun(aa + h)
        if (D2.mid <= D2.left &&
            D2.mid <= D2.right) {
          cand.min = c(cand.min, aa)
        }
      }
    }
  }
  
  # ------------------------------------------------------------
  # Step 4. Select the M1 threshold by the unrestricted M2 loss
  #
  # Importantly, the M2 coefficients are kept fixed here.
  # ------------------------------------------------------------
  if (length(cand.min) > 0) 
  {
    compare.mse = sapply(cand.min,
      function(aa)
        mse.loss(
          x, y, p, q,
          bL.2, bR.2, aa,
          Z, bZ.2
        )
    )
    a = cand.min[which.min(compare.mse)]
  } else {
    # Numerical safeguard:
    # if no interior local minimum is detected, retain the
    # unrestricted M2 threshold.
    a = M2.est$a
  }
  # ------------------------------------------------------------
  # Step 5. Refit the regression coefficients at the selected
  # M1 threshold
  # ------------------------------------------------------------
  final.b = b.est(x, y, p, q, a, a, Z)
  bL = final.b$bL
  bR = final.b$bR
  if (is.null(Z))
    bZ = NULL
  else
    bZ = final.b$bZ
  # ------------------------------------------------------------
  # Step 6. Enforce continuity at the selected threshold
  #
  # Adjust the right intercept so that
  # h(a; bL) = g(a; bR).
  # ------------------------------------------------------------
  rL = sum(bL * a^(0:p))
  rR = sum(bR * a^(0:q))
  
  bR[1] = bR[1] + (rL - rR)
  return(list(a = a,bL = bL,bR = bR,bZ = bZ,candidates = cand.min))
}
mean.model = function(x,y,p,q,lambda,aL,aR,Z=NULL) 
{
	n = length(x)
	M2.est = change.est.grid(x,y,aL,aR,p,q,Z)
	M1.est = break.est(x,y,p,q,M2.est,Z)
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
	M2.est = change.est.grid(x,y,aL,aR,p,q,Z.mean)
	M1.est = break.est(x,y,p,q,M2.est,Z.mean)
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


