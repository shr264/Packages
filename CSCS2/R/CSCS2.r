CSCS2 <- function( Y, lambda, L=NULL, maxitr=100, tol=1e-4, warmstart=FALSE) {
    
  require(lassoshooting)
  n = nrow(Y)
  p = ncol(Y)

  if (is.null(L)) L = diag(p)

  S = (t(Y)%*%Y)/n

  itr_log = eps_log = NULL

  L[1, 1] = 1/sqrt(S[1,1])
  for (k in 2:p){


    nuk_old = nuk_new = c(rep(0, k-1), 1) 
    r = 0

    repeat {

      r = r + 1
      
      km1_ = 1:(k-1)

      hk = lassoshooting(XtX    =  S[km1_, km1_, drop=FALSE], 
			 Xty    = -nuk_old[k] * S[km1_, k],
			 lambda =  0.5*lambda)
      nuk_new[km1_] = hk$coefficients

      sumterm = sum(nuk_new[km1_] * S[k, km1_])
      nuk_new[k] = (-sumterm + sqrt(sumterm^2 + 4*S[k, k]))/(2*S[k, k])

      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
	L[k, 1:k] = nuk_new
	eps_log = c(eps_log, maxdiff)
	itr_log = c(itr_log, r)
	break
      } else {
	nuk_old = nuk_new
      }
    }
  }

  list(L=L, itr=itr_log, eps=eps_log)

}




