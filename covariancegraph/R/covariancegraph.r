require(Rcpp)
require(fillgraph)

cmpGraph <- function(amat){
    g <- 1*!amat
    diag(g) <- 0
    g
}

icf <- function(amat, Y, start=NULL, tol = 1e-06){
  n = dim(Y)[1]
  p = dim(Y)[2]
  S = cov(Y)
  bi.graph = amat - diag(p)
  dimnames(bi.graph)[[1]] = as.list(seq(1:p))
  dimnames(bi.graph)[[2]] = as.list(seq(1:p))
  dimnames(S)[[1]] = as.list(seq(1:p))
  dimnames(S)[[2]] = as.list(seq(1:p))
  if(!is.matrix(S)){
    stop("Second argument is not a matrix!")
  }
  if(dim(S)[1]!=dim(S)[2]){
    stop("Second argument is not a square matrix!")
  }
  if(min(eigen(S)[[1]])<=0){
    stop("Second argument is not a positive definite matrix!")
  }

  p <- nrow(S)
  i <- 0

  ## prep spouses and non-spouses

  pa.each.node <-function(amat){
    ## List of the parents of each node.
    ## If amat is symmetric it returns the boundaries.
    p <- nrow(amat)
    b <- vector(p, mode="list")
    ip <- 1:p
    for(i in 1:p)
      b[[i]] <- ip[amat[,i]==1]
    b
  }
  spo <- pa.each.node(bi.graph)
  nsp <- pa.each.node(cmpGraph(bi.graph))
  number.spouses <- unlist(lapply(spo, length))
  no.spouses <- (1:p)[number.spouses==0]
  all.spouses <- (1:p)[number.spouses==(p-1)]
  nontrivial.vertices <- setdiff((1:p), no.spouses)

  if(length(nontrivial.vertices)==0){
    if(p==1){
      Sigma <- S
    }
    else{
      Sigma <- diag(diag(S))
      dimnames(Sigma) <- dimnames(S)
    }
    return(list(Sigmahat=Sigma, iterations=1))
  }

  if(is.null(start)){
    Sigma <- as.matrix(diag(diag(S))) # starting value
  }
  else{
    Sigma <- start
  }

  repeat{
    i <- i+1
    Sigma.old <- Sigma
    for(v in nontrivial.vertices){
      if(is.element(v, all.spouses)){
        B <- S[v,-v]%*%solve(S[-v,-v])
        lambda <- S[v,v]-B%*%S[-v,v]
      }
      else{
        B.spo.nsp <-
          Sigma[spo[[v]],nsp[[v]]]%*%solve(Sigma[nsp[[v]],nsp[[v]]])
        YZ <- S[v,spo[[v]]]-S[v,nsp[[v]]]%*%t(B.spo.nsp)
        B.spo <-
          YZ %*%
          solve( S[spo[[v]],spo[[v]]]
                 -S[spo[[v]],nsp[[v]]]%*%t(B.spo.nsp)
                 -B.spo.nsp%*%S[nsp[[v]],spo[[v]]]
                 +B.spo.nsp%*%S[nsp[[v]],nsp[[v]]]%*%t(B.spo.nsp) )
        lambda <- S[v,v]-B.spo%*%t(YZ)
        B.nsp <- -B.spo%*%B.spo.nsp
        B <- rep(0, p)
        B[spo[[v]]] <- B.spo
        B[nsp[[v]]] <- B.nsp
        B <- B[-v]
      }
      ## here I can improve by only using B[spo[[v]]]!
      Sigma[v,-v] <- B%*%Sigma[-v,-v]
      Sigma[v,nsp[[v]]] <- 0
      Sigma[-v,v] <- t(Sigma[v,-v])
      Sigma[v,v] <- lambda + B%*%Sigma[-v,v]
    }
    if(sum(abs(Sigma.old-Sigma)) < tol) break
  }
  dimnames(Sigma) <- dimnames(S)
  return(list(Shat=Sigma, iterations=i))
}

covariancegraph <- function(Y,G){

    n = dim(Y)[1]
    p = dim(Y)[2]
    S = cov(Y)

    out = chol(as.spam(S*(G) + p*diag(x = max(abs(S)), p)), pivot ="RCM")
    Q = ordering(out)
    S = S[Q,Q]
    G = G[Q,Q]

    Gfill = gNM1(G,p)

    L = t(chol(S))

    (E = matrix(0,0,2))
    for(i in 2:p){
        (zero = which(Gfill[i,1:i]==0))
        if(length(zero)>0){
            (flintrack = cbind(rep(i,length(zero)),zero))
            (E = rbind(E,flintrack))}
    }

    (a <- length(E[,1]))


    cat('...Algorithm 1')
    for(i in 1:a){
        if(E[i,2]>1){
            (L[E[i,1],E[i,2]] = -sum(L[E[i,1],1:(E[i,2]-1)]*L[E[i,2],1:(E[i,2]-1)]))
            (L[E[i,1],E[i,2]] = (L[E[i,1],E[i,2]])/(L[E[i,2],E[i,2]]))
        }
        else{
            L[E[i,1],E[i,2]] = 0}
    }


    (omegahat = L%*%t(L))
    omegahat = omegahat[invPerm(Q),invPerm(Q)]
    G = G[invPerm(Q),invPerm(Q)]
    Gfill = Gfill[invPerm(Q),invPerm(Q)]
    return(list(Shat=omegahat))
}



