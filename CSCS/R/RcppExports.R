CSCS <- function(Y = NULL, lambda = NULL, newL = NULL , maxitr = NULL, tol = NULL, returnOmega = FALSE) {
    if(is.null(Y)){stop('Must include data matrix')}
    if(is.null(lambda)){stop('Must include penalty parameter')}
    if(is.null(newL)){newL = diag(p)}
    if(is.null(maxitr)){maxitr = 100}
    if(is.null(tol)){tol = 10^(-4)}
    n = dim(Y)[1]
    p = dim(Y)[2]
    S = crossprod(Y)/n
    L = .Call('CSCS_CSCS', PACKAGE = 'CSCS', S, newL, lambda, maxitr, tol)
    if(returnOmega){out = crossprod(L)} else {out = L}
    out
}

