ccdr_cust <- function(X, B=NULL, lambda, maxitr=10, tol=10^(-4), verbose = FALSE) {
    n = dim(X)[1]
    p = dim(X)[2]
    S = (1/n)*crossprod(X)
    if(is.null(B)){
        B = Bold = matrix(1,p,p)
    } else {Bold = B}
    .Call('_partial_ccdr_cust', PACKAGE = 'partial', S, B, lambda, maxitr, tol, verbose)
}

partial_cust1 <- function(X, B=NULL, lambda, m, maxitr=10, tol=10^(-4),verbose = FALSE) {
    n = dim(X)[1]
    p = dim(X)[2]
    S = (1/n)*crossprod(X)
    if(is.null(B)){
        B = Bold = matrix(1,p,p)
    } else {Bold = B}
    .Call('_partial_partial_cust1', PACKAGE = 'partial', S, B, lambda, m, maxitr, tol,verbose)
}

partial_cust2 <- function(X, B=NULL, lambda, m, maxitr=10, tol=10^(-4),verbose = FALSE) {
    n = dim(X)[1]
    p = dim(X)[2]
    S = (1/n)*crossprod(X)
    if(is.null(B)){
        B = Bold = matrix(1,p,p)
    } else {Bold = B}
    .Call('_partial_partial_cust2', PACKAGE = 'partial', S, B, lambda, m, maxitr, tol,verbose)
}

