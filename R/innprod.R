## Evaluates all the inner-products ∫NᵢNⱼ of normalized B-spline
## basis functions of the supplied spline space
innprod <- function(x,deriv=0) {
  stopifnot(is(x,"spedecon_spline_basis"))
  stopifnot("Only implemented for deriv=0 and deriv=2"=
              (deriv == 0 || deriv == 2))
  stopifnot("Second-derivative inner-product matrix only impemented for order 4"=
              (deriv != 2 || x$order == 4))
  stopifnot("Only implemented for order 5 and smaller"=x$order < 6)
  stopifnot("Only implemented for evenly-spaced knots"=x$evenly_spaced)

  nb <- dimension(x)
  delta <- x$delta
  G <- matrix(0,nrow=nb,ncol=nb)
  cG <- .col(dim(G))
  rG <- .row(dim(G))

  if( deriv == 0 ) {
    if( x$order == 1 ) {
      G[cG == rG] <- delta
    } else if( x$order == 2 ) {
      G[abs(cG - rG) == 0] <- delta*4/6
      G[abs(cG - rG) == 1] <- delta*1/6
    } else if( x$order == 3 ) {
      G[abs(cG - rG) == 0] <- delta*66/120
      G[abs(cG - rG) == 1] <- delta*26/120
      G[abs(cG - rG) == 2] <- delta*1 /120
    } else if( x$order == 4 ) {
      G[abs(cG - rG) == 0] <- delta*2416/5040
      G[abs(cG - rG) == 1] <- delta*1191/5040
      G[abs(cG - rG) == 2] <- delta*120 /5040
      G[abs(cG - rG) == 3] <- delta*1   /5040
    } else if( x$order == 5 ) {
      G[abs(cG - rG) == 0] <- delta*156190/362880
      G[abs(cG - rG) == 1] <- delta*88234 /362880
      G[abs(cG - rG) == 2] <- delta*14608 /362880
      G[abs(cG - rG) == 3] <- delta*502   /362880
      G[abs(cG - rG) == 4] <- delta*1     /362880
    } else {
      stop("Shouldn't be possible to get here")
    }
  } else if( deriv == 2 ) {
    ## These are computed explicitly from the form of the canonical B-spline basis function.
    premult <- (1/delta^3)
    G[abs(cG - rG) == 0] <- premult*(8/3)
    G[abs(cG - rG) == 1] <- premult*(-3/2)
    G[abs(cG - rG) == 3] <- premult*(1/6)
  }

  return(G)
}

computeM <- function(x,gtwid) {
    stopifnot(is(x,"spedecon_spline_basis"))
    stopifnot(is(gtwid,"spedecon_gtwid"))

    nb <- dimension(x)
    M <- matrix(0,nrow=nb,ncol=nb)
    cM <- .col(dim(M))
    rM <- .row(dim(M))
    for(i in 0:(nb - 1)) {
      intgrnd <- function(t) Re(Mod(gtwid(t))^2*NjconjNktwid(t,x,0,i))
      M[abs(cM - rM) == i] <- integrate(intgrnd,-Inf,Inf)$value/(2*pi)
    }
    return(M)
}

computeE <- function(basis,gtwid,hn) {
  stopifnot(is(hn,"histogram"))
  stopifnot(is(basis,"spedecon_spline_basis"))
  stopifnot(is(gtwid,"spedecon_gtwid"))
  stopifnot(hn$equidist)

  k <- length(hn$mids)
  delta <- diff(hn$mids)[1]
  E <- matrix(NA,nrow=dimension(basis),ncol=k)

  for(j in 1:dimension(basis)) {
    for(ell in 1:k) {
      intgnd <- \(t) Re(gtwid(t)*Nitwid(t,basis,j)*exp(1i*t*hn$mids[ell])*sinc(delta*t/2)/(2*pi))
      E[j,ell] <- integrate( intgnd, -Inf, Inf)$value
    }
  }

  return(E)
}
