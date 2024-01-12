#' Pre-computations for [sped]
#'
#' `compute_ephemera()` does data-independent pre-computations for [sped] and can speed up repeated applications
#'
#' The computations in [sped] rely on several matrices and vectors that are determined by the error density, spline space, and histogram bins, but do not depend on the data.
#' Computing these is the most time-intensive element of the process, so if the estimator will be applied several times to different data, but the same error density, spline space, and histogram bins (likely in simulations), gains can be had by pre-computing those matrices and vectors just one time.
#'
#' For comparison, the [sped] function internally uses `padding = 0.4`, and `perknot = 2`.
#' 
#' @references
#' 
#' Kent D, Ruppert D (2023). “Smoothness-Penalized Deconvolution (SPeD) of a Density Estimate.” _Journal of the American Statistical Association_, to appear. ISSN 0162-1459, \doi{10.1080/01621459.2023.2259028}
#'
#' @param gtwid Object of class [spedecon_gtwid][new_spedecon_gtwid] describing the density of \eqn{Z} in the model \eqn{Y = X + Z}
#' @param padding Support of spline space is extended by `padding/2` beyond the data on each side
#' @param hn Object of class `histogram` holding any histogram with the desired bins. The bins must be equally-spaced, i.e. `hn$equidist` must be `TRUE`, but otherwise only `hn$breaks` and `hn$mids` are used.
#' @param spline_dim Numeric integer, dimension of spline space
#' @param perknot Number of positivity constraints per knot
#'
#' @return Object of class `spedecon_ephemera`, a list containing the pre-computed values.
#'
#' @examples
#' alpha <- 1e-3; n <- 1e3; s <- 0.3
#' Y <- rgamma(n,5,2) + rnorm(n,0,s)
#' gtwid <- gaussian_gtwid(sd=s)
#' hn <- hist(Y,breaks="FD",plot=FALSE)
#' ephemera <- compute_ephemera(gtwid=gtwid,hn=hn,padding=0.4,spline_dim=30,perknot=2)
#' sol1 <- sped(Y,gtwid,1e-3,ephemera=ephemera) # fast
#' sol2 <- sped(Y,gtwid,1e-3) # slow
#' attr(sol1,"coef") - attr(sol2,"coef")
#' @export
compute_ephemera <- function(
                             gtwid,
                             hn,
                             padding,
                             spline_dim,
                             perknot=2
                             ) {
  stopifnot(is(hn,"histogram"))
  stopifnot(is.numeric(spline_dim))
  stopifnot(is.numeric(padding))
  stopifnot(is.numeric(perknot))
  stopifnot(padding >= 0)
  spline_dim <- as.integer(spline_dim)
  perknot <- as.integer(perknot)

  hnrange <- range(hn$breaks)

  ## Hard-coded order 4 spline
  knots <- seq(hnrange[1] - diff(hnrange)*padding/2,
               hnrange[2] + diff(hnrange)*padding/2,
               l=spline_dim + 4)

  basis <- new_spedecon_spline_basis(knots=knots, order=4L,
                                 evenly_spaced=TRUE,
                                 delta=diff(knots[1:2]))

  ## The integrate-to-one constraint
  Aint <- rep(1,dimension(basis))
  bint <- 1/basis$delta
  
  ## Positivity constraint
  mima <- range(basis$knots)
  posat <- seq(mima[1],mima[2],l=(length(basis$knots) - 1)*perknot+1)
  Apos <- design(basis,x=posat) ## Design matrix for the basis, evaluated at posat
  bpos <- rep(0,nrow(Apos))
  
  A <- cbind(Aint,t(Apos))
  bvec <- c(bint,bpos)
  meq <- 1

  result <- new_spedecon_ephemera(
    basis=basis,
    G=innprod(basis),
    P=innprod(basis,deriv=2),
    M=computeM(basis,gtwid),
    E=computeE(basis,gtwid,hn),
    A=A, bvec=bvec, meq=meq,
    hnbreaks=hn$breaks
  )

  return(result)
}

## Takes a spline basis object x
Nitwid <- function(t,x,i) {
  x$delta*exp(-1i*t*(x$knots[i] + x$delta*x$order/2))*sinc(x$delta*t/2)^x$order
}
## Takes a spline basis object x
NjconjNktwid <- function(t,x,j,k) {
  x$delta^2*exp(1i*t*x$delta*(k-j))*sinc(x$delta*t/2)^(2*x$order)
}

## No checking beyond types of coef and knots. We'll check in the validate_ function
#' Creates object of class `spedecon_spline_sped_fit`
#'
#' Internal use only. Constructor for class `spedecon_spline_sped_fit`
#'
#' The basic type of an object of type `spedecon_spline_sped_fit` is a function; one can therefore evaluate, plot, etc. and ignore the other attributes if desired.
#' The function is represented as a spline, and has attributes `coef` and `basis`, which represents the coefficients and basis respectively.
#' `coef` is a numeric vector, while `basis` is an object of class `spedecon_spline_basis`, which is essentially just a list holding the knots and order of the spline space.
#' A `spedecon_spline_sped_fit` object also has attributes `alpha` and `constraint` which record the penalty parameter and constraint method used for the fit.
#' 
#' @param coef Numeric vector; the coefficients of the spline.
#' @param basis Object of class `spedecon_spline_basis`, representing a basis for the spline space.
#' @param alpha Positive numeric, the alpha that was used for the fit.
#' @param constraint The type of constraint used.
#'
#' @return Object of class `spedecon_spline_sped_fit`
new_spedecon_spline_sped_fit <- function(coef, basis, alpha, constraint) {
  stopifnot(is.numeric(coef))
  knots <- basis$knots
  order <- basis$order
  structure(\(x, deriv=0) drop(splines::splineDesign(x=x,knots=knots,ord=order,
                                                deriv=deriv, outer.ok=TRUE)%*%coef),
            coef=coef,
            basis=basis,
            alpha=alpha,
            constraint=constraint,
            class=c("spedecon_spline_sped_fit","spedecon_spline")
            )
}

## No checking beyond types of coef and knots. We'll check in the validate_ function
new_sped_spline <- function(coef, basis) {
  stopifnot(is.numeric(coef))
  knots <- basis$knots
  order <- basis$order
  structure(\(x, deriv=0) drop(splines::splineDesign(x=x,knots=knots,ord=order,
                                                deriv=deriv, outer.ok=TRUE)%*%coef),
            coef=coef,
            basis=basis,
            class="spedecon_spline"
            )
}

## Check if knots are sorted, the number of knots is higher than the
## order, etc.
validate_spedecon_spline <- function(x) {

}

#' @export
plot.spedecon_spline <- function(x,...,deriv=0) {
  stopifnot(is(x,"spedecon_spline"))
  spl <- \(t) x(t,deriv=deriv)
  curve(spl,min(attr(x,"basis")$knots),max(attr(x,"basis")$knots),...)
}

#' Creates object of class `spedecon_gtwid`
#'
#' Constructor for class `spedecon_gtwid`. Use helper functions [gaussian_gtwid()], [laplace_gtwid()], and [uniform_gtwid()] instead whenever possible.
#'
#' The `spedecon_gtwid` class is meant to represent the Fourier transform of a probability density.
#' The basic type is a function.
#' It also has a `family` attribute which can hold the name and parameters of the family of distributions.
#' 
#' @param gtwid Function representing the Fourier transform
#' @param family List with at least one entry `family[["family"]]` naming the family of distributions, and possibly other entries stating the values of the parameters in that family.
#'
#' @seealso Use [gaussian_gtwid()], [laplace_gtwid()], or [uniform_gtwid()] instead whenever possible.
#'
#' @return Object of class `spedecon_gtwid`
#' @export
new_spedecon_gtwid <- function(gtwid,family) {
  structure(\(t) gtwid(t),
            family=family,
            class="spedecon_gtwid"
            )
}

#' Fourier transform of Gaussian density
#'
#' Returns a `spedecon_gtwid` object representing the Fourier transform of a mean-zero Gaussian density
#'
#' @param sd Standard deviation
#'
#' @return Object of class `spedecon_gtwid`
#'
#' @examples
#' gtwid <- gaussian_gtwid(sd = 1)
#' @export
gaussian_gtwid <- function(sd) {
  new_spedecon_gtwid(gtwid=function(t) exp(-sd^2*t^2/2),
                 family=list(family="Gaussian",sd=sd))
}

## Scale b
#' Fourier transform of Laplace density
#'
#' Returns a `spedecon_gtwid` object representing the Fourier transform of a mean-zero Laplace density with scale `b`
#'
#' @param b Scale parameter
#'
#' @return Object of class `spedecon_gtwid`
#'
#' @examples
#' gtwid <- laplace_gtwid(b = 1)
#' @export
laplace_gtwid <- function(b) {
  new_spedecon_gtwid(gtwid=function(t) 1/(1+b^2*t^2),
                 family=list(family="Laplace",b=b))
}

## Half-width a
#' Fourier transform of Uniform density
#'
#' Returns a `spedecon_gtwid` object representing the Fourier transform of a Uniform\[-a,a\] density
#'
#' @param a Half-width
#'
#' @return Object of class `spedecon_gtwid`
#'
#' @examples
#' gtwid <- uniform_gtwid(a = 1)
#' @export
uniform_gtwid <- function(a) {
  new_spedecon_gtwid(gtwid=function(t) sinc(a*t),
                 family=list(family="Uniform",a=a))
}

## delta is the maximum knot spacing. If `evenly_spaced` is TRUE, then
## delta is also the spacing for all knots.
## Assume knots are sorted already.
## No non-trivial checking.
new_spedecon_spline_basis <- function(knots,order = 4,evenly_spaced,delta) {
  stopifnot(is.numeric(knots))
  stopifnot(is.integer(order))
  stopifnot(is.logical(evenly_spaced))
  if( missing(delta) )
    delta <- max(diff(knots))
  structure(list(knots=knots,
                 order=order,
                 delta=delta,
                 evenly_spaced=evenly_spaced),
            class="spedecon_spline_basis"
            )
}

## Plot the B-spline basis functions.
#' @export
plot.spedecon_spline_basis <- function(x,xmin,xmax,n=1000,...) {
  if( missing(xmin) ) xmin <- min(x$knots)
  if( missing(xmax) ) xmax <- max(x$knots)
  xs <- seq(xmin,xmax,l=n)
  matplot(xs,splines::splineDesign(knots=x$knots,x=xs,ord=x$order,outer.ok=TRUE),
          xlab="x", ylab=expression(N[i](x)),
          main="B-spline basis for spline space",
          type='l',
          lty=1,...
          )
}

dimension <- function(x) UseMethod("dimension")

## The dimension is the number of knots minus the order.
#' @export
dimension.spedecon_spline_basis <- function(x) {
  stopifnot(is(x,"spedecon_spline_basis"))
  return(length(x$knots) - x$ord)
}

design <- function(basis,x,deriv=0) UseMethod("design")

#' @export
design.spedecon_spline_basis <- function(basis,x,deriv=0) {
  stopifnot(is(basis,"spedecon_spline_basis"))
  stopifnot(is.numeric(x))
  splines::splineDesign(x=x,knots=basis$knots,ord=basis$order,
                        deriv=deriv, outer.ok=TRUE)
}

new_spedecon_ephemera <- function(G,P,M,E,A,bvec,meq,basis,hnbreaks) {
  stopifnot(is(basis,"spedecon_spline_basis"))
  stopifnot(is.matrix(G))
  stopifnot(is.matrix(P))
  stopifnot(is.matrix(M))
  stopifnot(is.matrix(E))
  stopifnot(is.matrix(A))
  stopifnot(is.numeric(bvec))
  stopifnot(is.numeric(meq))
  stopifnot(is.numeric(hnbreaks))

  structure(
    list(G=G,P=P,M=M,E=E,A=A,bvec=bvec,meq=meq,basis=basis),
    hnbreaks=hnbreaks,
    class="spedecon_ephemera"
  )
}

sinc <- function(x) ifelse(x!=0,1*sin(x)/x,1)
