#' Smoothness-Penalized Deconvolution
#'
#' `sped()` computes the Smoothness-Penalized Deconvolution estimate on the provided data and error distribution
#'
#' This function computes the "Smoothness-Penalized Deconvolution" (SPeD) estimate of a density under additive measurement error.
#' The essential inputs to the function are the data `Y`, the Fourier transform `gtwid` of the error density, and the penalty parameter `alpha`; more details follow here, but for a full description of the estimator please consult Kent and Ruppert (2023).
#' 
#' The data model is that we observe an iid sample distributed like \eqn{Y = X + Z}, with \eqn{Z} an error independent of \eqn{X}.
#' We wish to estimate the density \eqn{f(x)} of \eqn{X}.
#' It is assumed that we know the probability density of the errors \eqn{Z}, call it \eqn{g(z)}.
#' 
#' The estimator begins with a density estimate \eqn{h_n(y)}{hₙ(y)} of the density of \eqn{Y}, and minimizes the objective function
#' \deqn{\|g*v - h_n\|^2 + \alpha \|v^{(2)}\|^2}{‖g*v - hₙ‖² + α‖v⁽²⁾‖²}
#' in \eqn{v}, with \eqn{v} ranging over a space of cubic splines with equally-spaced knots; the dimension of this space can be adjusted with the argument `spline_dim`.
#' The SPeD estimate is not naturally a pdf, so it must be constrained.
#' When `constraint = "constrainedQP"`, the constraint is imposed directly into quadratic program minimizing the objective; when `constraint = "projection"`, the unconstrained estimate is computed and then projected onto the space of pdfs.
#' The preliminary density estimate \eqn{h_n}{hₙ} is computed internally as a histogram using Freedman-Diaconis choice of bin width, but a user-supplied histogram computed with [hist()] may be provided via the `hn` argument.
#' 
#' The computations require the *Fourier transform* \eqn{\tilde g(t)}{gtwid} of the probability density, and this must be supplied as an object of type `spedecon_gtwid`, which can be produced for common error densities using the helper functions [gaussian_gtwid()], [laplace_gtwid()], and [uniform_gtwid()].
#'
#' If the estimator will be re-computed many times for many realizations of data, substantial time can be saved by pre-computing all the auxiliary matrices and vectors one time, and supplying them through the `ephemera` argument.
#' This can be done whenever the repeated computations all use the same error density, same histogram bins, and same spline space, as those are what define the required matrices and vectors.
#' A helper function [compute_ephemera()] is provided to pre-compute these.
#'
#' @references
#' 
#' Kent D, Ruppert D (2023). “Smoothness-Penalized Deconvolution (SPeD) of a Density Estimate.” _Journal of the American Statistical Association_, to appear. ISSN 0162-1459, \doi{10.1080/01621459.2023.2259028}
#'
#' @param Y Numeric vector of data from the model \eqn{Y = X + Z}
#' @param gtwid Object of class [spedecon_gtwid][new_spedecon_gtwid] describing the density of \eqn{Z} in the model \eqn{Y = X + Z}. It should almost always be created by one of the helper functions [gaussian_gtwid()], [laplace_gtwid()], or [uniform_gtwid()].
#' @param alpha Positive numeric penalty parameter
#' @param constraint String, controls whether and how the solution is constrained to be a pdf. One of `"constrainedQP"`, `"projection"`, or `"unconstrained"` for constrained quadratic program, metric projection, or unconstrained, respectively
#' @param spline_dim Numeric integer, dimension of spline space
#' @param hn (optional) Object of class `histogram` holding pre-computed histogram computed from the data \eqn{Y}
#' @param ephemera (optional) Object of class `spedecon_ephemera` holding pre-computed computational bits
#' @param ... (optional) Other arguments
#'
#' @return Object of class [spedecon_spline_sped_fit][new_spedecon_spline_sped_fit]
#'
#' @examples
#' alpha <- 1e-3
#' n <- 1e3; s <- 0.3
#' Y <- rgamma(n,5,2) + rnorm(n,0,s) # Data, contaminated with Gaussian errors
#' sol <- sped(Y,gtwid=gaussian_gtwid(sd=s),1e-3)
#' plot(sol,n=1e3) # Plot the resulting estimate
#' curve(dgamma(x,5,2),col=2,n=1e3,add=TRUE) # The target density f() of X
#' sol(c(2,3,4)) # We can evaluate sol; it is a function
#' @export
sped <- function(
                 Y,
                 gtwid,
                 alpha,
                 constraint="constrainedQP",
                 spline_dim=30,
                 hn=NULL,
                 ephemera=NULL,
                 ...
                 ) {

  ## Process the arguments
  stopifnot("`Y` must be numeric"=is.numeric(Y))
  stopifnot("`alpha` must be positive"=alpha > 0)
  stopifnot("`gtwid` must be of class `spedecon_gtwid`. See helper function `?spedecon_gtwid`"=
              is(gtwid,"spedecon_gtwid"))
  stopifnot(
    '`constraint` must be one of "unconstrained", "projection", or "constrainedQP"'=
      (constraint %in% c("unconstrained", "projection", "constrainedQP"))
  )
  stopifnot(is.numeric(spline_dim))
  spline_dim <- as.integer(spline_dim)

  ## Compute the ephemera and other utility values
  if( is.null(hn) ) {
    hn <- hist(Y,breaks="FD",plot=FALSE)
    n <- length(Y)
  } else {
    # Make sure user-supplied hn is histogram
    stopifnot(is(hn,"histogram")) 
    n <- sum(hn$counts)
  }
  
  if( is.null(ephemera) ) {
    ephemera <- compute_ephemera(gtwid,hn,padding=0.4,spline_dim=spline_dim,perknot=2)
  } else {
    # Make sure user-supplied ephemera matches hn
    stopifnot("Ephemera must be computed using the same histogram bins as hn"=
                hn$breaks == ephemera$hnbreaks)
  } 
  
  # Compute the actual d vector
  d <- ephemera$E%*%(hn$counts/n)

  ## Compute the result
  if( constraint == "unconstrained" || constraint == "projection" ) {
    ## Need unconstrained theta for projection method too
    thetaunc <- solve(ephemera$M + alpha*ephemera$P, d)
    if( constraint == "projection" ) {
      theta <- quadprog::solve.QP(ephemera$G,
                                  t(thetaunc)%*%ephemera$G,
                                  Amat=ephemera$A,bvec=ephemera$bvec,
                                  meq=ephemera$meq)$solution
    } else {
      theta <- thetaunc
    }
  } else if( constraint == "constrainedQP" ) {
    theta <- quadprog::solve.QP(ephemera$M + alpha*ephemera$P,
                                d,
                                Amat=ephemera$A,bvec=ephemera$bvec,
                                meq=ephemera$meq)$solution
  }
  
  ## Construct spline function from coefficients
  fnalpha <- new_spedecon_spline_sped_fit(coef = theta,
                                      basis = ephemera$basis,
                                      alpha = alpha,
                                      constraint = constraint)

  

  return(fnalpha)
}
