#'  @title r-value computation [Ver 1.5]
#' 
#'  @description The function computes r-values given two vectors of p-values from primary and
#'    follow-up studies. r-values assess the False Discovery Rate (FDR) of repilcability
#'    claims across the primary and follow-up studies. 
#'  
#'  @note The function is also available as a web applet:  http://www.math.tau.ac.il/~ruheller/App.html
#'  
#'  @usage r.value(p1, p2, m, c2 = 0.5, l00= 0)
#'            
#'  @param p1,p2 Numeric vectors (of the same length!) of the p-values from the follow-up
#'    study (p2) and the corresponding p-values from the primary study (p1).
#'  @param m     Number of features examined in the primary study.
#'  @param c2    Parameter for relative boost to the p-values from the primary study.
#'    0.5 (default) is recommended, since was observed in simulations to yield
#'    similar power to procedure with the optimal value (which is unknown for real data).
#'  @param l00   Lower bound of the fraction of features (out of m) with true null hypotheses in both studies.
#'    For example, for GWAS on the whole genome, the choice of 0.8 is conservative
#'    in typical applications.
#'  @param variation When 'use.m.star' is selected m* is used. m* is defined as follows:
#'   \eqn{m^*=m\sum_{i=1}^{m}\frac{1}{i}}{m*=m*sum(1/i)}.
#'   When 'use.t' is selected c1 is computed given the threshold tt.
#'   Both variations guarantee that the procedure that decleares all r-values below  q as replicability claims,
#'   controls the FDR at level q, for any type of dependency of the p-values in the primary study.
#'   default is 'none'.
#'  @param tt  The selection rule threshold for p-values from the primary study. must be supplied when
#'   variation 'use.t' is selected.
#'    
#'  @return vector of length of p2, containing the r-values.
#'  
#'  @example
#'  pv <- read.csv("http://www.math.tau.ac.il/~ruheller/Software/CrohnExample.csv")
#'  rv <- r.value(p1=pv$p1,p2=pv$p2,m=635547,c2=0.5,l00=0.8)
#'  rv2 <- r.value(p1=pv$p1,p2=pv$p2,m=635547,c2=0.5,l00=0.8,variation="use.t",tt=1e-5)
#'  
#'  @authors Ruth Heller, Shay Yaacoby (shay66@gmail.com)

r.value <- function (p1, p2, m, c2 = 0.5, l00= 0, variation = c("none","use.m.star","use.t"), tt = NULL, Q = 0.05)
{
  variation <- match.arg(variation)
  
  if (variation != "use.t" & !is.null(tt)) 
    warning("threshold tt is ignored")
  if (variation == "use.t" & !is.null(tt))
  {
    if (tt <= (1-c2)/(1-l00*(1-c2*Q)) * Q/m )
    {
      warning("since t < c(q)q/m, no modification to the original r-value computation was necessary (see section Derivation & Properties in the article)")
      variation <- "none"
    }
    if (tt >= (1-c2)/(1-l00*(1-c2*Q))*Q/(1+sum(1/(1:(m-1)))))
      warning("for the selected threshold t, the 'use.t' variation won't lead\nto more discoveries than the 'use.m.star' variation.")
    
  }
  if (variation == "use.m.star") m <- m*sum(1L/(1L:m)) 
  k <- length(p1) 
  #------- input checks: -------
  if (k==0 | length(p2)==0)
    stop("p-value vectors cannot have length zero")
  if (k!=length(p2))
    stop("p-value vectors must be of equallengths")
  if ((1<=c2)|(c2<=0))
    stop("c2 value should be in the interval (0,1)")
  if ((1<=l00)|(l00<0))
    stop("l00 value should be in the interval [0,1)")
  if (m<k)
    stop("Number of features in the primary stage (m) must be at least as large as the number of features followed-up.")
  if (any(is.na(c(p1,p2))))
    stop("NA's are not allowed as p-values")
  if (any(c(p1,p2)>1) | any(c(p1,p2)<=0))
    stop("p-values must be in the interval (0,1]")
  if (variation == "use.t" & is.null(tt)) 
    stop("specify threshold tt")
    
  #------- function definition: compute r-value of given value of x: -------
  r.value.x <- function (x) {
    
    c1 <- switch(variation,
                 none       = (1-c2)/(1-l00*(1-c2*x)),
                 use.m.star = (1-c2)/(1-l00*(1-c2*x)),
                 use.t = {
                   if (tt <= (1-c2)/(1-l00*(1-c2*x)) * x/m)
                     (1-c2)/(1-l00*(1-c2*x))
                   else
                   {
                     lower <- 1e-6 ; upper <- (1-c2)/(1-l00*(1-c2*x))
                     f <- function(a)
                     {
                       if (tt*m/(a*x) < 10)
                       {
                         a*(1+sum(1/(1:(ceiling(tt*m/(a*x)-1))))) - (1-c2)/(1-l00*(1-c2*x))
                       }
                       else
                       {
                         a*(1-digamma(1)+1/(2*ceiling(tt*m/(a*x)-1))+log(ceiling(tt*m/(a*x)-1))) - (1-c2)/(1-l00*(1-c2*x))
                       }
                     }
                     
                     c1_sol1 <- uniroot(f,c(lower,upper))
                     next_step <- tt*m/(ceiling(tt*m/(c1_sol1$root*x))*x)
                     if (f(next_step)<=0)
                     {
                       c1_sol2 <- uniroot(f,c(next_step,upper))
                       c1_sol2$root
                     }
                     else
                     {
                       c1_sol1$root 
                     }
                   }
                 })
    
    E   <- pmax(m*p1/c1, k*p2/c2)
    oe  <- order(E, decreasing =TRUE)
    oer <- order(oe)
    r   <- cummin((E/rank(E, ties.method= "max"))[oe])[oer]
    r   <- pmin(r,1)
    return(r)
  }
  
  rv <- rep(NA,length(p2))
  tol = min(c(p1,p2)[c(p1,p2)!=0], 0.0001)
  for (i in 1:length(p2))
  {
    aux <- function(x) return(r.value.x(x)[i]-x)
    if (aux(1)>=0) rv[i] <- 1
    else 
    {
      sol <- uniroot(aux,c(tol,1),tol = tol)
      rv[i] <- sol$root
    }
  }  
  return(rv)
}