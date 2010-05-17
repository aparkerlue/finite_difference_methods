is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

## Rotate `x' counterclockwise by 90 degrees.
rotate <- function(x, digits = NULL) {
  if (is.null(digits)) apply(x, 1, rev)
  else round(apply(x, 1, rev), digits=digits)
}
