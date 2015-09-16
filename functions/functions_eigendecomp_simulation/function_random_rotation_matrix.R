# function to generate random roations of size n
f.random_rotation_incl_flip <- function(n){
  QR <- qr(matrix(rnorm(n^2), ncol=n)) # A = QR
  M <- qr.Q(QR) %*%  diag(sign(diag(qr.R(QR)))) # diag(R) > 0
  return(M)
}