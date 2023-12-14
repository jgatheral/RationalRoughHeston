#####################################################################################
# Jim Gatheral, 2023
# Various rational approximations to the rough Heston solution
#####################################################################################

########################################################################
# Pade approximations to h(a,t)
########################################################################

########################################################################

h.Pade22 <- function(params)function(a,t){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  b1 <- -a*(a+1i)/2 * 1/gamma(1+al)
  b2 <- -b1 * lamTilde * nu * gamma(1+al)/gamma(1+2*al)               
  
  g0 <- rm/nu
  g1 <- -1/aa*1/gamma(1-al)*g0/nu
  
  den <- g0^2 +b1*g1
  
  p1 <- b1
  p2 <- (b1^2*g0+b2*g0^2)/den
  q1 <- (b1*g0-b2*g1)/den
  q2 <- (b1^2+b2*g0)/den
  
  y <- t^al
  
  h.pade <- (p1*y + p2*y^2)/(1 + q1*y + q2*y^2)
  return(h.pade)
}

########################################################################

h.Pade33 <- function(params)function(a,x){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  b1 <- -a*(a+1i)/2 * 1/gamma(1+al)
  b2 <- -b1 * lamTilde * nu * gamma(1+al)/gamma(1+2*al)               
  b3 <- (-b2 * lamTilde * nu  + nu^2 * b1^2/2) * gamma(1+2*al)/gamma(1+3*al)
  b4 <- (-b3 * lamTilde *nu + nu^2 * b1*b2)* gamma(1+3*al)/gamma(1+4*al)
  
  g0 <- rm/nu
  g1 <- -1/(aa*nu)*1/gamma(1-al)*g0
  g2 <- -1/(aa*nu)*(gamma(1-al)/gamma(1-2*al)*g1 -1/2*nu^2*g1*g1)
  g3 <- -1/(aa*nu)*(gamma(1-2*al)/gamma(1-3*al)*g2 -nu^2 * g1*g2) 
  
  den <- g0^3 +2*b1*g0*g1-b2*g1^2+b1^2*g2+b2*g0*g2
  
  p1 <- b1
  p2 <- (b1^2*g0^2 + b2*g0^3 + b1^3*g1 + b1*b2*g0*g1 - b2^2*g1^2 +b1*b3*g1^2 +b2^2*g0*g2 - b1*b3*g0*g2)/den
  q1 <- (b1*g0^2 + b1^2*g1 - b2*g0*g1 + b3*g1^2 - b1*b2*g2 -b3*g0*g2)/den
  q2 <- (b1^2*g0 + b2*g0^2 - b1*b2*g1 - b3*g0*g1 + b2^2*g2 - b1*b3*g2)/den
  q3 <- (b1^3 + 2*b1*b2*g0 + b3*g0^2 -b2^2*g1 +b1*b3*g1 )/den
  p3 <- g0*q3
  
  y <- x^al
  
  h.pade <- (p1*y + p2*y^2 + p3*y^3)/(1 + q1*y + q2*y^2 + q3*y^3)
  return(h.pade)
}

########################################################################

h.Pade44 <- function (params)function(a, x) {
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  b1 <- -a*(a+1i)/2 * 1/gamma(1+al)
  b2 <- -b1 * lamTilde * nu * gamma(1+al)/gamma(1+2*al)               
  b3 <- (-b2 * lamTilde * nu  + nu^2 * b1^2/2) * gamma(1+2*al)/gamma(1+3*al)
  b4 <- (-b3 * lamTilde *nu + nu^2 * b1*b2)* gamma(1+3*al)/gamma(1+4*al)
  
  g0 <- rm/nu
  g1 <- -1/(aa*nu)*1/gamma(1-al)*g0
  g2 <- -1/(aa*nu)*(gamma(1-al)/gamma(1-2*al)*g1 -1/2*nu^2*g1*g1)
  g3 <- -1/(aa*nu)*(gamma(1-2*al)/gamma(1-3*al)*g2 -nu^2 * g1*g2) 
  
   den <- (g0^4 + 3 * b1 * g0^2 * g1 + b1^2 * g1^2 - 2 * b2 * 
            g0 * g1^2 + b3 * g1^3 + 2 * b1^2 * g0 * g2 + 2 * b2 * 
            g0^2 * g2 - 2 * b1 * b2 * g1 * g2 - 2 * b3 * g0 * g1 * 
            g2 + b2^2 * g2^2 - b1 * b3 * g2^2 + b1^3 * g3 + 2 * b1 * 
            b2 * g0 * g3 + b3 * g0^2 * g3 - b2^2 * g1 * g3 + b1 * 
            b3 * g1 * g3)
  p1 <- b1
  p2 <- (b1^2 * g0^3 + b2 * g0^4 + 2 * b1^3 * g0 * g1 + 2 * 
           b1 * b2 * g0^2 * g1 - b1^2 * b2 * g1^2 - 2 * b2^2 * g0 * 
           g1^2 + b1 * b3 * g0 * g1^2 + b2 * b3 * g1^3 - b1 * b4 * 
           g1^3 + b1^4 * g2 + 2 * b1^2 * b2 * g0 * g2 + 2 * b2^2 * 
           g0^2 * g2 - b1 * b3 * g0^2 * g2 - b1 * b2^2 * g1 * g2 + 
           b1^2 * b3 * g1 * g2 - 2 * b2 * b3 * g0 * g1 * g2 + 2 * 
           b1 * b4 * g0 * g1 * g2 + b2^3 * g2^2 - 2 * b1 * b2 * 
           b3 * g2^2 + b1^2 * b4 * g2^2 + b1 * b2^2 * g0 * g3 - 
           b1^2 * b3 * g0 * g3 + b2 * b3 * g0^2 * g3 - b1 * b4 * 
           g0^2 * g3 - b2^3 * g1 * g3 + 2 * b1 * b2 * b3 * g1 * 
           g3 - b1^2 * b4 * g1 * g3)/den
  p3 <- (b1^3 * g0^2 + 2 * b1 * b2 * g0^3 + b3 * g0^4 + b1^4 * 
           g1 + 2 * b1^2 * b2 * g0 * g1 - b2^2 * g0^2 * g1 + 2 * 
           b1 * b3 * g0^2 * g1 - 2 * b1 * b2^2 * g1^2 + 2 * b1^2 * 
           b3 * g1^2 - b2 * b3 * g0 * g1^2 + b1 * b4 * g0 * g1^2 + 
           b3^2 * g1^3 - b2 * b4 * g1^3 + b1 * b2^2 * g0 * g2 - 
           b1^2 * b3 * g0 * g2 + b2 * b3 * g0^2 * g2 - b1 * b4 * 
           g0^2 * g2 + b2^3 * g1 * g2 - 2 * b1 * b2 * b3 * g1 * 
           g2 + b1^2 * b4 * g1 * g2 - 2 * b3^2 * g0 * g1 * g2 + 
           2 * b2 * b4 * g0 * g1 * g2 - b2^3 * g0 * g3 + 2 * b1 * 
           b2 * b3 * g0 * g3 - b1^2 * b4 * g0 * g3 + b3^2 * g0^2 * 
           g3 - b2 * b4 * g0^2 * g3)/den
  q1 <- (b1 * g0^3 + 2 * b1^2 * g0 * g1 - b2 * g0^2 * g1 - 
           2 * b1 * b2 * g1^2 + b3 * g0 * g1^2 - b4 * g1^3 + b1^3 * 
           g2 - b3 * g0^2 * g2 + b2^2 * g1 * g2 + b1 * b3 * g1 * 
           g2 + 2 * b4 * g0 * g1 * g2 - b2 * b3 * g2^2 + b1 * b4 * 
           g2^2 - b1^2 * b2 * g3 - b2^2 * g0 * g3 - b1 * b3 * g0 * 
           g3 - b4 * g0^2 * g3 + b2 * b3 * g1 * g3 - b1 * b4 * g1 * 
           g3)/den
  q2 <- (b1^2 * g0^2 + b2 * g0^3 + b1^3 * g1 - b3 * g0^2 * 
           g1 + b1 * b3 * g1^2 + b4 * g0 * g1^2 - b1^2 * b2 * g2 + 
           b2^2 * g0 * g2 - 3 * b1 * b3 * g0 * g2 - b4 * g0^2 * 
           g2 - b2 * b3 * g1 * g2 + b1 * b4 * g1 * g2 + b3^2 * g2^2 - 
           b2 * b4 * g2^2 + b1 * b2^2 * g3 - b1^2 * b3 * g3 + b2 * 
           b3 * g0 * g3 - b1 * b4 * g0 * g3 - b3^2 * g1 * g3 + b2 * 
           b4 * g1 * g3)/den
  q3 <- (b1^3 * g0 + 2 * b1 * b2 * g0^2 + b3 * g0^3 - b1^2 * 
           b2 * g1 - 2 * b2^2 * g0 * g1 - b4 * g0^2 * g1 + b2 * 
           b3 * g1^2 - b1 * b4 * g1^2 + b1 * b2^2 * g2 - b1^2 * 
           b3 * g2 + b2 * b3 * g0 * g2 - b1 * b4 * g0 * g2 - b3^2 * 
           g1 * g2 + b2 * b4 * g1 * g2 - b2^3 * g3 + 2 * b1 * b2 * 
           b3 * g3 - b1^2 * b4 * g3 + b3^2 * g0 * g3 - b2 * b4 * 
           g0 * g3)/den
  q4 <- (b1^4 + 3 * b1^2 * b2 * g0 + b2^2 * g0^2 + 2 * b1 * 
           b3 * g0^2 + b4 * g0^3 - 2 * b1 * b2^2 * g1 + 2 * b1^2 * 
           b3 * g1 - 2 * b2 * b3 * g0 * g1 + 2 * b1 * b4 * g0 * 
           g1 + b3^2 * g1^2 - b2 * b4 * g1^2 + b2^3 * g2 - 2 * b1 * 
           b2 * b3 * g2 + b1^2 * b4 * g2 - b3^2 * g0 * g2 + b2 * 
           b4 * g0 * g2)/den
  p4 <- g0 * q4
  y <- x^al
  h.pade <- (p1 * y + p2 * y^2 + p3 * y^3 + p4 * y^4)/
                (1 + q1 * y + q2 * y^2 + q3 * y^3 + q4 * y^4)
  #    res <- 1/2 * (h.pade - rm) * (h.pade - rp)
  return(h.pade)
}

########################################################################

Power <- function(a,b){a^b}

h.Pade55 <- function (params)function(a, x) {
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  b1 <- -a*(a+1i)/2 * 1/gamma(1+al)
  b2 <- -b1 * lamTilde * nu * gamma(1+al)/gamma(1+2*al)               
  b3 <- (-b2 * lamTilde * nu  + nu^2 * b1^2/2) * gamma(1+2*al)/gamma(1+3*al)
  b4 <- (-b3 * lamTilde *nu + nu^2 * b1*b2)* gamma(1+3*al)/gamma(1+4*al)
  b5 <- (-b4 * lamTilde * nu + nu^2*(1/2*b2*b2 +b1*b3))* gamma(1+4*al)/gamma(1+5*al)
  
  g0 <- rm/nu
  g1 <- -1/(aa*nu)*1/gamma(1-al)*g0
  g2 <- -1/(aa*nu)*(gamma(1-al)/gamma(1-2*al)*g1 -1/2*nu^2*g1*g1)
  g3 <- -1/(aa*nu)*(gamma(1-2*al)/gamma(1-3*al)*g2 -nu^2 * g1*g2) 
  g4 <- -1/(aa*nu)*(gamma(1-3*al)/gamma(1-4*al)*g3 - nu^2*(1/2*g2*g2 + g1*g3)) 

  den <- -Power(g0,5) - 4*b1*Power(g0,3)*g1 - 3*Power(b1,2)*g0*Power(g1,2) + 3*b2*Power(g0,2)*Power(g1,2) + 2*b1*b2*Power(g1,3) - 
    2*b3*g0*Power(g1,3) + b4*Power(g1,4) - 3*Power(b1,2)*Power(g0,2)*g2 - 3*b2*Power(g0,3)*g2 - 2*Power(b1,3)*g1*g2 + 
    2*b1*b2*g0*g1*g2 + 4*b3*Power(g0,2)*g1*g2 - Power(b2,2)*Power(g1,2)*g2 - 2*b1*b3*Power(g1,2)*g2 - 3*b4*g0*Power(g1,2)*g2 + 
    Power(b1,2)*b2*Power(g2,2) - 2*Power(b2,2)*g0*Power(g2,2) + 4*b1*b3*g0*Power(g2,2) + b4*Power(g0,2)*Power(g2,2) + 
    2*b2*b3*g1*Power(g2,2) - 2*b1*b4*g1*Power(g2,2) - Power(b3,2)*Power(g2,3) + b2*b4*Power(g2,3) - 2*Power(b1,3)*g0*g3 - 
    4*b1*b2*Power(g0,2)*g3 - 2*b3*Power(g0,3)*g3 + 2*Power(b1,2)*b2*g1*g3 + 4*Power(b2,2)*g0*g1*g3 + 2*b4*Power(g0,2)*g1*g3 - 
    2*b2*b3*Power(g1,2)*g3 + 2*b1*b4*Power(g1,2)*g3 - 2*b1*Power(b2,2)*g2*g3 + 2*Power(b1,2)*b3*g2*g3 - 2*b2*b3*g0*g2*g3 + 
    2*b1*b4*g0*g2*g3 + 2*Power(b3,2)*g1*g2*g3 - 2*b2*b4*g1*g2*g3 + Power(b2,3)*Power(g3,2) - 2*b1*b2*b3*Power(g3,2) + 
    Power(b1,2)*b4*Power(g3,2) - Power(b3,2)*g0*Power(g3,2) + b2*b4*g0*Power(g3,2) - Power(b1,4)*g4 - 3*Power(b1,2)*b2*g0*g4 - 
    Power(b2,2)*Power(g0,2)*g4 - 2*b1*b3*Power(g0,2)*g4 - b4*Power(g0,3)*g4 + 2*b1*Power(b2,2)*g1*g4 - 2*Power(b1,2)*b3*g1*g4 + 
    2*b2*b3*g0*g1*g4 - 2*b1*b4*g0*g1*g4 - Power(b3,2)*Power(g1,2)*g4 + b2*b4*Power(g1,2)*g4 - Power(b2,3)*g2*g4 + 
    2*b1*b2*b3*g2*g4 - Power(b1,2)*b4*g2*g4 + Power(b3,2)*g0*g2*g4 - b2*b4*g0*g2*g4
  
  # Numerator
  
  p1 <- b1
  
  p2 <- (-(Power(b1,2)*Power(g0,4)) - b2*Power(g0,5) - 3*Power(b1,3)*Power(g0,2)*g1 - 3*b1*b2*Power(g0,3)*g1 - 
           Power(b1,4)*Power(g1,2) + Power(b1,2)*b2*g0*Power(g1,2) + 3*Power(b2,2)*Power(g0,2)*Power(g1,2) - 
           b1*b3*Power(g0,2)*Power(g1,2) + b1*Power(b2,2)*Power(g1,3) - 2*Power(b1,2)*b3*Power(g1,3) - 2*b2*b3*g0*Power(g1,3) + 
           b1*b4*g0*Power(g1,3) + b2*b4*Power(g1,4) - b1*b5*Power(g1,4) - 2*Power(b1,4)*g0*g2 - 4*Power(b1,2)*b2*Power(g0,2)*g2 - 
           3*Power(b2,2)*Power(g0,3)*g2 + b1*b3*Power(g0,3)*g2 + 2*Power(b1,3)*b2*g1*g2 + 2*b1*Power(b2,2)*g0*g1*g2 + 
           2*Power(b1,2)*b3*g0*g1*g2 + 4*b2*b3*Power(g0,2)*g1*g2 - 2*b1*b4*Power(g0,2)*g1*g2 - Power(b2,3)*Power(g1,2)*g2 + 
           Power(b1,2)*b4*Power(g1,2)*g2 - 3*b2*b4*g0*Power(g1,2)*g2 + 3*b1*b5*g0*Power(g1,2)*g2 - 
           Power(b1,2)*Power(b2,2)*Power(g2,2) + Power(b1,3)*b3*Power(g2,2) - 2*Power(b2,3)*g0*Power(g2,2) + 
           4*b1*b2*b3*g0*Power(g2,2) - 2*Power(b1,2)*b4*g0*Power(g2,2) + b2*b4*Power(g0,2)*Power(g2,2) - 
           b1*b5*Power(g0,2)*Power(g2,2) + 2*Power(b2,2)*b3*g1*Power(g2,2) - b1*Power(b3,2)*g1*Power(g2,2) - 
           3*b1*b2*b4*g1*Power(g2,2) + 2*Power(b1,2)*b5*g1*Power(g2,2) - b2*Power(b3,2)*Power(g2,3) + Power(b2,2)*b4*Power(g2,3) + 
           b1*b3*b4*Power(g2,3) - b1*b2*b5*Power(g2,3) - Power(b1,5)*g3 - 3*Power(b1,3)*b2*g0*g3 - 3*b1*Power(b2,2)*Power(g0,2)*g3 - 
           2*b2*b3*Power(g0,3)*g3 + b1*b4*Power(g0,3)*g3 + 2*Power(b1,2)*Power(b2,2)*g1*g3 - 2*Power(b1,3)*b3*g1*g3 + 
           4*Power(b2,3)*g0*g1*g3 - 4*b1*b2*b3*g0*g1*g3 + 2*b2*b4*Power(g0,2)*g1*g3 - 2*b1*b5*Power(g0,2)*g1*g3 - 
           2*Power(b2,2)*b3*Power(g1,2)*g3 + b1*Power(b3,2)*Power(g1,2)*g3 + 3*b1*b2*b4*Power(g1,2)*g3 - 
           2*Power(b1,2)*b5*Power(g1,2)*g3 - b1*Power(b2,3)*g2*g3 + 2*Power(b1,2)*b2*b3*g2*g3 - Power(b1,3)*b4*g2*g3 - 
           2*Power(b2,2)*b3*g0*g2*g3 + b1*Power(b3,2)*g0*g2*g3 + 3*b1*b2*b4*g0*g2*g3 - 2*Power(b1,2)*b5*g0*g2*g3 + 
           2*b2*Power(b3,2)*g1*g2*g3 - 2*Power(b2,2)*b4*g1*g2*g3 - 2*b1*b3*b4*g1*g2*g3 + 2*b1*b2*b5*g1*g2*g3 + 
           Power(b2,4)*Power(g3,2) - 3*b1*Power(b2,2)*b3*Power(g3,2) + Power(b1,2)*Power(b3,2)*Power(g3,2) + 
           2*Power(b1,2)*b2*b4*Power(g3,2) - Power(b1,3)*b5*Power(g3,2) - b2*Power(b3,2)*g0*Power(g3,2) + 
           Power(b2,2)*b4*g0*Power(g3,2) + b1*b3*b4*g0*Power(g3,2) - b1*b2*b5*g0*Power(g3,2) - Power(b1,2)*Power(b2,2)*g0*g4 + 
           Power(b1,3)*b3*g0*g4 - Power(b2,3)*Power(g0,2)*g4 + Power(b1,2)*b4*Power(g0,2)*g4 - b2*b4*Power(g0,3)*g4 + 
           b1*b5*Power(g0,3)*g4 + b1*Power(b2,3)*g1*g4 - 2*Power(b1,2)*b2*b3*g1*g4 + Power(b1,3)*b4*g1*g4 + 
           2*Power(b2,2)*b3*g0*g1*g4 - b1*Power(b3,2)*g0*g1*g4 - 3*b1*b2*b4*g0*g1*g4 + 2*Power(b1,2)*b5*g0*g1*g4 - 
           b2*Power(b3,2)*Power(g1,2)*g4 + Power(b2,2)*b4*Power(g1,2)*g4 + b1*b3*b4*Power(g1,2)*g4 - b1*b2*b5*Power(g1,2)*g4 - 
           Power(b2,4)*g2*g4 + 3*b1*Power(b2,2)*b3*g2*g4 - Power(b1,2)*Power(b3,2)*g2*g4 - 2*Power(b1,2)*b2*b4*g2*g4 + 
           Power(b1,3)*b5*g2*g4 + b2*Power(b3,2)*g0*g2*g4 - Power(b2,2)*b4*g0*g2*g4 - b1*b3*b4*g0*g2*g4 + b1*b2*b5*g0*g2*g4)/den
  
  p3 <- (-(Power(b1,3)*Power(g0,3)) - 2*b1*b2*Power(g0,4) - b3*Power(g0,5) - 2*Power(b1,4)*g0*g1 - 
           4*Power(b1,2)*b2*Power(g0,2)*g1 + Power(b2,2)*Power(g0,3)*g1 - 3*b1*b3*Power(g0,3)*g1 + Power(b1,3)*b2*Power(g1,2) + 
           5*b1*Power(b2,2)*g0*Power(g1,2) - 3*Power(b1,2)*b3*g0*Power(g1,2) + 2*b2*b3*Power(g0,2)*Power(g1,2) - 
           b1*b4*Power(g0,2)*Power(g1,2) - Power(b2,3)*Power(g1,3) + Power(b1,2)*b4*Power(g1,3) - 2*Power(b3,2)*g0*Power(g1,3) + 
           b2*b4*g0*Power(g1,3) + b1*b5*g0*Power(g1,3) + b3*b4*Power(g1,4) - b2*b5*Power(g1,4) - Power(b1,5)*g2 - 
           3*Power(b1,3)*b2*g0*g2 - 3*b1*Power(b2,2)*Power(g0,2)*g2 - 2*b2*b3*Power(g0,3)*g2 + b1*b4*Power(g0,3)*g2 + 
           2*Power(b1,2)*Power(b2,2)*g1*g2 - 2*Power(b1,3)*b3*g1*g2 + 4*b1*b2*b3*g0*g1*g2 - 4*Power(b1,2)*b4*g0*g1*g2 + 
           4*Power(b3,2)*Power(g0,2)*g1*g2 - 2*b2*b4*Power(g0,2)*g1*g2 - 2*b1*b5*Power(g0,2)*g1*g2 + Power(b2,2)*b3*Power(g1,2)*g2 - 
           2*b1*Power(b3,2)*Power(g1,2)*g2 + Power(b1,2)*b5*Power(g1,2)*g2 - 3*b3*b4*g0*Power(g1,2)*g2 + 3*b2*b5*g0*Power(g1,2)*g2 - 
           2*b1*Power(b2,3)*Power(g2,2) + 4*Power(b1,2)*b2*b3*Power(g2,2) - 2*Power(b1,3)*b4*Power(g2,2) - 
           2*Power(b2,2)*b3*g0*Power(g2,2) + 3*b1*Power(b3,2)*g0*Power(g2,2) + b1*b2*b4*g0*Power(g2,2) - 
           2*Power(b1,2)*b5*g0*Power(g2,2) + b3*b4*Power(g0,2)*Power(g2,2) - b2*b5*Power(g0,2)*Power(g2,2) + 
           b2*Power(b3,2)*g1*Power(g2,2) - Power(b2,2)*b4*g1*Power(g2,2) - b1*b3*b4*g1*Power(g2,2) + b1*b2*b5*g1*Power(g2,2) - 
           Power(b3,3)*Power(g2,3) + 2*b2*b3*b4*Power(g2,3) - b1*Power(b4,2)*Power(g2,3) - Power(b2,2)*b5*Power(g2,3) + 
           b1*b3*b5*Power(g2,3) - Power(b1,2)*Power(b2,2)*g0*g3 + Power(b1,3)*b3*g0*g3 + Power(b2,3)*Power(g0,2)*g3 - 
           4*b1*b2*b3*Power(g0,2)*g3 + 3*Power(b1,2)*b4*Power(g0,2)*g3 - 2*Power(b3,2)*Power(g0,3)*g3 + b2*b4*Power(g0,3)*g3 + 
           b1*b5*Power(g0,3)*g3 + b1*Power(b2,3)*g1*g3 - 2*Power(b1,2)*b2*b3*g1*g3 + Power(b1,3)*b4*g1*g3 + b1*Power(b3,2)*g0*g1*g3 - 
           b1*b2*b4*g0*g1*g3 + 2*b3*b4*Power(g0,2)*g1*g3 - 2*b2*b5*Power(g0,2)*g1*g3 - b2*Power(b3,2)*Power(g1,2)*g3 + 
           Power(b2,2)*b4*Power(g1,2)*g3 + b1*b3*b4*Power(g1,2)*g3 - b1*b2*b5*Power(g1,2)*g3 + Power(b2,4)*g2*g3 - 
           3*b1*Power(b2,2)*b3*g2*g3 + Power(b1,2)*Power(b3,2)*g2*g3 + 2*Power(b1,2)*b2*b4*g2*g3 - Power(b1,3)*b5*g2*g3 - 
           b2*Power(b3,2)*g0*g2*g3 + Power(b2,2)*b4*g0*g2*g3 + b1*b3*b4*g0*g2*g3 - b1*b2*b5*g0*g2*g3 + 2*Power(b3,3)*g1*g2*g3 - 
           4*b2*b3*b4*g1*g2*g3 + 2*b1*Power(b4,2)*g1*g2*g3 + 2*Power(b2,2)*b5*g1*g2*g3 - 2*b1*b3*b5*g1*g2*g3 - 
           Power(b3,3)*g0*Power(g3,2) + 2*b2*b3*b4*g0*Power(g3,2) - b1*Power(b4,2)*g0*Power(g3,2) - Power(b2,2)*b5*g0*Power(g3,2) + 
           b1*b3*b5*g0*Power(g3,2) + b1*Power(b2,3)*g0*g4 - 2*Power(b1,2)*b2*b3*g0*g4 + Power(b1,3)*b4*g0*g4 + 
           Power(b2,2)*b3*Power(g0,2)*g4 - 2*b1*Power(b3,2)*Power(g0,2)*g4 + Power(b1,2)*b5*Power(g0,2)*g4 - b3*b4*Power(g0,3)*g4 + 
           b2*b5*Power(g0,3)*g4 - Power(b2,4)*g1*g4 + 3*b1*Power(b2,2)*b3*g1*g4 - Power(b1,2)*Power(b3,2)*g1*g4 - 
           2*Power(b1,2)*b2*b4*g1*g4 + Power(b1,3)*b5*g1*g4 + b2*Power(b3,2)*g0*g1*g4 - Power(b2,2)*b4*g0*g1*g4 - b1*b3*b4*g0*g1*g4 + 
           b1*b2*b5*g0*g1*g4 - Power(b3,3)*Power(g1,2)*g4 + 2*b2*b3*b4*Power(g1,2)*g4 - b1*Power(b4,2)*Power(g1,2)*g4 - 
           Power(b2,2)*b5*Power(g1,2)*g4 + b1*b3*b5*Power(g1,2)*g4 + Power(b3,3)*g0*g2*g4 - 2*b2*b3*b4*g0*g2*g4 + 
           b1*Power(b4,2)*g0*g2*g4 + Power(b2,2)*b5*g0*g2*g4 - b1*b3*b5*g0*g2*g4)/den
  
  p4 <- (-(Power(b1,4)*Power(g0,2)) - 3*Power(b1,2)*b2*Power(g0,3) - Power(b2,2)*Power(g0,4) - 2*b1*b3*Power(g0,4) - 
           b4*Power(g0,5) - Power(b1,5)*g1 - 3*Power(b1,3)*b2*g0*g1 + b1*Power(b2,2)*Power(g0,2)*g1 - 
           4*Power(b1,2)*b3*Power(g0,2)*g1 + 2*b2*b3*Power(g0,3)*g1 - 3*b1*b4*Power(g0,3)*g1 + 
           3*Power(b1,2)*Power(b2,2)*Power(g1,2) - 3*Power(b1,3)*b3*Power(g1,2) + Power(b2,3)*g0*Power(g1,2) + 
           2*b1*b2*b3*g0*Power(g1,2) - 3*Power(b1,2)*b4*g0*Power(g1,2) - Power(b3,2)*Power(g0,2)*Power(g1,2) + 
           2*b2*b4*Power(g0,2)*Power(g1,2) - b1*b5*Power(g0,2)*Power(g1,2) - Power(b2,2)*b3*Power(g1,3) - 
           2*b1*Power(b3,2)*Power(g1,3) + 4*b1*b2*b4*Power(g1,3) - Power(b1,2)*b5*Power(g1,3) - b3*b4*g0*Power(g1,3) + 
           b2*b5*g0*Power(g1,3) + Power(b4,2)*Power(g1,4) - b3*b5*Power(g1,4) - Power(b1,2)*Power(b2,2)*g0*g2 + 
           Power(b1,3)*b3*g0*g2 - 2*Power(b2,3)*Power(g0,2)*g2 + 2*b1*b2*b3*Power(g0,2)*g2 + Power(b3,2)*Power(g0,3)*g2 - 
           2*b2*b4*Power(g0,3)*g2 + b1*b5*Power(g0,3)*g2 - 2*b1*Power(b2,3)*g1*g2 + 4*Power(b1,2)*b2*b3*g1*g2 - 
           2*Power(b1,3)*b4*g1*g2 + 4*b1*Power(b3,2)*g0*g1*g2 - 4*b1*b2*b4*g0*g1*g2 + 2*b3*b4*Power(g0,2)*g1*g2 - 
           2*b2*b5*Power(g0,2)*g1*g2 + 2*b2*Power(b3,2)*Power(g1,2)*g2 - 2*Power(b2,2)*b4*Power(g1,2)*g2 - 
           2*b1*b3*b4*Power(g1,2)*g2 + 2*b1*b2*b5*Power(g1,2)*g2 - 3*Power(b4,2)*g0*Power(g1,2)*g2 + 3*b3*b5*g0*Power(g1,2)*g2 - 
           b2*Power(b3,2)*g0*Power(g2,2) + Power(b2,2)*b4*g0*Power(g2,2) + b1*b3*b4*g0*Power(g2,2) - b1*b2*b5*g0*Power(g2,2) + 
           Power(b4,2)*Power(g0,2)*Power(g2,2) - b3*b5*Power(g0,2)*Power(g2,2) - Power(b3,3)*g1*Power(g2,2) + 
           2*b2*b3*b4*g1*Power(g2,2) - b1*Power(b4,2)*g1*Power(g2,2) - Power(b2,2)*b5*g1*Power(g2,2) + b1*b3*b5*g1*Power(g2,2) + 
           b1*Power(b2,3)*g0*g3 - 2*Power(b1,2)*b2*b3*g0*g3 + Power(b1,3)*b4*g0*g3 + Power(b2,2)*b3*Power(g0,2)*g3 - 
           2*b1*Power(b3,2)*Power(g0,2)*g3 + Power(b1,2)*b5*Power(g0,2)*g3 - b3*b4*Power(g0,3)*g3 + b2*b5*Power(g0,3)*g3 + 
           Power(b2,4)*g1*g3 - 3*b1*Power(b2,2)*b3*g1*g3 + Power(b1,2)*Power(b3,2)*g1*g3 + 2*Power(b1,2)*b2*b4*g1*g3 - 
           Power(b1,3)*b5*g1*g3 - 3*b2*Power(b3,2)*g0*g1*g3 + 3*Power(b2,2)*b4*g0*g1*g3 + 3*b1*b3*b4*g0*g1*g3 - 3*b1*b2*b5*g0*g1*g3 + 
           2*Power(b4,2)*Power(g0,2)*g1*g3 - 2*b3*b5*Power(g0,2)*g1*g3 + Power(b3,3)*Power(g1,2)*g3 - 2*b2*b3*b4*Power(g1,2)*g3 + 
           b1*Power(b4,2)*Power(g1,2)*g3 + Power(b2,2)*b5*Power(g1,2)*g3 - b1*b3*b5*Power(g1,2)*g3 + Power(b3,3)*g0*g2*g3 - 
           2*b2*b3*b4*g0*g2*g3 + b1*Power(b4,2)*g0*g2*g3 + Power(b2,2)*b5*g0*g2*g3 - b1*b3*b5*g0*g2*g3 - Power(b2,4)*g0*g4 + 
           3*b1*Power(b2,2)*b3*g0*g4 - Power(b1,2)*Power(b3,2)*g0*g4 - 2*Power(b1,2)*b2*b4*g0*g4 + Power(b1,3)*b5*g0*g4 + 
           2*b2*Power(b3,2)*Power(g0,2)*g4 - 2*Power(b2,2)*b4*Power(g0,2)*g4 - 2*b1*b3*b4*Power(g0,2)*g4 + 
           2*b1*b2*b5*Power(g0,2)*g4 - Power(b4,2)*Power(g0,3)*g4 + b3*b5*Power(g0,3)*g4 - Power(b3,3)*g0*g1*g4 + 
           2*b2*b3*b4*g0*g1*g4 - b1*Power(b4,2)*g0*g1*g4 - Power(b2,2)*b5*g0*g1*g4 + b1*b3*b5*g0*g1*g4)/den
  
  p5 <- -((g0*(Power(b1,5) + 4*Power(b1,3)*b2*g0 + 3*b1*Power(b2,2)*Power(g0,2) + 3*Power(b1,2)*b3*Power(g0,2) + 
                 2*b2*b3*Power(g0,3) + 2*b1*b4*Power(g0,3) + b5*Power(g0,4) - 3*Power(b1,2)*Power(b2,2)*g1 + 3*Power(b1,3)*b3*g1 - 
                 2*Power(b2,3)*g0*g1 - 2*b1*b2*b3*g0*g1 + 4*Power(b1,2)*b4*g0*g1 - Power(b3,2)*Power(g0,2)*g1 - 
                 2*b2*b4*Power(g0,2)*g1 + 3*b1*b5*Power(g0,2)*g1 + Power(b2,2)*b3*Power(g1,2) + 2*b1*Power(b3,2)*Power(g1,2) - 
                 4*b1*b2*b4*Power(g1,2) + Power(b1,2)*b5*Power(g1,2) + 2*b3*b4*g0*Power(g1,2) - 2*b2*b5*g0*Power(g1,2) - 
                 Power(b4,2)*Power(g1,3) + b3*b5*Power(g1,3) + 2*b1*Power(b2,3)*g2 - 4*Power(b1,2)*b2*b3*g2 + 2*Power(b1,3)*b4*g2 + 
                 2*Power(b2,2)*b3*g0*g2 - 4*b1*Power(b3,2)*g0*g2 + 2*Power(b1,2)*b5*g0*g2 - 2*b3*b4*Power(g0,2)*g2 + 
                 2*b2*b5*Power(g0,2)*g2 - 2*b2*Power(b3,2)*g1*g2 + 2*Power(b2,2)*b4*g1*g2 + 2*b1*b3*b4*g1*g2 - 2*b1*b2*b5*g1*g2 + 
                 2*Power(b4,2)*g0*g1*g2 - 2*b3*b5*g0*g1*g2 + Power(b3,3)*Power(g2,2) - 2*b2*b3*b4*Power(g2,2) + 
                 b1*Power(b4,2)*Power(g2,2) + Power(b2,2)*b5*Power(g2,2) - b1*b3*b5*Power(g2,2) - Power(b2,4)*g3 + 
                 3*b1*Power(b2,2)*b3*g3 - Power(b1,2)*Power(b3,2)*g3 - 2*Power(b1,2)*b2*b4*g3 + Power(b1,3)*b5*g3 + 
                 2*b2*Power(b3,2)*g0*g3 - 2*Power(b2,2)*b4*g0*g3 - 2*b1*b3*b4*g0*g3 + 2*b1*b2*b5*g0*g3 - Power(b4,2)*Power(g0,2)*g3 + 
                 b3*b5*Power(g0,2)*g3 - Power(b3,3)*g1*g3 + 2*b2*b3*b4*g1*g3 - b1*Power(b4,2)*g1*g3 - Power(b2,2)*b5*g1*g3 + 
                 b1*b3*b5*g1*g3))/den)
  
  # Denominator
  
  q1 <- (-(b1*Power(g0,4)) - 3*Power(b1,2)*Power(g0,2)*g1 + b2*Power(g0,3)*g1 - Power(b1,3)*Power(g1,2) + 
           4*b1*b2*g0*Power(g1,2) - b3*Power(g0,2)*Power(g1,2) - Power(b2,2)*Power(g1,3) - 2*b1*b3*Power(g1,3) + b4*g0*Power(g1,3) - 
           b5*Power(g1,4) - 2*Power(b1,3)*g0*g2 - b1*b2*Power(g0,2)*g2 + b3*Power(g0,3)*g2 + 4*Power(b1,2)*b2*g1*g2 + 
           2*b1*b3*g0*g1*g2 - 2*b4*Power(g0,2)*g1*g2 + 2*b2*b3*Power(g1,2)*g2 + b1*b4*Power(g1,2)*g2 + 3*b5*g0*Power(g1,2)*g2 - 
           2*b1*Power(b2,2)*Power(g2,2) + Power(b1,2)*b3*Power(g2,2) - 2*b1*b4*g0*Power(g2,2) - b5*Power(g0,2)*Power(g2,2) - 
           Power(b3,2)*g1*Power(g2,2) - b2*b4*g1*Power(g2,2) + 2*b1*b5*g1*Power(g2,2) + b3*b4*Power(g2,3) - b2*b5*Power(g2,3) - 
           Power(b1,4)*g3 - Power(b1,2)*b2*g0*g3 + Power(b2,2)*Power(g0,2)*g3 + b4*Power(g0,3)*g3 - 2*Power(b1,2)*b3*g1*g3 - 
           4*b2*b3*g0*g1*g3 - 2*b5*Power(g0,2)*g1*g3 + Power(b3,2)*Power(g1,2)*g3 + b2*b4*Power(g1,2)*g3 - 2*b1*b5*Power(g1,2)*g3 + 
           Power(b2,3)*g2*g3 - Power(b1,2)*b4*g2*g3 + Power(b3,2)*g0*g2*g3 + b2*b4*g0*g2*g3 - 2*b1*b5*g0*g2*g3 - 2*b3*b4*g1*g2*g3 + 
           2*b2*b5*g1*g2*g3 - Power(b2,2)*b3*Power(g3,2) + b1*Power(b3,2)*Power(g3,2) + b1*b2*b4*Power(g3,2) - 
           Power(b1,2)*b5*Power(g3,2) + b3*b4*g0*Power(g3,2) - b2*b5*g0*Power(g3,2) + Power(b1,3)*b2*g4 + 2*b1*Power(b2,2)*g0*g4 + 
           Power(b1,2)*b3*g0*g4 + 2*b2*b3*Power(g0,2)*g4 + b1*b4*Power(g0,2)*g4 + b5*Power(g0,3)*g4 - Power(b2,3)*g1*g4 + 
           Power(b1,2)*b4*g1*g4 - Power(b3,2)*g0*g1*g4 - b2*b4*g0*g1*g4 + 2*b1*b5*g0*g1*g4 + b3*b4*Power(g1,2)*g4 - 
           b2*b5*Power(g1,2)*g4 + Power(b2,2)*b3*g2*g4 - b1*Power(b3,2)*g2*g4 - b1*b2*b4*g2*g4 + Power(b1,2)*b5*g2*g4 - 
           b3*b4*g0*g2*g4 + b2*b5*g0*g2*g4)/den
  
  q2 <- (-(Power(b1,2)*Power(g0,3)) - b2*Power(g0,4) - 2*Power(b1,3)*g0*g1 - b1*b2*Power(g0,2)*g1 + b3*Power(g0,3)*g1 + 
           2*Power(b1,2)*b2*Power(g1,2) + Power(b2,2)*g0*Power(g1,2) - b4*Power(g0,2)*Power(g1,2) + b1*b4*Power(g1,3) + 
           b5*g0*Power(g1,3) - Power(b1,4)*g2 - Power(b1,2)*b2*g0*g2 - 2*Power(b2,2)*Power(g0,2)*g2 + 3*b1*b3*Power(g0,2)*g2 + 
           b4*Power(g0,3)*g2 - 2*b1*Power(b2,2)*g1*g2 - 4*b1*b4*g0*g1*g2 - 2*b5*Power(g0,2)*g1*g2 - b2*b4*Power(g1,2)*g2 + 
           b1*b5*Power(g1,2)*g2 + 2*b1*b2*b3*Power(g2,2) - 2*Power(b1,2)*b4*Power(g2,2) - Power(b3,2)*g0*Power(g2,2) + 
           3*b2*b4*g0*Power(g2,2) - 2*b1*b5*g0*Power(g2,2) + b3*b4*g1*Power(g2,2) - b2*b5*g1*Power(g2,2) - Power(b4,2)*Power(g2,3) + 
           b3*b5*Power(g2,3) + Power(b1,3)*b2*g3 + 3*Power(b1,2)*b3*g0*g3 + 3*b1*b4*Power(g0,2)*g3 + b5*Power(g0,3)*g3 + 
           Power(b2,3)*g1*g3 - 2*b1*b2*b3*g1*g3 + Power(b1,2)*b4*g1*g3 + Power(b3,2)*g0*g1*g3 - b2*b4*g0*g1*g3 - 
           b3*b4*Power(g1,2)*g3 + b2*b5*Power(g1,2)*g3 - Power(b2,2)*b3*g2*g3 - b1*Power(b3,2)*g2*g3 + 3*b1*b2*b4*g2*g3 - 
           Power(b1,2)*b5*g2*g3 - b3*b4*g0*g2*g3 + b2*b5*g0*g2*g3 + 2*Power(b4,2)*g1*g2*g3 - 2*b3*b5*g1*g2*g3 + 
           b2*Power(b3,2)*Power(g3,2) - Power(b2,2)*b4*Power(g3,2) - b1*b3*b4*Power(g3,2) + b1*b2*b5*Power(g3,2) - 
           Power(b4,2)*g0*Power(g3,2) + b3*b5*g0*Power(g3,2) - Power(b1,2)*Power(b2,2)*g4 + Power(b1,3)*b3*g4 - Power(b2,3)*g0*g4 + 
           Power(b1,2)*b4*g0*g4 - b2*b4*Power(g0,2)*g4 + b1*b5*Power(g0,2)*g4 + Power(b2,2)*b3*g1*g4 + b1*Power(b3,2)*g1*g4 - 
           3*b1*b2*b4*g1*g4 + Power(b1,2)*b5*g1*g4 + b3*b4*g0*g1*g4 - b2*b5*g0*g1*g4 - Power(b4,2)*Power(g1,2)*g4 + 
           b3*b5*Power(g1,2)*g4 - b2*Power(b3,2)*g2*g4 + Power(b2,2)*b4*g2*g4 + b1*b3*b4*g2*g4 - b1*b2*b5*g2*g4 + 
           Power(b4,2)*g0*g2*g4 - b3*b5*g0*g2*g4)/den
  
  q3 <- (-(Power(b1,3)*Power(g0,2)) - 2*b1*b2*Power(g0,3) - b3*Power(g0,4) - Power(b1,4)*g1 - Power(b1,2)*b2*g0*g1 + 
           2*Power(b2,2)*Power(g0,2)*g1 - b1*b3*Power(g0,2)*g1 + b4*Power(g0,3)*g1 + b1*Power(b2,2)*Power(g1,2) - 
           2*Power(b1,2)*b3*Power(g1,2) - 2*b2*b3*g0*Power(g1,2) - b5*Power(g0,2)*Power(g1,2) + b2*b4*Power(g1,3) - 
           b1*b5*Power(g1,3) + Power(b1,3)*b2*g2 + 3*Power(b1,2)*b3*g0*g2 + 3*b1*b4*Power(g0,2)*g2 + b5*Power(g0,3)*g2 + 
           2*Power(b3,2)*g0*g1*g2 - 2*b2*b4*g0*g1*g2 - b3*b4*Power(g1,2)*g2 + b2*b5*Power(g1,2)*g2 - b1*Power(b3,2)*Power(g2,2) + 
           b1*b2*b4*Power(g2,2) - b3*b4*g0*Power(g2,2) + b2*b5*g0*Power(g2,2) + Power(b4,2)*g1*Power(g2,2) - b3*b5*g1*Power(g2,2) - 
           Power(b1,2)*Power(b2,2)*g3 + Power(b1,3)*b3*g3 + Power(b2,3)*g0*g3 - 4*b1*b2*b3*g0*g3 + 3*Power(b1,2)*b4*g0*g3 - 
           2*Power(b3,2)*Power(g0,2)*g3 + b2*b4*Power(g0,2)*g3 + b1*b5*Power(g0,2)*g3 - Power(b2,2)*b3*g1*g3 + 
           3*b1*Power(b3,2)*g1*g3 - b1*b2*b4*g1*g3 - Power(b1,2)*b5*g1*g3 + 3*b3*b4*g0*g1*g3 - 3*b2*b5*g0*g1*g3 - 
           Power(b4,2)*Power(g1,2)*g3 + b3*b5*Power(g1,2)*g3 + b2*Power(b3,2)*g2*g3 - Power(b2,2)*b4*g2*g3 - b1*b3*b4*g2*g3 + 
           b1*b2*b5*g2*g3 - Power(b4,2)*g0*g2*g3 + b3*b5*g0*g2*g3 - Power(b3,3)*Power(g3,2) + 2*b2*b3*b4*Power(g3,2) - 
           b1*Power(b4,2)*Power(g3,2) - Power(b2,2)*b5*Power(g3,2) + b1*b3*b5*Power(g3,2) + b1*Power(b2,3)*g4 - 
           2*Power(b1,2)*b2*b3*g4 + Power(b1,3)*b4*g4 + Power(b2,2)*b3*g0*g4 - 2*b1*Power(b3,2)*g0*g4 + Power(b1,2)*b5*g0*g4 - 
           b3*b4*Power(g0,2)*g4 + b2*b5*Power(g0,2)*g4 - b2*Power(b3,2)*g1*g4 + Power(b2,2)*b4*g1*g4 + b1*b3*b4*g1*g4 - 
           b1*b2*b5*g1*g4 + Power(b4,2)*g0*g1*g4 - b3*b5*g0*g1*g4 + Power(b3,3)*g2*g4 - 2*b2*b3*b4*g2*g4 + b1*Power(b4,2)*g2*g4 + 
           Power(b2,2)*b5*g2*g4 - b1*b3*b5*g2*g4)/den
  
  q4 <- (-(Power(b1,4)*g0) - 3*Power(b1,2)*b2*Power(g0,2) - Power(b2,2)*Power(g0,3) - 2*b1*b3*Power(g0,3) - b4*Power(g0,4) + 
           Power(b1,3)*b2*g1 + 4*b1*Power(b2,2)*g0*g1 - Power(b1,2)*b3*g0*g1 + 4*b2*b3*Power(g0,2)*g1 - b1*b4*Power(g0,2)*g1 + 
           b5*Power(g0,3)*g1 - Power(b2,3)*Power(g1,2) + Power(b1,2)*b4*Power(g1,2) - 2*Power(b3,2)*g0*Power(g1,2) + 
           2*b1*b5*g0*Power(g1,2) + b3*b4*Power(g1,3) - b2*b5*Power(g1,3) - Power(b1,2)*Power(b2,2)*g2 + Power(b1,3)*b3*g2 - 
           2*Power(b2,3)*g0*g2 + 2*b1*b2*b3*g0*g2 + Power(b3,2)*Power(g0,2)*g2 - 2*b2*b4*Power(g0,2)*g2 + b1*b5*Power(g0,2)*g2 + 
           2*Power(b2,2)*b3*g1*g2 - 4*b1*b2*b4*g1*g2 + 2*Power(b1,2)*b5*g1*g2 - Power(b4,2)*Power(g1,2)*g2 + b3*b5*Power(g1,2)*g2 - 
           b2*Power(b3,2)*Power(g2,2) + Power(b2,2)*b4*Power(g2,2) + b1*b3*b4*Power(g2,2) - b1*b2*b5*Power(g2,2) + 
           Power(b4,2)*g0*Power(g2,2) - b3*b5*g0*Power(g2,2) + b1*Power(b2,3)*g3 - 2*Power(b1,2)*b2*b3*g3 + Power(b1,3)*b4*g3 + 
           Power(b2,2)*b3*g0*g3 - 2*b1*Power(b3,2)*g0*g3 + Power(b1,2)*b5*g0*g3 - b3*b4*Power(g0,2)*g3 + b2*b5*Power(g0,2)*g3 - 
           b2*Power(b3,2)*g1*g3 + Power(b2,2)*b4*g1*g3 + b1*b3*b4*g1*g3 - b1*b2*b5*g1*g3 + Power(b4,2)*g0*g1*g3 - b3*b5*g0*g1*g3 + 
           Power(b3,3)*g2*g3 - 2*b2*b3*b4*g2*g3 + b1*Power(b4,2)*g2*g3 + Power(b2,2)*b5*g2*g3 - b1*b3*b5*g2*g3 - Power(b2,4)*g4 + 
           3*b1*Power(b2,2)*b3*g4 - Power(b1,2)*Power(b3,2)*g4 - 2*Power(b1,2)*b2*b4*g4 + Power(b1,3)*b5*g4 + 
           2*b2*Power(b3,2)*g0*g4 - 2*Power(b2,2)*b4*g0*g4 - 2*b1*b3*b4*g0*g4 + 2*b1*b2*b5*g0*g4 - Power(b4,2)*Power(g0,2)*g4 + 
           b3*b5*Power(g0,2)*g4 - Power(b3,3)*g1*g4 + 2*b2*b3*b4*g1*g4 - b1*Power(b4,2)*g1*g4 - Power(b2,2)*b5*g1*g4 + b1*b3*b5*g1*g4)/den
  
  q5 <- (-Power(b1,5) - 4*Power(b1,3)*b2*g0 - 3*b1*Power(b2,2)*Power(g0,2) - 3*Power(b1,2)*b3*Power(g0,2) - 
           2*b2*b3*Power(g0,3) - 2*b1*b4*Power(g0,3) - b5*Power(g0,4) + 3*Power(b1,2)*Power(b2,2)*g1 - 3*Power(b1,3)*b3*g1 + 
           2*Power(b2,3)*g0*g1 + 2*b1*b2*b3*g0*g1 - 4*Power(b1,2)*b4*g0*g1 + Power(b3,2)*Power(g0,2)*g1 + 2*b2*b4*Power(g0,2)*g1 - 
           3*b1*b5*Power(g0,2)*g1 - Power(b2,2)*b3*Power(g1,2) - 2*b1*Power(b3,2)*Power(g1,2) + 4*b1*b2*b4*Power(g1,2) - 
           Power(b1,2)*b5*Power(g1,2) - 2*b3*b4*g0*Power(g1,2) + 2*b2*b5*g0*Power(g1,2) + Power(b4,2)*Power(g1,3) - 
           b3*b5*Power(g1,3) - 2*b1*Power(b2,3)*g2 + 4*Power(b1,2)*b2*b3*g2 - 2*Power(b1,3)*b4*g2 - 2*Power(b2,2)*b3*g0*g2 + 
           4*b1*Power(b3,2)*g0*g2 - 2*Power(b1,2)*b5*g0*g2 + 2*b3*b4*Power(g0,2)*g2 - 2*b2*b5*Power(g0,2)*g2 + 
           2*b2*Power(b3,2)*g1*g2 - 2*Power(b2,2)*b4*g1*g2 - 2*b1*b3*b4*g1*g2 + 2*b1*b2*b5*g1*g2 - 2*Power(b4,2)*g0*g1*g2 + 
           2*b3*b5*g0*g1*g2 - Power(b3,3)*Power(g2,2) + 2*b2*b3*b4*Power(g2,2) - b1*Power(b4,2)*Power(g2,2) - 
           Power(b2,2)*b5*Power(g2,2) + b1*b3*b5*Power(g2,2) + Power(b2,4)*g3 - 3*b1*Power(b2,2)*b3*g3 + Power(b1,2)*Power(b3,2)*g3 + 
           2*Power(b1,2)*b2*b4*g3 - Power(b1,3)*b5*g3 - 2*b2*Power(b3,2)*g0*g3 + 2*Power(b2,2)*b4*g0*g3 + 2*b1*b3*b4*g0*g3 - 
           2*b1*b2*b5*g0*g3 + Power(b4,2)*Power(g0,2)*g3 - b3*b5*Power(g0,2)*g3 + Power(b3,3)*g1*g3 - 2*b2*b3*b4*g1*g3 + 
           b1*Power(b4,2)*g1*g3 + Power(b2,2)*b5*g1*g3 - b1*b3*b5*g1*g3)/den
  
  
  
  y <- x^al
  
  h.pade <- (p1*y + p2*y^2 + p3*y^3 + p4*y^4 + p5*y^5)/(1 + q1*y + q2*y^2 + q3*y^3 +q4*y^4 + q5*y^5)
  return(h.pade)
}


########################################################################
# Pade approximations to D^\alpha h(a,t)
########################################################################

d.h.Pade22 <- function(params)function(a,t){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  h.pade <- h.Pade22(params)(a,t)
  
  res <- 1/2*(nu*h.pade-rm)*(nu*h.pade-rp)
  
  return(res)
}


########################################################################
d.h.Pade33 <- function(params)function(a,t){
  
  rho <- params$rho
  nu <- params$nu
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  h.pade <- h.Pade33(params)(a,t)
  
  res <- 1/2*(nu*h.pade-rm)*(nu*h.pade-rp)
  
  return(res)
}

########################################################################

d.h.Pade44 <- function(params)function(a,t){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  h.pade <- h.Pade44(params)(a,t)
  
  res <- 1/2*(nu*h.pade-rm)*(nu*h.pade-rp)
  
  return(res)
}

########################################################################

d.h.Pade55 <- function(params)function(a,t){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  h.pade <- h.Pade55(params)(a,t)
  
  res <- 1/2*(nu*h.pade-rm)*(nu*h.pade-rp)
  
  return(res)
}




########################################################################
#
# Characteristic function using Padé approximation
#
########################################################################

phiRoughHestonRational.raw <- function(params, xiCurve, h.approx, n=100) function(a, t) {
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+1/2
  lam <- params$lam
  
  lamp <- lam/nu
  
  lamTilde <-  lamp -(0+1i) * rho * a
  aa <- sqrt(a * (a + (0+1i)) + lamTilde^2)
  rm <- lamTilde - aa
  rp <- lamTilde + aa
  
  ti <- (0:n)/n * t
  h <- h.approx(params)(a,ti)
  dah <- 1/2*(nu*h-rm)*(nu*h-rp)
  g <- dah + lam * h
  xi <- xiCurve(ti)
  conv <- t(g) %*% rev(xi) * t
  cor <- (g[1]*xi[n+1]+g[n+1]*xi[1])/2*t
  psi <- (conv-cor)/n
  return(exp(psi))
  }

phiRoughHestonRational <- function(params, xiCurve, h.approx, n=100) function(a, t){
  phi1 <- function(u){ifelse(u==0,1,phiRoughHestonRational.raw(params, xiCurve, h.approx, n)(u,t))}
  return(sapply(a,phi1))
}
