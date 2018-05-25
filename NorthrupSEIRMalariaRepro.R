rm(list=ls())

#REPRODUCES FIGURE 1,3,4 FROM SEIR PAPER

require(deSolve)

XI <- function(Nh,Nv){
  num <- c_vh*a_v*mu_2h*(lambda_v - mu_v)*Nv
  denom <- mu_2v*mu_v(lambda_h - mu_h)*Nh
  ans <- num/denom
  return(ans)
}

SEIR<-function(t,y,p){
  u = y[1];
  w = y[2];
  R = y[3];
  x = y[4];
  z = y[5];
  #Nh = y[6];
  #Nv = y[7];
  with(as.list(p), {
    du.dtau = lambda*(1-u) + beta*R + r*w + gamma*w*u - xi*u*z;
    dw.dtau = nu*(1-u-R) + gamma*w^2 - (r+alpha + gamma + lambda + nu)*w;
    dR.dtau = alpha*w + gamma*w*R - (beta+lambda)*R;
    dx.dtau = a*(1-x) - b*x*w - c*x*R;
    dz.dtau = e*(1-x) - (a+e)*z;
    #dNh.dtau = (lambda - epsilon)*(1-Nh)*Nh - gamma*Nh*w;
    #dNv.dtau = (a-1)*(1-Nv)*Nv;
    return(list(c(du.dtau,dw.dtau,dR.dtau,dx.dtau,dz.dtau)));
  })
}

gamma = 0.001
xi = 4.066
alpha = 0.3
r = 0.2
beta = 0.35
lambda = 0.00184
nu = 2.0
b = 10
c = 1
epsilon = lambda - 0.005
a = 1.002
e = 2.4

u0 = 0.85
w0 = 0.05
R0 = 0.04
x0 = 0.999999
z0 = 0.000001


p = c(gamma = gamma,
      xi = xi,
      alpha = alpha,
      r = r,
      beta = beta,
      lambda = lambda,
      nu = nu,
      b = b,
      c = c,
      epsilon = epsilon,
      a = a,
      e = e)

y0 = c(u = u0,
       w = w0,
       R = R0,
       x = x0,
       z = z0)


step = .1; #it's useful to define the number of steps separately
t = seq(from=0,to=500,by=step);

out = ode(y=y0,times=t,func=SEIR,parms=p)

plot(out[,1],out[,2],type="l",lwd=2,xlab="Time",ylab="Fraction Sus.",main="Figure1'");
lines(out[,1],out[,3],type="l",lwd=2,col='red');
lines(out[,1],out[,4],type="l",lwd=2,col='blue');

#################
#################

#################
#################
end <- nrow(out)
y0 <- out[end,2:6]

XI <- function(Nh,Nv,xi){
  ans = xi*(Nv/Nh)
  return(ans)
}

SEIRXI<-function(t,y,p){
  u = y[1];
  w = y[2];
  R = y[3];
  x = y[4];
  z = y[5];
  Nh = y[6];
  Nv = y[7];
  with(as.list(p), {
    du.dtau = lambda*(1-u) + beta*R + r*w + gamma*w*u - XI(Nh,Nv,xi)*u*z;
    dw.dtau = nu*(1-u-R) + gamma*w^2 - (r+alpha + gamma + lambda + nu)*w;
    dR.dtau = alpha*w + gamma*w*R - (beta+lambda)*R;
    dx.dtau = a*(1-x) - b*x*w - c*x*R;
    dz.dtau = e*(1-x) - (a+e)*z;
    dNh.dtau = (lambda - epsilon)*(1-Nh)*Nh - gamma*Nh*w;
    dNv.dtau = (a-1)*(1-Nv)*Nv;
    return(list(c(du.dtau,dw.dtau,dR.dtau,dx.dtau,dz.dtau,dNh.dtau,dNv.dtau)));
  })
}

gamma = 0.001
xi = 4.066
alpha = 0.3
r = 0.2
beta = 0.35
lambda = 0.00184
nu = 2.0
b = 10
c = 1
epsilon = lambda - 0.005
a = 1.002
e = 2.4

u0 = y0[1]
w0 = y0[2]
R0 = y0[3]
x0 = y0[4]
z0 = y0[5]
Nh0 = 0.01
Nv0 = 0.01


p = c(gamma = gamma,
      xi = xi,
      alpha = alpha,
      r = r,
      beta = beta,
      lambda = lambda,
      nu = nu,
      b = b,
      c = c,
      epsilon = epsilon,
      a = a,
      e = e)

y0 = c(u = u0,
       w = w0,
       R = R0,
       x = x0,
       z = z0,
       Nh = Nh0,
       Nv = Nv0)


step = .1; #it's useful to define the number of steps separately
t = seq(from=0,to=20000,by=step);

out = ode(y=y0,times=t,func=SEIRXI,parms=p)

plot(out[,1],out[,7],type="l",lwd=2,col='green',xlab="Time",ylab="Fraction",main="Figure3");
lines(out[,1],out[,3],type="l",lwd=2,col='red');
lines(out[,1],out[,4],type="l",lwd=2,col='blue');
lines(out[,1],out[,2],type="l",lwd=2,col='black');

plot(out[,1],out[,8],type = "l",lwd=2,col='green',xlab="Time",ylab="Fraction",main="Figure4");
lines(out[,1],out[,5],type="l",lwd=2,col='red');
lines(out[,1],out[,6],type="l",lwd=2,col='blue');
