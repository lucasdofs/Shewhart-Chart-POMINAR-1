############### Author: LUCAS SALES ############
#Institution:Federal University of Rio Grande do Norte
#e-mail: ldo-sales@hotmail.com


################# Functions ###############################
##########################################################

#This function generates m observations from a POMINAR(1) process
#Theta is the vector of parameters(alpha,beta,lambda,p) of the process POMINAR(1) 
# m is the sample size.
rpominar=function(theta,m){
  alpha=theta[1] 
  beta=theta[2] 
  lambda=theta[3]
  p=theta[4]
  c1=p*alpha+beta*(1-p)
  y=epsilon=NULL
  y[1]= round(lambda/(1-c1),0) 
  for(i in 2:(m+300)){
    epsilon[i]= rpois(1,lambda) 
    u=runif(1) 
    if(u<p) {
      y[i]=rbinom(1,y[i-1],alpha) + epsilon[i] #INAR(1)
    }  else{
      y[i]=rpois(1,y[i-1]*beta) + epsilon[i] #INARCH(1)
    }
  }
  return(y[301:(300+m)])
}


#The conditional maximum log likelihood of POMINAR(1) process.
#This function provides the EML estimate for the parameters,
#the hessian matrix and the value of maximum log likelihood
#x is the data, initial.value is a vector with the initial values for the 
#parameters (alpha,beta,lambda,p)
EML_POMINAR <- function(x,initial.value){
  N<-length(x)
  LL<-function(y){
    alpha<- y[1]
    beta<- y[2]
    lambda<-y[3]
    p<-y[4]
    p1 = rep(NA, times = N)
    p1[1] = 0
    for(t in 2:N){
      soma <- 0
      soma1<-0
      for(k in 0:(min(x[t-1],x[t]))){
        soma = soma + choose(x[t-1],k)*alpha^k*(1-alpha)^(x[t-1]-k)
        *(exp(-lambda)*((lambda^(x[t] - k))))/(factorial(x[t] - k))
      } 
      soma1=soma1+p*soma+(1-p)*(exp(-(beta*x[t-1]+lambda))*
                    (beta*x[t-1]+lambda)^{x[t]})/(factorial(x[t]))
      p1[t] = log(soma1)   
    }
    -sum(p1)   
  }
  lw=c(0.01,0.01,0.01,0)#lw=-valorInf del parametro
  up = c(0.99, 0.99,Inf,1)###up=valorSup del parametro
  optim(initial.value,LL,lower=lw, upper=up,method="L-BFGS-B", hessian=TRUE)   
}


#This function provides the Upper control Limit (UCL) 
#for monitoring the average of POMINAR(1) process.
#b is the number of MC repetions (1000 is a good number in this case)
#theta is the vector of parameters (alpha,beta,lambda, p) of the process.
#m is the sample size to obtain the 0.9973th quantile (5000 is a good sample size) 
#and n is the size of rational subgroups.
UCLpominar=function(b,theta,m,n){
  q1=NULL
  for(i in 1:b){
    aux=matrix(rpominar(theta,n*m),nrow=m)
    dados=rowMeans(aux)
    q1[i]=quantile(dados,probs=c(0.9973)) 
  }
  (q=colMeans(q1))
  
}



#This function provides the ARLS for the propose method
#theta is the vector of parameters (alpha,beta,lambda, p) of the process.
#delta is the shift in mean of the process. 
#UCL is the Upper control limit obtained by the UCLpominar function
#m is the sample size for process monitoring  (2000 is a good sample size for
#delta<1.5, for delta>1.5 500 is a good sample size) 
#and n is the size of rational subgroups.
ARLpominar=function(theta,UCL,delta,m,n){
  alpha=theta[1] 
  beta=theta[2] 
  lambda=theta[3]
  p=theta[4]
  c1=p*alpha+beta*(1-p)
  c2=p*alpha^2+beta^2*(1-p)
  c3=p*alpha*(1-alpha)+beta*(1-p)
  c4=c3+2*lambda*c1
  #Theoretical variance
  sig2hat=(lambda^2*((1-c1)^2 -(1-c2)) + lambda*(1-c1)*(c4+1-c1))/ (( 1-c2)*(1-c1)^2)
  
  for(i in 1:10000){
    aux=matrix(rpominar(theta,n*m),nrow=m)
    dados[i,1:m]=rowMeans(aux)+delta*sqrt(sig2/n)
  }
  for(i in 1:10000){
    #print(i)
    nmaf[i]=which(dados[i,]>UCL)[1]
  }
  (arl=mean(nmaf,na.rm=T))
  return(arl)
}




