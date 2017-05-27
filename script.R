# Code to compare different resampling methods in the conext of regression
# Written by Paula Rodriguez y Felipe Gonzalez

#IMPORTS
library('rrcov') #for the minimum ellipsoid volume
library('MASS') #Multivariate normal
library('gridExtra') #Export tables
library('grid')
library('gtable')

# Funciton that returns the design matrix. Includes the the constant value (column of ones) untill the 
# cubic parameter (beta3)
# This code is based on the script received as part of the project's intruction
# PARAMETER: n (inetegr) the row size of the design matrix
# RETURNS: The desgin matrix of size n x 4
get_designMatrix = function(n)
{
    n=90
    unos=rep(1,n)
    x1=seq(0,10,length=n)
    x2=rep(c(-1,-.5,0,.5,1),length=n)
    tmp=seq(0,5,length=10)
    x3=rep(tmp,length=n)
    X=cbind(unos,x1,x2,x3)
    return(X)
}

# Function that returns a fixed value for the betas
# RETURNS: Vector of size 4, with the beta values
getFixedBetas = function()
{
  beta0=sqrt(2)
  beta1=1
  beta2=-2
  beta3=pi
  beta=c(beta0,beta1,beta2,beta3)
  return(beta)
}

# Function that gives the least square approximation solution of the system Xx = b.
# PARAMETER: X ( n x m matrix) design matrix.
#            b ( n x 1 vector) target or solution vector
#            A ( m x m matrix) corresponds to the inverse of (X^t)(X)
# RETURNS: An m vector to the least approximation solution of the mentioned system
# Recall that the least square approximation corresppons to coeeficinets the projection 
# of b onto the column space of X, that can be obtaioned solving (X^t)(X)x = (X^t)b
# This code is based on the script received as part of the project's intruction
estibeta=function(b,X,A = NULL){
  
  if(is.null(A))
  {
    A=solve(t(X)%*%X)
    print('calculating inverse: (X^t)(X)')
  }

  tmp=A%*%t(X)
  return(tmp%*%b)
} 


estiY=function(Y,X,A){
  n = nrow(X)
  betahat=estibeta(Y,X,A)
  Yhat=X%*%betahat
  resi=Y-Yhat
  tmp=X%*%A
  tmp=tmp%*%t(X)
  iden=diag(nrow=n)
  H=iden-tmp
  h=diag(H)
  resiadj=resi/sqrt(h)
  return(cbind(Yhat,resi,sqrt(h),resiadj))
}


#Function that constructs the target vector with homoscedastic error given the design matrix
# and the beta values
# PARAMETER: X ( n x m matrix) design matrix.
#          : beta ( m vector) vector with the values for betas in the design matrix
#          : sigma (number > 0) corresponds to the standard deviation of the normal samble
# RETURNS: an n vector corresponding to (X)beta + error
# This code is based on the script received as part of the project's intruction
getHomoscedasticError = function(X, beta, sigma = 0.5)
{
  n = nrow(X)
  Y=X%*%beta+rnorm(n, sd = sigma)
  return(Y)
}

#AUXILIARY FUNCTION
# Generates an heteroscedastic non gaussian exponential vecor error of size n
rdexp=function(n){
  y=rexp(n)
  u=runif(n)
  y[u>0.5]=-y[u>0.5]
  return(y)
} # fin de rdexp

#Function that constructs the target vector with heteroscedastic error given the design matrix
# and the beta values
# PARAMETER: X ( n x m matrix) design matrix.
#          : beta ( m vector) vector with the values for betas in the design matrix
#          : sigma (number > 0) corresponds to the standard deviation of the normal samble
# RETURNS: an n vector corresponding to (X)beta + error
# This code is based on the script received as part of the project's intruction
getHeteroscedasticError = function(X, beta, sigma = 0.5)
{
  n = nrow(X)
  Y = getHomoscedasticError(X,beta,sigma)
  Y=X%*%beta+0.025*Y*rdexp(n)
  
  return(Y)
}



# Function that excecutes the ADJUSTED RESIDUE RESAMPLING procedure for the given parameters.
# PARAMETER: X ( n x m matrix) design matrix.
#            Y ( n vector) target vector for the regression 
#            B (integer) size of the reampling
#            smooth_sample (boolean) if the bootstrap sampling should be done smoothly or not
# RETURNS: n vector with star{beta} - hat{beta}
# Given a design matrix X and a target vector Y (probably with error), the procedure first calculates
# hat{beta} as the least square approximation of the system Xx = Y. Then does B adjusted residue
# resampling procedures obtaining star{X} and star{beta} as the solution of (star{X})x = Y. The procedure
# returns the B observations of the vectors: star{beta} - hat{beta}
# This code is based on the script received as part of the project's intruction
adjResidueResampling = function(X,Y,B, smooth_sample = FALSE)
{
  
  n = nrow(X)
  p = 1
  # Formula based on A.W. Boman et al.
  hSmooth = (4/((2+p)*n))^(1/(4+p))
  
  A = solve(t(X)%*%X)
  betahat=estibeta(Y,X,A)
  tmp=estiY(Y,X,A)
  Yhat=tmp[,1];
  resiadj=tmp[,2]
  #output
  refer=matrix(0,B,4)
  for(b in 1:B){
    if(smooth_sample){
      resampledResidues = sample(resiadj,n,replace=T)  
      Ystar=X%*%betahat + rnorm(n, mean = resampledResidues, sd = hSmooth)  
    }else{
      Ystar=X%*%betahat + sample(resiadj,n,replace=T)  
    }
    
    betastar=estibeta(Ystar,X,A)
    refer[b,]=betastar-betahat
  }#end for  
  return(refer)
  
}

# Function that excecutes the PAIR RESAMPLING procedure for the given parameters.
# PARAMETER: X ( n x m matrix) design matrix.
#            Y ( n vector) target vector for the regression 
#            B (integer) size of the reampling
#            smooth_sample (boolean) if the bootstrap sampling should be done smoothly or not
# RETURNS: n vector with star{beta} - hat{beta}
# Given a design matrix X and a target vector Y (probably with error), the procedure first calculates
# hat{beta} as the least square approximation of the system Xx = Y. Then does B pair
# resampling procedures obtaining star{X} and star{beta} as the solution of (star{X})x = Y. The procedure
# returns the B observations of the vectors: star{beta} - hat{beta}
# This code is based on the script received as part of the project's intruction
pairResampling = function(X,Y,B, smooth_sample = FALSE)
{
  
  n = nrow(X)
  m = ncol(X)
  
  p = 1
  # Formula based on: A.W. Boman et al.
  hSmooth1 = (4/((2+p)*n))^(1/(4+p))
  
  p = 4
  # Formula based on: A.W. Boman et al.
  hSmooth4 = (4/((2+p)*n))^(1/(4+p))
  
  # Smoothing standard deviation
  sd =  hSmooth4*diag(rep(1,m))
  
  A = solve(t(X)%*%X)
  betahat=estibeta(Y,X,A)
  
  indi=1:n
  referpa=matrix(0,B,m)
  for(b in 1:B){
    indstar=sample(indi,n,replace=T)
    
    X_resampled = X[indstar,]
    Y_resampled = Y[indstar]
    if(smooth_sample){
      Xstar = t(apply(X_resampled,1,function(row)  mvrnorm(n = 1,mu = row, Sigma = sd)))
      Ystar = rnorm(n = n , mean = Y_resampled, sd = hSmooth1)  
    }else{
      Xstar = X[indstar,]
      Ystar = Y[indstar]   
    }
    
    Astar=solve(t(Xstar)%*%Xstar)
    betastar=estibeta(Ystar,Xstar,Astar)
    referpa[b,]=betastar-betahat
  } #end for  
  return(referpa)
  
}



# Function that excecutes the instructed experiments and stores its results in a data frame
excecuteExperiments = function(print_progress = TRUE,  num_ite = 10, B = 10000)
{
 
  #constants
  pairs = 'Parejas'
  residue = 'Residuos'
  hetero = 'Heterosedastico'
  homo = 'Homosedastico'
  
  
  #Resulting data frame
  result = data.frame(size = c(), typeOfError = c(), smooth = c(), devStandar = c(), typeOfReesampling = c(), mahaDistance = c(), betasInsideEllipse = c(), ellipseVolume = c() )
  
  size = c(90,270)
  devStandar = c(0.5,2)
  typeOfError = c(hetero,homo)
  smooth = c(TRUE, FALSE)
  typeOfReesampling = c(pairs,residue)

  
  variables = expand.grid(smooth,typeOfReesampling,devStandar,typeOfError,size)
  colnames(variables) = c('smooth','typeOfReesampling','devStandar','typeOfError','size')
  
  for(i in 1:nrow(variables))
  {
    vol = c()
    mahaDistance =  c()
    inside = c()
    row = variables[i,]
    for(j in 1:num_ite)
    {
    
      X = get_designMatrix(row$size)
      A = solve(t(X)%*%X)
      beta = getFixedBetas()
      if(row$typeOfError == hetero){
        Y = getHeteroscedasticError(X,beta, row$devStandar)
      }else if(row$typeOfError == homo){
        Y = getHomoscedasticError(X,beta, row$devStandar)
      }else{
        stop(paste('Type of error not supported: ', row$typeOfError))
      }
      
      if(row$typeOfReesampling == pairs){
        
        resampledResponse = pairResampling(X = X,Y = Y, B = B, smooth_sample = row$smooth)
        
      }else if(row$typeOfReesampling == residue){
        
        resampledResponse = adjResidueResampling(X = X,Y = Y, B = B, smooth_sample = row$smooth)
        
      }else{
        stop(paste('Type of reesampling not supported: ', row$typeOfReesampling))
      }
       
      #Gets the ellipse
      ellipse = CovMve(resampledResponse[,3:4],alpha=.95)
      vol = c(vol,getDet(ellipse))
      mahaDistance =  c(mahaDistance,(beta[3:4] - getCenter(ellipse))%*%getCov(ellipse)%*%(beta[3:4] - getCenter(ellipse)))
      inside = c(inside,as.numeric(mahaDistance < 1))
    
    }
    
    temp = data.frame(mahaDistance = mean(mahaDistance), betasInsideEllipse = mean(inside), ellipseVolume = mean(vol))
    
    #binds variables with results
    final_row = cbind(row,temp)
    #organizes the row
    final_row = final_row[,c('size','typeOfError','smooth','devStandar', 'typeOfReesampling','mahaDistance','betasInsideEllipse', 'ellipseVolume')]
    result = rbind(result,final_row )
    if(print_progress)
    {
      print(paste('Finished:', toString(i),'of',nrow(variables), sep = ' '))
      print(final_row)
    }
          
  }
  
  return(result)
  
}

a = excecuteExperiments(print_progress = TRUE,  num_ite = 200, B = 500)
a_temp = a[order(-a$betasInsideEllipse,a$ellipseVolume),]
colnames(a_temp) = c('Tam. Muestra','Tipo Error','Suavizado','Desv. Estan', 'Tipo Remuestreo','Dist. Mahalanobis','% Betas dentro Elipse', 'Volumen Elipse')
pdf('results.pdf', height=13, width=13)

t1 <- tableGrob(a_temp)
title <- textGrob('Results')
padding <- unit(5,"mm")

table <- gtable_add_rows(t1, 
                         heights = grobHeight(title) + padding,
                         pos = 0)
table <- gtable_add_grob(table, title, 1, 1, 1, ncol(table))

grid.newpage()
grid.draw(table)
dev.off()



