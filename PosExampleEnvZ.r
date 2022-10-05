
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met22.RData")
P=dim(Met)[2]
N=dim(Met)[1]
N11=3


Time = c(0.01, 2, 2.5, 2.9, 4, 4.98,5)
#Time=seq(0.5,5, by = .5)
S=length(Time)

# a- Replace the "0"s by a very small value
for (i in 1:N){
  for (j in 1:P){
    if (Met[i,j]==0){
      Met[i,j]=runif(1,0,0.1)
    }
  }
} 

NMet=matrix(0, 2*N,P)
NNMet1=matrix(0, 2*(N11),P)
#NNMet1=NMet[((2*N-5):(2*N)),]
NM1=Met[((N-2):N),]
nMet1=matrix(0,2,P)
a=1
x=c(1:3)[-a]
while (a<4){
  for (i in 1:2){
    nMet1[i,]=log(NM1[a,]/NM1[x[i],])
  }
  NNMet1[(((a-1)*2+1):(2*a)),]=nMet1
  #nMet1=NNMet1[(((a-1)*2+1):(2*a)),]
  a=a+1
  x=c(1:3)[-a]
  
}
NMet[((2*N-5):(2*N)),]=NNMet1

NNMet=matrix(0, 2*(N11),P)

NM=Met[1:N11,]
a1=1
while (a1<(S)) {
  
  b=1
  x=c(1:3)[-b]
  nMet=matrix(0, 2,P)
  while (b<4){
    for (i in 1:2){
      nMet[i,]=log(NM[b,]/NM[x[i],])
    }
    NNMet[(((b-1)*2+1):(2*b)),]=nMet
    #nMet=NNMet[(((b-1)*2+1):(2*b)),]
    b=b+1
    x=c(1:3)[-b]
    
  }
  NMet[(((a1-1)*6+1):(6*a1)),(1:P)]=NNMet
  #NNMet=NMet[(((a1-1)*6+1):(6*a1)),]
  NM=Met[(((a1)*3+1):(3*(a1+1))),]
  a1=a1+1
  
}

# 2- Let us standardize the obtained data

AnalMet=scale(NMet, center=TRUE, scale=TRUE)

# 3-Regression Analysis

N1=2*N11
x0=AnalMet[1:N1, ]
x1=AnalMet[(N1+1):(2*N1), ]
CD=c()
Save0=matrix(0,P-1,P)
Save1=matrix(0,P-1,P)
AllTP=matrix(0,S,P-1)
AllMod=matrix(0, P-1, P*S)
CoefDet=matrix(0,S,P)
Resid= c()
RSS=matrix(0,S,P)
#z1=sample(0:5,1)
#z2=sample(0:5,1)
# a-Now the compuation of coefficients at t0
library("chemometrics")
library("penalized")
library("survival")
for (c in 1:P){
  ny0=(-1)*x0[,c]
  nx0=x0[,-c]
  ny1=(-1)*x1[,c]
  nx1=x1[,-c]
  time0=data.frame(ny0,nx0)
  time1=data.frame(ny1,nx1)
  CV0=lassoCV(ny0~nx0, data=time0, K=3, fraction = seq(0.1, 0.5, by = 0.1))
  dev.off()
  #Pen=penalized(ny0, nx0, unpenalized=~0,lambda1=CV0$sopt, lambda2=0, positive=TRUE, data=time0, model="linear")
  Pen=penalized(ny0, nx0, lambda1=CV0$sopt , lambda2=0, positive=TRUE, data=time0, model="linear")
  #C0=lassocoef(ny0~nx0, data=time0, sopt=CV0$sopt, plot.opt=FALSE)
  Coef=coefficients(Pen,"penalized")
  Save0[,c]=Coef
  Rsq=var(nx0%*%Coef)/var(ny0)
  CD[[length(CD)+1]]=Rsq
  Resid[[length(Resid)+1]]=sum((ny0-(nx0%*%Coef))**2)
  
  # b-At time point 2(t1)
  # Compute the new variables
  
  Kappa=c(6,10,2,7,11,4,19,13,8,20,14,9,18,3,17,5,12,16,1,15)
  
  y11 = matrix(0, N1+P-1,1)
  X=matrix(0,(N1+P-1),2*(P-1))
  KapLam=matrix(0, length(Kappa), P+1)
  for (i in 1:N1){
    y11[i]=ny1[i]
  }
  
  for (i in 1:N1){
    for (j in 1:(P-1)){
      X[i,j]=nx1[i,j]
    }
  }
  for (i in 1:(P-1)){
    X[N1+i,P-1+i]=1
  }
  for (i in 1:(length(Kappa))){   
    
    # c-Enter the new Beta11, then the new data frame
    
    Beta00=matrix(0, (2*(P-1)), 1)
    for (j in P:(2*(P-1))){
      Beta00[j]=sqrt(Kappa[i])*(Coef[j-P+1])
    }
    Y=y11+X%*%Beta00
    ntime=data.frame(Y,X)      
    CV1=lassoCV(Y~X,data=ntime, K=10, fraction = seq(0.1, 0.5, by = 0.1))
    #Pen1=penalized(Y, X, unpenalized=~0,lambda1=CV1$sopt, lambda2=0, positive=TRUE, data=ntime, model="linear")
    Pen1=penalized(Y, X,lambda1=CV1$sopt, lambda2=0, positive=TRUE, data=ntime, model="linear")
    Coef1=coefficients(Pen1,"penalized")
    dev.off()
    for (l in 3:(P+1)){
      KapLam[i,1]=(CV1$sopt)*(1+sqrt(Kappa[i]))
      KapLam[i,2]=Kappa[i]
      KapLam[i,l]=Coef1[l-2]
    } 
  }
  CV11=lassoCV(ny1~nx1, data=time1, K=3, fraction = seq(0.1, 0.5, by = 0.1))
  dev.off()
  L=numeric(length(Kappa))
  for (n in 1:length(Kappa)){
    L[n]=abs(CV11$sopt - KapLam[n,2])
  }
  ind1=which.min(L)
  Save1[,c]=KapLam[ind1,][3:(P+1)]
}
AllMod[(1:(P-1)),(1:P)]=Save0
AllMod[(1:(P-1)),((P+1):(2*P))]=Save1

CoefDet[1,]=CD
CoefDet[2,]=CD
RSS[1,]=Resid
RSS[2,]=Resid
ind2=which.max(CD)
AllTP[1,]=Save0[,ind2]
AllTP[2,]=Save1[,ind2]
NM=P
SCD=sort(CD)
NMHCD=SCD[(length(SCD)-NM+1):(length(SCD))]
MHCD0=matrix(0, P+1, NM)
MHCD1=matrix(0, P+1, NM)
ALLTPMHCD=matrix(0, P+1,NM*(length(Time)) )
for (i in 1:length(NMHCD)){
  s=which(CD==(NMHCD[i]))
  MHCD0[1,i]=sample(s,1)
  MHCD1[1,i]=sample(s,1)
  MHCD1[2,i]=NMHCD[i]
  MHCD0[2,i]=NMHCD[i]
  MHCD0[(3:(P+1)),(1:(NM))][,i]=Save0[,sample(s,1)]
  MHCD1[(3:(P+1)),(1:(NM))][,i]=Save1[,sample(s,1)]
}
ALLTPMHCD[(1:(P+1)),(1:NM)]=MHCD0
ALLTPMHCD[(1:(P+1)),((NM+1):(2*NM))]=MHCD1

# Now, let us compute the coefficients for S-2 remaining time points

Incr=2
RespMet=c(ind2,ind2)
x2=AnalMet[(Incr*N1+1):((Incr +1)*N1),]
while (Incr<(S-1)){
  nCD=c()
  nResid=c()
  nSave1=matrix(0,P-1,P)
  nSave2=matrix(0,P-1,P)
  for (c in 1:P){
    nny1=(-1)*x1[,c]
    nnx1=x1[,-c]
    ny2=(-1)*x2[,c]
    nx2=x2[,-c]
    time11=data.frame(nny1,nnx1)
    time2=data.frame(ny2,nx2)
    nCV1=lassoCV(nny1~nnx1, data=time11, K=3, fraction = seq(0.1, 0.5, by = 0.1))
    dev.off()
    #nPen1=penalized(nny1, nnx1,unpenalized=~0, lambda1=nCV1$sopt, lambda2=0, positive=TRUE, data=time11, model="linear")
    nPen1=penalized(nny1, nnx1, lambda1=nCV1$sopt, lambda2=0, positive=TRUE, data=time11, model="linear")
    nCoef1=coefficients(nPen1,"penalized")
    
    #nSave1[,c]=nCoef1
    nRsq=var(nnx1%*%nCoef1)/var(nny1)
    nCD[[length(nCD)+1]]=nRsq
    nResid[[length(nResid)+1]]=sum((nny1-(nnx1%*%nCoef1))**2)
    # At the next time point 
    # Compute the new variables
    
    
    y22 = matrix(0, N1+P-1,1)
    nX=matrix(0,(N1+P-1),2*(P-1))
    nKapLam=matrix(0, length(Kappa), P+1)
    
    for (i in 1:N1){
      y22[i]=ny2[i]
    }
    
    for (i in 1:N1){
      for (j in 1:(P-1)){
        nX[i,j]=nx2[i,j]
      }
    }
    for (i in 1:(P-1)){
      nX[N1+i,P-1+i]=1
    }
    #pdf("test.pdf")
    for (i in 1:(length(Kappa))){  
      
      #Enter the new Beta11, then the new data frame
      
      Beta11=matrix(0, 2*(P-1), 1)
      for (j in P:(2*(P-1))){
        Beta11[j]=sqrt(Kappa[i])*(nCoef1[j-P+1])
      }
      nY=y22+nX%*%Beta11
      nntime=data.frame(nY,nX)      
      
      CV2=lassoCV(nY~nX,data=nntime, K=10, fraction = seq(0.1, 0.5, by = 0.1))
      dev.off()
      #nPen2=penalized(nY, nX,unpenalized=~0, lambda1=CV2$sopt, lambda2=0, positive=TRUE, data=nntime, model="linear")
      nPen2=penalized(nY, nX, lambda1=CV2$sopt, lambda2=0, positive=TRUE, data=nntime, model="linear")
      nCoef2=coefficients(nPen2,"penalized")
      
      
      for (l in 3:(P+1)){
        nKapLam[i,1]=(CV2$sopt)*(1+sqrt(Kappa[i]))
        nKapLam[i,2]=Kappa[i]
        nKapLam[i,l]=nCoef2[l-2]
      }
    }
    #dev.off()
    
    CV22=lassoCV(ny2~nx2, data=time2, K=3, fraction = seq(0.1, 0.5, by = 0.1))
    dev.off()
    nL=numeric(length(Kappa))
    for (n in 1:length(Kappa)){
      nL[n]=abs(CV22$sopt - nKapLam[n,2])
    }
    ind11=which.min(nL)
    nSave2[,c]=nKapLam[ind11,][3:(P+1)]
    
    
  }
  AllMod[(1:(P-1)),(Incr*P+1):((Incr+1)*P)]=nSave2
  nSCD=sort(nCD)
  nNMHCD=nSCD[(length(nSCD)-NM+1):(length(nSCD))]
  nMHCD=matrix(0, P+1, NM)
  for (i in 1:(NM)){
    ns=which(nCD==(nNMHCD[i]))
    nMHCD[1,i]=ns[1]
    nMHCD[2,i]=nNMHCD[i]
    nMHCD[(3:(P+1)),(1:(NM))][,i]=nSave2[,ns[1]]
  }
  ALLTPMHCD[(1:(P+1)),(((Incr)*NM+1):((Incr+1)*NM))]=nMHCD
  
  CoefDet[Incr+1,]=nCD
  RSS[Incr+1,]=nResid
  ind3=which.max(nCD)
  RespMet[length(RespMet)+1]=ind3
  AllTP[Incr+1,]=nSave2[,ind3]
  x1=x2
  Incr=Incr+1
  x2=AnalMet[((Incr*N1+1):((Incr+1)*N1)),]
}
x8=AnalMet[((N1*(S-2)+1):(N1*(S-1))), ]
x9=AnalMet[((N1*(S-1)+1):(N1*(S))), ]
NCD=c()
NResid=c()
nSave3=matrix(0, P-1, P)
#z6=sample(0:5,1)
for (c in 1:P){
  nny1=(-1)*x8[,c]
  nnx1=x8[,-c]
  ny2=(-1)*x9[,c]
  nx2=x9[,-c]
  time11=data.frame(nny1,nnx1)
  time2=data.frame(ny2,nx2)
  nCV1=lassoCV(nny1~nnx1, data=time11, K=3, fraction = seq(0.1, 0.5, by = 0.1))
  dev.off()
  #nPen1=penalized(nny1, nnx1, unpenalized=~0,lambda1=nCV1$sopt, lambda2=0, positive=TRUE, data=time11, model="linear")
  nPen1=penalized(nny1, nnx1, lambda1=nCV1$sopt, lambda2=0, positive=TRUE, data=time11, model="linear")
  nCoef1=coefficients(nPen1,"penalized")
  
  #nSave1[,c]=nCoef1
  nRsq=var(nnx1%*%nCoef1)/var(nny1)
  NCD[[length(NCD)+1]]=nRsq
  NResid[[length(NResid)+1]]=sum((nny1-(nnx1%*%nCoef1))**2)
  
  # At the next time point 
  # Compute the new variables
  
  
  y22 = matrix(0, N1+P-1,1)
  nX=matrix(0,(N1+P-1),2*(P-1))
  nnKapLam=matrix(0, length(Kappa), P+1)
  
  for (i in 1:N1){
    y22[i]=ny2[i]
  }
  
  for (i in 1:N1){
    for (j in 1:(P-1)){
      nX[i,j]=nx2[i,j]
    }
  }
  for (i in 1:(P-1)){
    nX[N1+i,P-1+i]=1
  }
  #pdf("test.pdf")
  for (i in 1:(length(Kappa))){  
    
    #Enter the new Beta11, then the new data frame
    
    Beta11=matrix(0, 2*(P-1), 1)
    for (j in P:(2*(P-1))){
      Beta11[j]=sqrt(Kappa[i])*(nCoef1[j-P+1])
    }
    nY=y22+nX%*%Beta11
    nntime=data.frame(nY,nX)      
    
    CV2=lassoCV(nY~nX,data=nntime, K=10, fraction = seq(0.1, 0.5, by = 0.1))
    dev.off()
    #nPen2=penalized(nY, nX, unpenalized=~0,lambda1=CV2$sopt, lambda2=0, positive=TRUE, data=nntime, model="linear")
    nPen2=penalized(nY, nX,lambda1=CV2$sopt, lambda2=0, positive=TRUE, data=nntime, model="linear")
    nCoef2=coefficients(nPen2,"penalized")
    
    
    for (l in 3:(P+1)){
      nnKapLam[i,1]=(CV2$sopt)*(1+sqrt(Kappa[i]))
      nnKapLam[i,2]=Kappa[i]
      nnKapLam[i,l]=nCoef2[l-2]
    }
  }
  #dev.off()
  
  CV22=lassoCV(ny2~nx2, data=time2, K=3, fraction = seq(0.1, 0.5, by = 0.1))
  dev.off()
  nL=numeric(length(Kappa))
  for (n in 1:length(Kappa)){
    nL[n]=abs(CV22$sopt - nnKapLam[n,2])
  }
  ind11=which.min(nL)
  nSave3[,c]=nnKapLam[ind11,][3:(P+1)]
  
  
}
AllMod[(1:(P-1)),((S-1)*P+1):(S*P)]=nSave3
nSCD1=sort(NCD)
nNMHCD1=nSCD1[(length(nSCD1)-NM+1):(length(nSCD1))]
nMHCD1=matrix(0, P+1, NM)
for (i in 1:(NM)){
  ns=which(NCD==(nNMHCD1[i]))
  nMHCD1[1,i]=ns[1]
  nMHCD1[2,i]=nNMHCD1[i]
  nMHCD1[(3:(P+1)),(1:(NM))][,i]=nSave3[,ns[1]]
}
ALLTPMHCD[(1:(P+1)),(((S-1)*NM+1):(S*NM))]=nMHCD1

CoefDet[S,]=NCD
RSS[S,]=NResid
ind3=which.max(NCD)
RespMet[length(RespMet)+1]=ind3
AllTP[S,]=nSave3[,ind3]
save(AllTP, file="PosAllTPEnvZ.RData")
save(RespMet, file="PosRespMetEnvZ.RData")
save(ALLTPMHCD, file="ALLTPMHCDEnvZ.RData")
save(CoefDet, file="CoefDetPosEnvZ.RData")
save(RSS, file="RSSPosEnvZ.RData")
save(AllMod, file="AllModEnvZ.RData")

# II. Profile of models with R>0.8

i=1
count=c()
while (i<(S+1)){
  V=CoefDet[i,]
  count[[length(count)+1]]=(length(which(V>0.8))) 
  i=i+1
}
plot(Time, count,"l",lwd=5, col = "blue", xlab="Time",ylab="Number of models with R²>0.8")
plot(Time,CoefDet[,5],"l",lwd=4,col = "blue", xlab="Time",ylab="Coeff det(R²) ")


# III. Plot of the coefficients of determination

hist(CoefDet, breaks = 50, col = rgb(1,0,1,0.5), xlab = "R²",ylab="Number of models", main = "")
hist(RSS, breaks = 50, col =rgb(1,0,1,0.5), xlab = "RSS ",ylab="Number of models", main = "")
