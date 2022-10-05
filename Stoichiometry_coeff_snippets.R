Stoichcoeff=function(Met, Time, NewAllMod, m, c1, c2, Incr){
  # The parameter "Met" corresponds to the N times P matrix of compounds concentrations at different time points given by the entries of the vector "Time".
  # n=S*x, where x corresponds to the number of initials conditions (replicates) for which the profiles where computed. S is the number of time points.
  # P is the number of compounds.
  # The parameter "NewAllMod" is the P times P*S matrix of all beta coefficients (including the responses with coefficients of 1) at all time points.
  # m is the number of reactions
  # c1 and c2 are as in the main paper
  # Incr is the corresponding time point. To be incremented for each time point
  n=dim(Met)[2]
  S=length(Time)
  
  # III. Computation of the coefficients
  
  Coeff=matrix(0, m,n)
  RespMet=c(1:n)                                                     
  FM=NewAllMod[(1:n),(((Incr-1)*n+1):(Incr*n))]                      
  Num_NZ_PerMod=c()                                                                                       
  ALL_NZ=c()                                                         
  NZ=c()
  
  # a. We identify the number of compounds (non-zero betas) per model for the supposed time point
  
  for (i in 1:n){
    S1=which(FM[,i]!=0)
    S2=length(S1)
    NZ=append(NZ,S1)
    ALL_NZ=unique(sort(NZ))
    Num_NZ_PerMod[[length(Num_NZ_PerMod)+1]]=S2
    
  }
  LL=length(ALL_NZ)
  Q=sum(Num_NZ_PerMod)
  
  # b. We enter the constrained matrix
  
  Amat=matrix(0, 2*m+2*LL+2*Q, m*LL+Q)
  A=matrix(0, m,m*LL)
  B=matrix(0,LL, m*LL)
  for (i in 1:LL){
    A[(1:m),(((i-1)*m +1):(i*m))]=diag(m)
  }
  Amat[(1:m),(1:(m*LL))]=A
  Amat[(m+1):(2*m),(1:(m*LL))]=-A
  for (j in 1:LL){
    B[j,][((j-1)*m+1):(j*m)]=rep(1,m)
    
  }
  Amat[(2*m+1):(2*m+LL),(1:(m*LL))]=B
  Amat[(2*m+LL+1):(2*m+2*LL),(1:(m*LL))]=-B
  Mat=matrix(0, 2*Q, m*LL+Q)
  
  s=1
  T2=FM[,s]
  a1=RespMet[s]
  v1=which(T2!=0)
  L=length(v1)
  Ind=c()
  for (k in 1:L){
    Ind=append(Ind, which(ALL_NZ==v1[k]))
  }
  e1=which(ALL_NZ==a1)
  AA=matrix(0, 2*L, m*LL+Q)
  D=matrix(0, L, m*LL+Q)
  for (t in 1:L){
    if (v1[t]!=a1){
      D[t,][((Ind[t]-1)*m+1):(Ind[t]*m)]=rep(1,m)
      D[t,][((e1[1]-1)*m+1):(e1[1]*m)]=rep(-T2[v1[t]],m)
    }
    else {
      D[t,][((e1[1]-1)*m+1):(e1[1]*m)]=rep((1-T2[v1[t]]),m)
    }
    
  }
  D[(1:L),((m*LL+1):(m*LL+L))]=-diag(L)
  AA[(1:L), (1:(m*LL+Q))]=D
  AA[((L+1):(2*L)), (1:(m*LL+Q))]=-D
  Mat[(1:(2*L)),(1:(m*LL+Q))]=AA 
  
  
  
  
  for (j in 2:n){
    T2=FM[,j]
    a1=RespMet[j]
    v1=which(T2!=0)
    L=length(v1)
    Ind=c()
    for (k in 1:L){
      Ind=append(Ind, which(ALL_NZ==v1[k]))
    }
    e1=which(ALL_NZ==a1)
    AA=matrix(0, 2*L, m*LL+Q)
    D=matrix(0, L, m*LL+Q)
    for (t in 1:L){
      if (v1[t]!=a1){
        D[t,][((Ind[t]-1)*m+1):(Ind[t]*m)]=rep(1,m)
        D[t,][((e1[1]-1)*m+1):(e1[1]*m)]=rep(-T2[v1[t]],m)
      }
      else {
        D[t,][((e1[1]-1)*m+1):(e1[1]*m)]=rep((1-T2[v1[t]]),m)
      }
      
    }
    x=sum(Num_NZ_PerMod[(1:(j-1))])
    D[(1:L),((m*LL+x+1):(m*LL+x+L))]=-diag(L)
    AA[(1:L), (1:(m*LL+Q))]=D
    AA[((L+1):(2*L)), (1:(m*LL+Q))]=-D
    Mat[((2*x+1):(2*x+2*L)),(1:(m*LL+Q))]=AA 
    
  }
  Amat[(2*m+2*LL+1):(2*m+2*LL+2*Q), (1:(m*LL+Q))]=Mat
  bvec=c(rep(c1,m), rep(0,m), rep(c2,LL), rep(-1,LL), rep(0, 2*Q))
  
  cvec=c(rep(1, m*LL),rep(1,Q))
   
  lb=c(rep(0,m*LL), rep(0,Q))
  ub=rep(2,m*LL+Q)
    
  Amat1=rbind(-Amat, diag(dim(Amat)[2]), -diag(dim(Amat)[2]))
  bvec1=c(-bvec, lb, -ub)
    
  library("Rcplex")
  Sol1=Rcplex(cvec,Amat,bvec, Qmat=NULL,lb,ub,objsense="min",sense="L")
  
  #c. The solution
  Coeff=matrix(0,m, n)
  for (t in 1:LL){
    Coeff[,ALL_NZ[t]]=Sol1$xopt[((t-1)*m+1):(t*m)]
  }
  return(Coeff)
}