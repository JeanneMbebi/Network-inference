# This code finds the reactant substrate complexes buy computing the stoichiometric coefficients of the compounds based on the information obtained from all the models.
# Application: Example of Envz OmR system.



load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met22.RData")

n=dim(Met)[2]
Time= c(0.01, 2, 2.5, 2.9, 4, 4.98,5)
S=length(Time)

# II. We include the responses in the data set 
NewRespMet=rep(c(1:n),S)
NewAllMod=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/AllmodelsEnvZ.csv")) 

# III. Computation of the coefficients

m=9
c1=2
c2=3
Coeff=matrix(0, m,n)
Incr=5

RespMet=c(1:n)                                                     # List of response metabolites for each model
FM=NewAllMod[(1:n),(((Incr-1)*n+1):(Incr*n))]                      # Matrix of the five best perfoming models
Num_NZ_PerMod=c()                                                  # List of non zero coefficients numbers per model                                             
ALL_NZ=c()                                                         # List of all compounds involved at the putative time point over all models
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
#x=sum(Num_NZ_PerMod[(1:(j-1))])
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

#cvec=rep(1,m*LL+Q)
cvec=c(rep(1, m*LL),rep(1,Q))
#reaction=2
# b=c(1,rep(0,m-1))
# for (i in 1:LL){
#   cvec[((i-1)*m+reaction):(m*i+(reaction-1))]=b
# }


lb=c(rep(0,m*LL), rep(0,Q))
ub=rep(2,m*LL+Q)
#save(cvec, file="cvec.RData")
write.csv(cvec, file = "cvec.csv")
#save(Amat, file="Amat.RData")
write.csv(Amat, file = "Amat.csv")
#save(bvec, file="bvec.RData")
write.csv(bvec, file = "bvec.csv")
#save(lb, file="lb.RData")
write.csv(lb, file = "lb.csv")
#save(ub, file="ub.RData")
write.csv(ub, file = "ub.csv")

Amat1=rbind(-Amat, diag(dim(Amat)[2]), -diag(dim(Amat)[2]))
bvec1=c(-bvec, lb, -ub)
#library("quadprog")
#Sol1 =solve.QP(Dmat=0*diag(dim(Amat)[2]), dvec=-cvec, Amat=t(Amat1),  bvec=t(bvec1))

library("Rcplex")
Sol1=Rcplex(cvec,Amat,bvec, Qmat=NULL,lb,ub,objsense="min",sense="L")

#c. The solution
Coeff=matrix(0,m, n)
for (t in 1:LL){
  Coeff[,ALL_NZ[t]]=Sol1$xopt[((t-1)*m+1):(t*m)]
}
EpsThetha=Sol1$xopt[(m*LL+1):(m*LL+Q)]

