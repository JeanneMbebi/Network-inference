# This code finds the reactant substrate complexes buy computing the stoichiometric coefficients of the compounds
# based on the information obtained from all the models.
# Application: Example, Glycine180, Glycine250


# I. We load the matrix containing the results obtained for all models over all time points
# 
# load("/home/onana/Documents/Rcodes/Projectcodes/AllModex.RData")
#load("/home/onana/Documents/Rcodes/Projectcodes/Met1.RData")

# n=dim(Met)[2]
# Time=c(0, 0.011, 2, 2.5, 2.9, 4)
# S=length(Time)
# 
# # II. We include the responses in the data set
# 
# NewRespMet=rep(c(1:n),S)
# NewAllMod=matrix(0,n, S*n)
# 
# for (l in 1:(S*n)){
#   T1=AllMod[,l]
#   Nzeros=which(T1!=0)
#   a=NewRespMet[l]
#   v=sort(c(Nzeros,a))
#   e=which(v==a)
#   f=which(v!=a)
#   L=length(v)
#   if (length(e)==1){
#     NewAllMod[a,l]=1
#     for (i in 1:(length(Nzeros))){
#       if (Nzeros[i]>a){
#         NewAllMod[(Nzeros[i]+1),l]=T1[Nzeros[i]]
#       }
#       if (Nzeros[i]<a){
#         NewAllMod[(Nzeros[i]),l]=T1[Nzeros[i]]
#       }
#     }
#   }
#   else{
#     for (k in 1:(length(Nzeros))){
#      NewAllMod[(Nzeros[k]+1),l]=T1[Nzeros[k]]
#       NewAllMod[a,l]=1
#     }
#   
#     for (i in 1:(length(f))){
#       if (f[i]>a){
#         NewAllMod[(v[f[i]]+1),l]=T1[v[f[i]]]
#       }
#       if (f[i]<a){
#         NewAllMod[(v[f[i]]),l]=T1[v[f[i]]]
#       }
#     }
#   }
# }
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met11.RData")

n=dim(Met)[2]
Time= c(0.01, 2, 2.5, 2.9, 4, 4.98,5)
S=length(Time)

# II. We include the responses in the data set

NewRespMet=rep(c(1:n),S)

NewAllMod=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels11.csv")) 



# III. Computation of the coefficients

m=10
c1=2
c2=7
Coeff=matrix(0, m,n)
Incr=1

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

############################# NOISE #####################################
NewAllMod1_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise1_005.csv")) 
NewAllMod2_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise2_005.csv"))
NewAllMod3_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise3_005.csv"))
NewAllMod4_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise4_005.csv"))
NewAllMod5_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise5_005.csv"))
NewAllMod6_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise6_005.csv"))
NewAllMod7_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise7_005.csv"))
NewAllMod8_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise8_005.csv"))
NewAllMod9_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise9_005.csv"))
NewAllMod10_005=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise10_005.csv"))
NewAllMod1_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise1_05.csv"))
NewAllMod2_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise2_05.csv"))
NewAllMod3_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise3_05.csv"))
NewAllMod4_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise4_05.csv"))
NewAllMod5_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise5_05.csv"))
NewAllMod6_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise6_05.csv"))
NewAllMod7_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise7_05.csv"))
NewAllMod8_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise8_05.csv"))
NewAllMod9_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise9_05.csv"))
NewAllMod10_05=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise10_05.csv"))
NewAllMod1_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise1_1.csv"))
NewAllMod2_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise2_1.csv"))
NewAllMod3_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise3_1.csv"))
NewAllMod4_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise4_1.csv"))
NewAllMod5_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise5_1.csv"))
NewAllMod6_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise6_1.csv"))
NewAllMod7_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise7_1.csv"))
NewAllMod8_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise8_1.csv"))
NewAllMod9_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise9_1.csv"))
NewAllMod10_1=as.matrix(read.table("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Allmodels_noise10_1.csv"))
NewAllMod_all=cbind(NewAllMod1_005,NewAllMod2_005, NewAllMod3_005, NewAllMod4_005, NewAllMod5_005, NewAllMod6_005, NewAllMod7_005, NewAllMod8_005, NewAllMod9_005, NewAllMod10_005, NewAllMod1_05, NewAllMod2_05, NewAllMod3_05, NewAllMod4_05, NewAllMod5_05, NewAllMod6_05, NewAllMod7_05, NewAllMod8_05, NewAllMod9_05, NewAllMod10_05, NewAllMod1_1, NewAllMod2_1, NewAllMod3_1, NewAllMod4_1, NewAllMod5_1, NewAllMod6_1, NewAllMod7_1, NewAllMod8_1, NewAllMod9_1, NewAllMod10_1 )


load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy1_005.RData") 
Met1_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy2_005.RData")
Met2_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy3_005.RData")
Met3_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy4_005.RData")
Met4_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy5_005.RData")
Met5_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy6_005.RData")
Met6_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy7_005.RData")
Met7_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy8_005.RData")
Met8_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy9_005.RData")
Met9_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy10_005.RData")
Met10_005= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy1_05.RData")
Met1_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy2_05.RData")
Met2_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy3_05.RData")
Met3_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy4_05.RData")
Met4_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy5_05.RData")
Met5_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy6_05.RData")
Met6_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy7_05.RData")
Met7_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy8_05.RData")
Met8_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy9_05.RData")
Met9_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy10_05.RData")
Met10_05= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy1_1.RData")
Met1_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy2_1.RData")
Met2_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy3_1.RData")
Met3_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy4_1.RData")
Met4_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy5_1.RData")
Met5_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy6_1.RData")
Met6_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy7_1.RData")
Met7_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy8_1.RData")
Met8_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy9_1.RData")
Met9_1= Met_noisy1
load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met_noisy10_1.RData")
Met10_1= Met_noisy1
Met_all=cbind(Met1_005, Met2_005, Met3_005, Met4_005, Met5_005, Met6_005, Met7_005, Met8_005, Met9_005, Met10_005, Met1_05, Met2_05, Met3_05, Met4_05, Met5_05, Met6_05, Met7_05, Met8_05, Met9_05, Met10_05, Met1_1, Met2_1, Met3_1, Met4_1, Met5_1, Met6_1, Met7_1, Met8_1, Met9_1, Met10_1 )
 







f =dim(NewAllMod9_1)[1]
g=dim(NewAllMod9_1)[2]


Met=Met1_005
d=dim(Met)[1]
e=dim(Met)[2]
  
#load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/BASF_update/Plos one/codes/Met11.RData")

NewAllMod = NewAllMod_all[1:f, 1:g]
c=1
Stoich_L2norm=c()
while (c<31){
  Coeff_noise=matrix(0, m,n)
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
  cvec=c(rep(1, m*LL),rep(1,Q))
   
  lb=c(rep(0,m*LL), rep(0,Q))
  ub=rep(2,m*LL+Q)
  Sol1=Rcplex(cvec,Amat,bvec, Qmat=NULL,lb,ub,objsense="min",sense="L")
  
  for (w in 1:LL){
    Coeff_noise[,ALL_NZ[t]]=Sol1$xopt[((w-1)*m+1):(w*m)]
  }
  Stoich_L2norm=c(Stoich_L2norm, norm(Coeff, Coeff_noise, mode = "L2"))
  
  NewAllMod=NewAllMod_all[1:f, ((c-1)*g+1):(c*g)]
  Met=Met_all[1:d, ((c-1)*e+1):(c*e)]
  c=c+1
  
}
noises=c(rep(0.05, 10), rep(0.5, 10), rep(1, 10))
data_noise=data.frame(Stoich_L2norm, noises)
boxplot(Stoich_L2norm~noises,data=data_noise, main=expression(paste("Effect of noise in stoichiometric matrix at ", t[1], "=0.01 ")), xlab="Noise", ylab=expression(paste(L[2], "-norm between stoichiometric matrix with and without noise")))

