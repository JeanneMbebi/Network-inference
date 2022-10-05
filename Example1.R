N=7
#IC=matrix(0,3, N)
 #IC[1,]=c(1, 0,1, 1,2, 2, 2)
 #IC[2,]=c(1, 0, 0, 1, 1, 0,0)
 #IC[3,]=c(1, 2, 0, 1, 2, 2,0)
# IC[2,6]=1
# IC[3,5]=1
# IC[3,3]=1
# IC[3,7]=1
# 
# 
for (t in 1:3){
  for (s in 1:N){
    IC[t,s]=sample(0:2,1)
  }
}
#load("/home/mpimp-golm.mpg.de/onana/Desktop/Documents/BASF_Project/Codes/IC1.RData")
#times=c(0.01, 2, 2.5, 2.9, 4, 4.98)
#times=seq(0.5,5, by = .5)
times=c(0.01, 2, 2.5, 2.9, 4, 4.98,5)
P=length(times)
Met=matrix(0, (P-1)*3, N)
Prof=matrix(0, N*3,(P-1))
library("deSolve")
for (i in 1:3){
RPMod = function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    dA= -2*k_2*A*A+(2*k_1+k_8)*B-(k_7+k_10)*A*C-(k_11+k_13)*A*D +k_9*C*C +k_12*E + k_14*f
    dB= k_2*A*A-(k_1+k_3+k_8)*B + k_4*C*C + k_7*A*C
    dC= (2*k_3+k_7)*B + 2*k_6*D-(k_8-k_10)*A*C - k_9*C*C -k_17*C*E+k_18*G*G
    dD= k_5*C*C -k_6*D-(k_11+k_13)*A*D +k_12*E +k_14*f
    dE= k_11*A*D-(k_12+k_15)*E + k_16*f-k_17*C*E +k_18*G*G
    df= k_13*A*D-(k_14+k_16)*f + k_15*E 
    dG= 2*k_17*C*E -2*k_18*G*G
    
    return(list(c(dA, dB, dC, dD, dE,df, dG)))#, dH, dI, dJ)))
    #return(list(c(dA, dB, dC, dD)))#, dH, dI, dJ)))
  })
}
Pars=c(k_1=1, k_2=1, k_3=1, k_4=1, k_5=1, k_6=1 ,k_7=1 ,k_8=1, k_9=1 , k_10=1, k_11 =2, k_12 =2,k_13 =2 ,k_14 =2, k_15 =2,k_16 =2,k_17 = 1,k_18 = 1 )

#Pars=c(k_1=1, k_2=1, k_3=1 
#c(k_1=1 , k_2=1, k_3=1, k_4=1, k_5=1, k_6=1, k_7=1, k_8=1, k_9=1, k_10=1, k_11=sample(1:3, 1), k_12=sample(1:3, 1), k_13=sample(1:3, 1), k_14=sample(1:3, 1), k_15=sample(1:3, 1), k_16=sample(1:3, 1), k_17=sample(1:3, 1), k_18=sample(1:3, 1))

#yini=c(A=IC[i,1], B=IC[i,2], C=IC[i,3], D=IC[i,4], E=IC[i,5], FF=IC[i,6], G=IC[i,7], H=IC[i,8], I=IC[i,9], J=IC[i,10])
yini=c(A=IC[i,1], B=IC[i,2], C=IC[i,3], D=IC[i,4], E=IC[i,5], f=IC[i,6],G=IC[i,7])

out=ode(yini, times, RPMod, Pars)
for (j in 1:(P-1)){
  Met[((j-1)*3+i),]=out[-1,-1][j,]
 
}
for (k in 1:N){
  Prof[((k-1)*3+i),]=out[-1,-1][,k]
}
}
save(Met, file="Met1.RData")
save(Prof, file="Prof1.RData")
save(IC, file="IC1.RData")
save(Pars, file="Reactrate1.RData")
Ti=seq(1:P)
a=1
opts = c(1,2,3,4,5,6,1,2,3,4) 
cols=c("red", "blue", "green", "yellow", "orange", "cyan", "magenta", "darkblue", "darkorange", "darkgreen")
#png('Profiles1.png')
plot(c(2,4.98),c(min(Met),max(Met)),type="n",xlab="Time",ylab="Compounds profile")
while (a<=N){
  R=Prof[((a-1)*3+1):(3*a),1:(P-1)]
  
  for (k in 1:3){
    lines(times[-1], R[k,],lty=opts[a],col=cols[a],lwd=2 )
  }
  a=a+1
}
legend(2,max(Met), c("A","B","C","D","E","F","G"),lty=c(1,1),lwd=3, col=c("red", "blue", "green", "yellow", "orange", "cyan", "magenta"), angle=c(90,0),cex=1, x.intersp=0.5,bty = "n")
#legend(2,6, c("A","B","C","D","E"),lty=c(1,1),lwd=3, col=c("red", "blue", "green", "yellow", "orange"), angle=c(90,0),cex=1, x.intersp=0.5,bty = "n")
#dev.off()
