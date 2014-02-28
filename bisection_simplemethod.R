
#Functions
bisect <- function(fn, lower,upper, zero=1.e-07) {
  f.lo <- fn(lower) #function must give negative result
  f.hi <- fn(upper) #function must give positive result

  
  if (f.lo * f.hi > 0) stop("interval borders wrong\n")
  chg <- upper - lower #narrowed interval
  
  while (abs(chg) > zero) {  #still not minimized
    x.new <- (lower + upper) / 2 #new input for function (half interval!)
    f.new <- fn(x.new) #solve function
    if (abs(f.new) <= zero)  
      break("solution found")#sufficiently close to zero=solution
    if (f.lo * f.new < 0) upper <- x.new #try harder
    if (f.hi * f.new < 0) lower <- x.new #try harder 
    chg <- upper - lower #narrowed interval go in loop again

  }
  list(solution = x.new, value = f.new)
}


H<-function(psimax,psimin,conc,MIC,kappa) {
((psimax-psimin)*(conc/MIC)^kappa)/((conc/MIC)^kappa-(psimin/psimax))
}

psi<-function(psimax,alpha,psiminA,A,MICA,KA,psiminB,B,MICB,KB){
psimax-H(psimax,psiminA,A,MICA,KA)-H(psimax,psiminB,B,MICB,KB)-alpha*H(psimax,psiminA,A,MICA,KA)*H(psimax,psiminB,B,MICB,KB)
}
#Parameters
MICA=1.27
MICB=0.38
doseA=seq(from=0,to=MICA,length.out=100)
pointsB<-c()
fctvalues<-c()
for (i in doseA ){
  psizero<-function(doseB){
    psi(psimax=0.04672,alpha=0,psiminA=-0.125,A=i,MICA=MICA,KA=0.925,psiminB=-0.145,B=doseB,MICB=MICB,KB=1.23)
  }
  fctvalue<-psizero(i)
  #fn=psi(psimax=0.04672,alpha=0.5,psiminA=-0.125,A=doseA,MICA=1.27,KA=0.925,psiminB=-0.145,B=0.1,MICB=0.38,KB=1.23)
  fctvalues<-c(fctvalues,fctvalue)
  pointsA=bisect(psizero,0,MICB*1e6)
  pointsB=c(pointsB,pointsA$solution)
}

#Plot Isoboles
plot(doseA,pointsB)
#print(pointsB)
print(fctvalues)
#psi(psimax=0.04672,alpha=0.5,psiminA=-0.125,A=0.1,MICA=1.27,KA=0.925,psiminB=-0.145,B=0.1,MICB=0.38,KB=1.23)
#solution1=bisect(fn1, 0.001,1)
#print("Antibiotic A")
#print(solution1)
#print("--------------------------------------------------")
#solution2=bisect(fn2, 1, 0.001)
#print("Antibiotic B")
#print(solution2)
#print("--------------------------------------------------")
#solution3=bisect(fn3, 3, 4)
#print(solution3)
