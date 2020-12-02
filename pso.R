######## Particle Swarm Optimization ########

# Test Function
eggholder<-function(arg){
  x<-arg[1]
  y<-arg[2]
  part.1<- -1*(y+47)*(sin(sqrt(abs((x/2)+y+47))))
  part.2<- -1*x*(sin(sqrt(abs(x-(y+47)))))
 
  return(part.2+part.1)
}
# Search Domain
def.range<-matrix(0,2,2)
# For booth function
def.range[1,]<-  512
def.range[2,]<- -512

update<-function(w,c.1,c.2,gbest,pbest,par.pos,movements,n.of.par,def.range){
  
  dim <-ncol(def.range)
  
  for(l in 1:n.of.par){
    for(m in 1:dim){
      part.1 <-(w*movements[l,m])
      part.2 <-c.1*(runif(1))*(pbest[l,m]-par.pos[l,m])
      part.3 <-c.2*(runif(1))*(par.pos[gbest,m]-par.pos[l,m])
      grad   <-part.1+part.2+part.3
      par.pos[l,m]   <- par.pos[l,m]+grad
      
      # If the new position is out of bounds
      if(par.pos[l,m]>max(def.range[,m]) || par.pos[l,m]<min(def.range[,m])){
        par.pos[l,m]   <- runif(1,min = min(def.range[,m]), max = max(def.range[,m]) )
        movements[l,m] <- abs(par.pos[l,m]-grad) 
      }else{
        movements[l,m] <- grad
      }
      
    }
  }
  result        <-list(par.pos,movements)
  names(result) <-c("par.pos","movements")
  return(result)
}

my.pso<-function(func,n.of.par=50,max.iter=1000,def.range,c.1=2,c.2=2,w.max=0.9,w.min=0.3){
  
  start.time <- Sys.time()
  
  dim       <-ncol(def.range)
  par.pos   <-matrix(0,n.of.par,dim)
  movements <-matrix(runif(n.of.par*dim),n.of.par,dim)
  cand.sol  <-matrix(0,n.of.par,1)
  
  # Random Positions
  for(i in 1:dim){
    par.pos[,i] <-runif(n.of.par, min = min(def.range[,i]), max = max(def.range[,i]))   
  }
  
  # First Solutions
  for(j in 1:n.of.par){
    arg<-c(par.pos[j,])
    cand.sol[j,1]<-func(arg)
  }
  
  # Global and Particle Bests for first solution
  gbest      <-which.min(cand.sol)
  gbest.val  <-cand.sol[gbest]
  pbest      <-cbind(par.pos,cand.sol)
  
  k  <-1
  
  w  <-w.max-(((w.max-w.min)*k)/max.iter)
  out       <-update(w,c.1,c.2,gbest,pbest,par.pos,movements,n.of.par,def.range)
  par.pos   <-out$par.pos
  movements <-out$movements
  
  for(k in 2:max.iter){
    w<-w.max-(((w.max-w.min)*k)/max.iter)
    
    for(j in 1:n.of.par){
      arg<-c(par.pos[j,])
      cand.sol[j,1]<-func(arg)
    }
    
    lbest<-which.min(cand.sol)
    
    if(cand.sol[lbest]<gbest.val){
      gbest.val<-cand.sol[lbest]
      gbest<-lbest
    }
    
    for(s in 1:n.of.par){
      if(pbest[s,3]>cand.sol[s]){
        pbest[s,]<-c(par.pos[s,],cand.sol[s])
      }
    }
    
    out       <-update(w,c.1,c.2,gbest,pbest,par.pos,movements,n.of.par,def.range)
    par.pos   <-out$par.pos
    movements <-out$movements
    
    
  }
  end.time  <-Sys.time()
  calc.time <-end.time - start.time
  
  message("->> Calculation Time     :  " ,calc.time , " Sec")
  message("->> Min Point Value      : " , gbest.val)
  message("->> Solution Coordinates : " , paste0(pbest[gbest,1:dim],collapse = ","))
  
  
}
  
set.seed(NULL)
result<-my.pso(func=eggholder,def.range = def.range,max.iter = 10000)













