library(irlba)
library(RSpectra)


#----------------------------initialization-------------------------------------
init.func <- function(N, D, Ti, node.exist, seed){
  set.seed(seed)
  alpha <- runif(N, -2, -1)
  beta <- runif(Ti, -2, 2)
  Z <- list()
  for (t in 1:Ti) {
    Nt <- node.exist.num[t]
    Jmat <- diag(rep(1,Nt)) - matrix(1, Nt, Nt)/Nt
    Z0 <- matrix(runif(Nt*D, -0.5, 0.5), Nt, D)
    Z0 <- Jmat%*%Z0
    Z0 <- qr.Q(qr(Z0))
    Z0 <- Z0*sqrt(Nt*D/5)/norm(Z0,'F')
    Z[[t]] <- Z0
  }
  return(list(Z=Z, alpha=alpha, beta=beta))
}



#--------------------------------optimization-----------------------------------
opt.func <- function(Z.ini, alpha.ini, beta.ini, A, X, node.exist, nu=0.5, step.size=0.2, epsi=0.0005, niter=500, trace=0){

  N <- length(alpha.ini)
  Ti <- length(beta.ini)
  D <- dim(Z.ini[[1]])[2]
  
  # step size
  step.size.z <- step.size/(Ti*N*D)
  step.size.alpha <- step.size/(Ti*N)
  step.size.beta <- step.size/sum(unlist(lapply(1:Ti, function(t) norm(X[[t]],'F')^2)))
  
  
  diff <- 1
  iter <- 0
  obj <- -Inf
  
  Z <- list()
  alpha <- rep(NA, N)
  beta <- rep(NA, Ti)
  
  while (diff>epsi & iter<niter) {
    Theta.hat <- lapply(1:Ti, function(t) alpha.ini[node.exist[[t]]] %*% t(rep(1,length(node.exist[[t]]))) + 
                          rep(1,length(node.exist[[t]])) %*% t(alpha.ini[node.exist[[t]]]) + 
                          Z.ini[[t]] %*% t(Z.ini[[t]]) + beta.ini[t]*X[[t]])
     
    obj1 <- lapply(1:Ti, function(t) (A[[t]] * Theta.hat[[t]] - exp(Theta.hat[[t]])) * 
                     (matrix(rep(1, length(node.exist[[t]])), length(node.exist[[t]]), length(node.exist[[t]])) - 
                        diag(1, length(node.exist[[t]]), length(node.exist[[t]]))))
    tmp.obj <- sum(unlist(lapply(1:Ti, function(t) sum(obj1[[t]]))))
    obj <- c(obj,tmp.obj)
    if(trace==1){
      print(iter)
      print(diff)
      print(tmp.obj)
    }
    
    if(obj[length(obj)] < obj[length(obj)-1]){
      step.size.z <- step.size.z * nu
      step.size.alpha <- step.size.alpha * nu
      step.size.beta <- step.size.beta * nu
    }
    
    alpha.grad <- rep(0, N)
    for (t in 1:Ti) {
      Nt <- length(node.exist[[t]])
      Theta.t <- alpha.ini[node.exist[[t]]] %*% t(rep(1,Nt)) + 
        rep(1,Nt) %*% t(alpha.ini[node.exist[[t]]]) + 
        Z.ini[[t]] %*% t(Z.ini[[t]]) + beta.ini[t]*X[[t]]
      
      Omega.t <- A[[t]] - exp(Theta.t)
      Omega.t[Omega.t>2] <- 2  
      Omega.t[Omega.t<-2] <- -2 
      diag(Omega.t) <- 0
      
      Z.grad <- 2 * Omega.t %*% Z.ini[[t]]
      alpha.grad[node.exist[[t]]] <- alpha.grad[node.exist[[t]]] + 2 * Omega.t %*% rep(1,Nt)
      beta.grad <- sum(Omega.t * X[[t]])
      
      
      Z[[t]] <- Z.ini[[t]] + step.size.z * Z.grad
      
      Jmat <- diag(rep(1,Nt)) - matrix(1, Nt, Nt)/Nt
      Z[[t]] <- Jmat %*% Z[[t]]
      Z[[t]] <- qr.Q(qr(Z[[t]]))
      Z[[t]] <- Z[[t]]*sqrt(Nt*D/5)/norm(Z[[t]],'F')
      
      beta[t] <- beta.ini[t] + step.size.beta * beta.grad
      
    }
    alpha <- alpha.ini + step.size.alpha * alpha.grad
    
    Z.diff.list <- c()
    for (t in 1:Ti) {
      dec <- svd(t(Z[[t]]) %*% Z.ini[[t]])
      O <- dec$v %*% t(dec$u)
      re0 <- norm(Z[[t]] - Z.ini[[t]] %*% O,'F')^2 / norm(Z.ini[[t]], 'F')^2
      Z.diff.list <- c(Z.diff.list, re0)
    }
    Z.diff <- mean(Z.diff.list)
    
    
    alpha.diff <- sum((alpha-alpha.ini)^2)/sum(alpha.ini^2)
    beta.diff <- sum((beta-beta.ini)^2)/sum(beta.ini^2)
    
    diff <- max(Z.diff, alpha.diff, beta.diff)
    
    iter <- iter+1
    
    # update
    Z.ini <- Z
    alpha.ini <- alpha
    beta.ini <- beta
  }
  
  Theta.hat <- lapply(1:Ti, function(t) alpha[node.exist[[t]]] %*% t(rep(1,length(node.exist[[t]]))) + rep(1,length(node.exist[[t]])) %*% t(alpha[node.exist[[t]]]) + Z[[t]] %*% t(Z[[t]]) + beta[t]*X[[t]])
  obj1 <- lapply(1:Ti, function(t) (A[[t]] * Theta.hat[[t]] - exp(Theta.hat[[t]])) * (matrix(rep(1, length(node.exist[[t]])*length(node.exist[[t]])), length(node.exist[[t]]), length(node.exist[[t]])) - diag(1, length(node.exist[[t]]), length(node.exist[[t]]))))
  tmp.obj <- sum(unlist(lapply(1:Ti, function(t) sum(obj1[[t]]))))
  obj <- c(obj, tmp.obj)
  
  
  if (iter==niter)
  {cat("Warning! May Not Convergence", "\n")}
  cat("Number of iterations=",iter,"\n")
  cat("The objective value=",tmp.obj,"\n")
  
  return(list(params=list(Z=Z, alpha=alpha, beta=beta), Theta=Theta.hat, obj=obj, iter=iter))                    
  
}



load('keyword_adj_20core.rda')

dim(year_keyword_adj_20core)
keywords <- rownames(year_keyword_adj_20core)
N <- dim(year_keyword_adj_20core)[1]
Ti <- dim(year_keyword_adj_20core)[3]

node_id <- seq(dim(year_keyword_adj_20core)[1])
rownames(year_keyword_adj_20core) <- node_id
colnames(year_keyword_adj_20core) <- node_id
rownames(isnew_keyword_adj_20core) <- node_id
colnames(isnew_keyword_adj_20core) <- node_id

node.exist <- list()
node.exist.num <- rep(NA, Ti)
for (t in 1:Ti) {
  na.num <- colSums(is.na(year_keyword_adj_20core[,,t]))
  node.exist[[t]] <- as.numeric(names(na.num[na.num<dim(year_keyword_adj_20core)[1]]))
  node.exist.num[t] <- length(node.exist[[t]])
}


A <- list()
X <- list()
for (t in 1:Ti) {
  A[[t]] <- year_keyword_adj_20core[,,t][node.exist[[t]], node.exist[[t]]]
  X[[t]] <- isnew_keyword_adj_20core[,,t][node.exist[[t]], node.exist[[t]]]
}


new.prop <- c()
for (t in 1:Ti) {
  new.prop[t] <- sum(colSums(X[[t]])==dim(X[[t]])[1]-1)/dim(X[[t]])[1]
}


para.ini <- init.func(N=N, D=2, Ti=Ti, node.exist=node.exist, seed=66) 
para.ini$beta
result <- opt.func(para.ini$Z, para.ini$alpha, sort(para.ini$beta), A=A, X=X, node.exist, nu=0.6, step.size=5, epsi=0.00001, niter=500, trace=1)

result$params$alpha <- array(result$params$alpha, dim=c(N,1))
rownames(result$params$alpha) <- keywords

for (t in 1:Ti) {
  rownames(result$params$Z[[t]]) <- keywords[node.exist[[t]]]
}
head(result$params$Z)
result$params$alpha
result$params$beta
result$obj[length(result$obj)]


