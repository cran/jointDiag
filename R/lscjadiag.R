##
lscjadiag <- function(C, B=NULL, eps=.Machine$double.eps, 
		itermax=200, keepTrace=FALSE) {
	# constants
	K <- dim(C)[1]
	n <- dim(C)[3]
	if (is.null(B)) B <- diag(K)
	# initialisation
	Ms <- matrix(0,K^2,K)
	L <- matrix(0,K,K)
	for (k in 1:K) {
		T <- sapply(1:n, function(naux) C[,,naux] %*% B[,k])
		Ms[,k] = matrix(T %*% t(T),K^2,1)
		L[,k] = matrix(Ms[,k],K,K) %*% B[,k]
	}
	s = sqrt(colSums(B*L))
	sqrts = sqrt(s)
	B = B/matrix(sqrts,K,K,byrow=T)	
	ssqrts = s*sqrts
	L = L/matrix(ssqrts,K,K,byrow=T)
 	M = matrix(rowSums(Ms/matrix(s,K^2,K,byrow=T)),K,K)
	gamma = colSums(B*(M %*% B))
	crit = sum(gamma) - K
	B_trace <- t(B)
	lcrit <- crit
	iter <- 0
	improve <- Inf
	# iteration
	while (improve>eps & iter<itermax) {
		R = chol(M)
		B0 = B
		aux <- solve(R)
		B = aux %*% (t(aux) %*% L)
		gamma = 1 - colSums(B*L)*gamma
		for (k in 1:K) {
			T <- sapply(1:n, function(naux) C[,,naux] %*% B[,k])
			Ms[,k] = matrix(T %*% t(T),K^2,1)
			L[,k] = matrix(Ms[,k],K,K) %*% B[,k]
		}
		s = sqrt(colSums(B*L))
		sqrts = sqrt(s)
		B = B/matrix(sqrts,K,K,byrow=T)	
		ssqrts = s*sqrts
		L = L/matrix(ssqrts,K,K,byrow=T)
		M = matrix(rowSums(Ms/matrix(s,K^2,K,byrow=T)),K,K)
		gamma = colSums(B*(M %*% B))
		crit = sum(gamma) - K
		improve <- rev(lcrit)[1]-crit
		lcrit <- c(lcrit,crit)
		if(keepTrace) B_trace <- cbind(B_trace,t(B))
		iter <- iter+1
	}
	if (iter == itermax) {
		warning("Convergence not reached")
	}
	if (keepTrace) {
		return(list(B=t(B),criter=lcrit,
			B_trace=array(B_trace,dim=c(K,K,iter+1))))
	}
	else {
		return(list(B=t(B),criter=lcrit,B_trace=NULL))
	}
}


















