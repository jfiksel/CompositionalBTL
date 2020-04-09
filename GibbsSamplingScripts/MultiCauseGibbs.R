### Function to sample T for M for r in L
sample.T <- function(M, C, G_L, A_L_int) {
    B <- array(NA, dim = c(C, C, nrow(A_L_int)))
    for(j in 1:C) {
        probs <- sweep(G_L, MARGIN=2, M[,j], `*`)
        for(r in 1:nrow(A_L_int)) {
            prob <- probs[r,] 
            if(sum(prob) == 0) {
                prob <- rep(.0001, C)
            }
            B[,j,r] <- rmultinom(n = 1, size = A_L_int[r, j], prob = prob) 
        } 
    }
    T.mat <- apply(B, 1:2, sum)
    return(T.mat)
}

compositional_multi_btl <- function(A_U, A_L, G_L = NULL, causes, ndraws = 10000,
                              burnin = 1000, thin = 5,
                              delta = 1, gamma = 1, power = 1,
                              epsilon = .001,
                              alpha = 5, beta = .5,
                              tau = .5, max.gamma = 1E3,
                              init.seeds = c(123,1234,12345)) {
    C <- length(causes)
    tau.vec <- rep(tau, C)
    ### Make pseudo data
    A_U <- t(apply(A_U, 1, function(x) round_preserve_sum(x, 2)))
    A_U_int <- A_U * 100
    
    A_L <- t(apply(A_L, 1, function(x) round_preserve_sum(x, 2)))
    A_L_int <- A_L * 100
    
    v <- colSums(A_U_int)
    posterior.list <- lapply(seq_along(init.seeds), function(chain){
        seed <- init.seeds[chain]
        post.samples <- vector("list", ndraws)
        post.samples <- vector("list", ndraws)
        post.samples[[1]]$M <- initialize.M(C)
        post.samples[[1]]$p <- as.vector(rdirichlet(1, rep(1, C)))
        names(post.samples[[1]]$p) <- causes
        ### sample B for U
        post.samples[[1]]$B <- sample.B(post.samples[[1]]$M, post.samples[[1]]$p, v)
        ### Sample B for L
        post.samples[[1]]$T.mat <- sample.T(post.samples[[1]]$M, C, G_L, A_L_int)
        
        gamma.init <- rgamma(C, alpha, beta)
        post.samples[[1]]$gamma <- sample.gamma2(gamma.init, epsilon,
                                                 alpha, beta, post.samples[[1]]$M,
                                                 tau.vec, max.gamma)
        for(i in 2:ndraws){
            post.samples[[i]]$M <- sample.M2(B = post.samples[[i-1]]$B,
                                             T.mat = post.samples[[i-1]]$T.mat,
                                             C = C,
                                             power = power,
                                             gamma.vec = post.samples[[i-1]]$gamma,
                                             epsilon = epsilon)
            post.samples[[i]]$p <- sample.p(post.samples[[i-1]]$B, power = power, delta = delta, C = C)
            #names(post.samples[[i]]$p) <- causes
            post.samples[[i]]$B <- sample.B(post.samples[[i]]$M, post.samples[[i]]$p, v= v)
            post.samples[[i]]$T.mat <- sample.T(post.samples[[i]]$M, C, G_L, A_L_int)
            post.samples[[i]]$gamma <- sample.gamma2(post.samples[[i-1]]$gamma,
                                                     epsilon, alpha, beta,
                                                     post.samples[[i]]$M, tau.vec,
                                                     max.gamma)
            if((i%%5000)==0) message(paste("Chain", chain, "Draw", i))
        }
        # return.index <- seq(burnin, ndraws, by = thin)
        ### Put everything into a matrix, to be converted into an mcmc object
        ### Number of params is 2 * C ^ 2 (M matrix and B matrix) + 2 * p (gamma vector & p vector)
        n.params <- C^2 + C * 2
        post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
        for(i in 1:nrow(post.samples.mat)){
            samp <- post.samples[[i]]
            p.vec <- unname(samp$p)
            M.vec <- as.vector(samp$M)
            gamma.vec <- samp$gamma
            post.samples.mat[i,] <- c(p.vec, M.vec, gamma.vec)
        }
        ### Column names is first cause names (with prefix p)
        ### then M (as.vector goes by column)
        cnames <- c(paste0("p[", 1:C, "]"),
                    paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]"),
                    paste0("gamma[", 1:C, "]"))
        
        colnames(post.samples.mat) <- cnames
        return(mcmc(post.samples.mat))
    })
    posterior.list <- mcmc.list(posterior.list)
    posterior.list <- window(posterior.list, start = burnin, thin = thin)
    return(posterior.list)
}

compositional_multi_ensemble_btl <- function(A_U, A_L, G_L = NULL, causes, ndraws = 10000,
                                             burnin = 1000, thin = 5,
                                             delta = 1, gamma = 1, power = 1,
                                             epsilon = .001,
                                             alpha = 5, beta = .5,
                                             tau = .5, max.gamma = 1E3,
                                             init.seeds = c(123,1234,12345)) {
    C <- length(causes)
    tau.vec <- rep(tau, C)
    ### Make pseudo data
    ### Make pseudo data
    K <- dim(A_U)[3]
    for(k in 1:K) {
        A_U[,,k] <- t(apply(A_U[,,k], 1, function(x) round_preserve_sum(x, 2)))
        A_L[,,k] <- t(apply(A_L[,,k], 1, function(x) round_preserve_sum(x, 2)))
    }
    A_U_int <- A_U * 100
    A_L_int <- A_L * 100
    v <- matrix(NA, nrow = C, ncol = K)
    for(k in 1:K) {
        v[,k] <- colSums(A_U_int[,,k])
    }
    
    posterior.list <- lapply(seq_along(init.seeds), function(chain){
        seed <- init.seeds[chain]
        post.samples <- vector("list", ndraws)
        post.samples <- vector("list", ndraws)
        M.array <- sapply(1:K, function(k) {
            M.mat <- initialize.M(C)
        }, simplify = 'array')
        
        post.samples[[1]]$M.array <- M.array
        init.p.list <- lapply(1:ncol(v), function(k) {
            as.vector(rdirichlet(1, rep(1, C)))
        })
        init.p <- Reduce("+", init.p.list) / length(init.p.list)
        init.p <- init.p / sum(init.p)
        post.samples[[1]]$p <- init.p
        #names(post.samples[[1]]$p) <- causes
        ### sample B for U
        B.array <- sapply(1:K, function(k) {
            sample.B(post.samples[[1]]$M.array[,,k], post.samples[[1]]$p, v[,k])    
        }, simplify="array")
        post.samples[[1]]$B.array <- B.array
        ### Sample B for L
        T.array <- sapply(1:K, function(k) {
            sample.T(post.samples[[1]]$M.array[,,k], C, G_L, A_L_int[,,k])  
        }, simplify = "array")
        post.samples[[1]]$T.array <- T.array
        
        gamma.init.mat <- matrix(NA, nrow = C, ncol = K)
        for(k in 1:K) {
            gamma.init.mat[,k] <- rgamma(C, alpha, beta) 
        }
        
        post.samples[[1]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
        for(k in 1:K) {
            M.mat <- post.samples[[1]]$M.array[,,k]
            post.samples[[1]]$gamma.mat[,k] <- sample.gamma2(gamma.init.mat[,k],
                                                             epsilon,
                                                             alpha, beta, M.mat, tau.vec,
                                                             max.gamma) 
        }
        for(i in 2:ndraws){
            post.samples[[i]]$M.array <- sapply(1:K, function(k) 
                sample.M2(post.samples[[i-1]]$B.array[,,k], 
                          T.mat = post.samples[[i-1]]$T.array[,,k],
                          C = C,
                          power = power,
                          gamma.vec = post.samples[[i-1]]$gamma.mat[,k],
                          epsilon = epsilon), simplify="array")
            post.samples[[i]]$p <- sample.p.ensemble.power(post.samples[[i-1]]$B,
                                                           power,
                                                           delta)
            #names(post.samples[[i]]$p) <- causes
            post.samples[[i]]$B.array <- sapply(1:K, function(k) {
                sample.B(post.samples[[i]]$M.array[,,k], post.samples[[i]]$p, v[,k])    
            }, simplify="array")
            post.samples[[i]]$T.array <- sapply(1:K, function(k) {
                sample.T(post.samples[[i]]$M.array[,,k], C, G_L, A_L_int[,,k])    
            }, simplify="array")
            post.samples[[i]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
            for(k in 1:K) {
                M.mat <- post.samples[[i]]$M.array[,,k]
                post.samples[[i]]$gamma.mat[,k] <- sample.gamma2(post.samples[[i-1]]$gamma.mat[,k],
                                                                 epsilon, alpha, beta, M.mat,
                                                                 tau.vec, max.gamma) 
            }
            if((i%%5000)==0) message(paste("Chain", chain, "Draw", i))
        }
        # return.index <- seq(burnin, ndraws, by = thin)
        ### Put everything into a matrix, to be converted into an mcmc object
        ### Number of params is 2 * C ^ 2 (M matrix and B matrix) + 2 * p (gamma vector & p vector)
        n.params <-  K * C^2 + C + K * C
        post.samples.mat <- matrix(nrow = length(post.samples), ncol = n.params)
        for(i in 1:nrow(post.samples.mat)){
            samp <- post.samples[[i]]
            p.vec <- samp$p
            M.vec <- unlist(lapply(1:K, function(k) as.vector(samp$M.array[,,k])))
            gamma.vec <- unlist(lapply(1:K, function(k) samp$gamma.mat[,k]))
            post.samples.mat[i,] <- c(p.vec, M.vec, gamma.vec)
        }
        ### Column names is first cause names (with prefix p)
        ### then M (as.vector goes by column)
        p.names <- paste0("p[", 1:C, "]")
        M.names <- unlist(lapply(1:K, function(k) {
            paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), ",", k, "]")
        }))
        gamma.names <- unlist(lapply(1:K, function(k) {
            paste0("gamma[", 1:C, ",", k, "]")
        }))
        cnames <- c(p.names,
                    M.names,
                    gamma.names)
        
        colnames(post.samples.mat) <- cnames
        return(mcmc(post.samples.mat))
    })
    posterior.list <- mcmc.list(posterior.list)
    posterior.list <- window(posterior.list, start = burnin, thin = thin)
    return(posterior.list)
}