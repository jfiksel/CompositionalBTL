sample.B <- function(M, p, v) {
    C <- length(p)
    B <- matrix(NA, ncol = C, nrow = C)
    for(j in 1:C){
        prob <- M[,j] * p
        if(sum(prob) == 0) {
            prob <- rep(.0001, length(p))
        }
        B[,j] <- rmultinom(n = 1, size = v[j], prob = prob)
    }
    return(B)
}

sample.M <- function(B, T.mat = NULL, power, gamma = 1, epsilon = .01) {
    C <- nrow(B)
    M <- matrix(NA, ncol = C, nrow = C)
    for(j in 1:nrow(M)){
        alpha <- (B[j,]  + T.mat[j,]) * power + gamma * epsilon 
        alpha[j] <- alpha[j] + gamma
        M[j,] <- rdirichlet(1, alpha)
    }
    return(M)
}

sample.M2 <- function(B, T.mat, C, gamma.vec, epsilon = .01, power) {
    M <- matrix(NA, ncol = C, nrow = C)
    for(i in 1:nrow(M)) {
        alpha <- (B[i,]  + T.mat[i,]) * power + gamma.vec * epsilon 
        alpha[i] <- alpha[i] + gamma.vec[i]
        M[i,] <- rdirichlet(1, alpha) 
    }
    return(M)
}



sample.p <- function(B, power, delta, C){
    alpha <- rep(NA, C)
    alpha <- (rowSums(B) * power) + delta
    p.samp <- rdirichlet(1, alpha)
    return(as.vector(p.samp))
}

sample.p.ensemble.power <- function(B.array, power, delta){
    alpha <- apply(B.array,1,sum) * power + delta
    p.samp <- rdirichlet(1, alpha)
    return(as.vector(p.samp))
}

initialize.M <- function(C) {
    M <- matrix(NA, nrow = C, ncol = C)
    for(i in 1:nrow(M)){
        M[i,] <- rdirichlet(1, rep(1, C))
    }
    return(M)
}

log.pi2 <- function(gamma.val, e, alpha, beta, M, M.row, tol=.00001) {
    C <- ncol(M)
    gammafn.q1 <- C * gamma.val * e + gamma.val + tol
    gammafn.q2 <- gamma.val * e + tol
    gammafn.q3 <- gamma.val * e + gamma.val + tol
    first.term <-
        lgamma(gammafn.q1) -
        (C - 1) * lgamma(gammafn.q2) -
        lgamma(gammafn.q3) 
    second.term <- (alpha - 1) * log(gamma.val + tol) - beta * gamma.val
    M.vec <- M[M.row,]
    M.gamma <- M.vec[M.row]
    M.notgamma <- M.vec[-M.row]
    third.term <-
        (gamma.val * e + gamma.val) * log(M.gamma + tol) +
        sum(gamma.val * e * log(M.notgamma + tol))
    return(first.term + second.term + third.term)
}

sample.gamma2 <- function(gamma.vec, epsilon, alpha, beta, M, tau.vec, max.gamma) {
    ### will use log scale for numeric stability
    # gamma.prop <- rgamma(1, shape = gamma.t * tau, rate = tau)
    # g.num <- dgamma(gamma.t, shape = gamma.prop * tau, rate = tau, log = TRUE)
    # g.denom <- dgamma(gamma.prop, shape = gamma.t * tau, rate = tau, log = TRUE)
    sampled.gammas <- sapply(seq_along(gamma.vec), function(i) {
        gamma.t <- gamma.vec[i]
        tau <- tau.vec[i]
        gamma.prop <- rlnorm(1, meanlog = log(gamma.t), sdlog = tau)
        while(gamma.prop >= max.gamma) {
            gamma.prop <- rlnorm(1, meanlog = log(gamma.t), sdlog = tau)
        }
        g.num <- dlnorm(gamma.t, meanlog = log(gamma.prop), sdlog = tau, log = TRUE)
        g.denom <- dlnorm(gamma.prop, meanlog = log(gamma.t), sdlog = tau, log = TRUE)
        pi.num <- log.pi2(gamma.prop, epsilon, alpha, beta, M, M.row = i)
        pi.denom <- log.pi2(gamma.t, epsilon, alpha, beta, M, M.row = i)
        log.alpha <- g.num + pi.num - g.denom - pi.denom
        accept.prob <- min(1, exp(log.alpha))
        u <- runif(1)
        if(u <= accept.prob) {
            return(gamma.prop)
        } else {
            return(gamma.t)
        }
    })
    return(sampled.gammas)
}

round_preserve_sum <- function(x, digits = 2) {
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    return(y / up)
}


compositional_btl <- function(A_U, A_L = NULL, G_L = NULL, causes, ndraws = 10000,
                                 burnin = 1000, thin = 5,
                                 delta = 1, power = 1,
                                 epsilon = .001,
                                 alpha = 5, beta = .5,
                                 tau = .5, max.gamma = 1E3,
                                 init.seeds = c(123,1234,12345)) {
    C <- length(causes)
    tau.vec <- rep(tau, C)
    ### Make pseudo data
    A_U <- t(apply(A_U, 1, function(x) round_preserve_sum(x, 2)))
    A_U_int <- A_U * 100
    d_U <- matrix(NA, nrow = nrow(A_U), ncol = 100)
    for(i in 1:nrow(d_U)){
        rep_factor <- as.integer(A_U_int[i,])
        rep_factor[rep_factor < 0] <- 0
        rep_factor[C] <- 100 - sum(rep_factor[1:(C-1)])
        d_U[i,] <- rep(1:C, rep_factor)
    }
    
    if(!is.null(A_L)) {
        A_L <- t(apply(A_L, 1, function(x) round_preserve_sum(x, 2)))
        A_L_int <- A_L * 100
        d_L <- matrix(NA, nrow = nrow(A_L), ncol = 100)
        for(i in 1:nrow(d_L)){
            rep_factor <- as.integer(A_L_int[i,])
            rep_factor[rep_factor < 0] <- 0
            rep_factor[C] <- 100 - sum(rep_factor[1:(C-1)])
            d_L[i,] <- rep(1:C, rep_factor)
        }
        
        ### Make pseudo latent variables for labeled set
        G_L <- t(apply(G_L, 1, function(x) round_preserve_sum(x, 2)))
        G_L_int <- G_L * 100
        z_L <- matrix(NA, nrow = nrow(G_L), ncol = 100)
        for(i in 1:nrow(z_L)){
            rep_factor <- as.integer(G_L_int[i,])
            rep_factor[rep_factor < 0] <- 0
            rep_factor[C] <- 100 - sum(rep_factor[1:(C-1)])
            z_L[i,] <- rep(1:C, rep_factor)
        } 
    }
    
    
    v <- sapply(1:C, function(c) sum(d_U==c))
    if(!is.null(G_L)) {
        T.mat <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:nrow(T.mat)) {
            for(j in 1:ncol(T.mat)) {
                T.mat[i,j] <- sum(z_L == i & d_L == j)
            }
        }  
    } else {
        T.mat <- matrix(0, nrow = C, ncol = C)
    }
    
    
    posterior.list <- lapply(seq_along(init.seeds), function(chain){
        seed <- init.seeds[chain]
        post.samples <- vector("list", ndraws)
        post.samples <- vector("list", ndraws)
        post.samples[[1]]$M <- initialize.M(C)
        post.samples[[1]]$p <- as.vector(rdirichlet(1, rep(1, C)))
        names(post.samples[[1]]$p) <- causes
        post.samples[[1]]$B <- sample.B(post.samples[[1]]$M, post.samples[[1]]$p, v)
        
        gamma.init <- rgamma(C, alpha, beta)
        post.samples[[1]]$gamma <- sample.gamma2(gamma.init, epsilon,
                                                 alpha, beta, post.samples[[1]]$M,
                                                 tau.vec, max.gamma)
        for(i in 2:ndraws){
            post.samples[[i]]$M <- sample.M2(B = post.samples[[i-1]]$B,
                                             T.mat = T.mat,
                                             C = C,
                                             power = power,
                                             gamma.vec = post.samples[[i-1]]$gamma,
                                             epsilon = epsilon)
            post.samples[[i]]$p <- sample.p(post.samples[[i-1]]$B, power = power, delta = delta, C = C)
            #names(post.samples[[i]]$p) <- causes
            post.samples[[i]]$B <- sample.B(post.samples[[i]]$M, post.samples[[i]]$p, v= v)
            post.samples[[i]]$gamma <- sample.gamma2(post.samples[[i-1]]$gamma,
                                                     epsilon, alpha, beta,
                                                     post.samples[[i]]$M, tau.vec,
                                                     max.gamma)
            if((i%%10000)==0) message(paste("Chain", chain, "Draw", i))
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
    post.samples.list  <- mcmc.list(posterior.list)
    post.samples.list <- window(post.samples.list, start = burnin, thin = thin)
    return(post.samples.list)
}


compositional_ensemble_btl <- function(A_U, A_L, G_L = NULL, causes, ndraws = 10000,
                                       burnin = 1000, thin = 5,
                                       delta = 1, power = 1,
                                       epsilon = .001,
                                       alpha = 5, beta = .5,
                                       tau = .5, max.gamma = 1E3,
                                       init.seeds = c(123,1234,12345)) {
    
    C <- length(causes)
    tau.vec <- rep(tau, C)
    ### Make pseudo data
    K <- dim(A_U)[3]
    d_U <- array(NA, dim= c(nrow(A_U), 100, K))
    d_L <- array(NA, dim= c(nrow(A_L), 100, K))
    for(k in 1:K) {
        A_U[,,k] <- t(apply(A_U[,,k], 1, function(x) round_preserve_sum(x, 2)))
        A_U_int <- A_U[,,k] * 100
        
        for(i in 1:nrow(d_U[,,k])){
            rep_factor <- as.integer(A_U_int[i,])
            rep_factor[rep_factor < 0] <- 0
            rep_factor[C] <- 100 - sum(rep_factor[1:(C-1)])
            d_U[i,,k] <- rep(1:C, rep_factor)
        }
        
        A_L[,,k] <- t(apply(A_L[,,k], 1, function(x) round_preserve_sum(x, 2)))
        A_L_int <- A_L[,,k] * 100
        
        for(i in 1:nrow(d_L[,,k])){
            rep_factor <- as.integer(A_L_int[i,])
            rep_factor[rep_factor < 0] <- 0
            rep_factor[C] <- 100 - sum(rep_factor[1:(C-1)])
            d_L[i,,k] <- rep(1:C, rep_factor)
        }
    }
    
    
    ### Make pseudo latent variables for labeled set
    G_L_int <- t(apply(G_L, 1, function(x) round_preserve_sum(x, 2)))
    G_L_int <- G_L_int * 100
    z_L <- matrix(NA, nrow = nrow(G_L), ncol = 100)
    for(i in 1:nrow(z_L)){
        rep_factor <- as.integer(G_L_int[i,])
        rep_factor[rep_factor < 0] <- 0
        rep_factor[C] <- 100 - sum(rep_factor[1:(C-1)])
        z_L[i,] <- rep(1:C, rep_factor)
    }
    
    v <- matrix(NA, nrow = C, ncol = K)
    for(k in 1:K) {
        v[,k] <- sapply(1:C, function(c) sum(d_U[,,k]==c))
    }
    T.array <- array(NA, dim = c(C, C, K))
    for(k in 1:K) {
        for(i in 1:C) {
            for(j in 1:C) {
                T.array[i,j,k] <- sum(z_L == i & d_L[,,k] == j)
            }
        }  
    }
    
    posterior.list <- lapply(seq_along(init.seeds), function(chain) {
        seed <- init.seeds[chain]
        set.seed(seed)
        ### Initialize array of M matrices
        M.array <- sapply(1:K, function(k) {
            T.mat <- T.array[,,k]
            M.mat <- initialize.M(C)
        }, simplify = 'array')
        
        post.samples <- vector("list", ndraws)
        
        post.samples[[1]]$M.array <- M.array
        # post.samples[[1]]$p <- rep(1 / length(causes), length(causes))
        init.p.list <- lapply(1:ncol(v), function(k) {
            as.vector(rdirichlet(1, rep(1, C)))
        })
        init.p <- Reduce("+", init.p.list) / length(init.p.list)
        init.p <- init.p / sum(init.p)
        post.samples[[1]]$p <- init.p
        # post.samples[[1]]$p <- sapply(causes, function(c) mean(test.cod.mat[,1] == c))
        # names(post.samples[[1]]$p) <- causes
        
        B.array=sapply(1:K, function(k) {
            sample.B(post.samples[[1]]$M.array[,,k], post.samples[[1]]$p, v[,k])    
        }, simplify="array")
        post.samples[[1]]$B.array <- B.array
        gamma.init.mat <- matrix(NA, nrow = C, ncol = K)
        for(k in 1:K) {
            gamma.init.mat[,k] <- rgamma(C, alpha, beta) 
        }
        
        post.samples[[1]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
        for(k in 1:K) {
            T.mat <- T.array[,,k]
            M.mat <- post.samples[[1]]$M.array[,,k]
            post.samples[[1]]$gamma.mat[,k] <- sample.gamma2(gamma.init.mat[,k],
                                                             epsilon,
                                                             alpha, beta, M.mat, tau.vec,
                                                             max.gamma) 
        }
        
        for(i in 2:ndraws){
            post.samples[[i]]$M.array <- sapply(1:K, function(k) 
                sample.M2(post.samples[[i-1]]$B.array[,,k], 
                          T.mat = T.array[,,k],
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
            
            post.samples[[i]]$gamma.mat <- matrix(NA, nrow = C, ncol = K)
            for(k in 1:K) {
                M.mat <- post.samples[[i]]$M.array[,,k]
                post.samples[[i]]$gamma.mat[,k] <- sample.gamma2(post.samples[[i-1]]$gamma.mat[,k],
                                                                 epsilon, alpha, beta, M.mat,
                                                                 tau.vec, max.gamma) 
            }
            #post.samples[[i]]$gamma <- gamma.init
            #if((i%%1000)==0) print(paste("Run", i, post.samples[[i]]$gamma))
            #print(paste("Draw", i))
            if((i%%10000)==0) message(paste("Chain", chain, "Draw", i))
        }
        #return(post.samples)
        ### Put everything into a matrix, to be converted into an mcmc object
        ### Number of params is 2 * k * C ^ 2 (M matrix and B matrix for each algorithm)
        ###                     + C (p vector) + K * C (gamma vector for each algorithm)
        n.params <- K * C^2 + C + K * C
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
    post.samples.list  <- mcmc.list(posterior.list)
    post.samples.list <- window(post.samples.list, start = burnin, thin = thin)
    return(post.samples.list)
}



getP <- function(post.samples) {
    p.list <- lapply(seq_along(post.samples), function(i){
        x <- post.samples[[i]]
        p <- x$p
        causes <- names(p)
        p <- unname(p)
        return(data.frame(p = p, cause = causes, sample = i))
    })
    p.df <- do.call(rbind, p.list)
    return(p.df)
}

getM <- function(post.samples) {
    M.list <- lapply(post.samples, function(x) x$M)
    M.mean <- Reduce("+", M.list) / length(M.list)
    return(M.mean)
}

getRhat <- function(post.samples) {
    rhat <- gelman.diag(post.samples, multivariate = F, autoburnin = F)$psrf
    rhat.df <- tibble(Parameter = rownames(rhat), Rhat = unname(rhat[,1]))
    return(rhat.df)
}
