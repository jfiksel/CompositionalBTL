#source("sampling_functions.R")

initialize.p <- function(v) {
    ### v is a vector with predicted counts for each COD
    #return(v / sum(v))
    #rep(1/length(v),length(v))
    ### First normalize v
    v.norm <- v / sum(v)
    ### Then add normalized version of v to 1
    alpha <- v.norm + 1
    as.vector(rdirichlet(1, alpha))
}

sample.B <- function(M, p, v) {
    C <- length(p)
    B <- matrix(NA, ncol = C, nrow = C)
    for(i in 1:C){
        prob <- M[,i] * p
        if(sum(prob) == 0) {
            prob <- rep(.0001, length(p))
        }
        B[,i] <- rmultinom(n = 1, size = v[i], prob = prob)
    }
    return(B)
}

sample.p.ensemble.power <- function(B.array, power, delta){
    alpha <- apply(B.array,1,sum) * power + delta
    p.samp <- rdirichlet(1, alpha)
    return(as.vector(p.samp))
}

round_preserve_sum <- function(x, digits = 2) {
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    return(y / up)
}

calibva.ensemble.power.sampler <- function(A_U, A_L, G_L = NULL, causes, ndraws = 10000,
                                           burnin = 1000, thin = 5,
                                           delta = 1, gamma = 1, power = 1,
                                           epsilon.M = .00001,
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
    G_L <- t(apply(G_L, 1, function(x) round_preserve_sum(x, 2)))
    G_L_int <- G_L * 100
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
        predictions.list <- lapply(1:length(post.samples), function(i) matrix(NA, dim(A_U)[1], ncol = C))
        for(i in 1:nrow(post.samples.mat)){
            samp <- post.samples[[i]]
            ### get predictions first
            B.array <- samp$B.array
            predictions <- lapply(1:K, function(k) {
                z_probs <- apply(B.array[,,k], 2, function(x) x / sum(x))
                A_U[,,k] %*% t(z_probs)
            })
            predictions.list[[i]] <- Reduce("+", predictions) / K
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
        colnames(post.samples.mat) <- cnames
        return(list(post.samples = mcmc(post.samples.mat), predictions = predictions.list))
    })
    post.samples.list <- lapply(1:length(posterior.list), function(i) posterior.list[[i]]$post.samples)
    samples_seq <- seq((burnin + 1), ndraws, by = thin)
    predictions.list <- lapply(1:length(posterior.list), function(i) {
        preds_list <- posterior.list[[i]]$predictions
        preds_list <- preds_list[samples_seq]
        preds_mean <- Reduce("+", preds_list) / length(preds_list)
        return(preds_mean)
    })
    predictions <-  Reduce("+", predictions.list) / length(predictions.list)
    predictions.df <- do.call(rbind, predictions.list)
    post.samples.list  <- mcmc.list(post.samples.list)
    post.samples.list <- window(post.samples.list, start = burnin, thin = thin)
    return(list(posterior = post.samples.list, predictions = predictions))
}
