
# cozigam(COZIGAM)

cozigam <- function(formula, constraint = "proportional", zero.delta = NULL,
    maxiter = 20, conv.crit.in = 1e-5, conv.crit.out = 1e-3, size = NULL,
    log.tran = FALSE, family, ...)
{

      if (is.character(family))
          fam <- eval(parse(text = family))
      if (is.function(family))
          fam <- family()

      if (constraint == "proportional") {
          if (fam$family == "gaussian" | fam$family == "Gamma") {
              cozigam.res <- PCOZIGAM.cts(formula, maxiter, conv.crit.in,
                  conv.crit.out, log.tran = log.tran, family = fam, ...)
              attr(cozigam.res, "family.type") <- "continuous"
          }
          else if (fam$family == "poisson" | fam$family == "binomial") {
              cozigam.res <- PCOZIGAM.dis(formula, maxiter, conv.crit.in,
                  conv.crit.out, size = size, family = fam, ...)
              attr(cozigam.res, "family.type") <- "discrete"
          }
          else stop("family not recognized")
          attr(cozigam.res, "constraint") <- "proportional"
      }

      else if (constraint == "component") {
          if (missing(zero.delta) | is.null(zero.delta))  stop("zero.delta not provided")
          if (fam$family == "gaussian" | fam$family == "Gamma") {
              cozigam.res <- COZIGAM.cts(formula, zero.delta, maxiter, conv.crit.in,
                  conv.crit.out, log.tran = log.tran, family = fam, ...)
              attr(cozigam.res, "family.type") <- "continuous"
          }
          else if (fam$family == "poisson" | fam$family == "binomial") {
              cozigam.res <- COZIGAM.dis(formula, zero.delta, maxiter, conv.crit.in,
                  conv.crit.out, size = size, family = fam, ...)
              attr(cozigam.res, "family.type") <- "discrete"
          }
          else stop("family not recognized")
          attr(cozigam.res, "constraint") <- "component"
      }
      else stop("constraint not recognized")

      invisible(cozigam.res)

}

# COZIGAM.cts(COZIGAM)
# Component-Specific-COnstrained Zero-Inflated GAM: Continuous Case

COZIGAM.cts <- function(formula, zero.delta, maxiter = 20, conv.crit.in = 1e-5,
    conv.crit.out = 1e-3, log.tran = FALSE, family = gaussian(), ...)
{
        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        n <- length(y)
        z <- (y!=0)
        if(log.tran==TRUE) y[z] <- log(y[z])

        if (is.character(family))
            family <- eval(parse(text = family))
        if (is.function(family))
            family <- family()
        if(family$family == "gaussian") {
            d.V <- function(mu) rep.int(0, length(mu))
            d.eta.mu <- function(mu) rep.int(0, length(mu))
            loglikfun <- function(y, mu, disp, z, p) {
                sum(z*(log(p)+dnorm(y,mu,sqrt(disp),log=TRUE))+(1-z)*log(1-p))
            }
        }
        else if(family$family == "Gamma") {
            d.V <- function(mu) 2*mu
            d.eta.mu <- function(mu) 2/(mu^3)
            loglikfun <- function(y, mu, disp, z, p) {
                sum(z*(log(p)+ifelse(y==0, 1, dgamma(y,1/disp,scale=mu*disp,log=TRUE)))+(1-z)*log(1-p))
            }
        }
        else
            stop("family not recognized")

        variance <- family$variance
        linkinv <- family$linkinv
        linkfun <- family$linkfun
        mu.eta <- family$mu.eta

        ## Estimate dispersion parameter using non-zero data...
        fm <- as.formula(sub(split$response,"y",deparse(formula)))
        G1 <- gam(fm, subset=z, family=family, fit=FALSE, ...)
        fit.nz <- gam(fm, subset=z, family=family, ...)
        disp <- fit.nz$sig2; est.disp <- TRUE

        ## Set initial values...
        mu <- pmax(y, 0.01)
        p <- rep(0.7, n)
        eta <- linkfun(mu) # g_mu(mu)
        gp <- eta # g_p(p)-alpha, equivalently alpha=0, delta=1

        ## Model setup
        G <- gam(fm, fit=FALSE, family=family, ...)

        np <- ncol(G$X) # length of parameter beta
        n.smooth <- length(G$smooth)
        if(n.smooth!=length(zero.delta))
            stop("length of 'zero.delta' does not match # of smooth terms")
        alpha <- 0
        zdelta <- !is.na(zero.delta)
        delta <- rep.int(1, n.smooth)
        delta[zdelta] <- 0 # fix zero delta's

        M <- list() # extract the design matrix of each smooth term
        for(k in 1:n.smooth) {
            M[[k]] <- G$X[, G$smooth[[k]]$first.para:G$smooth[[k]]$last.para]
        }
        X.par <- as.matrix(G$X[, 1:(G$smooth[[1]]$first.para-1)]) # design matrix of parametric part

        ## Begin outer loop: iteratively update beta and alpha & delta
        rep.out <- 1; norm.out <- 1
        while(norm.out>conv.crit.out & rep.out<maxiter) {

            b.old <- b <- rep.int(0, np)

            Xp <- matrix(0, ncol=ncol(X.par), nrow=nrow(X.par)) # the design matrix of the second SS
            for(k in 1:n.smooth) {
               Xp <- cbind(Xp, delta[k]*M[[k]])
            }

            ## inner loop: approx logLik by WSS
            ## update beta by P-IRLS
            rep.in <- 1; norm.in <- 1
            while(norm.in>conv.crit.in & rep.in<maxiter) {

                mu.eta.val <- mu.eta(eta)
                weights1 <- sqrt(z*(mu.eta.val^2)/variance(mu)/disp)
                weights2 <- sqrt(p*(1-p))

                good <- (weights1 > 0) & (weights2 > 0) & (mu<1000) & (mu.eta.val != 0)
                good[is.na(good)] <- FALSE
                if (all(!good)) {
                    converged <- FALSE
                    warning("No observations informative at iteration")
                    break
                }
                z1 <- eta[good] + (y - mu)[good]/mu.eta.val[good]
                z2 <- gp[good] + 1/(p*(1-p))[good]*(z-p)[good]

                y.c <- c(z1, z2)
                w.c <- c(weights1[good], weights2[good])
                X.c <- rbind(G$X[good,], Xp[good,])

                mgfit <- magic(y=y.c, X=X.c, sp=G$sp, S=G1$S, off=G$off, rank=G$rank, C=G$C, w=w.c, gcv=FALSE)

                b.old <- b
                b <- mgfit$b
                sp <- mgfit$sp
                eta <- G$X%*%b

                Eta <- NULL # the design matrix of binomial model
                bp <- list() # extract the paramters for each smooth term
                eta.p <- list() # record the estimate the each smooth function
                gp <- 0
                for(k in 1:n.smooth) {
                    bp[[k]] <- b[G$smooth[[k]]$first.para:G$smooth[[k]]$last.para]
                    eta.p[[k]] <- M[[k]]%*%bp[[k]]
                    gp <- gp + delta[k]*eta.p[[k]]
                    if(!zdelta[k]) Eta <- cbind(Eta, eta.p[[k]])
                }

                mu <- linkinv(eta)
                p <- .Call("logit_linkinv", alpha + gp, PACKAGE = "stats")
                norm.in <- sum((b-b.old)^2)
                rep.in <- rep.in + 1
            }

            ## Update alpha & delta
            ## Fit a generalized linear model given beta
            alpha.old <- alpha
            delta.old <- delta

            if(is.null(Eta)) alpha <- coef(glm(z ~ 1, family=binomial)) # Homogeneous ZI
            else {
                fit.glm <- glm(z ~ Eta, family=binomial) # z ~ eta.p[[1]] +...+ eta.p[[k]]
                alpha <- fit.glm$coef[1]
                delta[!zdelta] <- fit.glm$coef[-1]
            }

            ## Update p for new alpha & delta
            gp <- 0
            for(k in 1:n.smooth) {
                gp <- gp + delta[k]*eta.p[[k]]
            }
            p <- .Call("logit_linkinv", alpha + gp, PACKAGE = "stats")

            if(is.null(Eta)) norm.out <- abs(alpha-alpha.old) # actually only one outer iteration
            else norm.out <- max(abs(alpha-alpha.old), abs(delta-delta.old))
            rep.out <- rep.out + 1
            cat("iteration =", rep.out, "\t", "norm =", norm.out, "\n")

        }

        names(alpha) <- names(delta) <- NULL
        theta <- list(alpha=alpha, delta=delta, beta=b)
        if(rep.out==maxiter) converged <- FALSE
            else converged <- TRUE

        ## Compute observed information...
        Lambda <- matrix(0, np, np) # penalty matrix

        rho <- z-p; d.rho <- rep.int(-1, n)
        tau <- z*(y-mu)*mu.eta.val/disp/variance(mu)
        d.tau <- z*mu.eta.val^2*( -variance(mu)/mu.eta.val-(y-mu)*(variance(mu)*d.eta.mu(mu)+
            d.V(mu)/mu.eta.val) )/disp/(variance(mu)^2)


        Lam <- list(); n.S <- numeric(n.smooth);
        Xp <- matrix(0, ncol=ncol(X.par), nrow=nrow(X.par))
        P <- NULL
        for(k in 1:n.smooth) {
            n.S[k] <- length(G1$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp[k]*G1$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[j]*G1$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp[sum(n.S[1:(k-1)])+1]*G1$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[sum(n.S[1:(k-1)])+j]*G1$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G1$smooth[[k]]$first.para
            last <- G1$smooth[[k]]$last.para
            Lambda[first:last, first:last] <- Lam[[k]]
            Xp <- cbind(Xp, delta[k]*M[[k]])
            if(!zdelta[k]) {
                Xrho <- t(G$X)%*%rho; Xrho[-(G1$smooth[[k]]$first:G1$smooth[[k]]$last)] <- 0
                P <- cbind(P, Xrho) # used in computing partial derivatives
            }
        }

        ## Find the effective degrees of freedom
        FM <- solve(t(X.c)%*%diag(w.c^2)%*%X.c+Lambda)%*%t(X.c)%*%diag(w.c^2)%*%X.c

        edf.smooth <- numeric(n.smooth)
        edf <- diag(FM)
        for(k in 1:n.smooth) {
            first <- G$smooth[[k]]$first.para
            last <- G$smooth[[k]]$last.para
            edf.smooth[k] <- sum(edf[first:last])
        }

        if (G$nsdf > 0)
            term.names <- colnames(G$X)[1:G$nsdf]
        else term.names <- array("", 0)
        for (i in 1:n.smooth) {
            k <- 1
            for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) {
                term.names[j] <- paste(G$smooth[[i]]$label, ".",
                  as.character(k), sep = "")
                k <- k + 1
            }
        }
        names(theta$beta) <- term.names
        names(edf) <- term.names

        G.rho <- diag(as.vector(d.rho*p*(1-p)))
        G.tau <- diag(as.vector(d.tau*mu.eta.val))

        if(is.null(Eta)) {
            I.alpha <- sum(-d.rho*p*(1-p))
            I.beta <- -t(G$X)%*%G.tau%*%G$X - t(Xp)%*%G.rho%*%Xp + Lambda
            I.alpha.beta <- t(Xp)%*%(-d.rho*p*(1-p))

            n.theta <- 1 + np
            I.theta <- matrix(0, ncol=n.theta, nrow=n.theta)
            I.theta[1,1] <- I.alpha
            I.theta[1,2:n.theta] <- I.theta[2:n.theta,1] <- I.alpha.beta
            I.theta[2:n.theta,2:n.theta] <- I.beta
            colnames(I.theta) <- rownames(I.theta) <- NULL
            V.theta <- solve(I.theta)
            se.alpha <- sqrt(V.theta[1,1])
            se.delta <- rep.int(0, n.smooth)
            V.b <- V.theta[2:n.theta,2:n.theta]
        }

        else {
            I.alpha <- sum(-d.rho*p*(1-p))
            I.delta <- -t(Eta)%*%G.rho%*%Eta
            I.alpha.delta <- t(Eta)%*%(-d.rho*p*(1-p))
            I.beta <- -t(G$X)%*%G.tau%*%G$X - t(Xp)%*%G.rho%*%Xp + Lambda
            I.alpha.beta <- t(Xp)%*%(-d.rho*p*(1-p))
            I.delta.beta <- -t(Xp)%*%G.rho%*%Eta - P

            n.delta <- sum(!zdelta) # number of non-zero delta
            n.theta <- 1 + n.delta + np
            I.theta <- matrix(0, ncol=n.theta, nrow=n.theta)
            I.theta[1,1] <- I.alpha
            I.theta[1,2:(n.delta+1)] <- I.alpha.delta
            I.theta[1,(n.delta+2):n.theta] <- I.alpha.beta
            I.theta[2:(n.delta+1), 1] <- I.alpha.delta
            I.theta[2:(n.delta+1),2:(n.delta+1)] <- I.delta
            I.theta[2:(n.delta+1),(n.delta+2):n.theta] <- t(I.delta.beta)
            I.theta[(n.delta+2):n.theta,1] <- I.alpha.beta
            I.theta[(n.delta+2):n.theta,2:(n.delta+1)] <- I.delta.beta
            I.theta[(n.delta+2):n.theta,(n.delta+2):n.theta] <- I.beta
            colnames(I.theta) <- rownames(I.theta) <- NULL
            V.theta <- solve(I.theta)
            se.alpha <- sqrt(V.theta[1,1])
            se.delta <- rep.int(0, n.smooth)
            se.delta[!zdelta] <- sqrt(diag(V.theta)[2:(n.delta+1)])
            V.b <- V.theta[(n.delta+2):n.theta,(n.delta+2):n.theta]
        }

        loglik <- loglikfun(y, mu, disp, z, p)  # log-likelihood
        ploglik <- loglik - as.numeric(0.5*t(b)%*%Lambda%*%b)  # penalized log-likelihood
        DS <- diag(eigen(Lambda)$values[abs(eigen(Lambda)$values)>1e-10])
        # Laplace approximation of log marginal likelihood
        logE <- 0.5*determinant(DS)$modulus+ploglik+(n.theta-nrow(DS))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
        attr(logE, "logarithm") <- NULL

        cat("\n")
        cat("==========================================","\n")
        cat("estimated alpha =", alpha, "(", se.alpha, ")", "\n")
        for (k in 1:n.smooth) {
            cat(paste("estimated delta",k, sep=""),"=", delta[k], "(", se.delta[k], ")", "\n")
        }
        cat("==========================================","\n","\n")

        res <- list(coefficient=theta, V.theta=V.theta, converged=converged,
              mu=mu, linear.predictor=eta, dispersion=disp, formula=formula,
              fit.nonzero=fit.nz, score=mgfit$score, G=G, V.beta=V.b,
              y=y, p=p, loglik=loglik, ploglik=ploglik, family=family, edf=edf,
              est.disp=est.disp, edf.smooth=edf.smooth, sp=sp, se.delta=se.delta,
              se.alpha=se.alpha, logE=logE)
        class(res) <- "cozigam"
        res

}

# COZIGAM.dis(COZIGAM)
# Component-Specific-COnstrained Zero-Inflated GAM: Distrete Case

COZIGAM.dis <- function(formula, zero.delta, maxiter = 20, conv.crit.in = 1e-5,
    conv.crit.out = 1e-3, size = NULL, family = poisson(), ...)
{
        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        n <- length(y)

        if (is.character(family))
            family <- eval(parse(text = family))
        if (is.function(family))
            family <- family()
        if(family$family == "poisson") {
            d.V <- function(mu) rep.int(1, length(mu))
            d.eta.mu <- function(mu) -1/(mu^2)
            den <- dpois; disp <- 1; est.disp <- FALSE
            size <- rep.int(1, n)
            loglikfun <- function(y, mu, p) {
                e <- as.numeric(y!=0)
                sum((1-e)*log(1-p+p*dpois(y,mu,log=FALSE))+e*(log(p)+dpois(y,mu,log=TRUE)))
            }
        }
        else if(family$family == "binomial") {
            d.V <- function(mu) 1-2*mu
            d.eta.mu <- function(mu) (2*mu-1)/(mu^2*(1-mu)^2)
            den <- function(y, mu) dbinom(y*size, size, mu)
            disp <- 1; est.disp <- FALSE
            loglikfun <- function(y, mu, p) {
                e <- as.numeric(y!=0)
                sum((1-e)*log(1-p+p*dbinom(y*size,size,mu,log=FALSE))+e*(log(p)+dbinom(y*size,size,mu,log=TRUE)))
            }
        }
        else stop("family not recognized")

        variance <- family$variance
        linkinv <- family$linkinv
        linkfun <- family$linkfun
        mu.eta <- family$mu.eta

        ## Estimate dispersion parameter using non-zero data...
        fm <- as.formula(sub(split$response,"y",deparse(formula)))

        ## Set initial values...
        if(family$family == "binomial") mu <- pmin(pmax(y, 0.01), 0.99)
        else mu <- pmax(y, 0.01)
        p <- rep(0.7, n)
        eta <- linkfun(mu) # g_mu(mu)
        gp <- eta # g_p(p)-alpha, equivalently alpha=0, delta=1

        ## Model setup
        psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
        G <- gam(fm, fit=FALSE, family=family, ...)

        np <- ncol(G$X) # length of parameter beta
        n.smooth <- length(G$smooth) # number of smooth terms
        if(n.smooth!=length(zero.delta))
            stop("length of 'zero.delta' does not match # of smooth terms")
        alpha <- 0
        zdelta <- !is.na(zero.delta)
        delta <- rep.int(1, n.smooth)
        delta[zdelta] <- 0 # fix zero delta's

        M <- list() # extract the design matrix of each smooth term
        for(k in 1:n.smooth) {
            M[[k]] <- G$X[, G$smooth[[k]]$first.para:G$smooth[[k]]$last.para]
        }
        X.par <- as.matrix(G$X[, 1:(G$smooth[[1]]$first.para-1)]) # design matrix of parametric part

        ## Begin outer loop: EM algorithm
        rep.out <- 1; norm.out <- 1
        while(norm.out>conv.crit.out & rep.out<maxiter) {

            psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
            b.old <- b <- rep.int(0, np)

            Xp <- matrix(0, ncol=ncol(X.par), nrow=nrow(X.par)) # the design matrix of the second SS
            for(k in 1:n.smooth) {
               Xp <- cbind(Xp, delta[k]*M[[k]])
            }

            ## inner loop: approx logLik by WSS
            ## update beta by P-IRLS
            rep.in <- 1; norm.in <- 1
            while(norm.in>conv.crit.in & rep.in<maxiter) {

                mu.eta.val <- mu.eta(eta)
                weights1 <- sqrt(psi*size*(mu.eta.val^2)/variance(mu)/disp)
                weights2 <- sqrt(p*(1-p))

                good <- (weights1 > 0) & (weights2 > 0) & (mu<1000) & (mu.eta.val != 0)
                good[is.na(good)] <- FALSE
                if (all(!good)) {
                    converged <- FALSE
                    warning("No observations informative at iteration")
                    break
                }

                z1 <- eta[good] + (y - mu)[good]/mu.eta.val[good]
                z2 <- gp[good] + 1/(p*(1-p))[good]*(psi-p)[good]

                y.c <- c(z1, z2)
                w.c <- c(weights1[good], weights2[good])
                X.c <- rbind(G$X[good,], Xp[good,])

                mgfit <- magic(y=y.c, X=X.c, sp=G$sp, S=G$S, off=G$off, rank=G$rank, C=G$C, w=w.c, gcv=FALSE)

                b.old <- b
                b <- mgfit$b
                sp <- mgfit$sp
                eta <- G$X%*%b

                Eta <- NULL # the design matrix of binomial model
                bp <- list() # extract the paramters for each smooth term
                eta.p <- list() # record the estimate the each smooth function
                gp <- 0
                for(k in 1:n.smooth) {
                    bp[[k]] <- b[G$smooth[[k]]$first.para:G$smooth[[k]]$last.para]
                    eta.p[[k]] <- M[[k]]%*%bp[[k]]
                    gp <- gp + delta[k]*eta.p[[k]]
                    if(!zdelta[k]) Eta <- cbind(Eta, eta.p[[k]])
                }

                mu <- linkinv(eta)
                p <- .Call("logit_linkinv", alpha + gp, PACKAGE = "stats")
                norm.in <- sum((b-b.old)^2)
                rep.in <- rep.in + 1
            }

            ## Update alpha & delta
            ## Fit a generalized linear model given beta
            alpha.old <- alpha
            delta.old <- delta

            if(is.null(Eta)) alpha <- coef(glm(psi ~ 1, family=quasibinomial)) # Homogeneous ZI
            else {
                fit.glm <- glm(psi ~ Eta, family=quasibinomial) # z ~ eta.p[[1]] +...+ eta.p[[k]]
                alpha <- fit.glm$coef[1]
                delta[!zdelta] <- fit.glm$coef[-1]
            }

            ## Update p for new alpha & delta
            gp <- 0
            for(k in 1:n.smooth) {
                gp <- gp + delta[k]*eta.p[[k]]
            }
            p <- .Call("logit_linkinv", alpha + gp, PACKAGE = "stats")

            if(is.null(Eta)) norm.out <- abs(alpha-alpha.old)
            else norm.out <- max(abs(alpha-alpha.old), abs(delta-delta.old))
            rep.out <- rep.out + 1
            cat("iteration =", rep.out, "\t", "norm =", norm.out, "\n")

        }

        names(alpha) <- names(delta) <- NULL
        theta <- list(alpha=alpha, delta=delta, beta=b)
        if(rep.out==maxiter) converged <- FALSE
            else converged <- TRUE

        ## Compute observed information...
        Lambda <- matrix(0, np, np) # penalty matrix

        rho <- psi-p; Ed.rho <- rep.int(-1, n)
        tau <- size*psi*(y-mu)*mu.eta.val/disp/variance(mu)
        Ed.tau <- size*psi*mu.eta.val^2*( -variance(mu)/mu.eta.val-(y-mu)*(variance(mu)*d.eta.mu(mu)+
            d.V(mu)/mu.eta.val) )/disp/(variance(mu)^2)


        Lam <- list(); n.S <- numeric(n.smooth);
        Xp <- matrix(0, ncol=ncol(X.par), nrow=nrow(X.par))
        P <- NULL
        for(k in 1:n.smooth) {
            n.S[k] <- length(G$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp[k]*G$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[j]*G$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp[sum(n.S[1:(k-1)])+1]*G$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[sum(n.S[1:(k-1)])+j]*G$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G$smooth[[k]]$first.para
            last <- G$smooth[[k]]$last.para
            Lambda[first:last, first:last] <- Lam[[k]]
            Xp <- cbind(Xp, delta[k]*M[[k]])
            if(!zdelta[k]) {
                Xrho <- t(G$X)%*%rho; Xrho[-(G$smooth[[k]]$first:G$smooth[[k]]$last)] <- 0
                P <- cbind(P, Xrho) # used in computing partial derivatives
            }
        }

        ## Find the effective degrees of freedom
        FM <- solve(t(X.c)%*%diag(w.c^2)%*%X.c+Lambda)%*%t(X.c)%*%diag(w.c^2)%*%X.c

        edf.smooth <- numeric(n.smooth)
        edf <- diag(FM)
        for(k in 1:n.smooth) {
            first <- G$smooth[[k]]$first.para
            last <- G$smooth[[k]]$last.para
            edf.smooth[k] <- sum(edf[first:last])
        }

        if (G$nsdf > 0)
            term.names <- colnames(G$X)[1:G$nsdf]
        else term.names <- array("", 0)
        for (i in 1:n.smooth) {
            k <- 1
            for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) {
                term.names[j] <- paste(G$smooth[[i]]$label, ".",
                  as.character(k), sep = "")
                k <- k + 1
            }
        }
        names(theta$beta) <- term.names
        names(edf) <- term.names

        G.rho <- diag(as.vector(Ed.rho*p*(1-p)))
        G.tau <- diag(as.vector(Ed.tau*mu.eta.val))
        W.tau <- diag(as.vector(psi*(1-psi)*size^2*(y-mu)^2*mu.eta.val^2/(disp^2)/(variance(mu)^2)))
        W.rho <- diag(as.vector(psi*(1-psi)))
        W.tau.rho <- diag(as.vector(psi*(1-psi)*size*(y-mu)*mu.eta.val/disp/variance(mu)))

        if(is.null(Eta)) {  # Homogeneous ZI
            I.alpha <- sum(-Ed.rho*p*(1-p)-psi*(1-psi))
            I.beta <- -t(G$X)%*%(G.tau+W.tau)%*%G$X - t(Xp)%*%(G.rho+W.rho)%*%Xp
                      - t(G$X)%*%W.tau.rho%*%Xp - t(Xp)%*%W.tau.rho%*%G$X + Lambda
            I.alpha.beta <- t(Xp)%*%(-Ed.rho*p*(1-p)) - (t(Xp)%*%W.rho+t(G$X)%*%W.tau.rho)%*%rep.int(1,n)

            n.theta <- 1 + length(b)
            n.delta <- 0
            I.theta <- matrix(NA, ncol=n.theta, nrow=n.theta)
            I.theta[1,1] <- I.alpha
            I.theta[1,2:n.theta] <- I.alpha.beta
            I.theta[2:n.theta,1] <- I.alpha.beta
            I.theta[2:n.theta,2:n.theta] <- I.beta

            colnames(I.theta) <- rownames(I.theta) <- NULL
            V.theta <- solve(I.theta)
            se.alpha <- sqrt(V.theta[1,1])
            se.delta <- rep.int(0, n.smooth)
            V.b <- V.theta[-1,-1]
        }

        else {
            I.alpha <- sum(-Ed.rho*p*(1-p)-psi*(1-psi))
            I.delta <- -t(Eta)%*%(G.rho+W.rho)%*%Eta
            I.alpha.delta <- -t(Eta)%*%(Ed.rho*p*(1-p)+psi*(1-psi))

            I.beta <- -t(G$X)%*%(G.tau+W.tau)%*%G$X - t(Xp)%*%(G.rho+W.rho)%*%Xp
                      - t(G$X)%*%W.tau.rho%*%Xp - t(Xp)%*%W.tau.rho%*%G$X + Lambda
            I.alpha.beta <- t(Xp)%*%(-Ed.rho*p*(1-p)) - (t(Xp)%*%W.rho+t(G$X)%*%W.tau.rho)%*%rep.int(1,n)
            I.delta.beta <- -t(Xp)%*%G.rho%*%Eta - (t(Xp)%*%W.rho+t(G$X)%*%W.tau.rho)%*%Eta

            n.delta <- sum(!zdelta) # number of non-zero delta
            n.theta <- 1 + n.delta + length(b)
            I.theta <- matrix(NA, ncol=n.theta, nrow=n.theta)
            I.theta[1,1] <- I.alpha; I.theta[1,2:(n.delta+1)] <- I.alpha.delta
            I.theta[1,(n.delta+2):n.theta] <- I.alpha.beta
            I.theta[2:(n.delta+1), 1] <- I.alpha.delta
            I.theta[2:(n.delta+1),2:(n.delta+1)] <- I.delta
            I.theta[2:(n.delta+1),(n.delta+2):n.theta] <- t(I.delta.beta)
            I.theta[(n.delta+2):n.theta,1] <- I.alpha.beta
            I.theta[(n.delta+2):n.theta,2:(n.delta+1)] <- I.delta.beta
            I.theta[(n.delta+2):n.theta,(n.delta+2):n.theta] <- I.beta

            colnames(I.theta) <- rownames(I.theta) <- NULL
            V.theta <- solve(I.theta)
            se.alpha <- sqrt(V.theta[1,1])
            se.delta <- rep.int(0, n.smooth)
            se.delta[!zdelta] <- sqrt(diag(V.theta)[2:(n.delta+1)])
            V.b <- V.theta[-(1:(n.delta+1)),-(1:(n.delta+1))]
        }

        loglik <- loglikfun(y, mu, p)  # log-likelihood
        ploglik <- loglik - as.numeric(0.5*t(b)%*%Lambda%*%b)  # penalized log-likelihood

        # Laplace approximation of log marginal likelihood
        DS <- diag(eigen(Lambda)$values[abs(eigen(Lambda)$values)>1e-10])
        logE <- 0.5*determinant(DS)$modulus+ploglik+(np+1+n.delta-nrow(DS))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
        attr(logE, "logarithm") <- NULL

        cat("\n")
        cat("==========================================","\n")
        cat("estimated alpha =", alpha, "(", se.alpha, ")", "\n")
        for (k in 1:n.smooth) {
            cat(paste("estimated delta",k, sep=""),"=", delta[k], "(", se.delta[k], ")", "\n")
        }
        cat("==========================================","\n","\n")

        res <- list(coefficient=theta, V.theta=V.theta, converged=converged,
              mu=mu, linear.predictor=eta, dispersion=disp, formula=formula,
              score=mgfit$score, G=G, V.beta=V.b, y=y, p=p, psi=psi, loglik=loglik,
              ploglik=ploglik, family=family, edf=edf, est.disp=est.disp,
              edf.smooth=edf.smooth, sp=sp, se.delta=se.delta, se.alpha=se.alpha,
              logE=logE)
        class(res) <- "cozigam"
        res

}

# PCOZIGAM.cts(COZIGAM)
# Proportional COnstrained Zero-Inflated GAM: Continuous Case

PCOZIGAM.cts <- function(formula, maxiter = 20, conv.crit.in = 1e-5,
    conv.crit.out = 1e-3, log.tran = FALSE, family = gaussian(), ...)
{
        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        n <- length(y)
        z <- (y!=0)
        if(log.tran==TRUE) y[z] <- log(y[z])

        if (is.character(family))
            family <- eval(parse(text = family))
        if (is.function(family))
            family <- family()
        if(family$family == "gaussian") {
            d.V <- function(mu) rep.int(0, length(mu))
            d.eta.mu <- function(mu) rep.int(0, length(mu))
            loglikfun <- function(y, mu, disp, z, p) {
                sum(z*(log(p)+dnorm(y,mu,sqrt(disp),log=TRUE))+(1-z)*log(1-p))
            }
        }
        else if(family$family == "Gamma") {
            d.V <- function(mu) 2*mu
            d.eta.mu <- function(mu) 2/(mu^3)
            loglikfun <- function(y, mu, disp, z, p) {
                sum(z*(log(p)+ifelse(y==0, 1, dgamma(y,1/disp,scale=mu*disp,log=TRUE)))+(1-z)*log(1-p))
            }
        }
        else
            stop("family not recognized")

        variance <- family$variance
        linkinv <- family$linkinv
        linkfun <- family$linkfun
        mu.eta <- family$mu.eta

        fm <- as.formula(sub(split$response,"y",deparse(formula)))

        ## Estimate dispersion parameter using non-zero data...
        fit.nz <- gam(fm, subset=z, family=family, ...)
        G1 <- gam(fm, subset=z, family=family, fit=FALSE, ...)
        disp <- fit.nz$sig2; est.disp <- TRUE

        ## Set initial values...
        alpha <- 0; delta <- 1

        mu <- pmax(y, 0.01)
        p <- rep(0.7, n)
        eta1 <- linkfun(mu); eta2 <- delta*eta1

        ## Model setup
        G <- gam(fm, fit=FALSE, family=family, ...)

        ## Begin outer loop: iteratively update beta and alpha & delta
        rep.out <- 1; norm.out <- 1
        while(norm.out>conv.crit.out & rep.out<maxiter) {

            b.old <- b <- rep(0, ncol(G$X))

            ## inner loop: approx logLik by WSS
            ## update beta by P-IRLS
            rep.in <- 1; norm.in <- 1
            while(norm.in>conv.crit.in & rep.in<maxiter) {

                mu.eta.val <- mu.eta(eta1)
                weights1 <- sqrt(z*(mu.eta.val^2)/variance(mu)/disp)
                weights2 <- sqrt(p*(1-p))

                good <- (weights1 > 0) & (weights2 > 0) & (mu<1000) & (mu.eta.val != 0)
                good[is.na(good)] <- FALSE
                if (all(!good)) {
                    converged <- FALSE
                    warning("No observations informative at iteration")
                    break
                }
                z1 <- eta1[good] + (y - mu)[good]/mu.eta.val[good]
                z2 <- eta2[good] + 1/(p*(1-p))[good]*(z-p)[good]

                y.c <- c(z1, z2)
                w.c <- c(weights1[good], weights2[good])
                X.c <- rbind(G$X[good,], delta*G$X[good,])

                mgfit <- magic(y=y.c, X=X.c, sp=G$sp, S=G1$S, off=G$off, rank=G$rank, C=G$C, w=w.c, gcv=FALSE)

                b.old <- b
                b <- mgfit$b
                sp <- mgfit$sp
                eta1 <- G$X%*%b; eta2 <- delta*eta1
                mu <- linkinv(eta1)
                p <- .Call("logit_linkinv", eta2+alpha, PACKAGE = "stats")
                norm.in <- sum((b-b.old)^2)
                rep.in <- rep.in + 1
            }

            ## Update alpha & delta
            ## Fit a generalized linear model given eta1
            alpha.old <- alpha
            delta.old <- delta

            fit.glm <- glm(z ~ eta1, family=binomial)
            alpha <- fit.glm$coef[1]
            delta <- fit.glm$coef[2]

            ## Update p for new alpha & delta
            eta2 <- delta*eta1
            p <- .Call("logit_linkinv", eta2+alpha, PACKAGE = "stats")

            norm.out <- max(abs(alpha-alpha.old), abs(delta-delta.old))
            rep.out <- rep.out + 1
            cat("iteration =", rep.out, "\t", "norm =", norm.out, "\n")

        }

        names(alpha) <- names(delta) <- NULL
        theta <- list(alpha=alpha, delta=delta, beta=b)
        if(rep.out==maxiter) converged <- FALSE
            else converged <- TRUE

        ## Compute observed information...
        n.smooth <- length(G$smooth)
        np <- length(b)
            Lambda <- matrix(0, np, np)


        Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
        for(k in 1:n.smooth) {
            n.S[k] <- length(G1$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp[k]*G1$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[j]*G1$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp[sum(n.S[1:(k-1)])+1]*G1$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[sum(n.S[1:(k-1)])+j]*G1$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G1$smooth[[k]]$first.para
            last <- G1$smooth[[k]]$last.para
            Lambda[first:last, first:last] <- Lam[[k]]
        }

        ## Find the effective degrees of freedom
        FM <- solve(t(X.c)%*%diag(w.c^2)%*%X.c+Lambda)%*%t(X.c)%*%diag(w.c^2)%*%X.c

        edf.smooth <- numeric(n.smooth)
        edf <- diag(FM)
        for(k in 1:n.smooth) {
            first <- G$smooth[[k]]$first.para
            last <- G$smooth[[k]]$last.para
            edf.smooth[k] <- sum(edf[first:last])
        }

        if (G$nsdf > 0)
            term.names <- colnames(G$X)[1:G$nsdf]
        else term.names <- array("", 0)
        for (i in 1:n.smooth) {
            k <- 1
            for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) {
                term.names[j] <- paste(G$smooth[[i]]$label, ".",
                  as.character(k), sep = "")
                k <- k + 1
            }
        }
        names(theta$beta) <- term.names
        names(edf) <- term.names

        rho <- z-p; d.rho <- rep.int(-1, n)
        tau <- z*(y-mu)*mu.eta.val/disp/variance(mu)
        d.tau <- z*mu.eta.val^2*( -variance(mu)/mu.eta.val-(y-mu)*(variance(mu)*d.eta.mu(mu)+
            d.V(mu)/mu.eta.val) )/disp/(variance(mu)^2)

        G.rho <- diag(as.vector(d.rho*p*(1-p)))
        G.tau <- diag(as.vector(d.tau*mu.eta.val))
        eta1 <- as.vector(eta1)

        I.alpha <- sum(diag(G.rho))
        I.delta <- t(eta1)%*%G.rho%*%eta1
        I.alpha.delta <- t(rep.int(1, n))%*%G.rho%*%eta1
        I.beta <- t(G$X)%*%G.tau%*%G$X + delta^2*t(G$X)%*%G.rho%*%G$X - Lambda
        I.alpha.beta <- delta*t(G$X)%*%G.rho%*%rep.int(1, n)
        I.delta.beta <- delta*t(G$X)%*%G.rho%*%eta1 + t(G$X)%*%rho

        I.theta <- matrix(0, ncol=np+2, nrow=np+2)
        I.theta[1,1] <- I.alpha; I.theta[1,2] <- I.theta[2,1] <- I.alpha.delta
        I.theta[3:(np+2), 1] <- I.theta[1,3:(np+2)] <- I.alpha.beta
        I.theta[2,2] <- I.delta
        I.theta[3:(np+2), 2] <- I.theta[2,3:(np+2)] <- I.delta.beta
        I.theta[3:(np+2),3:(np+2)] <- I.beta
        I.theta <- -I.theta

        V.theta <- solve(I.theta)
        se.alpha <- sqrt(V.theta[1,1])
        se.delta <- sqrt(V.theta[2,2])
        V.beta <- V.theta[3:(np+2),3:(np+2)]

        loglik <- loglikfun(y, mu, disp, z, p) # log-likelihood
        ploglik <- loglik - as.numeric(0.5*t(b)%*%Lambda%*%b) # penalized log-likelihood

        # Model comparison/selection
        # Strictly positive eigenvalues of the penalty matrix Lambda
        DS <- diag(eigen(Lambda)$values[abs(eigen(Lambda)$values)>1e-10])
        # Laplace approximation of log marginal likelihood
        logE <- 0.5*determinant(DS)$modulus+ploglik+(np+2-nrow(DS))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
        attr(logE, "logarithm") <- NULL

        cat("\n")
        cat("==========================================","\n")
        cat("estimated alpha =", alpha, "(", se.alpha, ")", "\n")
        cat("estimated delta =", delta, "(", se.delta, ")", "\n")
        cat("==========================================","\n","\n")

        res <- list(coefficient=theta, V.theta=V.theta, converged=converged,
              mu=mu, linear.predictor=eta1, dispersion=disp, formula=formula,
              fit.nonzero=fit.nz, score=mgfit$score, G=G, V.beta=V.beta,
              y=y, p=p, loglik=loglik, ploglik=ploglik, family=family, edf=edf,
              est.disp=est.disp, edf.smooth=edf.smooth, sp=sp, se.delta=se.delta,
              se.alpha=se.alpha, logE=logE)
        class(res) <- "cozigam"
        res

}

# PCOZIGAM.dis(COZIGAM)
# Proportional COnstrained Zero-Inflated GAM: Distrete Case

PCOZIGAM.dis <- function(formula, maxiter = 30, conv.crit.in = 1e-5,
    conv.crit.out = 1e-3, size = NULL, family = poisson(), ...)
{
        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        n <- length(y)

        if (is.character(family))
            family <- eval(parse(text = family))
        if (is.function(family))
            family <- family()
        if(family$family == "poisson") {
            d.V <- function(mu) rep.int(1, length(mu))
            d.eta.mu <- function(mu) -1/(mu^2)
            den <- dpois; disp <- 1; est.disp <- FALSE
            size <- rep.int(1, n)
            loglikfun <- function(y, mu, p) {
                e <- as.numeric(y!=0)
                sum((1-e)*log(1-p+p*dpois(y,mu,log=FALSE))+e*(log(p)+dpois(y,mu,log=TRUE)))
            }
        }
        else if(family$family == "binomial") {
            d.V <- function(mu) 1-2*mu
            d.eta.mu <- function(mu) (2*mu-1)/(mu^2*(1-mu)^2)
            den <- function(y, mu) dbinom(y*size, size, mu)
            disp <- 1; est.disp <- FALSE
            loglikfun <- function(y, mu, p) {
                e <- as.numeric(y!=0)
                sum((1-e)*log(1-p+p*dbinom(y*size,size,mu,log=FALSE))+e*(log(p)+dbinom(y*size,size,mu,log=TRUE)))
            }
        }
        else stop("family not recognized")

        variance <- family$variance
        linkinv <- family$linkinv
        linkfun <- family$linkfun
        mu.eta <- family$mu.eta

        fm <- as.formula(sub(split$response,"y",deparse(formula)))

        ## Set initial values...
        alpha <- 0; delta <- 1

        if(family$family == "binomial") mu <- pmin(pmax(y, 0.01), 0.99)
        else mu <- pmax(y, 0.01)
        p <- rep(0.6, n)
        eta1 <- linkfun(mu); eta2 <- delta*eta1

        ## Model setup
        psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
        G <- gam(fm, fit=FALSE, family=family, ...)

        ## Begin outer loop: EM algorithm
        rep.out <- 1; norm.out <- 1
        while(norm.out>conv.crit.out & rep.out<maxiter) {

            psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
            b.old <- b <- rep(0, ncol(G$X))

            ## inner loop: approx logLik by WSS
            ## updating beta by P-IRLS
            rep.in <- 1; norm.in <- 1
            while(norm.in>conv.crit.in & rep.in<maxiter) {

                mu.eta.val <- mu.eta(eta1)
                weights1 <- sqrt(psi*size*(mu.eta.val^2)/variance(mu)/disp)
                weights2 <- sqrt(p*(1-p))

                good <- (weights1 > 0) & (weights2 > 0) & (mu<1000) & (mu.eta.val != 0)
                good[is.na(good)] <- FALSE
                if (all(!good)) {
                    converged <- FALSE
                    warning("No observations informative at iteration")
                    break
                }
                z1 <- eta1[good] + (y - mu)[good]/mu.eta.val[good]
                z2 <- eta2[good] + 1/(p*(1-p))[good]*(psi-p)[good]

                y.c <- c(z1, z2)
                w.c <- c(weights1[good], weights2[good])
                X.c <- rbind(G$X[good,], delta*G$X[good,])

                mgfit <- magic(y=y.c, X=X.c, sp=G$sp, S=G$S, off=G$off, rank=G$rank, C=G$C, w=w.c, gcv=FALSE)

                b.old <- b
                b <- mgfit$b
                sp <- mgfit$sp
                eta1 <- G$X%*%b; eta2 <- delta*eta1
                mu <- linkinv(eta1)
                p <- .Call("logit_linkinv", eta2+alpha, PACKAGE = "stats")
                norm.in <- sum((b-b.old)^2)
                rep.in <- rep.in + 1

            }

            ## Update alpha & delta
            ## Fit a generalized linear model given eta1
            alpha.old <- alpha
            delta.old <- delta

            fit.glm <- glm(psi ~ eta1, family=quasibinomial)
            alpha <- fit.glm$coef[1]
            delta <- fit.glm$coef[2]

            ## Update p for new alpha & delta
            eta2 <- delta*eta1
            p <- .Call("logit_linkinv", eta2+alpha, PACKAGE = "stats")

            norm.out <- max(abs(alpha-alpha.old), abs(delta-delta.old))
            rep.out <- rep.out + 1
            cat("iteration =", rep.out, "\t", "norm =", norm.out, "\n")

        }

        names(alpha) <- names(delta) <- NULL
        theta <- list(alpha=alpha, delta=delta, beta=b)
        if(rep.out==maxiter) converged <- FALSE
            else converged <- TRUE

        ## Compute observed information...
        n.smooth <- length(G$smooth) # penalty matrix
        np <- length(b)
        Lambda <- matrix(0, np, np)
        Lam <- list(); n.S <- numeric(n.smooth)
        for(k in 1:n.smooth) {
            n.S[k] <- length(G$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp[k]*G$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[j]*G$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp[sum(n.S[1:(k-1)])+1]*G$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[sum(n.S[1:(k-1)])+j]*G$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G$smooth[[k]]$first.para
            last <- G$smooth[[k]]$last.para
            Lambda[first:last, first:last] <- Lam[[k]]
        }

        ## Find the effective degrees of freedom
        FM <- solve(t(X.c)%*%diag(w.c^2)%*%X.c+Lambda)%*%t(X.c)%*%diag(w.c^2)%*%X.c

        edf.smooth <- numeric(n.smooth)
        edf <- diag(FM)
        for(k in 1:n.smooth) {
            first <- G$smooth[[k]]$first.para
            last <- G$smooth[[k]]$last.para
            edf.smooth[k] <- sum(edf[first:last])
        }

        if (G$nsdf > 0)
            term.names <- colnames(G$X)[1:G$nsdf]
        else term.names <- array("", 0)
        for (i in 1:n.smooth) {
            k <- 1
            for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) {
                term.names[j] <- paste(G$smooth[[i]]$label, ".",
                  as.character(k), sep = "")
                k <- k + 1
            }
        }
        names(theta$beta) <- term.names
        names(edf) <- term.names

        rho <- psi-p; Ed.rho <- rep.int(-1, n)
        tau <- size*psi*(y-mu)*mu.eta.val/disp/variance(mu)
        Ed.tau <- size*psi*mu.eta.val^2*( -variance(mu)/mu.eta.val-(y-mu)*(variance(mu)*d.eta.mu(mu)+
            d.V(mu)/mu.eta.val) )/disp/(variance(mu)^2)

        I.alpha <- sum(-Ed.rho*p*(1-p)-psi*(1-psi))
        I.delta <- sum(eta1^2*(-Ed.rho*p*(1-p)-psi*(1-psi)))
        I.alpha.delta <- sum(eta1*(-Ed.rho*p*(1-p)-psi*(1-psi)))
        I.beta <- t(G$X)%*%(diag(as.vector(-Ed.tau*mu.eta.val-psi*(1-psi)*size^2*(y-mu)^2*mu.eta.val^2/(disp^2)/(variance(mu)^2))))%*%G$X +
                  delta^2*t(G$X)%*%(diag(as.vector(-Ed.rho*p*(1-p)-psi*(1-psi))))%*%G$X -
                  2*delta*t(G$X)%*%(diag(as.vector(psi*(1-psi)*size*(y-mu)*mu.eta.val/disp/variance(mu))))%*%G$X + Lambda
        I.alpha.beta <- -delta*t(G$X)%*%(Ed.rho*p*(1-p)+psi*(1-psi)) -
                        t(G$X)%*%(psi*(1-psi)*size*(y-mu)*mu.eta.val/disp/variance(mu))
        I.delta.beta <- -delta*t(G$X)%*%(eta1*(Ed.rho*p*(1-p)+psi*(1-psi))) -
                        t(G$X)%*%(psi-p+eta1*psi*(1-psi)*size*(y-mu)*mu.eta.val/disp/variance(mu))

        I.theta <- cbind(rbind(I.alpha,I.alpha.delta,I.alpha.beta),
                         rbind(I.alpha.delta,I.delta,I.delta.beta),
                         rbind(t(I.alpha.beta),t(I.delta.beta),I.beta))

        colnames(I.theta) <- rownames(I.theta) <- NULL
        V.theta <- solve(I.theta)
        se.alpha <- sqrt(V.theta[1,1])
        se.delta <- sqrt(V.theta[2,2])
        V.b <- V.theta[-(1:2),-(1:2)]

        cat("\n")
        cat("==========================================","\n")
        cat("estimated alpha =", alpha, "(", se.alpha, ")", "\n")
        cat("estimated delta =", delta, "(", se.delta, ")", "\n")
        cat("==========================================","\n","\n")

        loglik <- loglikfun(y, mu, p) # log-likelihood
        ploglik <- loglik - as.numeric(0.5*t(b)%*%Lambda%*%b)  # penalized log-likelihood

        # Model comparison/selection
        # Strictly positive eigenvalues of the penalty matrix Lambda
        DS <- diag(eigen(Lambda)$values[abs(eigen(Lambda)$values)>1e-10])
        # Laplace approximation of log marginal likelihood
        logE <- 0.5*determinant(DS)$modulus+ploglik+(np+2-nrow(DS))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
        attr(logE, "logarithm") <- NULL

        res <- list(coefficient=theta, V.theta=V.theta, converged=converged,
              mu=mu, linear.predictor=eta1, dispersion=disp, formula=formula,
              score=mgfit$score, G=G, V.beta=V.b, y=y*size, psi=psi,
              p=p, family=family, edf=edf, loglik=loglik, ploglik=ploglik,
              est.disp=est.disp, edf.smooth=edf.smooth, se.delta=se.delta,
              se.alpha=se.alpha, logE=logE)
        class(res) <- "cozigam"
        res

}

# disgam(COZIGAM)
# Log marginal likelihood for discrete GAM

disgam <- function(formula, size=NULL, family = poisson(), ...)
{
        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        n <- length(y)

        if (is.character(family))
            family <- eval(parse(text = family))
        if (is.function(family))
            family <- family()
        if(family$family == "poisson") {
            d.V <- function(mu) rep.int(1, length(mu))
            d.eta.mu <- function(mu) -1/(mu^2)
            den <- dpois; disp <- 1; est.disp <- FALSE
            size <- rep.int(1, n)
            loglikfun <- function(y, mu) sum(dpois(y,mu,log=TRUE))
        }
        else if(family$family == "binomial") {
            d.V <- function(mu) 1-2*mu
            d.eta.mu <- function(mu) (2*mu-1)/(mu^2*(1-mu)^2)
            den <- function(y, mu) dbinom(y*size, size, mu)
            disp <- 1; est.disp <- FALSE
            loglikfun <- function(y, mu) sum(dbinom(y*size,size,mu,log=TRUE))
        }
        else stop("family not recognized")
        variance <- family$variance
        linkinv <- family$linkinv
        linkfun <- family$linkfun
        mu.eta <- family$mu.eta

        fm <- as.formula(sub(split$response,"y",deparse(formula)))
        G <- gam(fm, fit=FALSE, family=family, weights=size, ...)
        fit.gam <- gam(fm, family=family, weights=size, ...)
        b <- fit.gam$coefficients
        sp <- fit.gam$sp
        mu <- fit.gam$fitted.values

        n.smooth <- length(G$smooth) # penalty matrix
        np <- length(b)
        Lambda <- matrix(0, np, np)
        Lam <- list(); n.S <- numeric(n.smooth)
        for(k in 1:n.smooth) {
            n.S[k] <- length(G$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp[k]*G$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[j]*G$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp[sum(n.S[1:(k-1)])+1]*G$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp[sum(n.S[1:(k-1)])+j]*G$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G$smooth[[k]]$first.para
            last <- G$smooth[[k]]$last.para
            Lambda[first:last, first:last] <- Lam[[k]]
        }
        mu.eta.val <- mu.eta(fit.gam$linear.predictors)
        d.tau <- size*mu.eta.val^2*( -variance(mu)/mu.eta.val-(y-mu)*(variance(mu)*d.eta.mu(mu)+
            d.V(mu)/mu.eta.val) )/disp/(variance(mu)^2)
        I.beta <- t(G$X)%*%(diag(as.vector(-d.tau*mu.eta.val)))%*%G$X + Lambda
        V.beta <- solve(I.beta)
        loglik <- loglikfun(y, mu) # log-likelihood
        ploglik <- loglik - as.numeric(0.5*t(b)%*%Lambda%*%b)  # penalized log-likelihood

        DS <- diag(eigen(Lambda)$values[abs(eigen(Lambda)$values)>1e-10])
        # Laplace approximation of log marginal likelihood
        logE <- 0.5*determinant(DS)$modulus+ploglik+(np-nrow(DS))/2*log(2*pi)-0.5*determinant(I.beta)$modulus
        attr(logE, "logarithm") <- NULL

        res <- list(fit.gam=fit.gam, logE=logE, V.beta=V.beta, loglik=loglik, ploglik=ploglik,
          mu=mu, formula=formula, family=family)

}


# plot.cozigam(COZIGAM)

plot.cozigam <- function(x, plot.2d = "contour", too.far = 0,
    n.1d = 100, n.2d = 30, theta = 30, phi = 30, select = NULL, image.col = "topo",
    persp.col = "lightblue", contour.col = "red", n.Col = 100, shade.ci = FALSE,
    shade.col = "gray80", Rug = TRUE, xlab, ylab, ...)
{
      if (!inherits(x, "cozigam"))
          stop("use only with \"COZIGAM\" objects")
      require(mgcv)
      G <- x$G
      n.smooth <- length(G$smooth)
      b <- x$coefficient$beta; V.b <- x$V.beta

      opar <- par(ask=FALSE)
      on.exit(par(opar))
      par(ask = TRUE)
      if(is.null(select)) plot.seq <- 1:n.smooth
      else plot.seq <- select

      for(k in plot.seq) {

          first <- G$smooth[[k]]$first.para
          last <- G$smooth[[k]]$last.para

          if(G$smooth[[k]]$dim==1) {
                raw <- c(unlist(G$mf[G$smooth[[k]]$term]))
                xx <- seq(min(raw), max(raw), length=n.1d)
                if (G$smooth[[k]]$by != "NA") {
                      by <- rep(1, n.1d)
                      pr <- data.frame(xx=xx, by=by)
                      colnames(pr) <- c(G$smooth[[k]]$term, G$smooth[[k]]$by)
                }
                else {
                      pr <- data.frame(xx=xx)
                      colnames(pr) <- G$smooth[[k]]$term
                }
                X.pred <- PredictMat(G$smooth[[k]], pr)
                fit <- X.pred%*%b[first:last]
                se.fit <- sqrt(rowSums((X.pred%*%V.b[first:last,first:last])*X.pred))
                ci.upper <- fit + 2*se.fit
                ci.lower <- fit - 2*se.fit
                if(missing(xlab)) xlab <- colnames(pr)
                if(missing(ylab)) ylab <- NA

                if (shade.ci) {
                    plot(fit~xx, type="l", col="black", xlab=xlab, ylab=ylab,
                        ylim=c(min(ci.lower),max(ci.upper)), ...)
                    polygon(c(xx, xx[n.1d:1], xx[1]), c(ci.upper, ci.lower[n.1d:1], ci.upper[1]),
                        col = shade.col, border = NA)
                    lines(xx, fit)
                    if(Rug) {
                        rug(raw, col="black")
                    }
                }
                else {
                    plot(fit~xx, type="l", col="black", xlab=xlab, ylab=ylab,
                        ylim=c(min(ci.lower),max(ci.upper)), ...)
                    lines(ci.upper~xx, type="l", col="black", lty=2)
                    lines(ci.lower~xx, type="l", col="black", lty=2)
                    if(Rug) {
                        rug(raw, col="black")
                    }
                }
          }

          else if(G$smooth[[k]]$dim==2) {
                raw <- list(xx = c(unlist(G$mf[G$smooth[[k]]$term[1]])),
                            zz = c(unlist(G$mf[G$smooth[[k]]$term[2]])) )
                xs <- seq(min(raw$xx), max(raw$xx), length=n.2d)
                zs <- seq(min(raw$zz), max(raw$zz), length=n.2d)
                gx <- rep(xs,n.2d); gz <- rep(zs,rep(n.2d,n.2d))

                if (G$smooth[[k]]$by != "NA") {
                      by <- rep(1, n.2d*n.2d)
                      pr <- data.frame(xx=gx, zz=gz, by=by)
                      colnames(pr) <- c(G$smooth[[k]]$term, G$smooth[[k]]$by)
                }
                else {
                      pr <- data.frame(xx=gx, zz=gz)
                      colnames(pr) <- G$smooth[[k]]$term
                }
                X.pred <- PredictMat(G$smooth[[k]], pr)
                fit <- matrix(X.pred%*%b[first:last], n.2d, n.2d)
                if (too.far > 0)
                    exclude <- exclude.too.far(gx, gz, raw$xx, raw$zz, dist = too.far)
                else exclude <- rep(FALSE, n.2d*n.2d)
                fit[exclude] <- NA
                if(missing(xlab)) xlab <- G$smooth[[k]]$term[1]
                if(missing(ylab)) ylab <- G$smooth[[k]]$term[2]

                if(plot.2d=="contour") {
                      if(is.null(image.col)) contour(xs, zs, fit, col=contour.col, xlab=xlab, ylab=ylab, ...)
                      else {
                            if(image.col=="terrain") Col <- terrain.colors(n.Col)
                            else if(image.col=="topo") Col <- topo.colors(n.Col)
                            else if(image.col=="heat") {
                                  Col <- heat.colors(n.Col)
                                  contour.col <- "blue"
                            }
                            else if(image.col=="gray") {
                                  Col <- gray(seq(0.1, 0.9, length = n.Col))
                                  contour.col <- "black"
                            }
                            else stop("color scheme not recognised")

                            image(xs, zs, fit, col=Col, xlab=xlab, ylab=ylab, ...)
                            contour(xs, zs, fit, col=contour.col, add=TRUE)
                      }
                      if (Rug) {
                            if (is.null(list(...)[["pch"]]))
                                  points(raw$xx, raw$zz, pch = ".", ...)
                            else points(raw$xx, raw$zz, ...)
                      }
                }
                else if(plot.2d=="persp")
                      persp(xs, zs, fit, col=persp.col, theta=theta, phi=phi, xlab=xlab, ylab=ylab, ...)
                else stop("2-dim plot unavailable")
          }
	    else warning("dim>=3 cannot be visualized")
      }

}


# predict.cozigam(COZIGAM)

predict.cozigam <- function(object, newdata, type="link", se.fit=FALSE, ...)

{
      if (!inherits(object, "cozigam"))
          stop("use only with \"COZIGAM\" objects")
      if (type != "link" && type != "terms" && type != "response") {
          warning("Unknown type, reset to link.")
          type <- "link"
      }
      if (missing(newdata)) newdata <- NULL
      require(mgcv)
      G <- object$G

      if (missing(newdata) | is.null(newdata)) newdata <- G$mf
      else {
          Terms <- delete.response(terms(G))
          allNames <- all.vars(Terms)
          ff <- reformulate(allNames)
          if (sum(!(allNames %in% names(newdata)))) {
              stop("not all required variables have been supplied in  newdata!\n")
          }
          newdata <- eval(model.frame(ff, data = newdata), parent.frame())
      }

      linkinv <- G$family$linkinv
      dmu.deta <- G$family$mu.eta
      n.smooth <- length(G$smooth)
      n.obs <- nrow(newdata)
      b <- object$coefficient$beta; V.b <- object$V.beta
      alpha <- object$coefficient$alpha
      delta <- object$coefficient$delta

      se.terms <- fit <- matrix(0, nrow=n.obs, ncol=n.smooth)

      X.pred <- NULL

      for(k in 1:n.smooth) {

          first <- G$smooth[[k]]$first.para
          last <- G$smooth[[k]]$last.para

          pr <- newdata[G$smooth[[k]]$term]
          if (G$smooth[[k]]$by != "NA") {
                by <- rep(1, n.obs)
                pr$by <- by
                colnames(pr) <- c(G$smooth[[k]]$term, G$smooth[[k]]$by)
          }
          else colnames(pr) <- G$smooth[[k]]$term

          X.pred[[k]] <- PredictMat(G$smooth[[k]], pr)
          fit[, k] <- X.pred[[k]]%*%b[first:last]

          if(k==1) X.pd <- X.pred[[k]]
          else X.pd <- cbind(X.pd, X.pred[[k]])

          if(se.fit)
              se.terms[, k] <- sqrt(rowSums((X.pred[[k]]%*%V.b[first:last,first:last])*X.pred[[k]]))

      }

      if(attr(object, "constraint") == "proportional") {
          eta.fit <- apply(fit, 1, sum) + b[1]
          mu <- linkinv(eta.fit)
          eta.p <- alpha + delta*eta.fit
          p <- .Call("logit_linkinv", eta.p, PACKAGE = "stats")
      }
      else {
          eta.fit <- apply(fit, 1, sum) + b[1]
          mu <- linkinv(eta.fit)
          eta.p <- alpha + apply(t(t(fit)*delta), 1, sum)
          p <- .Call("logit_linkinv", eta.p, PACKAGE = "stats")
      }

      if(se.fit) {
          first <- G$smooth[[1]]$first.para
          last <- G$smooth[[n.smooth]]$last.para
          se <- sqrt(rowSums((X.pd%*%V.b[first:last,first:last])*X.pd))
      }

      if(type=="link") {
          if(se.fit) pred <- data.frame(fit=eta.fit, se=se, p=p)
          else pred <- data.frame(fit=eta.fit, p=p)
      }

      if(type=="terms") {
          if(se.fit) pred <- list(fit=fit, se=se.terms, p=p)
          else pred <- data.frame(fit=fit, p=p)
      }

      if(type=="response") {
          if(se.fit) pred <- data.frame(fit=mu, se=se*abs(dmu.deta(eta.fit)), p=p)
          else pred <- data.frame(fit=mu, p=p)
      }

      pred

}


# print.cozigam(COZIGAM)

print.cozigam <- function (x, ...)
{
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    alpha <- x$coef$alpha
    delta <- x$coef$delta
    se.alpha <- x$se.alpha
    se.delta <- x$se.delta
    n.smooth <- length(x$G$smooth)
    cat("\nCoefficients of Constraint","\n")
    cat(" alpha =", round(alpha,5), "(", round(se.alpha,5), ")", "\n")
    for (k in 1:n.smooth) {
        cat(paste(" delta",k, sep=""),"=", round(delta[k],5), "(", round(se.delta[k],5), ")", "\n")
    }

    if (n.smooth == 0)
        cat("Total model degrees of freedom", sum(x$edf), "\n")
    else {
        edf <- 0
        for (i in 1:n.smooth) edf[i] <- sum(x$edf[x$G$smooth[[i]]$first.para:x$G$smooth[[i]]$last.para])
        cat("\nEstimated degrees of freedom:\n", edf, "  total = ",
            sum(x$edf), "\n")
    }
}


# print.summary.cozigam(COZIGAM)

print.summary.cozigam <- function (x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    if (length(x$p.coeff) > 0) {
        cat("\nParametric coefficients:\n")
        printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars,
            na.print = "NA", ...)
    }
    cat("\n")
    if (x$m > 0) {
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars,
            has.Pvalue = TRUE, na.print = "NA", cs.ind = 1, ...)
    }
    cat("\nScale est. = ", formatC(x$dispersion, digits = 5, width = 8,
        flag = "-"), "  n = ", x$n, "\n", sep = "")
}


# summary.cozigam(COZIGAM)

summary.cozigam <- function(object, dispersion = NULL, ...)
{

    pinv <- function(V, M, rank.tol = 1e-06) {
        D <- La.svd(V)
        M1 <- length(D$d[D$d > rank.tol * D$d[1]])
        if (M > M1)
            M <- M1
        if (M + 1 <= length(D$d))
            D$d[(M + 1):length(D$d)] <- 1
        D$d <- 1/D$d
        if (M + 1 <= length(D$d))
            D$d[(M + 1):length(D$d)] <- 0
        res <- D$u %*% diag(D$d) %*% D$v
        attr(res, "rank") <- M
        res
    }
    p.table <- pTerms.table <- s.table <- NULL
    covmat <- object$V.beta
    covmat.unscaled <- covmat/object$dispersion
    est.disp <- object$est.disp
    if (!is.null(dispersion)) {
        covmat <- dispersion * covmat.unscaled
        est.disp <- FALSE
    }
    else dispersion <- object$dispersion
    se <- 0
    for (i in 1:length(object$coef$beta)) se[i] <- covmat[i,i]^0.5
    residual.df <- length(object$y) - sum(object$edf)
    alpha <- object$coef$alpha
    delta <- object$coef$delta
    se.alpha <- object$se.alpha
    se.delta <- object$se.delta
    n.delta <- length(se.delta)
    names.delta <- paste("delta", 1:n.delta, sep="")

    if (object$G$nsdf > 0) {
        p.coeff <- object$coef$beta[1:object$G$nsdf]
        p.coeff[length(p.coeff)+(1:(n.delta+1))] <- c(alpha, delta)
        names(p.coeff) <- c(names(object$coef$beta[1:object$G$nsdf]),"alpha",names.delta)
        p.se <- c(se[1:object$G$nsdf], se.alpha, se.delta)
        p.t <- p.coeff/p.se
        if (!est.disp) {
            p.pv <- 2 * pnorm(abs(p.t), lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate",
                "Std. Error", "z value", "Pr(>|z|)"))
        }
        else {
            p.pv <- 2 * pt(abs(p.t), df = residual.df, lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate",
                "Std. Error", "t value", "Pr(>|t|)"))
        }
    }
    else {
        p.coeff <- c(alpha, delta)
        names(p.coeff) <- c("alpha","delta")
        p.se <- c(se.alpha, se.delta)
        p.t <- p.coeff/p.se
        if (!est.disp) {
            p.pv <- 2 * pnorm(abs(p.t), lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate",
                "Std. Error", "z value", "Pr(>|z|)"))
        }
        else {
            p.pv <- 2 * pt(abs(p.t), df = residual.df, lower.tail = FALSE)
            p.table <- cbind(p.coeff, p.se, p.t, p.pv)
            dimnames(p.table) <- list(names(p.coeff), c("Estimate",
                "Std. Error", "t value", "Pr(>|t|)"))
        }
    }
    term.labels <- attr(object$G$pterms, "term.labels")
    nt <- length(term.labels)
    if (nt > 0) {
        np <- length(object$G$assign)
        Vb <- matrix(covmat[1:np, 1:np], np, np)
        bp <- array(object$coef$beta[1:np], np)
        pTerms.pv <- array(0, nt)
        attr(pTerms.pv, "names") <- term.labels
        pTerms.df <- pTerms.chi.sq <- pTerms.pv
        for (i in 1:nt) {
            ind <- object$G$assign == i
            b <- bp[ind]
            V <- Vb[ind, ind]
            if (length(b) == 1) {
                V <- 1/V
                pTerms.df[i] <- nb <- 1
                pTerms.chi.sq[i] <- V * b * b
            }
            else {
                V <- pinv(V, length(b), rank.tol = .Machine$double.eps^0.5)
                pTerms.df[i] <- nb <- attr(V, "rank")
                pTerms.chi.sq[i] <- t(b) %*% V %*% b
            }
            if (!est.disp)
                pTerms.pv[i] <- pchisq(pTerms.chi.sq[i], df = nb,
                  lower.tail = FALSE)
            else pTerms.pv[i] <- pf(pTerms.chi.sq[i]/nb, df1 = nb,
                df2 = residual.df, lower.tail = FALSE)
        }
        if (!est.disp) {
            pTerms.table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)
            dimnames(pTerms.table) <- list(term.labels, c("df",
                "Chi.sq", "p-value"))
        }
        else {
            pTerms.table <- cbind(pTerms.df, pTerms.chi.sq/pTerms.df,
                pTerms.pv)
            dimnames(pTerms.table) <- list(term.labels, c("df",
                "F", "p-value"))
        }
    }
    else {
        pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0, 0)
    }

    m <- length(object$G$smooth)
    df <- edf <- s.pv <- chi.sq <- array(0, m)
    if (m > 0) {
        for (i in 1:m) {
            start <- object$G$smooth[[i]]$first.para
            stop <- object$G$smooth[[i]]$last.para
            V <- covmat[start:stop, start:stop]
            p <- object$coef$beta[start:stop]
            M1 <- object$G$smooth[[i]]$df
            M <- min(M1, ceiling(2 * object$edf.smooth[i]))
            V <- pinv(V, M)
            chi.sq[i] <- t(p) %*% V %*% p
            er <- names(object$coef$beta)[start]
            er <- substring(er, 1, nchar(er) - 2)

            if (object$G$smooth[[i]]$by != "NA") {
                er <- paste(er, ":", object$G$smooth[[i]]$by, sep = "")
            }
            names(chi.sq)[i] <- er
            edf[i] <- sum(object$edf[start:stop])
            df[i] <- attr(V, "rank")
            if (!est.disp)
                s.pv[i] <- pchisq(chi.sq[i], df = df[i], lower.tail = FALSE)
            else s.pv[i] <- pf(chi.sq[i]/df[i], df1 = df[i],
                df2 = residual.df, lower.tail = FALSE)
        }
        if (!est.disp) {
            s.table <- cbind(edf, df, chi.sq, s.pv)
            dimnames(s.table) <- list(names(chi.sq), c("edf",
              "Est.rank", "Chi.sq", "p-value"))
        }

        else {
            s.table <- cbind(edf, df, chi.sq/df, s.pv)
            dimnames(s.table) <- list(names(chi.sq), c("edf",
              "Est.rank", "F", "p-value"))
        }
    }

    nobs <- length(object$y)

    ret <- list(p.coeff = p.coeff, se = se, p.t = p.t, p.pv = p.pv,
        residual.df = residual.df, m = m, chi.sq = chi.sq, s.pv = s.pv,
        family = object$family,
        formula = object$formula, n = nobs,
        edf = edf, dispersion = dispersion,
        cov.unscaled = covmat.unscaled, cov.scaled = covmat,
        p.table = p.table, s.table = s.table)

    class(ret) <- "summary.cozigam"
    ret

}


# test functions

f0 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f1 <- function(x) sin(pi*x)

test <- function(x, z, sx=0.3, sz=0.4)
{ (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
  0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))

}


# zigam(COZIGAM)

zigam <- function(formula, maxiter = 20, conv.crit = 1e-3,
    size = NULL, log.tran = FALSE, family, ...)
{

      if (is.character(family))
          fam <- eval(parse(text = family))
      if (is.function(family))
          fam <- family()
      if (fam$family == "gaussian" | fam$family == "Gamma") {
          zigam.res <- ZIGAM.cts(formula, log.tran = log.tran, family = fam, ...)
          attr(zigam.res, "family.type") <- "continuous"
      }
      else if (fam$family == "poisson" | fam$family == "binomial") {
          zigam.res <- ZIGAM.dis(formula, maxiter, conv.crit, size = size, family = fam, ...)
          attr(zigam.res, "family.type") <- "discrete"
      }
      else stop("family not recognized")
      attr(zigam.res, "constraint") <- "none"

      invisible(zigam.res)

}


# ZIGAM.cts(COZIGAM)

ZIGAM.cts <- function(formula, log.tran = FALSE, family = gaussian(), ...)
{

        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        z <- as.numeric(y!=0)
        n <- length(y)
        if(log.tran==TRUE) y[y!=0] <- log(y[y!=0])

        if (is.character(family))
                family <- eval(parse(text = family))
        if (is.function(family))
                family <- family()
        if(family$family == "gaussian") {
            d.V <- function(mu) rep.int(0, length(mu))
            d.eta.mu <- function(mu) rep.int(0, length(mu))
            loglikfun <- function(y, mu, disp, z, p) {
                sum(z*(log(p)+dnorm(y,mu,sqrt(disp),log=TRUE))+(1-z)*log(1-p))
            }
        }
        else if(family$family == "Gamma") {
            d.V <- function(mu) 2*mu
            d.eta.mu <- function(mu) 2/(mu^3)
            loglikfun <- function(y, mu, disp, z, p) {
                sum(z*(log(p)+ifelse(y==0, 1, dgamma(y,1/disp,scale=mu*disp,log=TRUE)))+(1-z)*log(1-p))
            }
        }
        else
            stop("family not recognized")

        variance <- family$variance
        linkinv <- family$linkinv
        linkfun <- family$linkfun
        mu.eta <- family$mu.eta

        fm1 <- as.formula(sub(split$response,"y",deparse(formula)))
        fm2 <- as.formula(sub(split$response,"z",deparse(formula)))

        G1 <- gam(fm1, subset=(y!=0), fit=FALSE, family=family, ...)
        G2 <- gam(fm2, family=binomial, fit=FALSE, ...)

        fit.gam <- gam(G = G1)
        fit.lr <- gam(G = G2)

        b1 <- fit.gam$coef; b2 <- fit.lr$coef
        mu <- rep.int(0,n)
        mu[y!=0] <- fitted(fit.gam)
        mu.eta.val <- mu.eta(mu)
        p <- fit.lr$fitted
        disp <- fit.gam$sig2

        np1 <- length(b1); np2 <- length(b2)
        sp1 <- fit.gam$sp; sp2 <- fit.lr$sp

        n.smooth <- length(G1$smooth)
        Lambda1 <- matrix(0, np1, np1)
            Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
            for(k in 1:n.smooth) {
                n.S[k] <- length(G1$smooth[[k]]$S)
                if(k==1) {
                    Lam[[k]] <- sp1[k]*G1$S[[k]]
                    if(n.S[k]>1) {
                        for(j in 2:n.S[k]) {
                            Lam[[k]] <- Lam[[k]]+sp1[j]*G1$S[[j]]
                        }
                    }
                }
                else {
                    Lam[[k]] <- sp1[sum(n.S[1:(k-1)])+1]*G1$S[[sum(n.S[1:(k-1)])+1]]
                    if(n.S[k]>1) {
                        for(j in 2:n.S[k]) {
                            Lam[[k]] <- Lam[[k]]+sp1[sum(n.S[1:(k-1)])+j]*G1$S[[sum(n.S[1:(k-1)])+j]]
                        }
                    }
                }
                first <- G1$smooth[[k]]$first.para
                last <- G1$smooth[[k]]$last.para
                Lambda1[first:last, first:last] <- Lam[[k]]
            }
        n.smooth <- length(G2$smooth)
        Lambda2 <- matrix(0, np2, np2)

        Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
        for(k in 1:n.smooth) {
            n.S[k] <- length(G2$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp2[k]*G2$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp2[j]*G2$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp2[sum(n.S[1:(k-1)])+1]*G2$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp2[sum(n.S[1:(k-1)])+j]*G2$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G2$smooth[[k]]$first.para
            last <- G2$smooth[[k]]$last.para
            Lambda2[first:last, first:last] <- Lam[[k]]
        }

        DS1 <- diag(eigen(Lambda1)$values[abs(eigen(Lambda1)$values)>1e-10])
        DS2 <- diag(eigen(Lambda2)$values[abs(eigen(Lambda2)$values)>1e-10])

        X1 <- matrix(0, ncol=ncol(G1$X), nrow=n)
        X1[y!=0,] <- G1$X
        X2 <- G2$X

        rho <- z-p; d.rho <- rep.int(-1, n)
        tau <- (z*(y-mu)*mu.eta.val/disp/variance(mu))
        d.tau <- (z*mu.eta.val^2*( -variance(mu)/mu.eta.val-(y-mu)*(variance(mu)*d.eta.mu(mu)+
          d.V(mu)/mu.eta.val) )/disp/(variance(mu)^2))

        I.beta <- t(X1)%*%(diag(as.vector(-d.tau*mu.eta.val)))%*%X1 + Lambda1
        I.gamma <- t(X2)%*%(diag(as.vector(-d.rho*p*(1-p))))%*%X2 + Lambda2

        I.theta <- matrix(0, ncol=np1+np2, nrow=np1+np2)
        I.theta[1:np1,1:np1] <- I.beta
        I.theta[(np1+1):(np1+np2),(np1+1):(np1+np2)] <- I.gamma

        loglik <- loglikfun(y, mu, disp, z, p)
        ploglik <- loglik - as.numeric(0.5*t(b1)%*%Lambda1%*%b1) - as.numeric(0.5*t(b2)%*%Lambda2%*%b2)

        logE <- 0.5*determinant(DS1)$modulus + 0.5*determinant(DS2)$modulus+ploglik +
          (np1+np2-(nrow(DS1)+nrow(DS2)))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
        attr(logE, "logarithm") <- NULL
        V.theta <- solve(I.theta)
        V.beta <- V.theta[1:np1, 1:np1]
        V.gamma <- V.theta[(np1+1):(np1+np2),(np1+1):(np1+np2)]

        fit.gam$Vp <- V.beta; fit.lr$Vp <- V.gamma

        res <- list(formula=formula, logE=logE, ploglik=ploglik, loglik=loglik, fit.gam=fit.gam, fit.lr=fit.lr,
          mu=mu, p=p, dispersion=disp, V.beta=V.beta, V.gamma=V.gamma, X1=X1, X2=X2, family=family)
        res

}


# ZIGAM.dis(COZIGAM)

ZIGAM.dis <- function(formula, maxiter = 20, conv.crit = 1e-3,
    size = NULL, family = poisson(), ...)
{

        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        n <- length(y)

        if (is.character(family))
            family <- eval(parse(text = family))
        if (is.function(family))
            family <- family()
        if(family$family == "poisson") {
            d.V <- function(mu) rep.int(1, length(mu))
            d.eta.mu <- function(mu) -1/(mu^2)
            d.f0 <- function(mu) -exp(-mu)
            den <- dpois; disp <- 1; est.disp <- FALSE
            size <- rep.int(1, n)
            loglikfun <- function(y, mu, p) {
                e <- as.numeric(y!=0)
                sum((1-e)*log(1-p+p*dpois(y,mu,log=FALSE))+e*(log(p)+dpois(y,mu,log=TRUE)))
            }
        }
        else if(family$family == "binomial") {
            d.V <- function(mu) 1-2*mu
            d.eta.mu <- function(mu) (2*mu-1)/(mu^2*(1-mu)^2)
            d.f0 <- function(mu) -n*(1-mu)^(n-1)
            den <- function(y, mu) dbinom(y*size, size, mu)
            disp <- 1; est.disp <- FALSE
            loglikfun <- function(y, mu, p) {
                e <- as.numeric(y!=0)
                sum((1-e)*log(1-p+p*dbinom(y*size,size,mu,log=FALSE))+e*(log(p)+dbinom(y*size,size,mu,log=TRUE)))
            }
        }
        else stop("family not recognized")

        variance <- family$variance
        linkinv <- family$linkinv
        linkfun <- family$linkfun
        mu.eta <- family$mu.eta

        fm1 <- as.formula(sub(split$response,"y",deparse(formula)))
        fm2 <- as.formula(sub(split$response,"psi",deparse(formula)))

        if(family$family == "binomial") mu <- pmin(pmax(y, 0.01), 0.99)
        else mu <- pmax(y, 0.01)
        p <- rep(0.7, n)
        psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
        norm <- 1; repli <- 0
        while( norm > conv.crit & repli < maxiter) {

            psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
            G1 <- gam(fm1, family=family, fit=FALSE, ...)
            G2 <- gam(fm2, family=quasibinomial, fit=FALSE, ...)
            G1$w <- psi*size
            fit.gam <- gam(G = G1)
            fit.lr <- gam(G = G2)
            b <- coef(fit.gam)
            g <- coef(fit.lr)

            mu.old <- mu; p.old <- p
            mu <- fit.gam$fitted
            p <- fit.lr$fitted
            norm <- max(abs(p-p.old), sum((mu-mu.old)^2))
            repli <- repli + 1
            cat("iteration =", repli, "\t", "norm =", norm, "\n")
        }

        b1 <- fit.gam$coef; b2 <- fit.lr$coef
        np1 <- length(b1); np2 <- length(b2)
        sp1 <- fit.gam$sp; sp2 <- fit.lr$sp
        mu.eta.val <- mu.eta(fit.gam$linear.predictor)

        n.smooth <- length(G1$smooth)
        Lambda1 <- matrix(0, np1, np1)
        Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
        for(k in 1:n.smooth) {
            n.S[k] <- length(G1$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp1[k]*G1$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp1[j]*G1$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp1[sum(n.S[1:(k-1)])+1]*G1$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp1[sum(n.S[1:(k-1)])+j]*G1$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G1$smooth[[k]]$first.para
            last <- G1$smooth[[k]]$last.para
            Lambda1[first:last, first:last] <- Lam[[k]]
        }

        n.smooth <- length(G2$smooth)
        Lambda2 <- matrix(0, np2, np2)
        Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
        for(k in 1:n.smooth) {
            n.S[k] <- length(G2$smooth[[k]]$S)
            if(k==1) {
                Lam[[k]] <- sp2[k]*G2$S[[k]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp2[j]*G2$S[[j]]
                    }
                }
            }
            else {
                Lam[[k]] <- sp2[sum(n.S[1:(k-1)])+1]*G2$S[[sum(n.S[1:(k-1)])+1]]
                if(n.S[k]>1) {
                    for(j in 2:n.S[k]) {
                        Lam[[k]] <- Lam[[k]]+sp2[sum(n.S[1:(k-1)])+j]*G2$S[[sum(n.S[1:(k-1)])+j]]
                    }
                }
            }
            first <- G2$smooth[[k]]$first.para
            last <- G2$smooth[[k]]$last.para
            Lambda2[first:last, first:last] <- Lam[[k]]
        }

        DS1 <- diag(eigen(Lambda1)$values[abs(eigen(Lambda1)$values)>1e-10])
        DS2 <- diag(eigen(Lambda2)$values[abs(eigen(Lambda2)$values)>1e-10])

        X1 <- G1$X; X2 <- G2$X

        loglik <- loglikfun(y, mu, p) # log-likelihood
        ploglik <- loglik - as.numeric(0.5*t(b1)%*%Lambda1%*%b1) -  as.numeric(0.5*t(b2)%*%Lambda2%*%b2)

        # Model selection criterion
        I.theta <- matrix(0, ncol=np1+np2, nrow=np1+np2)
        tau.mu <- -size
        rho.p <- rep.int(-1, n)
        a <- 1-p+p*den(0, mu)
        tau.mu[y==0] <- ( -p*mu.eta.val/disp/variance(mu)/a*(den(0,mu)-size*mu*den(0,mu)*
          (variance(mu)*d.eta.mu(mu)+d.V(mu)/mu.eta.val)*mu.eta.val/variance(mu)+size*(1-p)*mu*d.f0(mu)/a) )[y==0]
        rho.p[y==0] <- ( (1-den(0,mu))/(a^2)*(-a*(1-2*p)-(1-den(0,mu))*p*(1-p)) )[y==0]

        G.tau.mu <- diag(as.vector(mu.eta.val*tau.mu))
        G.rho.p <- diag(as.vector(p*(1-p)*rho.p))

        I.theta[1:np1,1:np1] <- t(X1) %*% G.tau.mu %*% X1 - Lambda1
        I.theta[(np1+1):(np1+np2),(np1+1):(np1+np2)] <- t(X2) %*% G.rho.p %*% X2 - Lambda2
        I.theta <- -I.theta

        logE <- 0.5*determinant(DS1)$modulus + 0.5*determinant(DS2)$modulus+ploglik +
          (np1+np2-(nrow(DS1)+nrow(DS2)))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
        attr(logE, "logarithm") <- NULL

        V.theta <- solve(I.theta)
        V.beta <- V.theta[1:np1, 1:np1]
        V.gamma <- V.theta[(np1+1):(np1+np2),(np1+1):(np1+np2)]

        fit.gam$Vp <- V.beta; fit.lr$Vp <- V.gamma

        res <- list(formula=formula, logE=logE, ploglik=ploglik, loglik=loglik, fit.gam=fit.gam, fit.lr=fit.lr,
            mu=mu, p=p, psi=psi, dispersion=disp, V.beta=V.beta, V.gamma=V.gamma, X1=X1, X2=X2, family=family)
        res

}


