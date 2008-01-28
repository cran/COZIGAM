## cozigam.R
## Hai Liu
## Jan 2008


# cozigam(COZIGAM)

cozigam <- function(formula, maxiter = 30, conv.crit.in = 1e-5, 
    conv.crit.out = 1e-4, size = NULL, log.tran = FALSE, family, ...)
{

      if (is.character(family)) 
          fam <- eval(parse(text = family))
      if (is.function(family)) 
          fam <- family()
      if (fam$family == "gaussian" | fam$family == "Gamma") {
          cozigam.res <- COZIGAM.cts(formula, maxiter, conv.crit.in, 
              conv.crit.out, log.tran = log.tran, family = fam, ...)
          attr(cozigam.res, "family.type") <- "continuous"
      }
      else if (fam$family == "poisson" | fam$family == "binomial") {
          cozigam.res <- COZIGAM.dis(formula, maxiter, conv.crit.in, 
              conv.crit.out, size = size, family = fam, ...)
          attr(cozigam.res, "family.type") <- "discrete"
      }
      else stop("family not recognized")

      invisible(cozigam.res)

}


# COZIGAM.cts(COZIGAM)

COZIGAM.cts <- function(formula, maxiter = 20, conv.crit.in = 1e-5, 
    conv.crit.out = 1e-3, log.tran = FALSE, family = gaussian(), ...)
{
        require(mgcv)
        split <- interpret.gam(formula)
        y <- eval(parse(text=split$response))
        n <- length(y)
        z <- (y!=0)
        if(log.tran) y[z] <- log(y[z])
        
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

        ## estimating dispersion using non-zero data...
        fit.nz <- gam(fm1, subset=z, family=family, ...)
        disp <- fit.nz$sig2; est.disp <- TRUE
        
        ## setting initial values...
        alpha <- 0; delta <- 1

        mu <- pmax(y, 0.01)
        p <- rep(0.7, n)
        eta1 <- linkfun(mu); eta2 <- delta*eta1
         
        ## getting design matrices & penalty matrix
        G1 <- gam(fm1, fit=FALSE, family=family, ...)
        G2 <- gam(fm2, fit=FALSE, family=quasi(variance="mu(1-mu)",link="logit"), ...)

        rep.out <- 1; norm.out <- 1
        while(norm.out>conv.crit.out & rep.out<maxiter) {

            b.old <- b <- rep(0, ncol(G1$X))
        
            ## loop: approx logLik by WSS
            ## updating beta
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
                X.c <- rbind(G1$X[good,], delta*G2$X[good,])

                mgfit <- magic(y.c, X.c, G1$sp, G1$S, G1$off, G1$rank, C=G1$C, w=w.c, gcv=FALSE)

                b.old <- b
                b <- mgfit$b
                sp <- mgfit$sp
                eta1 <- G1$X%*%b; eta2 <- delta*eta1
                mu <- linkinv(eta1)
                p <- .Call("logit_linkinv", eta2+alpha, PACKAGE = "stats")
                norm.in <- sum((b-b.old)^2)
                rep.in <- rep.in + 1
                #cat("iteration =", rep.in, "\t", "norm =", norm.in, "\n")
            }   
 
            ## updating alpha & delta
            ## fit a generalized linear model given eta1
            alpha.old <- alpha
            delta.old <- delta
       
            fit.glm <- glm(z ~ eta1, family=binomial)
            alpha <- fit.glm$coef[1]
            delta <- fit.glm$coef[2]

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

        ## observed information...
        n.smooth <- length(G1$smooth)
        np <- length(b)
            Lambda <- matrix(0, np, np)


        Lam <- list(); n.S <- numeric(n.smooth)
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

        FM <- solve(t(X.c)%*%diag(w.c^2)%*%X.c+Lambda)%*%t(X.c)%*%diag(w.c^2)%*%X.c

        edf.smooth <- numeric(n.smooth)
        edf <- diag(FM)
        for(k in 1:n.smooth) {
            first <- G1$smooth[[k]]$first.para
            last <- G1$smooth[[k]]$last.para
            edf.smooth[k] <- sum(edf[first:last])
        }

        if (G1$nsdf > 0) 
            term.names <- colnames(G1$X)[1:G1$nsdf]
        else term.names <- array("", 0)
        for (i in 1:n.smooth) {
            k <- 1
            for (j in G1$smooth[[i]]$first.para:G1$smooth[[i]]$last.para) {
                term.names[j] <- paste(G1$smooth[[i]]$label, ".", 
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

        I.alpha <- sum(-d.rho*p*(1-p))
        I.delta <- sum(-p*(1-p)*eta1^2*d.rho)
        I.alpha.delta <- sum(-p*(1-p)*eta1*d.rho)
        I.beta <- (t(G1$X)%*%(diag(as.vector(-d.tau*mu.eta.val)))%*%G1$X) +
                  delta^2*t(G1$X)%*%(diag(as.vector(-d.rho*p*(1-p))))%*%G1$X + Lambda
        I.alpha.beta <- delta*t(G1$X)%*%(-d.rho*p*(1-p)) 
        I.delta.beta <- delta*t(G1$X)%*%(-d.rho*p*(1-p)*eta1) 
 
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

        loglik <- loglikfun(y, mu, disp, z, p)
        ploglik <- loglik - as.numeric(0.5*t(b)%*%Lambda%*%b)
                
        res <- list(coefficient=theta, V.theta=V.theta, converged=converged, 
              mu=mu, linear.predictor=eta1, dispersion=disp, formula=formula, 
              fit.nonzero=fit.nz, score=mgfit$score, G=G1, beta=b, V.beta=V.b, 
              y=y, p=p, loglik=loglik, ploglik=ploglik, family=family, edf=edf,
              est.disp=est.disp, edf.smooth=edf.smooth)
        class(res) <- c("cozigam", "gam", "glm", "lm")
        res
}


# COZIGAM.dis(COZIGAM)

COZIGAM.dis <- function(formula, maxiter = 30, conv.crit.in = 1e-5, 
    conv.crit.out = 1e-4, size = NULL, family = poisson(), ...)
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

        fm1 <- as.formula(sub(split$response,"y",deparse(formula)))
        fm2 <- as.formula(sub(split$response,"psi",deparse(formula)))

        ## setting initial values...
        alpha <- 0; delta <- 1

        if(family$family == "binomial") mu <- pmin(pmax(y, 0.01), 0.99)
        else mu <- pmax(y, 0.01)
        p <- rep(0.6, n)
        eta1 <- linkfun(mu); eta2 <- delta*eta1
         
        ## getting design matrices & penalty matrix
        psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
        G1 <- gam(fm1, fit=FALSE, family=family, ...)
        G2 <- gam(fm2, fit=FALSE, family=quasi(variance="mu(1-mu)",link="logit"), ...)

        ## outer loop: EM algorithm
        rep.out <- 1; norm.out <- 1
        while(norm.out>conv.crit.out & rep.out<maxiter) {

            psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
            b.old <- b <- rep(0, ncol(G1$X))
        
            ## inner loop: approx logLik by WSS
            ## updating beta
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
                X.c <- rbind(G1$X[good,], delta*G2$X[good,])

                mgfit <- magic(y.c, X.c, G1$sp, G1$S, G1$off, G1$rank, C=G1$C, w=w.c, gcv=FALSE)

                b.old <- b
                b <- mgfit$b
                sp <- mgfit$sp
                eta1 <- G1$X%*%b; eta2 <- delta*eta1
                mu <- linkinv(eta1)
                p <- .Call("logit_linkinv", eta2+alpha, PACKAGE = "stats")
                norm.in <- sum((b-b.old)^2)
                rep.in <- rep.in + 1
            
            }   
 
            ## updating alpha & delta
            ## fit a generalized linear model given eta1
            alpha.old <- alpha
            delta.old <- delta
       
            fit.glm <- glm(psi ~ eta1, family=quasibinomial)
            alpha <- fit.glm$coef[1]
            delta <- fit.glm$coef[2]

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

        ## observed information...
        n.smooth <- length(G1$smooth)
        np <- length(b)
        Lambda <- matrix(0, np, np)
        Lam <- list(); n.S <- numeric(n.smooth)
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

        FM <- solve(t(X.c)%*%diag(w.c^2)%*%X.c+Lambda)%*%t(X.c)%*%diag(w.c^2)%*%X.c

        edf.smooth <- numeric(n.smooth)
        edf <- diag(FM)
        for(k in 1:n.smooth) {
            first <- G1$smooth[[k]]$first.para
            last <- G1$smooth[[k]]$last.para
            edf.smooth[k] <- sum(edf[first:last])
        }

        if (G1$nsdf > 0) 
            term.names <- colnames(G1$X)[1:G1$nsdf]
        else term.names <- array("", 0)
        for (i in 1:n.smooth) {
            k <- 1
            for (j in G1$smooth[[i]]$first.para:G1$smooth[[i]]$last.para) {
                term.names[j] <- paste(G1$smooth[[i]]$label, ".", 
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
        I.beta <- t(G1$X)%*%(diag(as.vector(-Ed.tau*mu.eta.val-psi*(1-psi)*size^2*(y-mu)^2*mu.eta.val^2/(disp^2)/(variance(mu)^2))))%*%G1$X +
                  delta^2*t(G1$X)%*%(diag(as.vector(-Ed.rho*p*(1-p)-psi*(1-psi))))%*%G1$X -
                  2*delta*t(G1$X)%*%(diag(as.vector(psi*(1-psi)*size*(y-mu)*mu.eta.val/disp/variance(mu))))%*%G1$X + Lambda
        I.alpha.beta <- -delta*t(G1$X)%*%(Ed.rho*p*(1-p)+psi*(1-psi)) - 
                        t(G1$X)%*%(psi*(1-psi)*size*(y-mu)*mu.eta.val/disp/variance(mu)) 
        I.delta.beta <- -delta*t(G1$X)%*%(eta1*(Ed.rho*p*(1-p)+psi*(1-psi))) - 
                        t(G1$X)%*%(psi-p+eta1*psi*(1-psi)*size*(y-mu)*mu.eta.val/disp/variance(mu)) 
 
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

        loglik <- loglikfun(y, mu, p)
        ploglik <- loglik - as.numeric(0.5*t(b)%*%Lambda%*%b)

        res <- list(coefficient=theta, V.theta=V.theta, converged=converged, 
              mu=mu, linear.predictor=eta1, dispersion=disp, formula=formula, 
              score=mgfit$score, G=G1, beta=b, V.beta=V.b, y=y*size, psi=psi, 
              p=p, family=family, edf=edf, loglik=loglik, ploglik=ploglik,
              est.disp=est.disp, edf.smooth=edf.smooth)
              
        class(res) <- c("cozigam", "gam", "glm", "lm")
        res

}


# plot.cozigam(COZIGAM)

plot.cozigam <- function(x, plot.2d = "contour", too.far = 0, 
    n.1d = 100, n.2d = 30, theta = 30, phi = 30, select = NULL, image.col = "topo", 
    persp.col = "lightblue", contour.col = "red", n.Col = 100, shade.ci = FALSE,
    shade.col = "gray80", Rug = TRUE, ...)
{
      if (!inherits(x, "cozigam")) 
          stop("use only with \"COZIGAM\" objects")
      require(mgcv)
      G <- x$G
      n.smooth <- length(G$smooth)
      b <- x$beta; V.b <- x$V.beta

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

                if (shade.ci) {
                    plot(fit~xx, type="l", col="black", xlab=colnames(pr), 
                        ylim=c(min(ci.lower),max(ci.upper)), ...)
                    polygon(c(xx, xx[n.1d:1], xx[1]), c(ci.upper, ci.lower[n.1d:1], ci.upper[1]), 
                        col = shade.col, border = NA)
                    lines(xx, fit)
                    if(Rug) {
                        rug(raw, col="black")
                    }
                }
                else {
                    plot(fit~xx, type="l", col="black", xlab=colnames(pr), 
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
                if(too.far>0) {
                      exclude <- exclude.too.far(gx, gz, raw$xx, raw$zz, dist = too.far)
                      fit[exclude] <- NA
                }
                if(plot.2d=="contour") {
                      if(is.null(image.col)) contour(xs, zs, fit, col=contour.col)
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

                            image(xs, zs, fit, col = Col, xlab = G$smooth[[k]]$term[1], 
                                  ylab = G$smooth[[k]]$term[2],...)
                            contour(xs, zs, fit, col=contour.col, add=TRUE)
                      }
                      if (Rug) {
                            if (is.null(list(...)[["pch"]])) 
                                  points(raw$xx, raw$zz, pch = ".", ...)
                            else points(raw$xx, raw$zz, ...)
                      }
                }
                else if(plot.2d=="persp") 
                      persp(xs, zs, fit, col=persp.col, theta=theta, phi=phi, xlab = G$smooth[[k]]$term[1],
                          ylab = G$smooth[[k]]$term[2], ...)
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
      b <- object$beta; V.b <- object$V.beta
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
      
      eta.fit <- apply(fit, 1, sum) + b[1]
      mu <- linkinv(eta.fit)
      eta.p <- alpha + delta*eta.fit
      p <- .Call("logit_linkinv", eta.p, PACKAGE = "stats")
      
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

# test functions

f0 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f1 <- function(x) sin(pi*x)

testfn <- function(x, z, sx=0.3, sz=0.4) 
{ (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
  0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))

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
    se.alpha <- sqrt(object$V.theta[1,1])
    se.delta <- sqrt(object$V.theta[2,2])

 
    if (object$G$nsdf > 0) {
        p.coeff <- object$coef$beta[1:object$G$nsdf]
        p.coeff[length(p.coeff)+(1:2)] <- c(alpha, delta)
        names(p.coeff) <- c(names(object$coef$beta[1:object$G$nsdf]),"alpha","delta")
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


# print.cozigam(COZIGAM)

print.cozigam <- function (x, ...) 
{
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    alpha <- x$coef$alpha
    delta <- x$coef$delta
    se.alpha <- sqrt(x$V.theta[1,1])
    se.delta <- sqrt(x$V.theta[2,2])
    cat("\nCoefficients of Constraint","\n")
    cat(" alpha =", round(alpha,5), "(", round(se.alpha,5), ")", "\n")
    cat(" delta =", round(delta,5), "(", round(se.delta,5), ")", "\n")
    n.smooth <- length(x$G$smooth)
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

