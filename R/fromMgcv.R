gam.fit2 <- 
function (x, y, sp, S = list(), rS = list(), off, H = NULL, weights = rep(1, 
    nobs), start = NULL, etastart = NULL, mustart = NULL, offset = rep(0, 
    nobs), family = gaussian(), control = gam.control(), intercept = TRUE, 
    deriv = TRUE, gamma = 1, scale = 1, pearson = FALSE, printWarn = TRUE, 
    ...) 
{
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    n <- nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    nSp <- length(S)
    if (nSp == 0) 
        deriv <- FALSE
    St <- totalPenalty(S, H, off, sp, ncol(x))
    Sr <- mroot(St)
    z1 <- w1 <- matrix(0, nobs, nSp)
    beta1old <- matrix(0, ncol(x), nSp)
    upe <- list(beta = rep(0, ncol(x)), beta1 = beta1old, trA = 0, 
        trA1 = rep(0, nSp))
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("Invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
        V <- variance(mu)
        if (pearson) {
            alpha1 <- alpha <- sum((y - mu)^2/V)
        }
        else {
            alpha1 <- dev
        }
        trA1 <- trA <- 0
        if (deriv) 
            GCV1 <- UBRE1 <- trA1 <- alpha1 <- rep(0, nSp)
        else GCV1 <- UBRE1 <- trA1 <- alpha1 <- NULL
        GCV <- nobs * alpha/(nobs - gamma * trA)^2
        UBRE <- alpha/nobs - scale + 2 * gamma/n * trA
        scale.est <- alpha/(nobs - trA)
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop("Length of start should equal ", nvars, 
                  " and correspond to initial coefs for ", deparse(xnames))
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
        boundary <- conv <- FALSE
        rV = matrix(0, ncol(x), ncol(x))
        old.pdev <- 0
        for (iter in 1:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("No observations informative at iteration ", 
                  iter)
                break
            }
            mevg <- mu.eta.val[good]
            mug <- mu[good]
            yg <- y[good]
            weg <- weights[good]
            z <- (eta - offset)[good] + (yg - mug)/mevg
            var.mug <- variance(mug)
            w <- sqrt((weg * mevg^2)/var.mug)
            if (deriv && iter > 1) {
                d2g <- family$d2link(mug)
                dV <- family$dvar(mug)
                eta1 <- (x %*% upe$beta1)[good, ]
                z1 <- as.vector((yg - mug) * d2g * mevg) * eta1
                w1 <- as.vector(-0.5 * w^3/weg * (dV/mevg + 2 * 
                  var.mug * d2g)) * eta1
            }
            if (sum(good) < ncol(x)) 
                stop("Not enough informative observations.")
            oo <- .C(C_update_beta, as.double(x[good, ]), as.double(Sr), 
                as.double(unlist(rS)), as.double(sp), as.double(w), 
                as.double(w1), as.double(z), as.double(z1), as.integer(ncol(Sr)), 
                rSncol = as.integer(unlist(lapply(rS, ncol))), 
                m = as.integer(length(rS)), n = as.integer(sum(good)), 
                q = as.integer(ncol(x)), get.trA = as.integer(0), 
                as.integer(deriv), rank.tol = as.double(.Machine$double.eps), 
                beta = as.double(upe$beta), trA = as.double(upe$trA), 
                beta1 = as.double(upe$beta1), trA1 = as.double(upe$trA1), 
                rV = as.double(rV), rank = as.integer(1))
            upe$beta <- oo$beta
            upe$beta1 <- matrix(oo$beta1, oo$q, oo$m)
            if (any(!is.finite(upe$beta))) {
                conv <- FALSE
                warning("Non-finite coefficients at iteration ", 
                  iter)
                break
            }
            start <- upe$beta
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
                warning("Step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  if (deriv) 
                    upe$beta1 <- (upe$beta1 + beta1old)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                warning("Step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  if (deriv) 
                    upe$beta1 <- (upe$beta1 + beta1old)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            pdev <- dev + t(start) %*% St %*% start
            if (iter > 1 && pdev > old.pdev) {
                ii <- 1
                while (pdev - old.pdev > (0.1 + abs(old.pdev)) * 
                  .Machine$double.eps * 10) {
                  if (ii > 200) 
                    stop("inner loop 3; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  if (deriv) 
                    upe$beta1 <- (upe$beta1 + beta1old)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                  pdev <- dev + t(start) %*% St %*% start
                }
            }
            if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                old.pdev <- pdev
                coef <- coefold <- start
                if (deriv) 
                  beta1old <- upe$beta1
            }
        }
        oo <- .C(C_update_beta, as.double(x[good, ]), as.double(Sr), 
            as.double(unlist(rS)), as.double(sp), as.double(w), 
            as.double(w1), as.double(z), as.double(z1), as.integer(ncol(Sr)), 
            rSncol = as.integer(unlist(lapply(rS, ncol))), m = as.integer(length(rS)), 
            n = as.integer(sum(good)), q = as.integer(ncol(x)), 
            get.trA = as.integer(1), as.integer(deriv), rank.tol = as.double(.Machine$double.eps), 
            beta = as.double(upe$beta), trA = as.double(upe$trA), 
            beta1 = as.double(upe$beta1), trA1 = as.double(upe$trA1), 
            rV = as.double(rV), rank = as.integer(1))
        rV <- matrix(oo$rV, ncol(x), ncol(x))
        upe$beta <- oo$beta
        upe$trA <- oo$trA
        upe$beta1 <- matrix(oo$beta1, oo$q, oo$m)
        upe$trA1 <- oo$trA1
        V <- variance(mug)
        if (pearson) 
            alpha1 <- alpha <- sum(weights[good] * (yg - mug)^2/V)
        else {
            dev <- sum(dev.resids(y, mu, weights))
            alpha1 <- alpha <- dev
        }
        trA1 <- trA <- upe$trA
        GCV <- nobs * alpha/(nobs - gamma * trA)^2
        UBRE <- alpha/nobs + 2 * gamma * trA * scale/nobs - scale
        scale.est <- alpha/(length(mug) - trA)
        if (deriv) {
            mu.eta.val <- mu.eta(eta[good])
            if (pearson) {
                d2g <- family$d2link(mug)
                dV <- family$dvar(mug)
                eta1 <- (x %*% upe$beta1)[good, ]
                temp <- as.vector(-dV/V^2 * mu.eta.val * (yg - 
                  mug)^2 - 2/V * (yg - mug) * mu.eta.val) * eta1
                temp <- as.matrix(temp * as.vector(weights[good]))
                alpha1 <- colSums(temp)
            }
            else {
                temp <- weights[good] * (yg - mug) * mu.eta.val/V
                temp <- t(temp) %*% x[good, ]
                alpha1 <- -2 * as.numeric(temp %*% upe$beta1)
            }
            trA1 <- upe$trA1
            GCV1 <- nobs * alpha1/(nobs - gamma * trA)^2 + 2 * 
                nobs * alpha * gamma * trA1/(nobs - gamma * trA)^3
            UBRE1 <- alpha1/nobs + 2 * gamma * scale/nobs * trA1
        }
        else UBRE1 <- GCV1 <- NULL
        if (!conv && printWarn) 
            warning("Algorithm did not converge")
        if (printWarn && boundary) 
            warning("Algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (printWarn && family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (printWarn && family$family == "poisson") {
            if (any(mu < eps)) 
                warning("fitted rates numerically 0 occurred")
        }
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
        names(coef) <- xnames
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    aic.model <- aic(y, n, mu, weights, dev)
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        family = family, linear.predictors = eta, deviance = dev, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.null = nulldf, y = y, converged = conv, boundary = boundary, 
        alpha = alpha, alpha1 = alpha1, trA = trA, trA1 = trA1, pdev = pdev,
        GCV = GCV, GCV1 = GCV1, UBRE = UBRE, UBRE1 = UBRE1, rV = rV, 
        scale.est = scale.est, aic = aic.model, rank = oo$rank)
}

gam.fit3 <- 
function (x, y, sp, S = list(), rS = list(), off, H = NULL, weights = rep(1, 
    nobs), start = NULL, etastart = NULL, mustart = NULL, offset = rep(0, 
    nobs), family = gaussian(), control = gam.control(), intercept = TRUE, 
    deriv = 2, use.svd = TRUE, gamma = 1, scale = 1, printWarn = TRUE, 
    ...) 
{
    scale <- abs(scale)
    if (!deriv %in% c(0, 1, 2)) 
        stop("unsupported order of differentiation requested of gam.fit3")
    x <- as.matrix(x)
    iter <- 0
    coef <- rep(0, ncol(x))
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    n <- nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    nSp <- length(S)
    if (nSp == 0) 
        deriv <- FALSE
    St <- totalPenalty(S, H, off, sp, ncol(x))
    Sr <- mroot(St)
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("Invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
        V <- variance(mu)
        alpha <- dev
        trA2 <- trA1 <- trA <- 0
        if (deriv) 
            GCV2 <- GCV1 <- UBRE2 <- UBRE1 <- trA1 <- rep(0, 
                nSp)
        else GCV2 <- GCV1 <- UBRE2 <- UBRE1 <- trA2 <- trA1 <- NULL
        GCV <- nobs * alpha/(nobs - gamma * trA)^2
        UBRE <- alpha/nobs - scale + 2 * gamma/n * trA
        scale.est <- alpha/(nobs - trA)
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop("Length of start should equal ", nvars, 
                  " and correspond to initial coefs for ", deparse(xnames))
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        etaold <- eta
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
        boundary <- conv <- FALSE
        rV = matrix(0, ncol(x), ncol(x))
        old.pdev <- 0
        for (iter in 1:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("No observations informative at iteration ", 
                  iter)
                break
            }
            mevg <- mu.eta.val[good]
            mug <- mu[good]
            yg <- y[good]
            weg <- weights[good]
            z <- (eta - offset)[good] + (yg - mug)/mevg
            var.mug <- variance(mug)
            w <- sqrt((weg * mevg^2)/var.mug)
            if (sum(good) < ncol(x)) 
                stop("Not enough informative observations.")
            oo <- .C(C_pls_fit, y = as.double(z), as.double(x[good, 
                ]), as.double(w), as.double(Sr), as.integer(sum(good)), 
                as.integer(ncol(x)), as.integer(ncol(Sr)), eta = as.double(z), 
                penalty = as.double(1), as.double(.Machine$double.eps * 
                  100))
            start <- oo$y[1:ncol(x)]
            penalty <- oo$penalty
            eta <- drop(x %*% start)
            if (any(!is.finite(start))) {
                conv <- FALSE
                warning("Non-finite coefficients at iteration ", 
                  iter)
                break
            }
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found:please supply starting values", 
                    call. = FALSE)
                warning("Step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                warning("Step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            pdev <- dev + penalty
            if (control$trace) 
                cat("penalized deviance =", pdev, "\n")
            div.thresh <- 10 * (0.1 + abs(old.pdev)) * .Machine$double.eps^0.5
            if (iter > 1 && (pdev - old.pdev > div.thresh)) {
                ii <- 1
                while (pdev - old.pdev > div.thresh) {
                  if (ii > 200) 
                    stop("inner loop 3; can't correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- (eta + etaold)/2
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                  pdev <- dev + t(start) %*% St %*% start
                  if (control$trace) 
                    cat("Step halved: new penalized deviance =", 
                      pdev, "\n")
                }
            }
            if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) {
                if (max(abs(start - coefold)) > control$epsilon * 
                  max(abs(start + coefold))/2) {
                  old.pdev <- pdev
                  coef <- coefold <- start
                  etaold <- eta
                }
                else {
                  conv <- TRUE
                  coef <- start
                  break
                }
            }
            else {
                old.pdev <- pdev
                coef <- coefold <- start
                etaold <- eta
            }
        }
        dev <- sum(dev.resids(y, mu, weights))
        good <- weights > 0
        varmu <- variance(mu)[good]
        if (any(is.na(varmu))) 
            stop("NAs in V(mu)")
        if (any(varmu == 0)) 
            stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good]))) 
            stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0)
        mevg <- mu.eta.val[good]
        mug <- mu[good]
        yg <- y[good]
        weg <- weights[good]
        etag <- eta[good]
        z <- (eta - offset)[good] + (yg - mug)/mevg
        var.mug <- variance(mug)
        w <- sqrt((weg * mevg^2)/var.mug)
        g1 <- 1/mevg
        g2 <- family$d2link(mug)
        g3 <- family$d3link(mug)
        V <- family$variance(mug)
        V1 <- family$dvar(mug)
        V2 <- family$d2var(mug)
        P1 <- D1 <- array(0, nSp)
        P2 <- D2 <- matrix(0, nSp, nSp)
        trA1 <- array(0, nSp)
        trA2 <- matrix(0, nSp, nSp)
        rV = matrix(0, ncol(x), ncol(x))
        dum <- 1
        if (control$trace) 
            cat("calling gdi...")
        oo <- .C(C_gdi, X = as.double(x[good, ]), E = as.double(Sr), 
            rS = as.double(unlist(rS)), sp = as.double(exp(sp)), 
            z = as.double(z), w = as.double(w), mu = as.double(mug), 
            eta = as.double(etag), y = as.double(yg), p.weights = as.double(weg), 
            g1 = as.double(g1), g2 = as.double(g2), g3 = as.double(g3), 
            V0 = as.double(V), V1 = as.double(V1), V2 = as.double(V2), 
            beta = as.double(coef), D1 = as.double(D1), D2 = as.double(D2), 
            P = as.double(dum), P1 = as.double(P1), P2 = as.double(P2), 
            trA = as.double(dum), trA1 = as.double(trA1), trA2 = as.double(trA2), 
            rV = as.double(rV), rank.tol = as.double(.Machine$double.eps * 
                100), conv.tol = as.double(control$epsilon), 
            rank.est = as.integer(1), n = as.integer(length(z)), 
            p = as.integer(ncol(x)), M = as.integer(nSp), Encol = as.integer(ncol(Sr)), 
            rSncol = as.integer(unlist(lapply(rS, ncol))), deriv = as.integer(deriv), 
            use.svd = as.integer(use.svd))
        if (control$trace) 
            cat("done!\n")
        rV <- matrix(oo$rV, ncol(x), ncol(x))
        coef <- oo$beta
        trA <- oo$trA
        P <- oo$P
        delta <- nobs - gamma * trA
        delta.2 <- delta * delta
        GCV <- nobs * dev/delta.2
        GACV <- dev/nobs + P * 2 * gamma * trA/(delta * nobs)
        UBRE <- dev/nobs - 2 * delta * scale/nobs + scale
        scale.est <- dev/(nobs - trA)
        if (deriv) {
            trA1 <- oo$trA1
            D1 <- oo$D1
            P1 <- oo$P1
            delta.3 <- delta * delta.2
            GCV1 <- nobs * D1/delta.2 + 2 * nobs * dev * trA1 * 
                gamma/delta.3
            GACV1 <- D1/nobs + 2 * P/delta.2 * trA1 + 2 * gamma * 
                trA * P1/(delta * nobs)
            UBRE1 <- D1/nobs + gamma * trA1 * 2 * scale/nobs
            if (deriv == 2) {
                trA2 <- matrix(oo$trA2, nSp, nSp)
                D2 <- matrix(oo$D2, nSp, nSp)
                P2 <- matrix(oo$P2, nSp, nSp)
                GCV2 <- outer(trA1, D1)
                GCV2 <- (GCV2 + t(GCV2)) * gamma * 2 * nobs/delta.3 + 
                  6 * nobs * dev * outer(trA1, trA1) * gamma * 
                    gamma/(delta.2 * delta.2) + nobs * D2/delta.2 + 
                  2 * nobs * dev * gamma * trA2/delta.3
                GACV2 <- D2/nobs + outer(trA1, trA1) * 4 * P/(delta.3) + 
                  2 * P * trA2/delta.2 + 2 * outer(trA1, P1)/delta.2 + 
                  2 * outer(P1, trA1) * (1/(delta * nobs) + trA/(nobs * 
                    delta.2)) + 2 * trA * P2/(delta * nobs)
                GACV2 <- (GACV2 + t(GACV2)) * 0.5
                UBRE2 <- D2/nobs + 2 * gamma * trA2 * scale/nobs
            }
            else {
                GACV2 <- P2 <- trA2 <- D2 <- UBRE2 <- GCV2 <- NULL
            }
        }
        else {
            GACV1 <- P1 <- P2 <- trA1 <- trA2 <- D1 <- D2 <- GACV2 <- UBRE2 <- GCV2 <- UBRE1 <- GCV1 <- NULL
        }
        if (!conv && printWarn) 
            warning("Algorithm did not converge")
        if (printWarn && boundary) 
            warning("Algorithm stopped at boundary value")
        eps <- 10 * .Machine$double.eps
        if (printWarn && family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("fitted probabilities numerically 0 or 1 occurred")
        }
        if (printWarn && family$family == "poisson") {
            if (any(mu < eps)) 
                warning("fitted rates numerically 0 occurred")
        }
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
        names(coef) <- xnames
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    aic.model <- aic(y, n, mu, weights, dev)
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        family = family, linear.predictors = eta, deviance = dev, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.null = nulldf, y = y, converged = conv, boundary = boundary, 
        D1 = D1, D2 = D2, P = P, P1 = P1, P2 = P2, trA = trA, pdev = pdev,
        trA1 = trA1, trA2 = trA2, GCV = GCV, GCV1 = GCV1, GCV2 = GCV2, 
        GACV = GACV, GACV1 = GACV1, GACV2 = GACV2, UBRE = UBRE, 
        UBRE1 = UBRE1, UBRE2 = UBRE2, rV = rV, scale.est = scale.est, 
        aic = aic.model, rank = oo$rank.est)
}
