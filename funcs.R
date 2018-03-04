do.gw.lass <- function(X,Y,W.i,i) {
		x <- X * W.i
		cvfit = cv.glmnet(x, Y)
		coef.fit <- as.matrix(coef(cvfit, s = "lambda.min"))
		r2 <- cvfit$glmnet.fit$dev.ratio[
			which(cvfit$glmnet.fit$lambda == cvfit$lambda.min)]
		yhat <- predict(cvfit,newx=x,s="lambda.min")
		return(list(coef.fit, r2, yhat[i]))
}

gw.lass.cv <- function (bw, X, Y, kernel = "bisquare", adaptive = FALSE, dp.locat, 
    p = 2, theta = 0, longlat = F, dMat, verbose = T) 
{
    dp.n <- length(dp.locat[, 1])
    if (is.null(dMat)) 
        DM.given <- F
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    CV <- numeric(dp.n)
    for (i in 1:dp.n) {
        if (DM.given) 
            dist.vi <- dMat[, i]
        else {
            dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                p = p, theta = theta, longlat = longlat)
        }
        W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
        W.i[i] <- 0
        # gw.resi <- try(gw_reg(X, Y, W.i, FALSE, i))
        #print(head(X))
        gw.resi <- try(do.gw.lass(X[,-1],Y,W.i))
        
        #cat(head(X))
  
        if (!inherits(gw.resi, "try-error")) {
            yhat.noi <- X[i, ] %*% gw.resi[[1]]
            CV[i] <- Y[i] - yhat.noi
        }
        else {
            CV[i] <- Inf
            break
        }
    }
    if (!any(is.infinite(CV))) 
        CV.score <- t(CV) %*% CV
    else {
        CV.score <- Inf
    }
    cat("Adaptive bandwidth:", bw, "CV score:", CV.score, "\n")
	#df.all <<- rbind(df.all, c(bw, CV.score))
    
    CV.score
}

bw.lass <- function (formula, data, approach = "CV", kernel = "bisquare", 
    adaptive = FALSE, p = 2, theta = 0, longlat = F, dMat) 
{
    if (is(data, "Spatial")) {
        dp.locat <- coordinates(data)
        data <- as(data, "data.frame")
    }
    else {
        stop("Given regression data must be Spatial*DataFrame")
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    dp.n <- nrow(data)
    if (dp.n > 1500) {
        cat("Take a cup of tea and have a break, it will take a few minutes.\n")
        cat("          -----A kind suggestion from GWmodel development group\n")
    }
    if (missing(dMat)) {
        DM.given <- F
        if (dp.n + dp.n <= 10000) {
            dMat <- gw.dist(dp.locat = dp.locat, rp.locat = dp.locat, 
                p = p, theta = theta, longlat = longlat)
            DM.given <- T
        }
    }
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    if (adaptive) {
        upper <- dp.n
        lower <- 20
    }
    else {
        if (DM.given) {
            upper <- range(dMat)[2]
            lower <- upper/5000
        }
        else {
            dMat <- NULL
            if (p == 2) {
                b.box <- bbox(dp.locat)
                upper <- sqrt((b.box[1, 2] - b.box[1, 1])^2 + 
                  (b.box[2, 2] - b.box[2, 1])^2)
                lower <- upper/5000
            }
            else {
                upper <- 0
                for (i in 1:dp.n) {
                  dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                    p = p, theta = theta, longlat = longlat)
                  upper <- max(upper, range(dist.vi)[2])
                }
                lower <- upper/5000
            }
        }
    }
    bw <- NA
    #df.i <- matrix(nrow = 0, ncol = 2)
    if (approach == "cv" || approach == "CV") 
        bw <- gold(gw.lass.cv, lower, upper, adapt.bw = adaptive, 
            x, y, kernel, adaptive, dp.locat, p, theta, longlat, 
            dMat)
    else if (approach == "aic" || approach == "AIC" || approach == 
        "AICc") 
        bw <- gold(gwr.aic, lower, upper, adapt.bw = adaptive, 
            x, y, kernel, adaptive, dp.locat, p, theta, longlat, 
            dMat)
    as.vector(bw)
}


gwr.lass <- function (formula, data, regression.points, bw, kernel = "bisquare", 
    adaptive = FALSE, p = 2, theta = 0, longlat = F, dMat, F123.test = F, 
    cv = T, W.vect = NULL) 
{
    timings <- list()
    timings[["start"]] <- Sys.time()
    this.call <- match.call()
    p4s <- as.character(NA)
    if (missing(regression.points)) {
        rp.given <- FALSE
        regression.points <- data
        rp.locat <- coordinates(data)
        hatmatrix <- T
    }
    else {
        rp.given <- TRUE
        hatmatrix <- F
        if (is(regression.points, "Spatial")) {
            rp.locat <- coordinates(regression.points)
        }
        else if (is.numeric(regression.points) && dim(regression.points)[2] == 
            2) 
            rp.locat <- regression.points
        else {
            warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
            rp.locat <- dp.locat
        }
    }
    griddedObj <- F
    if (is(regression.points, "Spatial")) {
        if (is(regression.points, "SpatialPolygonsDataFrame")) 
            polygons <- polygons(regression.points)
        else griddedObj <- gridded(regression.points)
    }
    if (is(data, "Spatial")) {
        p4s <- proj4string(data)
        dp.locat <- coordinates(data)
        data <- as(data, "data.frame")
    }
    else {
        stop("Given regression data must be Spatial*DataFrame")
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    var.n <- ncol(x)
    rp.n <- nrow(rp.locat)
    dp.n <- nrow(data)
    betas <- matrix(nrow = rp.n, ncol = var.n)
    betas.SE <- matrix(nrow = rp.n, ncol = var.n)
    betas.TV <- matrix(nrow = rp.n, ncol = var.n)
    S <- matrix(nrow = dp.n, ncol = dp.n)
    idx1 <- match("(Intercept)", colnames(x))
    if (!is.na(idx1)) 
        colnames(x)[idx1] <- "Intercept"
    colnames(betas) <- colnames(x)
    lms <- lm(formula, data = data)
    lms$x <- x
    lms$y <- y
    gTSS <- c(cov.wt(matrix(y, ncol = 1), wt = rep(as.numeric(1), 
        dp.n), method = "ML")$cov * dp.n)
    if (missing(dMat)) {
        DM.given <- F
        DM1.given <- F
        if (dp.n + rp.n <= 10000) {
            dMat <- gw.dist(dp.locat = dp.locat, rp.locat = rp.locat, 
                p = p, theta = theta, longlat = longlat)
            DM.given <- T
        }
    }
    else {
        DM.given <- T
        DM1.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != rp.n) 
            stop("Dimensions of dMat are not correct")
    }
    W <- matrix(nrow = dp.n, ncol = rp.n)
    lass.yhat <- vector(length = (dp.n))
    for (i in 1:rp.n) {
        if (DM.given) 
            dist.vi <- dMat[, i]
        else {
            if (rp.given) 
                dist.vi <- gw.dist(dp.locat, rp.locat, focus = i, 
                  p, theta, longlat)
            else dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                p = p, theta = theta, longlat = longlat)
        }
        W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
        W[, i] <- W.i
        if (!is.null(W.vect)) 
            W.i <- W.i * W.vect
        #gwsi <- gw_reg(x, y, W.i, hatmatrix, i)
 		#cat(colnames(x))
 		gwsi <- do.gw.lass(x[,-1],y,W.i,i)
 
        #gwsi.out <<- gwsi
        #cat(dim(betas))
        betas[i, ] <- (gwsi[[1]])
        lass.yhat[i] <- gwsi[[3]]
        
        if (hatmatrix) {
            S[i, ] <- gwsi[[2]]
            Ci <- gwsi[[3]]
            betas.SE[i, ] <- diag(Ci %*% t(Ci))
        }
    }
    GW.diagnostic <- NA
    Ftests <- list()
    if (hatmatrix) {
        diags <- gwr_diag(y, x, betas, S)
        tr.S <- diags[8]
        tr.StS <- diags[9]
        Q <- t(diag(dp.n) - S) %*% (diag(dp.n) - S)
        RSS.gw <- diags[5]
        yhat <- gw.fitted(x, betas)
        residual <- y - yhat
        CV <- numeric(dp.n)
        local.R2 <- numeric(dp.n)
        if (cv) 
            CV <- gwr.cv.contrib(bw, x, y, kernel, adaptive, 
                dp.locat, p, theta, longlat, dMat)
        sigma.hat1 <- RSS.gw/(dp.n - 2 * tr.S + tr.StS)
        Stud_residual <- residual
        q.diag <- diag(Q)
        for (i in 1:dp.n) {
            Stud_residual[i] <- residual[i]/sqrt(sigma.hat1 * 
                q.diag[i])
            betas.SE[i, ] <- sqrt(sigma.hat1 * betas.SE[i, ])
            betas.TV[i, ] <- betas[i, ]/betas.SE[i, ]
            W.i <- W[, i]
            if (!is.null(W.vect)) 
                W.i <- W.i * W.vect
            TSSw <- t((y - mean(y)) * W.i) %*% (y - mean(y))
            RSSw <- t((y - yhat) * W.i) %*% (y - yhat)
            local.R2[i] <- (TSSw - RSSw)/TSSw
        }
        AIC <- diags[1]
        AICc <- diags[2]
        edf <- diags[3]
        enp <- diags[4]
        gw.R2 <- diags[6]
        gwR2.adj <- diags[7]
        GW.diagnostic <- list(RSS.gw = RSS.gw, AIC = AIC, AICc = AICc, 
            enp = enp, edf = edf, gw.R2 = gw.R2, gwR2.adj = gwR2.adj)
        Ftests <- list()
        if (F123.test) {
            F.test.parameters <- list(dp.n = dp.n, var.n = var.n, 
                dMat = dMat, x = x, bw = bw, adaptive = adaptive, 
                kernel = kernel, betas = betas, RSS.lm = sum(lms$residuals^2), 
                DF.lm = lms$df.residual, RSS.gw = RSS.gw, tr.S = tr.S, 
                tr.StS = tr.StS, Q = Q)
            Ftests <- F1234.test(F.test.parameters)
        }
    }
    GW.arguments <- list(formula = formula, rp.given = rp.given, 
        hatmatrix = hatmatrix, bw = bw, kernel = kernel, adaptive = adaptive, 
        p = p, theta = theta, longlat = longlat, DM.given = DM1.given, 
        F123.test = F123.test)
    if (hatmatrix) {
        if (is.null(W.vect)) {
            gwres.df <- data.frame(betas, y, yhat, residual, 
                CV, Stud_residual, betas.SE, betas.TV, local.R2)
            colnames(gwres.df) <- c(c(c(colnames(betas), c("y", 
                "yhat", "residual", "CV_Score", "Stud_residual")), 
                paste(colnames(betas), "SE", sep = "_")), paste(colnames(betas), 
                "TV", sep = "_"), "Local_R2")
        }
        else {
            gwres.df <- data.frame(betas, y, yhat, residual, 
                CV, Stud_residual, betas.SE, betas.TV, W.vect, 
                local.R2)
            colnames(gwres.df) <- c(c(c(colnames(betas), c("y", 
                "yhat", "residual", "CV_Score", "Stud_residual")), 
                paste(colnames(betas), "SE", sep = "_")), paste(colnames(betas), 
                "TV", sep = "_"), "E_weigts", "Local_R2")
        }
    }
    else {
        if (is.null(W.vect)) 
            gwres.df <- data.frame(betas)
        else {
            gwres.df <- data.frame(betas, W.vect)
            colnames(gwres.df) <- c(colnames(betas), "E_weigts")
        }
    }
    rownames(rp.locat) <- rownames(gwres.df)
    if (is(regression.points, "SpatialPolygonsDataFrame")) {
        polygons <- polygons(regression.points)
        rownames(gwres.df) <- sapply(slot(polygons, "polygons"), 
            function(i) slot(i, "ID"))
        SDF <- SpatialPolygonsDataFrame(Sr = polygons, data = gwres.df, 
            match.ID = F)
    }
    else {
        SDF <- SpatialPointsDataFrame(coords = rp.locat, data = gwres.df, 
            proj4string = CRS(p4s), match.ID = F)
        if (griddedObj) 
            gridded(SDF) <- T
    }
    timings[["stop"]] <- Sys.time()
    res <- list(GW.arguments = GW.arguments, GW.diagnostic = GW.diagnostic, 
        lm = lms, SDF = SDF, timings = timings, this.call = this.call, 
        Ftests = Ftests, yhat = lass.yhat)
    class(res) <- "gwrm"
    invisible(res)
}
 