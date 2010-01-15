joint.glm <- function(form.mean, form.disp, data,
                      family.mean = gaussian, family.disp = Gamma(link='log'),
                      eps = 10^(-6), iter.max = 100, maxit = 50, reml = TRUE){

  ##This function will fit the Joint Model through the framework describe
  ##in [McCullagh1989]
  ##By default, it fit a gaussian distribution for the mean component
  ##and a Gamma distribution with log link function for the dispersion
  ##component

  ##form.mean     : a R ``formula'' which specify the linear model for the
  ##                mean
  ##form.disp     : a R ``formula'' which specify the linear model for the
  ##                dispersion
  ##data          : a ``data.frame'' corresponding to the data.
  ##                The first column correspond to the response
  ##                observations, while others to the explicative
  ##                variables.
  ##family.mean   : The ``family'' for the mean component.
  ##family.disp   : The ``family'' for the dispersion component.
  ##eps           : Optional Numeric. It specifies the precision
  ##                for the convergence test.
  ##iter.max      : The maximum number of iteration to fit the Joint
  ##                Model.
  ##maxit         : The maximum number of iteration to fit a generic
  ##                GLM - see function ``glm.fit''
  ##reml          : Logical. If TRUE, the Restricted Maximum Likelihood
  ##                Estimation is used.

  ##Convert formulas as real R formulas
  form.mean <- as.formula(form.mean)
  form.disp <- as.formula(form.disp)
  
  ##Set the initial weigth vector
  weights <- rep(1, length(data[,1]))
  ##Set an arbitrary value for the EQL beacause
  ##difference between two iterations will be a criteria
  ##for our test of convergence.
  dev.new <- Inf
  
  ##Introduce a logical object ``flag'' which specifie
  ##if convergence is reached or not.
  ##   TRUE    ->   convergence test OK
  ##   FALSE   ->   convergence test non OK
  ##Initialize it to FALSE to enter in the ``while'' loop
  flag <- FALSE

  ##Introduce also a counter for the number of
  ##iterations
  n.iter <- 0

  ##First, fit the GLM associated to the mean
  mod.mean <- glm(form.mean, family = family.mean, data = data,
                  maxit = maxit, weights = weights)
  
  while ( !flag ){

    ##We can then compute the dispersion response
    ##i.e. the deviance contribution
    mod.mean.alias <- mod.mean          #craation of an alias of
                                        #``mod.mean'' object
    mod.mean.alias$prior.weights <- rep(1, length(data[,1]))
                                        #Set the prior weights
                                        #to  1 to evaluate the deviance
                                        #contribution
    d <- residuals(mod.mean.alias, 'deviance')^2
    
    if (reml)
      d <- d / (1 - hatvalues(mod.mean))
    
    data.disp <- cbind(d = d, data[,-1])
    data.disp <- as.data.frame(data.disp)
    names(data.disp) <- c('d', names(data)[-1])

    weights.disp <- rep(1, length(d))
    if (reml)
      weights.disp <- 1 - hatvalues(mod.mean)
    
    mod.disp <- glm(form.disp, family = family.disp,
                    data = data.disp, weights = weights.disp,
                    maxit = maxit)
    
    gam <- mod.disp$coefficients
    phi <- mod.disp$fitted

    ##Update for the mean weights
    weights <- 1 / phi
    mod.mean <- glm(form.mean, family = family.mean, data = data,
                    maxit = maxit, weights = weights)

    dev.old <- dev.new
            
    dev.new <- eql(mod.mean,mod.disp)["eql"]
    cat("Percentage of EQL variation is : ",
        abs( (dev.old - dev.new) / dev.old * 100),"\n")
    flag <- ( iter.max < n.iter) || abs( (dev.old - dev.new) / dev.new * 100) < eps
   
    n.iter <- n.iter + 1
    
  }

  return(list(mod.mean = mod.mean, mod.disp = mod.disp,
              iterations = n.iter, eql = dev.new))
}

eql <- function(mod.mean, mod.disp, df.adj = TRUE){
  ##A function which compute the Extended Quasi-Deviance
  ##(EQD) for the Joint Model

  ##mod.mean : Object of class ``glm'' corresponding to
  ##           the mean component
  ##mod.disp : Object of class ``glm'' corresponding to
  ##           the dispersion component
  ##df.adj   : Logique. If ``TRUE'' adjustment on degree
  ##           of freedom is applied

  dispersion <- mod.disp$fitted              #dispersion parameter
                                             #estimate
  dev.res <- mod.disp$y
  
  fct.var <- mod.mean$family$variance        #Variance function

  if(df.adj)
    eqllik <- - dev.res / (2 * dispersion)  -
      1/2 * mod.mean$df.residual / length(dispersion) *
        log(2*pi*dispersion*fct.var(mod.mean$y))
  else
    eqllik <- -dev.res / (2 * dispersion) -
      1/2 * log(2*pi*dispersion*fct.var(mod.mean$y))
  
  eql <- sum(eqllik)
  eqd <- -2 * sum(dispersion * eqllik)

  return(c(eql = eql, eqd = eqd))
}
