# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# Helper Functions --------------------------------------------------------

.weight.calc = function(X,tau,init=FALSE){
  # helper function to either calculate weights based on Reconstruction OR init randomly (based on dim of the data matrix X)

  # get dimensions of argument
  n = nrow(X)
  p = ncol(X)

  if(init) { # if used for random init of weights
    weights = matrix(runif(n=n*p),nrow=n,ncol=p) # output random weigths at desired levels
  } else {
    weights = ifelse(X>=0,tau,1-tau) # calculate weights of argument
  }
  return(weights) # returns weigths as matrix of same size
}

.check.loss = function(X, tau) { # calculates scaled quantile check loss (used in cross validation)
  return(sum(.weight.calc(X,tau)*abs(X))/length(X))
}


# Internal minimization function -----------------------------------------

.pqc.optim = function(X,Z,basis,tau, lambda = 1, const=FALSE) { # This is the function containing the Quadratic Programming
  not.na = !is.na(X) # X is a vector, so check how many entries are not empty
  X = X[not.na] # remove NA entries from vector
  Z = Z[not.na] # same for the weight vector


  n = length(X) # remaining number of observations in data column/row

  k = ncol(basis) # target rank, taken from the basis

  if (!is.null(k)) { # excpetion for k=1 (then basis is not a matrix)
    basis = matrix(basis[not.na,],ncol=k) # make sure that basis is a matrix
  } else { # if not then assign k=1 and convert basis to vector
    k = 1
    basis = matrix(basis[not.na],ncol=1)
  }

  s.mat = matrix(0, nrow=n,ncol=2*n) # this forms the last columns of the constraint matrix (called B in document)

  for (i in 1:n) { # adds indicators for every positive negative residual (could have been done by joining indicators)
    s.mat[i,i] = 1
    s.mat[i,i+n] = -1
  }
  A = cbind(basis/sqrt(matrix(Z, nrow=n,ncol=k,byrow=FALSE)), s.mat) # first columns of constraint matrix, basis divided by weights

  u.mat = diag(x = lambda,nrow = k,ncol=k) # regularization factors corresponding to the v or u to be estimated

  H = matrix(0, nrow = k+2*n,ncol=k+2*n) # matrix in the objective function
  H[1:k,1:k] = u.mat # regulization factors
  H[(k+1):(n+k),(k+1):(n+k)] = tau^2*diag(ncol=n,nrow=n)*2 # factors for positive residuals
  H[(n+k+1):(2*n+k),(n+k+1):(2*n+k)] = (1-tau)^2*diag(ncol=n,nrow=n)*2 # factors for negative residuals
  b = as.matrix(X/sqrt(Z))

  d = rep(0,k+2*n) # make sure no order 1 part in objective function

  alpha = quadprog::solve.QP(Dmat = H,dvec = d,Amat = t(A),bvec = b,meq=n) # call quadprog, meq makes sure everything is handled as equality
  out = alpha$solution[1:k] # we only need the first k results
  return(out)
}


#' Wrapper for the optimization function
#' @importFrom foreach foreach %dopar%
.pqc.next = function (X, Z, basis, tau, lambda, doPar, muEst) {
  n = ncol(X) # get data matrix dimensions
  p = nrow(X)
  k = ncol(basis) # get desired rank

  if (muEst) { # if mu is to be estimated (only for calc of U)
    X = X - basis[,1] # substract mu
    k = k-1 # reduce the rank in the function
    basis = basis[,-1] # remove mu from basis
  }

  out = matrix(NA, nrow = p, ncol = k)# initialize output matrix

  if (doPar) { # if parallel is activated us this
    out = foreach(i = 1:p, .combine = 'rbind') %dopar% {   # this will calculate every row of U/V at once
      .pqc.optim(X[i,, drop=FALSE], Z[i,,drop=FALSE], basis, tau, lambda) # call to quadratic programming wrapper (above)
    }
  } else { # if not parallel
    for(i in 1:p) { # loop for every row in U/V
      out[i,] = .pqc.optim(X[i,, drop=FALSE], Z[i,,drop=FALSE], basis, tau, lambda) # call to quadratic programming wrapper (above)
    }
  }

  if (muEst) { # now adjust U in case mu was estimated (add all-one column)
    out = cbind(rep(1,p), out)
  }
  return(out)
}

# Main PQC algorithm -----------------------------------------------------

#' Main quantile components algorithm
#'
#' @param data data matrix with rows as data entries and columns as variables
#' @param projDim no. of desired principal components
#' @param tau regulization parameter
#' @param lambda regulization parameter
#' @param muEst calculate a constant mu
#' @param epsilon approxiamtion parameter
#' @param iterTol no. of max iterataions
#' @param convTol set algorithm to stop if weigths did not change for convTol no. of consecutive iterations (deactivated for tau=0.5)
#' @param preOut continues based on previous output (if doSeq is TRUE, the components of a previos PCA method are sufficient)
#' @param progBar optional progressbar
#' @param doPar use parallel backend (*nix systems recommended),
#' @param doSeq run via sequential optimization
#' @param randomSeed select a random seed as in set.seed for fixed initialization
#'
#' @return list object containing components and loadings
#'
#' @export
pqcomp = function(data, # data matrix with rows as data entries and columns as variables
                  projDim, # no. of desired principal components
                  tau, # asymmetry parameter (b'w 0 and 1)
                  lambda, # regulization parameter
                  muEst = TRUE,# calculate a constant
                  epsilon = 1/log10(nrow(data)), # approxiamtion parameter
                  iterTol = 200, # no. of max iterataions
                  # Optionals
                  convTol = NA, # set algorithm to stop if weigths did not change for convTol no. of consecutive iterations (deactivated for tau=0.5)
                  preOut = NA, # continues based on previous output (if doSeq is TRUE, the components of a previos PCA method are sufficient)
                  progBar = TRUE, # optional progressbar
                  doPar =TRUE, # use parallel backend (*nix systems recommended),
                  doSeq = FALSE, # run via sequential optimization
                  randomSeed = NA # select a random seed as in set.seed for fixed initialization
)
{
  ## this is entirely optional but will show a progress bar, loss value, and weight convergence indicator while the algorithm runs
  if (progBar) {
    pb = progress::progress_bar$new( # creates template for progress bar
      format = "  estimating tau: :tau [:bar] :percent loss: :approx :wconv",
      total = iterTol, clear = FALSE, width= 70)
    pb$tick(0)
  }
  if (doPar) { # installs required parallel packages
    cores = parallel::detectCores() - 1 # sets no. of cores to number of available cores (minus 1)
    cl = parallel::makeForkCluster(cores)
    doParallel::registerDoParallel(cl)
  }

  if(!is.na(convTol) && tau==0.5) { # exception not to allow early stopping if tau=0.5
    warning("Cannot use timeout if tau=0.5, setting convTol=NA")
    convTol=NA
  }

  # get data dimensions and rename rank
  n = nrow(data)
  p = ncol(data)
  k = projDim

  loss = c() # init container for loss values

  conv = 0 # convergence break counter if convTol!=NA

  X = t(data) # transposes data matrix
  na = is.na(X) # save indeces of NA values

  # Init U and V with random values
  if (!is.na(randomSeed)) { # exception if seed was given
    set.seed(randomSeed)
  }
  if (!is.na(preOut)) { # check for previous output
    V = preOut$components # take components from previous output
    Z = matrix(1,nrow=p,ncol=n) # initialize weigths at 1 (no weight)
    if (!doSeq) { # sequential algo only requires a base
      U = preOut$loadings # take loadings from previous output
      Z = preOut$weights # take weights from previous output
    }
  } else { # regular case: random initialization
    Z = matrix(1,nrow=p,ncol=n) # init weights
    Z[na] = NA # also set NA values were X are NA (avoid division by NA later)
    U = matrix(rnorm(n*k, mean=0, sd=1), nrow=n, ncol=k) # init U randomly
    V = matrix(rnorm(p*k, mean=0, sd=1), nrow=p,ncol=k) # init V randomly

    if (muEst) { # if mu is estimated, do an initial calculation of mu
      mu = .pqc.next(X, Z, basis = matrix(1, ncol=1, nrow=n), tau=tau, lambda = lambda, doPar=doPar, muEst=FALSE) # this will calculate baseline mu

      V = cbind(mu, V) # adjoin V with mu
      U = cbind(rep(1,n),U) # adjoin U with all-one column
    }
  }

  R.old = matrix(0, nrow=p, ncol=n) # init Residiual matrix for previous iteration

  # Set-up loop requirements
  continue=TRUE
  iter = 0

  while (continue) { # loop until continue=F
    # calculate loadings
    U.next = .pqc.next(X=t(X), Z=t(Z), basis = V, tau = tau, lambda = lambda, doPar=doPar, muEst = muEst) # here muEst tells the minimization wrapper to demean X and add the all-on column back later
    if (doSeq) {
      U = U.next # give loadings directly to the components calculation
    }
    V.next = .pqc.next(X=X, Z=Z, basis = U, tau = tau, lambda = lambda, doPar=doPar, muEst = FALSE) # calculate the components matrix (here muEst says not to demean X)

    # Calculate residuals
    R = X - V.next%*%t(U.next)

    # Calculate the quantile check loss
    diff = sum(abs(.weight.calc(R[!na],tau)*R[!na]))

    # Check wether the weights from the previous iteration are the same
    weight.conv = all(.weight.calc(R[!na],tau)==(.weight.calc(R.old[!na],tau)))

    # Save check loss for later information
    loss = append(loss,diff)

    # Calculate weigths for the next iteration
    Z.next = sqrt((.weight.calc(R,tau)^2)*(R^2)+epsilon^2)

    # Prepare for next iteration, rename all matrices
    R.old = R
    V = V.next
    U = U.next
    Z = Z.next

    iter = iter+1 # increase iteration number

    if (iter>iterTol) { # break if max iterations are reached
      continue=FALSE
    } else { # if not
      if (progBar) { # show progress if progBar is activated
        pb$tick(tokens = list(tau = tau, approx = floor(diff), wconv = ifelse(weight.conv==TRUE, "[c]", "[.]")))
      } else{ # print simple indicator that iteration has ended
        cat(".",file=stderr())
      }
    }
    if (weight.conv==TRUE) { # additional exception if algo is set to timeout when weights converged
      if (!is.na(convTol)) { # increase a counter everytime the weights are the same
        conv = conv + 1
        continue = !(conv>convTol) # set algo to break if weights stayed the same for as long as specified
      }
    } else {
      conv = 0 # reset counter everytime weights change
    }
  }
  print("Finished!")
  if (doPar) { # parallel needs to stop the processes
    parallel::stopCluster(cl)
  }

  # Finally name the components matrix accordingly
  if (muEst) {
    colnames = c("constant")
  } else {
    colnames = c()
  }
  for(i in 1:k) {colnames = append(colnames,paste("PC",i,sep=""))}

  colnames(V) = colnames
  rownames(V) = colnames(data)

  # Output as a list containing some additional information
  return(list(components=V,loadings=U, loss = loss, converged = weight.conv))
}


# Cross Validation --------------------------------------------------------

.createNA <- function (x, shareNA) { # routine to create random NA values in a dataset according to % size
  n = nrow(x) # dimensions
  p = ncol(x)
  NAloc = rep(FALSE, n * p) # init indeces as all not NA
  NAloc[sample.int(n * p, floor(n * p * shareNA))] = TRUE # randomly add TRUE for holdout
  x[matrix(NAloc, nrow = n, ncol = p)] = NA # set data to NA at the TRUE locations
  return(x)
}

#' Cross validation wrapper for pqcomp(): Select ranges for rank and lambda.
#' Results are calculated on a random hold out on the supplied dataset.
#' See Madeleine Udell, Corinne Horn, Reza Zadeh and Stephen Boyd (2016), "Generalized Low Rank Models", Foundations and Trends® in Machine Learning: Vol. 9: No. 1, pp 1-118.
#'
#' @param data data matrix with rows as data entries and columns as variables
#' @param projDimRange # range of desired principal components
#' @param tau asymmetry parameter (b'w 0 and 1)
#' @param lambdaRange regulization parameter
#' @param muEst calculate a constant mu
#' @param epsilon approxiamtion parameter
#' @param shareNA size of holdout (0<shareNA<1)
#' @param iterTol no. of max iterataions
#' @param convTol set algorithm to stop if weigths did not change for convTol no. of consecutive iterations (deactivated for tau=0.5)
#' @param preOut continues based on previous output (if doSeq is TRUE, the components of a previos PCA method are sufficient)
#' @param progBar optional progressbar
#' @param doPar use parallel backend (*nix systems recommended),
#' @param doSeq run via sequential optimization
#'
#' @return list containing pqcomp() results
#'
#' @export
cv.pqc = function(data,
                  projDimRange,
                  tau,
                  lambdaRange,
                  muEst = TRUE,
                  epsilon = 1/log10(nrow(data)),
                  shareNA = 0.1,
                  iterTol = 200,
                  convTol = NA,
                  progBar = TRUE,
                  doPar = TRUE,
                  doSeq = FALSE) {
  num = length(lambdaRange)*length(projDimRange) # how many parameter combinations
  Y.na = .createNA(data, shareNA) # create NA values in data
  na.loc = is.na(Y.na) # save location of true parameters
  Y.true = data[na.loc] # save true values

  results = list() # init result list

  i = 1 # init counter
  for (lambda in lambdaRange) { # loop for every lambda
    for (k in projDimRange) { # loop for every projDim
      print(paste(i,num,sep="/")) # just print at which iteration we are

      # regular PQC call with parameter combination
      pqc = pqcomp(Y.na, projDim = k,tau = tau,lambda = lambda, muEst=muEst, epsilon = epsilon, iterTol = iterTol, convTol = convTol, progBar = progBar, doPar = doPar,  doSeq=doSeq)

      na.pred = (pqc$loadings%*%t(pqc$components))[na.loc] # save predictions

      result = list(lambda = lambda, k=k, na.pred = na.pred, na.true = Y.true,loss=pqc$loss[length(pqc$loss)], converged = pqc$converged) # create list of result infos
      results = append(results, list(result)) # adjoin to results
      i = i+1 # increment counter
    }
  }
  return(results)
}

#' Calculate cross validation scores via quantile check loss on the result of cv.pqc().
#'
#' @param cv.result list output of cv.pqc()
#' @param tau asymmetry parameter (0<tau<1)
#'
#' @return list containing pqcomp() results
#'
#' @export
cv.scores = function(cv.result, tau) { # This will calculate the CV score using the object created by cv.pqc
  scores = lapply(cv.result,function(x) data.frame(k=x$k, lambda = x$lambda, loss = x$loss, cv.score = .check.loss(x$na.true-x$na.pred,tau))) # calculate the CV score based on true and predicted values for every result
  score.df = do.call("rbind", scores) # put scores into dataframe
  score.df["lambda"] = factor(score.df$lambda) # put lambda in the column into factors
  score.df["cv.score.resc"] = rescale(score.df$cv.score, to = c(0,1)) # add a rescaled score for comparison
  return(score.df)
}
