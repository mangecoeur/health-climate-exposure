import numpy as np
from numpy import sum, mean, nansum, floor, round, log, isnan, median, var, Inf


def acs(execesses):
    """Compute average cluster size for a given three-dimensional array of excesses

    Description:

          Returns average cluster size

    Usage: acs(execesses)

    Args:

          excesses: p (space) x q (space) x n (time) array of excesses

    Returns:

          A p x q matrix containing the average cluster size

    Authors:

         Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 28 Feb 2006
         Chris Ferro <c.a.t.ferro@reading.ac.uk>
    """




    def avenumbexceed(y):
        x = matrix(y, ncol=3, nrow=len(y) / 3)
        nansum(x > 0) / nansum(apply(x > 0, 1, any, na_rm = T))


    apply(execesses, c(1, 2), avenumbexceed)


def anomalyfield(data, firstmonth, firstyear, calc="mac"):
    """
    Calculate anomalies by subtracting out mean annual cycle or long-term mean of
    a three-dimensional p x q x n array of monthly mean values. The structure of this
    array must be as follows: - first dimension p is longitude points
                              - second dimension q is latitude points
                              - third dimension n is time (months of the year in
                                ascending order, e.g. Jan 1948, Feb 1948, **kwargs, Sep 2005)

    Description:

         Returns a p x q x n array containing anomalies

    Usage:

         anomalyfield(data,firstmonth,firstyear,calc="mac")

    Args:

         data:       p x q x n array (see data structure description above)
         firstmonth: number of the first month containing data in the
                     data array (e.g. 1 (for January))
         firstyear:  first year of data in the data array (e.g. 1948)
         calc:       string of characters that can be "mac" (Default) to  calculate
                     anomalies by subtracting the mean annual cycle, or "ltm" to
                     calculate anomalies by subtracting the long-term mean.

    Returns:

         .anomaly:   p x q x n array containing anomalies

    Authors:

         Caio Coelho <c.a.d.s.coelho@reading.ac.uk>  28 October 2005
         Chris Ferro <c.a.t.ferro@reading.ac.uk>

    Examples:

         x <- seq(-20, 20, 5)
         y <- seq(30, 60, 5)
         dim <- c(len(x), len(y), 100)
         data <- array(rnorm(prod(dim)), dim)
         out<-anomalyfield(data,1,1948)
         out<-anomalyfield(data,5,1970)
         out<-anomalyfield(data,5,1970,calc="ltm")
    Returns
    """

    if data.shape[2] < 12:
        raise RuntimeError("Not possible to  compute either mean annual cycle or long-term mean. Less than a year of data available")

    if (calc == "mac"):
        def mean_month (data):
            apply(matrix(data, nrow=12), 1, mean)

    anomaly = data

    if (firstmonth == 1):
        initialmonth = firstmonth
        finalmonth = floor(dim(data)[3] / 12) * 12
    elif (firstmonth != 1):
        initialmonth = 12-firstmonth+2
        finalmonth = (floor((dim(data)[3]-(12-firstmonth+1)) / 12) * 12)+initialmonth-1


    m = aperm(apply(data[, , initialmonth:finalmonth], c(1, 2), mean.month), c(2, 3, 1))
    anom = data[, , initialmonth:finalmonth]
    for i in range(1, ((finalmonth - initialmonth + 1) / 12)):
        if (firstmonth == 1):
            anom[, , (((i - 1) * 12) + 1):(12 * i)] = data[, , (((i - 1) * 12) + 1):(12 * i)] - m
        else:
            anom[, , (((i - 1) * 12) + 1):(12 * i)] = data[, , (((i - 1) * 12) + initialmonth):((12 * i) + initialmonth - 1)] - m


    if firstmonth == 1 & finalmonth == dim(data)[3]:
        anomaly = anom


    if firstmonth == 1 & finalmonth < dim(data)[3]:
        anomaly[:,: , firstmonth:finalmonth] = anom
        j = 1
        for i in range((finalmonth + 1), dim(data)[3]):
            anomaly[:,:, i] = data[:,:, i] - m[:,:, j]
            j = j + 1

    if firstmonth != 1 & finalmonth == dim(data)[3]:
        j = 1
        for i in range(firstmonth,12):
             anomaly[:,:, j] = data[:,:, j] - m[:,:, i]
            j = j + 1

    anomaly[, , initialmonth:finalmonth] = anom


    if firstmonth != 1 & finalmonth < dim(data)[3]:
        j = 1
        for i in range(firstmonth,12):
            anomaly[:,:, j] = data[:,:, j] - m[:,:, i]
            j = j + 1

    anomaly[:,: , initialmonth:finalmonth] = anom
    j = 1
    for i in range(finalmonth + 1):
        dim(data)[3]
        anomaly[:,:, i] = data[:,:, i] - m[:,:, j]
        j = j + 1




    if (firstmonth == 1):
        y1 = firstyear
    if (firstmonth != 1):
        y1 = firstyear+1

    y2 = y1 +((dim(anom)[3]) / 12)-1

    print("Mean annual cycle was computed with", dim(anom)[3], "months of data")
    print("from Jan", y1, "to Dec", y2, "(i.e.", ((dim(anom)[3]) / 12), "years)")



    if (calc == "ltm"):
        anomaly = data

    m =apply(data, c(1, 2), mean)
    for i in range(1, dim(data)[3])):
        anomaly[:,:, i] = data[:,:, i] - m

    print("Long-term mean was computed with", dim(data)[3], "months of data",)


    dict(anomaly=anomaly)

def momentskew(z):
    mean(((z - mean(z, na_rm = T)) / sd(z, na_rm=True)) ^ 3, na_rm = T)
def ykskew (z):
    q =quantile(z, na_rm=True)
    (q[2] - 2 * q[3] + q[4]) / (q[4] - q[2])
def quantil (w, p=0.5):
    quantile(w, p, na_rm=True)


def boundexcesses(paretoparam):
    # Compute upper bound of excesses for a given array of Generalized Pareto
    # distribution parameters
    #
    #  Description:
    #
    #       Returns upper bound of excesses
    #
    # Usage: boundexcesses(paretoparam)
    #
    # Input:
    #
    #       paretoparam: 2 x p x q array with scale (first level of the array) and
    #              shape (second level of the array) parameters
    #
    # Output:
    #
    #       A matrix containing the upper bound of excesses
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 28 Feb 2006


    paretoparam[2,,][paretoparam[2,, ] >= -0.008] =0
    b =-(paretoparam[1,,] / paretoparam[2,,])
    b[b == -Inf] =10000
    return b


def corfield (array3d, ts, **kwargs):
    # Compute correlation between each point of a three-dimensional array, with
    # first two space dimensions and third time dimension, and a given vector of
    # the same len as the time dimension of the three-dimentional array.
    #
    # Description:
    #
    #      Returns a space matrix with the correlations
    #
    # Usage:
    #
    #      corfield(array3d,ts, **kwargs)
    #
    # Input:
    #
    #      array3d: three-dimensional array with first two space dimensions
    #               (e.g. longitude and latitude) and third time dimension
    #
    #      ts: vector containing the time series to be correlated with the
    #          three-dimensional array
    #
    #     **kwargs: Additional argument passed to the correlation function
    #          (e.g. method = "pearson", method = "kendall", method = "spearman",
    #           use = "all.obs", use = "complete.obs", use = "pairwise.complete.obs")
    #          Default computes pearson correlation coefficient (method = "pearson")
    #          using all data available (use = "all.obs")
    #
    # Output:
    #
    #          Matrix of correlations
    #
    # Authors:
    #
    #      Dag Johan Steinskog <dag.johan.steinskog@nersc.no> 9 June 2005
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 100)
    #      data <- array(rnorm(prod(dim)), dim)
    #      timeseries<-rnorm(100)
    #      corfield(data,timeseries)
    #      corfield(data,timeseries,method = "pearson")
    #      corfield(data,timeseries,method = "kendall")
    #      corfield(data,timeseries,method = "spearman")


    def correl(z):
        cor(z, ts, use="complete.obs", **kwargs)
    apply(array3d, c(1, 2), correl)


def covfield(array3d, ts, **kwargs):
    # Compute covariance between each point of a three-dimensional array, with
    # first two space dimensions and third time dimension, and a given vector of
    # the same len as the time dimension of the three-dimentional array.
    #
    # Description:
    #
    #      Returns a space matrix with the covariances
    #
    # Usage:
    #
    #      covfield(array3d,ts, **kwargs)
    #
    # Input:
    #
    #      array3d: three-dimensional array with first two space dimensions
    #               (e.g. longitude and latitude) and third time dimension
    #
    #      ts: a vector containing the time series to be covariate with the
    #          three-dimensional array
    #
    #     **kwargs: Additional argument passed to the covariance function
    #          (e.g. method = "pearson", method = "kendall", method = "spearman",
    #           use = "all.obs", use = "complete.obs", use = "pairwise.complete.obs")
    #          Default computes pearson covariance (method = "pearson")
    #          using all data available (use = "all.obs")
    #
    # Output:
    #
    #      Matrix of covariances
    #
    # Authors:
    #
    #      Dag Johan Steinskog <dag.johan.steinskog@nersc.no> 25 October 2005
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 100)
    #      data <- array(rnorm(prod(dim)), dim)
    #      timeseries<-rnorm(100)
    #      covfield(data,timeseries)
    #      covfield(data,timeseries,method = "pearson")
    #      covfield(data,timeseries,method = "kendall")
    #      covfield(data,timeseries,method = "spearman")

    def  covar (z):
        cov(z, ts, use="complete.obs", **kwargs)
    apply(array3d, c(1, 2), covar)



def detrend(y, smooth=False, **kwargs):
    # Description:
    #
    #   Removes time trend from a time series
    #
    # Usage:
    #
    #   detrend(y,smooth=False,**kwargs)
    #
    # Arguments:
    #
    #   y:      a vector containing a time series of data
    #
    #   smooth: Logical. If False (the default) removes linear time trend by
    #           fitting ordinary least squares regression to the time series y.
    #           If True, then time trend is removed by fitting a local weighted
    #           regression (lowess) to the time series y
    #
    #   **kwargs:    Additional arguments passed to lowess
    #
    # Value:
    #
    #   Vector with the detrended data.
    #
    # Author:
    #
    #   Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 9 Jan 2006
    #
    # Examples:
    #
    #      ts <- rnorm(100,25,2)
    #      detrend(ts)
    #      detrend(ts,smooth=True)
    #      detrend(ts,smooth=True,f=1/3)

    if (not smooth):
        x = seq(np.arange(1, len(y))
        lm(y ~ x).res

    else:
        y - lowess(y, **kwargs).y




def detrendfield(array3d, smooth=F, **kwargs):


    # Description:
    #
    #   Removes time trend at each grid point of a three-dimensional array of data
    #   with first two space dimensions (e.g. longitude and latitude) and
    #   third time dimension
    #
    # Usage:
    #
    #   detrendfield(array3d,smooth=False,**kwargs)
    #
    # Arguments:
    #
    #   array3d: a three-dimensional array, with first two space dimensions (e.g.
    #            longitude and latitude) and third time dimension
    #
    #   smooth: Logical. If False (the default) removes linear time trend by
    #           fitting ordinary least squares regression at each grid. If True,
    #           then time trend is removed by fitting a local weighted regression
    #           (lowess) at each grid point
    #
    #   **kwargs:    Additional arguments passed to lowess
    #
    # Value:
    #
    #   Three-dimensional array with the detrended data.
    #
    # Author:
    #
    #   Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 9 Jan 2006
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 100)
    #      data <- array(rnorm(prod(dim),25,2), dim)
    #      detrendfield(data)
    #      detrendfield(data,smooth=True)
    #      detrendfield(data,smooth=True,f=1/3)


    aperm(apply(array3d, c(1, 2), detrend, smooth, **kwargs), c(2, 3, 1))


def eof(lon, lat, array3d, **kwargs):
    # Compute EOFs and PCs of a three-dimensional p x q x n array, with first
    # two space dimensions p (longitude) and q (latitude), and third time
    # dimension n.
    #
    # Description:
    #
    #      Returns principal component time series, spatial patterns (loadings),
    #      standard deviation of each principal component, percentage of total
    #      variance accounted by each principal component, a vector of longitudes
    #      and a vector of latitudes
    #
    # Usage:
    #
    #      eofs(lon,lat,array3d)
    #
    # Input:
    #
    #      array3d: a three-demensional array with p longitude points and q latitude
    #               points as the first two dimensions and n as the third time
    #               dimension
    #
    #      lon: vector of longitudes
    #
    #      lat: vector of latitudes
    #
    #      **kwargs: Additional arguments to be passed for prcomp function such as,
    #           retx:   logical value indicating whether the rotated variables
    #                   should be returned.
    #           center: logical value indicating whether the variables should
    #                   be zero centered. Alternately, a vector of
    #                   len equal the number of columns of x can be supplied.
    #                   The value is passed to scale.
    #           scale.: logical value indicating whether the variables should
    #                   be scaled to have unit variance before the analysis takes
    #                   place. The default is False for consistency with S, but
    #                   in general scaling is advisable. Alternatively, a vector
    #                   of len equal the number of columns of x can be supplied.
    #                   The value is passed to scale.
    #              tol: value indicating the magnitude below which components
    #                   should be omitted. (Components are omitted if their
    #                   standard deviations are less than or equal to tol times the
    #                   standard deviation of the first component). With the default
    #                   null setting, no components are omitted. Other settings for
    #                   tol could be tol = 0 or tol = sqrt(.Machine.double.eps),
    #                   which would omit essentially constant components.
    # Output:
    #
    #               .pcs: matrix of principal components (first column contains
    #                     first principal component and so on)
    #             .stdev: vector containing the standard deviation of each principal
    #                     component, starting from the first.
    #               .eof: three-dimensional array containing the loadings.
    #      .accountedvar: vector containing the amount of variance accounted for
    #                     each principal component
    #               .lon: vector of longitudes
    #               .lat: vector of latitudes
    #
    #
    # Authors:
    #
    #      Dag Johan Steinskog <dag.johan.steinskog@nersc.no> 10 June 2005
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #
    # Examples:
    #
    #      lon <- seq(50,70,20)
    #      lat <- seq(-10,10,20)
    #      t <- 20
    #      x <- array(rnorm(len(lon)*len(lat)*t,sd=10),dim=c(2,2,20))
    #      eofs(lon,lat,x)

    data = reshapefield(lon, lat, array3d).out

    aux = prcomp(data, **kwargs)

    eofs =reshapefield(lon, lat, t(aux.r)).out

    varsquared = (aux.sdev) ^ 2
    sumvarsquared = sum(varsquared)

    accountedvar = round((varsquared / sumvarsquared) * 100, 2)

    return dict(pcs=aux.x, stdev = aux.sdev, eofs = eofs, accountedvar = accountedvar, lon = lon, lat = lat)



def fa(x, y, reduction=F, nmodes=3):
    # Bayesian forecast assimilation for cross-validated calibration and combination
    # of forecasts
    #
    # Description: This function calls functions faprior, falikelihood,
    # and faprediction and estimates the mean (ya) and the
    # square root of the diagonal elements of the
    # covariance matrix (d) of the posterior distribution
    # y|x ~ N(ya,d)
    #
    # Usage: fa(x,y,reduction,nmodes)
    #
    # Arguments:
    # x:    n x p data matrix containing time series of p forecast variables
    # y:    n x q data matrix containing time series of q observable
    #   variables
    # reduction: Logical. If True performs maximum covariance analysis
    #            on the cross-covariance matrix (y'*x) and uses the
    #            selected number of modes ('nmodes', see below) to model the
    #            likelihood function. Default is False.
    # nmodes: number of modes used in the likelihood modelling. Default
    #         is 3. Argument only used if 'reduction' is True.
    #
    # Outputs:
    # .ya: q x n matrix containing q predicted observable state time series
    # .d : q x q matrix containing predicted observable error variance
    #      (i.e. the square root of the diagonal elementes of observable
    #       error covariance matrix d.)
    #
    # Authors: Caio Coelho <c.a.d.s.coelho@reading.ac.uk>       13 Feb 2006
    #          David Stephenson <d.b.stephenson@reading.ac.uk>
    #
    # Reference: Stephenson D. B., C.A.S. Coelho, F. J. Doblas-Reyes and
    #            M. Balmaseda, 2005: Forecast Assimilation: A unified
    #            framework for the combination of multi-model weather and
    #            climate predictions.  Tellus A . Vol.  57, 253-264.
    #
    # Example 1: Calibration and combination of forecasts for a single location
    #
    #          xfa<-read.table("xfa.txt") # xfa.txt contains 3 time series of forecasts
    #          yfa<-read.table("yfa.txt") # yfa.txt contains 1 time series of observations
    #          out<-fa(xfa,yfa)
    #          out.ya # combined and calibrated forecast mean
    #          out.d  # combined and calibrated forecast standard deviation
    #
    # Example 2: Calibration and combination of grided coupled model
    #            forecasts
    #
    #          xfa<-read.table("xfasst.txt") # xfa.txt contains 168 time series of forecasts
    #          # Columns 1-56 of xfa.txt are forecasts produced by model A
    #          # Columns 57-112 of xfa.txt are forecasts produced by model B
    #          # Columns 113-168 of xfa.txt are forecasts produced by model C
    #          yfa<-read.table("yfasst.txt") # yfa.txt contains 56 time series of observations
    #          out<-fa(xfa,yfa,reduction=T,3)
    #          out.ya # combined and calibrated forecast mean
    #          out.d  # combined and calibrated forecast standard deviation


    y = as.matrix(y)
    x = as.matrix(x)
    n = dim(x)[1]
    ya =matrix(nrow=n, ncol=ncol(y))
    d =matrix(nrow=n, ncol=ncol(y))
    if (not reduction):
        # Cross-validation loop
        nn = 1
        while (nn < n + 1):
            xcross = as.matrix(x[-nn, ])
            ycross = as.matrix(y[-nn, ])
            # Estimate prior parameters yb (1 x q) and c (q x q)
            prior = faprior(ycross)
            # Estimate likelihood parameters g (p x q), yzerogt (n-1 x p) and
            # s (p x q)
            lik = falikelihood(xcross, ycross)
            # Estimate posterior parameters ya (1 x q) and d (q x q)
            out = faprediction(x[nn, ], ycross, prior.yb, prior.c, lik.g, lik.yzerogt,
            lik.s)
            ya[nn, ] =out.ya
            d[nn, ] =sqrt(diag(out.d))
            nn = nn + 1
          # end of cross-validation loop

    else:
        # Cross-validation loop
        nn = 1
        while (nn < n + 1):
            xcross = as.matrix(x[-nn, ])
            ycross = as.matrix(y[-nn, ])
            mca = svd(t(ycross) % * % xcross)
            xtrunc = xcross % * % mca.v[, np.arange(1,nmodes)]
            ytrunc = ycross % * % mca.u[, np.arange(1,nmodes)]
            # Estimate prior parameters yb (1 x q) and c (q x q)
            prior = faprior(ytrunc)
            # Estimate likelihood parameters g (p x q), yzerogt (n-1 x p) and
            # s (p x q)
            lik = falikelihood(xtrunc, ytrunc)
            # Estimate posterior parameters ya (1 x q) and d (q x q)
            out = faprediction(x[nn,] % * % mca.v[, np.arange(1,nmodes)], ytrunc, prior.yb, prior.c, lik.g, lik.yzerogt, lik.s)
            ya[nn,] = out.ya % * % t(mca.u[, np.arange(1,nmodes)])
            d[nn,] = t(sqrt(diag(mca.u[, np.arange(1,nmodes)] % * % out.d % * % t(mca.u[, np.arange(1,nmodes)]))))
            nn = nn + 1
          # end of cross-validation loop

    return dict(ya=ya, d=d)


def falikelihood(x, y):
    # This function estimates the parameters G, yoG' and S
    # of the likelihood distribution x|y ~ N([y-yo]G',S),
    # where the symbol ' denotes transpose matrix.
    #
    # NOTE: in the code all parameters appear
    # in lower case letters.
    #
    # USAGE: falikelihood(x,y)
    #
    # INPUTS:
    # x  -> n x p data matrix containing time series of p forecast variables
    # y  -> n x q data matrix containing time series of q observable var.
    #
    # OUTPUTS:
    # g  -> p x q matrix containing the calibration operator
    # yzerogt -> 1 x q matrix containing the intercept of the
    #        regression between x and y
    # s  -> p x p forecast error covariance matrix
    #
    # Authors: Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #          David Stephenson <d.b.stephenson@reading.ac.uk>
    #


    #
    # find out p, q and n dimentions
    #
    p = ncol(x)
    q = ncol(y)
    n = nrow(y)
    #
    # Define g, yzerogt and s
    #
    g = matrix(nrow=p, ncol=q)
    yzerogt = matrix(nrow=n, ncol=p)
    s = matrix(nrow=p, ncol=p)
    #
    # Calculate column means of x and y
    #
    # xmean is a 1 x p matrix
    #
    xmean = t(apply(x, 2, mean))
    #
    # ymean is a 1 x q matrix
    #
    ymean = t(apply(y, 2, mean))
    #
    # Estimate covariance matrices
    #
    Sxx = var(x) * ((n - 1) / n)
    Syy = var(y) * ((n - 1) / n)
    Sxy = var(x, y) * ((n - 1) / n)
    Syx = t(Sxy)
    #
    # Inverte Syy
    #
    Syyi = solve(Syy)
    #
    # Estimate forecast error covariance (s)
    #
    s = Sxx - Sxy % * % Syyi % * % Syx
    #
    # Estimate intercept (yzerogt) and slope (g) from
    # the regression between x and y
    #
    # g is a p x q matrix
    #
    g = Sxy % * % Syyi
    #
    # yzerogt is a 1 x q matrix.
    #
    yzerogt =  - (xmean - ymean % * % t(g))
    return dict(g=g, yzerogt=yzerogt, s=s)

def faprediction(x, y, yb, c, g, yzerogt, s):
    # This function estimates the mean ya and the covariance d
    # of the posterior distribution y|x ~ N(ya,d),
    # where ya = yb + [x - yb g' + yo g'] L'
    #       d  = (I - Lg) c
    #       L' = (gcg' +s)^(-1) gc
    #
    # The symbol ' indicates transpose matrix
    #
    # USAGE: faprediction(x, y, yb, c, g, yzerogt, s)
    #
    # INPUTS:
    # x  -> n x p data matrix containing time series of p forecast variables
    # y  -> n x q data matrix containing time series of q observable vars.
    # yb -> 1 x q matrix containing column means of y
    # c  -> q x q observable state covariance matrix
    # g  -> p x q matrix containing the calibration operator
    # yzerogt -> 1 x q matrix containing the intercept of the
    #            regression between x and y
    # s  -> p x q forecast error covariance matrix
    #
    # OUTPUTS:
    #
    # ya -> 1 x q matrix containing q predicted observable state time series
    # d  -> q x q matrix containing predicted observable error covaraiance
    #
    # Authors: Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #          David Stephenson <d.b.stephenson@reading.ac.uk>
    #


    # Find out p, q and n
    #
    p = ncol(x)
    q = ncol(y)
    n = nrow(y)
    #
    # Define matrices ya and d
    #
    ya = matrix(nrow=n, ncol=q)
    d = matrix(nrow=q, ncol=q)
    #
    # Transpose g
    #
    gt = t(g)
    #
    # Calculate gain/weight matrix transpose (Lt)
    #
    Lt = solve(g % * % c % * % gt + s) % * % g % * % c
    #
    # Calculate posterior mean (ya)
    #
    ya = yb + (x - yb % * % gt + yzerogt) % * % Lt
    #
    # Calculate the posterior covariance matrix (d)
    #
    d = (diag(q) - t(Lt) % * % g) % * % c
    #
    list(ya=ya, d=d)

def faprior(y):


    # This function estimates the mean (yb) and the background
    # observable error covariance (c) of the prior distribution
    # of y given by y~N(yb,c)
    #
    # USAGE: faprior(y)
    #
    # INPUT:
    # y  -> n x q data matrix containing q observable time series
    #
    # OUTPUTS:
    # yb -> 1 x q matrix containing the column means of y
    # c  -> q x q observable state covariance matrix
    #
    # Authors: Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #          David Stephenson <d.b.stephenson@reading.ac.uk>
    #



    # Calculate column means of y and store in yb (1 x q)
    yb = t(apply(y, 2, mean))
    #
    #
    # Estimates observable state covariance matrix (c) - biased estimate 1/N
    #
    c = var(y) * (nrow(y) - 1) / nrow(y)
    list(yb=yb, c=c)


def mygpd_fit (xdat, threshold, npy=365, ydat=None, sigl=None,
                      shl=None, siglink=identity, shlink=identity, show=True,
                      method="Nelder-Mead", maxit=10000, **kwargs):


    # Same as gpd.fit function from ismev package but with standard error
    # calculation disactivated

    z = list()
    npsc = len(sigl) + 1
    npsh = len(shl) + 1
    n = len(xdat)
    z.trans = False
    if ( is_function(threshold)):
        raise RuntimeError("`threshold' cannot be a function")
    u = np.repeat(threshold, len.out = n)
    if (len(unique(u)) > 1):
        z.trans = True
    xdatu = xdat[xdat > u]
    xind = np.arange(1,n)[xdat > u]
    u = u[xind]
    in2 = sqrt(6 * var(xdat)) / pi
    in1 = mean(xdat, na_rm=True) - 0.57722 * in2
    if sigl is None:
        sigmat = as.matrix(np.repeat(1, len(xdatu)))
        siginit = in2

    else:
        z.trans = True
        sigmat = cbind(np.repeat(1, len(xdatu)), ydat[xind, sigl])
        siginit = c(in2, np.repeat(0, len(sigl)))

    if shl is None:
        shmat = as.matrix(np.repeat(1, len(xdatu)))
        shinit = 0.1

    else:
        z.trans = True
        shmat = cbind(np.repeat(1, len(xdatu)), ydat[xind, shl])
        shinit = c(0.1, np.repeat(0, len(shl)))

    init = c(siginit, shinit)
    z.model = list(sigl, shl)
    z.link = deparse(substitute(c(siglink, shlink)))
    z.threshold = threshold
    z.nexc = len(xdatu)
    z.data = xdatu
    gpd.lik = function(a)

    sc = siglink(sigmat % * % (a[seq(1, len=npsc)]))
    xi = shlink(shmat % * % (a[seq(npsc + 1, len=npsh)]))
    y = (xdatu - u) / sc
    y = 1 + xi * y
    if (min(sc) <= 0):
        l = 10 ^ 6
    else:
        if (min(y) <= 0):
            l = 10 ^ 6
        else:
            l = sum(log(sc)) + sum(log(y) * (1 / xi + 1))


    l

    x = optim(init, gpd.lik, hessian=True, method=method,
    control = list(maxit=maxit, **kwargs))
    sc = siglink(sigmat % * % (x.par[seq(1, len=npsc)]))
    xi = shlink(shmat % * % (x.par[seq(npsc + 1, len=npsh)]))
    z.conv = x.convergence
    z.nllh = x.value
    z.vals = cbind(sc, xi, u)
    if (z.trans)
    z.data = -log(as.vector((1 + (xi * (xdatu - u)) / sc) ^ (-1 / xi)))

    z.mle = x.par
    z.rate = len(xdatu) / n
    # z.cov <- solve(x.hessian)
    # z.se <- sqrt(diag(z.cov))
    z.n = n
    z.npy = npy
    z.xdata = xdat
    if (show):
        if (z.trans):
            print(z[c(2, 3)])
        if (len(z[[4]]) == 1):
            print(z[4])
        print(z[c(5, 7)])
        if (~z.conv):
            print(z[c(8, 10, 11, 13)])

    return z

def returnperiod(excess, parameter):

    # Compute return period for a given matrix of excesses and a given array of
    # Generalized Pareto distribution parameters
    #
    #  Description:
    #
    #       Returns return period of excesses
    #
    # Usage: returnperiod(excess,parameter)
    #
    # Input:
    #
    #       excess: p x q matrix of excesses
    #
    #       parameter: 2 x p x q array with scale (first level of the array) and
    #                  shape (second level of the array) parameters
    #
    # Output:
    #
    #       A matrix containing the return period of excesses
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 28 Feb 2006


    rp =excess

    for i in range(1,dim(excess)[1]):
        for j in range(1, dim(excess)[2]):
            if (any( isnan(c(excess[i, j], parameter[1, i, j], parameter[2, i, j])))):
                rp[i, j] =np.nan
            else:
                rp[i, j] =1 / (1-pgpd(excess[i, j], loc=0, scale=parameter[1, i, j], shape=parameter[2, i, j]))

    rp

def tvt(y, span, index, prob):

    # Compute time-varying threshold for a given monthly time series
    #
    #  Description:
    #
    #       Returns smooth long term mean, time-varying threshold and fitted values
    #       given by long term mean plus mean seasonal cycle
    #
    # Usage: svt(y,span,index,prob)
    #
    # Input:
    #
    #       y: vector containing a monthly time series of data
    #
    #    span: fraction of the total number of points of 'y' to be used to compute
    #          the long term mean
    #
    #   index: index vector composed of logicals (True and False) of the same
    #          len as 'y' defining the months to be considered when computing the
    #          time-varying threshold
    #
    #    prob: a value between 0 and 1 that informs the percentage of point
    #          that will be left above the time-varying threshold
    #
    # Outputs
    #
    #    .mfit: fitted data (i.e. long term mean + mean annual cycle)
    #
    #     .ltm: long term mean
    #
    # .theshold: time-varying threshold
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 28 Feb 2006
    #      Chris Ferro <c.a.t.ferro@reading.ac.uk>


    x = np.arange(1,len(y))
    lo =loess(y~x, span = span)

    resid =np.repeat(np.nan, len(y))
    resid[~ isnan(y)] =lo.residuals

    meanresid =np.repeat(np.nan, 12)
    for (i in 1: 12)
        meanresid[i] =mean(resid[seq(i, len(y), 12)], na_rm=True)


    fitted =np.repeat(np.nan, len(y))
    fitted[~ isnan(y)] =lo.fitted

    mfit =fitted + np.repeat(meanresid, len(y) / 12)

    step = sort(y[index] - mfit[index])
    step = unique(step[step > 0])
    th = mfit[index]
    p = mean(y[index] > th, na_rm=True)
    k = 0
    while (p > prob):
        k = k + 1
        th = mfit[index] + step[k]
        p = mean(y[index] > th, na_rm=True)


    dict(threshold=th, ltm=fitted, mfit=mfit)

def xdependence (array3d, ts, u, fun):
    # Compute extreme dependence measures between a given three-dimensional array
    # with first two space dimensions (e.g. longitude and latitude) and third
    # time dimension, and a given time series with the same len as the time
    # dimension of the three dimensional array
    #
    # Description:
    #
    #      Returns a map (matrix) of the extreme dependence measures described in
    #      Coles et al. (1999) and in section 8.4 of Coles (2001)
    #
    # Usage:
    #
    #      xdependence(array3d,ts,u,fun)
    #
    # Input:
    #
    #      array3d: three-dimensional array with p longitude points and q latitude
    #               points as the first two dimensions and n as the third time
    #               dimension
    #
    #          ts: time series of data with the same len as the time dimension
    #              of array3d for which extreme dependence with each grid point of
    #              array3d will be computed
    #
    #           u: high threshold between 0 and 1 that will define the extreme
    #              level to compute dependence. Must be a high quantile (e.g. 0.95)
    #
    #         fun: String of characters that defines which extreme dependence
    #              measure is computed (e.g. "chi", "chibar" or "chicount").
    #              If fun is "chi", computes chi(u) as in Eq.(3.2) of
    #              Coles et al. 1999 and in section 8.4 of Coles (2001).
    #              If fun is "chibar", computes chibar(u) as in section 3.3.2 of
    #              Coles et al. 1999 and in section 8.4 of Coles (2001).
    #              If fun is "chicount", computes the conditional probability of
    #              each grid point variable of array3d to be large (above a large
    #              threshold u) conditional on (given that) the variable ts is
    #              also large (with value above the same large threshold u)
    #              by counting values above threshold and computing relative
    #              frequencies
    #
    # Output:
    #
    #      .out: an map (matrix) of the extreme dependence measure for each grid
    #            point of array3d.
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 21 Dec 2005
    #      Christopher Ferro <c.a.t.ferro@reading.ac.uk>
    #
    # References: Coles, S., Heffernan, J., and Tawn, J., 1999: Dependence measures
    #             of extreme value analyses. Extremes 2:4, 339-365.
    #
    #             Coles, S., 2001: An introduction to statistical modeling of
    #             extreme values. Spring series in statistics. 208pp.
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 1000)
    #      data <- array(rnorm(prod(dim)), dim)
    #      xdependence(data,data[4,5,],u=0.95,fun="chi").out
    #      xdependence(data,data[4,5,],u=0.95,fun="chibar").out
    #      xdependence(data,data[4,5,],u=0.95,fun="chicount").out


    out1 =aperm(apply(array3d, c(1, 2), rank), c(2, 3, 1)) / dim(array3d)[3]

    if (fun == "chi"):
        ts =rank(ts) / len(ts)
        def chi (y):
            2-(log((sum(ts <= u & y <= u)+0.5) / (len(ts)+1)) / (log((sum(ts <= u)+0.5) / (len(ts)+1))))

        out =apply(out1, c(1, 2), chi)
        out[out < 0] =0


    if (fun == "chibar"):
        ts =rank(ts) / len(ts)
        def  chibar (y):
            ((2 * log((sum(ts >= u)+0.5) / (len(ts)+1))) / log((sum(ts >= u & y >= u)+0.5) / (len(ts)+1)))-1

        out = apply(out1, c(1, 2), chibar)


    if (fun == "chicount"):
        def    chicount (y):
            len(y[ts >= quantile(ts, u)][y[ts >= quantile(ts, u)] >= quantile(y, u)]) / len(y[ts >= quantile(ts, u)])

        out =apply(array3d, c(1, 2), chicount)


    dict(out=out)



def xdependence1(lon, lat, array3d, ts, u, fun, nonmissing=0.5):
    # Compute extreme dependence measures between a given three-dimensional array
    # with first two space dimensions (e.g. longitude and latitude) and third
    # time dimension, and a given time series with the same len as the time
    # dimension of the three dimensional array
    #
    # Note: This function allows the user to specify through the parameter
    #       nonmissing (defaul is 0.5) the acceptable percentage of nonmissing
    #       values in the time series of each grid point of the three-dimensional
    #       array. This is the main difference between this function and
    #       xdependence.r
    #
    # Description:
    #
    #      Returns a map (matrix) of the extreme dependence measures described in
    #      Coles et al. (1999) and in section 8.4 of Coles (2001)
    #
    # Usage:
    #
    #      xdependence1(lon,lat,array3d,ts,u,fun,nonmissing)
    #
    # Input:
    #
    #          lon: vector with p longitude values
    #
    #          lat: vector with q latitude values
    #
    #      array3d: three-dimensional array with p longitude points and q latitude
    #               points as the first two dimensions and n as the third time
    #               dimension
    #
    #          ts: time series of data with the same len as the time dimension
    #              of array3d for which extreme dependence with each grid point of
    #              array3d will be computed
    #
    #           u: high threshold between 0 and 1 that will define the extreme
    #              level to compute dependence. Must be a high quantile (e.g. 0.95)
    #
    #         fun: String of characters that defines which extreme dependence
    #              measure is computed (e.g. "chi", "chibar" or "chicount").
    #              If fun is "chi", computes chi(u) as in Eq.(3.2) of
    #              Coles et al. 1999 and in section 8.4 of Coles (2001).
    #              If fun is "chibar", computes chibar(u) as in section 3.3.2 of
    #              Coles et al. 1999 and in section 8.4 of Coles (2001).
    #              If fun is "chicount", computes the conditional probability of
    #              each grid point variable of array3d to be large (above a large
    #              threshold u) conditional on (given that) the variable ts is
    #              also large (with value above the same large threshold u)
    #              by counting values above threshold and computing relative
    #              frequencies
    #
    # nonmissing: Only grid points with fraction given by this fraction
    #            (between 0 and 1) of non-missing values are used to compute
    #            the statistics specified in fun. Default is 0.5.
    #
    # Output:
    #
    #      .out: an map (matrix) of the extreme dependence measure for each grid
    #            point of array3d.
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 21 Dec 2005
    #      Christopher Ferro <c.a.t.ferro@reading.ac.uk>
    #
    # References: Coles, S., Heffernan, J., and Tawn, J., 1999: Dependence measures
    #             of extreme value analyses. Extremes 2:4, 339-365.
    #
    #             Coles, S., 2001: An introduction to statistical modeling of
    #             extreme values. Spring series in statistics. 208pp.
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 1000)
    #      data <- array(rnorm(prod(dim)), dim)
    #      xdependence1(x,y,data,data[4,5,],u=0.95,fun="chi").out
    #      xdependence1(x,y,data,data[4,5,],u=0.95,fun="chibar").out
    #      xdependence1(x,y,data,data[4,5,],u=0.95,fun="chicount").out


    if (all( isnan(ts))):
        raise RuntimeError("all missing values for ts")

    data = reshapefield(lon, lat, array3d).out

    # compute percentage of non-missing values at each grid point
    aux = apply(data, 2, function(x)
    sum(~ isnan(x)) / (len(x)))

    # identify grid points with more than 50% missing values
    index = (1: len(aux))[aux < nonmissing]

    data[, index] =np.nan

    array3d =reshapefield(lon, lat,as.matrix(data)).out
    out1 =aperm(apply(array3d, c(1, 2), rank, na.last = "keep"), c(2, 3, 1))
    aux =apply(array3d, c(1, 2), function(x)
    sum(~ isnan(x)))
    aux1 =array(np.nan, c(dim(out1)))

    for i in range(1, (dim(array3d)[3])):
       aux1[:,:, i] = out1[:,:, i] / aux


    if (fun == "chi"):
        ts =rank(ts, na.last = "keep") / sum(~ isnan(ts))
        def chi(y):
            index = ~( isnan(y) | isnan(ts))
            if (sum(index) == 0):
                return (np.nan)
            2 - (log((sum(ts[index] <= u & y[index] <= u) + 0.5) / (len(ts[index]) + 1)) / (
            return log((sum(ts[index] <= u) + 0.5) / (len(ts[index]) + 1))))

    out =apply(aux1, c(1, 2), chi)
    out[out < 0] =0


    if (fun == "chibar"):
        ts = rank(ts, na.last = "keep") / sum(~ isnan(ts))
        def chibar (y):
            index = ~( isnan(y) | isnan(ts))
            if (sum(index) == 0):
                return (np.nan)
            return ((2 * log((sum(ts[index] >= u) + 0.5) / (len(ts[index]) + 1))) / log(
                (sum(ts[index] >= u & y[index] >= u) + 0.5) / (len(ts[index]) + 1))) - 1

    out =apply(aux1, c(1, 2), chibar)


    if (fun == "chicount"):
        def    chicount (y):
        index = ~( isnan(y) | isnan(ts))
        if (sum(index) == 0):
            return (np.nan)
        len(y[index][ts[index] >= quantile(ts[index], u)][
        return y[index][ts[index] >= quantile(ts[index], u)] >= quantile(y[index], u)]) / len(y[index][ts[index] >= quantile(ts[index], u)])

    out =apply(array3d, c(1, 2), chicount)
    dict(out=out)


def xexcess(array3d, fun='meanexcess', p=0.9, threshold=False, upper=True):
    # Compute mean excess and variance of excess for a given three-dimensional
    # array with first two space dimensions and third time dimension.
    #
    # Description:
    #
    #      Returns a matrix with the mean excess or the variance of excess for a
    #      given quantile or user defined threshold.
    #
    # Usage:
    #
    #      xexcess(array3d,fun,p=0.9,threshold = False, upper=True)
    #
    # Input:
    #
    #      array3d: three-dimensional array with p longitude points and
    #               q latitude points as the first 2 dimensions and n as the
    #               third time dimension
    #
    #      fun: String containing the name of the statistics to be computed
    #           that must be one of the following options "meanexcess",
    #           "varexcess", or "medianexcess".
    #
    #      p: The quantile to be computed (Default is 0.9) or a particular threshold
    #         chosen by the user. If p is a threshold then the argument threshold
    #         (see below) must be set to True.
    #
    #      threshold: Logic argument. Default is False. If True p must be a given
    #                 threshold value.
    #
    #      upper: Logic argument. Defalts is True (examines upper tail of the
    #             distribution).If False lower tail is examined.
    #
    # Output:
    #
    #      .out: a matrix of the computed statistics
    #
    # Authors:
    #
    #      Dag Johan Steinskog <dag.johan.steinskog@nersc.no> 14 June 2005
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #      Christopher Ferro <c.a.t.ferro@reading.ac.uk>
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 100)
    #      data <- array(rnorm(prod(dim)), dim)
    #      xexcess(data,"meanexcess").out
    #      xexcess(data,"varexcess").out
    #      xexcess(data,"meanexcess",p=1.5,threshold =True).out
    #      xexcess(data,"meanexcess",p=1.5,threshold =True,upper=False).out
    #

    if (fun == "meanexcess"):
        def mean_excess (y, p):

            if (threshold):
                u = p
            else:
                u =quantile(y, p, na_rm=True)

            if (upper) mean(y[y > u]-u, na_rm=True)
            else mean(y[y < u]-u, na_rm=True)


    out =apply(array3d, c(1, 2), mean.excess, p)


    if (fun == "varexcess"):
        def  var_excess (y, p):
            if (threshold):
                u = p
            else:
                u =quantile(y, p, na_rm=True)

            if (upper):
                var(y[y > u]-u, na_rm=True)
            else:
                var(y[y < u]-u, na_rm=True)



    out =apply(array3d, c(1, 2), var.excess, p)



    if (fun == "medianexcess"):
        def median_excess (y, p):
            if (threshold):
                u = p
            else:
                u =quantile(y, p, na_rm=True)

            if (upper):
                median(y[y > u]-u, na_rm=True)
            else:
                median(y[y < u]-u, na_rm=True)


    out =apply(array3d, c(1, 2), median_excess, p)


    dict(out=out)



def xgev(array3d, upper=True, n=365):
    # Compute location, shape and scale parameters of a Generalized Extreme Value
    # Distribution for block annual maxima or minima of a given three-dimensional
    # array with first two space dimensions and third time dimension.
    #
    # Description:
    #
    #      Returns an array with the location, shape and scale parameters of
    #      a Generalized Extreme Value Distribution for block annual maxima or
    #      minima of a given three-dimensional array with first two space
    #      dimensions and third time dimension.
    #
    # Usage:
    #
    #      xgev(array3d,upper=True,n=12)
    #
    # Input:
    #
    #      array3d: a 3-dimensional array with p longitude points and q latitude
    #               points as the first 2 dimensions and n as the third time
    #               dimension
    #
    #      upper: Logic argument. Default is True (examines block annual maxima).
    #             If False block annual minima is examined.
    #
    #      n: Length of the time block. Default is 365 (e.g for daily data
    #         annual maxima/minima)
    #

    # Output:
    #
    #      .out: an array of the computed statistics. First dimention contains
    #            the estimated parameters. (e.g. .out[1,,] is the location parameter
    #            .out[2,,] is the scale parameterand and .out[3,,] is the shape
    #            parameter)
    #
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 14 Dec 2005
    #      Christopher Ferro <c.a.t.ferro@reading.ac.uk>
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 1000)
    #      data <- array(rnorm(prod(dim)), dim)
    #      xgev(data)
    #      xgev(data,n=50)
    #      xgev(data,upper=False)
    #


    d =dim(array3d)

    def block_max(x):
        apply(matrix(x, n, d[3] / n), 2, max)


        bmax =aperm(apply(array3d, c(1, 2), block.max), c(2, 3, 1))

    def estim_par(y):
        if (upper):
            fgev(y, std_err = False).est
        else:
            fgev(-y, std_err = False).est



    out =apply(bmax, c(1, 2), estim_par)
    dict(out=out)


def xindex (y, u, threshold=False):
    # Description:
    #
    #   Evaluates the intervals estimator for the extremal index, an index for
    #   time clusters, for a given time series and threshold.
    #
    # Usage:
    #
    #   xindex(y,u,threshold=False)
    #
    # Arguments:
    #
    #   y: a time series for which the extremal index will be computed
    #
    #   u: Threshold. Can be either a quantile (e.g. 0.8) or a threshold value
    #      (see definition of theshold below).
    #
    #   threshold: Logical. If False (the default) u must be a quantile value.
    #              If True, then u must be a threshold value.
    #
    # Value:
    #
    #   Estimate of the extremal index.
    #
    # Details:
    #
    # Warning:
    #
    #   The estimator is not constrained to lie in [0,1] and a
    #   default value of 1 is returned if there are fewer than two
    #   extreme values.
    #
    # Authors:
    #
    #   Christopher Ferro <c.a.t.ferro@reading.ac.uk> 21 Dec 2005
    #   Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #
    # References:
    #
    #   Ferro CAT & Segers J (2003) Inference for clusters of extreme
    #   values. Journal of the Royal Statistical Society B 65, 545-556.
    #
    # See Also:
    #
    # Examples:
    #
    #   y<-rnorm(1000,25,2)
    #   xindex(y,0.9)
    #   xindex(y,28,threshold=True)

    # generates a logical vector indicating which positions correspond to
    # extreme values.
    if (threshold):
        z = (y >= u)
    else:
        z = (y >= quantile(y, u, na_rm=True))

    if (sum(z, na_rm=True) <= 1):
        warning("estimator undefined: too few exceedances")
        return (1)
    else:
        nz = len(z)  # len of sequence
        s = c(1: nz)[z]  # exceedance times
        t = diff(s)  # interexceedance times
        if (max(t, na_rm=True) <= 2):
            t1 = mean(t, na_rm=True)
            t2 = mean(t ^ 2, na_rm=True)
        else:
            t1 = mean(t-1, na_rm=True)
            t2 = mean((t-1) * (t-2), na_rm=True)


    min(1, 2 * (t1 ^ 2) / t2)

def xindexfield (y, u, threshold=False):
    # Description:
    #
    #   Computes the intervals estimator for the extremal index, an index for
    #   time clusters, at each grid point of a three-dimensional array of data
    #   with first two space dimensions (e.g. longitude and latitude) and
    #   third time dimension
    #
    # Usage:
    #
    #   xindexfield(y,u,threshold=False)
    #
    # Arguments:
    #
    #   y: a three-dimensional array, with first two space dimensions (e.g.
    #      longitude and latitude) and third time dimension
    #
    #   u: Threshold. Can be either a quantile (e.g. 0.8) or a threshold value
    #      (see definition of theshold below).
    #
    #   threshold: Logical. If False (the default) u must be a quantile value.
    #              If True, then u must be a threshold value.
    #
    # Value:
    #
    #   Estimate of the extremal index at each grid point.
    #
    # Details:
    #
    # Warning:
    #
    #   The estimator is not constrained to lie in [0,1] and a
    #   default value of 1 is returned if there are fewer than two
    #   extreme values.
    #
    # Authors:
    #
    #   Christopher Ferro <c.a.t.ferro@reading.ac.uk> 3 Jan 2006
    #   Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #
    # References:
    #
    #   Ferro CAT & Segers J (2003) Inference for clusters of extreme
    #   values. Journal of the Royal Statistical Society B 65, 545-556.
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 100)
    #      data <- array(rnorm(prod(dim),25,2), dim)
    #      xindexfield(data,0.9)
    #      xindexfield(data,28,threshold=True)

    apply(y, c(1, 2), xindex, u, threshold)


def xpareto (array3d, p=0.9, upper=True):
    # Compute shape and scale parameters of a Generalized Pareto Distribution for
    # a given three-dimensional array with first two space dimensions and third
    # time dimension.
    #
    # Description:
    #
    #      Returns an array with the shape and scale parameters of a Generalized
    #      Pareto Distribution for the exceedances above a given quantile that
    #      defines the threshold.
    #
    # Usage:
    #
    #      xpareto(array3d,p=0.9, upper=True)
    #
    # Input:
    #
    #      array3d: three-dimensional array with p longitude points and q latitude
    #               points as the first 2 dimensions and n as the third time
    #               dimension
    #
    #      p:     quantile to be computed (Default is 0.9).
    #
    #      upper: Logic argument. Default is True (examines upper tail of the
    #             distribution). If False lower tail is examined.
    #
    # Output:
    #
    #      .out: an array of the computed statistics. First dimension contains
    #            the estimated parameters. (e.g. .out[1,,] is the scale parameter
    #            and .out[2,,] is the shape parameter). Second and third dimensions
    #            are the same space dimensions as of array3d.
    #
    #
    # Authors:
    #
    #      Dag Johan Steinskog <dag.johan.steinskog@nersc.no> 16 June 2005
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk>
    #      Christopher Ferro <c.a.t.ferro@reading.ac.uk>
    #
    # Examples:
    #
    #      x <- seq(-20, 20, 5)
    #      y <- seq(30, 60, 5)
    #      dim <- c(len(x), len(y), 1000)
    #      data <- array(rnorm(prod(dim)), dim)
    #      xpareto(data)
    #      xpareto(data,p=0.95)
    #      xpareto(data,p=0.1,upper=False)
    #


    def estim_par (y, p):
        if (upper):
            fpot(y, quantile(y, p, na_rm=True), "gpd", std_err=False).est
        else:
            fpot(-y, quantile(-y, 1-p, na_rm=True), "gpd", std_err=False).est



    out =apply(array3d, c(1, 2), estim_par, p)
    dict(out=out)


def xparetotvt (lon, lat, array3d, span, p=0.9, index, nonmissing=0.5):

    # Fit Generalized Pareto distribution with time-varying threshold at each grid
    # point for a given three-dimensional array of montly data with first two space
    # dimensions (e.g. longitude and latitude) and third time dimension.
    #
    # Description:
    #
    #      Returns maps of KS test p-values, fitted parameters (constant scale and
    #      constant shape), map of fraction of points above the time-varying
    #      threshold, and the time varying threshold.
    #
    # Usage:
    #
    #      xparetotvt(lon,lat,array3d,span,p,index,nonmissing)
    #
    # Input:
    #
    #       lon: vector with p longitude values
    #
    #       lat: vector with q latitude values
    #
    #   array3d: a three-dimensional array of monthly data with p longitude points
    #            and q latitude points as the first two dimensions and n as the
    #            third time dimension
    #
    #      span: fraction of the total number of points 'n' to be used to compute
    #            the long term mean
    #
    #         p: a value between 0 and 1 that informs the percentage of point
    #            that will be left below the time-varying threshold
    #
    #     index: index vector of len n composed of logicals (True and False)
    #            defining the months to be considered when computing the
    #            time-varying threshold
    #
    # nonmissing: Only grid points with fraction given by 'nonmissing'
    #            (between 0 and 1) of non-missing values are used to estimate the
    #            Generalized Pareto distribution parameters. Default is 0.5.
    #
    # Output:
    #
    # .pvaloutput: p x q map of KS test p-values
    #
    #     .output: 2 x p x q array with scale (first level of the array) and
    #              shape (second level of the array) parameters
    #
    #       .frac: p x q map with fraction of points above the time-varying threshold
    #
    #         .th: p x q x n' array containing the time varying threshold, where
    #              n' is the number of months considered when computing the
    #              time-varying threshold
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 28 Feb 2006
    #      Chris Ferro <c.a.t.ferro@reading.ac.uk>


    # reshape three-dimensional array into a data matrix with
    # time as first dimension and space as sencond dimension
    data = reshapefield(lon, lat, array3d).out

    # compute percentage of non-missing values at each grid point
    aux = apply(data[index,], 2, function(x)
    sum(~ isnan(x)) / (len(x)))

    # identify grid points with less than 50% missing values
    indexgrid = (1: len(aux))[aux >= nonmissing]

    out1 = np.repeat(np.nan, dim(data)[2])
    out2 = np.repeat(np.nan, dim(data)[2])
    pval1 = np.repeat(np.nan, dim(data)[2])
    frac1 = np.repeat(np.nan, dim(data)[2])
    mthresh = matrix(nrow=(dim(data)[1]) / (12 / sum(index[1:12])), ncol=dim(data)[2])

    # function to perform KS test
    def    ff (x, loc, scale, shape):
        if (all( isnan(c(scale, shape)))):
            return (np.nan)

        ks_test(x, "pgpd", loc, scale, shape).p



    # Estimate Generalized Pareto distribution parameters (scale and shape) for
    # exceedances above the threshold p
    def estim_par (y, p, index):
        threshold =tvt(y, span, index, 1 - p).threshold

        yy =y[index]

        fraction =sum(~ isnan(yy[yy > threshold])) / sum(~ isnan(yy))

        params =gpd.fit(yy[~ isnan(yy)], threshold[~ isnan(yy)], show = False).mle

        pval =ff((yy - threshold)[yy > threshold], 0, params[1], params[2])

        list(params=params, pval=pval, fraction=fraction, threshold=threshold)



    out = apply(data[, indexgrid], 2, estim_par, p, index)

    for i in range(1, len(indexgrid)):
        # scale parameter
        out1[indexgrid[i]] =out[[i]].params[1]

    mthresh[, indexgrid[i]] =out[[i]].threshold

    # shape parameter
    out2[indexgrid[i]] =out[[i]].params[2]

    # p-value
    pval1[indexgrid[i]] =out[[i]].pval

    frac1[indexgrid[i]] =out[[i]].fraction



    aux =reshapefield(lon, lat, t(as.matrix(pval1))).out
    pvaloutput =aux[,, 1, drop = T]

    aux1 =reshapefield(lon, lat, t(as.matrix(frac1))).out
    frac =aux1[,, 1, drop = T]

    th =reshapefield(lon, lat, mthresh).out

    output =array(np.nan, c(2, dim(array3d)[1], dim(array3d)[2]))
    output[1,,] =reshapefield(lon, lat, t(as.matrix(out1))).out
    output[2,,] =reshapefield(lon, lat, t(as.matrix(out2))).out

    dict(pvaloutput=pvaloutput, output=output, frac=frac, th=th)


def xparetotvtcov(lon, lat, array3d, span, p=0.9, index, covariates=None, sigl=None, shl=None,
                              siglink=identity, shlink=identity, nonmissing=0.5):
    # Fit Generalized Pareto distribution with time-varying threshold at each grid
    # point for a given three-dimensional array of montly data with first two space
    # dimensions (e.g. longitude and latitude) and third time dimension. Allows
    # linear modelling of the paramters.
    #
    # Description:
    #
    #      Returns fitted parameters.
    #
    # Usage:
    #
    #      xparetotvtcov(lon,lat,array3d,span,p,index,covariates,sigl,shl,siglink,shlink,nonmissing)
    #
    # Input:
    #
    #       lon: vector with p longitude values
    #
    #       lat: vector with q latitude values
    #
    #   array3d: a three-dimensional array of monthly data with p longitude points
    #            and q latitude points as the first two dimensions and n as the
    #            third time dimension
    #
    #      span: fraction of the total number of points 'n' to be used to compute
    #            the long term mean
    #
    #         p: a value between 0 and 1 that informs the percentage of point
    #            that will be left below the time-varying threshold
    #
    #     index: index vector of len n composed of logicals (True and False)
    #            defining the months to be considered when computing the
    #            time-varying threshold
    #
    # covariates: Matrix of covariates for generalized linear modelling of
    #            the parameters. The number of rows should be 'n' (i.e. the same as
    #        the time dimension of 'array3d'.
    #
    # sigl, shl: Numeric vectors of integers, giving the columns of 'ydat'
    #            that contain covariates for generalized linear modelling of
    #            the scale and shape parameters repectively .
    #
    # siglink, shlink: Inverse link functions for generalized linear
    #          modelling of the scale and shape parameters repectively.
    #
    # nonmissing: Only grid points with fraction given by 'nonmissing'
    #            (between 0 and 1) of non-missing values are used to estimate the
    #            Generalized Pareto distribution parameters. Default is 0.5.
    #
    # Output:
    #
    #     .output: n' x p x q array with estimates parameters, where n' is
    #              the total number of estimated parameters
    #
    # Authors:
    #
    #      Caio Coelho <c.a.d.s.coelho@reading.ac.uk> 28 Feb 2006
    #      Chris Ferro <c.a.t.ferro@reading.ac.uk>


    # reshape three-dimensional array into a data matrix with
    # time as first dimension and space as sencond dimension
    data = reshapefield(lon, lat, array3d).out

    # compute percentage of non-missing values at each grid point
    aux = apply(data[index,], 2, function(x)
    sum( np.isfinite(x)) / (len(x)))

    # identify grid points with less than 50% missing values
    indexgrid = np.arange(1, len(aux))[aux >= nonmissing]

    outaux =matrix(ncol=dim(data)[2], nrow=(len(sigl) + len(shl) + 2))

    # Estimate Generalized Pareto distribution parameters (scale and shape) for
    # exceedances above the threshold p
    def estim_par (y, p, index, covariates, sigl, shl, siglink, shlink):

        threshold =tvt(y, span, index, 1 - p).threshold

    yy =y[index]
    covs =as.matrix(covariates[index,])

    index1 = ~( isnan(yy) | apply(t(apply(covs, 1, isnan)), 1, any))

    params =mygpd.fit(yy[index1], threshold[index1], ydat=as_matrix(
        covs[index1,]), sigl = sigl, shl = shl, siglink = siglink, shlink = shlink, show = False).mle

    dict(params=params)



    out =apply(data[:, indexgrid], 2, estim.par, p, index, covariates, sigl, shl, siglink, shlink)

    for j in range(1, (len(sigl) + len(shl) + 2)):
        for i in range(1, len(indexgrid)):
            outaux[j, indexgrid[i]] =out[[i]].params[j]



    output =aperm(reshapefield(lon, lat, outaux).out, c(3, 1, 2))

    dict(output=output)

