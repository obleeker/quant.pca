# Principal Quantile Components Analysis
## Introduction
Principal Quantile Components Analysis is an analogue to PCA in an asymmetric L1 norm. 
For the theory behind see Tran and Osipenko (2016).
This method is based on the robust quantile matrix factorization algorithm developed by Zhu et al. (2017).
It approaches PCA in an asymmetric linear norm via smooth approximation of the objective function.
The objective function is regularized and sequentially optimized via quadratic programming. 

Setting this method to the symmetric case (`tau=0.5`) this is equal to performing L1-PCA. 

This implementation can be used in conjunction with the functional data analysis (`fda`) package for analyis of multivariate time series.   

The code was developed as part of this authors master thesis. 

## Installation
You can install quant.pca from github with:
```
devtools::install_github("obleeker/qant.pca")
```

## Example

```
# generate data
n = 100
X = data.frame(cbind(rnorm(n),rnorm(n),rnorm(n)))

# running PQC on the 0.9-quantile 
pqcomp(data = X, projDim = 2, tau =  0.9, lambda = 0.1, muEst = T)
```
Optimal selection of the regularization `lambda` can be found via the supplied cross-validation method (`cv.pqc()`).

Notes: Parallel execution is not currently supported on Windows.

## References 
[1]: https://ssrn.com/abstract=2854822  "Tran et. al. 2016"
[2]: https://ieeexplore.ieee.org/document/8057117 "Zhu et. al. 2017"

(1) [Tran, N. M., Burdejová, P., Ospienko, M. and Härdle, W. K. (2019), Principal component
analysis in an asymmetric norm.][1]

(2) [Zhu, R., Niu, D. and Li, Z. (2017), Robust web service recommendation via quantile matrix
factorization.][2]


