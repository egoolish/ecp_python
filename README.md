Overview
--------------------------------------------------

**ecp** is the Python equivalent of the ecp R package containing various procedures for finding multiple change-points. Two methods make use of dynamic programming and pruning, with no distributional assumptions other than the existence of certain absolute moments in one method. Hierarchical and exact search methods are included. All methods return the set of estimated change-points as well as other summary information.

Please see [package manual](https://cran.r-project.org/web/packages/ecp/ecp.pdf) for more details.


Datasets
--------------------------------------------------

The package contains two datasets
- ACGH: bladder tumor micro-array data
- DJIA: Dow Jones industrial average index


Change-point methods
--------------------------------------------------

The package contains algorithms for detecting multiple change-points. A list containing brief descriptives of these algorithms are as follows:
- Agglomerative hierarchical estimation
- Divisive hierarchical estimation
- Estimation by pruned objectives via energy statistics
- Estimation by pruned objectives via Kolmogorov-Smirnov statistics



Installation
--------------------------------------------------

``` python
pip install git+https://github.com/egoolish/ecp_python
```

Example
--------------------------------------------------

The example.py file contains various basic examples demonstrating how to use the ecp package, as well as the similarities to the R version.
Note that this file does require the package **rpy2** as well as its dependencies. 