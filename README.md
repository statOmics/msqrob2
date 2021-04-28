# R package: msqrob2
#

## Implementation of the MSqRob analysis of differentially expressed proteins using the Features infrastructure

The `msqrob2` package ports and modernises the method presented in
[`MSqRob`](https://github.com/statOmics/MSqRob) and
[`MSqRobSum`](https://github.com/statOmics/MSqRobSum) to use the
[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/articles/Features.html)
class infrastructure.

## Installation

To install the current version of *msqrob2*, run.

```
if(!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager")
}
BiocManager::install("statOmics/msqrob2")
```

The dependencies of the package are listed in the DESCRIPTION file of the package.


## Issues and bug reports

Please use https://github.com/statOmics/msqrob2/issues to submit issues, bug reports, and comments.

## Usage

See vignettes on [msqrob2Examples](https://statomics.github.io/msqrob2Examples)
