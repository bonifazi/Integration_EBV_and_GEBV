# Functions for integrating international (G)EBVs at the national level

## Description
This repository contains `R` functions for integrating (genomic) estimated breeding values ((G)EBVs) into national evaluations following [Bonifazi et al., 2023, GSE](https://doi.org/10.1186/s12711-023-00813-2). These functions compute Effective Records Contribution (ERC), de-regressed proof (DRP), adjusted de-regressed ERC (dERC*) and adjusted DRP (DRP*).

## Requirements
Two R packages are required: [`dplyr`](https://dplyr.tidyverse.org/) (for data manipulation) and [`assertthat`](https://github.com/hadley/assertthat) (for asserting arguments' type). For installation, please refer to the respective packages' links.

## Installation
Either download [All_functions.R](https://github.com/bonifazi/Integration_EBV_and_GEBV/blob/main/All_functions.R) or copy-paste the code contained in it. This file contains all R functions.  
Load the R functions into your R session either by including them in your scripts (copy-paste and run them before using them) or by sourcing the downloaded R file into your R code as
```r
source("my_directory/All_functions.R")
```
## Usage and documentation
See the documentation of each R function for a description of its input, output, and usage. To view the help documentation of any function, use the [`docstring`](https://github.com/Dasonk/docstring) R package as `docstring(fun = <function_name>)`. For instance:
```r
#install.packages("docstring") # install docstring pkg
library(docstring) # load docstring pkg
# view the documentation associated with the function compute_ERC_and_DRP
docstring(fun = compute_ERC_and_DRP)
```

## Licence
This project is licensed under the MIT licence - see the [LICENSE](https://github.com/bonifazi/Integration_EBV_and_GEBV/blob/main/LICENSE) file for details.

## Citation
If you use this code in your research or find it helpful, please consider citing our paper:
_Bonifazi, R., Calus, M.P.L., ten Napel, J. et al. Integration of beef cattle international pedigree and genomic estimated breeding values into national evaluations, with an application to the Italian Limousin population. Genet Sel Evol 55, 41 (2023)._ https://doi.org/10.1186/s12711-023-00813-2

```bibtex
@article{Bonifazi2023,
  title={Integration of beef cattle international pedigree and genomic estimated breeding values into national evaluations, with an application to the Italian Limousin population},
  author={Bonifazi, R. and Calus, M.P.L. and ten Napel, J. and others},
  journal={Genetics Selection Evolution},
  year={2023},
  volume={55},
  pages={41},
  doi={10.1186/s12711-023-00813-2}
}
```

## Contact
For questions or support, please contact renzo.bonifazi@wur.nl or renzo.bonifazi@outlook.it.

### Schematic overview of the integration procedure
<div style="text-align:center">
  <p align="center">
    <img src="https://github.com/bonifazi/Integration_EBV_and_GEBV/assets/74569672/d9bf25b7-8d2e-4b89-877d-c26979f97433" alt="Schematic overview of the integration procedure">
  </p>
   <p align="center">
    Figure 1: Schematic overview of the procedure for the integration of international (G)EBVs at the national level.<sup>1</sup>
</div>
<sup>1: Figure from Bonifazi et al., 2023, GSE. See the publication for more details.</sup>
