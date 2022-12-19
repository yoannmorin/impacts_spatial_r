# Computing impact measures in R

Compute impact measures for SDM models in R, and simulations to compute dipersion measures and significance tests for impact measures.

## Why this code ?
If you're familiar with R and spatial econometrics you know that the SPDEP package provides code to estimate spatial econometrics models and their impacts.

A few years ago there was a case non covered by this package: if you wanted to estimate impacts for a SDM but don't want the spatial lag of each variable (now the "durbin=" argument solved this)


We can write a traditionnal "full" SDM model as :

$$Y = \rho WY + X\beta + W X\delta + \varepsilon$$


The model for which we want to be able to compute impacts measures is :

$$Y = \rho WY + X\beta + W \tilde{X} \gamma + \varepsilon $$

where $\tilde{X}$ is a subset of $X$.



Now in cross-section estimation SPDEP provides a way to estimate impacts for those models, but in the panel context there are still no way to estimates impacts for "full" SDM or "non-full" SDM.
Examples provided here allow the computation of impact measures and z-test for cross-section and panel models.



## Description of the files
There are 5 R files in this repository :
- impacts_sar_cs_compare.R
- impacts_sdm_cs_compare.R
- Impacts_SDM_cs.R
- impacts_sar_panel_compare.R
- impacts_sdm_panel.R


The first two files compare the results obtained using SPDEP package and the code provided here for the SAR and the "full" SDM (in cross-section)
Results are exactly the same.

The third file compute impact measures for "non-full" SDM (in cross-section), and simulations to compute dipersion measures and significance tests for impact measures.

The fourth file compare the results of impacts using the "impacts" R function and my implementation for a panel SAR model, results are exactly the same.

The Fifth file compute impacts measure their standard-deviation z-tests and p-values for a Panel SDM model (a "non-full" SDM but it can be easily adapted).


## How to use this code

- To limit the number of files to download, I used datasets from R packages for the examples.
- Because I'm not an expert in R the code provided needs to be manually edited to be adapted your data (no functions), I don't plan to make it in function forms for now as it would be to time consumming.
- The code should have enough notes for you to be able to understand it and modify it for your data, in case of trouble you can contact me at yoann.morin@protonmail.com
- Note : the code is not as efficient as spdep implementation.



## Methodology

### Computing impact measures
Impact measures are computed following LeSage and Pace (2009), no real difficulty here.

### Computing dispersion measures, z-tests and p-values.
This part was a bit more complicated because there was no clear indication on the methodology to simulates impact measures in LeSage and Pace (2009) or in SPDEP package documentation.

Elhorst (2014) was the most complete explanation I found on the methodology to simulate impact measures.

Steps to follow :

- Estimate parameters of the model (rho, sigma, betas and their variance-covariance matrix)
- Compute impact measures (direct impact, total impact and then indirect impact) following LeSage and Pace (2009)
- Simulates R times parameters values for a normal law centered aronnd the true parameter value, with variance and covariance being specified with the variance-covariance matrix estimated previously.
- Remove rows of simulated observations for rho outside of the [-1;1[ range.
- Compute impacts measures (direct, total and indirect) for each simulated value.
- Now you can compute means and standard deviation for each variable and impact measure.
- Apply a simple t-test/z-test to test for the significativity of the impact measure.



## References

- Impact measures

LeSage, J., & Pace, R. K. (2009). *Introduction to spatial econometrics.* Chapman and Hall/CRC.

- Details to compute dispersion measures and significance tests for impact measures

Elhorst, J. P. (2014). *Matlab software for spatial panels.* International Regional Science Review, 37(3), 389-405.
