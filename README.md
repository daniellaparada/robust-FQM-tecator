The following real data example of the implementation of robust
estimators for functional quadratic regression models. These robust
estimators of the principal directions with robust regression estimators
based on a bounded loss function and a preliminary residual scale
estimator. This is part of a work in progress done in collaboration with
Prof. Dr. Graciela Boente.

Let’s first load some <code>R</code> packages.

    library('fda')          
    library('robustbase')   # to use lmrob 
    library(fda.usc)        # tecator dataset and fda tools
    library(gdata)          # to use uppertriangle
    library(lattice)        # to plot

Let’s now load some custom <code>R</code> functions.

    source('funciones-auxiliares.R') 
    source('descomponer3.R')         # for covariance decomposition
    source('minimizar.R')            # minimization menu

We will use the Tecator data set available in <code>R</code>
(<http://lib.stat.cmu.edu/datasets/tecator>). Each observation consists
of a spectrometric curve that corresponds to the absorbance measured on
an equally spaced grid of 100 wavelengths between 850 and 1050 nm. The
contents of fat protein and moisture were also recorded through analytic
chemistry methods. The goal of the analysis is to predict the fat
content (*y*) using some characteristics of the spectrometric curve. In
particular, we consider the first derivative (*X*) of the spectrometric
curve.

    datos <- data(tecator)
    absorp <- tecator$absorp.fdata
    absorp1 <- fdata.deriv(absorp,nderiv = 1)  #  computes the first derivative

Plots of the spectrometric data and it’s first derivative are displayed
below.

![](README_files/figure-markdown_strict/initial%20plots-1.png)![](README_files/figure-markdown_strict/initial%20plots-2.png)![](README_files/figure-markdown_strict/initial%20plots-3.png)

Functional boxplots of the spectrometric data and it’s first derivative
are displayed below.

![](README_files/figure-markdown_strict/functional%20boxplots-1.png)![](README_files/figure-markdown_strict/functional%20boxplots-2.png)

# Robust Estimator under a Functional Linear/Quadratic Model

FLM:
*y* = *α*<sub>0</sub> + ⟨*X*, *β*<sub>0</sub>⟩ + *ϵ*

FQM:
*y* = *α*<sub>0</sub> + ⟨*X*, *β*<sub>0</sub>⟩ + *ϵ*

In both cases, we choose 4 principal directions which explain more than
97% of the total variability.

    ## [1] 0.9841754

    ## [1] 0.9841754

![](README_files/figure-markdown_strict/robust%20est-1.png)![](README_files/figure-markdown_strict/robust%20est-2.png)

    ## [1] "The outliers from the robust linear fit are:"

    ##  [1]   7  34  35  43  44  45 118 119 121 122 126 127 128 129 130 140 143 168 171
    ## [20] 172 185 186 215

    ## [1] "The common outliers with the functional boxplot are"

    ## [1]  35 140

    ## [1] "The outliers from the robust quadratic fit are:"

    ##  [1]   4   8  10  20  31  34  35  38  40  43  44  45 102 108 117 121 122 123 125
    ## [20] 126 127 129 131 140 171 172 174 175 183 186 211 215

    ## [1] "The common outliers with the functional boxplot are"

    ## [1]  35 140

Curves corresponding to atypical residuals from FQM.

![](README_files/figure-markdown_strict/robust%20est2-1.png)

Residuals boxplots from FLM and FQM.
![](README_files/figure-markdown_strict/robust%20est3-1.png)![](README_files/figure-markdown_strict/robust%20est3-2.png)

Beta from FLM fit.
![](README_files/figure-markdown_strict/robust%20est4-1.png)![](README_files/figure-markdown_strict/robust%20est4-2.png)

    ## [1] -403

![](README_files/figure-markdown_strict/robust%20est4-3.png)![](README_files/figure-markdown_strict/robust%20est4-4.png)![](README_files/figure-markdown_strict/robust%20est4-5.png)![](README_files/figure-markdown_strict/robust%20est4-6.png)![](README_files/figure-markdown_strict/robust%20est4-7.png)

Beta and gamma from FQM fit.
![](README_files/figure-markdown_strict/robust%20est5-1.png)

    ## [1] -403

![](README_files/figure-markdown_strict/robust%20est5-2.png)

Residuals vs predicted.
![](README_files/figure-markdown_strict/robust%20est6-1.png)![](README_files/figure-markdown_strict/robust%20est6-2.png)![](README_files/figure-markdown_strict/robust%20est6-3.png)![](README_files/figure-markdown_strict/robust%20est6-4.png)
