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

Let’s now load some custom <code>R</code> functions that are needed to
compute robust estimators.

    source('funciones-auxiliares.R') 
    source('descomponer3.R')         # for covariance decomposition
    source('minimizar.R')            # minimization menu

We will use the Tecator data set available in <code>R</code>
(<http://lib.stat.cmu.edu/datasets/tecator>). Each observation consists
of a spectrometric curve that corresponds to the absorbance measured on
an equally spaced grid of 100 wavelengths between 850 and 1050 nm. The
contents of fat protein and moisture were also recorded through analytic
chemistry methods. The goal of the analysis is to predict the fat
content (*y*) using some characteristics of the spectrometric curve.

    datos <- data(tecator)
    absorp <- tecator$absorp.fdata
    absorp1 <- fdata.deriv(absorp,nderiv = 1)  #  computes the first derivative

Plots of the spectrometric data and it’s first derivative are displayed
below.

![](README_files/figure-markdown_strict/initial%20plots-1.png)

Functional boxplots of the spectrometric data and it’s first derivative
are displayed below.

![](README_files/figure-markdown_strict/functional%20boxplots-1.png)

In particular, we will now consider the first derivative of the
spectrometric curve which we denote *X*.

![](README_files/figure-markdown_strict/X(t)-1.png)

# Robust Estimators under a Functional Linear/Quadratic Model

FLM:
*y* = *α*<sub>0</sub> + ⟨*X*, *β*<sub>0</sub>⟩ + *ϵ*

FQM:
*y* = *α*<sub>0</sub> + ⟨*X*, *β*<sub>0</sub>⟩ + ⟨*X*, *Υ*<sub>0</sub>*X*⟩ + *ϵ*

In both cases, we choose 4 principal directions which explain more than
98% of the total variability as seen below.

    ## [1] 0.9841754

    ## [1] 0.9841754

Residuals boxplots both for FLM and FQM are displayed below.

![](README_files/figure-markdown_strict/robust%20est%20resbox-1.png)

In both cases, we show the common outliers between residuals boxplots
and functional boxplots.

    ## [1] "The outliers from the robust linear fit are:"

    ##  [1]   7  34  35  43  44  45 118 119 121 122 126 127 128 129 130 140 143 168 171
    ## [20] 172 185 186 215

    ## [1] "The common outliers with the functional boxplot are:"

    ## [1]  35 140

    ## [1] "The outliers from the robust quadratic fit are:"

    ##  [1]   4   8  10  20  31  34  35  38  40  43  44  45 102 108 117 121 122 123 125
    ## [20] 126 127 129 131 140 171 172 174 175 183 186 211 215

    ## [1] "The common outliers with the functional boxplot are:"

    ## [1]  35 140

Curves corresponding to atypical residuals from FQM are displayed below
in red.

![](README_files/figure-markdown_strict/robust%20est%20atypical%20curves-1.png)

*β* curve estimated from FLM fit is displayed below.

![](README_files/figure-markdown_strict/robust%20est%20beta%20FLM-1.png)

*β* curve and *υ* surface estimated from FQM fit are displayed below.

![](README_files/figure-markdown_strict/robust%20est%20beta%20gamma%20FQM-1.png)![](README_files/figure-markdown_strict/robust%20est%20beta%20gamma%20FQM-2.png)

Estimates comparison from FLM and FQM.

![](README_files/figure-markdown_strict/robust%20est%20betas-1.png)

![](README_files/figure-markdown_strict/robust%20est%20betas2-1.png)

Residuals vs predicted, both for FLM and FQM are displayed below.

![](README_files/figure-markdown_strict/robust%20est6-1.png)

# Classic Estimators under a Functional Linear/Quadratic Model

# Comparison between Classical and Robust
