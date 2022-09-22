The following is a real data example of the implementation of robust
estimators for functional quadratic regression models, compared to the
least squares approach to which will refer as “classic”. This robust
proposal involves robust estimators of the principal directions with
robust regression estimators based on a bounded loss function and a
preliminary residual scale estimator. This is part of a work in progress
done in collaboration with Prof. Dra. Graciela Boente.

Let’s first load some <code>R</code> packages.

    library('fda')          # fda tools
    library('robustbase')   # lmrob 
    library(fda.usc)        # tecator dataset and fda tools
    library(gdata)          # uppertriangle
    library(lattice)        # plot

Let’s now load some custom <code>R</code> functions that are needed to
compute robust estimators.

    source('funciones-auxiliares.R') # aux functions
    source('descomponer3.R')         # for covariance decomposition
    source('minimizar.R')            # minimization menu

We will use the Tecator data set available in the <code>fda.usc</code>
library from <code>R</code>
(<http://lib.stat.cmu.edu/datasets/tecator>). Each observation from this
dataset consists of a spectrometric curve that corresponds to the
absorbance measured on an equally spaced grid of 100 wavelengths between
850 and 1050 nm. The contents of fat protein and moisture were also
recorded through analytic chemistry methods. The goal of the analysis is
to predict the fat content (*y*) using some characteristics of the
spectrometric curve.

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

For the functional regression model, we consider two relationships
between the functional explanatory variable *X* and the scalar response
*y*.

Linear (FLM):
*y* = *α*<sub>0</sub> + ⟨*X*, *β*<sub>0</sub>⟩ + *ϵ*

Quadratic (FQM):
*y* = *α*<sub>0</sub> + ⟨*X*, *β*<sub>0</sub>⟩ + ⟨*X*, *Υ*<sub>0</sub>*X*⟩ + *ϵ*

In both cases, we choose 4 principal directions which explain more than
98% of the total variability as seen below.

    ## [1] 0.9841754

    ## [1] 0.9841754

Residuals boxplots both for FLM and FQM are displayed below.

![](README_files/figure-markdown_strict/robust%20est%20resbox-1.png)

In both cases, we show the common outliers between residuals boxplots
and functional boxplots. As seen below, only observations 35 and 140 are
common outliers with the functional boxplot for both linear (FLM) and
quadratic (FQM) fit.

    ## [1] "The outliers from the robust FLM fit are:"

    ##  [1]   7  34  35  43  44  45 118 119 121 122 126 127 128 129 130 140 143 168 171
    ## [20] 172 185 186 215

    ## [1] "The common outliers with the functional boxplot are:"

    ## [1]  35 140

    ## [1] "The outliers from the robust FQM fit are:"

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

Estimates comparison *β* from FLM and FQM in the robust fit. The shape
is similar, but for the FQM robust fit, the variation range of *β* in
*y*−axis is wider.

![](README_files/figure-markdown_strict/robust%20est%20betas-1.png)

![](README_files/figure-markdown_strict/robust%20est%20betas2-1.png)

Residuals vs predicted, both for FLM and FQM in the robust fit are
displayed below.

![](README_files/figure-markdown_strict/robust%20est6-1.png)

# Classic Estimators under a Functional Linear/Quadratic Model and comparison between their robust counterparts

In both cases, FLM and FQM, we choose 4 principal directions for the
classic fit which explain more than 97% of the total variability as seen
below.

    ## [1] 0.9791543

*β* curve estimated from FLM classic fit is displayed below.

![](README_files/figure-markdown_strict/classic%20est%20beta%20FLM-1.png)

We can now compare *β* estimates in FLM between classic and robust fit.

![](README_files/figure-markdown_strict/robvsclas%20est%20beta%20FLM-1.png)
*β* curve and *υ* surface estimated from FQM classic fit are displayed
below.

![](README_files/figure-markdown_strict/classic%20est%20beta%20gamma%20FQM-1.png)![](README_files/figure-markdown_strict/classic%20est%20beta%20gamma%20FQM-2.png)

We can now compare *β* estimates in FQM between classic (CL) and robust
(ROB) fit.

![](README_files/figure-markdown_strict/robvsclas%20est%20beta%20FQM-1.png)
Residuals vs predicted, both for FLM and FQM in the classic fit are
displayed below.

![](README_files/figure-markdown_strict/classic%20est6-1.png)

# Analysis without outliers

We now repeat the classical analysis, both for FLM and FQM, without the
outliers that have been identified in the residual boxplot from the
robust FQM fit. We will refer it as CL woithout out. We choose 4
principal directions which explain more than 97% of the variability as
seen below.

    ## [1] 0.9774696

Predicted vs residuals in the classic (CL) FLM and FQM fit, both without
outliers, are shown in the first pair of plots. Estimated *β* parameters
are compared among classic (CL), robust (ROB) and classic after outlier
removal (CL without out) in the second pair of plots, both for FLM and
FQM fit.

![](README_files/figure-markdown_strict/plotscompar-1.png)![](README_files/figure-markdown_strict/plotscompar-2.png)

We now show surface differences between classic and robust fit, and
between classic without outliers and robust fit, both for FQM. As seen
below, robust fit behaves similarly to the classical one when atypical
observations are removed.

![](README_files/figure-markdown_strict/dif%20suf%20CL%20ROB-1.png)![](README_files/figure-markdown_strict/dif%20suf%20CL%20ROB-2.png)![](README_files/figure-markdown_strict/dif%20suf%20CL%20ROB-3.png)
