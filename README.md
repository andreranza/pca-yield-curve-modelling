Yield Curve Modelling with PCA for Market Risk assessment
================
Andrea Ranzato
2020-04-06

## Abstract

This dissertation illustrates how *principal components* computed on a
time series of US yield curves identify the most common kind of
movements occurring in interest rates of different maturities. In
addition, it shows that the significant PCs can be employed to build a
*risk factor model* for an interest rate sensitive portfolio comprised
of US Government bonds which might be used to achive interest risk
immunization. Conversely to traditional price sensitivity measures of
bond securities, such as *duration*, these models take into account *non
parallel shifts* of the yield curve, usually known as *tilt* and
*curvature*. Furthermore, we investigate the empirical distribution of
the eigenvalues and eigenvectors resulting from the singular value
decomposition of the centered and scaled absolute interest rate changes
using the *bootstrap*.

## Objectives

1.  Understanding the most common movements of the U.S. yield curve from
    a multivariate time series of interest rates of different
    maturities.
2.  Building a linear risk factor model using the most significant
    principal components for a portfolio of U.S. Government bonds.
3.  Obtaining a “probabilistic” view of the obtained results using the
    bootstrap.

## Seminal paper: Common Factors Affecting Bond Returns (1991)

  - Bond portfolios are affected by yield curve risk.
  - Yield curve risk: movements of interest rates of different
    maturities.
  - How to manage *multidimensional* interest rate risk?
  - Traditional duration analysis identifies only *parallel shifts* of
    the yield curve.
  - Statistical approach can identify not only parallel shifts, but also
    *tilts* and *curvature* movements.

![seminal-paper](/Users/andrearanzato/github/pca-yield-curve-modelling/report/report-thesis_files/figure-gfm/cover-seminal-paper.png){width=50%}

## PCA in Market Risk Analysis

### The Yield Curve

A spot interest rate is the interest rate having validity from now, time
zero, until time *t* in the future.

The concept of the *spot yield curve* is introduced since the assumption
of a flat rate used to discount cash flows at different maturity is not
realistic. Conversely, one can observe patterns in the required yields
as a function of the maturity of the bond. Usually, the higher the
maturity of a bond, the higher is the expected return reflected by the
yield as it is possible to see in the figure below. However, in certain
circumstances this might not be case.

For the purpose of the current analysis, it will be relevant not much
the particular shape assumed by the yield curve on a certain day,
instead, it will be of interest understanding its absolute *change*
between two consecutive points in time. The notion of *yield curve
change* is discussed extensively by Jones (1991) who discusses how the
knowledge of the most common kind of movements having occurred
historically proves useful not only for interest risk analysis, but also
to figurate out portfolio allocations that might capitalize on such
movements.

<div class="figure" style="text-align: center">

<img src="/Users/andrearanzato/github/pca-yield-curve-modelling/report/report-thesis_files/figure-gfm/yieldcurve-1.png" alt="Spot Yield Curve on March 10, 2020"  />

<p class="caption">

Spot Yield Curve on March 10, 2020

</p>

</div>

### PCA on the Interest Rates Term Structure

Fixed-income portfolios are exposed to interest rate risk. In
particular, financial institutions are willing to quantify the exposure
of bond portfolios to *unequal* fluctuations in interest rates of
different maturities. However, standard measures of bond price
volatility\[1\], such as duration and convexity, do not provide a good
estimate of such changes. In fact, their explanation power is limited
to:

1.  variations of the required yield, disregarding the existence of
    interest rates of different maturities.
2.  *parallel* shifts of the yield.

In other words, these measures rely on the simplifying assumption of a
*flat* yield curve (Fabozzi 2013).

For instance, “consider the situation of a U.S. government bond trader.
The trader’s portfolio is likely to consist of many bonds with different
maturities. \[Hence\], there is an exposure to movements in the one-year
rate, the two-year rate, the three-year rate, and so on. \[…\] He or she
must be concerned with all the different ways in which the U.S. Treasury
yield curve changes its shape through time” (Hull 2018, 185). As a
result, bonds portfolio are sensitive to several sources of risk,
potentially as many interest rates impact its value. For this reason,
there is the need to reduce the number comprising the set of risk
factors to a manageable one, usually three or four.

Hence, a statistical approach\[2\] requires to collect historical
observations (see fig. below) of interest rates at different maturities
to infer what the future changes might be on the basis of those observed
in the past\[3\]. This would imply to estimate a model which we require
to achieve the following objectives:

  - Identifying the fundamental yield curve movements.  
  - Being able to replicate as faithfully as possible the covariation of
    the overall system of interest rates in a compact way.

The first task is accomplished by *principal component analysis*,
whereas the second by a *linear factor model*\[4\], whose factors are
selected among the PCs with higher explanatory power. Hence, those
principal components represent a reduced set of “basis” that combined
linearly with the eigenvectors are able to replicate with accuracy the
past behaviour of the interest rates in a compact manner, and more
importantly to identify those common movements.

<div class="figure" style="text-align: center">

<img src="report-thesis_files/figure-gfm/termstructure-1.png" alt="Figure: US Treasury Interest Rates 2006-2020"  />

<p class="caption">

Figure: US Treasury Interest Rates 2006-2020

</p>

</div>

| DATE       | MAT1MO | MAT3MO | MAT6MO | MAT1YR | MAT2YR | MAT3YR | MAT5YR | MAT7YR | MAT10YR | MAT20YR | MAT30YR |
| :--------- | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | ------: | ------: | ------: |
| 2006-02-09 |   4,32 |   4,52 |   4,67 |   4,66 |   4,66 |   4,62 |   4,55 |   4,55 |    4,54 |    4,72 |    4,51 |
| 2006-02-10 |   4,36 |   4,53 |   4,70 |   4,70 |   4,69 |   4,67 |   4,59 |   4,59 |    4,59 |    4,76 |    4,55 |
| 2006-02-13 |   4,38 |   4,55 |   4,71 |   4,70 |   4,68 |   4,66 |   4,58 |   4,58 |    4,58 |    4,76 |    4,56 |
| 2006-02-14 |   4,42 |   4,55 |   4,72 |   4,71 |   4,69 |   4,68 |   4,61 |   4,61 |    4,62 |    4,80 |    4,60 |
| 2006-02-15 |   4,39 |   4,55 |   4,70 |   4,70 |   4,71 |   4,68 |   4,60 |   4,60 |    4,61 |    4,78 |    4,58 |
| 2006-02-16 |   4,38 |   4,55 |   4,69 |   4,69 |   4,69 |   4,67 |   4,59 |   4,59 |    4,59 |    4,77 |    4,57 |

Table: Sample of the data set of the US Term Structure.

### Singular Value Decomposition of the term structure

To illustrate how PCA can be used as a dimension reduction tool for
interest risk analysis, we consider the multivariate time series of
daily [U.S. spot
rates](https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield)
for a period of 3521 trading days, from Feb 9, 2006 to Mar 10, 2020.

We observe that the term structure is highly correlated, in particular,
among rates of close maturity. The figures in the table below confirm
this fact, since the correlation between yields is higher around the
main diagonal, meaning that interest rates of similar maturity move
closely. Conversely, the correlation is weaker between interest rates of
different maturity.

|         | MAT1MO | MAT3MO | MAT6MO | MAT1YR | MAT2YR | MAT3YR | MAT5YR | MAT7YR | MAT10YR | MAT20YR | MAT30YR |
| :------ | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | ------: | ------: | ------: |
| MAT1MO  |  1,000 |  0,998 |  0,994 |  0,988 |  0,972 |  0,950 |  0,885 |  0,806 |   0,709 |   0,541 |   0,481 |
| MAT3MO  |  0,998 |  1,000 |  0,998 |  0,994 |  0,980 |  0,959 |  0,894 |  0,816 |   0,718 |   0,549 |   0,487 |
| MAT6MO  |  0,994 |  0,998 |  1,000 |  0,998 |  0,987 |  0,968 |  0,906 |  0,828 |   0,731 |   0,561 |   0,496 |
| MAT1YR  |  0,988 |  0,994 |  0,998 |  1,000 |  0,994 |  0,977 |  0,918 |  0,841 |   0,743 |   0,571 |   0,505 |
| MAT2YR  |  0,972 |  0,980 |  0,987 |  0,994 |  1,000 |  0,994 |  0,950 |  0,882 |   0,785 |   0,615 |   0,547 |
| MAT3YR  |  0,950 |  0,959 |  0,968 |  0,977 |  0,994 |  1,000 |  0,976 |  0,923 |   0,835 |   0,673 |   0,606 |
| MAT5YR  |  0,885 |  0,894 |  0,906 |  0,918 |  0,950 |  0,976 |  1,000 |  0,983 |   0,927 |   0,799 |   0,743 |
| MAT7YR  |  0,806 |  0,816 |  0,828 |  0,841 |  0,882 |  0,923 |  0,983 |  1,000 |   0,977 |   0,889 |   0,845 |
| MAT10YR |  0,709 |  0,718 |  0,731 |  0,743 |  0,785 |  0,835 |  0,927 |  0,977 |   1,000 |   0,964 |   0,935 |
| MAT20YR |  0,541 |  0,549 |  0,561 |  0,571 |  0,615 |  0,673 |  0,799 |  0,889 |   0,964 |   1,000 |   0,989 |
| MAT30YR |  0,481 |  0,487 |  0,496 |  0,505 |  0,547 |  0,606 |  0,743 |  0,845 |   0,935 |   0,989 |   1,000 |

Table: Correlation Matrix of the US Term Structure.

As a result, figures in the table above demonstrate that “treasury
yields do not move around in a completely uncorrelated fashion. If they
did, it would be impossible to analyze the interest rates risk of a bond
portfolio in any meaningful way; even the notion of portfolio duration
would be meaningless” (Fabozzi 2012, 797).

Hence, to assess the yield curve risk affecting a fixed-income
portfolio, it would be useful to understand, at least historically, what
are the most common types of shifts that have occurred in the yield
curve.

As stated previously, principal component analysis performed on the
*interest rates changes* is capable of detecting them, in the form of
principal components. Usually\[5\], the first principal component
records an almost *parallel shift* of the yield curve, the second one a
change in the slope (*tilt*), and the third one a change located in the
middle of the term structure (*curvature* or *convexity*). The first
degree of intuition for this representation is provided by the figure
below which illustrates the first three *eigenvectors* resulting from
the singular value decomposition of the interest rates changes. In red,
the first eigenvector is approximately a parallel line since it takes
similar values across the entire spectrum of maturities. For this
reason, it captures parallel movements of the yield curve. Subsequently,
the eigenvector in green is almost increasing, hence it explains
movements which are downward in nature on early maturities and upward on
later ones. Ultimately, the third eigenvector in blue, is decreasing at
the beginning and increasing at the end. Therefore, it describes
inverted “bumps” of the yield curve (Alexander 2008a)\[6\].

<div class="figure" style="text-align: center">

<img src="report-thesis_files/figure-gfm/eigenvectors-1.png" alt="Eigenvectors"  />

<p class="caption">

Eigenvectors

</p>

</div>

The first table below shows the actual values of the eigenvectors, also
known as *loadings*, whereas the second one illustrates the
corresponding eigenvalues in decreasing order of explanatory power.
Since the correlation matrix is positive definite, the eigenvalues are
all positive. In addition, we know that each principal component
contributes to explain an amount of variance corresponding to its
associated eigenvalue. For instance, the first principal component PC1
explains 62.97% of the total covariation between changes in the US
interest rates. Instead, the first three PCs, considered jointly,
capture 91.05% of the total covariation in the system.

| Maturity |      w1 |      w2 |      w3 |      w4 |      w5 |
| :------- | ------: | ------: | ------: | ------: | ------: |
| MAT1MO   | \-0,123 | \-0,454 |   0,537 |   0,625 |   0,291 |
| MAT3MO   | \-0,164 | \-0,514 |   0,239 | \-0,216 | \-0,722 |
| MAT6MO   | \-0,221 | \-0,451 | \-0,126 | \-0,495 |   0,217 |
| MAT1YR   | \-0,283 | \-0,336 | \-0,302 | \-0,144 |   0,481 |
| MAT2YR   | \-0,335 | \-0,045 | \-0,412 |   0,311 | \-0,078 |
| MAT3YR   | \-0,354 |   0,026 | \-0,306 |   0,256 | \-0,149 |

Table: Slice of the first five Eigenvectors with corresponding loadings

| Eigenvalue Id. | Value |     % | Cumulative |
| -------------: | ----: | ----: | ---------: |
|              1 | 2,632 | 0,630 |      0,630 |
|              2 | 1,534 | 0,214 |      0,844 |
|              3 | 0,858 | 0,067 |      0,911 |
|              4 | 0,673 | 0,041 |      0,952 |
|              5 | 0,479 | 0,021 |      0,973 |

Table: First five eigenvalues.

<img src="report-thesis_files/figure-gfm/eigenvalues-histogram-1.png" style="display: block; margin: auto;" />

Once the eigenvalues and eigenvectors are found using the singular value
decomposition, the principal components are computed using equation
![\\underset{(3520\\times 11)}Z=\\underset{(3520\\times 11)}X\\underset{(11\\times 11)}A](https://latex.codecogs.com/png.latex?%5Cunderset%7B%283520%5Ctimes%2011%29%7DZ%3D%5Cunderset%7B%283520%5Ctimes%2011%29%7DX%5Cunderset%7B%2811%5Ctimes%2011%29%7DA
"\\underset{(3520\\times 11)}Z=\\underset{(3520\\times 11)}X\\underset{(11\\times 11)}A"),
in which:

  - ![X](https://latex.codecogs.com/png.latex?X "X") contains
    ![p](https://latex.codecogs.com/png.latex?p "p") time series of
    interest rates changes measured in basis points at the 11 different
    maturities
  - ![A](https://latex.codecogs.com/png.latex?A "A") is the orthogonal
    matrix of eigenvectors
  - ![Z](https://latex.codecogs.com/png.latex?Z "Z") collects the
    resulting variables of the transformation performed by matrix
    ![A](https://latex.codecogs.com/png.latex?A "A") on the original set
    of variables in ![X](https://latex.codecogs.com/png.latex?X "X")

In other words, each *k-th* principal component in
![Z](https://latex.codecogs.com/png.latex?Z "Z"), is produced by
“weighting” *all* the columns in
![X](https://latex.codecogs.com/png.latex?X "X") with the coefficients
in the *k-th* column vector of
![A](https://latex.codecogs.com/png.latex?A "A"). For this reason, each
of them is a *linear combination* of the original variables in
![X](https://latex.codecogs.com/png.latex?X "X").

The role of principal components analysis consists in finding those
weights (the columns of ![A](https://latex.codecogs.com/png.latex?A
"A")), that lead the derived variables in
![Z](https://latex.codecogs.com/png.latex?Z "Z") to have desirable
properties:

  - highest possible variance
  - zero correlation between each others

The figure below shows a small slice of the *first three* PCs which are
the result of the weighted sum of the interest rates changes with the
first three eigenvectors
![w\_{1}](https://latex.codecogs.com/png.latex?w_%7B1%7D "w_{1}"),
![w\_{2}](https://latex.codecogs.com/png.latex?w_%7B2%7D "w_{2}") and
![w\_{3}](https://latex.codecogs.com/png.latex?w_%7B3%7D "w_{3}") as
weights. Therefore, we might say informally that PC1 “incorporates” the
“information” of the first eigenvector which captures a parallel shift
of the yield, and so on.

<img src="report-thesis_files/figure-gfm/three-pcs-graph-1.png" style="display: block; margin: auto;" />

Furthermore, it is useful to recall that the PCs are *orthogonal*, thus
their correlation is zero. This property is particularly useful in
market risk models because it allows to handle *uncorrelated* risk
factors.

In conclusion: “the first principal component captures a *common trend*
\[…\] \[in\] interest rates \[changes\]. That is, if the first principal
component changes at a time when the other components are fixed, then
\[the interest rates\] all move by roughly the same amount. For this
reason we often called the first component the *trend component*. \[…\]
Then the second principal component usually captures a change in slope
of the term structure. \[…\] For this reason we often called the third
component the *curvature* or *convexity* component” (Alexander 2008a).

### Linear Factor Model

Factor models are conceptually independent from principal component
analysis and they constitute an independent field of research in
statistics. Nevertheless, PCA provides a feasible estimation strategy
for their “factors” . This should not be a surprise, given that
“analysis of principal components are more of a means to an end rather
than an end in themselves, because they frequently serve as intermediate
steps in much larger investigations” (Johnson and Wichern 2014).

We require from a linear factor model\[7\] to approximate with the
highest possible accuracy the covariation of the observed interest rates
changes, using just a small set of risk factors.

In particular, we approximate the random vector of interest rates
changes   
![({{\\Delta}{{R}^{(m1)}}}, {{\\Delta}{{R}^{(m3)}}}, \\ldots,
{{\\Delta}{{R}^{(y30)}}})](https://latex.codecogs.com/png.latex?%28%7B%7B%5CDelta%7D%7B%7BR%7D%5E%7B%28m1%29%7D%7D%7D%2C%20%7B%7B%5CDelta%7D%7B%7BR%7D%5E%7B%28m3%29%7D%7D%7D%2C%20%5Cldots%2C%20%7B%7B%5CDelta%7D%7B%7BR%7D%5E%7B%28y30%29%7D%7D%7D%29
"({{\\Delta}{{R}^{(m1)}}}, {{\\Delta}{{R}^{(m3)}}}, \\ldots, {{\\Delta}{{R}^{(y30)}}})")  
In particular, we approximate the actual interest rates changes by means
of three risk factors represented by the first three principal
components, and factor weights the corresponding eigenvectors. We call
*principal component representation* at time *t* the daily interest
rates changes provided by the following linear factor model:   
![\\underset{(3520 \\times 11)}{{\\Delta}{R}}\\approx\\underset{(n
\\times 3)}{Z}\\underset{(3 \\times 11)}{W^{T}}
](https://latex.codecogs.com/png.latex?%5Cunderset%7B%283520%20%5Ctimes%2011%29%7D%7B%7B%5CDelta%7D%7BR%7D%7D%5Capprox%5Cunderset%7B%28n%20%5Ctimes%203%29%7D%7BZ%7D%5Cunderset%7B%283%20%5Ctimes%2011%29%7D%7BW%5E%7BT%7D%7D%20
"\\underset{(3520 \\times 11)}{{\\Delta}{R}}\\approx\\underset{(n \\times 3)}{Z}\\underset{(3 \\times 11)}{W^{T}} ")  

We shall see the above matrix multiplication as follows: each *m-th*
column of
![\\Delta{R}](https://latex.codecogs.com/png.latex?%5CDelta%7BR%7D
"\\Delta{R}") is computed by “weighting” *all* the three PCs in
![Z](https://latex.codecogs.com/png.latex?Z "Z"), with the coefficients
of the *m-th column* of
![W^T](https://latex.codecogs.com/png.latex?W%5ET "W^T").

The linear factor model estimated on our data is able to explain about
91.05% of the total covariation of the interest rates changes given that
we have used only the first three components. A better approximation can
always be achieved by adding more PCs at the cost of increasing the
dimensionality of the model. Still, the model demonstrates the power of
principal component analysis, since we are capable of explaining almost
the entirety of the variance in the system exploiting only three
components, instead of the original entire set of variables. This model
will be employed later to describe the *profit and loss* of a portfolio
composed of U.S. Government bonds.

The following two plots show the actual interest rates changes and the
principal component approximation. On the top, it is shown the linear
factor model estimated by means of principal component which makes a
pretty good job in replicating the overall covariation of the actual
system depicted on the bottom.

<img src="report-thesis_files/figure-gfm/representation-vs-actual-1.png" style="display: block; margin: auto;" /><img src="report-thesis_files/figure-gfm/representation-vs-actual-2.png" style="display: block; margin: auto;" />

## Bootstrap: Eigenvalues and Eigenvectors

In this section, we apply the *bootstrap* to approximate the theoretical
distribution of the eigenvalues and eigenvector loadings. In other
words, we would like to have a *probabilistic* representation of the
results obtained in previous sections. However, an important remark is
needed from the very first. The following figures should be interpreted
taking in consideration that each bootstrap sample is truly a *random
sample*\[8\], whereas the figures obtained in the previous section are
affected by either the *autocorrelation* existing within a single
interest rate and the *cross-correlation* subsisting among them. In
particular, the proportion of variance explained by the first three
components obtained with the bootstrap will substantially outperform
those obtained considering consecutive time frames. In fact, in the
latter case, the estimated eigenvalues and eigenvectors are inevitably
affected by the temporal dependence that underlies two consecutive
observations of interest rates changes. For this reason, this
application highlights the benefits of handling samples comprised of
independent and identically distributed observations, compared with time
series samples which might not be stationary in nature.

A potential future analysis of the bootstrap to this particular
application should incorporate the estimated autocorrelation and
cross-correlation of the interest rates changes, in order to obtain more
representative “artificial samples” . Consequently, the estimated
eigenvectors and eigenvalues on each sample will be more accurate.

To conclude the premise, we might reasonably state that the results
deriving from the non-parametric bootstrap performed in this particular
case should be considered significant either in case of stationary
interest rates or *asymptotically*, meaning that these results would be
attained only if we had an infinite amount of past observations such
that any temporal dependency is zeroed.

After this needed specification, we briefly introduce the bootstrap.

The bootstrap\[9\] is a computer intensive resampling technique based on
the simple but powerful idea of *repeated sampling* from a collection of
available observations, with the objective of evaluating the uncertainty
surrounding a parameter of interest. In our case, the original data set
is made up of 3520 daily observation across eleven interest rates. Then,
the procedure is conducted as follows. We ask R to compose 10.000
distinct cross-sectional samples made of 587 observations each, by
drawing randomly\[10\] from the actual data set. Subsequently, the
eigenvectors and eigenvalues are computed on each of the ten thousand
samples. As a result, we are capable of obtaining an estimate of their
associated variability, under the *i.i.d.* assumption.

The results of such process are illustrated in the following figures and
tables.

<img src="report-thesis_files/figure-gfm/eigenvalues-densities-1.png" style="display: block; margin: auto;" />
<img src="report-thesis_files/figure-gfm/eigenvalues-boxplots-1.png" style="display: block; margin: auto;" />

The two figures above illustrate the empirical distributions and box
plots of the first three eigenvalues resulting from 10.000 bootstrap
samples. As we expected, the first eigenvalue consistently attains
higher values compared to the second and third one, meaning that, under
stationary conditions of interest rates changes, the first component
contributes significantly to explain most of the covariation in the
system, confirming, with the specifications mentioned above, the figures
obtained earlier. Furthermore, we notice a similar degree of variability
associated to
![\\lambda\_1](https://latex.codecogs.com/png.latex?%5Clambda_1
"\\lambda_1") and
![\\lambda\_2](https://latex.codecogs.com/png.latex?%5Clambda_2
"\\lambda_2") which is significantly higher than the one attained by
![\\lambda\_3](https://latex.codecogs.com/png.latex?%5Clambda_3
"\\lambda_3").

Thus, we can conclude, with reasonable confidence, that the third
component contributes to a much lesser extent in explaining the overall
covariation of the interest rates changes. This is true, because the
estimated density is highly concentrated around the measures of central
tendencies, whilst the other two exhibit more spread.

After having considered the eigenvalues singularly, we evaluate the
variability of the cumulative variance explained by the first three
components. The results shown in the figure and in the table below
demonstrate that the first three components are capable of explaining an
amount of variance equal to 0.916 on average.

|  Min. | 1st Quantile | Median |  Mean | 3rd Quantile |  Max. |
| ----: | -----------: | -----: | ----: | -----------: | ----: |
| 0,876 |        0,909 |  0,916 | 0,916 |        0,923 | 0,952 |

Table: Descriptive statistics resulting from 10.000 bootstrap samples of
the cumulative variance explained by the first three components

<img src="report-thesis_files/figure-gfm/density-var-explained-first-three-evals-1.png" style="display: block; margin: auto;" />

For inferential purposes, it is useful to consider the 95% estimated
confidence interval for the four statistics obtained from 10.000
bootstrap samples. Therefore, if we were to repeat the estimation
process of the statistics above one hundred times, we are confident that
they would be included in the provided intervals 95 times out of 100.

Finally, the last plot represents the “probabilistic” counterpart of the
typical structure embodied by the three eigenvectors typical structure
which identifies a *parallel shift*, *tilt* and *curvature* of the yield
curve, respectively. The figure should be read from the top to the
bottom. Each quadrant illustrates the values taken by the corresponding
loading resulting from 10.000 bootstrap samples. Even if the typical
pattern of the first three eigenvectors can be recognized, we should
notice that the symmetry existing along the three eigenvectors tends to
hide it partially. This kind of symmetry is due to the centering
transformation performed on the interest rates. In fact, if we were to
obtaining the eigenvectors from the spectral decomposition of the
correlation matrix rather than singualar value decomposition, the
typical pattern would be even more manifested. Nevertheless, the first
eigenvector depicted in red, assumes with higher probability similar
values across the eleven maturities, confirming the idea that it
captures a parallel shift of the yield curve.

In the same way we can observe that the second eigenvector most of the
times represent a tilt of the yield curve, meaning that yields on
shorter maturities witness a positive change, whereas the yields on
longer ones change negatively. We might notice also that, with lower
probability, the second eigenvector captures an inverted behaviour in
which the longer yields maturities have positive changes.

Finally, the third eigenvector almost always shows the typical curvature
component.
<img src="report-thesis_files/figure-gfm/bootstrap-loadings-1.png" style="display: block; margin: auto;" />

## Bond Portfolio

This section demonstrates how to build a principal component factor
model for interest rate sensitive portfolios.

The portfolio chosen for this application comprises 26 U.S. Government
bonds expiring in a range of time that goes from a few months to 30
years (see Table below for all the details). The market prices of the
securities refers to March 10, 2020, as provided by [Business
Insider](https://markets.businessinsider.com/bonds/finder). Bonds data
were retrieved by means of a Python script which have semi-automated the
process of data collection. If run on a terminal, the program keeps
asking for the bond data which can be copied and pasted from the site.
Afterwards, the figures are stored in a dictionary which is then
converted into a *.csv* file. Subsequently, that file has been read
using R in order to perform the analysis.

For the sake of simplicity, we assume to have acquired from the market
one unit of each security at the market price listed on March 10, 2020
which we set to be our reference time point *t*.

| NAME                   | MRKT\_PRICE | CP\_RATE | POSTED\_YTM | ISN          | ISSUE\_PRICE | ISSUE\_DATE | FACE\_VAL | MATURITY\_DATE | COUPON\_PYMT\_DATE | NUMB\_PAYMENTS | STARTCOUPON\_DATE | FINALCOUPON\_DATE | PURCHASE\_DATE | DAYS\_TO\_COUPON | YTM2\_POSTED | SEMI\_COUPON\_AMOUNT | N\_CFs | SEMI\_ANN\_YTM | ANNUAL\_YTM | DELTA\_YTM | GROUP |
| :--------------------- | ----------: | -------: | ----------: | :----------- | -----------: | :---------- | --------: | :------------- | :----------------- | -------------: | :---------------- | :---------------- | :------------- | :--------------- | -----------: | -------------------: | -----: | -------------: | ----------: | ---------: | :---- |
| US TREASURY 2047       |      137,04 |    2,750 |        2,22 | US912810RZ30 |        98,97 | 2017-11-15  |       100 | 2047-11-15     | 2020-05-15         |              2 | 2018-05-15        | 2047-11-14        | 2020-03-10     | 66 days          |        0,794 |                1,375 |     56 |          0,012 |       0,024 |      0,175 | II    |
| US TREASURY 2020       |       99,30 |    1,500 |        1,76 | US9128282Q23 |        99,94 | 2017-08-15  |       100 | 2020-08-15     | 2020-08-15         |              2 | 2018-02-15        | 2020-08-14        | 2020-03-10     | 158 days         |        0,661 |                0,750 |      1 |          0,029 |       0,059 |      4,166 | VI    |
| US TREASURY 2043       |      131,44 |    2,875 |        2,20 | US912810RB61 |        97,93 | 2013-05-15  |       100 | 2043-05-15     | 2020-05-15         |              2 | 2013-11-15        | 2043-05-14        | 2020-03-10     | 66 days          |        0,789 |                1,438 |     47 |          0,013 |       0,026 |      0,448 | II    |
| US TREASURY 2036       |      155,39 |    4,500 |        1,91 | US912810FT08 |        99,51 | 2006-02-15  |       100 | 2036-02-15     | 2020-08-15         |              2 | 2006-08-15        | 2036-02-14        | 2020-03-10     | 158 days         |        0,706 |                2,250 |     32 |          0,008 |       0,016 |    \-0,296 | VI    |
| US TREASURY 2026 15.02 |      127,43 |    6,000 |        1,51 | US912810EW46 |        98,37 | 1996-02-15  |       100 | 2026-02-15     | 2020-08-15         |              2 | 1996-08-15        | 2026-02-14        | 2020-03-10     | 158 days         |        0,584 |                3,000 |     12 |          0,012 |       0,025 |      0,989 | VI    |
| US TREASURY 2028       |      105,56 |    0,521 |        0,15 | US9128283R96 |        99,54 | 2018-01-15  |       100 | 2028-01-15     | 2020-07-15         |              2 | 2018-07-15        | 2028-01-14        | 2020-03-10     | 127 days         |        0,072 |                0,261 |     16 |        \-0,002 |     \-0,003 |    \-0,487 | IV    |

Table: A sample of the Bond Portfolio as of March 10, 2020

The portfolio cash flow along the term structure is illustrated in the
figure below wherein the gray lines represent the eleven *vertices* of
the risk factors for which the historical observations are available. We
recall that the value of portfolio is sensitive to movements in those
key rates. However, since most of the projected cash flows at time *t*
occur between those vertices, we employ Svensson’s model\[11\] in order
to obtain an *interpolated* value for the spot rates at the maturities
dictated by each of the portfolio cash flows. This is possible since,
each cash flow can be regarded as a zero-coupon bond whose corresponding
yield is the one approximated by means of Svensson’s model.

<img src="report-thesis_files/figure-gfm/cash-flows-figure-1.png" style="display: block; margin: auto;" />

Figure below shows the interpolation performed by Svensson’s model\[12\]
on the yield curve observed on March 10, 2020 at the requested cash flow
maturities measured in years. In other words, since the portfolio cash
flow is a vector having dimension (1, 170), Svensson’s model adds 159
points to the yield curve at the constant desired time points indicated
by the cash flow maturities.

Afterwards, the process shown for the daily yield curve of March 10 is
repeated for each of the past curves, starting from February 09, 2006 to
March 10, 2020, in order to reconstruct the cross-sectional term
structure. The final result is shown in bottom figure. In practical
terms, applying Svensson’s model to the historical yield curves grows
the dataset of interest rates from an original dimension of
![(3521\\times 11)](https://latex.codecogs.com/png.latex?%283521%5Ctimes%2011%29
"(3521\\times 11)") to a size of
![(3521\\times 170)](https://latex.codecogs.com/png.latex?%283521%5Ctimes%20170%29
"(3521\\times 170)"). In other words, we use Svensson’s model to
reconstruct the past term structure as if we were capable of observing
each day the yield curve at the maturities dictated by the portfolio
cash flows.

<img src="report-thesis_files/figure-gfm/svensson-yield-curve-figure-1.png" style="display: block; margin: auto;" />
<img src="report-thesis_files/figure-gfm/svensson-term-structure-1.png" style="display: block; margin: auto;" />

The functional form for the estimated spot yield curve by Svensson’s
model, implemented in the *YieldCurve* R package can be seen in its
documentation.

The approximated yield curve on March 10, 2020 is then employed to
compute the vector ![p^{T}=(PV01^{(1)},...,
PV01^{(170)})](https://latex.codecogs.com/png.latex?p%5E%7BT%7D%3D%28PV01%5E%7B%281%29%7D%2C...%2C%20PV01%5E%7B%28170%29%7D%29
"p^{T}=(PV01^{(1)},..., PV01^{(170)})") at time *t* of present value of
basis point move associated to each of the portfolio cash flow. As a
result we obtain an *exact* representation of what would happen at time
*t* to the present value of each
![C^{(i)}](https://latex.codecogs.com/png.latex?C%5E%7B%28i%29%7D
"C^{(i)}") if the yield curve was to move downward by 0.01% at each of
the *n* maturities.

| Maturity |     C($) | R10032020 |     PV($) | R10032020-0.01% |     PV($) |    PV01 |
| -------: | -------: | --------: | --------: | --------------: | --------: | ------: |
|  0,05833 |   1,3125 |   0,00570 |   1,31206 |         0,00560 |   1,31207 | 0,00001 |
|  0,18333 |  17,4375 |   0,00501 |  17,42152 |         0,00491 |  17,42184 | 0,00032 |
|  0,22778 |   1,0625 |   0,00482 |   1,06134 |         0,00472 |   1,06136 | 0,00002 |
|  0,35278 |   0,6175 |   0,00443 |   0,61654 |         0,00433 |   0,61656 | 0,00002 |
|  0,39722 |   0,9375 |   0,00433 |   0,93589 |         0,00423 |   0,93593 | 0,00004 |
|  0,43889 | 123,0000 |   0,00425 | 122,77138 |         0,00415 | 122,77675 | 0,00537 |
|  0,48333 |   1,3125 |   0,00418 |   1,30986 |         0,00408 |   1,30992 | 0,00006 |
|  0,56667 |   1,3125 |   0,00409 |   1,30947 |         0,00399 |   1,30954 | 0,00007 |
|  0,69444 |  17,4375 |   0,00403 |  17,38887 |         0,00393 |  17,39007 | 0,00120 |
|  0,73611 |   1,0625 |   0,00403 |   1,05936 |         0,00393 |   1,05944 | 0,00008 |
|  0,86389 |   0,6175 |   0,00406 |   0,61534 |         0,00396 |   0,61540 | 0,00005 |
|  0,90833 |   0,9375 |   0,00408 |   0,93404 |         0,00398 |   0,93413 | 0,00008 |
|  0,95000 |  22,2500 |   0,00410 |  22,16365 |         0,00400 |  22,16575 | 0,00210 |
|  0,98611 |   1,3125 |   0,00413 |   1,30718 |         0,00403 |   1,30731 | 0,00013 |
|  1,07222 |   1,3125 |   0,00419 |   1,30663 |         0,00409 |   1,30677 | 0,00014 |
|  1,19722 |  17,4375 |   0,00431 |  17,34804 |         0,00421 |  17,35011 | 0,00207 |
|  1,24167 |   1,0625 |   0,00435 |   1,05679 |         0,00425 |   1,05692 | 0,00013 |
|  1,36667 | 100,6175 |   0,00448 | 100,00500 |         0,00438 | 100,01861 | 0,01361 |
|  1,41111 |   0,9375 |   0,00452 |   0,93155 |         0,00442 |   0,93168 | 0,00013 |
|  1,45278 |  22,2500 |   0,00457 |  22,10312 |         0,00447 |  22,10632 | 0,00320 |

Table: Risk factors sensitivites

Prior to proceed the analysis, we introduce the *profit and loss* of the
portfolio which is defined as the change in value between two
consecutive days.   
![\\Delta{P\_{t}}=P\_{t} -
P\_{t-1}](https://latex.codecogs.com/png.latex?%5CDelta%7BP_%7Bt%7D%7D%3DP_%7Bt%7D%20-%20P_%7Bt-1%7D
"\\Delta{P_{t}}=P_{t} - P_{t-1}")  
Then, “the \[profit and loss\] on the portfolio is approximated as a
weighted sum of the changes in the interest rate risk factors with
weights given by the present values of a basis point at the maturity
corresponding to the interest rate” (Alexander 2008a). Therefore, we
shall represent the profit and loss at time *t* in the following way:   
![P\_{t} - P\_{t-1} = -
\\sum\_{i=1}^{170}{C^{(i)}\\delta01\_{t}}(R\_{t}^{(i)}-R\_{t-1}^{(i)})](https://latex.codecogs.com/png.latex?P_%7Bt%7D%20-%20P_%7Bt-1%7D%20%3D%20-%20%5Csum_%7Bi%3D1%7D%5E%7B170%7D%7BC%5E%7B%28i%29%7D%5Cdelta01_%7Bt%7D%7D%28R_%7Bt%7D%5E%7B%28i%29%7D-R_%7Bt-1%7D%5E%7B%28i%29%7D%29
"P_{t} - P_{t-1} = - \\sum_{i=1}^{170}{C^{(i)}\\delta01_{t}}(R_{t}^{(i)}-R_{t-1}^{(i)})")  
where the vector ![{\\Delta{R\_{t}}^{T}}=(\\Delta{R\_{t}^{(1)}},
\\ldots,
\\Delta{R\_{t}^{(170)}})](https://latex.codecogs.com/png.latex?%7B%5CDelta%7BR_%7Bt%7D%7D%5E%7BT%7D%7D%3D%28%5CDelta%7BR_%7Bt%7D%5E%7B%281%29%7D%7D%2C%20%5Cldots%2C%20%5CDelta%7BR_%7Bt%7D%5E%7B%28170%29%7D%7D%29
"{\\Delta{R_{t}}^{T}}=(\\Delta{R_{t}^{(1)}}, \\ldots, \\Delta{R_{t}^{(170)}})")
contains the *changes* between two consecutive yield curves estimated
using Svensson. The minus sign is due to the convention of representing
the losses as positive quantities. “The \[…\] vector
\[![p^{T}](https://latex.codecogs.com/png.latex?p%5E%7BT%7D "p^{T}")\]
is held fixed at its current value so that we are measuring the interest
rate risk of the *current* portfolio” (Alexander 2008a).

At this point, if principal component were not known, we would be
constrained to model the joint behaviour of the entire set of 159
interest rates changes. This would imply to consider 159 variances and
25281 covariances, for a total number of 14535 parameters which is
clearly unfeasible. Conversely, using the principal component
approximation derived previously, we are able to represent the interest
rates changes using just three principal components, and at the same
time being still capable of explaining most of the covariation in the
system.

Then, the *principal component representation* using the first three
components is obtained from the singular value decomposition of the
Svensson interest rates changes.   
![\\Delta{R\_{t}^{(i)}} \\approx w\_{1}^{(i)}{z\_{t}^{(1)}} +
w\_{2}^{(i)}{z\_{t}^{(2)}} +
w\_{3}^{(i)}{z\_{t}^{(3)}}](https://latex.codecogs.com/png.latex?%5CDelta%7BR_%7Bt%7D%5E%7B%28i%29%7D%7D%20%5Capprox%20w_%7B1%7D%5E%7B%28i%29%7D%7Bz_%7Bt%7D%5E%7B%281%29%7D%7D%20%2B%20w_%7B2%7D%5E%7B%28i%29%7D%7Bz_%7Bt%7D%5E%7B%282%29%7D%7D%20%2B%20w_%7B3%7D%5E%7B%28i%29%7D%7Bz_%7Bt%7D%5E%7B%283%29%7D%7D
"\\Delta{R_{t}^{(i)}} \\approx w_{1}^{(i)}{z_{t}^{(1)}} + w_{2}^{(i)}{z_{t}^{(2)}} + w_{3}^{(i)}{z_{t}^{(3)}}")  
Using matrix notation the PCA approximation for the entire set of
observations would be as follows:   
![\\underset{(3520\\times 170)}{{\\Delta}{R}}\\approx\\underset{(3520
\\times 3)}{Z}\\underset{(3\\times 170)}{W}](https://latex.codecogs.com/png.latex?%5Cunderset%7B%283520%5Ctimes%20170%29%7D%7B%7B%5CDelta%7D%7BR%7D%7D%5Capprox%5Cunderset%7B%283520%20%5Ctimes%203%29%7D%7BZ%7D%5Cunderset%7B%283%5Ctimes%20170%29%7D%7BW%7D
"\\underset{(3520\\times 170)}{{\\Delta}{R}}\\approx\\underset{(3520 \\times 3)}{Z}\\underset{(3\\times 170)}{W}")  
Then, the *j-th* principal component *risk factor sensitivity* of the
portfolio is computed as:   
![k\_{j} = - \\sum\_{i=1}^{170}
{PV01^{(i)}}{w\_{j}^{(i)}}](https://latex.codecogs.com/png.latex?k_%7Bj%7D%20%3D%20-%20%5Csum_%7Bi%3D1%7D%5E%7B170%7D%20%7BPV01%5E%7B%28i%29%7D%7D%7Bw_%7Bj%7D%5E%7B%28i%29%7D%7D
"k_{j} = - \\sum_{i=1}^{170} {PV01^{(i)}}{w_{j}^{(i)}}")  
or, equivalently,   
![k\_{j}= -
{p^{T}}w\_{j}](https://latex.codecogs.com/png.latex?k_%7Bj%7D%3D%20-%20%7Bp%5E%7BT%7D%7Dw_%7Bj%7D
"k_{j}= - {p^{T}}w_{j}")  
which measures the change in the portfolio value when the principal
component risk factor changes, keeping all the other principal component
risk factors constant (Alexander 2008b, 1:32).

The entire set of risk factor sensitivities are computed as follows:   
![\\underset{(1\\times 3)}{k^{T}}=p^{T}{W}](https://latex.codecogs.com/png.latex?%5Cunderset%7B%281%5Ctimes%203%29%7D%7Bk%5E%7BT%7D%7D%3Dp%5E%7BT%7D%7BW%7D
"\\underset{(1\\times 3)}{k^{T}}=p^{T}{W}")  
Eventually, *the principal component factor model representation* of the
portfolio profit and loss is:   
![\\Delta{P\_{t}}\\approx
k^{T}{z\_{t}}](https://latex.codecogs.com/png.latex?%5CDelta%7BP_%7Bt%7D%7D%5Capprox%20k%5E%7BT%7D%7Bz_%7Bt%7D%7D
"\\Delta{P_{t}}\\approx k^{T}{z_{t}}")  
where ![{z}\_{t}^{T} = (z\_{t}^{(1)},
z\_{t}^{(2)},z\_{t}^{(3)})](https://latex.codecogs.com/png.latex?%7Bz%7D_%7Bt%7D%5E%7BT%7D%20%3D%20%28z_%7Bt%7D%5E%7B%281%29%7D%2C%20z_%7Bt%7D%5E%7B%282%29%7D%2Cz_%7Bt%7D%5E%7B%283%29%7D%29
"{z}_{t}^{T} = (z_{t}^{(1)}, z_{t}^{(2)},z_{t}^{(3)})") and
![k^{T}=(k\_1,k\_2,k\_3)](https://latex.codecogs.com/png.latex?k%5E%7BT%7D%3D%28k_1%2Ck_2%2Ck_3%29
"k^{T}=(k_1,k_2,k_3)") denote the principal component risk factor at
time *t*, and their (constant) risk factor sensitivities.

Using PCA we reduced the number of risk factor from 159 to
![k=3](https://latex.codecogs.com/png.latex?k%3D3 "k=3").

The estimated risk factor sensitivities are shown in the following
table:

|      w1 |    w2 |    w3 |
| ------: | ----: | ----: |
| \-0,355 | 0,205 | 0,052 |

Table: Risk factors sensitivites

Therefore, the resulting PCA factor model for this example is:   
![-p\\\&l\_{t}\\approx 0.3547\\times z\_{t}^{(1)}-0.2048\\times
z\_{t}^{(2)}-0.0521\\times
z\_{t}^{(3)}](https://latex.codecogs.com/png.latex?-p%5C%26l_%7Bt%7D%5Capprox%200.3547%5Ctimes%20z_%7Bt%7D%5E%7B%281%29%7D-0.2048%5Ctimes%20z_%7Bt%7D%5E%7B%282%29%7D-0.0521%5Ctimes%20z_%7Bt%7D%5E%7B%283%29%7D
"-p\\&l_{t}\\approx 0.3547\\times z_{t}^{(1)}-0.2048\\times z_{t}^{(2)}-0.0521\\times z_{t}^{(3)}")  
which can be used to immunize the portfolio against the most frequent
movements of the yield curve.

<img src="report-thesis_files/figure-gfm/p&l-graph-1.png" style="display: block; margin: auto;" />

## Conclusions

This text examined the US spot term structure observed between 2006 and
2020. In particular, the empirical analysis was conducted on a
multivariate time series comprised of eleven key interest rates provided
by the US Treasury which showed high correlation, in particular between
yields of closer maturity. As a consequence, we were able to approximate
the dynamics of the interest rates changes with an accuracy of almost
91% using only the first three principal components computed on those
interest rates. In addition, we confirmed the existence of the
traditional structure for the eigenvectors of the sample correlation
matrix identifying a parallel *shift* a *tilt* and *curvature* as main
yield curve movements occurring between consecutive time instances.
However, as regard to stability, this was not entirely true when the
analysis was performed on shorter successive time windows, even if some
regularities were still evident, in particular as far as the first
eigenvector was concerned. This might suggests that on shorter time
windows principal component analysis is affected by short term
volatilities of the interest rates.

Further, we showed the approximated empirical distribution associated to
the first three eigenvalues and eigenvectors resulting from 10.000
bootstrap samples. We were immediately able to conclude that even if the
first two eigenvalues attained consistently higher values compared to
the third one, there were a greater amount of uncertainty underlying
their distributions. Moreover, we were able to show that the first three
components were able to explain about 92% of the variation in the
interest rates changes within a 95% confidence interval, under the
*i.i.d.* assumption which is not necessarily true in the context of
financial time series.

In the last section we included in the analysis a portfolio made of
bonds expiring at different maturities. We then succeeded to approximate
the historical yield curves using Svensson’s model to interpolate the
historical daily yield curves, in order to obtain a series of risk
factors associated to portfolio cash flows projected from March 10,
2020. Finally we employed the principal component representation to
approximate the profit and loss of the portfolio which is of interest in
risk management applications.

To conclude, it seems useful to provide some ideas for potential future
analysis:

  - Coding a bootstrap function that might take into consideration the
    autocorrelation of the interest rates.
  - Investigating more deeply the differences arising from applying
    spectral decomposition on the correlation matrix instead of the
    covariance.
  - Applying the technique of cash flow mapping and investigating the
    possible differences with the results obtained with Svensson’s
    model.
  - Making further steps in the understanding of the uses of the profit
    and loss obtained by means of PCA.

## Bibliography

<div id="refs" class="references">

<div id="ref-carol2">

Alexander, Carol. 2008a. *Practical Financial Econometrics*. Vol. 2.
Market Risk Analysis. Wiley.

</div>

<div id="ref-carol1">

———. 2008b. *Quantitative Methods in Finance*. Vol. 1. Market Risk
Analysis. Wiley.

</div>

<div id="ref-yield">

Choudhry, Moorad. 2004. *Analysing and Interpreting the Yield Curve*.
Sixth. Wiley.

</div>

<div id="ref-connor">

Connor, Gregory. 1995. “The Three Types of Factor Models: A Comparison
of Their Explanatory Power.” *Financial Analysts Journal* 51 (3): 41–43.

</div>

<div id="ref-embrechts">

Embrechts, Paul, Alexander J. McNeil, and Rudiger Frey. 2015.
*Quantitative Risk Management: Concepts, Techniques and Tools*.
Princeton University Press.

</div>

<div id="ref-fabozzihandbook">

Fabozzi, Frank J. 2012. *The Handbook of Fixed Income Securities*.
Eighth. McGraw-Hill.

</div>

<div id="ref-fabozzi">

———. 2013. *Bond Markets, Analysis, and Strategies*. Eight. Pearson.

</div>

<div id="ref-hull">

Hull, John C. 2018. *Risk Management and Financial Institutions*. Fifth.
Wiley.

</div>

<div id="ref-ams">

Johnson, Richard A., and Dean W. Wichern. 2014. *Applied Multivariate
Statistical Analysis*. Sixth. Pearson.

</div>

<div id="ref-jolly">

Jolliffe, I. T. 2002. *Principal Component Analysis*. Second. Springer.

</div>

<div id="ref-yieldcurvestrategies">

Jones, Frank J. 1991. “Yield Curve Strategies.” *The Journal of Fixed
Income* 1 (2): 43–48.

</div>

<div id="ref-rao">

Rao, Radhakrishna. 1964. “The Use and Interpretation of Principal
Component Analysis in Applied Research.” *The Indian Journal of
Statistics* 26 (4): 329–58.

</div>

<div id="ref-moody">

Redfern, David, and Douglas McLean. 2004. “Principal Component Analysis
for Yield Curve Modelling.” Moody’s Analytics Research.

</div>

<div id="ref-sironiresti">

Sironi, Andrea, and Andrea Resti. 2007. *Risk Management and
Shareholders’ Value in Banking*. Wiley.

</div>

<div id="ref-tsay">

Tsay, Ruey S. 2010. *Analysis of Financial Time Series*. Third. Wiley.

</div>

</div>

1.  Standard measures and models to manage interest rate risk such as
    the *Duration Gap* are discussed by (Sironi and Resti 2007).
    Particular attention is devoted to the asset-liabilities mismatch.

2.  The random process generating the behaviour charachterizing the
    yield curve is also modeled by means of stochastic differential
    equations. Those models are reviewed in (Choudhry 2004). In
    contrast, the approach adopted here is statistical in nature. The
    difference between the two frameworks can be well appreciated in
    Redfern and McLean (2004)

3.  Yet, this approach is subject to another kind of risk known as
    *historical bias*. It might not be the case that the past will
    repeat itself in the future. In trying to reduce the historical bias
    risk, one might plug into the model subjective information adopting
    a *bayesian* approach. *Model risk* should also be taken into
    account.

4.  *Multi-factor models* are extensively used in finance. Their aim is
    to explain the covariance structure of portfolios’ asset returns, as
    a function of *p* underlying factors. In the literature, one can
    identify three types of factors models depending, on the estimation
    strategy adopted. A comparison between their explanatory power on
    U.S. equities is provided by (Connor 1995). On the other hand, the
    estimation process of the three models, - *macroeconomic*,
    *fundamental* and *statistical* - is discussed in (Embrechts,
    McNeil, and Frey 2015) and, more extensively, in (Tsay 2010).
    Generally speaking, both the macroeconomic and fundamental require
    the researcher to provide the observations on a factor, for example
    in the form of a economic indicator or a financial index, whereas
    the latter allow the analyst to estimate directly from the data the
    main factors that might drive the overall risk of the portfolio.

5.  Fabozzi (2012) reports the academic studies which applied PCA on the
    term structure of different countries. In particular: Robert
    Litterman and Jose Scheinkman, “Common Factors Affecting Bond
    Returns”, Journal of Fixed Income (September 1991), pp. 54-61;
    Alfred Buhler and Heinz Zimmermann, “A Statistical Analysis of the
    Term Structure of Interest Rates in Switzerland and Germany”,
    Journal of Fixed Income (December 1996), pp. 55-67; Joel R. Barber
    and Mark L. Copper, “Immunization Using Principal Component
    Analysis”, Journal of Portfolio Management (Fall 1996),
    pp. 99-105; Rita L. D’Ecclesia and Stavros Zenios, “Risk Factor
    Analysis and Portfolio Immunization in the Italian Bond Market”,
    Journal of Fixed Income (September 1994), pp. 51-58; Bennett W.
    Golub and Leo M. Tilman, “Measuring Yield Curve Risk Using Principal
    Components Analysis, Value at Risk, and Key Rate Durations”, Journal
    of Portfolio Management (Summer 1997), pp. 72-84; Lionel Martellini
    and Philippe Priaulet, Fixed-Income Securities: Dynamic Methods for
    Interest Rate Risk Pricing and Hedging (Chichester, England: Wiley,
    2000); Sandrine Lardic, Philippe Priaulet, and Stephane Priaulet,
    “PCA of Yield Curve Dynamics: Questions of Methodologies”, Journal
    of Bond Trading and Management (April 2003), pp. 327-349. Fabozzi
    (2012) reports also the number of factors (PCs) employed and the
    percentage of variance explained. In all of the cases it is greater
    than 80%.

6.  It is worth noting that, within this case, the eigenvectors, and
    consequently the principal components, have straightforward
    interpretation since the system of interest rates is *ordered* and
    comprised of variables of the same kind. Differently, the principal
    components would be the result of linear combinations of variables
    having independent meaning. Therefore, they would be still valid
    from a mathematical point of view, yet their interpretation would be
    more obscure. These issues are discussed by Rao (1964) in a very
    technical paper, and more accessibly by (Jolliffe 2002) in Chapter
    11.

7.  The general formulation would be
    ![R=XW+\\Psi](https://latex.codecogs.com/png.latex?R%3DXW%2B%5CPsi
    "R=XW+\\Psi")

8.  The *bootstrap* produces approximately samples with i.i.d.
    observations, notwithstanding technicalities attributable to the
    process with which a computer generates random numbers.

9.  This technique was conceived by renowned statistician Bradley Efron
    as an improvement over the *jackniffe*.

10. The random draws should be performed with replacement.

11. An alternative strategy relies on *Cash Flow Mapping* which consists
    in assigning each cash flow falling on non key maturities to the
    vertices by keeping the main financial characteristics of the
    portfolio invariant. Since this process would have been quite
    laborious to carry out in R, we have decide to embrace Svensson’s
    model which is capable of approximating the spot rates at the
    requested non standard maturities.

12. In R the *YieldCurve* package implements Nelson-Siegel, Diebold-Li
    and Svensson models.
