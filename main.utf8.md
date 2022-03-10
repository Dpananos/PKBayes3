---
title: "Paper 3"
author:
  - "A. Demetri Pananos"
  - "Daniel J. Lizotte"
  - "Simon J. Bonner"
output: pdf_document
date: "March 10, 2022"
bibliography: references.bib  
---







# Introduction

One goal of personalized medicine is optimized dosing of drugs for individuals [-@morse2015personalized].  When considering optimal doses, a thorough understanding pharmacokinetic (what the body does to the drug) and/or pharmacodynamic (what the drug does to the body) effects are crucial.  To this end, models describing the mediation of pharmacokinetic/pharmacodynamic effects via clinical, genetic, and lifestyle factors have an important role in deciding which patients should get what dose and are sometimes published by research teams collaborating with drug manufacturers using data from clinical trials.

Independent investigators can find themselves in a situation in which data collection from a particular population of interest is achievable. If the data come from practice (e.g. a personalized medical clinic), there may be questions about how new variables not previously studied in clinical trials effect the pharmacokientics/pharamcodynamics of a particular drug.  Running large studies in order to examine the effects of these new variables, or discover effects of other variables, may be unrealistic due to a variety of constraints.  Consequently, investigators must think about how best to model the pharmacokinetics, for use in decision making *and* exploration, using the data available to them.

The oral anti-coagulant apixaban provides an illustrative example.  Pharmacokinetic models have been previously published [-@cirincione2018population, -@ueshima2018population] in collaboration with the drug's manufacturer using data from clinical trials.  These studies identified age, sex, body weight, renal function,  patient race, and CYP3A4 inhibitors as modulators of apixaban pharmacokinetics [-@cirincione2018population], though according to authors the effects of some of these variables were not large enough to require clinical dose adjustment. However, even after adjusting for the aforementioned factors, concentrations of apixaban in real life applications have been observed to be larger than what was reported in clinical trials [-@sukumar2019apixaban], raising questions as to the optimal dosing of apixaban for patients outside these trials. Additionally, recent research has indicated appropriateness of dose adjustment criteria are unclear [-@vu2021critical], citing there is no reduction in safety in patients above 75 years of age, below 60kg of body weight, and egfr lower than 50 mL/min. The uncertainty regarding dosing criteria and additional variability in concentrations in day to day use suggest that, while previously published models may be internally valid, these models may not be representative of all populations in which apixaban is to be applied.  That is to say, the models may lack a degree of external validity, thus supporting the idea that pharmacokinetic models may need to be tailored for specific populations of interest. When viewed through a Bayesian lens, the previous modeling work can act as an informative prior on various pharmacokinetic/pharmacodynamic measures.  Creating new models for populations of interest is then more of a "fine tuning" than an all together new approach.  Pharmacokinetic models for use in a specific population may then have two goals:  to adjust dosing criteria to be specific for a population, and/or to explore how additional variables (for example, concomitant medications) not included in the previous studies effect particular parts of the pharmacokinetics of apixaban.  

This study seeks to demonstrate how investigators can fit similar models to their pharmacokinetic data with the aim of accomplishing the goals of accurate modeling of pharmacokientics and exploration of effects of new variables.  We use apixaban as a specific example, but our methodology can be generalized to other drugs where pharmacokientics are of interest.  Importantly, we only focus on pharmacokientics since blood plasma levels correlated closely with the pharmacodynamic effect of apixaban [-@byon2019apixaban]. Our approach leverages a Bayesian methodology to building pharmacokinetic models so that we may incorporate prior information from previous studies. Additionally, we describe how investigators can use *all* relevant data available to them to fit these models and make inferences, even if that data come from controlled studies.  Finally, we show how sparsity inducing priors can be applied to new variables in order to explore how those variables may effect apixaban pharmacokinetics, encouraging negligible effect sizes but allowing for large effects to be detected. We present a small simulation study to demonstrate how smallest meaningful effects can be detected through these priors as a function of sample size. Finally, we use an open source Bayesian language to develop our models, making our code freely available.  Previous models are constructed in a proprietary software tool set, which can come at a high cost. Creation of these models in a free tool removes a barrier to research, making these methods more widely available.



# Background

## Apixaban

Our study uses clinical and experimental plasma concentrations from patients who were prescribed apixaban. Apixaban is a direct acting oral anti-coagulant often prescribed for prevention of stroke and systemic embolism in patients with atrial fibrillation (AF) [-@BMSmonograph; -@byon2019apixaban].  Studies as recent as 2019 have reported excess variability in observed apixaban plasma concentrations in patients with AF [-@sukumar2019apixaban]. Since apixaban plasma concentrations correlate closely with anti-coagulation[-@upreti2013effect,-@frost2013safety,-@frost2013apixaban], excess variability in these concentrations may mean increased risk of bleeding. These findings have raised questions towards the optimal dosing of apixaban in older adults with AF encountered outside of clinical trials.

Additional research into determining factors which explain this excess variability beyond known clinical factors [-@gulilat2020drug] has consequently begun.

## Variable Importance & Simulation Study

Existing studies often use variable selection methods (e.g. variants of stepwise selection, including fitting all submodels [*CITE*]) when faced with the determining which variables effect the phenomenon under study. Many studies have noted that these techniques result in bias away from the null, exaggerated precision, inaccurate or uninterpretable p values due to inability to properly incorporate uncertainty in the selection process, and can fail to select the "true" model with high confidence even when modeling assumptions are consistent with the true data generating process [*CITE MANY PAPERS*].  Hence, even in the best case scenario where the selection procedure identifies the correct variables, the resulting estimates may not be reliable.  Results from simulation studies on selection methods make a convincing argument to avoid selection methods all together.

Selection methods are intended to answer the question "which variables are important in modeling the outcome", and although simulation studies have demonstrated deficiencies with variable selection, they often do not provide an alternative answer.  From a Bayesian perspective, selection to include a variable in or out of a model defines a sort of prior on the parameter value; there is a strong preference for a null effect estimate unless the data provide sufficient evidence for free estimation of that effect.  Efforts to operationalize this prior structure in terms of Bayesian inference have lead to a wide variety of sparsity inducing priors, which include spike and slab priors [*CITE*], and horseshoe and Finish horseshoe priors [*CITE*]. These approaches admit that while unlikely that the effects of unimportant variables are exactly 0, they may be small enough to be negligible.  These priors place the majority of their probability mass near 0, encouraging small effects to be estimated as something negligibly small, but allow for large effects to be identified (with perhaps a small amount of bias depending on the prior hyperparameters and prior structure).

Bias towards a null effect  can be acceptable when the goal is exploration and prediction.  The bias can act as a regularization for predictions hence combating overfitting, and can hedge estimates of novel effects when they are exaggerated due to high variance.  From our perspective, discovery is not about getting the estimate right the first time, its about making progress and identifying directions for further investigation.  Bias in the estimates from sparsity inducing priors puts investigators on the right path with the hopes that additional studies will provide more precise and potentially unbiased estimates.

Biasing the effect downward is likely favorable in scenarios where the covariate of interest is only recorded from one dataset available to investigators, much like the scenarios we describe here.  Densely sampled data may come from highly controlled studies, with very explicit inclusion/exclusion criteria.  Unless the covariate is of primary interest in those studies, it may be the case that subjects are highly homogeneous in many respects (e.g. All being healthy young adults, with no concomitant medications).  In these cases, it may be unlikely that the covariate of interest was recorded, if subjects with that covariate were to be included at all.  It is more probable that novel covariates are collected from observational data (e.g. from a personalized medicine clinic which sees patients irregularly).  The observational nature of the data allows a wider collection of patients to be observed, hence making it more likely that the novel covariate is observed in subjects.  However, observational data are limited in so far as there is likely one measurement per subject, likely due to the high patient burden of dense sampling.  This means that, unless the pharmacokientics are very well understood and all relevant covariates are measured, effects of novel covariates may be confounded due to omitted variable bias.  If the bias is away from the null effect or if the estimate is highly variable due to the covariate not being sufficiently variable within the sample, regularization via informative priors can combat this. If the bias is towards the null effect, regularization makes this worse (assuming the effect would be detected at all in the absence of regularization and that the bias due to regularization is sufficiently large to hide this effect).  

To this end, we present a simulation study in which we use a sparsity inducing prior to estimate the effect of a concomitant medication on apixaban pharmacokinetics.  In particular, the medication is assumed to inhibit a particular gene important in the elimination of apixaban, making the bioavailability larger.  We place a double exponential (or Laplace) prior on the effect of the concomitant medication, as well as a prior on the parameter for the Laplace distribution.  This is similar to putting a LASSO penalty on the effect as well as a prior on the LASSO penalty strength [*CITE?*].  Although our simulation only has a single variable of interest, many variables can be used with this prior structure.


<!-- $$ p(\beta_F) = \int \mbox{normal}(\beta_k \vert 0, \sigma) \,  \mbox{exponential}\left(\sigma^2 \Bigg\vert \dfrac{1}{2\tau^2} \right) \,d\sigma = \mbox{Laplace}(\beta_k \vert \tau) $$ -->

For our simulation, we generate data from the posterior of a previously fit model [*CITE SECOND PAPER*]. We simulate 10 datasets from a pre-specified number of densely sampled patients (we examine 5, 10, 20, 30, 40, and 50 densely sampled patients) haven taken their first dose of the drug with a pre-specified and fixed effect of a concomitant medication on the bio-availability of the drug.  We assume that investigators can sparsely sample patients more easily, and so we simulate 10x more sparsely sampled patients who have already achieved steady state. We do this so as to more closely resemble real life scenarios in which patients come into a clinic for a plasma measurement having already been on the drug for sometime. We examine effects of 0, 0.125, 0.25, 0.5, 1.0, and 1.5 on the logit scale (we use the logit scale since bioavailability is constrained to be between 0 and 1). 

## "You May Think Method X might work, and you'd be wrong"

In this paper, we propose a pooling of both sparsely and densely sampled data in a single model. Pooling information is not a new approach, and reasonable arguments could be made to use simpler models.  After all, if the sparsely sampled data models a continuous outcome as a function of covariates, why would investigators use a complex model when something simple like linear regression (or linear regression on log concentrations) may be sufficient.  While simpler approaches and criticisms of using unnecessarily complex models are valid, both linear modeling and mixed effects models for pooling suffer from important drawbacks in the case when attempting to combine sparsely sampled and densely sampled data from different studies.  We examine those drawbacks below.

Linear regression can be, and has been [*CITE*], used to model apixaban concentrations as a function of time and other covariates using sparsely sampled data.  When certain criteria are met, there is good reason to do so.  The concentration profile, $y(t)$, from a first order absorption with linear elimination pharmacokinetic model looks like

$$ y(t) = \frac{F \cdot D}{C l} \frac{k_{e} \cdot k_{a}}{k_{e}-k_{a}}\left(e^{-k_{a}t}-e^{-k_{e}t}\right) $$
The elimination phase occurs when $t$ is sufficiently large, resulting in $y(t)$ being approximately exponential and $\log(y(t))$ being linear in time with slope $-k_e$.  Assuming measurement error is additive on the log scale facilitates use of linear regression.

This approach is common in pharmacokinetics when estimating the elimination rate but suffers from three important drawbacks generally. First, the elimination rate is not allowed to vary as a function of known factors which effect elimination rate, such as kidney function.  This can be ameliorated by specifying an interaction between time and those covariates known to effect elimination rate (though this has not been done in all papers [*CITE*]).  Second, an exponential approximation is only appropriate when time is sufficiently large.  Clearly, the exponential approximation breaks down near $t=t_{max} = \log(k_a/k_e)/(k_a-k_e)$ and is completely inappropriate in the absorption phase when $\partial y(t) / \partial t >0$.  This effects estimates of max concentration in an appreciable way, resulting in an upward bias of $C_{\max}$ (the bias is proportional to $-k_a t_{\max}$).  Additionally, because $t_{\max}$ is not modeled per individual, estimates of $C_{\max}$ must rely on a point estimate of $t_{\max}$.  This results in uncertainty estimates of $C_{\max}$ which may be too narrow for a given individual.  Finally, the effects of covariates on other aspects of the pharmacokinetcs are undetermined. Assuming a linear model is used to model concentrations on the log scale, we find

$$ \log(y(t)) \approx \log(D) + \log(F) - \log(Cl) + \log(k_e) + \log(k_a) - \log(k_e-k_a) - k_et = \log(D) + \beta_0 + \beta_1t \>. $$
Here, $\beta_0 =  \log(F) - \log(Cl) + \log(k_e) + \log(k_a) - \log(k_e-k_a)$.  If covariates are included in the model, then although changes in log concentration may be accurate (in so far as the magnitude of the increase or decrease in log concentration is concerned), *where that change occurs is under determined*.  Did concentration increase because bioavailability ($F$ in the log linear model) increased, or was it because the clearance rate ($Cl$ in the log linear model) decreased?  We can't say for certain from this model. In order to determine if a change in concentration was due to an increase/decrease in a pharmacokinetic parameter, each pharamacokinetic parameter must be modeled as functions of covariates. How salient these drawbacks are is up to the investigator, but if any of them are important to decision making for personalized medicine then a linear model may not be appropriate.

Mixed effects models can be used to pool information from many datasets.  Meta-analysis is perhaps the most prevalent example of this approach.  A typical example may be pooling data from studies conducted using similar protocols across multiple centers.  Ideally, the data are collected under similar protocols, making assumptions regarding the likelihood and exchangeability of appropriate units tenable.  In the scenario we describe, where information from at least two studies with different protocols are to be pooled, we believe a mixed effect model specifying between study variation is *not* appropriate due to subjects not being exchangeable between studies.  Recall, a sequence of random variables $\theta_1, \dots, \theta_n$ is said to be exchangeable in their joint density if $p(\theta_1, \dots, \theta_n)$ is invariant permutations of the indicies $(1, \dots, n)$ [*CITE BDA3*].  If no other information, other than observed data, is available to distinguish any of the $\theta_j$ from any others, and no ordering or grouping of parameters can be made, one must assume exchageability of the $\theta$.  In the scenarios we describe, we do have additional information which can be used to distinguish the $\theta$.  In particular, sparsely sampled data will have a larger estimated residual error than densely sampled data.  This is because the residual variance is a combination of within and between subject variation. There are then 2 residual errors to be estimated: one for the densely sampled data and one for the sparsely sampled data.  When pooling sparsely sampled and densely sampled data together, individuals within subjects are exchangeable because of the common residual variance within study. However, subjects are not exchangeable between studies because permutations of the subject indicies fail to account for which subject should be associated with which residual variance.


# Methods

By now, we have established a few points which are worth reiterating:

* In cases where specific populations of interest display sufficient differences as compared to data presented in clinical trials, the tailoring of a pharmacokinetic model to the intended population for use in personalized medicine applications may be desirable.  Apixaban is one such drug which displays such differences, and we use data from apixaban as a case study.


* Models developed for specific populations may serve two purposes: prediction and exploration of novel effects.  Previous studies rely on variable selection procedures to identify novel predictors of pharmacokinetics.  As of late, many simulation studies have demonstrated deficiencies in variable selection procedures.  We propose the use of sparsity inducing priors to regularize negligible effects towards 0 while keeping those variables in the models.

* Investigators may want to make use of all pharmacokinetic data available to them, be they sparsely sampled and of an observational nature, or densely sampled from well controlled clinical studies.  Typical methods for pooling this information (e.g. mixed effect models, or simply concatenating datasets) are not universally appropriate in cases when data come from studies with vastly different protocols, thereby violating exchangeability.

In what follows, we present a Bayesian model which addresses all three points above.  We present a model of apixaban pharmacokinetics which combines data from two studies with different protocols.  We demonstrate how sparse priors can be used to estimate the effect of a potentially novel predictor, and present a simulation study to investigate how relative sample sizes between the two studies and effect size of the predictor effect estimation.

## Bayesian Model

Our model specifies a population level effect of covariates (age, sex, weight (kg), serum creatinine $\mu \mbox{mol}$) on patient clearance, time to max concentration, and the ratio between absorption and elimination rates (a unitless parameter we refer to as $\alpha$). These effects are shared between all populations, allowing information from one dataset to partially inform model fit on the other.  We also include a population level effect of concomitant amidarone on bioavailability of apixaban.  

We fit our model using Stan [-@gelman2015stan], an open source probabilistic programming language with interfaces to Python, R, Stata, Matlab, and more.  Fitting two datasets jointly is formally equivalent to fitting one dataset first and passing the posterior as a prior for a model for the other dataset.  However, such an approach requires  of the posterior, which may result in loss of information (such as covariance between draws, unless explicitly modeled).  The recommended approach is then to fit datasets jointly, as done in [-@fallingBetancourt].



### Dense Sampling (Dataset 1) Model

Since patients are observed multiple times in these data, this offers the opportunity to estimate random effects for Clearence $Cl$, time to mac concentration $t_{\max}$, ratio between elimination and absorption rates $\alpha = k_e/k_a$.

Let $X$ be a matrix of mean centered and standardized covariates for dataset 1.  For patient $j$, we model the pharmacokinetic parameters as

$$ \log(Cl_{j}) = \mu_{Cl} + X\beta_{Cl} + z_{Cl, j} \sigma_{Cl}  $$
$$ \log(t_{\max,j}) = \mu_{t_{\max}} + X\beta_{t_{max}} + z_{t_{\max, j}}\sigma_{t_{\max}}  $$

$$ \operatorname{logit}(\alpha_{j})  =  \mu_{\alpha} + X\beta_{\alpha} + z_{\alpha, j}\sigma_{\alpha}  $$

Here, the $\mu$ are the population level means for the indicated pharmacokinetic parameters, the $\beta$ are the regression coefficients, the $z$ are standard normal random variables to account for random effects, and the $\sigma$ are the standard deviations of the population distribution for the indicated pharmacokinetic parameters. Additionally, we model the population mean for the bioavailability as $F = 1/(1 + e^{\mu_F}$, with a prior on $\mu_F$.  Both the $\mu$ and the $\beta$ are shared between datasets.

For dataset 1, we also model a delay between ingestion and absorption of apixaban.  The delay is modeled as 

$$ \delta_j = 0.5 \times b  $$

Where $b$ is a beta distributed random variable with parameters learned from the data.  The factor of 0.5 is used to ensure that at $t=0.5$ hours after ingestion, the predicted to be non-zero.

We use a one compartment pharmacokinetic model with first order elimination as our conditional mean

$$  C_{j}(t)= \begin{cases}\frac{F \cdot D}{C l_j} \frac{k_{e, j} \cdot k_{a, j}}{k_{e, j}-k_{a, j}}\left(e^{-k_{a, j}(t-\delta_j)}-e^{-k_{e, j}(t-\delta_j)}\right) & \delta_j \leq t \\ 0 & \text { else }\end{cases} \>. $$
Here, all parameters are estimated from data, and we have used the facts that

$$ t_{\max }=\frac{\ln \left(k_{a}\right)-\ln \left(k_{e}\right)}{k_{a}-k_{e}} $$
$$ \alpha = \dfrac{k_e}{k_a} $$

in order to solve for $k_e$ and $k_a$ for use in our PK model.  Finally, we specify a lognormal likelihood for dataset 1

$$ y_j \sim \mbox{Lognormal}(C_j(t), \sigma_1) \>. $$

For information of prior distributions, see our supplement.

### Real Life Data (Dataset 2) Model

Much of the structure from the previous model is translated to dataset 2.  However, there are a few differences.

Because patients in this dataset are not measured multiple times we do not estimate random effects or a time delay.  Hence, we model

$$ \log(Cl) = \mu_{Cl} + X\beta_{Cl} $$

$$ \log(t_{\max}) = \mu_{t_{\max}} + X\beta_{t_{\max}} $$

$$\operatorname{logit}(a) = \mu_\alpha + X\beta_{\alpha} $$

Where the $mu$ and the $\beta$ are shared between datasets.  Additionally, we model the bioavailbility as 

$$ \operatorname{logit}(F) = \mu_F   + \beta_{amio} \mbox{amiodarone}$$

as dataset 2 has information regarding concomitant amiodarone.  The prior for the effect of concomitant amiodarone is a sparsity inducing prior, meaning it encourages negligible effects (near 0) but will allow for large effects to be identified.  The effect of amiodarone is not included in the model for dataset 1 as the inclusion/exclusion criteria specified that subjects were not using CYP3A4 inhibitors.

Finally, data from dataset 2 comes from patients who have been taking apixaban twice daily and are assumed to be at steady state.  Hence, their initial plasma concentration on ingestion is not 0, but can be modeled using the pharmacokientics none the less.  Assuming the patients have been taking apixaban twice a day, 12 hours apart, for the last 5 days, the initial concentration can be shown to be

$$ c_0 = \sum_{j=1}^{10} C(12j) \>. $$

Here, $C(t)$ is the pharmacokintic profile with estimated coefficients for that patient.  See our supplement for a proof of this proposition.

Our pharmacokinetic profile is again provided by a one compartment first order elimination 

$$  C(t)= c_0 + \frac{F \cdot D}{C l} \frac{k_{e} \cdot k_{a}}{k_{e}-k_{a}}\left(e^{-k_{a}(t)}-e^{-k_{e}(t)}\right)  $$
and we assume a lognormal likelihood

$$ y_j \sim \mbox{Lognormal}(C(t), \sigma_2) \>. $$
Note the likelihood for dataset 2 has a different observational noise component ($\sigma_2$) as compared to dataset 1 ($\sigma_1$).  This is because random effects can not be estimated from dataset 2, hence the residual variance is part observational noise and part between subject variability conditional on the subject covariates.

# Results

The results from our simulation study are shown below. The precision of the estimate of effect of concomitant drug use increases as the number of densely sampled (and sparsely sampled) patients increases.  Show in read are the sample means of the 10 runs (black dots).  On average we see a small amount of bias in the estimates.  This is expected since the sparsity inducing priors have the majority of their density in a small neighborhood of 0, regularizing effects towards 0.  For purposes of discovery, these biases may be acceptable.


When using real data, our model can accurately predict both densely sampled and sparsely sampled data.  Shown in figure x is a log-log plot of predicted and actual concentrations for both datasets.  The model makes more accurate predictions for densely sampled patients (because it is able to estimate the random effect in each pharmacokinetic parameter).  The apparent increase in prediction error for the sparsely sampled can be explained by the absence of random effects for each patient.  The within and between patient variation manifests as measurement error solely, thus leading to lower predictive ability.  This perspective is supported when examining the measurement error posterior distributions for both dense/sparse patients.

With a model for the pharmacokinetics of apixaban in hand, estimates of salient pharmacokinetic phenomena can be easily obtained.  In figure y, we use our model to estimate the max concentration for the reference patient under different doses of amiodarone.  Through our model, we estimate concomitant amiodarone increases bioavailability, which in turn increases max concentration.  Shown in black is the expected max concentration conditioned on concomitant amiodarone dose, as well as 95% equal tailed posterior credible intervals.


\begin{center}\includegraphics{main_files/figure-latex/unnamed-chunk-1-1} \end{center}








\begin{center}\includegraphics{main_files/figure-latex/plot-model-predictions-1} \end{center}

```
## [1] "figures/prediction_v_actual.png"
```


\begin{center}\includegraphics{main_files/figure-latex/unnamed-chunk-2-1} \end{center}


# References
