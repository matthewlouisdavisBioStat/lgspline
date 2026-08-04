# Introduction

This R package implements Lagrangian multiplier smoothing splines, which reformulate smoothing splines through constrained optimization. This approach provides direct access to predictor-response relationships through interpretable coefficients, unlike other formulations that require post-fitting algebraic manipulation.

Additionally, this package allows for fitting survival models, GLMs, MMRMs, and many other models sunder arbitrary linear equality and inequality constraints upon coefficients.

I have ensured that asymptotic confidence interval coverage remains at the nominal 95% level for all models contained in this package, including for survival models and models with marginal correlation structures.

I will not be continuing to pursue publication for this idea, due to personal life, work, and difficulty in finding an appropriate journal. For those interested, an unpublished manuscript with reproducible code appears in the "Article" folder, which justifies and explains the proposed method. 

# Installation 

install.packages('lgspline')
devtools::install_github("matthewlouisdavisBioStat/lgspline")

# Citation 

If you use this package in your research, please cite:

Davis, M. (2025). Lagrangian Multiplier Smoothing Splines. https://github.com/matthewlouisdavisBioStat/lgspline/

# Contact

For questions or feedback, please open an issue on GitHub.
