======================================================


       GRAM - Gaussian Rating Aggregation Model

 
 Jacob Whitehill, Cecilia Aguerrebere and Benjamin Hylak

 Do Learners Know What’s Good for Them? Crowdsourcing 
 Subjective Ratings of OERs to Predict Learning Gains
 
 Int. Conf. on Educational Data Mining, 2019

 Version 1.0 - April, 2019
 by Cecilia Aguerrebere <caguerrebere@gmail.com>

======================================================


 Requirements
-------------------

- R environment for statistical computing <https://www.r-project.org/>
- data.table package <https://cran.r-project.org/web/packages/data.table/index.html>


  Running the code
-------------------

The script run_experiments.R runs all the experiments presented in
sections 6.1 and 6.2 of the article.

The function EM in the lib.R file implements the GRAM model, with 
parameters:

  - d_train <data.frame>      Data frame containing the training set.
  - Niter <numeric>           Number of iterations to run the EM algorithm.
  - V <numeric>               Standard deviation of the prior distribution 
			      of the educational resources quality.
  - U <numeric>               Mean of the prior distribution of the educational 
			      resources quality.
  - EPS <numeric>	      Small additive factor to avoid numerical problems.
  - useMu <numeric>	      Indicates how to estimate the mean in the EM 
			      algorithm:
				-- 1, as defined by the EM algorithm
				-- else, mu is set to 0 
  - sigmaType <string>	      Indicates how to estimate the standard deviation 
			      in the EM algorithm: 
				-- "sigma", as defined by the EM algorithm 
				-- "pretest", using pretest scores only
				-- else, standard deviation is set to 1
  - verbose <numeric>         Indicates whether to print (1) or not (0) the 
			      Spearman and Pearson correlation at each EM iteration.
  - verbose_pval <numeric>    Indicates whether to print (1) or not (0) the 
			      p-values of the correlations. 


-------------------------------------------------------------------------------
 Copyright and License
-----------------------

Copyright (c)

This code is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

------------------------------------------------------------------------------
 Feedback
----------

Do not hesitate to contact me if you have any comments, especially about errors,
bugs, or strange results. <caguerrebere@gmail.com>
