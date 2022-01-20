% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a normal distribution for the continuous prior and an
% exponential distribution for the likelihood. The parameter describing the
% likelihood is lambda which is the reciprocal of the scale factor by convention
% Also returns the values of the continuous prior and continuous posterior
% and likelihood at a specificed set of points to facilitate plotting.


function [obs_increments, prior_pts, post_pts, like_pts, err] =  obs_increment_normal_exp(ensemble, like_lambda, y)

% Does a case where the prior is normal and the observation is exponential
% Same idea could be used for the other way around where the prior was exponential and the
% likelihood was normal. Recall that the exponential is just a gamma with k = 1.

% Set error return to default successful
err = 0;
ens_size = size(ensemble, 2);

% Get the sample statistics from this ensemble
prior_mean = mean(ensemble);
prior_var = var(ensemble);
prior_std = sqrt(prior_var);

% Will use gamma to do the stats for likelihood which requires the scale
like_scale = 1/like_lambda;

% Formula for posterior normal distribution
post_mean = prior_mean - like_lambda * prior_var;
post_var = prior_var;
post_std = sqrt(post_var);

% BUT this posterior is only for y>0, so need to normalize
trimmed_cdf = normcdf(0, post_mean, post_std);
post_wt = 1 / (1 - trimmed_cdf);

% Get the quantiles for the prior ensemble
q = normcdf(ensemble, prior_mean, prior_std);

% Get the posterior values by inverting this but remember that it's trucated
for i = 1: ens_size
   q_needed = q(i) / post_wt + trimmed_cdf; 
   updated_ensemble(i) = norminv(q_needed, post_mean, post_std);
end

% Compute the increments
obs_increments = updated_ensemble - ensemble;

% Compute the point values of the prior mean for plotting
prior_pts = normpdf(y, prior_mean, prior_std);
like_pts = gampdf(y, 1, like_scale);

for i = 1:size(y, 2)
   if(y(i) < 0) 
      post_pts(i) = 0;
   else
      post_pts(i) = normpdf(y(i), post_mean, post_std) * post_wt;
   end
end

