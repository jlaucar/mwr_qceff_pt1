% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a normal distribution for the continuous prior and a
% weighted sum of two normal distributions for the likelihood. The parameters describing the
% likelihood are the observation and the observation error variance of a standard
% normal. The variance of a second normal is 10 times the variance of the standard normal.
% Also returns the values of the continuous prior and continuous posterior
% and likelihood at a specificed set of points to facilitate plotting.

function [obs_increments, prior_pts, like_pts, post_pts, err] =  ...
   obs_increment_double_like(ensemble, observation, obs_error_var, y_pts)


% Specify the probability of the obs error coming from the broader normal
like_weight(2) = 0.1;
like_weight(1) = 1 - like_weight(2);

like_var(1) = obs_error_var;
like_var(2) = obs_error_var * 10;
like_sd = sqrt(like_var);

% Get the value of the likelihood for the requested points
like_pts = like_weight(1) * normpdf(y_pts, observation, like_sd(1));
like_pts = like_pts + like_weight(2) * normpdf(y_pts, observation, like_sd(2));

% Set error return to default successful
err = 0;

ens_size = size(ensemble, 2);

% Compute prior ensemble mean and variance
prior_mean = mean(ensemble);
prior_var = var(ensemble);
prior_sd = sqrt(prior_var);

% Figure out the position in the prior CDF of the original ensemble members
pos = normcdf(ensemble, prior_mean, prior_sd);

% Compute the point values of the prior mean for plotting
prior_pts = normpdf(y_pts, prior_mean, prior_sd);

% Compute the posterior mean and variance
% Use product of gaussians
% Compute the posterior variance
for i = 1:2
   post_var(i) = 1 / (1 / prior_var + 1 / like_var(i));
   post_sd(i) = sqrt(post_var(i));

   % Compute posterior mean
   post_mean(i) = post_var(i) * (prior_mean / prior_var + observation / like_var(i));

   % Also need a weight in this case since we have two gaussians in likelihood
   weight(i) =  exp(-0.5 * (prior_mean^2 / prior_var + ... 
      observation^2 / like_var(i) - post_mean(i)^2 / post_var(i)));
end

% Need to add the specified weight fraction
weight(1:2) = like_weight .* weight(1:2);

% Normalize the weights
weight_sum = sum(weight);
weight = weight / weight_sum;

% Find the x-value for each of the original cdfs
for i = 1:ens_size
   ens_max = max(ensemble);
   ens_min = min(ensemble);
   % Make the tolerance a fraction of the ensemble range
   tol = (ens_max - ens_min) * 0.001;
   post_ens(1, i) = cdf_search_gaussians(pos(i), ens_min, ens_max, 2, post_mean, post_sd, weight, tol);
end

% Compute the increments
obs_increments = post_ens - ensemble;

% Get the value of the posterior for the requested points
post_pts = weight(1) * normpdf(y_pts, post_mean(1), post_sd(1));
post_pts = post_pts + weight(2) * normpdf(y_pts, post_mean(2), post_sd(2));
