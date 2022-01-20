% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a binormal distribution for the continuous prior and a
% normal distribution for the likelihood. The parameters describing the
% likelihood are the observation and the observation error variance.
% Also returns the values of the continuous prior and continuous posterior
% at a specificed set of points to facilitate plotting.

function [obs_increments, prior_pts, post_pts, err] =  obs_increment_binormal(ensemble, observation, obs_error_var, y_pts)

% For this demonstration, there is no clustering done. The clusters are hard code specified here and
% are assumed to be in order.

% Set error return to default successful
err = 0;

ens_size = size(ensemble, 2);
cluster_size(1) = 4;
cluster_size(2) = ens_size - cluster_size(1); 

% Compute prior ensemble mean and variance for the standard gaussian part of the likelihood
prior_mean(1) = mean(ensemble(1:cluster_size(1)));
prior_var(1) = var(ensemble(1:cluster_size(1)));
prior_mean(2) = mean(ensemble(cluster_size+1:ens_size));
prior_var(2) = var(ensemble(cluster_size+1:ens_size));
prior_sd(1:2) = sqrt(prior_var(1:2));

% Prior weights are relative to number of ensemble members 
prior_weight(1:2) = cluster_size(1:2) ./ ens_size;

% Figure out the quantile in the prior CDF of the original ensemble members
pos = prior_weight(1) * normcdf(ensemble, prior_mean(1), prior_sd(1));
pos = pos + prior_weight(2) * normcdf(ensemble, prior_mean(2), prior_sd(2));

% Get the points for the prior distribution
prior_pts = prior_weight(1) * normpdf(y_pts, prior_mean(1), prior_sd(1));
prior_pts = prior_pts + prior_weight(2) * normpdf(y_pts, prior_mean(2), prior_sd(2));

% Use product of gaussians
% Compute the posterior variance
for i = 1:2
   post_var(i) = 1 / (1 / prior_var(i) + 1 / obs_error_var);
end

% Mean and weight can be different for each ensemble posterior kernel
for i = 1:2
   % Compute posterior mean
   post_mean(i) = post_var(i) * (prior_mean(i) / prior_var(i) + observation / obs_error_var);

   % Also need a weight in this case since we have two gaussians in likelihood
   weight(i) =  exp(-0.5 * (prior_mean(i)^2 / prior_var(i) + ... 
      observation^2 / obs_error_var - post_mean(i)^2 / post_var(i)));

   % Multiply by the prior relative weights
   weight(i) = weight(i) * prior_weight(i);
end

% Normalize the weights
weight_sum = sum(weight);
weight = weight / weight_sum;

% Find the value of these positions in the posterior continuous distribution
post_sd(1:2) = sqrt(post_var(1:2));

% Find the x-value for each of the original cdfs
for i = 1:ens_size
   min_ens = min(ensemble);
   max_ens = max(ensemble);
   tol = (max_ens - min_ens) * 0.0001;
   [post_ens(1, i), approx_cdf] = cdf_search_gaussians(pos(i), min_ens, max_ens, 2, post_mean, post_sd, weight, tol);
end

% Compute the increments
obs_increments = post_ens- ensemble;

% Get the value of the posterior for the requested points
post_pts = weight(1) * normpdf(y_pts, post_mean(1), post_sd(1));
post_pts = post_pts + weight(2) * normpdf(y_pts, post_mean(2), post_sd(2));
