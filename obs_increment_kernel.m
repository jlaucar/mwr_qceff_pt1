% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes observation increments for a Gaussian kernel QCEF given the prior
% ensemble, the observation, and the observation error variance. Also returns
% the value of the piecewise continuous prior and posterior distributions 
% and the piecewise continuous approximate likelihood at a specified set of points
% to facilitate plotting.


function [obs_increments, prior_pts, post_pts, err] =  obs_increment_kernel(ensemble, observation, obs_error_var, y_pts)

% Case where continuous prior is represented as sum of gaussian kernels centered around
% the prior ensemble members. Note that there is a stochastic version of this in DART that has
% been there forever.

% Set error return to default successful
err = 0;

ens_size = size(ensemble, 2);

% Compute prior ensemble mean and variance for the standard gaussian part of the likelihood
prior_var = var(ensemble);

% Narrow the prior variance since these are kernels; Need to think about this factor for practical application
prior_var = prior_var / 10;

% If both prior and observation error variance are return error
if(prior_var <= 0 && obs_error_var <= 0)
   err = 1;
   return;
end

% Compute the posterior variance; Worry about obscure error conditions for practical application
if(prior_var == 0)
   post_var(1:ens_size) = 0;
elseif(obs_error_var == 0)
% If obs_error_var is 0, posterior mean is observation and variance is 0
   post_var(1:ens_size) = 0;
else
% Use product of gaussians
   % Compute the posterior variance
   post_var(1:ens_size) = 1 / (1 / prior_var + 1 / obs_error_var);
end

% Mean and weight can be different for each ensemble posterior kernel
for i = 1:ens_size
   % Compute posterior mean
   post_mean(i) = post_var(1) * (ensemble(i) / prior_var + observation / obs_error_var);

   % Also need a weight in this case since we have multiple gaussians
   weight(i) =  exp(-0.5 * (ensemble(i)^2 / prior_var + ... 
      observation^2 / obs_error_var - post_mean(i)^2 / post_var(1)));
end

% Figure out how the position in the prior CDF of the original ensemble members
sum_pos(1, 1:ens_size) = 0;
for i = 1:ens_size
   pos = normcdf(ensemble, ensemble(i), sqrt(prior_var));
   sum_pos(1:ens_size) = sum_pos(1:ens_size) + pos;
end
sum_pos = sum_pos ./ ens_size;

% Find the value of these positions in the posterior continuous distribution
post_sd = sqrt(post_var);

% Normalize the weights
weight_sum = sum(weight);
weight = weight / weight_sum;

% Find the x-value for each of the original cdfs
for i = 1:ens_size
   min_ens = min(ensemble);
   max_ens = max(ensemble);
   tol = (max_ens - min_ens) * 0.0001;
   [post_ens(1, i), approx_cdf] = cdf_search_gaussians(sum_pos(i), min_ens, max_ens, ens_size, post_mean, post_sd, weight, tol);
end

% Compute the increments
obs_increments = post_ens- ensemble;

% Get the value of the posterior for the requested points
prior_pts = zeros(size(y_pts));
post_pts = prior_pts;
for i = 1:ens_size
   prior_pts = prior_pts + (1 / ens_size) * normpdf(y_pts, ensemble(i), sqrt(prior_var));
   post_pts = post_pts + weight(i) * normpdf(y_pts, post_mean(i), sqrt(post_var(i)));
end
