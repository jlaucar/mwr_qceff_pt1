% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a lognormal distribution for the continuous prior and a
% lognormal distribution for the likelihood. The parameters describing the
% likelihood are the observation and the observation error variance in log space.
% Also returns the values of the continuous prior and continuous posterior
% at a specificed set of points to facilitate plotting.


function [obs_increments, prior_pts, post_pts, err] =  obs_increment_eakf(ensemble, observation, obs_error_var, y)

% Obs increment for a lognormal/lognormal 
% The observations passed in are for the distribution in log space

% Set error return to default successful
err = 0;

ens_size = size(ensemble, 2);

% Convert the ensemble to log space
log_ens = log(ensemble);

% Compute prior ensemble mean and variance
prior_mean = mean(log_ens);
prior_var = var(log_ens);
prior_sd = sqrt(prior_var);

[prior_mean, prior_var, prior_sd];
[observation, obs_error_var, sqrt(obs_error_var)];

post_var = 1 / (1 / prior_var + 1 / obs_error_var);
post_sd = sqrt(post_var);

% Compute posterior mean
post_mean = post_var * (prior_mean / prior_var + observation / obs_error_var);

% Compute the point values of the prior mean for plotting up higher
for i = 1:size(y, 2);
   if(y(i) <= 0)
      prior_pts(i) = 0;
   else
      prior_pts(i) = normpdf(log(y(i)), prior_mean, prior_sd);
   end
end

% Shift the prior ensemble to have the posterior mean
log_updated_ens = log_ens - prior_mean + post_mean;

% Contract the ensemble to have the posterior_variance
var_ratio = post_var / prior_var;
log_updated_ens = sqrt(var_ratio) * (log_updated_ens - post_mean) + post_mean;
updated_ensemble = exp(log_updated_ens);

% Compute the increments
obs_increments = updated_ensemble - ensemble;

% Compute the point values of the post mean for plotting
for i = 1:size(y, 2);
   if(y(i) <= 0)
      post_pts(i) = 0;
   else
      post_pts(i) = normpdf(log(y(i)), post_mean, post_sd);
   end
end
