% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a normal distribution for the continuous prior and a
% normal distribution for the likelihood. The parameters describing the
% likelihood are the observation and the observation error variance.
% Also returns the values of the continuous prior and continuous posterior
% at a specificed set of points to facilitate plotting.
%
% This is a modified version of the script that is part of the DART_LAB tutorial
% in DART. Note that the EAKF here is computed using the traditional algorithm
% rather than the equivalent quantile conserving method.



function [obs_increments, prior_pts, post_pts, err] =  obs_increment_eakf(ensemble, observation, obs_error_var, y)

% Set error return to default successful
err = 0;

% Compute prior ensemble mean and variance
prior_mean = mean(ensemble);
prior_var = var(ensemble);

% Compute the point values of the prior mean for plotting up higher
prior_pts = normpdf(y, prior_mean, sqrt(prior_var));

% If both prior and observation error variance are non-positive return error
if(prior_var <= 0 && obs_error_var <= 0)
   err = 1;
   return;
end

% Compute the posterior mean and variance
% If prior variance is 0, posterior mean is prior_mean and variance is 0
if(prior_var == 0)
   post_mean = prior_mean;
   post_var = 0;
elseif(obs_error_var == 0)
% If obs_error_var is 0, posterior mean is observation and variance is 0
   post_mean = observation;
   post_var = 0;
else
% Use product of gaussians
   % Compute the posterior variance
   post_var = 1 / (1 / prior_var + 1 / obs_error_var);

   % Compute posterior mean
   post_mean = post_var * (prior_mean / prior_var + observation / obs_error_var);
end

% Shift the prior ensemble to have the posterior mean
updated_ensemble = ensemble - prior_mean + post_mean;

% Contract the ensemble to have the posterior_variance
var_ratio = post_var / prior_var;
updated_ensemble = sqrt(var_ratio) * (updated_ensemble - post_mean) + post_mean;

% Compute the increments
obs_increments = updated_ensemble - ensemble;

% Compute the point values of the prior mean for plotting
post_pts = normpdf(y, post_mean, sqrt(post_var));
