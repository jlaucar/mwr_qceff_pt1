% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a bounded normal distribution for the continuous prior and a
% normal distribution for the likelihood. The parameters describing the
% likelihood are the observation and the observation error variance.
% Also returns the values of the continuous prior and continuous posterior
% at a specificed set of points to facilitate plotting.


function [obs_increments, prior_pts, post_pts, err] =  obs_increment_pos_eakf(ensemble, observation, obs_error_var, y)

% Obs increment for an eakf for a positive definite quantity. 
% Fits a normal but then makes the lower bound 0. 
% Does a standard continuous update

% Set error return to default successful
err = 0;

ens_size = size(ensemble, 2);

% Compute prior ensemble mean and variance
prior_mean = mean(ensemble);
prior_var = var(ensemble);
prior_sd = sqrt(prior_var);

post_var = 1 / (1 / prior_var + 1 / obs_error_var);
post_sd = sqrt(post_var);

% Compute posterior mean
post_mean = post_var * (prior_mean / prior_var + observation / obs_error_var);

% Find the quantile for each prior ensemble member
base_q = normcdf(0, prior_mean, prior_sd);
% Normalization factor to make the bounded prior a pdf
prior_wt = 1 / (1 - base_q);

% quantile for prior in the bounded prior
for i = 1:ens_size
   raw_q = normcdf(ensemble(i), prior_mean, prior_sd);
   q(i) = prior_wt * (raw_q - base_q);
end

% Compute the point values of the prior mean for plotting
prior_pts = prior_wt * normpdf(y, prior_mean, sqrt(prior_var));
for i = 1:size(y, 2)
   if(y(i) < 0) prior_pts(i) = 0; end
end

% Need to get the appropriate weight for the bounded posterior to be a pdf
post_base_q = normcdf(0, post_mean, post_sd);
% Normalization factor to make the bounded prior a pdf
post_wt = 1 / (1 - post_base_q);

% Need to get posterior x
% q = w * [F(x) - F(b))
% x = F-1[q/w + F(b)]
for i = 1:ens_size
   updated_ensemble(i) = norminv(q(i) / post_wt + post_base_q, post_mean, post_sd);
end

% Compute the increments
obs_increments = updated_ensemble - ensemble;

% Compute the point values of the prior mean for plotting
post_pts = post_wt * normpdf(y, post_mean, sqrt(post_var));
for i = 1:size(y, 2)
   if(y(i) < 0) post_pts(i) = 0; end
end
