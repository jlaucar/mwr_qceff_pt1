% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a beta distribution for the continuous prior and a
% beta distribution for the likelihood. The parameters describing the
% likelihood are the alpha and beta.
% Also returns the values of the continuous prior and continuous posterior
% and likelihood at a specificed set of points to facilitate plotting.


function [obs_increments, prior_pts, post_pts, like_pts, err] = ...
     obs_increment_beta_beta(ensemble, like_alpha, like_beta, y)

% Does a case where the prior is beta and the likelihood is beta
% Note that the beta includes any power of x, for instance a line or parabola
% Also note this works on interval [0 1], but this could be mapped to any bounded interval

% Set error return to default successful
err = 0;

% Get the sample statistics from this ensemble
prior_mean = mean(ensemble);
prior_var = var(ensemble);

% Get the alpha and beta for the prior (see wikipedia)
nu = prior_mean * (1 - prior_mean) / prior_var - 1;
prior_alpha = prior_mean * nu;
prior_beta = (1 - prior_mean) * nu;

% Compute the posterior alpha and beta
post_alpha = prior_alpha + like_alpha - 1;
post_beta = prior_beta + like_beta - 1;

% Get the quantiles for the prior ensemble
q = betacdf(ensemble, prior_alpha, prior_beta);

% Get the posterior values by inverting this but remember that it's trucated
updated_ensemble = betainv(q, post_alpha, post_beta);

% Compute the increments
obs_increments = updated_ensemble - ensemble;

% Compute the point values of the prior mean for plotting
prior_pts = betapdf(y, prior_alpha, prior_beta);
like_pts =  betapdf(y, like_alpha,  like_beta);
post_pts =  betapdf(y, post_alpha,  post_beta);
