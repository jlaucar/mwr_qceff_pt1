% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a gamma distribution for the continuous prior and a
% gamma distribution for the likelihood. The parameters describing the
% likelihood are the shape and the scale.
% Also returns the values of the continuous prior and continuous posterior
% at a specificed set of points to facilitate plotting.


function [obs_increments, prior_pts, post_pts, err] =  obs_increment_gamma_gamma(ensemble, like_shape, like_scale, y)

% Gamma prior, gamma likelihood update. This is analogous to Craig Bishop's GIG filter

% Set error return to default successful
err = 0;

% Get the sample statistics from this ensemble, again following Bishop
prior_mean = mean(ensemble);
prior_type1_var = var(ensemble) / prior_mean^2;
prior_shape = 1 / prior_type1_var;
prior_scale = prior_mean * prior_type1_var;

% Compute the posterior shape and scale
post_shape = prior_shape + like_shape - 1;
post_scale = (prior_scale * like_scale) / (prior_scale + like_scale);

% Get the quantiles for the prior ensemble
q = gamcdf(ensemble, prior_shape, prior_scale);

% Get the posterior values by inverting this 
updated_ensemble = gaminv(q, post_shape, post_scale);

% Compute the increments
obs_increments = updated_ensemble - ensemble;

% Compute the point values of the prior mean for plotting
prior_pts = gampdf(y, prior_shape, prior_scale);
post_pts = gampdf(y, post_shape, post_scale);
