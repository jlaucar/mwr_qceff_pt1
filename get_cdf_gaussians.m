% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Evaluates the CDF at x for a PDF that is a weighted sum of normals with different means
% and standard deviations.

function[cdf] = get_cdf_gaussians(x, num, mu, sigma, wt)

cdf = 0;
for i = 1:num
   cdf = cdf + wt(i) * normcdf(x, mu(i), sigma(i));
end
