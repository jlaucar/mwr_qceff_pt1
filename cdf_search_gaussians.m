% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes the inverse of the CDF (quantile function) corresponding to a PDF that is 
% a weighted sum of num_g normals with specified covariance and means.
% Does a search to find a good estimate of the location of a given cdf value.
% Takes two initial guesses that may bound the search interval, the value of the cdf to be found
% and eventually info on things like stopping criteria. Returns the value of x with approximately
% the correct CDF. Calls a function to compute the CDF. For now, assume that the parameters for that
% function are set some other place.


function[x_val, approx_cdf] = cdf_search_gaussians(cdf_val, x_lo, x_hi, num_g, mu, sigma, wt, tol)  

% Begin by computing the value of the cdf at the suggested end points
cdf_lo = get_cdf_gaussians(x_lo, num_g, mu, sigma, wt);
cdf_hi = get_cdf_gaussians(x_hi, num_g, mu, sigma, wt);

% Now make sure we actually have a lower bound
% Eventually may want to have an ability to assert this
while(cdf_lo > cdf_val)
   x_lo = x_hi - 2*(x_hi - x_lo);
   cdf_lo = get_cdf_gaussians(x_lo, num_g, mu, sigma, wt); 
end

% Make sure we have an upper bound
while(cdf_hi < cdf_val)
   x_hi = x_lo + 2*(x_hi - x_lo);
   cdf_hi = get_cdf_gaussians(x_hi, num_g, mu, sigma, wt); 
end

% Now do a binary search
% Note that since these things are smooth we can do better than a binary search, can use gradient
cdf_mid = 100.0;
while(abs(cdf_mid - cdf_val) > tol)
   x_mid = x_lo + (x_hi - x_lo) / 2;
   cdf_mid = get_cdf_gaussians(x_mid, num_g, mu, sigma, wt);
   if(cdf_mid > cdf_val) 
      x_hi = x_mid;
      cdf_hi = cdf_mid;
   else
      x_lo = x_mid;
      cdf_lo = cdf_mid;
   end
end

x_val = x_mid;
approx_cdf = cdf_mid;

end
