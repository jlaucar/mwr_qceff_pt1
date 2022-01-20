% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.

% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a particle filter for the 'continuous' prior and a
% normal distribution for the likelihood. The parameters describing the
% likelihood are the observation and the observation error variance.
% In essence, this is a straightforward way to do a deterministic analogy of the
% stochastic particle filter.

function [obs_increments, err] =  obs_increment_particle(ensemble, observation, obs_error_var, y)

ens_size = size(ensemble, 2);
% Set error return to default successful
err = 0;

% Sort the ensemble members
[s_ens, s_indx] = sort(ensemble);

% Compute the values of the likelihood at the ensemble members
weight = normpdf(s_ens, observation, sqrt(obs_error_var));

% Normalize the weights
weight_sum = sum(weight);
weight = weight / weight_sum;

% The CDF has discrete jumps
post_cdf(1) = weight(1);
for i = 2:ens_size
   post_cdf(i) = post_cdf(i-1) + weight(i);
end

% For now, assume the prior has ranks that are i / (n+1);
% Obviously this can be done much more quickly
for i = 1:ens_size
   quant(i) = i / (ens_size + 1);
   for j = 1: ens_size
      if(quant(i) < post_cdf(j))
         s_updated_ens(i) = s_ens(j);
         break
      end
   end
end

% Now convert back into the unsorted order
for i = 1:ens_size
   updated_ens(s_indx(i)) = s_updated_ens(i);
end

% Compute the increments
obs_increments = updated_ens - ensemble;
