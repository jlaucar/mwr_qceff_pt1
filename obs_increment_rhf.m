% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% Computes Quantile Conserving Ensemble Filter increments for an input prior 
% ensemble using a rank histogram distribution for the peicewise continuous prior and a
% piecewise continuous approximation of a normal distribution for the likelihood. 
% The parameters describing the
% likelihood are the observation and the observation error variance.
% Also returns the values of the piecewise continuous prior and posterior
% and the piecewise continuous likelihood 
% at a specificed set of points to facilitate plotting.

% This is a modified version of the script that is part of the DART_LAB tutorial
% in DART.


function [obs_increments, prior_pts, post_pts, like_pts, err] =  obs_increment_rhf(ensemble, observation, obs_error_var, y_pts)
% obs_increment_rhf Computes increments for a rank histogram filter
% For now this implements the default options

% Set error return to default successful
err = 0;

% Get the ensemble size
ens_size  = size(ensemble, 2);
prior_sd  = std(ensemble);
prior_var = prior_sd^2;

% Allocate space for some of the density calculation bits
like_dense = zeros(1,ens_size); % likelihood density
mass       = zeros(1,ens_size);
height     = zeros(1,ens_size);

% Sort the ensemble members and keep the indices
[x, e_ind] = sort(ensemble);

% Compute the likelihood of each member given the observation
like = exp(-1 * (x - observation).^2 / (2 * obs_error_var));

% Compute the mean likelihood density in each interior bin
for i = 2:ens_size
   like_dense(i) = (like(i - 1) + like(i)) / 2;
end

% For unit normal, find distance from mean to where cdf is 1/(n+1)
dist_for_unit_sd =  -1 * weighted_norm_inv(1, 0, 1, 1/(ens_size + 1));

% Variance of tails is just sample prior variance
% Mean is adjusted so that 1/(ens_size + 1) is outside
left_mean = x(1) + dist_for_unit_sd * prior_sd;
left_var  = prior_var;
left_sd   = prior_sd;

% Same for the right tail
right_mean = x(ens_size) - dist_for_unit_sd * prior_sd;
right_var  = prior_var;
right_sd   = prior_sd;

% Eventually want to support options, for now
gaussian_likelihood_tails = false;

if(gaussian_likelihood_tails)
   % Need to fill this in
else
   % Block to do flat tails for likelihood follows
   % This removes assumptions about likelihood and cuts cost
   new_sd_left      = left_sd;
   new_mean_left    = left_mean;
   prod_weight_left = like(1);
   mass(1)          = like(1) / (ens_size + 1);

   % Same for right tail
   new_sd_right       = right_sd;
   new_mean_right     = right_mean;
   prod_weight_right  = like(ens_size);
   mass(ens_size + 1) = like(ens_size) / (ens_size + 1);
   % End block for flat tail likelihood
end

% The mass in each interior box is the height times the width
% The height of the likelihood is like_dense
% For the prior, mass is 1 / ((n+1) width), and mass = height x width so
% The height of the prior is 1 / ((n+1) width); Multiplying by width leaves 1/(n+1)

% In prior, have 1/(n+1) mass in each bin, multiply by mean likelihood
% density to get approximate mass in updated bin

for i = 2:ens_size
   mass(i) = like_dense(i) / (ens_size + 1);
   % Height of prior in this bin is mass/width; Only needed for trapezoidal
   % If two ensemble members are the same, set height to -1 as flag
   if(x(i) == x(i-1))
      height(i) = -1;
   else
      height(i) = 1 / ((ens_size + 1) * (x(i) - x(i-1)));
   end
end

% Now normalize the mass in the different bins to get a pdf
mass_sum = sum(mass);
nmass = mass / mass_sum;

% Get the weight for the final normalized tail gaussians
% This is the same as left_amp=(ens_size + 1)*nmass(1)
left_amp = prod_weight_left / mass_sum;
% This is the same as right_amp=(ens_size + 1)*nmass(ens_size+1)
right_amp = prod_weight_right / mass_sum;

% Find cumulative mass at each box boundary and middle boundary
cumul_mass = zeros(1,ens_size+1);
for i = 1:ens_size + 1
   cumul_mass(i+1) = cumul_mass(i) + nmass(i);
end

% Begin internal box search at bottom of lowest box
lowest_box = 1;

new_ens        = zeros(1,ens_size);
obs_increments = zeros(1,ens_size);
% Find each new ensemble member's location
for i = 1:ens_size
   % Each update ensemble member has 1/(n+1) mass before it
   umass = i / (ens_size + 1);

   % If it is in the inner or outer range have to use normal
   if(umass < cumul_mass(2))
      % It is in the left tail
      % Get position of x in weighted gaussian where the cdf has value umass
      new_ens(i) = weighted_norm_inv(left_amp, new_mean_left, ...
         new_sd_left, umass);
   elseif (umass > cumul_mass(ens_size + 1))
      % It's in the right tail
      % Get the position of x in weighted gaussian where the cdf has value umass
      new_ens(i) = weighted_norm_inv(right_amp, new_mean_right, ...
         new_sd_right, 1 - umass);
      % Coming in from the right, use symmetry after pretending it's on left
      new_ens(i) = new_mean_right + (new_mean_right - new_ens(i));
   else
      % In one of the inner boxes
      for j = lowest_box:ens_size - 1
         % Find the box that this mass is in
         if(umass >= cumul_mass(j+1) && umass <= cumul_mass(j+2))

            % Only rectangular is implemented for now
            rectangular_quadrature = true;

            if(rectangular_quadrature)
               % Rectangular quadrature block first
               % Linearly interpolate in mass
               new_ens(i) = x(j) + ((umass - cumul_mass(j+1)) / ...
                  (cumul_mass(j+2) - cumul_mass(j+1))) * (x(j+1) - x(j));
            else
               % Trapezoidal interpolation block goes here
            end

            % Don't need to search lower boxes again
            lowest_box = j;
            break
         end
      end
   end
end

% Convert to increments for unsorted
for i = 1:ens_size
   obs_increments(e_ind(i)) = new_ens(i) - x(i);
end

prior_pts(1:size(y_pts, 2)) = 0;
% Compute the values of the prior and posterior at each of the requested y points for plotting
for i = 1: size(y_pts, 2)
   % This is expensive, could be optimized
   if(y_pts(i) < x(1)) 
      y_bin = 1;
   elseif(y_pts(i) > x(ens_size))
      y_bin = ens_size + 1;
   else
      for j = 1:ens_size - 1
         if(y_pts(i) >= x(j) && y_pts(i) <= x(j+1))
            y_bin = j + 1;
         end
      end
   end

   if(y_bin > 1 && y_bin < ens_size + 1)
      prior_pts(i) = height(y_bin);
      post_pts(i) = nmass(y_bin) / (x(y_bin) - x(y_bin - 1)); 
      like_pts(i) = like_dense(y_bin);
   elseif(y_bin == 1)
      prior_pts(i) = normpdf(y_pts(i), left_mean, left_sd); 
      post_pts(i) = left_amp * normpdf(y_pts(i), new_mean_left, new_sd_left); 
      like_pts(i) = like_dense(2) / 2; 
   elseif(y_bin == ens_size + 1)
      prior_pts(i) = normpdf(y_pts(i), right_mean, right_sd); 
      post_pts(i) = right_amp * normpdf(y_pts(i), new_mean_right, new_sd_right); 
      like_pts(i) = like_dense(ens_size) / 2;
   end

end


%-----------------------------------------------


function [x] = weighted_norm_inv(alpha, mean, sd, p)

% Find the value of x for which the cdf of a N(mean, sd)
% multiplied times alpha has value p.

% Can search in a standard normal, then multiply by sd at end and add mean
% Divide p by alpha to get the right place for weighted normal
np = p / alpha;

% Find spot in standard normal
x = norm_inv(np);

% Add in the mean and normalize by sd
x = mean + x * sd;

