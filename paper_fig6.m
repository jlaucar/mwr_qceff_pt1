% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% This script generates Figure 6.

% The ensemble prior has 4 members drawn from a normal distribution with mean
% -2 and standard deviation 0.2 and 6 members drawn from a normal distribution
% with mean 2 and standard deviation 0.2. The likelihood is a normal with
% mean 2 and variance 1. Three continuous priors and a deterministic QCEF type
% particle filter are applied. The continuous priors are a standard EAKF (normal),
% a binormal, and a Gaussian kernel filter. 

ens_size = 10;

% Create the random draw for the prior ensemble
r_seed = 28;
rng(r_seed);
y_prior(1:4) = -2 + 0.2*randn(1, 4);
y_prior(5:10) = 2 + 0.2*randn(1, 6);
% Sort the prior ensemble for use below
sort_ens = sort(y_prior);

% Specifiy the likelihood
observation = 2;
obs_error_var = 1;

% Set of uniformly spaced horizontal points for doing plots
y = -3.5:0.05:3.5;

% Get the EAKF ensemble increments and prior and posterior points for plotting
[eakf_incs, eakf_prior_pts, eakf_post_pts, err] = obs_increment_eakf(y_prior, observation, obs_error_var, y);

% Get the gaussian kernel ensemble increments and prior and posterior points for plotting
[kernel_incs, kernel_prior_pts, kernel_post_pts, err] = obs_increment_kernel(y_prior, observation, obs_error_var, y);

% Get the binormal ensemble increments and prior and posterior points for plotting
[binorm_incs, binorm_prior_pts, binorm_post_pts, err] = obs_increment_binormal(y_prior, observation, obs_error_var, y);

% Get the particle increments 
[particle_incs, err] = obs_increment_particle(y_prior, observation, obs_error_var, y);

% Evaluate the normal likelihood at points for plotting
like = normpdf(y, observation, sqrt(obs_error_var));

% Establish colors for the three cases
colormap('parula');
my_map = colormap('parula');
c1 = my_map(1, :);
c2 = my_map(22, :);
c3 = my_map(43, :);

% Put on the zero lines below the different curves
bx = [min(y), max(y)];
by = [1 1];
plot(bx, by * 0, 'k');
hold on
plot(bx, by * -1, 'k');
plot(bx, by * -3, 'k');

% Plot the continuous priors at the top on their own axes
l_wid = 3;
ast_width = 1.5;
plot(y, eakf_prior_pts, '--', 'color', c1, 'linewidth', l_wid);
hold on
plot(y, binorm_prior_pts, 'color', c2, 'linewidth', l_wid);
plot(y, kernel_prior_pts', 'color', c3, 'linewidth', l_wid);

% Plot the likelihood
offset = 1;
plot(y, like - offset, 'color', c1, 'linewidth', l_wid) 

% Plot the continous analyses
offset = 3;
h_leg(1) = plot(y, eakf_post_pts - offset, '--', 'color', c1, 'linewidth', l_wid);
h_leg(2) = plot(y, binorm_post_pts - offset, 'color', c2, 'linewidth', l_wid);
h_leg(3) = plot(y, kernel_post_pts - offset', 'color', c3, 'linewidth', l_wid);

% Plot the prior ensemble just below the prior curves
ay = ones(size(y_prior));
plot(y_prior, ay * -.1, 'k*', 'linewidth', ast_width);

% Plot the various analysis ensembles
% Note that increments are added onto prior here for this example
plot(y_prior + eakf_incs, ay * -3.2, '*', 'color', c1,  'linewidth', ast_width); 
plot(y_prior + binorm_incs, ay * -3.4, '*', 'color', c2,  'linewidth', ast_width); 
plot(y_prior + kernel_incs, ay * -3.6, '*', 'color', c3,  'linewidth', ast_width); 
plot(y_prior + particle_incs, ay * -3.8, 'k*', 'linewidth', ast_width); 

% Plot a number below the posterior particles
% This is not automated if the ensemble is changed!!!
% Members 5(1), 6(2), 7(2), 8(2), 9(2), 10(1)
num_x = [1.2 1.5 1.8 2.1 2.4 2.7];
num_y = ones(size(num_x)) * -4.2;
number = [1 2 2 2 2 1];
for i = 1:6
   text(num_x(i), num_y(i), num2str(number(i)), 'fontsize', 14);
   % Draw a line from the number to the ensemble member
   cx = [num_x(i) + 0.02, sort_ens(i + 4)];
   cy = [num_y(i) + 0.1, -3.85];
   plot(cx, cy, 'k');
end

% Make things easily visible
pbaspect([1 1 1]);
set(gca, 'fontsize', 16, 'linewidth', 2);

set(gca, 'YTick', [-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1.0 1.5]);

tick_lab = get(gca, 'YTickLabel');
tick_lab = {0 0.5 1 1.5 0 0.5 0 0.5 1.0 1.5};
set(gca, 'YTickLabel', tick_lab);

axis([min(y), max(y), -4.5, 1.5]);
xlabel 'Observation';
ylabel 'Probability';

% Label the panels
text(-3.3, 1.2, 'Prior', 'Fontsize', 16);
text(-3.3, -0.8, 'Likelihood', 'Fontsize', 16);
text(-3.3, -2.8, 'Analysis', 'Fontsize', 16);

% Place a legend
h = legend(h_leg, 'EAKF', 'BiNormal', 'Kernel', 'location', 'w');
set(h, 'fontsize', 12)
