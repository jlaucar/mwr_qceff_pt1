% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% This script generates Figure 2.

% The ensemble prior is a draw from a normal distribution
% Three different continuous priors are used: Normal, a Rank Histogram Filter, 
% and a sum of 10 Gaussian kernels

ens_size = 10;

% Create the random draw for the prior ensemble
r_seed = 28;
rng(r_seed);
y_prior = randn(1, ens_size);

% Specify an observed value and observation error variance
observation = 1;
obs_error_var = 1;

% Set of uniformaly spaced horizontal points for doing plots
y = -2.5:0.005:3.5;

% Get the EAKF ensemble increments and prior and posterior points for plotting
[eakf_incs, eakf_prior_pts, eakf_post_pts, err] = obs_increment_eakf(y_prior, observation, obs_error_var, y);
% Posterior ensemble is prior plus increments
eakf_post = y_prior + eakf_incs;

% Get the RHF ensemble increments and prior and posterior points for plotting
% Also returns the piecewise constant approximate RHF likelihood
[rhf_incs, rhf_prior_pts, rhf_post_pts, rhf_like_pts, err] = obs_increment_rhf(y_prior, observation, obs_error_var, y);
% Posterior ensemble is prior plus increments
rhf_post = y_prior + rhf_incs;

% Get the gaussian kernel ensemble increments and prior and posterior points for plotting
[kernel_incs, kernel_prior_pts, kernel_post_pts, err] = obs_increment_kernel(y_prior, observation, obs_error_var, y);
% Posterior ensemble is prior plus increments
kernel_post = y_prior + kernel_incs;

% Evaluate the normal likelihood at points for plotting
like = normpdf(y, observation, sqrt(obs_error_var));
% Rescale the piecewise constant RHF likelihood for plotting
rhf_like_pts = rhf_like_pts / max(rhf_like_pts) * max(like);

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
plot(bx, by * -2, 'k');

% Plot the continuous priors at the top on their own axes
l_wid = 3;
l_wid2 = 3;
ast_width = 1.5;
plot(y, eakf_prior_pts, '--', 'color', c1, 'linewidth', l_wid);
hold on
plot(y, rhf_prior_pts, 'color', c2, 'linewidth', l_wid);
plot(y, kernel_prior_pts', 'color', c3, 'linewidth', l_wid2);

% Plot the continuous and piecewise continuous likelihoods
offset = 1;
plot(y, like - offset, '--', 'color', c1, 'linewidth', l_wid) 
plot(y, rhf_like_pts - offset, 'color', c2, 'linewidth', l_wid);

% Plot the continuous analyses
offset = 2;
h_leg(1) = plot(y, eakf_post_pts - offset, '--', 'color', c1, 'linewidth', l_wid);
h_leg(2) = plot(y, rhf_post_pts - offset, 'color', c2, 'linewidth', l_wid);
h_leg(3) = plot(y, kernel_post_pts - offset', 'color', c3, 'linewidth', l_wid2);

% Plot the prior ensemble just below the continuous prior 
ay = ones(size(y_prior));
plot(y_prior, ay * -.1, 'k*', 'linewidth', ast_width);

% Plot the various analysis ensembles
plot(y_prior + eakf_incs, ay * -2.1, '*', 'color', c1, 'linewidth', ast_width); 
plot(y_prior + rhf_incs, ay * -2.2, '*', 'color', c2, 'linewidth', ast_width); 
plot(y_prior + kernel_incs, ay * -2.3, '*', 'color', c3, 'linewidth', ast_width); 

% Make things easily visible
pbaspect([1 1 1]);
set(gca, 'fontsize', 16, 'linewidth', 2);

tick_lab = get(gca, 'YTickLabel');
tick_lab = {0 0.5 0 0.5 0 0.5 1};
set(gca, 'YTickLabel', tick_lab);

axis([-2, 3, -2.4, 1]);
xlabel 'Observation';
ylabel 'Probability';

% Label the panels
text(-1.9, 0.2, 'Prior', 'Fontsize', 16);
text(-1.9, -0.8, 'Likelihood', 'Fontsize', 16);
text(-1.9, -1.8, 'Analysis', 'Fontsize', 16);

% Place a legend
h = legend(h_leg, 'EAKF', 'RHF', 'Kernel', 'location', 'nw');
set(h, 'fontsize', 12)
