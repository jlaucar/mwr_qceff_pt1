% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% This script generates Figure 5.

% The ensemble prior is a draw from a beta distribution, the continuous 
% prior is beta and the likelihood is beta.

ens_size = 10;

% Create the random draw for the prior ensemble
r_seed = 14;
rng(r_seed);
prior_alpha = 0.5;
prior_beta = 0.5;
y_prior = betarnd(prior_alpha, prior_beta, 1, ens_size);

% Define the likelihood parameters
like_alpha = 2;
like_beta = 5;

% Set of uniformly spaced horizontal points for doing plots
y = 0.001:0.001:0.999;

% Get the beta/beta ensemble increments and prior and posterior points for plotting
% Also returns the points for plotting the likelihood
[bb_incs, bb_prior_pts, bb_post_pts, bb_like_pts, err] = ...
   obs_increment_beta_beta(y_prior, like_alpha, like_beta, y);
% Posterior is prior plus increments
bb_post = y_prior + bb_incs;

% Put on the zero lines below the different curves
bx = [min(y), max(y)];
by = [1 1];
plot(bx, by * 0, 'k');
hold on
plot(bx, by * -3, 'k');
plot(bx, by * -6, 'k');

% Plot the continuous prior at the top on its own axes
l_wid = 3;
ast_width = 1.5;
plot(y, bb_prior_pts, 'b', 'linewidth', l_wid);
hold on

% Plot the continuous likelihood
offset = 3;
plot(y, bb_like_pts - offset, 'b', 'linewidth', l_wid) 

% Plot the continuous analysis
offset = 6;
plot(y, bb_post_pts - offset, 'b', 'linewidth', l_wid);

% Plot the prior ensemble just below the continuous prior
ay = ones(size(y_prior));
plot(y_prior, ay * -.1, 'k*', 'linewidth', ast_width);

% Plot the analysis ensemble
plot(y_prior + bb_incs, ay * -6.1, 'b*', 'linewidth', ast_width); 

% Make things easily visible
pbaspect([1.5 1 1]);
set(gca, 'fontsize', 16, 'linewidth', 2);

set(gca, 'YTick', [-6 -5 -4 -3 -2 -1 0 1 2 3 4]);

tick_lab = get(gca, 'YTickLabel');
tick_lab = {0 1 2 0 1 2  0 1 2 3 4};
set(gca, 'YTickLabel', tick_lab);

axis([0, 1, -6.5, 4.1]);
xlabel 'Observation';
ylabel 'Probability';

% Label the panels
text(0.8, 3.5, 'Prior', 'Fontsize', 16);
text(0.8, -2, 'Likelihood', 'Fontsize', 16);
text(0.8, -5, 'Analysis', 'Fontsize', 16);
