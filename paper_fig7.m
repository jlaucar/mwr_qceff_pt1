% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% This script generates Figure 7.

% The continuous prior is a normal while the likelihood is an exponential.

ens_size = 10;

% Create the random draw for the prior ensemble
r_seed = 12;
rng(r_seed);
y_prior = normrnd(2, 1, 1, ens_size);

% Define the likelihood, a single parameter lambda for the exponential
like_lambda = 1;

% Set of uniformly spaced horizontal points for doing plots
y = -3:0.001:10;

% Get the Normal/Exponential ensemble increments and prior and posterior points for plotting
[ne_incs, ne_prior_pts, ne_post_pts, ne_like_pts, err] = obs_increment_normal_exp(y_prior, like_lambda, y);
% Posterior ensemble is prior plus increments
ne_post = y_prior + ne_incs;

% Put on the zero lines below the different curves
bx = [min(y), max(y)];
by = [1 1];
plot(bx, by * 0, 'k');
hold on
plot(bx, by * -1, 'k');
plot(bx, by * -1.75, 'k');

% Plot the continuous prior at the top on its own axes
l_wid = 3;
ast_width = 1.5;
plot(y, ne_prior_pts, 'b', 'linewidth', l_wid);
hold on

% Plot the likelihood
offset = 1;
plot(y, ne_like_pts - offset, 'b', 'linewidth', l_wid) 

% Plot the continuous analysis
offset = 1.75;
plot(y, ne_post_pts - offset, 'b', 'linewidth', l_wid);

% Plot the prior ensemble just below the continuous prior
ay = ones(size(y_prior));
plot(y_prior, ay * -.1, 'k*', 'linewidth', ast_width);

% Plot the posterior ensemble
plot(y_prior + ne_incs, ay * -1.9, 'b*', 'linewidth', ast_width); 

% Make things easily visible
pbaspect([1.5 1 1]);
set(gca, 'fontsize', 16, 'linewidth', 2);

set(gca, 'YTick', [-1.75 -1.25 -1 -0.5 0 0.5]);

tick_lab = get(gca, 'YTickLabel');
tick_lab = {0 0.5 0 0.5 0 0.5};
set(gca, 'YTickLabel', tick_lab);

axis([-1, 6, -2, 0.5]);
xlabel 'Observation';
ylabel 'Probability';

% Label the panels
text(5, 0.3, 'Prior', 'Fontsize', 16);
text(4.5, -0.6, 'Likelihood', 'Fontsize', 16);
text(4.5, -1.4, 'Analysis', 'Fontsize', 16);
