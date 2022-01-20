% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% This script generates Figure 3.

% Computes QCEF for a continuous gamma prior and gamma likelihood

% This is the second case from Bishop 2016 (see paper text for details).

% The exact prior is gamma with prior_mean = 1; type1 variance = 1;
% This means shape k = 1 and scale theta = 1
true_prior_shape = 1;
true_prior_scale = 1;

ens_size = 10;

% Create the random draw for the prior ensemble
r_seed = 12;
rng(r_seed);
y_prior = gamrnd(1, 1, 1, ens_size);

% Compute the sample statistics for this ensemble, following Bishop
prior_mean = mean(y_prior);
prior_type1_var = var(y_prior) / prior_mean^2;
prior_shape = 1 / prior_type1_var;
prior_scale = prior_mean * prior_type1_var;

% Define the observation as per Bishop
% Start with his specifications for the inverse gamma observation error PDF
% This is observation 3, type_1 variance 1/4, type_2 variance 1/5
% This leads to alpha = 6, gamma = 15
% The shape = alpha+1; the scale = obs * t2_var
like_shape = 7;
like_scale = 0.6;

% Compute the continuous posterior shape and scale (Appendix a)
post_shape = prior_shape + like_shape - 1;
post_scale = (prior_scale * like_scale) / (prior_scale + like_scale);

% Set of uniformly spaced horizontal points for doing plots
y = 0:0.001:10;

% Get the gamma/gamma ensemble increments and posterior points for plotting
[gg_incs, gg_prior_pts, gg_post_pts, err] = obs_increment_gamma_gamma(y_prior, like_shape, like_scale, y);
% Posterior ensemble is prior plus increments
gg_post = y_prior + gg_incs;

% Evaluate the continuous gamma prior and likelihood at points for plotting
gg_prior_pts = gampdf(y, prior_shape, prior_scale);
gg_post_pts = gampdf(y, post_shape, post_scale);
like = gampdf(y, like_shape, like_scale);

% Put on the zero lines below the different curves
bx = [min(y), max(y)];
by = [1 1];
plot(bx, by * 0, 'k');
hold on
plot(bx, by * -0.5, 'k');
plot(bx, by * -1.25, 'k');

% Plot the continuous prior at the top on its own axis
l_wid = 3;
ast_width = 1.5;
plot(y, gg_prior_pts, 'b', 'linewidth', l_wid);
hold on

% Plot the continuous likelihood
offset = 0.5;
plot(y, like - offset, 'b', 'linewidth', l_wid) 

% Plot the continuous analysis
offset = 1.25;
plot(y, gg_post_pts - offset, 'b', 'linewidth', l_wid);

% Plot the prior ensemble just below the continuous prior
ay = ones(size(y_prior));
plot(y_prior, ay * -.1, 'k*', 'linewidth', ast_width);

% Plot the analysis ensemble
plot(y_prior + gg_incs, ay * -1.3, 'b*', 'linewidth', ast_width); 

% Make things easily visible
pbaspect([1 1 1]);
set(gca, 'fontsize', 16, 'linewidth', 2);

set(gca, 'YTick', [-1.25 -0.75 -0.5 -0.25 0 0.5 1.0]);


tick_lab = get(gca, 'YTickLabel');
tick_lab = {0 0.5 0 0.25 0 0.5 1.0};
set(gca, 'YTickLabel', tick_lab);

axis([min(y), 8, -1.4, 1.1]);
xlabel 'Observation';
ylabel 'Probability';


% Label the panels
text(7, 0.9, 'Prior', 'Fontsize', 16);
text(6, -0.2, 'Likelihood', 'Fontsize', 16);
text(6, -0.9, 'Analysis', 'Fontsize', 16);
