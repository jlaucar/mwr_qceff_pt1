% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% This script generates Figure 8.

% The continuous prior is a normal fit to the prior ensemble.
% Four different likelihoods are compared. 
% 1. Normal(1, 1)
% 2. Normal(7, 1)
% 3. 0.9*Normal(1, 1) + 0.1*Normal(1, sqrt(10))
% 4. 0.9*Normal(7, 1) + 0.1*Normal(7, sqrt(10))

ens_size = 10;

% Create the random draw for the prior ensemble
r_seed = 28;
rng(r_seed);
y_prior = randn(1, 10);

% Specify observed values and observation error variance
observation(1) = 1;
observation(2) = 7;
obs_error_var = 1;

% Set of uniformly spaced horizontal points for doing plots
y = -5:0.005:10;

% Get the EAKF ensemble increments and prior and posterior points for non-outlier observation
[eakf_incs, eakf_prior_pts, eakf_post_pts, err] = obs_increment_eakf(y_prior, observation(1), obs_error_var, y);
% Posterior ensemble is prior plus incrments
eakf_post = y_prior + eakf_incs;

% EAKF with outlier observation
[oeakf_incs, oeakf_prior_pts, oeakf_post_pts, err] = obs_increment_eakf(y_prior, observation(2), obs_error_var, y);
oeakf_post = y_prior + oeakf_incs;

% Fat tail sum of two normals likelihood with non-outlier observation
[fat_incs, fat_prior_pts, fat_like_pts, fat_post_pts, err] = ...
   obs_increment_double_like(y_prior, observation(1), obs_error_var, y);
fat_post = y_prior + fat_incs;

% Fat tail sum of two normals likelihood with outlier observation
[ofat_incs, ofat_prior_pts, ofat_like_pts, ofat_post_pts, err] = ...
   obs_increment_double_like(y_prior, observation(2), obs_error_var, y);
ofat_post = y_prior + ofat_incs;

% Evaluate the standard normal likelihoods at points for plotting
like = normpdf(y, observation(1), sqrt(obs_error_var));
olike = normpdf(y, observation(2), sqrt(obs_error_var));

% Establish colors for the cases
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
ast_width = 1.5;
plot(y, eakf_prior_pts, 'color', c1, 'linewidth', l_wid);
hold on

%  Plot the four likelihoods
offset = 1;
plot(y, like - offset, 'color', c2, 'linewidth', l_wid) 
plot(y, olike - offset, 'color', c3, 'linewidth', l_wid);
plot(y, fat_like_pts - offset, '--', 'color', c2, 'linewidth', l_wid) 
plot(y, ofat_like_pts - offset, '--', 'color', c3, 'linewidth', l_wid);

% Plot the four continuous analyses
offset = 2;
h_leg(1) = plot(y, eakf_post_pts - offset, 'color', c2, 'linewidth', l_wid);
h_leg(2) = plot(y, oeakf_post_pts - offset, 'color', c3, 'linewidth', l_wid);
h_leg(3) = plot(y, fat_post_pts - offset, '--', 'color', c2, 'linewidth', l_wid);
h_leg(4) = plot(y, ofat_post_pts - offset, '--', 'color', c3, 'linewidth', l_wid);

% Plot the prior ensemble just below the prior curves
ay = ones(size(y_prior));
plot(y_prior, ay * -.1, 'k*', 'linewidth', ast_width);

% Plot the various analysis ensembles
plot(y_prior + eakf_incs, ay * -2.1, '*', 'color', c2, 'linewidth', ast_width); 
plot(y_prior + oeakf_incs, ay * -2.2, '*', 'color', c3, 'linewidth', ast_width); 
plot(y_prior + fat_incs, ay * -2.3, '*', 'color', c2, 'linewidth', ast_width); 
plot(y_prior + ofat_incs, ay * -2.4, '*', 'color', c3, 'linewidth', ast_width); 

% Make things easily visible
pbaspect([1 1 1]);
set(gca, 'fontsize', 16, 'linewidth', 2);

tick_lab = get(gca, 'YTickLabel');
tick_lab = {0 0.5 0 0.5 0 0.5 1};
set(gca, 'YTickLabel', tick_lab);

axis([-2, max(y), -2.5, 0.5]);
xlabel 'Observation';
ylabel 'Probability';

% Label the panels
text(8, 0.2, 'Prior', 'Fontsize', 16);
text(7, -0.4, 'Likelihood', 'Fontsize', 16);
text(7.5, -1.8, 'Analysis', 'Fontsize', 16);

% Place a legend
h = legend(h_leg, 'Not Outlier', 'Outlier', 'Not Outlier, heavy tail', 'Outlier, heavy tail');
set(h, 'fontsize', 12)






