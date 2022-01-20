% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download

% This script was used with Matlab 2016b to generate figures for
% A Quantile Conserving Ensemble Filter Framework. Part I: Updating an Observed Variable
% by Jeffrey Anderson
% which was submitted to Monthly Weather Review.


% This script generates Figure 3 which compares 3 QCEF filters applied to the
% same prior ensemble:
% 1. An EAKF where continuous prior and likelihood are normal
% 2. A truncated normal prior and normal likelihood
% 3. A log-normal prior and log-normal likelihood

ens_size = 10;

% Create the random draw for the prior ensemble
r_seed = 12;
rng(r_seed);
log_y_prior = randn(1, ens_size);
% Make it a draw from a log-normal distribution
y_prior = exp(log_y_prior);

% Specify an observed value and observation errof variance
log_observation = -2;
observation = exp(log_observation);
log_obs_error_var = 4;
log_obs_error_sd = sqrt(log_obs_error_var);
% Specify observation error variance for the normal likelihood approximation
obs_error_sd = 1.0;
obs_error_var = obs_error_sd^2;

% Set of uniformly spaced horizontal points for doing plots
y = -1:0.001:6.2;

% Get the EAKF ensemble increments and prior and posterior points for plotting
[eakf_incs, eakf_prior_pts, eakf_post_pts, err] = obs_increment_eakf(y_prior, observation, obs_error_var, y);
% Posterior ensemble is prior plus increments
eakf_post = y_prior + eakf_incs;

% Get the bounded positive normal ensemble increments and prior and posterior points for plotting
[pos_eakf_incs, pos_eakf_prior_pts, pos_eakf_post_pts, err] = obs_increment_pos_eakf(y_prior, observation, obs_error_var, y);
% Posterior ensemble is prior plus increments
pos_eakf_post = y_prior + pos_eakf_incs;

% Get the lognormal prior, lognormal likelihood ensemble increments 
% and prior and posterior points for plotting
[ln_eakf_incs, ln_eakf_prior_pts, ln_eakf_post_pts, err] = obs_increment_ln_eakf(y_prior, log_observation, log_obs_error_var, y);
% Posterior ensemble is prior plus increments
ln_eakf_post = y_prior + ln_eakf_incs;

% Evaluate the normal likelihood at points for plotting
like = normpdf(y, observation, sqrt(obs_error_var));

% Evaluate the lognormal likelihood at points for plotting
for i = 1:size(y, 2)
   if(y(i) <= 0) 
      ln_like_pts(i) = 0;
   else
      ln_like_pts(i) = normpdf(log(y(i)), log_observation, log_obs_error_sd);
   end
end

% Normalize the likelihood for display
like_area = sum(ln_like_pts) * (y(2) - y(1));
ln_like_pts = ln_like_pts / like_area;

% Normalize the log likelihood for display
prior_area = sum(ln_eakf_prior_pts) * (y(2) - y(1));
ln_eakf_prior_pts = ln_eakf_prior_pts / prior_area;
post_area = sum(ln_eakf_post_pts) * (y(2) - y(1));
ln_eakf_post_pts = ln_eakf_post_pts / post_area;

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
plot(bx, by * -0.75, 'k');
plot(bx, by * -1.5, 'k');

% Plot the continuous priors at the top on their own axes
l_wid = 3;
ast_width = 1.5;
plot(y, eakf_prior_pts, '--', 'color', c1, 'linewidth', l_wid);
hold on
plot(y, pos_eakf_prior_pts, 'color', c2, 'linewidth', l_wid);
plot(y, ln_eakf_prior_pts', 'color', c3, 'linewidth', l_wid);

% Plot the continuous likelihoods
offset = 0.75;
plot(y, like - offset, '--', 'color', c1, 'linewidth', l_wid);
plot(y, ln_like_pts - offset, 'color', c3, 'linewidth', l_wid);

% Plot the continuous analysis
offset = 1.5;
h_leg(1) = plot(y, eakf_post_pts - offset, '--', 'color', c1, 'linewidth', l_wid);
h_leg(2) = plot(y, pos_eakf_post_pts - offset, 'color', c2, 'linewidth', l_wid);
h_leg(3) = plot(y, ln_eakf_post_pts - offset, 'color', c3, 'linewidth', l_wid);

% Plot the prior ensemble just below the continuous prior
ay = ones(size(y_prior));
plot(y_prior, ay * -.1, 'k*', 'linewidth', ast_width);

% Plot the various analysis ensembles
plot(y_prior + eakf_incs, ay * -1.6, '*', 'color', c1, 'linewidth', ast_width); 
plot(y_prior + pos_eakf_incs, ay * -1.7, '*', 'color', c2, 'linewidth', ast_width); 
plot(y_prior + ln_eakf_incs, ay * -1.8, '*', 'color', c3, 'linewidth', ast_width); 

% Make things easily visible
pbaspect([1 1 1]);
set(gca, 'fontsize', 16, 'linewidth', 2);

set(gca, 'YTick', [-1.5 -1 -0.75 -0.25 0 0.5]);

tick_lab = get(gca, 'YTickLabel');
tick_lab = {0 0.5 0 0.5 0 0.5};
set(gca, 'YTickLabel', tick_lab);

axis([min(y), max(y), -1.9, 0.5]);
xlabel 'Observation';
ylabel 'Probability';

% Label the panels
text(5, 0.4, 'Prior', 'Fontsize', 16);
text(4.5, -0.4, 'Likelihood', 'Fontsize', 16);
text(4.5, -0.95, 'Analysis', 'Fontsize', 16);

% Place a legend
h = legend(h_leg, 'EAKF', 'Truncated Normal', 'Log-Normal', 'location', 'nw');
set(h, 'fontsize', 12)
