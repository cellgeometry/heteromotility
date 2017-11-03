function [div_stats] = pfa_divergence(experiment, experiment_name, state_vectors, binsize)

%% Make x and y matrices of appropriate size

m_x = zeros(binsize(1), binsize(2));
m_y = zeros(binsize(1), binsize(2));
m_u = zeros(binsize(1), binsize(2));
m_v = zeros(binsize(1), binsize(2));
m_c = zeros(binsize(1), binsize(2));

for i = 1:binsize(2);
    m_x(:,i) = i;
end

for i = 1:binsize(1);
    m_y(i,:) = i;
end



%% Place vector values into appropriate matrices

for row = 1:length(state_vectors);
    m_u(state_vectors(row,2), state_vectors(row,1)) = state_vectors(row,3);
end

for row = 1:length(state_vectors);
    m_v(state_vectors(row,2), state_vectors(row,1)) = state_vectors(row,4);
end

for row = 1:length(state_vectors);
    m_c(state_vectors(row,2), state_vectors(row,1)) = state_vectors(row,5);
end

%% Calculate divergence, plot field and divergence

div = divergence(m_x, m_y, m_u, m_v);

f1 = figure;
colormap(redblue(256));
imagesc(div);
h = colorbar;
set(h,'fontsize',24);
hold on;
quiver(m_x, m_y, m_u, m_v, 'k');
caxis([-1.5,1.5]);
hold off;
saveas(f1, strcat('divergence_figs/', experiment, '_divergence.fig'), 'fig');
saveas(f1, strcat('divergence_figs/', experiment, '_divergence.png'), 'png');

f2 = figure;
surf(div);
hold on;
colormap(redblue(256));
caxis([-1.5,1.5]);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'ZTickLabel', []);
xlabel('PC1', 'FontSize', 26);
ylabel('PC2', 'FontSize', 26);
zlabel('Divergence', 'FontSize', 26);
h = colorbar;
set(h,'fontsize',24);
quiver(m_x, m_y, m_u, m_v, 'k');
hold off;
saveas(f2, strcat('divergence_figs/', experiment, '_divergence_surface.fig'), 'fig');

f3 = figure;
[xx, yy] = meshgrid(0:0.1:15);
surf(interp2(div, xx, yy));
shading interp;
grid off;
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'ZTickLabel', [])
colormap('jet');
caxis([-1.5,1.5]);
h = colorbar;
set(h,'fontsize',24);
xlabel('PC1', 'FontSize', 26);
ylabel('PC2', 'FontSize', 26);
zlabel('Divergene', 'FontSize', 26);
saveas(f3, strcat('divergence_figs/', experiment, '_divergence_surface_interp.fig'), 'fig');


% Save a surface map of the number of cells in each state bin
f4 = figure;
% flipping handles different defaults for quiver vs surf axes ordering
surf(interp2(fliplr(log(m_c+1)), fliplr(xx), flipud(yy))); 
shading interp;
grid off;
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'ZTickLabel', [])
colormap('jet');
h = colorbar;
set(h,'fontsize',24);
xlabel('PC1', 'FontSize', 26);
ylabel('PC2', 'FontSize', 26);
zlabel('log(Cell Number)', 'FontSize', 26);
saveas(f4, strcat('divergence_figs/', experiment, '_cell_occupancy_surface_interp.fig'), 'fig');

%% Find metastable states, write divergence stats


metastable = find( div < 0); % negative divergence ~ stable
prop_stable = length(metastable) / binsize(1)^2;

[stable, I] = min(div);
[s_i, s_j] = ind2sub(size(div), I);
m_s = [s_i, s_j];


% Count metastable states

num_ms = length(metastable);
avg_div = mean(div(:));

% Calculate mean transition magnitude

v_mags = sqrt(m_u.^2 + m_v.^2);
gt0_vmags = v_mags(v_mags>0);
avg_v_mag = mean(gt0_vmags);

%% Calculate Quantitiative metrics of 2D PFA

% Mean Vector Magnitude (for unit vector seperated displacements!!) 
m_w = m_c.*v_mags; % find weighted vector magnitudes
mean_mag = sum(m_w(:))/sum(m_c(:)); % only consider bins containing cells at t0

% Proportion of cells in metastable states at t0
p_metastable_t0 = sum(m_c(metastable))/sum(m_c(:));

% Directedness of Vector Field

dir_v = [ mean(m_u(m_c > 0)), mean(m_v(m_c > 0))];
mag_dir = sqrt(sum(dir_v.^2));

% Write out m_u m_v for statistical testing
mag_mat = sqrt(m_u.^2 + m_v.^2); % make bin x bin matrix of vector mags
mag_matw = mag_mat .* m_c; % weight by number of cells
mag_matw = mag_matw ./ sum(m_c(:)); % sum of mag_matw is now the mean!
mag_vec = mag_mat(m_c > 0); % linearize for output, only save >0 bins
mag_vecw = mag_matw(m_c > 0);
cells_vec = m_c(m_c > 0);
dir_ds = dataset({[mag_vec mag_vecw cells_vec] 'bin_vmag', 'bin_vmag_w', 'cells'});
export(dir_ds, 'file', strcat('divergence_figs/', experiment, '_bin_vmag.csv'), 'delim', ',');

%% Write out divergence field statistics

% Mean Divergence, # of MS states, prop of cells in MS state at t0,
% magnitude of vector sum of rate constants
div_stats = [avg_div, num_ms, p_metastable_t0, mag_dir];
ds = dataset({div_stats 'avg_div' 'num_ms' 'p_ms_cells_t0' 'mag_dir'});
export(ds, 'file', strcat('divergence_figs/', experiment, '_div_stats.csv'), 'delim', ',');

%% Determine if vector magnitudes are increased farther from MS states
% Checks every point in the divergence matrix, finds distance from 
% nearest metastable state and magnitude of flux in that state
% plots distance to metastable vs flux mag 
% performs lin regression to describe relationship


dvm_all = zeros(length(div(:)), 2);
for k = 1:length(div(:));
    [p_i, p_j] = ind2sub(size(div), k);
    p = [p_i, p_j];
    res = point_distmag(p, m_u, m_v, metastable);
    dvm_all(k,:) = res;
end

gt0_ind = dvm_all(:,2) > 0;
dvm_gt0 = dvm_all(gt0_ind,:);
sorted_dvm = sortrows(dvm_gt0, 1);

p = polyfit(sorted_dvm(:,1), sorted_dvm(:,2), 1); 
% p(1) slope, p(2) intercept
yfit = polyval(p,sorted_dvm(:,1));
yresid = sorted_dvm(:,2) - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(sorted_dvm(:,2))-1) * var(sorted_dvm(:,2));
rsq = 1 - SSresid/SStotal;
lin_result = [p(1), p(2), rsq];
% Writes out slope, intercept, rsq of linear fit
csvwrite( strcat('divergence_figs/', experiment, '_distvmag_fit.csv'), lin_result);

f4 = figure;
plot(sorted_dvm(:,1), sorted_dvm(:,2));
hold on;
plot(sorted_dvm(:,1), yfit);
title(strcat(experiment_name, ' - ', 'Flux Magnitude vs Distance from Stable State'));
xlabel('Distance from Metastable State [state units]');
ylabel('Flux Magnitude');
saveas(f4, strcat('divergence_figs/', experiment, 'dist_v_mag.fig'), 'fig');
hold off;

%% Close all plots

close all;

end