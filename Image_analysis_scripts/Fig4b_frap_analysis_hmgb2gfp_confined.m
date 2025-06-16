clear all
close all

% FRAP analysis - HMGB2-GFP in confined cells

data_dir = '/Volumes/whitelab/Lab Members/MirandaHunter/Microscopy/LSM880/2023/230526_frap_HMGB2gfp/for_analysis/';

tp_before_frap = 1; % number of images acquired before FRAP

t_interval = 0.2; % in seconds

bleach_cutoff = 0.25; % percent bleaching cutoff, ie throw out all cells that weren't bleached below this threshold


%% Import data

cd(data_dir);
files = dir;
filenames = {files.name};
sample_data = filenames(endsWith(filenames, 'txt'));

n_cells = length(sample_data);


for ii = 1:n_cells

    recovery_data = readtable(sample_data{ii});
    sample_name = append(extractBetween(sample_data{ii}, 'GFP_', '.txt'), '_', extractBefore(sample_data{ii}, '_FRAP'));
    recovery_data.Sample_name(:) = sample_name;
    replicate = extractAfter(sample_name, '_n');
    recovery_data.Replicate(:) = extractAfter(replicate, '_');

    % calculate percent bleach
    percent_bleach = recovery_data.Intensity(tp_before_frap+1) ./ recovery_data.Intensity(1);

    if percent_bleach > bleach_cutoff
        continue
    else
        recovery_data.Bleach_per(:) = percent_bleach*100;

        % normalize to pre-bleach fluorescence
        recovery_data.Intensity_norm = recovery_data.Intensity ./ recovery_data.Intensity(1);

        % normalize to pre AND post bleach fluorescence (the way done in Anna's paper)
        intensity_pre = recovery_data.Intensity(recovery_data.Time == tp_before_frap);
        intensity_post = recovery_data.Intensity(recovery_data.Time == (tp_before_frap+1));

        recovery_data.Intensity_norm_both = ((recovery_data.Intensity - intensity_post) ./ (intensity_pre - intensity_post)) *100;

        if ii == 1
            frap_results_table = recovery_data(1:100,:); % remove the last row which for some reason is all nans
        else
            frap_results_table = [frap_results_table; recovery_data(1:100,:)];
        end
    end

end

frap_results_table.Time_s_corr = (frap_results_table.Time-2)*t_interval;
frap_results_table.Group = extractBefore(frap_results_table.Sample_name, '_n');

%% Average across groups

n_timepoints = max(frap_results_table.Time);

groups = unique(extractBefore(frap_results_table.Sample_name, '_n'));
replicates = unique(frap_results_table.Replicate);

recovery_m = nan(n_timepoints, length(groups));
recovery_s = nan(n_timepoints, length(groups));
recovery_e = nan(n_timepoints, length(groups));

for ii = 1:length(groups)

    for jj = 1:n_timepoints
        
        tp_data = frap_results_table(frap_results_table.Time == jj & startsWith(frap_results_table.Sample_name, groups{ii}),:);
        [recovery_m(jj,ii), recovery_s(jj,ii), recovery_e(jj,ii)] = rmsnan(tp_data.Intensity_norm_both);

    end
end

%% Fit exponential

timeconst1 = nan(100,2);
timeconst2 = nan(100,2);
A1 = nan(100,2);
A2 = nan(100,2);
A1_per = nan(100,2);
A2_per = nan(100,2);
Rsq_all = nan(100,2);
mobile_fraction_exp = nan(100,2);

show_fits = 0;

for ii = 1:length(groups)

    group_data = frap_results_table(startsWith(frap_results_table.Sample_name, groups{ii}),:);
    cell_names = unique(group_data.Sample_name);

    %figure; tiledlayout(6,9);
    for jj = 1:length(cell_names)

        % find intensity at last TP
        recovery_data = group_data((strcmp(group_data.Sample_name, cell_names{jj}) & group_data.Time > 1),:);

        [f,results] = fit(recovery_data.Time_s_corr, recovery_data.Intensity_norm_both, 'exp2');
        rsq = results.rsquare;

        if show_fits
            figure;
            %nexttile;
            plot(f, recovery_data.Time_s_corr, recovery_data.Intensity_norm_both);
            set(gca, 'FontSize', 20);
            yline(f.a, 'k--', 'DisplayName', 'fast/slow diffusing cutoff');
            %set(gca, 'XTick', [], 'XTickLabels', [], 'YTick', [], 'YTickLabels', [], 'Xcolor', 'none', 'YColor', 'none');
            legend('location', 'southeast')
            xlabel('Time (s)')
            ylabel('% recovery')
            subtitle(append('Rsq = ', string(rsq)));
        end
        
        mf1 = recovery_data(recovery_data.Time == n_timepoints,:).Intensity_norm_both;
        mf2 = recovery_data(recovery_data.Time == (n_timepoints-1),:).Intensity_norm_both;
        mf_exp = (mf1+mf2) ./2;
        mobile_fraction_exp(jj,ii) = mf_exp;


        % b and d are the time constants (b = fast diffusion, d = slow diffusion)
        timeconst1(jj,ii) = f.b;
        timeconst2(jj,ii) = f.d;

        A1(jj,ii) = f.a;
        A2(jj,ii) = mf_exp - f.a;

        A1_per(jj,ii) = (f.a / mf_exp) *100;
        A2_per(jj,ii) = 100 - ((f.a / mf_exp) *100);
        Rsq_all(jj,ii) = rsq;

    end
end

ratio = A1 ./ A2;

cols = [136 136 136; 154 41 134]./256;


% plot % slow diffusing
figure;
boxplotMVH(flip(A2_per,2), cols, flip(groups))
ylabel({'% slow diffusing HMGB2'})
ylim([0 15])

% % plot as bar graph
% [m_A2 s_A2 e_A2] = rmsnan(flip(A2_per,2));
% [m_A1 s_A1 e_A1] = rmsnan(flip(A1_per,2));
% 
% figure;
% hold on
% bar([m_A1(1), m_A2(1); m_A1(2), m_A2(2)], 'stacked');



%% Plots

% individual curves for each group

n_cells_per_group = sum(~isnan(mobile_fraction_exp));

for ii = 1:length(groups)

    group_data = frap_results_table(startsWith(frap_results_table.Sample_name, groups{ii}),:);
    cell_names = unique(group_data.Sample_name);

    switch ii
        case 1 % confined
            cols = flip(magma(n_cells_per_group(ii)));

            figure;
            tiledlayout(1,3);
            nexttile;
            hold on;

        case 2 % unconfined
            cols = flip(viridis(n_cells_per_group(ii)));
            nexttile;
            hold on;
    end

    for jj = 1:length(cell_names)

        % find intensity at last TP
        recovery_data = group_data(strcmp(group_data.Sample_name, cell_names{jj}),:);

        plot([0:99.5]*0.2, recovery_data.Intensity_norm_both, '.-',...
            'MarkerFaceColor', cols(jj,:),...
            'MarkerEdgeColor', cols(jj,:),...
            'Color', cols(jj,:),...
            'LineWidth',2);

    end
    xlim([-0.5 20]);
    ylim([0 103]);
    set(gca, 'Box', 'off', 'FontSize', 26, 'LineWidth', 3);
    xlabel('time (s)');
    ylabel('% recovery');
    title(groups{ii});
end


% plot recovery over time - normalized (SEM)
%figure;
nexttile;
hold on
lineProps1.col{1} = [136 136 136]./256; 
lineProps1.width = 5;

lineProps2.col{1} = [154 41 134]./256; 
lineProps2.width = 5;

mseb([-1:98.5]*t_interval, recovery_m(:,2)', recovery_e(:,2)', lineProps1, 0); % unconfined
mseb([-1:98.5]*t_interval, recovery_m(:,1)', recovery_e(:,1)', lineProps2, 0); % confined

legend({'unconfined', 'confined'}, 'location', 'southeast', 'fontsize', 24, 'box', 'off')
xlim([-1 20]);
ylim([0 103]);
set(gca, 'Box', 'off', 'FontSize', 26, 'LineWidth', 3);
xlabel('time (s)');
ylabel('% recovery');

figure;
mseb([-1:98.5]*t_interval, recovery_m(:,2)', recovery_e(:,2)', lineProps1, 1); % unconfined
xlim([-1 20]);

figure;
mseb([-1:98.5]*t_interval, recovery_m(:,1)', recovery_e(:,1)', lineProps2, 1); % confined
xlim([-1 20]);


% recovery over time - normalized (SD)
figure;
hold on
lineProps1.col{1} = [136 136 136]./256; 
lineProps1.width = 5;

lineProps2.col{1} = [154 41 134]./256; 
lineProps2.width = 5;

mseb([-1:98.5]*t_interval, recovery_m(:,2)', recovery_s(:,2)', lineProps1, 1); % unconfined
mseb([-1:98.5]*t_interval, recovery_m(:,1)', recovery_s(:,1)', lineProps2, 1); % confined

legend({'unconfined', 'confined'}, 'location', 'southeast', 'fontsize', 20, 'box', 'off')
xlim([-1 20]);
ylim([0 103]);
set(gca, 'Box', 'off', 'FontSize', 26, 'LineWidth', 3);
xlabel('time (s)');
ylabel('% recovery');
    

% % stats
% [~,P_mf] = ttest2(mobile_fraction_all(:,1), mobile_fraction_all(:,2));
% 
% figure;
% boxplotMVH(flip(mobile_fraction_all,2), [136 136 136; 154 41 134]./256, flip(groups))
% ylim([0 100])
% ylabel('mobile fraction (%)')
% title('mobile fraction')
% subtitle(append('P = ', string(P_mf)))



%% HALF-TIME

% t1/2 with last time point

n_cells_per_group = sum(~isnan(mobile_fraction_exp));

htime_end = nan(100,2);

for ii = 1:length(groups)

    group_data = frap_results_table(startsWith(frap_results_table.Sample_name, groups{ii}),:);
    cell_names = unique(group_data.Sample_name);

    for jj = 1:length(cell_names)

        % find intensity at last TP
        recovery_data = group_data(strcmp(group_data.Sample_name, cell_names{jj}),:);
        intensity_lasttp = max(recovery_data.Intensity);
        intensity_firsttp = min(recovery_data.Intensity); % since we are not bleaching down to 0, use min intensity to calculate t1/2
        hmax = (intensity_lasttp - intensity_firsttp)/2 + intensity_firsttp; % half max intensity

        ind = recovery_data((recovery_data.Intensity >= hmax & recovery_data.Time > 1),:); % exclude pre-bleach tp
        htime_end(jj,ii) = min(ind.Time);
    end
end    

[~, P_t] = ttest2(htime_end(:,1), htime_end(:,2));

% plot
figure;
boxplotMVH(flip(htime_end,2)*0.2, [136 136 136; 154 41 134]./256, flip(groups))
ylim([0 3])
ylabel('half-time (s)')
title('half-time')
subtitle(append('P = ', string(P_t)))


%% Plot per replicate

for kk = 1:length(replicates)

    figure;
    tiledlayout(1,2);

    for ii = 1:length(groups)

        group_data = frap_results_table(startsWith(frap_results_table.Sample_name, groups{ii}) & strcmp(frap_results_table.Replicate, replicates{kk}),:);
        cell_names = unique(group_data.Sample_name);

        switch ii
            case 1 % confined
                cols = flip(magma(n_cells_per_group(ii)));
                nexttile;
                hold on;

            case 2 % unconfined
                cols = flip(viridis(n_cells_per_group(ii)));
                nexttile;
                hold on;
        end

        for jj = 1:length(cell_names)

            % find intensity at last TP
            recovery_data = group_data(strcmp(group_data.Sample_name, cell_names{jj}),:);

            plot([0:99.5]*0.2, recovery_data.Intensity_norm_both, '.-',...
                'MarkerFaceColor', cols(jj,:),...
                'MarkerEdgeColor', cols(jj,:),...
                'Color', cols(jj,:),...
                'LineWidth',2);

        end
        xlim([-0.5 20]);
        ylim([0 103]);
        set(gca, 'Box', 'off', 'FontSize', 26, 'LineWidth', 3);
        xlabel('time (s)');
        ylabel('% recovery');
        title(groups{ii});
    end

end


%% Make plots of fit for figure

% unconfined
ii = 2;
group_data = frap_results_table(startsWith(frap_results_table.Sample_name, groups{ii}),:);
cell_names = unique(group_data.Sample_name);

jj = 45;
% find intensity at last TP
recovery_data = group_data((strcmp(group_data.Sample_name, cell_names{jj}) & group_data.Time > 1),:);

[f,results] = fit(recovery_data.Time_s_corr, recovery_data.Intensity_norm_both, 'exp2');
rsq = results.rsquare;
figure; tiledlayout(1,2);
nexttile;
plot(f, recovery_data.Time_s_corr, recovery_data.Intensity_norm_both);
set(gca, 'FontSize', 20);
yline(f.a, 'k--', 'DisplayName', 'fast/slow diffusing cutoff');
box off;
%set(gca, 'XTick', [], 'XTickLabels', [], 'YTick', [], 'YTickLabels', [], 'Xcolor', 'none', 'YColor', 'none');
legend('location', 'southeast')
xlabel('Time (s)')
ylabel('% recovery')
subtitle(append('Rsq = ', string(rsq)));


% confined
ii = 1;
group_data = frap_results_table(startsWith(frap_results_table.Sample_name, groups{ii}),:);
cell_names = unique(group_data.Sample_name);

jj = 2;
recovery_data = group_data((strcmp(group_data.Sample_name, cell_names{jj}) & group_data.Time > 1),:);

[f,results] = fit(recovery_data.Time_s_corr, recovery_data.Intensity_norm_both, 'exp2');
rsq = results.rsquare;

nexttile;
plot(f, recovery_data.Time_s_corr, recovery_data.Intensity_norm_both);
set(gca, 'FontSize', 20);
yline(f.a, 'k--', 'DisplayName', 'fast/slow diffusing cutoff');
box off;
%set(gca, 'XTick', [], 'XTickLabels', [], 'YTick', [], 'YTickLabels', [], 'Xcolor', 'none', 'YColor', 'none');
legend('location', 'southeast')
xlabel('Time (s)')
ylabel('% recovery')
subtitle(append('Rsq = ', string(rsq)));






