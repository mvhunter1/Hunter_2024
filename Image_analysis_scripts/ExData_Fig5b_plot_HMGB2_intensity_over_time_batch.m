clear all
close all


%% Import trackmate files

trackmate_res_folder = '/Volumes/whitelab/Lab Members/MirandaHunter/Microscopy/LSM880/2023/230412_timelapse_129_198/batch_movies/';
cd(trackmate_res_folder);

trackmate_files = dir;
trackmate_file_names = {trackmate_files.name};

trackmate_results = trackmate_file_names(contains(trackmate_file_names, 'spots.csv'));

for ii = 1:length(trackmate_results)

    trackmate_csv = readtable(trackmate_results{ii});

    % remove any nans
    trackmate_csv = trackmate_csv(~isnan(trackmate_csv.ID),:);
    
    % need to append TRACK_ID to make unique across movies
    trackmate_csv.TRACK_ID = append(string(trackmate_csv.TRACK_ID), "_", string(ii));

    if ii == 1
        trackmate_results_table = trackmate_csv;
    else
        trackmate_results_table = [trackmate_results_table; trackmate_csv];
    end   
end


%% Remove anything that doesn't start at first TP

all_tracks = unique(trackmate_results_table.TRACK_ID);

% automatically determine the length of the longest tracks
track_lengths = [];
for jj = 1:length(all_tracks)

   [n_row,~] = size(trackmate_results_table(strcmp(trackmate_results_table.TRACK_ID, all_tracks{jj}),:));
   track_lengths(jj,1) = n_row;

end
n_tp_max = max(track_lengths);

% now filter out anything that doesn't have the max # of TPs
tracks_to_keep = all_tracks(track_lengths == n_tp_max);

trackmate_results_table_filt = trackmate_results_table(ismember(trackmate_results_table.TRACK_ID, tracks_to_keep),:);


%% Normalize to first TP


for jj = 1:length(tracks_to_keep)

    % find TP 0 data
    tp_0_intensity = trackmate_results_table_filt.MEAN_INTENSITY_CH1(trackmate_results_table_filt.POSITION_T == 0 & strcmp(trackmate_results_table_filt.TRACK_ID, tracks_to_keep{jj}),:);
    
    track_data = trackmate_results_table_filt(strcmp(trackmate_results_table_filt.TRACK_ID, tracks_to_keep{jj}),:);
    track_data.NORM_INTENSITY_CH1 = track_data.MEAN_INTENSITY_CH1 ./ tp_0_intensity;

    if jj == 1
        trackmate_results_table_filt_norm = track_data;
    else
        trackmate_results_table_filt_norm = [trackmate_results_table_filt_norm; track_data];
    end

end

% convert to hours
trackmate_results_table_filt_norm.POSITION_T_HRS = trackmate_results_table_filt_norm.POSITION_T ./60 ./60;

%% Fit a line to each curve and determine fit

cols = viridis(length(tracks_to_keep));

fit_error_threshold = 0.2; % remove any curve where the difference between the data and the line of best fit is bigger than this

show_fit_plots = 0;
slopes = [];

for ii = 1:length(tracks_to_keep)

    cell = tracks_to_keep{ii};
    plot_data = trackmate_results_table_filt_norm(strcmp(trackmate_results_table_filt_norm.TRACK_ID, cell),:);
    plot_data = sortrows(plot_data, 'POSITION_T_HRS', 'ascend');

    if show_fit_plots
        figure;
        hold on
        plot(plot_data.POSITION_T_HRS, plot_data.NORM_INTENSITY_CH1, '.-',...
            'MarkerFaceColor', 'k',...
            'MarkerEdgeColor', 'k',...
            'Color', 'k',...
            'LineWidth', 3);
    end

    % fit line
    [c,S] = polyfit(plot_data.POSITION_T_HRS, plot_data.NORM_INTENSITY_CH1,1);
    [y_est, delta] = polyval(c, plot_data.POSITION_T_HRS, S);

    % evaluate fit
    T = table(plot_data.POSITION_T_HRS, plot_data.NORM_INTENSITY_CH1, y_est, plot_data.NORM_INTENSITY_CH1-y_est, 'VariableNames', {'X', 'Y', 'Fit', 'FitError'});

    if (sum(abs(T.FitError) > fit_error_threshold) > 4) % if the difference between the data and the line of best fit is off by more than 0.1, discard that sample
        tracks_to_keep{ii} = []; % remove that track from analysis
        slopes(ii,1) = nan;
        continue
    else
        slopes(ii,1) = c(1);
        % plot
        if show_fit_plots
            plot(plot_data.POSITION_T_HRS, y_est, 'r--', 'LineWidth', 3);
        end
    end

end

tracks_to_keep = rmmissing(tracks_to_keep);
slopes = rmmissing(slopes);


%% Plot

figure;
hold on
for ii = 1:length(tracks_to_keep)
    cell = tracks_to_keep{ii};
    plot_data = trackmate_results_table_filt_norm(strcmp(trackmate_results_table_filt_norm.TRACK_ID, cell),:);
    plot_data = sortrows(plot_data, 'POSITION_T_HRS', 'ascend');
    %scatter(plot_data.POSITION_T_HRS, plot_data.NORM_INTENSITY_CH1, 50, cols(ii,:), 'filled');
    plot(plot_data.POSITION_T_HRS, plot_data.NORM_INTENSITY_CH1, '.-',...
        'MarkerFaceColor', cols(ii,:),...
        'MarkerEdgeColor', cols(ii,:),...
        'Color', cols(ii,:),...
        'LineWidth',4)
end
xlim([-0.5 17]);
%ylim([0 4]);
set(gca, 'Box', 'off', 'FontSize', 26, 'LineWidth', 3);
xlabel('time (hours)');
ylabel({'normalized HMGB2 intensity'});

% calculate mean intensity at final tp
intensity_final_tp = trackmate_results_table_filt_norm.NORM_INTENSITY_CH1((trackmate_results_table_filt_norm.POSITION_T_HRS == max(trackmate_results_table_filt_norm.POSITION_T_HRS)),:);
[m_end, s_end, e_end] = rmsnan(intensity_final_tp)

% calculate mean slope per min
[m_slope, s_slope, e_slope] = rmsnan(slopes); % slope is per 10 min, convert to per min
m_slope = m_slope *6 % convert to per hour
e_slope = e_slope /10


% export raw data for paper
raw_data = trackmate_results_table_filt_norm(contains(trackmate_results_table_filt_norm.TRACK_ID, tracks_to_keep),:);



