clear all
close all

% required to merge GFP and tdT images into a greyscale time stack (convert RGB to 8bit greyscale) to enable segmentation of nuclei at initial TP

%% Read in files

im_folder = '/Volumes/whitelab/Lab Members/MirandaHunter/Microscopy/LSM880/2023/230623_timelapse_fucci/for_quant/';
cd(im_folder)

merged_ims_names = dir(fullfile(im_folder, '*greyscale.tif'));
merged_ims_names = {merged_ims_names.name};

merged_ims = {};
gfp_ims = {};
tdT_ims = {};


for ii = 1:length(merged_ims_names)

    fprintf(append('Reading in movie ', string(ii), ' of ', string(length(merged_ims_names)), '...\n'))

    merged_ims{ii} = tiffreadVolume(merged_ims_names{ii});
    scene = extractBefore(merged_ims_names{ii}, '_xyCorrected');
    gfp_ims{ii} = tiffreadVolume(append(scene, '_ch1_xyCorrected.tif'));
    tdT_ims{ii} = tiffreadVolume(append(scene, '_ch2_xyCorrected.tif'));
end
fprintf('Done!\n')


time_res = 10; % in mins
n_movies = length(gfp_ims);

[n_px_x, n_px_y, n_tp, ~] = size(gfp_ims{1});


%% Segment nuclei in first TP

nuclei_masks = {};

show_segmentation = 1;

for ii = 1:n_movies

    im_movie = merged_ims{ii};

    first_im = im_movie(:,:,1);
    im_mean = mean(first_im(:));

    % set image mean as background
    im_seg = first_im > (im_mean);
    if show_segmentation; figure; imshow(im_seg); end

    % remove edge objects
    im_seg = imclearborder(im_seg);

    % fill holes
    im_seg = bwareafilt(im_seg, [1000 Inf]);
    %if show_segmentation; figure; imshow(im_seg); end

    im_seg = imfill(im_seg, 'holes');
    %if show_segmentation; figure; imshow(im_seg); end

    SE1 = strel('disk', 10);
    im_seg = imopen(im_seg, SE1);
    %if show_segmentation; figure; imshow(im_seg); end

    % final area filter to remove small junk
    im_seg = bwareafilt(im_seg, [1000 Inf]);
    if show_segmentation; figure; imshow(im_seg); end

    index = length(nuclei_masks) + 1;
    nuclei_masks{index} = im_seg;

end


%% Label nuclei and quantify fluorescence over time

gfp_fluor = [];
tdT_fluor = [];

for jj = 1:n_movies

    im_seg = nuclei_masks{jj};
    im_label = bwconncomp(im_seg,8);

    n_cells = im_label.NumObjects;

    for ii = 1:n_cells

        nuc_im = false(n_px_x, n_px_y); 
        nuc_im(im_label.PixelIdxList{ii}) = true;
        %figure; imshow(nuc_im)

        gfp_ims_nuc = uint16(double(gfp_ims{jj}) .* double(nuc_im));
        tdT_ims_nuc = uint16(double(tdT_ims{jj}) .* double(nuc_im));

        [~, n_col] = size(gfp_fluor);
        index = n_col+1;

        % calculate mean of each TP
        gfp_fluor(1:n_tp,index) = mean(mean(gfp_ims_nuc,1),2);
        tdT_fluor(1:n_tp,index) = mean(mean(tdT_ims_nuc,1),2);

    end
end

% Normalize to first TP
gfp_fluor_norm = gfp_fluor ./ gfp_fluor(1,:);
tdT_fluor_norm = tdT_fluor ./ tdT_fluor(1,:);


%% Plot

[n_tp, n_cells] = size(gfp_fluor);

tp_plot = [0:(n_tp-1)] ./ (60/time_res);

cols_g = viridis(n_cells);
cols_r = magma(n_cells);


% % plot per cell
% for ii = 1:n_cells
%     figure;
%     hold on
%     plot(1:n_tp, gfp_fluor(:,ii), 'Color', 'green')
%     plot(1:n_tp, tdT_fluor(:,ii), 'Color', 'magenta')
% end


% average
[m_gfp, s_gfp, e_gfp] = rmsnan(gfp_fluor');
[m_tdT, s_tdT, e_tdT] = rmsnan(tdT_fluor');

% set colors
lineProps2.col{1} = [160 64 152] ./256;
lineProps1.col{1} = [17 121 41] ./256;

% plot normalized
figure; hold on;
for ii = 1:n_cells
    plot(1:n_tp, gfp_fluor_norm(:,ii), 'Color', cols_g(ii,:))
end

figure; hold on;
for ii = 1:n_cells
    plot(1:n_tp, tdT_fluor_norm(:,ii), 'Color', cols_r(ii,:))
end

% % plot per cell
% for ii = 1:n_cells
%     figure;
%     hold on
%     plot(1:n_tp, gfp_fluor_norm(:,ii), 'Color', 'green')
%     plot(1:n_tp, tdT_fluor_norm(:,ii), 'Color', 'magenta')
% end


% average
[m_gfp_norm, s_gfp_norm, e_gfp_norm] = rmsnan(gfp_fluor_norm');
[m_tdT_norm, s_tdT_norm, e_tdT_norm] = rmsnan(tdT_fluor_norm');

figure;
hold on
mseb(tp_plot, m_tdT_norm, e_tdT_norm, lineProps2, 0);
mseb(tp_plot, m_gfp_norm, e_gfp_norm, lineProps1, 0);
xlim([-0.25 16.25])
set(gca, 'linewidth', 2, 'FontSize', 24);
xlabel({'confined time (hours)'});
ylabel({'normalized intensity (a.u)'})
legend('G1', 'G2M/S')





