clear all
close all

im_folder = ('/Volumes/whitelab/Lab Members/MirandaHunter/Microscopy/V16/2022/hmgb2_crispr/');

% check for connection to server
try
    cd(im_folder);
catch
    error('Not connected to labshare server.')
end

import_czis = 1;
import_mats = 0;

flip_ims = 1; % set to 1 if intestinal region is in bottom right (i.e head pointing towards the right of the image)

day_of_teaz = 'Dec 14'

%% Import all CZIs across all conditions and weeks for a particular set (which_set)

% Import CZIs
if import_czis
    
    EV_BF = {};
    EV_GFP = {};
    EV_tdT = {};
    EV_weeks = {};
    EV_pxsize = [];
    
    HMG_BF = {};
    HMG_GFP = {};
    HMG_tdT = {};
    HMG_weeks = {};
    HMG_pxsize = [];
    
    % each week is stored as a folder.
    weeks = dir;
    weeks = weeks([weeks(:).isdir]);
    weeks = {weeks.name};
    weeks = sort(weeks(~contains(weeks, '.'))); % make sure they are in ascending order, otherwise it messes up how the images are stored!
    
    for jj = 1:length(weeks)
        
        week = weeks{jj};
        cd(append(im_folder, week))
        fprintf('Importing image CZIs...\n')
        folder_info = dir;
        folder_names = {folder_info.name};
        folder_names = folder_names(~contains(folder_names, '.'));
        
        groups = folder_names;
        
        % find all the folders for that set of images
        for ii = 1:length(folder_names)
            
            cd(append(im_folder, week))
            group = folder_names{ii};
            
            cd(group)
            ims_folder = dir;
            im_names = {ims_folder.name};
            im_names = im_names(contains(im_names, '.czi'));
            
            count = 0;
            
            % check for duplicated images in the folder. Since i took 2 images per fish, there should be an even number of images in the folder. Throw an error otherwise.
            if (rem(length(im_names), 2))
                message = append('Odd number of images in ', group, ' ', week, '! Check for duplicate images.');
                error(message)
            end
            
            % import images
            for mm = 1:length(im_names)
                
                % odd number indexed images are 3.5X magnification
                % even number index images are 10X magnification
                % can probably just use 10X images for now...
                
                if (rem(mm,2))
                    continue
                else
                    count = count+1;
                    %fprintf('Importing image %d/%d...\n',  count, length(im_names)/2)
                    %im_path = fullfile(folder, group, im_names{mm});
                    data = bfopen(im_names{mm});
                    ims_all = data{1,1};
                    
                    BF_im = ims_all{1,1};
                    GFP_im = ims_all{4,1};
                    tdT_im = ims_all{9,1};
                    
                    % import metadata
                    omeMeta = data{1,4};
                    pixelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value();
                    pixelSizeX_num = pixelSizeX.doubleValue();
                    pixelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value();
                    pixelSizeY_num = pixelSizeY.doubleValue();
                    pixelSize_um2 = pixelSizeX_num * pixelSizeY_num;
                    
                    if flip_ims
                        % need to flip images to match previous orientation
                        BF_im = flipdim(BF_im,2);
                        GFP_im = flipdim(GFP_im,2);
                        tdT_im = flipdim(tdT_im,2);
                        %figure; subplot(1,3,1); imshow(imadjust(BF_im)); subplot(1,3,2); imshow(imadjust(GFP_im)); subplot(1,3,3); imshow(imadjust(tdT_im));
                    end
                    
                    switch group
                        case 'NT'
                            EV_BF{jj,count} = BF_im; % brightfield
                            EV_GFP{jj,count} = GFP_im; % GFP3
                            EV_tdT{jj,count} = tdT_im; % tdT4
                            EV_pxsize(jj,count) = pixelSize_um2;
                        case 'hmgb2_gRNA'
                            HMG_BF{jj,count} = BF_im; % brightfield
                            HMG_GFP{jj,count} = GFP_im; % GFP3
                            HMG_tdT{jj,count} = tdT_im; % tdT4
                            HMG_pxsize(jj,count) = pixelSize_um2;
                    end
                end
            end
            cd(append(im_folder, week, '/', group))
        end
    end
    cd(im_folder)
end

%     cd(im_folder)
%     save_folder = char('TEAZ_quants_replicate2');
%     mkdir(save_folder)
%     cd(save_folder)
%
%     % save .mat files
%     save('sgMH3_BF.mat', 'sgMH3_BF');
%     save('sgMH3_GFP.mat', 'sgMH3_GFP');
%     save('sgMH3_tdT.mat', 'sgMH3_tdT');
%     save('sgMH3_weeks.mat', 'sgMH3_weeks');
%
%     save('sgNT_BF.mat', 'sgNT_BF');
%     save('sgNT_GFP.mat', 'sgNT_GFP');
%     save('sgNT_tdT.mat', 'sgNT_tdT');
%     save('sgNT_weeks.mat', 'sgNT_weeks');
%
%     % also save a backup
%     backup_folder = char(append('image_mats_', date));
%     mkdir(backup_folder)
%     cd(backup_folder)
%     save('sgMH3_BF.mat', 'sgMH3_BF');
%     save('sgMH3_GFP.mat', 'sgMH3_GFP');
%     save('sgMH3_tdT.mat', 'sgMH3_tdT');
%     save('sgMH3_weeks.mat', 'sgMH3_weeks');
%
%     save('sgNT_BF.mat', 'sgNT_BF');
%     save('sgNT_GFP.mat', 'sgNT_GFP');
%     save('sgNT_tdT.mat', 'sgNT_tdT');
%     save('sgNT_weeks.mat', 'sgNT_weeks');

cd(im_folder);

clearvars -except EV_BF EV_GFP EV_tdT HMG_BF HMG_GFP HMG_tdT EV_pxsize HMG_pxsize groups weeks im_folder


%% Threshold GFP signal

EV_GFP_thresh = {};
HMG_GFP_thresh = {};


for ii = 1:length(groups)
    
    group = groups{ii};
    
    switch group
        case 'NT'
            ims_all = EV_GFP;
        case 'hmgb2_gRNA'
            ims_all = HMG_GFP;
    end
    
    [n_weeks,n_ims] = size(ims_all);
    
    for jj = 1:n_weeks
        
        for kk = 1:n_ims
            
            im = ims_all{jj,kk};
            
            %figure; imagesc(im); axis off;
            
            im_mean = mean(im(:));
            im_std = std(double(im(:)));
            im_bsub = im - im_mean;
            
            im_thresh = im > (im_mean+(im_std*0.5));
            %figure; imshow(im_thresh) % basic thresholded image
            
            % dilate the tumor mass a ton and use that as a mask to get rid of random stuff closer to the edges
            %im_thresh = imclearborder(im_thresh);
            SE1 = strel('disk', 30); 
            mask = imclose(im_thresh, SE1);
            %figure; imshow(mask)
            mask = imfill(mask, 'holes');
            %figure; imshow(mask)
            mask = imclearborder(mask);
            %figure; imshow(mask)
            
            %                 mask = mask .* im_thresh;
            %                 figure; imshow(mask)
            %
            %                 SE2 = strel('disk', 30);
            %                 mask = imdilate(mask, SE2);
            %                 mask = imfill(mask, 'holes');
            %                 figure; imshow(mask)
            %                 % %mask = bwareafilt(mask, [10000 1000000]);
            %                 %figure; imshow(mask);
            
            im_thresh = im_thresh .* mask;
            %figure; imshow(im_thresh); title('thresholded tumor without edge objects')
            
            switch group
                case 'NT'
                    EV_GFP_thresh{jj,kk} = im_thresh;
                case 'hmgb2_gRNA'
                    HMG_GFP_thresh{jj,kk} = im_thresh;
            end
        end
    end
end



%% Brightfield segmentation (for pigmented tumors)

fprintf('Thresholding brightfield images...\n');

EV_BF_thresh = {};
HMG_BF_thresh = {};

for kk = 1:length(groups)
    group = groups{kk};
    
    switch group
        case 'NT'
            BF_ims = EV_BF;
        case 'hmgb2_gRNA'
            BF_ims = HMG_BF;
    end
    
    [n_weeks,n_ims] = size(BF_ims);
    
    for jj = 1:n_weeks
        week = weeks{jj};
        
        for ii = 1:length(BF_ims)
            
            %fprintf(append('Thresholding ', which_set, ' ', group, ' brightfield image ', string(ii), '...\n'));
            
            BF_im = BF_ims{jj,ii};
            if isempty(BF_im)
                continue
            end
            
            % automatically determine where the peaks are in the histogram
            [counts,edges] = histcounts(BF_im);
            [~,locs] = findpeaks(counts, 'NPeaks', 3); % find the 3 tallest peaks in the histogram
            binwidth = edges(2) - edges(1);
            im_peaks = edges(locs) - (binwidth/2);
            
            % set threshold
            im_peaks = sort(im_peaks, 'ascend');
            %thresh = sum(im_peaks(1:2))/2;
            thresh = im_peaks(1) + 30; 
            im_thresh = BF_im < thresh;
            im_thresh = imclearborder(im_thresh);
            SE1 = strel('disk', 2); % originally 2
            im_thresh = imopen(im_thresh, SE1);
            %figure; subplot(1,2,1); imagesc(BF_im); axis off; subplot(1,2,2); imshow(im_thresh);
            
            switch group
                case 'NT'
                    EV_BF_thresh{jj,ii} = im_thresh;
                case 'hmgb2_gRNA'
                    HMG_BF_thresh{jj,ii} = im_thresh;
                    
            end
        end
    end
end

%% GFP and tdTomato segmentation

% fprintf('Thresholding GFP and tdTomato images...\n')
% 
% show_ims = 1;
% 
% sgNT_mask_int = {};
% sgMH3_mask_int = {};
% 
% sgNT_GFP_thresh = {};
% sgMH3_GFP_thresh = {};
% 
% sgNT_tdT_thresh = {};
% sgMH3_tdT_thresh = {};
% 
% for kk = 1:length(groups)
%     
%     group = groups{kk};
%     
%     switch group
%         case 'sgNT'
%             GFP_ims = sgNT_GFP;
%             tdT_ims = sgNT_tdT;
%             weeks = sgNT_weeks;
%         case 'sgMH3'
%             GFP_ims = sgMH3_GFP;
%             tdT_ims = sgMH3_tdT;
%             weeks = sgMH3_weeks;
%     end
%     
%     
%     for jj = 1:length(weeks)
%         
%         for ii = 1:length(GFP_ims)
%             
%             GFP_im = GFP_ims{jj,ii};
%             tdT_im = tdT_ims{jj,ii};
%             
%             if isempty(GFP_im)
%                 continue
%             else
%                 % GFP
%                 GFP_mean = mean(GFP_im(:));
%                 GFP_std = std(double(GFP_im(:)));
%                 GFP_im_bsub = GFP_im - GFP_mean;
%                 
%                 GFP_im_thresh = GFP_im > (GFP_mean + GFP_std);
%                 %figure; imshow(GFP_im_thresh) % basic thresholded image
%                 
%                 % tdT
%                 tdT_mean = mean(tdT_im(:));
%                 tdT_std = std(double(tdT_im(:)));
%                 tdT_im_bsub = tdT_im - tdT_mean;
%                 
%                 tdT_im_thresh = tdT_im > (tdT_mean+tdT_std*1.5);
%                 %figure; imshow(tdT_im_thresh) % basic thresholded image
%                 
%                 % dilate the tumor mass a ton and use that as a mask to get rid of random stuff closer to the edges
%                 % GFP
%                 SE1 = strel('disk', 2);
%                 GFP_mask = imdilate(GFP_im_thresh, SE1);
%                 GFP_mask = imfill(GFP_mask, 'holes');
%                 %figure; imshow(GFP_mask)
%                 
%                 GFP_mask_intestinal = bwareafilt(GFP_mask, [5000 1000000]); % this mask is the autofluorescent intestinal region
%                 %figure; imshow(GFP_mask_intestinal)
%                 
%                 % check for presence of an object in the bottom left corner of the thresholded image (some images don't have it for some reason)
%                 [rows, columns] = size(GFP_im);
%                 row2 = floor(rows/2);
%                 row3 = row2 + 1;
%                 lowerleft_im = imcrop(GFP_mask_intestinal, [1 row3 100 row2]);
%                 if nnz(lowerleft_im) > 1 % meaning there is an autofluorescent intestinal region in the bottom left corner
%                     GFP_mask_intestinal = GFP_mask_intestinal .* ~imclearborder(GFP_mask_intestinal);
%                     SE2 = strel('disk', 25);
%                     GFP_mask_intestinal = imdilate(GFP_mask_intestinal, SE2);
%                 else
%                     GFP_mask_intestinal = [];
%                 end
%                 
%                 % tdT
%                 SE2 = strel('disk', 2);
%                 tdT_mask = imdilate(tdT_im_thresh, SE1);
%                 tdT_mask = imfill(tdT_mask, 'holes');
%                 tdT_mask_intestinal = bwareafilt(tdT_mask, [10000 1000000]); % this mask is the autofluorescent intestinal region
%                 SE3 = strel('disk', 30);
%                 tdT_mask_intestinal = imdilate(tdT_mask_intestinal, SE3);
%                 % find the edge object
%                 tdT_mask_intestinal = tdT_mask_intestinal .* ~imclearborder(tdT_mask_intestinal);
%                 
%                 % combine both
%                 if ~isempty(GFP_mask_intestinal)
%                     mask_intestinal = GFP_mask_intestinal + tdT_mask_intestinal;
%                 else
%                     mask_intestinal = tdT_mask_intestinal;
%                 end
%                 
%                 % remove intestinal signal
%                 GFP_im_thresh = GFP_im_thresh .* ~mask_intestinal;
%                 tdT_im_thresh = tdT_im_thresh .* ~mask_intestinal;
%                 
%                 % remove noise from tdT image
%                 SE3 = strel('disk', 2);
%                 tdT_im_thresh = imopen(tdT_im_thresh, SE3);
%                 
%                 if show_ims
%                     figure; subplot(1,3,1); imshow(mask_intestinal); subplot(1,3,2); imshow(GFP_im_thresh); subplot(1,3,3); imshow(tdT_im_thresh);
%                     figure; subplot(1,2,1); imagesc(GFP_im); axis off; subplot(1,2,2); imagesc(tdT_im); axis off;
%                 end
%                 
%                 switch group
%                     case 'sgNT'
%                         sgNT_GFP_thresh{jj,ii} = GFP_im_thresh;
%                         sgNT_tdT_thresh{jj,ii} = tdT_im_thresh;
%                         sgNT_mask_int{jj,ii} = mask_intestinal;
%                     case 'sgMH3'
%                         sgMH3_GFP_thresh{jj,ii} = GFP_im_thresh;
%                         sgMH3_tdT_thresh{jj,ii} = tdT_im_thresh;
%                         sgMH3_mask_int{jj,ii} = mask_intestinal;
%                 end
%             end
%         end
%     end
% end


%% Display BF and GFP segmentation
% [n_weeks, n_im] = size(HMG_GFP);
% %for jj = 1:n_weeks
% jj = 3;
%     for ii = 1:n_im
%         figure;
%         subplot(2,2,1); imagesc(HMG_GFP{jj,ii}); axis off
%         subplot(2,2,2); imshow(HMG_GFP_thresh{jj,ii});
%         subplot(2,2,3); imagesc(HMG_BF{jj,ii}); axis off
%         subplot(2,2,4); imshow(HMG_BF_thresh{jj,ii});
%     end
%end


%% Combine BF and GFP masks.

fprintf('Merging GFP and brightfield segmentation...\n')

EV_mask_tumor = cellfun(@(GFPmask, BFmask) GFPmask + BFmask, EV_GFP_thresh, EV_BF_thresh, 'UniformOutput', false);
HMG_mask_tumor = cellfun(@(GFPmask, BFmask) GFPmask + BFmask, HMG_GFP_thresh, HMG_BF_thresh, 'UniformOutput', false);


% plot all HMG images
[n_weeks,n_im] = size(HMG_mask_tumor);
for jj = 1:n_weeks
    for ii = 1:n_im
        if isempty(HMG_GFP{jj,ii})
            continue
        else 
            figure;
            subplot(2,2,1); imshow(imadjust(HMG_GFP{jj,ii}));
            subplot(2,2,2); imshow(imadjust(HMG_tdT{jj,ii}));
            subplot(2,2,3); imagesc(imadjust(HMG_BF{jj,ii})); axis off;
            subplot(2,2,4); imshow(HMG_mask_tumor{jj,ii});
        end
    end
end

% plot all NT images
[n_weeks,n_im] = size(EV_mask_tumor);
for jj = 1:n_weeks
    for ii = 1:n_im
        if isempty(EV_GFP{jj,ii})
            continue
        else
            figure;
            subplot(2,2,1); imshow(imadjust(EV_GFP{jj,ii}));
            subplot(2,2,2); imshow(imadjust(EV_tdT{jj,ii}));
            subplot(2,2,3); imagesc(imadjust(EV_BF{jj,ii})); axis off;
            subplot(2,2,4); imshow(EV_mask_tumor{jj,ii});
        end
    end
end



%% Quantify tumor area

fprintf('Calculating tumor area...\n');

EV_tumor_area = cellfun(@nnz, EV_mask_tumor);
HMG_tumor_area = cellfun(@nnz, HMG_mask_tumor);

% convert to um2
EV_tumor_area_um = EV_tumor_area ./ EV_pxsize;
HMG_tumor_area_um = HMG_tumor_area ./ HMG_pxsize;

% Manual correction
% HMG_tumor_area_um(1,2:end) = 0;
% HMG_tumor_area_um(2,2:3) = 0;
% HMG_tumor_area_um(2,5:end) = 0;


% % For 2 weeks:
% %% manual correction for HMG: set all the areas except the first image to 0. the first image is the only one that has a tumor, the rest are the thresholding picking up background
% HMG_tumor_area_um(2:end) = 0;

% For 4 weeks:
% manual correction for HMG: set all the areas except 1 and 4 to 0.
% HMG_tumor_area_um(2:3) = 0;
% HMG_tumor_area_um(5:end) = 0;

[EV_m, EV_s, EV_e] = rmsnan(EV_tumor_area_um');
[HMG_m, HMG_s, HMG_e] = rmsnan(HMG_tumor_area_um');




%% Plotting

%groups_plot = {'sgMH1', 'sgMH2', 'sgMH3', 'sgMH4', 'sgMH5', 'sgMH6', 'sgMH7', 'sgMH9', 'sgNT'};
groups_plot = {'NT', 'hmgb2a/hmgb2b gRNA'};

% area over time
weeks_num = [3,4,6,8,10];
cols = turbo(2);
pt_size = 200;
line_width = 4;
set(0, 'DefaultFigureRenderer', 'painters');

figure;
subplot(1,2,1)
hold on
lineProps1.col{1} = cols(2,:); % sgMH3
lineProps2.col{1} = cols(1,:); % sgNT

lineProps1.width = line_width;
lineProps2.width = line_width;

mseb(weeks_num, EV_m, EV_e, lineProps2, 0);
ax_2 = scatter(weeks_num, EV_m, pt_size, cols(1,:), 'filled');
mseb(weeks_num, HMG_m, HMG_e, lineProps1, 0); % for transparent error bars: last variable = 1
ax_1 = scatter(weeks_num, HMG_m, pt_size, cols(2,:), 'filled');


handles_for_legend = [ax_1(1) ax_2(1)];

box off
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 2;
xlim([weeks_num(1)-0.25 weeks_num(end)+0.25]);
xlabel('weeks post-TEAZ');
ylabel('tumor area (um2)');
set(gca, 'XTick', [0:1:40]);
legend(handles_for_legend, flip(groups_plot), 'location', 'northwest', 'box', 'off')
%title('mean area')

% compare areas under curve
% % NOTE: this doesn't work bc we don't know which fish corresponds to what over time
% EV_area_data = [EV_tumor_area_um(1,:); EV_tumor_area_um(last_week,:)]';
% EV_aocs = trapz(test, EV_area_data);
% 
% HMG_area_data = [areas_2weeks(1:10,2), areas_4weeks(1:10,2)];
% HMG_aocs = trapz(weeks_num, HMG_area_data');
% 
% [~,p_aoc] = ttest2(EV_aocs, HMG_aocs);


% plot at final week mark
plot_data = nan(20,2);
n_ims_HMG = nnz(~isnan(HMG_tumor_area_um(length(weeks),:)));
n_ims_EV = nnz(~isnan(EV_tumor_area_um(length(weeks),:)));

plot_data(1:n_ims_EV,1) = EV_tumor_area_um(length(weeks), 1:n_ims_EV);
plot_data(1:n_ims_HMG,2) = HMG_tumor_area_um(length(weeks), 1:n_ims_HMG);


%figure;
subplot(1,2,2)
hold on
boxplot(plot_data, 'labels', groups_plot, 'Colors', 'k')
lines = findobj(gcf, 'type', 'line');
set(lines, 'linewidth', 2)
plotSpread(plot_data, 'distributionColors', turbo(length(groups_plot)), 'xNames', groups_plot, 'spreadWidth', 1)
box off
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 2;
ylabel('tumor area (um2)')
xlim([0.5 2.5])
%ylim([0 750])

[~,p] = ttest2(plot_data(:,1), plot_data(:,2))



