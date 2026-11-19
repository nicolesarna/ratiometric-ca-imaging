% RatiometricCaImagingSingleTimepoint
%
% Author:       Nicole Suzanne Sarna
% Affiliation:  Department of Bioengineering, Rice University
%               (Laboratory of Dr. Michael R. King)
% Date Created: April 30th, 2026
% Last Updated: May 20th, 2026
% Version:      1.0.1
%
% Adapted from: TCellAnalysis (N. S. Sarna)
%               https://github.com/nicolesarna/TCellAnalysis
% Method described in:
%   Sarna NS, Desai SH, Kaufman BG, Curry NM, Hanna AM, King MR.
%   Enhanced and sustained T cell activation in response to fluid shear
%   stress. iScience. 2024;27(6):109999. doi:10.1016/j.isci.2024.109999
%
% Please retain this authorship header in any copy or modified version of
% this file, and credit the author in any thesis, presentation, or publication
% reporting results produced with it.
%
% Archived at:  https://doi.org/10.5281/zenodo.24380711
% Source:       https://github.com/nicolesarna/ratiometric-ca-imaging
%
% If you use this software, or data produced with it, please cite:
%   Sarna, N. S. (2026). RatiometricCaImagingSingleTimepoint: MATLAB analysis
%   program for ratiometric Fluo-4/Fura-Red calcium imaging (Version 1.0.1)
%   [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24380711
%
% Use for ratiometric calcium imaging with single-timepoint snaps
% Code expects .czi file with 2 channels: Fluo-4 ch1 and FuraRed ch2 
% Segments cells in each image and assigns a CellID per object
% Output table represents one row per cell
% 
% Code used to process confocal images acquired by Ismael Ortiz, (Advisor: Dr.
% Cynthia A. Reinhart-King at Rice University)

clc; clear

maxObjects = 500000; % pre-allocate rows

% Initialize result arrays
imagenum       = zeros(1, maxObjects);
cellID_all     = zeros(1, maxObjects);
allarea        = zeros(1, maxObjects);
allCirc        = zeros(1, maxObjects);
allEccen       = zeros(1, maxObjects);
allDia         = zeros(1, maxObjects);
allSolid       = zeros(1, maxObjects);
gCTCFall       = zeros(1, maxObjects);
rCTCFall       = zeros(1, maxObjects);
gBCmeanall     = zeros(1, maxObjects);
rBCmeanall     = zeros(1, maxObjects);
ratioBCmeanall = zeros(1, maxObjects);
allpxArea      = zeros(1, maxObjects);
allpxDia       = zeros(1, maxObjects);

currentIdx = 1;

%%%%%%%%%%%%%%%%%% GUI boxes %%%%%%%%%%%%%%%%%%%%%%%%
waitfor(msgbox({'Please select single-timepoint .czi files to process.' ...
    'Use shift + left mouse to select more than one file.' ...
    'Files may only be selected from one directory at a time.'}, ...
    'Select Image Help'));

[fname1, pname1] = uigetfile( ...
    {'*.czi'; '*.*'}, 'Select Image(s)', 'MultiSelect', 'on');
fname1 = cellstr(fname1);

saveImages = 0;
choice = questdlg('Save QC segmentation images? (recommended)', 'Save Images?', 'Yes', 'No', 'No');
switch choice
    case 'Yes'; saveImages = 1;
    case 'No';  saveImages = 0;
end

saveData = 0;
choice = questdlg('Export data to Excel?', 'Export Data?', 'Yes', 'No', 'Yes');
switch choice
    case 'Yes'; saveData = 1;
    case 'No';  saveData = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Structuring elements (defined once, shared across all images)
se0   = strel('disk', 3);
se1   = strel('disk', 1);
se2   = strel('disk', 9);
bkgse = strel('disk', 15);

for p = 1:length(fname1)

    filename     = fullfile(pname1, fname1{p});
    [~, name, ~] = fileparts(filename);
    ID           = name;


    r    = bfGetReader(filename);
    numC = r.getSizeC();

    %%%%%%%%%%%%%%%%%% Read OME metadata %%%%%%%%%%%%%%%%%%
    omeMeta = r.getMetadataStore();

    % Spatial calibration — fall back to hardcoded value if metadata unavailable
    pxSize_um  = 0.085166;
    pxArea_um2 = pxSize_um^2;

    try
        % Search file metadata to obtain pixel size
        physX = omeMeta.getPixelsPhysicalSizeX(0);
        if ~isempty(physX)
            pxSize_um  = double(physX.value());
            pxArea_um2 = pxSize_um^2;
        else
            fprintf('  Warning: pixel size not found in metadata; using fallback %.6f um/px.\n', pxSize_um);
        end
    catch ME
        fprintf('  Warning: could not read pixel size (%s); using fallback %.6f um/px.\n', ME.message, pxSize_um);
    end

    %%%%%%%%%%%%%%%%%% Read single-timepoint planes %%%%%%%%%%%%%%%%%%
    Fluo4Idx       = r.getIndex(0, 0, 0) + 1;   % channel 1 = Fluo-4
    FuraRedIdx     = r.getIndex(0, 1, 0) + 1;   % channel 2 = Fura-Red

    FuraRedchannel = bfGetPlane(r, FuraRedIdx);
    Fluo4channel   = bfGetPlane(r, Fluo4Idx);

    r.close();

    %%%%%%%%%%%% Image enhancement and filtering %%%%%%%%%%%
    % NOTE: Enhancement / filtering / thresholding is performed to create
    % binary masks to segment cells appropriately. Fluorescence
    % quantification is only performed on the RAW image channels
    % (Fluo4channel and FuraRedchannel) 
   
    R1 = imadjust(FuraRedchannel);
    G1 = imadjust(Fluo4channel);

    R1filt = medfilt2(R1);
    G1filt = medfilt2(G1);

    visual = imadd(R1filt, G1filt);  % composite for QC overlay

    %%%%%%%%%%%%%%%%%%%%% Cell Mask Image Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Values adjusted May 20th, 2026 using representative image
    bwR1 = imbinarize(R1filt, 0.2);
    bwG1 = imbinarize(G1filt, 0.1);

    allbw = bwR1 | bwG1;

    % Foreground processing
    x1 = imdilate(allbw, se1);
    x2 = imopen(x1, se2);
    x3 = bwareafilt(x2, [1000 100000000]);

    % Background mask
    bkg1 = imdilate(x2, se2);
    bkg  = ~bkg1;

    % Watershed segmentation
    D    = -bwdist(~x3);
    maskD = imextendedmin(D, 7);
    D2   = imimposemin(D, maskD);
    Ld   = watershed(D2);
    segmentedCells = x3;
    segmentedCells(Ld == 0) = 0;

    % Clear border of objects and remove objects <1000 pixels^2
    y2 = imclearborder(segmentedCells);
    y3 = bwareafilt(y2, [1000 100000000]);

    if ~any(y3(:))
        fprintf('  "%s": no cells found, skipping.\n', name);
        continue
    end

    % Background quantification
    BkgInds = find(bkg);
    if isempty(BkgInds)
        fprintf('  "%s": no background region found, skipping.\n', name);
        continue
    end
    
    AvgGreenBkg = mean(double(Fluo4channel(BkgInds)));
    AvgRedBkg   = mean(double(FuraRedchannel(BkgInds)));

    % Regionprops
    objects    = regionprops(y3, 'Area', 'Centroid', 'Circularity', ...
        'Eccentricity', 'EquivDiameter', 'PixelList', 'Solidity');
    numObjects = length(objects);

    gAvg    = zeros(1, numObjects);
    rAvg    = zeros(1, numObjects);
    gIntDen = zeros(1, numObjects);
    rIntDen = zeros(1, numObjects);
    pxArea  = zeros(1, numObjects);
    calArea = zeros(1, numObjects);
    pxDia   = zeros(1, numObjects);
    calDia  = zeros(1, numObjects);
    cellcirc  = zeros(1, numObjects);
    celleccen = zeros(1, numObjects);
    cellSolid = zeros(1, numObjects);

    for n = 1:numObjects
        pixelList = objects(n).PixelList;
        ind = sub2ind(size(Fluo4channel), pixelList(:,2), pixelList(:,1));

        gAvg(n)    = mean(double(Fluo4channel(ind)));
        rAvg(n)    = mean(double(FuraRedchannel(ind)));
        gIntDen(n) = sum(double(Fluo4channel(ind)));
        rIntDen(n) = sum(double(FuraRedchannel(ind)));

        pxArea(n)  = objects(n).Area;
        calArea(n) = objects(n).Area * pxArea_um2;
        pxDia(n)   = objects(n).EquivDiameter;
        calDia(n)  = objects(n).EquivDiameter * pxSize_um;
        cellcirc(n)  = objects(n).Circularity;
        celleccen(n) = objects(n).Eccentricity;
        cellSolid(n) = objects(n).Solidity;
    end

    %%%%%%%%%%%%% Intensity calculations %%%%%%%%%%%%%%%%%%%%
    %   IntDen  = sum of raw pixel values
    %   CTCF    = IntDen - (area_px * mean_background)
    %   BCmean  = CTCF / area_px  =  mean - AvgBkg

    % BCmean (Background Corrected Mean) used for data processing

    GreenCTCF   = gIntDen - (pxArea .* AvgGreenBkg);
    gBCmean     = GreenCTCF ./ pxArea;

    RedCTCF     = rIntDen - (pxArea .* AvgRedBkg);
    rBCmean     = RedCTCF ./ pxArea;

    ratioBCmean = gBCmean ./ rBCmean;

    % Store results
    endIdx = currentIdx + numObjects - 1;
    imagenum(currentIdx:endIdx)        = p;
    cellID_all(currentIdx:endIdx)      = 1:numObjects;
    allpxArea(currentIdx:endIdx)       = pxArea;
    allarea(currentIdx:endIdx)         = calArea;
    allpxDia(currentIdx:endIdx)        = pxDia;
    allDia(currentIdx:endIdx)          = calDia;
    allCirc(currentIdx:endIdx)         = cellcirc;
    allEccen(currentIdx:endIdx)        = celleccen;
    allSolid(currentIdx:endIdx)        = cellSolid;
    gCTCFall(currentIdx:endIdx)        = GreenCTCF;
    rCTCFall(currentIdx:endIdx)        = RedCTCF;
    gBCmeanall(currentIdx:endIdx)      = gBCmean;
    rBCmeanall(currentIdx:endIdx)      = rBCmean;
    ratioBCmeanall(currentIdx:endIdx)  = ratioBCmean;
    currentIdx = currentIdx + numObjects;

    if saveImages == 1
        %%%% 1. Segmentation overlay (composite + magenta outlines + cell ID labels) 
        bw_perim2 = imdilate(bwperim(y3), se1);
        overlay1  = imoverlay(visual, bw_perim2, 'magenta');

        fg = figure('Visible', 'off');
        imshow(overlay1);
        for k = 1:numObjects
            c = objects(k).Centroid;
            text(c(1), c(2), sprintf('%d', k), ...
                'Color', [1 0.2 0.2], 'FontSize', 14);
        end
        title(sprintf('%s  —  %d cells', ID, numObjects));
        saveas(gcf, fullfile(pname1, sprintf('%s_segmentation.png', ID)), 'png');
        close(fg);

        %%%% 2. FuraRed - Fluo-4 - Composite with 20 um scale bar
        [imgH, imgW] = size(R1filt);
        scaleBarPx   = round(20 / pxSize_um);

        % White padding around the whole figure and between panels (pixels)
        pad = 40;
        %%%%%%%%%%%%%%%%%%%%%% EDIT label font size here
        labelFontSz = 40;
        %%%%%%%%%%%%%%%%%%%%%%% EDIT scale bar font size here
        sbarFontSz  = 30;
        %%%%%%%%%%%%%%%%%%%%%%% EDIT scale bar thickness here
        sbarLineW   = 10;

        % Reserve space above each image for the channel label
        labelAreaH = labelFontSz * 2;   % rough pixel height for label row

        Rnorm = mat2gray(double(R1filt)) * 0.8;   % 20% dimmer
        Gnorm = mat2gray(double(G1filt)) * 1;   % Original - brightness not adjusted

        % Build RGB panels: FuraRed in red, Fluo-4 in green, composite in red+green
        panelRed   = cat(3, Rnorm,             zeros(imgH,imgW), zeros(imgH,imgW));
        panelGreen = cat(3, zeros(imgH,imgW),  Gnorm,            zeros(imgH,imgW));
        panelComp  = cat(3, Rnorm,             Gnorm,            zeros(imgH,imgW));

        % Total figure size: 3 images + 4 vertical pads; label area + image + 2 vertical pads
        figW = imgW*3 + pad*4;
        figH = labelAreaH + imgH + pad*2;

        fm = figure('Visible', 'off', 'Color', 'w', ...
            'Units', 'pixels', 'Position', [100 100 figW figH]);

        panels = {panelRed, panelGreen, panelComp};
        labels = {'FuraRed', 'Fluo-4', 'Composite'};

        for col = 1:3
            % Normalized position of each image axes within the figure
            axL = (pad + (col-1)*(imgW + pad)) / figW;
            axB = pad / figH;
            axW = imgW / figW;
            axH = imgH / figH;

            ax = axes(fm, 'Position', [axL, axB, axW, axH], 'Color', 'k');
            image(ax, panels{col});
            set(ax, 'XLim', [0.5 imgW+0.5], 'YLim', [0.5 imgH+0.5], ...
                'XTick', [], 'YTick', [], 'YDir', 'reverse');
            axis(ax, 'off');
            hold(ax, 'on');

            % Channel label centered above panel in black text (sits in white padding area)
            % Annotation coords are in figure-normalized units
            annotation(fm, 'textbox', [axL, axB+axH, axW, labelAreaH/figH], ...
                'String', labels{col}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', labelFontSz, 'FontWeight', 'bold', ...
                'Color', 'k', 'EdgeColor', 'none', 'BackgroundColor', 'none');

            % 20 um scale bar, bottom-right of each image
            sbMargin = 20;
            barX1 = imgW - sbMargin - scaleBarPx;
            barX2 = imgW - sbMargin;
            barY  = imgH - sbMargin;

            line(ax, [barX1 barX2], [barY barY], 'Color', 'w', 'LineWidth', sbarLineW);
            % text(ax, (barX1+barX2)/2, barY - 6, ['20 ' char(956) 'm'], ...
            %     'Color', 'w', 'FontSize', sbarFontSz, ...
            %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end

        drawnow;
        exportgraphics(fm, fullfile(pname1, sprintf('%s_montage.png', ID)), ...
            'BackgroundColor', 'white', 'ContentType', 'image');
        close(fm);
    end

end % end file loop

% Trim arrays to actual data size
imagenum       = imagenum(1:currentIdx-1);
cellID_all     = cellID_all(1:currentIdx-1);
allpxArea      = allpxArea(1:currentIdx-1);
allarea        = allarea(1:currentIdx-1);
allpxDia       = allpxDia(1:currentIdx-1);
allDia         = allDia(1:currentIdx-1);
allCirc        = allCirc(1:currentIdx-1);
allEccen       = allEccen(1:currentIdx-1);
allSolid       = allSolid(1:currentIdx-1);
gCTCFall       = gCTCFall(1:currentIdx-1);
rCTCFall       = rCTCFall(1:currentIdx-1);
gBCmeanall     = gBCmeanall(1:currentIdx-1);
rBCmeanall     = rBCmeanall(1:currentIdx-1);
ratioBCmeanall = ratioBCmeanall(1:currentIdx-1);

fileNameList = reshape(fname1(imagenum), [], 1);

% Build output table — one row per cell
datamat = horzcat( ...
    imagenum', cellID_all', ...
    allpxArea', allarea', allpxDia', allDia', ...
    allCirc', allEccen', allSolid', ...
    rCTCFall', gCTCFall', rBCmeanall', gBCmeanall', ...
    ratioBCmeanall');

T = array2table(datamat, 'VariableNames', { ...
    'File', 'CellID', ...
    'Area_px', 'Area_um2', 'Diameter_px', 'Diameter_um', ...
    'Circularity', 'Eccentricity', 'Solidity', ...
    'RedCTCF', 'GreenCTCF', 'RedBCmean', 'GreenBCmean', ...
    'RatioBCmean_Fluo4overFuraRed'});

T = [table(fileNameList, 'VariableNames', {'FileName'}), T];
T = sortrows(T, {'File', 'CellID'});

if saveData == 1
    fullFileName = fullfile(pname1, 'SingleTimepointCalciumAnalysis.xls');
    writetable(T, fullFileName);
    fprintf('\nDone! Results saved to:\n  %s\n', fullFileName);
else
    fprintf('\nDone! Data export skipped.\n');
end
