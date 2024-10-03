%% MTEX TrueEBSD for WC Contiguity calculation
%
% authors: Vivian Tong, National Physical Laboratory, Teddington, UK; 
% Stefan Olovsjö, Seco Tools AB, R&D Materials and Technology, 737 82
% Fagersta, Sweden;
% Contact: vivian.tong@npl.co.uk
%
% Description:  
% Example script to run trueEBSD workflow
% MATLAB R2024a and mtex version forked from feature/grain3d, approx mtex6.0.beta3
%
% Inputs:  
% mtexdata trueEbsdWCCo
%
% Outputs: 
% Published html file containing code and outputs
%
% Version control
% 20241001 - create TrueEBSD example script using data from SECOvisit_020_1_site1

clear; close all; home;

% TrueEBSD version ID 
vId = '20240916 / app version 1.2.1';

%% Add trueEBSD related MATLAB paths 
addpath(genpath(cd));

%% Data Import
% Begin by loading an EBSD map with a list of images we want to use
% together with the EBSD map data. This contains an EBSD map of a WC-Co 
% composite acquired at 20 kV accelerating voltage, and four SEM images of 
% the same sample area within ebsd.opt.trueEbsdImgs: 
%
% # Band contrast (|ebsd.bc|) is used as the image for the EBSD map.
%
% # |fsdB3| is a colour image from the three FSD detectors mounted at
% the bottom of the EBSD camera, and the EBSD camera is retracted by 20 mm
% relative to the EBSD map acquisitiion position;
%
% # |fsdT3| is a greyscale image from the same beam scan as fsdB3 and 
% the FSD detectors at the top of the EBSD camera;
%
% # |fsdT1| is a greyscale image from the FSD detectors at the top of the
% EBSD camera at the EBSD map acquisitiion position;
%
% # |fsdT10| is a greyscale image from the FSD detectors at the top of the
% EBSD camera, in EBSD map acquisitiion position, and the electron beam 
% accelerating voltage lowered to 10 kV.
%
% * |ebsd.opt.trueEbsdImgs.pixSzImg| is the image pixel size in microns for
% all four images.

tic
mtexdata trueEbsdWCCo

display(ebsd);
display(ebsd.opt.trueEbsdImgs);

%% Set up TrueEBSD job
% @distortedImg imgList{:} is a TrueEBSD class containing information 
% about an image or EBSD map and its distortion types within the 
% TrueEBSD workflow.
%
% job is a @trueEbsd object containing a sequence of @distortedImg
% images.

% Construct distortedImg list and set up trueEBSD job
dataName = 'trueEbsdWCCo';

% Do some simple image denoising
ebsd.opt.trueEbsdImgs.fsdB3 = rescale(imboxfilt(ebsd.opt.trueEbsdImgs.fsdB3,3));
ebsd.opt.trueEbsdImgs.fsdT3 = rescale(imboxfilt(ebsd.opt.trueEbsdImgs.fsdT3,3));
ebsd.opt.trueEbsdImgs.fsdT1 = rescale(imboxfilt(ebsd.opt.trueEbsdImgs.fsdT1,3));
ebsd.opt.trueEbsdImgs.fsdT10 = rescale(imboxfilt(ebsd.opt.trueEbsdImgs.fsdT10,3));

% Construct @distortedImg imgList{:} 
imgList=cell(1,5);
imgList{1} = distortedImg('bc','drift-shift', ebsd, 'mapplottingConvention', ebsd.plottingConvention, 'highContrast',1,'edgePadWidth',3);
imgList{2} = distortedImg(ebsd.opt.trueEbsdImgs.fsdB3,'true', 'dxy', ebsd.opt.trueEbsdImgs.pixSzImg, 'highContrast',1,'edgePadWidth',5);
imgList{3} = distortedImg(ebsd.opt.trueEbsdImgs.fsdT3,'shift', 'dxy', ebsd.opt.trueEbsdImgs.pixSzImg, 'highContrast',1,'edgePadWidth',5);
imgList{4} = distortedImg(ebsd.opt.trueEbsdImgs.fsdT1,'tilt', 'dxy', ebsd.opt.trueEbsdImgs.pixSzImg, 'highContrast',1,'edgePadWidth',5);
imgList{5} = distortedImg(ebsd.opt.trueEbsdImgs.fsdT10,'true', 'dxy', ebsd.opt.trueEbsdImgs.pixSzImg, 'highContrast',1,'edgePadWidth',3);

% @trueEbsd job is a TrueEBSD class.
% The starting data for the TrueEBSD workflow are stored in job.imgList.
job = trueEbsd(imgList{:});

%%%
% Plot as-imported image sequence to check they are all of similar regions
% on the sample, but the image contrasts look quite different.

figure('WindowState', 'maximized'); 
t=tiledlayout('flow','TileSpacing','tight','Padding','tight');
title(t,'TrueEBSD starting image sequence');
for n=1:numel(imgList)
    nexttile; 
    imagesc('XData',imgList{n}.dx.*(1:size(imgList{n}.img,2)),...
        'YData',imgList{n}.dy*(1:size(imgList{n}.img,1)),...
        'CData',imgList{n}.img);
    colormap gray; axis image on ij;
end
linkaxes;

t1  = toc;
disp(['Finished set up trueEBSD job for ' dataName ' in ' num2str(t1,'%.1f') ' seconds']);

%% Resize images to match pixel size and FOV
% The EBSD map and images in job.imgList{:} are of the same sample area but
% have different pixel sizes. Here, we match up the pixel positions of the 
% the image sequence in job.imgList{:}.
%
% Inputs - distorted image sequence job.imgList{:}, target pixel size pixSzIn
%
% Outputs - distorted image sequence on a common pixel grid job.resizedList{:}

pixSzIn = 0; % target pixel length in microns, or 0 to default to smallest common pixel size
job = pixelSizeMatch(job,pixSzIn);

%%%
% Now job has a new property job.resizedList, which is where the outputs of
% pixelSizeMatch are stored.
display(job);

t1  = toc;
disp(['Finished resize images to match pixel size and FOV for ' dataName ' in ' num2str(t1,'%.1f') ' seconds']);

%% [Optional] Change cross-correlation function (XCF) ROI settings
% TrueEBSD image registration computes the cross-correlation function (XCF)
% between pairs of regions of interest (ROI) in sequential images. The ROI
% size and spacing within each image pair are tunable parameters in
% TrueEBSD.
%
% The pixelSizeMatch function automatically guesses some normally sensible
% parameters for a polycrystal EBSD map, but you can also set custom values.
% You can define one XCF setting in job.resizedList{n}.setXCF{:}
% When selecting ROI size, a good rule of thumb is an ROI at least 4 times
% wider than the measured local image shifts. It also needs to be a power of 2
% for the XCF to work properly. 
%
% You can set a custom XCF for each individual distortion model in the distorted
% image (job.resizedList{n}.distortionModel).
%
% Here, we deliberately misjudge and choose an ROI box that is too small
% for the EBSD map job.resizedList{1}, which is used to correct EBSD map drift
% (linear interpolation between rigid EBSD map rows).

customSetXCF1.ROISize=2^round(log2(32));
customSetXCF1.NumROI=struct;
customSetXCF1.NumROI.x = 40; % good rule of thumb: as many ROI as grains in FOV
customSetXCF1.NumROI.y = round(customSetXCF1.NumROI.x * size(job.resizedList{1}.img,1)/size(job.resizedList{1}.img,2)); % follow image aspect ratio
customSetXCF1.XCFMesh=250; % correlation peak upsampling, default 250
customSetXCF1.xcfImg = 'edge'; %choose whether to correlate edge transforms or images

customSetXCF2 = customSetXCF1;
customSetXCF2.ROISize=2^round(log2(128));

% assign customSetXCF
job.resizedList{1}.setXCF{2} = customSetXCF1;
job.resizedList{3}.setXCF{1} = customSetXCF2;
% or just rewrite individual properties
job.resizedList{1}.setXCF{1}.ROISize = 2^round(log2(64));
job.resizedList{3}.setXCF{1}.xcfImg = 'img';
job.resizedList{4}.setXCF{1}.xcfImg = 'img';
job.resizedList{5}.setXCF{1}.xcfImg = 'img';



%% Calculate local image shifts and fit to a distortion model
%
% These are the image pairs that will be used for cross-correlation. 
figure('WindowState', 'maximized'); 
t=tiledlayout('flow','TileSpacing','tight','Padding','tight');
title(t,'TrueEBSD image sequence for cross-correlation');
for n=1:numel(job.resizedList)
    nexttile; 
    imagesc('XData',job.resizedList{n}.dx.*(1:size(job.resizedList{n}.img,2)),...
        'YData',job.resizedList{n}.dy*(1:size(job.resizedList{n}.img,1)),...
        'CData',job.resizedList{n}.(job.resizedList{n}.setXCF{1}.xcfImg));
    colormap gray; axis image on ij;
end
linkaxes;

%% Compute image shifts
% Now we compute local image ROI shifts and fit them to distortion models. 
% After each image correction step, the average ROI shifts (X, Y and length
% components) are printed to the command window. 
%
% The 'fitErr' flag means that residual local image shifts are recomputed 
% after image correction but not included in the final result. If this
% number is small (around 1 pixel or less) then most likely the image 
% registration was successful. 
%
% Just now, we intentionally set the ROI box size too small for the EBSD
% map drift correction step. Therefore, the average residual shift is > 2
% pixels long. calcShifts tries to fix this by doubling the ROI size and
% retrying the image registration step. It will try to do this until either
% the ROI are too big to fit into the image, or the residual shifts are < 2
% pixels.
%
% The only exception to this is where the distortion name is 'true', such
% as between images 3 and 2 in this dataset. For this case,TrueEBSD assumes 
% that all the shifts between this image pair are zeros, and ignores the
% residual shifts, even if they are greater than 2 pixels.

job = calcShifts(job,'fitErr');


%%%
% Now job has a new property job.shifts, which is where the outputs of
% calcShifts are stored.

display(job);

t1  = toc;
disp(['Finished calculate image shifts and fit distortion models for ' dataName ' in ' num2str(t1,'%.1f') ' seconds']);

%% Undistort images
% This applies the image shifts between each image pair in job.shifts to
% the data in job.resizedList, and outputs a new property
% job.undistortedList which contains aligned image data. Now all pixels in
% this image sequence can be directly overlaid.

job = undistort(job);

%%% Plot images after distortion correction

figure('WindowState', 'maximized'); 
t=tiledlayout('flow','TileSpacing','tight','Padding','tight');
title(t,'TrueEBSD image sequence after alignment');
for n=1:numel(job.undistortedList)
    nexttile; 
    imagesc('XData',job.undistortedList{n}.dx.*(1:size(job.undistortedList{n}.img,2)),...
        'YData',job.undistortedList{n}.dy*(1:size(job.undistortedList{n}.img,1)),...
        'CData',job.undistortedList{n}.img);
    colormap gray; axis image on ij;
end
linkaxes;

t1  = toc;
disp(['Finished remove image distortions for ' dataName ' in ' num2str(t1,'%.1f') ' seconds']);


%%% Plot data as MTEX EBSD maps
% We can also plot all images as MTEX EBSD maps. This is a good check to
% make sure images are not indexed 'upside down' relative to the EBSD map.
% Since images are usually stored and read by MATLAB using the 'axis ij'
% convention, but EBSD maps can have other kinds of plotting convention defined
% in ebsd.plottingConvention, we need the ij2EbsdSquare helper function to
% rotate the image data into the ebsd map plottingConvention.

figure;
nextAxis;
plot(job.undistortedList{1}.ebsd('W C'), job.undistortedList{1}.ebsd('W C').orientations, ...
    job.undistortedList{1}.ebsd.plottingConvention,'coordinates','on');
title('Undistorted MTEX EBSD map (WC IPF-out of screen)','Color','k');
for n=1:numel(job.undistortedList)
    nextAxis;
    plot(job.undistortedList{1}.ebsd, ...
        ij2EbsdSquare(job.undistortedList{1}.ebsd,job.undistortedList{n}.img), ...
        job.undistortedList{1}.ebsd.plottingConvention,'coordinates','on');
    mtexColorMap gray;
    title(['Undistorted MTEX image ' num2str(n)],'Color','k');
end


%%% Finish
% This is the end of the TrueEBSD distortion correction workflow. 
% 
% You can save your data here, or do any further data analysis that you
% would on normal MTEX EBSD maps.

t1  = toc;
disp(['Finished TrueEBSD workflow for ' dataName ' in ' num2str(t1,'%.1f') ' seconds.']);

%% Further Analysis
% For this dataset, we want to measure the contiguity of the WC grains in
% this EBSD map. That will be covered in the next example script
% example_WCCo_contiguity.


disp('Exiting program now.');


