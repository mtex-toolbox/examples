%% MTEX TrueEBSD for analysing grain boundary voids in a copper polycrystal
%
% authors: Vivian Tong, National Physical Laboratory, Teddington, UK; 
% Contact: vivian.tong@npl.co.uk
%
% Description:  
% Example script to run trueEBSD workflow
% MATLAB R2024a and mtex version forked from feature/grain3d, approx mtex6.0.beta3
%
% Inputs:  
% 
%
% Outputs: 
% Published html file containing code and outputs
%
% Version control
% 20241001 - create TrueEBSD example script using ID29 data (not
% released).

clear; close all; home;

% TrueEBSD version ID 
vId = '20240916 / app version 1.2.1';

%% Add trueEBSD related MATLAB paths 
addpath(genpath(cd));

%% Data Import
% Begin by loading the EBSD map and associated SEM images from the Oxford
% Instruments *.h5oina file.
%
% The SEM image pixel size is stored in the same *.h5oina file, in the 
% Image header data.
%
% # Band contrast (|ebsd.bc|) is used as the image for the EBSD map.
%
% # |fsd1B| is a colour image from the three FSD detectors mounted at
% the bottom of the EBSD camera, and the EBSD camera is retracted by 40 mm
% relative to the EBSD map acquisition position;
%
% # |bse1| is a greyscale image acquired using the annular backscatter (ABS)
% detector at 10 kV and 0 degrees sample tilt. We can filter this in two ways to enhance 
% different types of image contrast: |bse1a| has gamma compression applied to
% enhance grain boundary contrast, and |bse1b| uses a moving median filter 
% to preserve edge contrast of voids in the denoised image.

tic

ebsd = gridify(rotate(...
    loadEBSD_h5oina("copper29.h5oina"),...
    reflection(xvector),'keepEuler'));
ebsd.plottingConvention = plottingConvention(vector3d.Z,-vector3d.X);

display(ebsd);
display(ebsd.opt.Images);

%%%
% Now we extract a reference image which has good contrast for finding the
% grain boundary voids, and an intermediate image to help with matching to
% the EBSD map features.

fsd1B = rescale(im2double(cat(3,ebsd.opt.Images.Lower_Centre_19, ...
    ebsd.opt.Images.Lower_Left_19, ...
    ebsd.opt.Images.Lower_Right_19))); 
bse1 = rescale(im2double(ebsd.opt.Images.ABSinner_0deg)); 



%%
% Do some image pre-processing to help us with 

fsd1a = imboxfilt(fsd1B,5);
bse1a = imboxfilt(nthroot(bse1,0.1),5);
bse1b = medfilt2(bse1,[3 3],'symmetric');


%% Set up TrueEBSD job
% @distortedImg imgList{n} is a TrueEBSD class containing information 
% about an image or EBSD map and its distortion type relative to the next
% [(n+1)th] image in the sequence.
% The final reference image bse1b has very low grain boundary contrast,
% so we set the 'highContrast' flag to 0.


% Construct distortedImg list and set up trueEBSD job
dataName = 'trueEbsdCopper';

% Construct @distortedImg imgList{:} 
imgList=cell(1,4);
imgList{1} = distortedImg('bc','drift-shift', ebsd, 'mapPlottingConvention', ebsd.plottingConvention, 'highContrast',1,'edgePadWidth',3);
imgList{2} = distortedImg(fsd1a,'tilt', 'dxy', double(ebsd.opt.Images.Header.X_Step), 'highContrast',1,'edgePadWidth',3); % 
imgList{3} = distortedImg(bse1a,'true', 'dxy', double(ebsd.opt.Images.Header.X_Step), 'highContrast',1,'edgePadWidth',3); % 
imgList{4} = distortedImg(bse1b,'true', 'dxy', double(ebsd.opt.Images.Header.X_Step), 'highContrast',0,'edgePadWidth',1); % BSE but with pores only

%%
% job is a |@trueEbsd| object containing a sequence of |@distortedImg|
% images.
%
% The starting data for the TrueEBSD workflow are stored in |job.imgList|.

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
% The EBSD map and images in |job.imgList{:}| are of the same sample area but
% have different pixel sizes. Here, we match up the pixel positions of the 
% the image sequence in |job.imgList{:}|.
%
% Inputs - distorted image sequence |job.imgList{:}|, target pixel size pixSzIn
%
% Outputs - distorted image sequence on a common pixel grid
% |job.resizedList{:}|.

pixSzIn = 0; % target pixel length in microns, or 0 to default to smallest common pixel size
job = pixelSizeMatch(job,pixSzIn);

%%%
% Now |job| has a new property |job.resizedList|, which is where the outputs of
% pixelSizeMatch are stored.
display(job);

t1  = toc;
disp(['Finished resize images to match pixel size and FOV for ' dataName ' in ' num2str(t1,'%.1f') ' seconds']);

%% Calculate local image shifts and fit to a distortion model
% Now we compute local image ROI shifts and fit them to distortion models. 
% After each image correction step, the average ROI shifts (X, Y and length
% components) are printed to the command window. 
%
% Images 4 and 3 are the same SEM image filtered using two methods, so there 
% is no distortion between them -- we had set this up in the input data as 
% |job.imgList{1}.distortionName = 'true'|. In this case, TrueEBSD assumes 
% that all the shifts between this image pair are zero and does not
% compute them from ROI cross-correlation.
%
% The 'fitErr' flag means that residual local image shifts are recomputed 
% after image correction but not included in the final result. If this
% number is small (around 1 pixel or less) then most likely the image 
% registration was successful. 
%
%%%
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

%%% Compute image shifts

job = calcShifts(job,'fitErr');

%%%
% Now job has a new property job.shifts, which is where the outputs of
% calcShifts are stored.

display(job);

t1  = toc;
disp(['Finished calculate image shifts and fit distortion models for ' dataName ' in ' num2str(t1,'%.1f') ' seconds']);

%% Undistort images
% This applies the image shifts between each image pair in |job.shifts| to
% the data in |job.resizedList|, and outputs a new property
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
plot(job.undistortedList{1}.ebsd('indexed'), job.undistortedList{1}.ebsd('indexed').orientations, ...
    job.undistortedList{1}.ebsd.plottingConvention,'coordinates','on');
title('Undistorted MTEX EBSD map (Copper IPF-out of screen)','Color','k');
for n=1:numel(job.undistortedList)
    nextAxis;
    plot(job.undistortedList{1}.ebsd, ...
        ij2EbsdSquare(job.undistortedList{1}.ebsd,job.undistortedList{n}.img), ...
        job.undistortedList{1}.ebsd.plottingConvention,'coordinates','on');
    mtexColorMap gray;
    title(['Undistorted MTEX image ' num2str(n)],'Color','k');
end

%%%
% This is the end of the TrueEBSD distortion correction workflow. 
t1  = toc;
disp(['Finished TrueEBSD workflow for ' dataName ' in ' num2str(t1,'%.1f') ' seconds.']);


%% Futher analysis
% For this dataset, we want to find out if there are preferred grain
% boundary nucleation sites.

disp('Continuing onto grain boundary voids analysis.');

%%%
% First we extract EBSD map and BSE image with voids contrast from the
% TrueEBSD processed data.
ebsd = job.undistortedList{1}.ebsd;
bse = job.undistortedList{4}.img; 

%%%
% We need crop both the EBSD map and BSE image to the largest rectangle that 
% contains the EBSD map information. 
[~,~,~,keepGrid] = FindLargestRectangles(~isnan(job.undistortedList{1}.img));
ebsd = gridify(ebsd(ij2EbsdSquare(ebsd,keepGrid)));
% ebsd forgets its plottingConvention after gridify so we remind it here
ebsd.plottingConvention = job.undistortedList{1}.mapPlottingConvention; 
%repeat for binarised reference image
bse = reshape(bse(keepGrid),size(ebsd));

%% Update EBSD phase assignment according to BSE image

phasesBse = ij2EbsdSquare(ebsd,(bse<0.8)); %voids = 1, copper = 0.

voidPhase = crystalSymmetry('1','mineral','voids','color',str2rgb('DarkBlue'));

ebsd(phasesBse).rotations = rotation('euler',0,0,0);
ebsd(phasesBse).CS = voidPhase;
ebsd = gridify(ebsd);
ebsd.plottingConvention = job.undistortedList{1}.mapPlottingConvention; 

display(ebsd);

%% Calculate grains and boundaries
% Now we calculate grains and grain boundaries in the copper grains only.
% By excluding the 'voids' phase in calcGrains, the grain boundaries will
% be drawn through the voids instead of around them. 

[grains,ebsd('Copper').grainId] = calcGrains(ebsd('Copper'),'angle',10*degree);
%redraw boundaries using new contour algorithm
grains = boundaryContours(grains);
gBs = grains.boundary('Copper','Copper'); %explicitly specify copper-copper boundaries to exclude map border
tPs = grains.triplePoints;

%%%
% plot ebsd map with voids overlaid
figure; newMtexFigure('layout',[2,1]);
nextAxis;
plot(ebsd,ebsd.bc,ebsd.plottingConvention,'micronbar','off'); 
mtexColorMap gray; hold on; 
plot(ebsd('indexed'),ebsd.plottingConvention,'FaceAlpha',0.7);
mtexTitle('Band Contrast and Phases');
nextAxis;
plot(ebsd,ones(size(ebsd)),ebsd.plottingConvention); 
mtexColorMap gray; hold on; %force white background
plot(ebsd('Copper'),ebsd('Copper').orientations,'FaceAlpha',0.5,...
    ebsd.plottingConvention,'micronbar','on'); hold on;
plot(gBs,ebsd.plottingConvention,'linecolor',str2rgb('gray'));
plot(ebsd('voids'),zeros(size(ebsd('voids'))),ebsd.plottingConvention);
mtexTitle('Copper Orientations (IPF-out of screen), Grain Boundaries and Voids');

%% Find nearest g.b. per copper void
% We want to find out whether certain grain boundaries or triple junctions
% are more or less susceptible to void nucleation. 
%
% Every grain boundary segment in |gBs| runs between two neighbour EBSD map
% points stored in |gBs.ebsdId|. Let's construct an EBSD map |gbPosMap| that
% shows all EBSD map points that border a grain boundary segment by
% its index to |gBs|.
%
gbPosMap = zeros(size(ebsd));
[ebsdIdList,ia,ic]=unique(gBs.ebsdId,'stable');
% ia are linear indices to unique values of gBs.ebsdId
% gBIdList are row indices to gBs
gBIdList = repmat([1:length(gBs)]',[1,2]);

% gbPosMap is an ebsd map, nonzero values are indices to gBs
gbPosMap(id2ind(ebsd,ebsdIdList))=gBIdList(ia);

%%%
% Now we can find the grain boundary nearest to each void pixel
%
% 1. Compute distance of all EBSD map positions from nearest grain boundary
[gbDist, gbNearest]=bwdist(gbPosMap); % gbNearest are linear indices to gbPosMap
% 2. Construct an EBSD map |voidsMap| showing the locations of all voids,
% where voids = 1 and not-voids = 0. 
voidsMap = zeros(size(ebsd));
voidsMap(ebsd.phase==ebsd('voids').phase(1))=1;
% 3. Find closest g.b. to this void
% we only find intersections between voids and gbs, because the voids
% cover up the copper-copper boundaries
voidsList=gbPosMap(gbNearest(logical(voidsMap))); %indices to gBs
voidsDist=gbDist(logical(voidsMap)); %distance from void to nearest gB in pixels

% separate out voids on the boundary vs close to the boundary (could be
% not on the boundary, could be a image matching error)

% guess a good threshold based on the TrueEBSD tilt fit residuals (95th
% percentile - arbitrary sensible guess?)
voidsList_threshPix = prctile(sqrt(job.shifts{2}{1}.fitError.ROI.Shift_X_1.^2 + job.shifts{2}{1}.fitError.ROI.Shift_Y_1.^2),95); 
%split voidsList into voids on a boundary and voids near the closest gB
voidsList_on = voidsList(voidsDist<=sqrt(2)); %pixel next to copper void
voidsList_near = voidsList(voidsDist>sqrt(2) & voidsDist<=voidsList_threshPix);
voidsList_notNear = voidsList(voidsDist>voidsList_threshPix);

%sort voidsList to prefer gb points on a void
voidsList= [voidsList_on;voidsList_near;voidsList_notNear];

% find nearest gB misorientation
[voidLocs_gIdPairs,ia2,ic2] = unique(gBs(voidsList).grainId,'rows','stable');

% find g.b. misorientation
mdf_voids = calcDensity(gBs(voidsList).misorientation);
mdf_all = calcDensity(gBs.misorientation);


% plot EBSD map with g.b. annotations
figure; newMtexFigure('figSize','large');
plot(ebsd('Copper'),ebsd('Copper').orientations,'FaceAlpha',0.3,ebsd.plottingConvention); hold on;
plot(gBs,ebsd.plottingConvention,'linecolor',str2rgb('gray'));
plot(ebsd('voids'),zeros(size(ebsd('voids'))),ebsd.plottingConvention); colormap gray; clim([0 1]);
plot(gBs(voidsList_on),ebsd.plottingConvention,'linecolor',str2rgb('DarkRed'),'linewidth',3);
plot(gBs(voidsList_near),ebsd.plottingConvention,'linecolor',str2rgb('Red'),'linewidth',3);
plot(gBs(voidsList_notNear),ebsd.plottingConvention,'linecolor',str2rgb('magenta'),'linewidth',2);



% plot misorientation distributions 
figure; 
newMtexFigure('figSize','tiny','outerplotspacing',30);
plotAngleDistribution(mdf_all,'DisplayName','All GBs');hold on
plotAngleDistribution(mdf_voids,'DisplayName','Void GBs');
plotAngleDistribution(ebsd('Copper').CS,ebsd('Copper').CS,'antipodal','DisplayName','Uniform MDF');
legend('show','Location','northwest'); 
xlabel('Misorientation angle / degrees');
ylabel('Frequency / mrd');


figure; newMtexFigure('layout',[2,2],'figSize','large','outerplotspacing',30,'innerplotspacing',50);
nextAxis(1,2); plotAxisDistribution(mdf_all,'colorRange','equal'); mtexTitle('All GBs');
nextAxis(2,2); plotAxisDistribution(mdf_voids,'colorRange','equal'); mtexTitle('Void GBs');
nextAxis(2,1); plotAxisDistribution(ebsd('Copper').CS,ebsd('Copper').CS,'antipodal','colorRange','equal'); mtexTitle('Uniform MDF');
mtexColorbar;


