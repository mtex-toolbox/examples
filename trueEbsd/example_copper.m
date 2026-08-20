%% TrueEBSD grain boundary voids in a copper polycrystal
%
% Authors: Vivian Tong. EBSD data from "Void-Microstructure Correlation in
% Thin Film Copper Power Semiconductor Metallization using MTEX", Matthias
% Grabner, Master's Thesis, Graz University of Technology, 2023.
% Contact: vivian.tong@extern.tu-freiberg.de
%
% Thin film copper metallization develops voids as it ages. This script
% first aligns an EBSD map with the SEM images of the same area using
% TrueEBSD, and then asks the question that alignment makes answerable: do
% the voids sit preferentially on grain boundaries or triple junctions, and
% are some boundary types more resistant than others?
%
% The alignment half of the workflow is explained step by step in
% <example_WCCo.html TrueEBSD on a WC-Co composite> — read that one first if
% you are new to TrueEBSD. This page keeps the alignment brief and spends
% its time on the voids analysis.
%
% This script needs the TrueEBSD toolbox, which is *not* part of MTEX and
% is distributed separately under Apache-2.0:
% <https://github.com/vtvivian/mtex-trueebsd>
%
%   addpath(genpath('<path to mtex-trueEbsd>'))
%
% Use the |mtex7-compat| branch with MTEX 7. Beyond MTEX it requires MATLAB
% R2024a or newer and the Image Processing, Curve Fitting, and Statistics
% and Machine Learning toolboxes. Runtime is minutes, not seconds.
%
%% Data Import
% The data set is a single Oxford Instruments |.h5oina| file holding both
% the EBSD map and the SEM images of the same area. It is 164 MB and is not
% shipped with MTEX, so it is downloaded on first use from Zenodo:
% <https://zenodo.org/records/16902083 zenodo.org/records/16902083>.

fName = fullfile(mtexDataPath,'EBSD','copper29.h5oina');

if ~isfile(fName)
  websave(fName,'https://zenodo.org/records/16902083/files/copper29.h5oina');
end

%%
% The map is stored mirrored with respect to its Euler angles, so it is
% reflected back with |'keepEuler'| — the orientations are the ground truth
% here, the coordinates are what needs correcting.

ebsd = gridify(rotate(EBSD.load(fName),reflection(xvector),'keepEuler'));
ebsd.how2plot.east = -xvector;
ebsd.how2plot.outOfScreen = zvector;

display(ebsd)

%%
% The images that ship inside the |.h5oina| container come back in
% |ebsd.opt.electron_image|, one field per detector, plus a |Header| that
% states their pixel size.
%
% # Band contrast (|ebsd.bc|) is the image belonging to the EBSD map.
%
% # |fsd1B| is a colour image from the three FSD detectors mounted at the
% bottom of the EBSD camera, with the camera retracted by 40 mm relative to
% the EBSD map acquisition position.
%
% # |bse1| is a greyscale image from the annular backscatter (ABS) detector
% at 10 kV and 0 degrees sample tilt.

semImgs = ebsd.opt.electron_image;

display(semImgs)

%%
% The backscatter image is filtered two ways, because the two things we
% need from it respond to opposite treatments. |bse1a| gets gamma
% compression, which brings up grain boundary contrast and is what the
% registration needs. |bse1b| gets a moving median, which preserves the
% edges of the voids and is what the analysis needs.

fsd1B = rescale(im2double(cat(3,semImgs.Lower_Centre_19, ...
  semImgs.Lower_Left_19, ...
  semImgs.Lower_Right_19)));
bse1  = rescale(im2double(semImgs.ABSinner_0deg));

fsd1a = imboxfilt(fsd1B,5);
bse1a = imboxfilt(nthroot(bse1,0.1),5);
bse1b = medfilt2(bse1,[3 3],'symmetric');

%% Set up the TrueEBSD job
% A <distortedImg.distortedImg.html distortedImg> holds one image or EBSD
% map, its pixel size, its plotting convention, and the name of the
% distortion separating it from the *next* entry in the sequence. The
% sequence runs from most distorted to ground truth, and the reference
% image carries |'true'|.
%
% Images 3 and 4 are the same backscatter image filtered two ways, so
% nothing separates them and both are |'true'|. The last one has very low
% grain boundary contrast by design, so it is flagged |'highContrast',0|.

dxyImg = double(semImgs.Header.X_Step);

imgList = createArray(4,1,'distortedImg');
imgList(1) = distortedImg('bc','shift-drift', ebsd, 'how2plot', ebsd.how2plot, ...
  'highContrast',1, 'edgePadWidth',3);
imgList(2) = distortedImg(fsd1a,'tilt', 'dxy', dxyImg, 'highContrast',1, 'edgePadWidth',3);
imgList(3) = distortedImg(bse1a,'true', 'dxy', dxyImg, 'highContrast',1, 'edgePadWidth',3);
imgList(4) = distortedImg(bse1b,'true', 'dxy', dxyImg, 'highContrast',0, 'edgePadWidth',1);

job = trueEbsd(imgList)

%%
% Plot the as-imported sequence to check that the maps cover similar
% regions of the sample.

plotImgList(imgList,'TrueEBSD starting image sequence')

%% Align the sequence
% The three alignment steps are the same as in the WC-Co example:
% resample everything onto one common pixel grid, cross-correlate
% regions of interest between successive pairs and fit the named
% distortion model, then apply the accumulated shifts.

pixSzIn = 0; % target pixel length in microns, or 0 for the smallest present
job = pixelSizeMatch(job,pixSzIn);

%%
% These are the images that will actually be cross-correlated — the edge
% transform where |xcfImg| is |'edge'|, the raw values where it is |'img'|.

plotImgList(job.resizedList,'TrueEBSD image sequence for cross-correlation','xcf')

%%
% |'fitErr'| re-measures the residual local shifts after correction without
% including them in the result. A residual around one pixel or less means
% the registration most likely succeeded. Those residuals are kept in
% |job.fitError| and are used further down to decide how close a void has
% to be to a boundary before we call it "on" the boundary.

job = calcShifts(job,'fitErr');

job = undistort(job);

display(job)

plotImgList(job.undistortedList,'TrueEBSD image sequence after alignment')

%%
% Plotting the images back onto the EBSD map is a good check that nothing
% is indexed upside down relative to the map. Images are stored and read by
% MATLAB in the |axis ij| convention, whereas an EBSD map carries whatever
% convention |ebsd.how2plot| says, so |ij2EbsdSquare| is needed to rotate
% the image data into the map's plotting convention.

ebsdOut = job.undistortedList(1).ebsd;

figure
nextAxis
plot(ebsdOut('indexed'), ebsdOut('indexed').orientations, ebsdOut.how2plot, 'coordinates','on')
title('Undistorted MTEX EBSD map (Copper IPF out of screen)','Color','k')

for n = 1:numel(job.undistortedList)

  im = job.undistortedList(n).img;
  if size(im,3) > 1, im = mean(im,3); end

  nextAxis
  plot(ebsdOut, ij2EbsdSquare(ebsdOut,im), ebsdOut.how2plot, 'coordinates','on')
  mtexColorMap gray
  title(['Undistorted MTEX image ' num2str(n)],'Color','k')
end

%% Turn the voids into a phase
% That is the end of the distortion correction. Everything below is
% ordinary MTEX on an EBSD map that now overlays its images pixel for
% pixel.
%
% Undistorting leaves a ragged border of pixels that only one of the maps
% covered, so both the map and the image are first cropped to the largest
% rectangle that is fully inside the EBSD data.

ebsd = job.undistortedList(1).ebsd;
bse  = job.undistortedList(4).img;

[~,~,~,keepGrid] = FindLargestRectangles(~isnan(job.undistortedList(1).img));
ebsd = gridify(ebsd(ij2EbsdSquare(ebsd,keepGrid)));
ebsd.how2plot = job.undistortedList(1).how2plot; % gridify forgets it
bse = reshape(bse(keepGrid),size(ebsd));

%%
% A void is a hole, so EBSD has nothing to index there and the aligned
% backscatter image is the only evidence of where it is. Thresholding that
% image gives a mask, and the mask becomes a phase of its own — which is
% what makes the voids available to |calcGrains|, |grains.boundary| and
% every other MTEX tool below.

phasesBse = ij2EbsdSquare(ebsd,(bse<0.8)); % voids = 1, copper = 0

voidPhase = crystalSymmetry('1','mineral','voids','color',str2rgb('DarkBlue'));

ebsd(phasesBse).rotations = rotation('euler',0,0,0);
ebsd(phasesBse).CS = voidPhase;
ebsd = gridify(ebsd);
ebsd.how2plot = job.undistortedList(1).how2plot;

display(ebsd)

%% Void size distribution
% Most voids are small — of the order of ten pixels.

grainsVoids = calcGrains(ebsd('indexed'),'angle',10*degree);

display(grainsVoids('voids'))

figure
histogram(grainsVoids('voids'),grainsVoids('voids').area,50);
xlabel('void area ({\mu}m^2)');

figure
histogram(grainsVoids('voids'),grainsVoids('voids').numPixel,50);
xlabel('void area (pixels)');

figure
histogram(grainsVoids('voids'),grainsVoids('voids').diameter/ebsd.dPos,50);
xlabel('void diameter (pixels)');

%% Copper grains and boundaries
% The copper grains are reconstructed from the copper phase alone. Leaving
% the voids out of |calcGrains| is deliberate: the boundaries are then
% drawn straight *through* the voids rather than around them, which is what
% lets a void be attributed to the boundary it sits on.
%
% The missing points are filled in, carrying |grainId| along, because that
% is what makes each boundary segment identifiable later by the EBSD map
% points on either side of it (|gBs.ebsdId|).

[~,ebsd('Copper').grainId] = calcGrains(ebsd('Copper'),'angle',10*degree);
ebsdCopper = gridify(smooth(ebsd('Copper'),'fill'));
ebsdCopper.how2plot = ebsd.how2plot;

[grains,ebsdCopper('Copper').grainId] = calcGrains(ebsdCopper('Copper'),'angle',10*degree);

% naming both phases excludes the map border
gBs = grains.boundary('Copper','Copper');

% triple point segment triplets rather than triplePoints, so that they can
% be treated the same way as gBs below
tPs = grains.triplePoints('Copper','Copper','Copper');
tPGbs = grains.boundary(tPs.boundaryId);

%%
% The map with the voids overlaid.

figure; newMtexFigure('layout',[2,1]);
nextAxis
plot(ebsd,ebsd.bc,ebsd.how2plot,'micronbar','off');
mtexColorMap gray; hold on
plot(ebsd('indexed'),ebsd.how2plot,'FaceAlpha',0.7);
mtexTitle('Band Contrast and Phases');
nextAxis
plot(ebsd('Copper'),ebsd('Copper').orientations,'FaceAlpha',0.5,...
  ebsd.how2plot,'micronbar','on'); hold on
plot(gBs,ebsd.how2plot,'linewidth',1,'linecolor','g');
plot(tPGbs,ebsd.how2plot,'linewidth',2,'linecolor','m');
plot(ebsd('voids'),zeros(size(ebsd('voids'))),ebsd.how2plot);
mtexTitle('Copper Orientations (IPF out of screen), Grain Boundaries and Voids');

%% Find the nearest boundary to each void
% Every boundary segment in |gBs| runs between the two neighbouring EBSD
% map points recorded in |gBs.ebsdId|. Painting those points with the index
% of their segment gives |gbPosMap| — an image of the boundary network in
% which every boundary pixel knows which segment it belongs to.

gbPosMap = zeros(size(ebsdCopper));
[gbEbsdIdList,ia] = unique(gBs.ebsdId,'stable');

% ia indexes the unique values of gBs.ebsdId, gBIdList holds row indices to gBs
gBIdList = repmat((1:length(gBs))',[1,2]);
gbPosMap(id2ind(ebsdCopper,gbEbsdIdList)) = gBIdList(ia);

% and the same for the triple point segments
tpPosMap = zeros(size(ebsdCopper));
[tpEbsdIdList,ia] = unique(tPGbs.ebsdId,'stable');
tpIdList = repmat((1:length(tPGbs))',[1,2]);
tpPosMap(id2ind(ebsdCopper,tpEbsdIdList)) = tpIdList(ia);

%%
% Indexing gets confusing here, so plot as we go.

figure
nextAxis
plot(ebsdCopper,gbPosMap,ebsdCopper.how2plot); colormap gray; hold on; mtexColorbar
plot(gBs,ebsdCopper.how2plot,'lineColor','g');
nextAxis
plot(ebsdCopper,tpPosMap,ebsdCopper.how2plot); colormap gray; hold on; mtexColorbar
plot(gBs,ebsdCopper.how2plot,'lineColor','g');
plot(tPGbs,ebsdCopper.how2plot,'lineColor','m','lineWidth',1);

%%
% A Euclidean distance transform of that image answers both questions at
% once: |gbDist| is how far each map point is from the nearest boundary
% pixel, and |gbNearest| is which segment that nearest pixel belongs to.

[gbDist, ix] = bwdist(gbPosMap);
gbNearest = gbPosMap(ix);

[tpDist, ix] = bwdist(tpPosMap);
tpNearest = tpPosMap(ix);

figure; newMtexFigure('layout',[2,2]);
nextAxis
plot(ebsdCopper,gbDist); colormap gray; hold on; mtexColorbar
plot(gBs,'lineColor','g');
nextAxis
plot(ebsdCopper,gbNearest); colormap gray; hold on; mtexColorbar
plot(gBs,'lineColor','g');
nextAxis
plot(ebsdCopper,tpDist); colormap gray; hold on; mtexColorbar
plot(tPGbs,'lineColor','m','lineWidth',1);
nextAxis
plot(ebsdCopper,tpNearest); colormap gray; hold on; mtexColorbar
plot(tPGbs,'lineColor','m','lineWidth',1);

%%
% Restricting those two maps to the void pixels leaves, for every void
% pixel, the boundary segment nearest to it and how far away that is.

voidsMapgb = nan(size(ebsd));
voidsMapgb(ebsd.phase==ebsd('voids').phase(1)) = gbNearest(ebsd.phase==ebsd('voids').phase(1));
voidsDistgb = (~isnan(voidsMapgb)) .* gbDist;

voidsMaptp = nan(size(ebsd));
voidsMaptp(ebsd.phase==ebsd('voids').phase(1)) = tpNearest(ebsd.phase==ebsd('voids').phase(1));
voidsDisttp = (~isnan(voidsMaptp)) .* tpDist;

voidsListgb = unique(voidsMapgb(~isnan(voidsMapgb)));

% a triple point is a triplet of segments, so reunite each match with its
% two partner segments
[r,~] = ind2sub(size(tPs.boundaryId),unique(voidsMaptp(~isnan(voidsMaptp))));
voidsListtp = sub2ind(size(tPs.boundaryId),repmat(r(:),[1 3]),repmat(1:3,[numel(r) 1]));

figure; newMtexFigure('layout',[2 2]);
nextAxis
plot(ebsdCopper,voidsMapgb); colormap gray; hold on; mtexColorbar
plot(gBs(voidsListgb),'lineColor','g');
nextAxis
plot(ebsdCopper,voidsDistgb); colormap gray; hold on; mtexColorbar
plot(gBs(voidsListgb),'lineColor','g');
nextAxis
plot(ebsdCopper,voidsMaptp); colormap gray; hold on; mtexColorbar
plot(tPGbs(voidsListtp),'lineColor','m','lineWidth',1);
nextAxis
plot(ebsdCopper,voidsDisttp); colormap gray; hold on; mtexColorbar
plot(tPGbs(voidsListtp),'lineColor','m','lineWidth',1);

%% How close is "on the boundary"?
% A void one pixel away from a boundary may be genuinely off it, or the
% alignment may be off by a pixel. The honest threshold is therefore the
% registration's own accuracy: the 95th percentile of the residual shift
% left over from each of the three correction steps, summed.

voidsList_threshPix = ...
    prctile(sqrt((job.fitError(1).xShiftsXcf/job.fitError(1).dx).^2 + ...
                 (job.fitError(1).yShiftsXcf/job.fitError(1).dy).^2),95) + ...
    prctile(sqrt((job.fitError(2).xShiftsXcf/job.fitError(2).dx).^2 + ...
                 (job.fitError(2).yShiftsXcf/job.fitError(2).dy).^2),95) + ...
    prctile(sqrt((job.fitError(3).xShiftsXcf/job.fitError(3).dx).^2 + ...
                 (job.fitError(3).yShiftsXcf/job.fitError(3).dy).^2),95);

disp(['Threshold distance from g.b. (pixels): ' num2str(voidsList_threshPix)]);

%%
% Split the void pixels into those on a boundary, near one, and away from
% one. Note this counts void *pixels*, not voids.

voidsListGb_on      = voidsMapgb(voidsDistgb<=1 & ~isnan(voidsMapgb));
voidsListGb_near    = voidsMapgb(voidsDistgb>1 & voidsDistgb<=voidsList_threshPix & ~isnan(voidsMapgb));
voidsListGb_notNear = voidsMapgb(voidsDistgb>voidsList_threshPix & ~isnan(voidsMapgb));

numVoidPix = nnz(~isnan(voidsMapgb));

% inconsistent void pixel counts between the two maps would be a bug
assert(nnz(~isnan(voidsMapgb))==nnz(~isnan(voidsMaptp)));

disp([num2str((numel(voidsListGb_near)+numel(voidsListGb_on))/numVoidPix*100) ...
  ' % of void pixels are on or near a GB (including GBs attached to TPs).']);
disp([num2str(numel(voidsListGb_notNear)/numVoidPix*100) ...
  ' % of void pixels are far from a GB.']);

%%
% The same for triple junctions, with the extra step of reuniting each
% matched segment with the other two of its triplet. Splitting "on" from
% "near" is not meaningful here, since only pixels directly intersecting
% the triple point count as "on".

t1 = voidsMaptp(voidsDisttp<=1 & ~isnan(voidsMaptp));
[r,~] = ind2sub(size(tPs.boundaryId),t1);
voidsListTp_on = sub2ind(size(tPs.boundaryId),repmat(r(:),[1 3]),repmat(1:3,[numel(r) 1]));

t1 = voidsMaptp(voidsDisttp>1 & voidsDisttp<=voidsList_threshPix & ~isnan(voidsMaptp));
[r,~] = ind2sub(size(tPs.boundaryId),t1);
voidsListTp_near = sub2ind(size(tPs.boundaryId),repmat(r(:),[1 3]),repmat(1:3,[numel(r) 1]));

voidsListTp_notNear = voidsMaptp(voidsDisttp>voidsList_threshPix & ~isnan(voidsMaptp));

disp([num2str((size(voidsListTp_near,1)+size(voidsListTp_on,1))/numVoidPix*100) ...
  ' % of void pixels are on or near a TP.']);
disp([num2str(size(voidsListTp_notNear,1)/numVoidPix*100) ...
  ' % of void pixels are far from a TP.']);

%%
% Drop the repeats and everything not near a void.

voidsListgb = unique([voidsListGb_on;voidsListGb_near]);
voidsListtp = unique([voidsListTp_on;voidsListTp_near]);

%% Which boundaries resist voids?
% Comparing the misorientation distribution of the boundaries that carry
% voids against the distribution over all boundaries. This is counted per
% boundary segment on or near a void rather than per void, so a large void
% crossing many segments weighs more than a small one.

mdf_voidsGb = calcDensity(gBs(voidsListgb).misorientation);
mdf_voidsTp = calcDensity(gBs(voidsListtp).misorientation);
mdf_all = calcDensity(gBs.misorientation);

figure; newMtexFigure;
plot(ebsd('Copper'),ebsd('Copper').orientations,'FaceAlpha',0.3,ebsd.how2plot); hold on
plot(gBs,ebsd.how2plot,'linecolor',str2rgb('gray'));
plot(ebsd('voids'),zeros(size(ebsd('voids'))),ebsd.how2plot); colormap gray; clim([0 1]);
plot(gBs(voidsListGb_near),ebsd.how2plot,'linecolor',str2rgb('DarkGreen'),'linewidth',3);
plot(gBs(voidsListGb_on),ebsd.how2plot,'linecolor',str2rgb('LightGreen'),'linewidth',3);
plot(tPGbs(voidsListTp_near),ebsd.how2plot,'linecolor',str2rgb('DarkRed'),'linewidth',3);
plot(tPGbs(voidsListTp_on),ebsd.how2plot,'linecolor','m','linewidth',3);

%%
% The answer is in the next two plots: boundaries around 60 degrees about
% [111] are markedly under-represented among the void boundaries. Those are
% the sigma-3 twin boundaries of FCC copper, and they resist void
% formation. No triple junction type stands out the same way in this
% material.

figure
newMtexFigure('figSize','tiny','outerplotspacing',30);
plotAngleDistribution(mdf_all,'DisplayName','All GBs'); hold on
plotAngleDistribution(mdf_voidsGb,'DisplayName','Void GBs');
plotAngleDistribution(mdf_voidsTp,'DisplayName','Void TPs');
plotAngleDistribution(ebsd('Copper').CS,ebsd('Copper').CS,'antipodal','DisplayName','Uniform MDF');
legend('show','Location','northwest');
xlabel('Misorientation angle / degrees');
ylabel('Frequency / mrd');

figure; newMtexFigure('layout',[2,2],'figSize','large','outerplotspacing',30,'innerplotspacing',50);
nextAxis(1,1); plotAxisDistribution(mdf_all,'colorRange','equal'); mtexTitle('All GBs');
nextAxis(1,2); plotAxisDistribution(mdf_voidsGb,'colorRange','equal'); mtexTitle('Void GBs');
nextAxis(2,1); plotAxisDistribution(mdf_voidsTp,'colorRange','equal'); mtexTitle('Void TPs');
nextAxis(2,2); plotAxisDistribution(ebsd('Copper').CS,ebsd('Copper').CS,'antipodal','colorRange','equal');
mtexTitle('Uniform MDF');
mtexColorbar;

%% Helper
% Tile a distortedImg sequence, plotted in microns on a shared axis.

function plotImgList(list,ttl,which)

if nargin < 3, which = 'img'; end

figure('WindowState','maximized');
t = tiledlayout('flow','TileSpacing','tight','Padding','tight');
title(t,ttl);

for n = 1:numel(list)

  if strcmp(which,'xcf')
    cData = list(n).(list(n).setXCF(1).xcfImg);
  else
    cData = list(n).img;
  end

  nexttile
  imagesc('XData',list(n).dx .* (1:size(list(n).img,2)),...
    'YData',list(n).dy .* (1:size(list(n).img,1)),...
    'CData',cData);
  colormap gray; axis image on ij
end
linkaxes

end
