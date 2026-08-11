%% TrueEBSD distortion correction on a WC-Co composite
%
% Authors: Vivian Tong; Stefan Olovsjö, Seco Tools AB, R&D Materials and
% Technology, 737 82 Fagersta, Sweden.
% Contact: vivian.tong@extern.tu-freiberg.de
%
% TrueEBSD spatially aligns an EBSD map and SEM images of the same sample
% area, correcting the distortions between them so that every pixel
% overlays. The method is described in Tong et al.,
% <http://arxiv.org/abs/2605.00703 arXiv 2605.00703>.
%
% This script needs the TrueEBSD toolbox, which is *not* part of MTEX and
% is distributed separately under Apache-2.0:
%
%   <https://github.com/vtvivian/mtex-trueebsd>
%
%   addpath(genpath('<path to mtex-trueEbsd>'))
%
% Use the |mtex7-compat| branch with MTEX 7. Beyond MTEX it requires MATLAB
% R2024a or newer and the Image Processing, Curve Fitting, and Statistics
% and Machine Learning toolboxes. Runtime is minutes, not seconds.

%% Data Import
% Begin by loading an EBSD map together with the list of images we want to
% use alongside it. This dataset is an EBSD map of a WC-Co composite
% acquired at 20 kV accelerating voltage, plus four SEM images of the same
% sample area stored in |ebsd.opt.trueEbsdImgs|:
%
% # Band contrast (|ebsd.bc|) is used as the image for the EBSD map.
%
% # |fsdB3| is a colour image from the three FSD detectors mounted at the
% bottom of the EBSD camera, with the camera retracted by 20 mm relative to
% the EBSD map acquisition position;
%
% # |fsdT3| is a greyscale image from the same beam scan as |fsdB3|, from
% the FSD detectors at the top of the EBSD camera;
%
% # |fsdT1| is a greyscale image from the FSD detectors at the top of the
% EBSD camera, at the EBSD map acquisition position;
%
% # |fsdT10| is a greyscale image from the FSD detectors at the top of the
% EBSD camera, in EBSD map acquisition position, with the electron beam
% accelerating voltage lowered to 10 kV.
%
% * |ebsd.opt.trueEbsdImgs.pixSzImg| is the image pixel size in microns,
% the same for all four images.

mtexdata trueEbsdWCCo

display(ebsd)
display(ebsd.opt.trueEbsdImgs)

%% Set up the TrueEBSD job
% A <distortedImg.distortedImg.html distortedImg> holds one image or EBSD
% map, its pixel size, its plotting convention, and the name of the
% distortion separating it from the *next* entry in the sequence. The
% sequence runs from most distorted to ground truth, and the reference
% image carries |'true'|.

% some simple image denoising first
img = ebsd.opt.trueEbsdImgs;
img.fsdB3  = rescale(imboxfilt(img.fsdB3,3));
img.fsdT3  = rescale(imboxfilt(img.fsdT3,3));
img.fsdT1  = rescale(imboxfilt(img.fsdT1,3));
img.fsdT10 = rescale(imboxfilt(img.fsdT10,3));

imgList = createArray(5,1,'distortedImg');
imgList(1) = distortedImg('bc','shift-drift', ebsd, 'how2plot', ebsd.how2plot, ...
  'highContrast',1, 'edgePadWidth',3);
imgList(2) = distortedImg(img.fsdB3, 'true',  'dxy', img.pixSzImg, 'highContrast',1, 'edgePadWidth',5);
imgList(3) = distortedImg(img.fsdT3, 'shift', 'dxy', img.pixSzImg, 'highContrast',1, 'edgePadWidth',5);
imgList(4) = distortedImg(img.fsdT1, 'tilt',  'dxy', img.pixSzImg, 'highContrast',1, 'edgePadWidth',5);
imgList(5) = distortedImg(img.fsdT10,'true',  'dxy', img.pixSzImg, 'highContrast',1, 'edgePadWidth',3);

%%%
% The job is a <trueEbsd.trueEbsd.html trueEbsd> object built from that one
% sequence. It is a *value* class, so every workflow method returns the job
% and must be reassigned.

job = trueEbsd(imgList)

%%%
% Plot the as-imported sequence to check that the maps cover similar
% regions of the sample. Note how different the image contrasts look — this
% is why registration is done on edge transforms rather than raw values.

plotImgList(imgList,'TrueEBSD starting image sequence')

%% Resize images to match pixel size and field of view
% The EBSD map and the images cover the same sample area but have different
% pixel sizes. |pixelSizeMatch| resamples everything onto one common grid,
% so that pixel (i,j) means roughly the same place in each. Images are
% resampled by linear interpolation; EBSD data by nearest neighbour,
% because orientations and phase labels have no meaningful in-between.
%
% This is bookkeeping only — no distortion has been corrected yet.

pixSzIn = 0; % target pixel length in microns, or 0 for the smallest present
job = pixelSizeMatch(job,pixSzIn);

%%%
% The job now has a new property |job.resizedList| holding the output.

display(job)

%% [Optional] Change the cross-correlation ROI settings
% TrueEBSD registration cross-correlates pairs of regions of interest (ROI)
% between sequential images. ROI size and spacing are the tunable
% parameters, held in |job.resizedList(n).setXCF(m)| — one entry per
% distortion-model stage of that hop.
%
% |pixelSizeMatch| guesses values that are usually sensible for a
% polycrystal EBSD map, so most users never touch this. If you do, note the
% settings are in *pixels*, so they must be written after |pixelSizeMatch|
% (which creates the grid) and before |calcShifts|. A good rule of thumb is
% an ROI at least four times wider than the local shifts you expect, and it
% must be a power of two for the cross-correlation to work properly.
%
% Here we deliberately misjudge and choose an ROI box that is far too small
% for the EBSD map drift correction, to demonstrate the automatic retry
% below.

customSetXCF1.ROISize = 2^round(log2(32));   % deliberately too small
customSetXCF1.NumROI = struct;
customSetXCF1.NumROI.x = 40;                 % rule of thumb: as many ROI as grains across
customSetXCF1.NumROI.y = round(customSetXCF1.NumROI.x * ...
  size(job.resizedList(1).img,1)/size(job.resizedList(1).img,2)); % follow the aspect ratio
customSetXCF1.xcfImg = 'edge';               % correlate edge transforms, or 'img' for raw values

customSetXCF2 = customSetXCF1;
customSetXCF2.ROISize = 2^round(log2(128));

% assign a whole settings struct
job.resizedList(1).setXCF(2) = customSetXCF1;
job.resizedList(3).setXCF(1) = customSetXCF2;

% or overwrite individual properties
job.resizedList(1).setXCF(1).ROISize = 2^round(log2(64));
job.resizedList(3).setXCF(1).xcfImg = 'img';
job.resizedList(4).setXCF(1).xcfImg = 'img';
job.resizedList(5).setXCF(1).xcfImg = 'img';

%% Calculate local image shifts and fit a distortion model
% These are the images that will actually be cross-correlated — the edge
% transform where |xcfImg| is |'edge'|, the raw values where it is |'img'|.
% Edge transforms are what make an EBSD band contrast map comparable to a
% backscatter image at all.

plotImgList(job.resizedList,'TrueEBSD image sequence for cross-correlation','xcf')

%%%
% Now compute the local ROI shifts and fit them to the distortion model
% named on each hop. After each correction step the average ROI shifts (X,
% Y and length components) are printed to the command window.
%
% The |'fitErr'| flag means residual local shifts are re-measured after
% correction, but not included in the result. If that residual is small —
% around one pixel or less — the registration most likely succeeded.
%
% We just set the ROI box too small for the EBSD map drift correction, so
% its average residual shift comes out greater than 2 pixels. |calcShifts|
% responds by doubling the ROI size and retrying, and keeps doing so until
% either the residual drops below 2 pixels or the ROI outgrows the image.
%
% The exception is a hop named |'true'|, such as images 2 and 3 here, where
% nothing separates the pair. TrueEBSD takes those shifts to be identically
% zero and ignores the residual even when it is large.

job = calcShifts(job,'fitErr');

%%%
% The job now has |job.shifts| — a cell array, one entry per hop, each
% holding one <pairShifts.pairShifts.html pairShifts> per distortion-model
% stage — and |job.fitError|, the residuals measured after correction.

display(job)

%% Undistort
% This accumulates the shifts in reverse — the first map receives every
% hop's shift, the reference none — and resamples with nearest-neighbour
% interpolation, so no orientation or phase label is ever invented by
% averaging two real ones. The result is |job.undistortedList|, in which
% every pixel of every map can be directly overlaid.

job = undistort(job);

display(job)

plotImgList(job.undistortedList,'TrueEBSD image sequence after alignment')

%% Plot the aligned data as MTEX EBSD maps
% Plotting the images back onto the EBSD map is a good check that nothing
% is indexed upside down relative to the map. Images are stored and read by
% MATLAB in the |axis ij| convention, whereas an EBSD map carries whatever
% convention |ebsd.how2plot| says, so |ij2EbsdSquare| is needed to rotate
% the image data into the map's plotting convention.
%
% Note |fsdB3| is a three-channel colour image, and TrueEBSD carries all
% its channels through the workflow. Plotting values onto an EBSD map needs
% one value per pixel, so multi-channel images are averaged down to one
% channel here.

ebsdOut = job.undistortedList(1).ebsd;

figure
nextAxis
plot(ebsdOut('W C'), ebsdOut('W C').orientations, ebsdOut.how2plot, 'coordinates','on')
title('Undistorted MTEX EBSD map (WC IPF out of screen)','Color','k')

for n = 1:numel(job.undistortedList)

  im = job.undistortedList(n).img;
  if size(im,3) > 1, im = mean(im,3); end

  nextAxis
  plot(ebsdOut, ij2EbsdSquare(ebsdOut,im), ebsdOut.how2plot, 'coordinates','on')
  mtexColorMap gray
  title(['Undistorted MTEX image ' num2str(n)],'Color','k')
end

%% Finish
% That is the end of the distortion correction workflow. Every image and
% the EBSD map now overlay pixel for pixel, and |ebsdOut| is an ordinary
% MTEX EBSD map — anything you would normally do with one works from here,
% including using an aligned image as a per-pixel property or as a phase.
%
% For this dataset the intended follow-on is measuring the contiguity of
% the WC grains.

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
