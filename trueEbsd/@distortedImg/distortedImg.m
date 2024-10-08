classdef distortedImg
    % Class to store distorted images + metadata for trueEBSD
    %
    % Constructor syntax
    %   disImg = distortedImg(img,distortionName)
    %   disImg = distortedImg("bc",distortionName,ebsd)
    %   disImg = distortedImg(img,distortionName,"dxy",[dx
    %   dy]);
    %   disImg = distortedImg(img,distortionName,"dxy",[dx])
    %
    %
    %  Inputs
    %
    % img – r*c*z numeric array
    % Input image values, where r = number of rows, c = number of columns,
    % z = image signal channels, e.g. scalar (z=1) for a greyscale image,
    % RGB (z=3) for a colour image. z is not restricted to 1 (greyscale) or
    % 3 (colour); for example, some SEM images have 5 signal channels from
    % 5 detector diodes, so z=5 in this case.
    %
    % distortionName – @char array
    % •	Name of the distortion type between the current image (which
    % contains this type of distortion) and the image after it (without
    % this distortion). The ‘distortionName’ value for the reference image
    % is ‘true’. •	Current available options: ‘shift’, ‘drift’, ‘tilt’,
    % ‘true’
    %
    % varargin{:} –
    % Optional inputs:
    % ebsd – @EBSD or @EBSDSquare MTEX object, same pixel positions as the image
    %  leave any image input field empty [] to load from UI
    %
    % Options as Name, Value pairs:
    %   'dxy', <pixel size>
    % - Pixel size in μm, numeric scalar, optional if you can read this
    %     from the EBSD object
    %
    % 'mapplottingConvention', <@plottingConvention>
    % @plottingConvention is a MTEX class
    % - Defaults in @distortedImg constructor to axis ij convention (+X
    %     points East, +Y points South, +Z points into screen) which is always
    %     correct for image file formats 
    % - needs to be specified for EBSD maps, convention depends on system
    %
    % 'highContrast', <1 or 0>
    % qualitative and subjective flag for ‘does this image have high edge contrast?’ 
    % - This affects how calcShifts() will work later: 
    %   (1) % cross-correlation of [n+1,n] image ROI pairs will always try to start
    %   on an image with high edge contrast, and 
    %   (2) if the n+1th image has low edge contrast, image intensities from the nth image will be
    %   remapped onto the spatial coordinates of the n+1th image so that the
    %   next iteration ([n+2,n+1]) of cross-correlation has high contrast
    %   features for image registration. 
    % 
    % 'edgePadWidth',<integer number of pixels>
    %  - edge transform kernel radius in pixel units.
    %
    % ‘pixelTime’, <number>
    %  - dwell time per pixel in milliseconds
    %  - Not currently used by any functions, but useful if you need to
    %      calculate drift rate
    %
    % ‘setXCF’, <setXCF @struct>
    % @struct describing region of interest (ROI) size and spacing for image cross-correlation in calcShifts().
    % - The default values are calculated in pixelSizeMatch()  --
    %     reasonable values for a typical microstructure map, but you can set
    %     this property yourself in job.resizedList{n}.setXCF.
    % - Set this after running pixelSizeMatch(job) and before
    %     calcShifts(job), because the ROI are defined in pixel units, so it
    %     should be performed on the re-sized/scaled data.
    % - Some guidance for selecting the ROI size:
    %       - A Hann windowing function is applied to each ROI to eliminate
    %           edge artefacts during cross-correlation. This tapers the image
    %           intensity at the edges of each ROI, so your ROI size should be
    %           selected to account for this (it needs to be bigger than you
    %           think).
    %      - ROI should be big enough to contain multiple non-collinear
    %           high contrast features in the image edge transform. This is
    %           needed for good image registration.
    %      - ROI should be big enough to contain multiple overlapping features
    %           in pairs of ROI in the nth and (n+1)th image, so that the local
    %           image shifts can be measured from the image cross-correlation peak
    %           position. A good rule of thumb is that the ROI should be about 4
    %           times as wide than the image shift you want to measure. 
    %       - ROI should be small enough that the differences between ROI can be
    %           approximate by rigid body shifts, i.e. the degree of image warping
    %           is small.
    %
    %

    properties %properties we import from the constructor function
        img = [] % from img
        distortionName = '' % @char array containing name of the distortion type
        % between the current image (which contains this type of
        % distortion) and the image after it (without this distortion). The
        % distortionName value for the reference image is 'true'.
        dx = nan % from option 'dxy', <pixel size> - scalar um, 
        % optional if you can read this from the EBSD object
        dy = nan % same as 'dx' 
        ebsd = EBSD % @EBSD or @EBSDSquare MTEX object, same pixel positions 
        % as disImg.img
        pixelTime = 0 % EBSD exposure time or image pixel dwell time in ms
        mapPlottingConvention = plottingConvention(-vector3d.Z,vector3d.X) %from 
        % option 'mapplottingConvention', <MTEX @plottingConvention> - defaults 
        % to 'axis image'
        highContrast = nan %from option 'highContrast', <1 or 0>, scalar/logical 
        % 1 = good edge contrast, 0 = poor edge contrast
        edgePadWidth = 1 %from option 'edgePadWidth',<integer number of pixels> 
        % used in get.edge. Increase this to enhance edge contrast in images with blurred boundaries
    end

    properties  %properties related to trueEbsd algorithm
        pos = vector3d % vector3d position of every pixel, computed later in pixelSizeMatch
        setXCF = {} % settings for image registration, used in calcShifts, defaults computed in constructor function
    end

    properties (Dependent=true) %properties we get
        distortionModel %char that goes into fRunDIC input
        edge
    end

    methods
        function disImg = distortedImg(img,distortionName,varargin)
            % Construct an instance of this class
            % inputs / outputs described in script header

            %empty object
            if nargin == 0, return; end

            [img,distortionName,varargin{:}] = convertStringsToChars(img,distortionName,varargin{:});
            disImg.distortionName = argin_check(distortionName,{'char'});


            if isempty(img)
                disp(['Load image file:']);
                [f1,p1]=uigetfile('*.*');
                disImg.img = im2double(imread(fullfile(p1,f1)));
            elseif isa(img,'numeric')
                disImg.img = im2double(img);
            end

            % handle optional inputs
            disImg.pixelTime = get_option(varargin,'pixelTime',0,{'double';'single';'uint8';'uint16';'uint32'});
            disImg.mapPlottingConvention  = get_option(varargin,'mapPlottingConvention',plottingConvention(-vector3d.Z,vector3d.X),{'plottingConvention'});

            % import EBSD object
            [disImg.ebsd,varargin] = getClass(varargin,'EBSD');  % includes EBSDSquare and EBSDHex
            if ~isempty(disImg.ebsd)
                if isa(disImg.ebsd,'EBSDhex')
                    error('Square grid EBSD maps only');
                elseif isa(disImg.ebsd,'EBSDsquare')
                    %do nothing
                else
                    disImg.ebsd = gridify(disImg.ebsd);
                    if isa(disImg.ebsd,'EBSDhex')
                        error('Square grid EBSD maps only');
                    end
                end
                disImg.dx = disImg.ebsd.dx;
                disImg.dy = disImg.ebsd.dy;

                % extract img if required
                if isa(img,'char')
                    disImg.img = im2double(ebsdSquare2ij(disImg.ebsd,img,disImg.mapPlottingConvention));
                end
            else
                dxy = get_option(varargin,'dxy',0,{'double';'single';'uint8';'uint16';'uint32'});
                disImg.dx = dxy(1);
                disImg.dy = dxy(end);
            end

            % default XCF settings
            % can be set using options setXCF, <setXCF @struct> otherwise
            % use default below
            %
            setXCF.ROISize=0; %this will be changed anyway during pixelSizeMatch
            %estimate this for 'typical' EBSD map - maybe 25 grains across?
            setXCF.NumROI.x = 24; % guess ~1 ROI per grain in FOV?
            setXCF.NumROI.y = round(setXCF.NumROI.x * size(disImg.img,1)/size(disImg.img,2)); % follow SEM aspect ratio x:y
            setXCF.XCFMesh=250; % correlation peak upsampling, default 250
            setXCF.xcfImg='edge'; % decide whether to use edge transform or image intensities for cross-correlation
            % % or load a custom one (but ROISize will always be rewritten
            % % when resizing images)

            disImg.setXCF = get_option(varargin,'setXCF',{setXCF},{'cell'});
            disImg.highContrast = get_option(varargin,'highContrast',0,{'double';'single';'uint8';'uint16';'uint32';'logical'});
            disImg.edgePadWidth = get_option(varargin,'edgePadWidth',1,{'double';'single';'uint8';'uint16';'uint32';'logical'});

        end % end constructor function


        function out = get.edge(disImg)
            % generate edge transform of image
            % input - @distortedImg properties
            %         - img = image values
            %         - edgePadWidth = edge transform kernel radius in pixel units.
            %           use option 'edgePadWidth',<integer number of
            %           pixels> in constructor
            % output - @distortedImg disImg.edge = edge transform
            out = GenBoundaryMap(disImg.img,disImg.edgePadWidth);
        end

        function out = get.distortionModel(disImg)
            % convert distortionName into inputs that fRunDIC understands
            % used to fit a continous function (evaluated per pixel) to
            % the measured image shifts (evalulated per region of interest)
            % input = @distortedImg property distortionName - @char array
            % output =  @distortedImg property distortionModel - @char array
            switch disImg.distortionName
                case 'shift'
                    out = {'poly11'};
                case 'shift-interp'
                    out = {'poly11','interpolate'};
                case 'drift'
                    out = {'linearinterp'};
                case 'drift-shift'
                    out = {'poly11','linearinterp'};
                case 'drift-interp'
                    out = {'linearinterp','interpolate'};
                case 'tilt'
                    out = {'projective','poly11','poly22'};
                case 'tilt-interp'
                    out = {'projective','poly11','poly22','interpolate'};
                case 'true'
                    out = {''};
                otherwise
                    error('uknown distortion model');
            end
        end
    end
end