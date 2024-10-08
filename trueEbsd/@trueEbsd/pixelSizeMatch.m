function job = pixelSizeMatch(job, varargin)
% @trueEbsd method
% resize all job.interImgs images to the same pixel size
%
% Syntax
% job = pixelSizeMatch(job, pixelsize)
% job = pixelSizeMatch(job) % default pixel size = smallest of job.imgList{:}.dx
% job = pixelSizeMatch(job, 0) % pixelsize = 0 - use default pixel size
% job = pixelSizeMatch(job,pixelsize,'offsetMatch','centre','extentMatch','largest');
% job = pixelSizeMatch(job,pixelsize,'offsetMatch','centre','extentMatch','smallest');
% job = pixelSizeMatch(job,pixelsize,'offsetMatch','topLeft','extentMatch','largest');
% job = pixelSizeMatch(job,pixelsize,'offsetMatch','topLeft','extentMatch','smallest');
% 
% 
% Inputs
% job = @trueEbsd object
% job.imgList contains input data as @distortedImg cell array
% 
% Optional inputs
% pixelsize = numeric scalar (default = smallest)
%
%  Name, Value pair options:
% 'offsetMatch', 'centre' or 'topLeft' (default = centre)
%    How to overlay the images if the fields of view are different sizes.
%    'centre' means the images will start with the central pixel aligned.
%    'topLeft' means the images will start with the top left corner aligned.
%    You should select the option that minimises the total rigid body shift
%    between images.
% 'extentMatch', 'largest' or 'smallest' (default = largest) 
%    Whether to crop images to the smallest image field of view (‘smallest’), or
%    extrapolate images to the largest image field of view ( largest’) by
%    repeating the edge and corner pixels to fill the space.
%

% default params
temp = [job.imgList{:}];
pixelsize = min([temp.dx]);
offsetMatch = 'centre'; %default
extentMatch = 'largest'; %default
if nargin>1
    [pixelsize,varargin] = getClass(varargin,'numeric',pixelsize);
    offsetMatch = get_option(varargin,'offsetMatch',offsetMatch,{'char'});
    extentMatch = get_option(varargin,'extentMatch',extentMatch,{'char'});
end
if ~(pixelsize>0) %0 or nan    
    pixelsize = min([temp.dx]);
    disp(['using default pixel size of ' num2str(pixelsize) ' um, minimum from imgList']);
else

    disp(['using user-defined pixel size of ' num2str(pixelsize) ' um']);
end


%preallocate variables
intpi=cell(1,numel(job.imgList)); % image interpolant
% interpolant function for old image
% with offsetMatch to set zero position
% @griddedInterpolant for image values (numeric r*c*z rray)

extentMax = zeros(numel(job.imgList),2); extentMin=extentMax;
v=cell(numel(job.imgList),1); %xyz pixel position vectors

%% Step 1. calculate original pixel position grids and  set up interpolants
% main outputs in this section:
% intpi{n} = griddedInterpolant object describing xyz pixel positions for image job.imgList{n}.img
% v{n} = grid vectors defining pixel positions of job.resizedList{n}.img
% v{1} = 1*r numeric vector containing row(Y) positions
% v{2} = 1*c numeric vector containing col(X) positions 
% v{3} = [1:z] for z image (colour) channels
% uses ndgrid convention in MATLAB (so [r,c,z] not [x,y,z])


for n = 1:numel(job.imgList) %loop through all images in workflow
    %old pixel size
    pxn=job.imgList{n}.dx;

    %preallocate variables
    v{n} = cell(1,ndims(job.imgList{n}.img)); % cell array to hold xy pixel position lists for each image
    g=v{n};
    for nd = 1:ndims(job.imgList{n}.img)
        if nd<3
            %create xy vector
            %initialise as [0,0] = centre of top left pixel
            v{n}{nd} = pxn*(0:size(job.imgList{n}.img,nd)-1);
            % handle offset between images in job.imgList{:}
            switch offsetMatch
                case 'centre'
                    %[0,0] = centre of image FOV
                    v{n}{nd} = v{n}{nd} - (max(v{n}{nd})/2);
                case 'topLeft'
                    %[0,0] = top left corner of image FOV -->
                    %centre of top left pixel = [pxn/2,pxn/2]
                    v{n}{nd} = v{n}{nd} + (pxn/2);
            end
        else
            %handle 3rd image as colour channels (pixel depth = 1 and don't
            %match anything) -- but don't assume it's always RGB 3
            %channels - e.g. forescatter images have 5 channels
            v{n}{nd} = 1:size(job.imgList{n}.img,nd);
        end
    end
    [g{:}] = ndgrid(v{n}{:}); % pixel positions
    job.imgList{n}.pos = vector3d(g{2}(:,:,1),g{1}(:,:,1),zeros(size(g{1}(:,:,1))));
    % intp = interpolant function
    % intpi = interpolant for image
    intpi{n}=griddedInterpolant(g{:},job.imgList{n}.img,'linear','none');

    % calculate exent of common pixel position grid
    for nd2=1:2 %only do this for xy coordinates
        extentMax(n,nd2)= max(intpi{n}.GridVectors{nd2} + pxn/2);
        extentMin(n,nd2)= min(intpi{n}.GridVectors{nd2} - pxn/2);
    end


end

%%  Step 2. calculate common pixel position grid
% with optional extra dimensions
%
% main output in this section:
% vc = grid vectors defining pixel positions of common image grid in job.resizedList{:}
% vc{nd2} = row-col pixel positions of in length units (default um)
% vc{1} = 1*r numeric vector containing row(Y) positions
% vc{2} = 1*c numeric vector containing col(X) positions 
% uses ndgrid convention in MATLAB (so [r,c] not [x,y])
% colour channels are left unchanged so numel(vc) = 2
% eventually goes into job.resizedList{n}.pos (via gn{:} - defined in Step 4)
% 

%1. first figure out how many pixels you need
%take into account half-pxn step at image edges
vc=cell(1,2); %common xy position vectors
offsetPos = zeros(1,2); % xy (east/south) position offset
for nd2=1:2 %only do this for xy coordinates
    switch extentMatch
        %you could do min:pixelsize:max but this method avoids rounding instabilities
        case 'largest'
            % vc{nd2} = min(extentMin(:,nd2)) : pixelsize : max(extentMax(:,nd2)) + eps*10000;
            npix = floor(eps*10000 + (max(extentMax(:,nd2))- min(extentMin(:,nd2)))/pixelsize);
        case 'smallest'
            % vc{nd2} = max(extentMin(:,nd2)) : pixelsize : min(extentMax(:,nd2)) + eps*10000;
            npix = floor(eps*10000 + (min(extentMax(:,nd2))- max(extentMin(:,nd2)))/pixelsize);            
        otherwise
            error("invalid extentMatch option: use either 'largest' or 'smallest'");
    end

    switch offsetMatch
        case 'centre' %[0,0] = centre of image FOV
            % 2. create xy pos vector with the new length (npix) and
            % pixelsize, then ...
            % 3. adjust vector offset to put [0,0] back in the middle
            offsetPos(nd2) = (-(npix+1)/2)*pixelsize;
           
            % sanity check:
            % if npix is even, [0,0] is at a pixel edge
            % but if npix is odd, [0,0] position is a pixel centre

        case 'topLeft' %[0,0] = top left corner of image FOV -->
            %centre of top left pixel = [pixelsize/2,pixelsize/2]
            offsetPos(nd2) = (-0.5)*pixelsize;
        otherwise
            error("invalid offsetMatch option: use either 'centre' or 'topLeft'");
    end
    %convert to length units (um) and apply position offset between images
    vc{nd2} = pixelsize*(1:npix) + offsetPos(nd2); 
end

%%  Step 3. estimate a sensible cross-correlation ROI size for the common grid
% this eventually goes into @distortedImg property setXCF.ROISize
% job.resizedList{n}.setXCF{1}.ROISize = roiSizePix
% each ROI should contain at least ~10 visible boundary segments - 3 grains
% across
% but only the central quarter of the ROI is visible because of the
% windowing function, so it needs to be around 6 grains across
% guess that a 'typical' EBSD map is maybe 25 grains across, so each ROI
% should be at least 1/4th of the map width numel(vc{1})

roiSizePix=2^(ceil(log2(numel(vc{1})/4))); 

%%  Step 4. resample images and edge transforms on common grid 
% Main inputs in this section:
% intpi{n} = griddedInterpolant object describing image values
%   for image job.imgList{n}.img from Step 1
% v{n} = original (old) image grid position vectors for each image in 
%   job.imgList{n} from Step 1
% vc = common (new) image grid row-col position vectors for all images in 
%   job.resizedList{:} from Step 2
%
% Calculated helper variables
% vn = common (new) image grid position vectors for all images in 
%   job.resizedList{:} (including colour channels)
% gn{:} = the same information as vn{n} but formatted as ndgrid

% Main outputs in this section:
% job.resizedList{n}.img = image values resampled on the common pixel position 
%   grid gn{:} using intpi{n} griddedInterpolant 
% job.resizedList{n}.pos = @vector3d of pixel positions from gn{:}
% job.resizedList{n}.ebsd = (if it exists) @EBSD where ebsd.id is resampled
% on the common pixel position grid gn{:} using intpi{n}
% griddedInterpolant using nearest neighbour interpolation (i.e.
% orientations are not actually interpolated)


% write outputs to job.resizedList
%initialise job.resizedList{n} as copy of job.imgList{n}
job.resizedList = job.imgList;
for n = 1:numel(job.imgList)

    vn=v{n}; %new vector
    gn=cell(size(v{n})); % ndgrid positions for image job.resized{n}.img, constructed from vn
    vn(1:2) = vc; %replace xy with common vectors but leave colour channels alone
    [gn{:}] = ndgrid(vn{:}); %new grid

    % resample image values on new grid
    job.resizedList{n}.img = intpi{n}(gn{:});    
    job.resizedList{n}.dx = pixelsize;
    job.resizedList{n}.dy = pixelsize;
    job.resizedList{n}.pos = vector3d(gn{2}(:,:,1),gn{1}(:,:,1),zeros(size(gn{1}(:,:,1))));
     

    % handle ebsd variables too
    if ~isempty(job.resizedList{n}.ebsd)
        %construct new EBSD position vectors
        % use function tools/ij2EbsdSquare to handle conversion between ij indexing
        % positions (MATLAB convention) and EBSD xyz positions (convention
        % depends on system)
        [posEbsdX, posEbsdY] = nwse2EbsdPos(job.imgList{n}.ebsd,job.resizedList{n}.pos.x, job.resizedList{n}.pos.y);
        posEbsdX=ij2EbsdSquare(job.imgList{n}.ebsd, posEbsdX);
        posEbsdY=ij2EbsdSquare(job.imgList{n}.ebsd, posEbsdY);
        

        %transform image ij coordinates back into ebsd.pos xyz convention
        % function tools/ebsdSquare2ij is the reverse of ij2EbsdSquare 
        % find the nearest ebsd point on sampled grid
        intpi{n}.Values = ebsdSquare2ij(job.imgList{n}.ebsd, 'id');
        intpi{n}.Method = 'nearest';
        % resample EBSD map on the original common sampling grid
        %remember to convert back into EBSDSquare indexing!
        ebsdNewId = ij2EbsdSquare(job.imgList{n}.ebsd, intpi{n}(gn{1}(:,:,1),gn{2}(:,:,1)));
        %remove nan and zeros (only positive integers)
        ix = ~isnan(ebsdNewId) & ebsdNewId>0;

        % resize / resample EBSD map
        % don't use EBSDsquare/interp because that assumes everything is on a
        % regular grid, which is sometimes ok but grabbing grabbing ebsd.id from the
        %griddedInterpolant intpi (i.e. the method in undistort(job)) is
        %more general.
        %
        %create EBSD fields with just 0 values, then repopulate
        %with the correct fields
        ebsdMap0 = zeros(size(ebsdNewId)); %zeros matrix
        prop1 = job.imgList{n}.ebsd.prop; %struct containing map props
        for fn = fieldnames(prop1)'
            prop1.(char(fn))= ebsdMap0;
            prop1.(char(fn))(ix)= job.imgList{n}.ebsd(ebsdNewId(ix)).prop.(char(fn));
        end
        rot1 = rotation.nan(size(ebsdMap0));
        rot1(ix) = job.imgList{n}.ebsd(ebsdNewId(ix)).rotations;
        phase1 = ebsdMap0;
        phase1(ix) = job.imgList{n}.ebsd(ebsdNewId(ix)).phase;
        % recreate EBSD object
        ebsd1 = EBSD(vector3d(posEbsdX, posEbsdY,ebsdMap0), rot1, phase1, ...
            job.resizedList{n}.ebsd.CSList, prop1);
        ebsd1.plottingConvention = job.imgList{n}.mapPlottingConvention;

        % write to output
        job.resizedList{n}.ebsd = gridify(ebsd1);
        % plottingConvention isn't passed on automatically in gridify so reassign this
        job.resizedList{n}.ebsd.plottingConvention = job.imgList{n}.mapPlottingConvention;

        %NOTE: don't try to use @EBSDsquare/interp because it only handles indexed EBSD
        %points when interpolating the map, this leads to e.g. ebsd.bc disappearing from unindexed
        %points
        % this will not work: job.resizedList{n}.ebsd = gridify(interp(job.resizedList{n}.ebsd,ebsdPosX, ebsdPosY));


        % translate ebsd object to match map image offset
        % use function tools/ebsdMapOffset to convert between map dirs (NWSE)
        % and positions vectors (XYZ)
        % can't use ebsd +- directly because that's for xyz not ij map
        % directions
        % offsetPos(1) = row shift, i.e. south
        % offsetPos(2) = column shift. i.e. east
        job.resizedList{n}.ebsd = ebsdMapOffset(job.resizedList{n}.ebsd,offsetPos(2),offsetPos(1));
    end

    % update ROI size in setXCF (from step 3), except leave reference image as it has no
    % 'next image' to correlate shifts to.
    if  n < numel(job.imgList)
        job.resizedList{n}.setXCF{1}.ROISize = roiSizePix;
    end

end