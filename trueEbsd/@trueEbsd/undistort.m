function job = undistort(job)
% 
% Syntax
% job = undistort(job)
% use this after running calcShifts(job)
%
% Inputs
% job.resizedList{1:n};
% job.shifts{1:n-1};
%
% Outputs 
% job.undistortedList{1:n}
% 
% % Description 
% This shifts the position of each pixel for every image in
% job.resizedList{n}.img by the shifts calculated in calcShifts, and
% resamples the shifted image at the resizedList image positions using
% nearest-neighbour interpolation, so that the final output image is on a
% regular square grid. 
% In theory you could modify this to use linear or
% cubic interpolation schemes but nearest-neighbour interpolation is least
% likely to produce image intensity artefacts.
% 
% @EBSD orientation maps are handled differently to normal images because
% the image pixel intensities are not directly related to the orientations.
% Instead, the @EBSD ‘id’ property is resampled, also using
% nearest-neighbour interpolation, to construct the output @EBSD object. In
% this case, you must use nearest neighbour interpolation because you are
% resampling EBSD map pixel labels which do not have ‘in-between’ values.
% 
% Since shifts are calculated for [n, n+1] image pairs, they are also
% stacked sequentially. So 
% -	job.undistortedList{1} is shifted by everything in job.shifts{:} 
% -	job.undistortedList{2} is shifted by job.shifts{2:end} 
% -	job.undistortedLIst{end} is not shifted at all, it is the reference image.


%initialise undistorted images as copy of resizedList
job.undistortedList = job.resizedList;

%initialise undistorted xy image positions as originals
% this needs to be outside the 1:nn for-loop because we are stacking
% distortions in the image list from n:1
xTrue = job.resizedList{end}.pos.x;
yTrue = job.resizedList{end}.pos.y;

%add image shifts
% need to stack them in reverse order (most distorted has n=1, reference
%image has n=end)
for n=numel(job.resizedList):-1:1
    if n==numel(job.resizedList)
        % no shifts for reference image
        % xTrue = xTrue;
        % yTrue = yTrue;
    else
        for m = 1:numel(job.shifts{n})
            xTrue = xTrue + job.shifts{n}{m}.x;
            yTrue = yTrue + job.shifts{n}{m}.y;
        end
    end
    %% create scattered interpolant function and resample on original image grid
    if ndims(job.resizedList{n}.img) > 3 % can deal with 2D images or colour channels ('extruded-3D') but not larger dimensions
        error("can't handle 4D images")
    end
    %repeat for every image channel
    for p = 1:size((job.resizedList{n}.img),3) % 2D image (no colours)
        v1 = job.resizedList{n}.img(:,:,p); %extract 1 colour channel
        if p==1
            intpiNew = scatteredInterpolant([xTrue(:) yTrue(:)],v1(:),'nearest','none');%relatively fast: 6s for 1.8M pixels
        else
            %  faster to update intpiNew with new values than recalculate interpolant
            intpiNew.Values = v1(:);
        end
        v1New = intpiNew(job.resizedList{n}.pos.x,job.resizedList{n}.pos.y);% very slow step: 74s for 1.8M pixels
        job.undistortedList{n}.img(:,:,p) = v1New;
    end

    %the last image should always be the same as this is the reference image
    % i.e. job.undistortedList{end}==job.resizedList{end}
    % which is why numel(job.shifts)+1==numel(job.resizedList) if everything is
    % running correctly

    %% handle EBSD data interpolation
    if ~isempty(job.undistortedList{n}.ebsd)
      
        % don't try to mess with ebsd.pos xyz positions, just find the
        % nearest ebsd point on sampled grid
        intpiNew.Values = reshape(...
            ebsdSquare2ij(job.resizedList{n}.ebsd, 'id'),...
            [],1);

        % resample EBSD map on the original common sampling grid
        %remember to convert back into EBSDSquare indexing!
        ebsdNewId = ij2EbsdSquare(job.resizedList{n}.ebsd, ...
            intpiNew(job.resizedList{n}.pos.x,job.resizedList{n}.pos.y));
        %remove nan and zeros (only positive integers)
        ix = ~isnan(ebsdNewId) & ebsdNewId>0;
        % MTEX interp doesn't work on non-gridded data and
        % uses ebsdsquare.dx in its calculation - basically won't work for a
        % scatter point cloud

        %create EBSD fields with just 0 values, then repopulate
        %with the correct fields
        ebsdMap0 = zeros(job.resizedList{n}.ebsd.size);
        prop1 = job.resizedList{n}.ebsd.prop;
        for fn = fieldnames(prop1)'
            prop1.(char(fn))= ebsdMap0;
            prop1.(char(fn))(ix)= job.resizedList{n}.ebsd(ebsdNewId(ix)).prop.(char(fn));
        end
        rot1 = rotation.nan(job.resizedList{n}.ebsd.size);
        rot1(ix) = job.resizedList{n}.ebsd(ebsdNewId(ix)).rotations;
        phase1 = ebsdMap0;
        phase1(ix) = job.resizedList{n}.ebsd(ebsdNewId(ix)).phase;
        % recreate EBSD object
        ebsd1 = EBSD(job.resizedList{n}.ebsd.pos, rot1, phase1, ...
            job.resizedList{n}.ebsd.CSList, prop1);

        % write to output
        job.undistortedList{n}.ebsd = gridify(ebsd1);
        job.undistortedList{n}.ebsd.plottingConvention = job.resizedList{n}.ebsd.plottingConvention;

    end
end




