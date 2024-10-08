function job = calcShifts(job, varargin)
% calculate ROI shifts from @trueEbsd variable
% basically works as a fRunDIC wrapper
% outputs job.shifts{imgList#}{distortionModel#}.shifts
%
% Syntax:
% job = calcShifts(job);
% job = calcShifts(job.'fitErr');
%
% Description:
%   1. Get job.resizedList{n}.edge for image cross-correlation
%   2. Decide image correlation sequence:
%     a. n = nStart: start with first image with high edge contrast (if no images have high edge contrast then nStart=1)
%     b. Call WP1_TrueEBSD/funcsV0/fRunDIC() to calculate XY shifts between [n+1, n] image pairs
%     c. Iterate n++ and repeat until n+1 = reference image (i.e. job.resizedList{end})
%     d. Whenever job.resizedList{n} has 'highContrast' = false, remap intensities using shifts measured from the previous iteration using the remapImage() which is nested inside this function
%     e. If nStart > 1, loop back round to n=nStart-1 and repeat for [nStart, nStart-1] image pairs
%         and repeat steps c. and d. with n-- until you reach n=1
%         The for-loop order runs the other way so that one image in the pair always
%         has 'highContrast' = true, which improves the correlation
%         quality.
%
% Inputs
%   job = @trueEbsd object. This function uses properties:
%   -	test image - job.resizedList{n}.img or .edge - set this in job.resizedList{n}.setXCF{1}.xcfImg, default is 'edge'
%   -	reference image - job.resizedList{n+1}.img or .edge
%   -	job.resizedList{n}.setXCF - structure containing ROI size and spacing settings,
%           description in @distortedImg/distortedImg and @trueEbsd/pixelSizeMatch
%   -	job.resizedList{n}.distortionModel{:} (description in @distortedImg)
%
% Optional input flags
%  'fitErr' - calculate residual shifts after fitting each image pair and
%  print the average shift length (in pixels) to the command line. This
%  gives an idea of the image registration quality,  but comes at a
%  computational cost of an extra image registration step that is not used
%  in the final output.
%  'retryMax', retryMax - max number of retry attempts by doubling the ROI
%  size to make the residual shifts smaller (default 100, which in practice
%  is infinite)
%
%
% Outputs
% job.shifts - (n,1)(m,1) nested cell array with where
%       n = (number of images in job - 1) = numel(job.resizedList)-1
%       m = number of distortionModels linked to
%           job.resizedList{n}.distortionName - find in current image
%           numel(job.resizedList{n}.distortionModel).
%           e.g. let job.resizedList{1}.distortionName = ‘tilt’;
%           so job.resizedList{1}.distortionModel = {'projective','poly11','poly22'};
%           therefore m=3.
% 	    job.shifts{n}{m} is a @struct containing these fields:
%       -	x = image x-shifts in μm, after fitting to the distortion model.
%       2D numeric array same size as image row/columns
%       -	y = same as x but for y-shifts
%       -	xshiftsROI = image x-shifts in μm, as-measured from ROI cross-correlation.
%           Numeric array or vector same number of elements as number of ROI in image.
%           If the ROI are positioned in a regular 2D grid, then this output is an array,
%           otherwise it is a vector.
%       -	yshiftsROI = same as xshiftsROI but for y-shifts.
%       -	ROI = another structure containing information about the ROI
%               size and spacing. This is only used internally by fRunDIC() and the
%               functions that it calls in funcsV0, so we don’t describe it fully
%               here.
%
%  print to the command window: mean ROI shifts [x; y; sqrt(x^2+y^2)] in
%  pixel units after each cross-correlation step
%
%
% Requirements:
%  - Uses MATLAB’s Curve Fitting Toolbox - fit()
%
% Version control
% 20240912 - option to rerun with 2x ROI size if residuals shift length is
% larger than 2 pixels
% 20240914 - 'retryMax' option bug fix -- reset ref and test properly for
% distortions with multiple distortionModels.


if ~isempty(varargin)
    %check for fitErr flag
    fitErr1 = check_option(varargin,'fitErr');
    %this means you need to calculate residual shifts

    % retry1 means you rerun the XCF with a 2x window size
    % use a random big number (100) means you retry until it works or you
    % run out of space
    % true unless this flag is on, OR if you don't calculate residual shifts (fitErr1 is false)
    retryMax = get_option(varargin,'retryMax',100,{'double';'single';'uint8';'uint16';'uint32'});
    if ~fitErr1 retryMax=0; end % if you don't calculate residual shifts you can't use this

end

% preallocate output array
job.shifts = cell(1,numel(job.resizedList)-1);

% decide where to start - nstart = smallest n with a
% high-contrast image, go to n = max, then go backwards from nstart to
% n = 1
hCList = cat(1,cat(1,job.resizedList{:}).highContrast);
nStart = find(hCList,1,'first');

if numel(job.resizedList) > nStart
    % if the only high contrast image is the last one (ground truth ref image)
    %skip the first for-loop
    %% loop through from most distorted (n=nstart) to least distorted (n=max)
    % Described in Steps 2a-c of the function header
    for n=nStart:numel(job.resizedList)-1
        retryInc = 1; % first try = 1, second try = 2 etc
        xcf1 = job.resizedList{n}.setXCF;
        while retryInc<=retryMax % this is a while loop for retry1 but test this later
            %initialise variables in for-loop
            ref = job.resizedList{n+1};
            test = job.resizedList{n};
            disp([newline 'Calculating shifts between images ' num2str(n+1) ' and ' num2str(n) ' (' test.distortionName '):']);
            imRef = ref.(ref.setXCF{1}.xcfImg);
            % if imTestNew was created at the end of the last loop iteration,
            % read imTest from that
            if exist('imTestNew','var'), imTest = imTestNew; else, imTest = test.(test.setXCF{1}.xcfImg); end; clear imTestNew


            % if no distortion (i.e. the distortion name is 'true')
            if isempty(test.distortionModel{1})
                %skip DIC and write zero shifts
                m=1;
                disp(['Mean X-shift length 0 pixels' newline 'Mean Y-shift length 0 pixels'  newline 'Mean shift length 0 pixels']);
                job.shifts{n}{m}.x = zeros(size(imRef));
                job.shifts{n}{m}.y = zeros(size(imRef));
                retryInc=Inf; %don't try to make ROI bigger because you don't use xcf anyway
            else
                %loop and stack shifts from multiple distortion models
                for m = 1:numel(test.distortionModel)
                    % divide image into ROIs, run cross-correlation on each
                    % ROI to calculate local xy shifts, output dicOut as structure array
                    dicOut = fRunDIC(imRef,imTest,xcf1{min(numel(xcf1),m)},test.distortionModel{m});

                    %unpack outputs from dicOut
                    % output is in pixels - convert to length units
                    pix2um = @(dxy,outPix) dxy*outPix;
                    job.shifts{n}{m}.x = pix2um(test.dx, dicOut.x);
                    job.shifts{n}{m}.y = pix2um(test.dy, dicOut.y);
                    job.shifts{n}{m}.xshiftsROI = pix2um(test.dx, dicOut.xshiftsROI);
                    job.shifts{n}{m}.yshiftsROI = pix2um(test.dy, dicOut.yshiftsROI);

                    % leave ROI outputs in pixels -- mainly used for debugging
                    job.shifts{n}{m}.ROI = dicOut.ROI;

                    %TODO - figure out how to include fieldnames in anonymous function
                    %to avoid repeated code

                    % update imTest to stack distortion models
                    if numel(test.distortionModel)>1 && m<numel(test.distortionModel)
                        imTest= remapImage(test.pos,job.shifts{n}{m},imTest,'test2ref');
                    end
                end
            end

            % if imRef has poor contrast (highContrast = 0), remap intensity
            % values of previous imTest onto imRef, which becomes the next
            % imTest
            if ~ref.highContrast
                imTestNew = remapImage(test.pos,job.shifts{n}{end},imTest,'test2ref');
            end

            % residual shifts - values calculated and displayed/stored but not used
            % in subsequent steps, it is a measure of correlation quality of
            % the previous image correlation step
            if fitErr1
                disp(['Residual shifts / pixels between images ' num2str(n+1) ' and ' num2str(n) ' (' test.distortionName ')']);
                job.shifts{n}{1}.fitError = fRunDIC(imRef,remapImage(test.pos,job.shifts{n}{end},imTest,'test2ref'),xcf1{end},'poly11');

                residShiftPix = mean(sqrt(abs(job.shifts{n}{1}.fitError.ROI.Shift_X_1(:).^2+abs(job.shifts{n}{1}.fitError.ROI.Shift_Y_1(:).^2))));
                if residShiftPix > 2  && retryInc<=retryMax  %if residual shift length is larger than 2 pixels and we said we want to retry
                    retryInc=retryInc+1; % go back and try the first pass again with a double sized ROI
                    disp(['Residual shifts length was greater than 2 pixels between images ' num2str(n+1) ' and ' num2str(n) ' (' test.distortionName ')']);
                    for m = 1:numel(test.distortionModel) %do this for all the xcf in this distortion model
                        newRoiSize = 2*xcf1{min(numel(xcf1),m)}.ROISize;
                        maxRoiAllowedSize = min((size(test.img,2)- xcf1{min(numel(xcf1),m)}.NumROI.x), (size(test.img,1)- xcf1{min(numel(xcf1),m)}.NumROI.y));
                        if newRoiSize < maxRoiAllowedSize
                            xcf1{min(numel(xcf1),m)}.ROISize = newRoiSize;
                            disp(['Retrying calculation with double ROI size of ' num2str(newRoiSize) ' pixels for distortionModel ' num2str(m) '...']);
                        else                         % unless it's too big
                            retryInc=Inf; % don't retry
                            disp(['Can''t retry calculation with double ROI size of ' num2str(newRoiSize) ' pixels for distortionModel ' num2str(m) ' because the image is too small.']);
                            break
                        end
                    end
                else
                    retryInc=inf; %don't retry
                end %end update ROI size for retry
            end %end fitErr - print xcf residuals

        end % end while retry
    end %end first for-loop
end %end if-check for whether or not to skip the first for-loop

%% and then loop backwards from n=nstart to n=1
% this code repeats the previous section but works backwards from the first
% high-contrast image back to the test image
% description in Steps 2e in this function header
% this is needed if the first (test) image has the highContrast property set to 0

if nStart>1
    % now work backwards
    for n=nStart-1:-1:1 %start at nStart-1, because nStart is the first ref
        %in both for-loops, n refers to the index of the test image
        retryInc = 1; % first try = 1, second try = 2 etc
        xcf1 = job.resizedList{n}.setXCF;
        while retryInc<=retryMax % this is a while loop for retry1 but test this later
             %initialise variables in for-loop
        ref = job.resizedList{n+1};
        test = job.resizedList{n};

        disp([newline 'Calculating shifts between images ' num2str(n+1) ' and ' num2str(n) ' (' test.distortionName '):']);

        imTest = test.(test.setXCF{1}.xcfImg);
        % if imRefNew was created at the end of the last loop iteration,
        % read imRef from that
        if exist('imRefNew','var'), imRef = imRefNew; else, imRef = ref.(ref.setXCF{1}.xcfImg); end; clear imRefNew


            % TODO - put this in a nested function to avoid repeating
            % code lines
            % % %
            %loop and stack shifts from multiple distortion models
            if isempty(test.distortionModel{1})
                %skip DIC and write zero shifts
                % m=1;
                disp(['Mean X-shift length 0 pixels' newline 'Mean Y-shift length 0 pixels'  newline 'Mean shift length 0 pixels']);
                job.shifts{n}{1}.x = zeros(size(imRef));
                job.shifts{n}{1}.y = zeros(size(imRef));
            else
                for m = 1:numel(test.distortionModel)
                    disp(['Distortion model ' test.distortionModel{m} ':']);
                    dicOut = fRunDIC(imRef,imTest,xcf1{min(numel(xcf1),m)},test.distortionModel{m});

                    % output is in pixels - convert to length units
                    pix2um = @(dxy,outPix) dxy*outPix;
                    job.shifts{n}{m}.x = pix2um(test.dx, dicOut.x);
                    job.shifts{n}{m}.y = pix2um(test.dy, dicOut.y);
                    job.shifts{n}{m}.xshiftsROI = pix2um(test.dx, dicOut.xshiftsROI);
                    job.shifts{n}{m}.yshiftsROI = pix2um(test.dy, dicOut.yshiftsROI);

                    % leave ROI outputs in pixels -- mainly used for debugging
                    job.shifts{n}{m}.ROI = dicOut.ROI;

                    %TODO - figure out how to include fieldnames in anonymous function
                    %to avoid repeated code

                    % update imTest to stack distortion models
                    if numel(test.distortionModel)>1 && m<numel(test.distortionModel)
                        imTest= remapImage(test.pos,job.shifts{n}{m},imTest,'test2ref');
                    end
                end
            end
            % % %

            % imTest should always have poor contrast (because of how nStart is
            % calculated, but check anyway). Remap intensity
            % values of imRef (which should always have high contrast) onto
            % imTest, and make this the next imRef
            if ~test.highContrast
                imRefNew = remapImage(ref.pos,job.shifts{n}{end},imRef,'ref2test');
            end

            % residual shifts - values calculated and displayed/stored but not used
            % in subsequent steps, it is a measure of correlation quality of
            % the previous image correlation step
            if fitErr1
                disp(['Residual shifts / pixels between images ' num2str(n+1) ' and ' num2str(n) ' (' test.distortionName ')']);
                job.shifts{n}{1}.fitError = fRunDIC(imRef,remapImage(test.pos,job.shifts{n}{end},imTest,'test2ref'),xcf1{end},'poly11');

                residShiftPix = mean(sqrt(abs(job.shifts{n}{1}.fitError.ROI.Shift_X_1(:).^2+abs(job.shifts{n}{1}.fitError.ROI.Shift_Y_1(:).^2))));
                if residShiftPix > 2  && retryInc<=retryMax  %if residual shift length is larger than 2 pixels and we said we want to retry
                    retryInc=retryInc+1; % go back and try the first pass again with a double sized ROI
                    disp(['Residual shifts length was greater than 2 pixels between images ' num2str(n+1) ' and ' num2str(n) ' (' test.distortionName ')']);
                    for m = 1:numel(test.distortionModel) %do this for all the xcf in this distortion model
                        newRoiSize = 2*xcf1{min(numel(xcf1),m)}.ROISize;
                        maxRoiAllowedSize = min((size(test.img,2)- xcf1{min(numel(xcf1),m)}.NumROI.x), (size(test.img,1)- xcf1{min(numel(xcf1),m)}.NumROI.y));
                        if newRoiSize < maxRoiAllowedSize
                            xcf1{min(numel(xcf1),m)}.ROISize = newRoiSize;
                            disp(['Retrying calculation with double ROI size of ' num2str(newRoiSize) ' pixels for distortionModel ' num2str(m) '...']);
                        else                         % unless it's too big
                            retryInc=Inf; % don't retry
                            disp(['Can''t retry calculation with double ROI size of ' num2str(newRoiSize) ' pixels for distortionModel ' num2str(m) ' because the image is too small.']);
                            break
                        end
                    end
                else
                    retryInc=inf; %don't retry
                end %end update ROI size for retry
            end %end fitErr - print xcf residuals

        end % end while retry


    end % end second for-loop (working backwards from nStart to 1)
end  %end if-check for whether or not to skip the second for-loop
%% subfunctions
    function imgOut = remapImage(posGrid,shifts,img,shiftdir)
        % Function to remap image into shifted image frame by sampling image interpolant
        % Described in Step 2d of the function header
        % on shifted XY positions between cross-correlation steps
        %
        % Inputs:
        % posGrid = @vector3d array of image pixel positions
        % shifts = job.shifts{n}{end} @struct containing xy shifts in um
        % img = r*c*z array, image values
        % shiftdir = 'test2ref' or 'ref2test' - decide which way to do the
        % remapping
        %
        % Outputs:
        % imgOut = image resampled at shifted xy positions, same size as posGrid

        switch shiftdir % decide whether to add or subtract image shifts
            case 'test2ref'
                sgn = 1; %add
            case 'ref2test'
                sgn = -1; %subtract
        end
        xRemap = posGrid.x + sgn*shifts.x;
        yRemap = posGrid.y + sgn*shifts.y;
        %relatively fast: 6s for 1.8M pixels
        newImgInterpolant = scatteredInterpolant(xRemap(:),yRemap(:),img(:),'nearest','none');
        % very slow step: 74s for 1.8M pixels
        imgOut = newImgInterpolant(posGrid.x,posGrid.y);
    end

end
