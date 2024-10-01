function [RegOutput]= fRunDIC(Image_ref,Image_test,setXCF,fitfunc)
%% Unpack inputs
ROI.size_pass_1 =setXCF.ROISize;
numROI = setXCF.NumROI;
XCF_mesh = setXCF.XCFMesh; % XCF mesh size default value is 250

%% run  main functions
%set up ROIs
% set up the filter windows, boundary
% filters = round([log2(ROI.size_pass_1)/2,log2(ROI.size_pass_1)/4,2*log2(ROI.size_pass_1),log2(ROI.size_pass_1)]);
Filters_setting = [0;0;round(ROI.size_pass_1/2);round(ROI.size_pass_1/4)];
[rowvals,colvals]=find(~(isnan(Image_ref)).*(~isnan(Image_test))==1);
xmin=min(colvals);
ymin=min(rowvals);
xmax=max(colvals);
ymax=max(rowvals);
boundary = [xmin+1,ymin+1,xmax-xmin-1,ymax-ymin-1]; %[left edge top edge width height], like imcrop.
[FFTfilter,hfilter] = fFilters(ROI.size_pass_1,Filters_setting);
% spread out the subregions in the defined ROI and set up the position of ROIs
[ROI.position_X_pass_1, ROI.position_Y_pass_1,ROI.num_x_pass_1,ROI.num_y_pass_1, ROI.coordinator_pass_1, ROI.num_pass_1] = fDIC_ROI_position(ROI.size_pass_1,numROI,boundary);

% perform the XCF and determine shift in x  shift in y and peak height
[ROI.Shift_X_1,ROI.Shift_Y_1,CCmax_1] = fDIC_xcf_mat_mex(Image_ref,Image_test,ROI,Filters_setting,XCF_mesh,hfilter,FFTfilter);

%% fit a surface to the shift vector
% figure, imagesc(ROI.Shift_X_1); axis equal; caxis([-15 -5]);
% figure, imagesc(ROI.Shift_Y_1); axis equal; caxis([-5 5]);

% ROI.position_X_pass_1 and ROI.position_Y_pass_1 are the ROI centre locations
switch fitfunc
    case 'poly11'
    %remove values with NaN XCF height
    ROI.Shift_X_1=ROI.Shift_X_1(~isnan(CCmax_1));
    ROI.Shift_Y_1=ROI.Shift_Y_1(~isnan(CCmax_1));
    ROI.position_X_pass_1=ROI.position_X_pass_1(~isnan(CCmax_1));
    ROI.position_Y_pass_1=ROI.position_Y_pass_1(~isnan(CCmax_1));
    CCmax=CCmax_1(~isnan(CCmax_1));

    % calculate best fit
    [xsurf, ~]=fit([ROI.position_X_pass_1(:),ROI.position_Y_pass_1(:)],ROI.Shift_X_1(:),fitfunc,'Robust','Bisquare','Weights',CCmax);
    [ysurf, ~]=fit([ROI.position_X_pass_1(:),ROI.position_Y_pass_1(:)],ROI.Shift_Y_1(:),fitfunc,'Robust','Bisquare','Weights',CCmax);
    RegOutput.transMat = [xsurf.p10+1,xsurf.p01,xsurf.p00;
        ysurf.p10,ysurf.p01+1,ysurf.p00;
        0,0,1];
    
    %for quiver plot
    xshiftsROI=xsurf(ROI.position_X_pass_1,ROI.position_Y_pass_1);
    yshiftsROI=ysurf(ROI.position_X_pass_1,ROI.position_Y_pass_1);
    
    % for entire image
    [Xpixels,Ypixels]=meshgrid(1:size(Image_ref,2),1:size(Image_ref,1)); %pixel grid in reference frame
    xshifts=xsurf(Xpixels,Ypixels); %pixel shifts in reference frame
    yshifts=ysurf(Xpixels,Ypixels);
    
    case 'poly22'
    %remove values with NaN XCF height
    ROI.Shift_X_1=ROI.Shift_X_1(~isnan(CCmax_1));
    ROI.Shift_Y_1=ROI.Shift_Y_1(~isnan(CCmax_1));
    ROI.position_X_pass_1=ROI.position_X_pass_1(~isnan(CCmax_1));
    ROI.position_Y_pass_1=ROI.position_Y_pass_1(~isnan(CCmax_1));
    CCmax=CCmax_1(~isnan(CCmax_1));

    % calculate best fit
    [xsurf, ~]=fit([ROI.position_X_pass_1(:),ROI.position_Y_pass_1(:)],ROI.Shift_X_1(:),fitfunc,'Robust','Bisquare','Weights',CCmax);
    [ysurf, ~]=fit([ROI.position_X_pass_1(:),ROI.position_Y_pass_1(:)],ROI.Shift_Y_1(:),fitfunc,'Robust','Bisquare','Weights',CCmax);
    
    
    %for quiver plot
    xshiftsROI=xsurf(ROI.position_X_pass_1,ROI.position_Y_pass_1);
    yshiftsROI=ysurf(ROI.position_X_pass_1,ROI.position_Y_pass_1);
    
    % for entire image
    [Xpixels,Ypixels]=meshgrid(1:size(Image_ref,2),1:size(Image_ref,1)); %pixel grid in reference frame
    xshifts=xsurf(Xpixels,Ypixels); %pixel shifts in reference frame
    yshifts=ysurf(Xpixels,Ypixels);
    
    
    case 'linearinterp'
    CCmax=reshape(CCmax_1,size(ROI.coordinator_pass_1));
    shiftX_temp=nan(size(CCmax));
    shiftY_temp=shiftX_temp;
    shiftX_temp(~isnan(CCmax))=ROI.Shift_X_1(~isnan(CCmax));
    shiftY_temp(~isnan(CCmax))=ROI.Shift_Y_1(~isnan(CCmax));
    
    % calculate best fit
    shiftX_rows=median(shiftX_temp,2,'omitnan');
    shiftY_rows=median(shiftY_temp,2,'omitnan');
    
    nancheck=isnan(shiftX_rows) | isnan(shiftY_rows);
    
    ROI.position_Y_pass_1(nancheck,:)=[];
    ROI.position_X_pass_1(nancheck,:)=[];
    shiftX_rows(nancheck)=[];
    shiftY_rows(nancheck)=[];
    ROI.Shift_X_1(nancheck,:)=[];
    ROI.Shift_Y_1(nancheck,:)=[];

    
    
    [xfit, ~]=fit(ROI.position_Y_pass_1(:,1),shiftX_rows(:),fitfunc);
    [yfit, ~]=fit(ROI.position_Y_pass_1(:,1),shiftY_rows(:),fitfunc);
    
    %for quiver plot
    xshifts_lineROI=xfit(ROI.position_Y_pass_1(:,1));
    yshifts_lineROI= yfit(ROI.position_Y_pass_1(:,1));
    xshiftsROI=repmat(xshifts_lineROI(:),1,size(ROI.position_Y_pass_1,2));
    yshiftsROI=repmat(yshifts_lineROI(:),1,size(ROI.position_Y_pass_1,2));
    
    
    % for entire image
%     [Xpixels,Ypixels]=meshgrid(1:size(Image_ref,2),1:size(Image_ref,1)); %pixel grid in reference frame
    Ypixels_line=1:size(Image_ref,1);
    xshifts_line=xfit(Ypixels_line);
    yshifts_line= yfit(Ypixels_line);
    
    xshifts = repmat(xshifts_line(:),1,size(Image_ref,2));
    yshifts = repmat(yshifts_line(:),1,size(Image_ref,2));
    
    
    case 'projective'
    %remove values with NaN XCF height
    ROI.Shift_X_1=ROI.Shift_X_1(~isnan(CCmax_1));
    ROI.Shift_Y_1=ROI.Shift_Y_1(~isnan(CCmax_1));
%     CCmax=CCmax_1(~isnan(CCmax_1));
    ROI.position_X_pass_1=ROI.position_X_pass_1(~isnan(CCmax_1));
    ROI.position_Y_pass_1=ROI.position_Y_pass_1(~isnan(CCmax_1));
    
    % calculate best fit
    zeromat=zeros(numel(ROI.position_X_pass_1),1);
    onesmat=zeromat+1;
    
    lhs_mat=[ROI.position_X_pass_1(:), ROI.position_Y_pass_1(:),onesmat(:), ...
        zeromat(:), zeromat(:), zeromat(:), ...
        -(ROI.position_X_pass_1(:).*(ROI.position_X_pass_1(:)-ROI.Shift_X_1(:))) ...
        -(ROI.position_Y_pass_1(:).*(ROI.position_X_pass_1(:)-ROI.Shift_X_1(:)));
        
        zeromat(:), zeromat(:), zeromat(:),...
        ROI.position_X_pass_1(:), ROI.position_Y_pass_1(:), onesmat(:), ...
        -(ROI.position_X_pass_1(:).*(ROI.position_Y_pass_1(:)-ROI.Shift_Y_1(:))) ...
        -(ROI.position_Y_pass_1(:).*(ROI.position_Y_pass_1(:)-ROI.Shift_Y_1(:)))];
    
    rhs_mat=[(ROI.position_X_pass_1(:)-ROI.Shift_X_1(:));
        (ROI.position_Y_pass_1(:)-ROI.Shift_Y_1(:))];
    
    [rottrans_mat,~] = robustfit(lhs_mat,rhs_mat,'bisquare',4.685,0); %like wmd solution
    rottrans_mat=reshape([rottrans_mat(:);1],3,3)';
    
    %for quiver plot
    onesmatgrid=ones(size(ROI.position_X_pass_1));
    transgridROI=rottrans_mat*([ROI.position_X_pass_1(:) ROI.position_Y_pass_1(:) onesmatgrid(:)])';
    transgridxROI=reshape(transgridROI(1,:),size(ROI.position_X_pass_1,1),size(ROI.position_X_pass_1,2));
    transgridyROI=reshape(transgridROI(2,:),size(ROI.position_X_pass_1,1),size(ROI.position_X_pass_1,2));
    xshiftsROI=ROI.position_X_pass_1-transgridxROI;
    yshiftsROI=ROI.position_Y_pass_1-transgridyROI;
    
    % for entire image
    [Xpixels,Ypixels]=meshgrid(1:size(Image_ref,2),1:size(Image_ref,1)); %Data grid
    onesmatgrid=ones(size(Image_ref));
    transgrid=rottrans_mat*([Xpixels(:) Ypixels(:) onesmatgrid(:)])';
    transgridx=reshape(transgrid(1,:),size(Image_ref,1),size(Image_ref,2));
    transgridy=reshape(transgrid(2,:),size(Image_ref,1),size(Image_ref,2));
    xshifts=Xpixels-transgridx;
    yshifts=Ypixels-transgridy;


    case 'interpolate'

    %remove values with NaN XCF height
    ROI.Shift_X_1=ROI.Shift_X_1(~isnan(CCmax_1));
    ROI.Shift_Y_1=ROI.Shift_Y_1(~isnan(CCmax_1));
    ROI.position_X_pass_1=ROI.position_X_pass_1(~isnan(CCmax_1));
    ROI.position_Y_pass_1=ROI.position_Y_pass_1(~isnan(CCmax_1));
    CCmax=CCmax_1(~isnan(CCmax_1));

    % normalise data to avoid scaling artefacts in scatteredInterpolant
    [xPosNorm,cXP,sXP] = normalize(ROI.position_X_pass_1(:));
    [yPosNorm,cYP,sYP] = normalize(ROI.position_Y_pass_1(:));
    [xShiftNorm,cXS,sXS] = normalize(ROI.Shift_X_1(:));
    [yShiftNorm,cYS,sYS] = normalize(ROI.Shift_Y_1(:));
    % scale/centring values to undo the normalisation later:
    % N = (A - C) ./ S
    % A = N .* S + C

    % create scattered interpolant object
    xsurf = scatteredInterpolant([xPosNorm,yPosNorm],xShiftNorm,'natural','linear');
    % no need to compute this twice at the same grid points -- copy and
    % overwrite with Shift_X_1 values instead
    % % ysurf = scatteredInterpolant([ROI.position_X_pass_1(:),ROI.position_Y_pass_1(:)],ROI.Shift_Y_1(:),'natural','linear');%don't do this
    ysurf = xsurf; ysurf.Values=yShiftNorm;

    %for quiver plot
    xshiftsROI=(xsurf(xPosNorm,yPosNorm) .* sXS) + cXS;
    yshiftsROI=(ysurf(xPosNorm,yPosNorm) .* sYS) + cYS;
    
    % for entire image
     %pixel grid in reference frame
    [Xpixels,Ypixels]=meshgrid(1:size(Image_ref,2),1:size(Image_ref,1));
    % normalise by scaling/centring values for x and y positions
    % before sampling interpolant object
    %pixel shifts in reference frame
    xshifts=(xsurf(...
        (Xpixels - cXP)./sXP,(Ypixels-cYP)./sYP)...
        .* sXS) + cXS;
    yshifts=(ysurf(...
        (Xpixels-cXP)./sXP,(Ypixels-cYP)./sYP)...
        .* sYS) + cYS;
    % undo normalisation to get this back in pixel units
    
    otherwise
    %remove values with NaN XCF height
    ROI.Shift_X_1=ROI.Shift_X_1(~isnan(CCmax_1));
    ROI.Shift_Y_1=ROI.Shift_Y_1(~isnan(CCmax_1));
%     CCmax=CCmax_1(~isnan(CCmax_1));
    ROI.position_X_pass_1=ROI.position_X_pass_1(~isnan(CCmax_1));
    ROI.position_Y_pass_1=ROI.position_Y_pass_1(~isnan(CCmax_1));
    
    % calculate best fit
    [xsurf, ~]=fit([ROI.position_X_pass_1(:),ROI.position_Y_pass_1(:)],ROI.Shift_X_1(:),fitfunc);
    [ysurf, ~]=fit([ROI.position_X_pass_1(:),ROI.position_Y_pass_1(:)],ROI.Shift_Y_1(:),fitfunc);
    
    %for quiver plot
    xshiftsROI=xsurf(ROI.position_X_pass_1,ROI.position_Y_pass_1);
    yshiftsROI=ysurf(ROI.position_X_pass_1,ROI.position_Y_pass_1);
    
    % for entire image
    [Xpixels,Ypixels]=meshgrid(1:size(Image_ref,2),1:size(Image_ref,1)); %pixel grid in reference frame
    xshifts=xsurf(Xpixels,Ypixels); %pixel shifts in reference frame
    yshifts=ysurf(Xpixels,Ypixels);
    
end

%% Pack outputs
disp(['Mean X-shift length ' num2str(mean(abs(ROI.Shift_X_1(:)))) ' pixels']);
disp(['Mean Y-shift length ' num2str(mean(abs(ROI.Shift_Y_1(:)))) ' pixels']);
disp(['Mean shift length ' num2str(mean(sqrt(abs(ROI.Shift_X_1(:).^2+abs(ROI.Shift_Y_1(:).^2))))) ' pixels']);

RegOutput.x= xshifts;
RegOutput.y = yshifts;
RegOutput.xshiftsROI= xshiftsROI;
RegOutput.yshiftsROI = yshiftsROI;
RegOutput.ROI = ROI;


