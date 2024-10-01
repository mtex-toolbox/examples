
function [Shift_X,Shift_Y,CCmax] = fDIC_xcf_mat(Image_ref,Image_test,ROI,Filters_setting,XCF_mesh,hfilter,FFTfilter)
% cross correlation is performed and shift X, shift Y and peak height are
% determined.
%#codegen
% Image_size = size(Image_ref);
Shift_X_temp=zeros(ROI.num_pass_1,1);
Shift_Y_temp=zeros(ROI.num_pass_1,1);
CCmax=zeros(ROI.num_pass_1,1);

parfor j=1:ROI.num_pass_1
        ROI_ref_1 = Image_ref((ROI.coordinator_pass_1{j}(1,2)-round(ROI.size_pass_1/2)):(ROI.coordinator_pass_1{j}(1,2)+round(ROI.size_pass_1/2)-1),(ROI.coordinator_pass_1{j}(1,1)-round(ROI.size_pass_1/2)):(ROI.coordinator_pass_1{j}(1,1)+round(ROI.size_pass_1/2)-1));
        % zero mean and normalise standard deviation
        ROI_ref_1 = (ROI_ref_1- mean(ROI_ref_1(:)))./std(ROI_ref_1(:));
        ROI_ref_1 = ROI_ref_1.*hfilter; % han filtering
        ROI_ref_1 = fft2(ROI_ref_1); % 2D fast fourier transform
        data_fill =[1:(Filters_setting(3)+Filters_setting(4)),ROI.size_pass_1-(Filters_setting(3)+Filters_setting(4)-1):ROI.size_pass_1];
        data_fill = data_fill';
        ROI_ref = FFTfilter(data_fill,data_fill).*ROI_ref_1(data_fill,data_fill); % apply high and low frequence filter
        
        ROI_test_1 = Image_test((ROI.coordinator_pass_1{j}(1,2)-round(ROI.size_pass_1/2)):(ROI.coordinator_pass_1{j}(1,2)+round(ROI.size_pass_1/2)-1),(ROI.coordinator_pass_1{j}(1,1)-round(ROI.size_pass_1/2)):(ROI.coordinator_pass_1{j}(1,1)+round(ROI.size_pass_1/2)-1));
        % zero mean and normalise standard deviation
        ROI_test_1 = (ROI_test_1- mean(ROI_test_1(:)))./std(ROI_test_1(:));
        ROI_test_1 = ROI_test_1.*hfilter; % han filtering
        ROI_test_1 = fft2(ROI_test_1); % 2D fast fourier transform
        ROI_test = FFTfilter(data_fill,data_fill).*ROI_test_1(data_fill,data_fill);

        [Shift_X_temp(j),Shift_Y_temp(j),CCmax(j)] = fReg(ROI_ref,ROI_test, ROI.size_pass_1,XCF_mesh, data_fill);
end
Shift_X = reshape(Shift_X_temp,size(ROI.position_X_pass_1));
Shift_Y = reshape(Shift_Y_temp,size(ROI.position_X_pass_1));

