function diff_totmag = GenBoundaryMap(testimg_noisy1,varargin)

if nargin==1
    padwidth=1;
else
    padwidth= varargin{1};
end

%% import image
%resize image
% testimg_noisy_rs=imresize(testimg_noisy,[rpts cpts],'nearest');

cpts=size(testimg_noisy1,2); %columns
rpts=size(testimg_noisy1,1); %rows

% make pretend-RGB 3 layer stack for greyscale images
if length(size(testimg_noisy1))== 2
    testimg_noisy1=repmat(testimg_noisy1,[1 1 3]);
end

%% normalise image contrast
testimg_noisy=zeros(size(testimg_noisy1));
for c=1:3
    temp=testimg_noisy1(:,:,c);
    valmax=prctile(temp(:),98);
    valmin=prctile(temp(:),2);
    testimg_noisy(:,:,c)=(temp-valmin)/(valmax-valmin);
end

%% RGB difference vectors calculation - precalcs
%Pad data with NaN
% possibly pad with more than 1 row/col of NaN (grain size / pixel size)


testimg_R=zeros(rpts+padwidth*2,cpts+padwidth*2,size(testimg_noisy,3)); %[padded rdata, padded cdata, RGB data]
testimg_L=testimg_R;
testimg_B=testimg_R;
testimg_T=testimg_R;
testimg_RT=testimg_R;
testimg_LT=testimg_R;
testimg_RB=testimg_R;
testimg_LB=testimg_R;
testimg_noisy_pad=testimg_R;

testimg_R(1+padwidth:end-padwidth,1+padwidth*2:end,:)=testimg_noisy; %right
testimg_L(1+padwidth:end-padwidth,1:end-padwidth*2,:)=testimg_noisy; %shifted to left
testimg_B(1+padwidth*2:end,1+padwidth:end-padwidth,:)=testimg_noisy; %bottom
testimg_T(1:end-padwidth*2,1+padwidth:end-padwidth,:)=testimg_noisy; %shifted to top

testimg_RT(1:end-padwidth*2,1+padwidth*2:end,:)=testimg_noisy;
testimg_LT(1:end-padwidth*2,1:end-padwidth*2,:)=testimg_noisy; 
testimg_RB(1+padwidth*2:end,1+padwidth*2:end,:)=testimg_noisy;
testimg_LB(1+padwidth*2:end,1:end-padwidth*2,:)=testimg_noisy; 

testimg_noisy_pad(1+padwidth:end-padwidth,1+padwidth:end-padwidth,:)=testimg_noisy; %not shifted


%% RGB difference vectors calculation - maps

diff_L=testimg_L-testimg_noisy_pad; %difference
diff_R=testimg_R-testimg_noisy_pad; %difference
diff_B=testimg_B-testimg_noisy_pad; %difference
diff_T=testimg_T-testimg_noisy_pad; %difference
diff_RT=testimg_RT-testimg_noisy_pad; %difference
diff_LT=testimg_LT-testimg_noisy_pad; %difference
diff_RB=testimg_RB-testimg_noisy_pad; %difference
diff_LB=testimg_LB-testimg_noisy_pad; %difference

diff_L_mag=sqrt(diff_L(:,:,1).^2+diff_L(:,:,2).^2+diff_L(:,:,3).^2);
diff_R_mag=sqrt(diff_R(:,:,1).^2+diff_R(:,:,2).^2+diff_R(:,:,3).^2);
diff_B_mag=sqrt(diff_B(:,:,1).^2+diff_B(:,:,2).^2+diff_B(:,:,3).^2);
diff_T_mag=sqrt(diff_T(:,:,1).^2+diff_T(:,:,2).^2+diff_T(:,:,3).^2);
diff_RT_mag=sqrt(diff_RT(:,:,1).^2+diff_RT(:,:,2).^2+diff_RT(:,:,3).^2);
diff_LT_mag=sqrt(diff_LT(:,:,1).^2+diff_LT(:,:,2).^2+diff_LT(:,:,3).^2);
diff_RB_mag=sqrt(diff_RB(:,:,1).^2+diff_RB(:,:,2).^2+diff_RB(:,:,3).^2);
diff_LB_mag=sqrt(diff_LB(:,:,1).^2+diff_LB(:,:,2).^2+diff_LB(:,:,3).^2);

diff_totmag1=diff_L_mag+diff_R_mag+diff_T_mag+diff_B_mag+diff_RT_mag+diff_LT_mag+diff_RB_mag+diff_LB_mag;
diff_totmag=diff_totmag1(1+padwidth:end-padwidth,1+padwidth:end-padwidth);


