function t = plot(job, varargin)

if nargin > 1 && iscell(varargin{1})
  imgList = varargin{1};
  varargin(1) = [];
else
  imgList = job.imgList;
end

figure('WindowState', 'maximized'); 
t = tiledlayout('flow','TileSpacing','tight','Padding','tight');
for n = 1:numel(imgList)
  nexttile;
  plot(imgList{n})   
end
linkaxes;

end