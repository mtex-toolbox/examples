function [] = saveFigs(figHandle,imName,pName,varargin)
% function wrapper to save figures from matlab
% V1 - Vivian Tong NPL March 2020
% 2020-08-11
% enable optional inputs (varargin) into export_fig
% 2020-08-25
% set white figure background

if ~isfolder(pName)
   mkdir(pName);
end

saveLoc = fullfile(pName,imName);

try figHandle.Color = 'w'; end %sometimes doesn't work with mtex figures, not dealbreaker

export_fig(append(char(saveLoc), '.png'),'-dpng','-r300','-a1',varargin{:},figHandle);
% savefig(figHandle,saveLoc);


