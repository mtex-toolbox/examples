function [] = saveVarsAndScript(savepname, dataName, scriptName)
% Vivian Tong NPL
% saves variables and script to disk for data backup purposes
%%%% Version Log %%%%
% V1 - Sept 2019
% V2 - 27 April 2020: 
% MAJOR BUG FIX! ln 24-26 evalin/assignin to save vars outside local workspace
% v3 - 23 Aug 2021
% list and save dependent functions as well
% v4 - 21 Sep 2021
% remove timestamps - use gitlab for version control instead
% v4 - 08 Nov 2021
% replace illegal characters in dataName before saving .mat file
% deal with windows path character limit in MAT file name
% v5 - 01 Jul 2022
% bugfix for regexp to replace illegal characters - use '\W?' instead of
% '\W*' so that it works for consecutive illegal characters too
% 18 Jul 2022
% remove m script copying, use git for version logging instead
%%%%

%% Save vraiables
%save directory
if ~exist (savepname, 'dir')
    mkdir(savepname);
end


%time now is
% timestamp = char(datestr(now,'yymmdd_HHMM'));

%save mat files (omit file extensions)
dataName2= removeweirdchars(dataName);%first remove illegal characters
saveFile = fullfile(savepname,dataName2); %mat file
if length(saveFile)>260 %windows path character limit
    saveFile(260:end)=[];
end

% copy mscript, use git for version control instead!
% %{
saveCodePath = fullfile(savepname,'src'); %m files
if ~exist (saveCodePath, 'dir')
    mkdir(saveCodePath);
end

% %list and save scripts
% [fList, ~] = matlab.codetools.requiredFilesAndProducts(scriptName);
%save main script
for n=1:numel(scriptName)
    [~,name1,ext1]=fileparts(scriptName{n});
    copyfile(scriptName{n}, fullfile(saveCodePath, [name1 ext1]));
end

close all; %close figures to avoid saving figure handles too
assignin('base','temp1',saveFile);
evalin('base','save([temp1 ''.mat''], ''-v7.3'')');