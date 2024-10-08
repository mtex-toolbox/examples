function [stringOut] =  removeweirdchars(stringIn)
%remove non AZ09 characters for saving file names 
stringOut=stringIn;
stringOut(regexp(stringIn,'\W?','start')) = ''; 