classdef trueEbsd 
    % @trueEbsd class constructor 
    % class guiding through the trueEBSD image correction process
    %
    % Constructor syntax
    % job = trueEbsd(disImgTEST, disImg2, disImg3, disImgREF)
    % job = trueEbsd(disImgEBSD, disImgNF, disImgFF, trueImg)    %
    % 
    % Input
    % cell array of @distortedImg in order from most distorted {1} to reference image {end}
    %  imgList{:} - cell array of @distortedImg objects in the correct order for image
    %  registration
    %  starting with most distorted imgList{1} to ground truth imgList{end}
    %
    % Output
    % job = @trueEbsd object
    %
    % Class Properties
    % imgList – cell array of input @distortedImg
    % resizedList – same size cell array of empty @distortedImg, will be used to store output from pixelSizeMatch()
    % undistortedList – same size cell array of empty @distortedImg, will be used to store output from undistort()
    % shifts – empty cell, will be used to store output from calcShifts()
    %
    % See also
    % distortedImg parentGrainReconstructor 
    %

    properties %props we define at the start
        imgList = {distortedImg} %cell array of distortedImg objects        
        resizedList = {distortedImg} %pixel sizes and FOV matched
        undistortedList = {distortedImg} %undistorted
        shifts = {}
    end %properties

    methods
        % constructor
        function job = trueEbsd(varargin)
            nImgs=numel(varargin);
            job.imgList=cell(nImgs,1);    
             for n=1:nImgs
                % check they're the correct type
                job.imgList{n} = argin_check(varargin{n},{'distortedImg'});       
            end
        end % constructor function

    end %methods

end %classdef