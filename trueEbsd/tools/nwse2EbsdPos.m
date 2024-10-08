function [xEbsdNew, yEbsdNew] = nwse2EbsdPos(ebsd,imgEast,imgSouth,varargin)
% transform image ij map positions (NWSE) into ebsd.pos (XYZ) 
% this is nothing to do with the order of matrix indexing, but what the actual
% numbers are, so imgEast and imgSouth can be any size (inc column vectors)

% by default, images are always plotted in matlab with +x = east, +y = down
% in a right handed set, this means that +z must point into the screeen
e = ebsd.plottingConvention.east;
n = ebsd.plottingConvention.north;
o = ebsd.plottingConvention.outOfScreen; 

% handle each case of axis transformations explicitly
% the screen nwse must be aligned to xyz vectors.
if e==vector3d.X && n==-vector3d.Y && o==-vector3d.Z %1
    % do nothing
    [xEbsdNew, yEbsdNew] = deal(imgEast,imgSouth);   

elseif e==vector3d.Y && n==-vector3d.X && o==vector3d.Z %2
    [xEbsdNew, yEbsdNew] = deal(imgSouth,imgEast);

elseif e==vector3d.X && n==vector3d.Y && o==vector3d.Z %3
    [xEbsdNew, yEbsdNew] = deal(imgEast,-imgSouth);

elseif e==vector3d.Y && n==vector3d.X && o==-vector3d.Z %4
    [xEbsdNew, yEbsdNew] = deal(-imgSouth,imgEast);

elseif e==-vector3d.X && n==vector3d.Y && o==-vector3d.Z %5
    [xEbsdNew, yEbsdNew] = deal(-imgEast,-imgSouth);

elseif e==-vector3d.Y && n==vector3d.X && o==vector3d.Z %6
    [xEbsdNew, yEbsdNew] = deal(-imgSouth,-imgEast);

elseif e==-vector3d.X && n==-vector3d.Y && o==vector3d.Z %7
    [xEbsdNew, yEbsdNew] = deal(-imgEast,imgSouth);

elseif e==-vector3d.Y && n==-vector3d.X && o==-vector3d.Z %8
    [xEbsdNew, yEbsdNew] = deal(imgSouth,-imgEast);

else
    error("Can't determine screen axes for EBSD map image");
end

%}