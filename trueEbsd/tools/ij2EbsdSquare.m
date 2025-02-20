function ebsdImg = ij2EbsdSquare(ebsd,ebsdImg,varargin)
% translate and rotate/flip gridified ebsd map to match axis ij indexing
% this is for indexing ebsd.pos in the correct order, but also works for
% any @EBSDsquare map property
% this is the reverse operation of ebsdPos2Img
% optional input - disImg.ebsdPlottingConvention


% by default, images are always plotted in matlab with +x = east, +y = down
% in a right handed set, this means that +z must point into the screeen
e = ebsd.how2plot.east;
n = ebsd.how2plot.north;
o = ebsd.how2plot.outOfScreen; 

% EBSDSquare indexing already partially compensates for 


% handle each case of axis transformations explicitly
% the screen nwse must be aligned to xyz vectors.
if e==vector3d.X && n==-vector3d.Y && o==-vector3d.Z %1
    % do nothing

elseif e==vector3d.Y && n==-vector3d.X && o==vector3d.Z %2
    
    % reverse order
    % flipud/lr
    ebsdImg = flipud(ebsdImg);
    % rotate
    ebsdImg = rot90(ebsdImg,-1);
    

elseif e==vector3d.X && n==vector3d.Y && o==vector3d.Z %3
    % rotate
    % ebsdImg = rot90(ebsdImg,);
    % flipud/lr
    ebsdImg = flipud(ebsdImg);

elseif e==vector3d.Y && n==vector3d.X && o==-vector3d.Z %4
    % rotate
    ebsdImg = rot90(ebsdImg,-3);
    % flipud/lr
    % ebsdImg = flip(ebsdImg);

elseif e==-vector3d.X && n==vector3d.Y && o==-vector3d.Z %5
    % rotate
    ebsdImg = rot90(ebsdImg,-2);
    % flipud/lr
    % ebsdImg = flip(ebsdImg);

elseif e==-vector3d.Y && n==vector3d.X && o==vector3d.Z %6
    %reverse order of operations
    % flipud/lr
    ebsdImg = fliplr(ebsdImg);
    % rotate
    ebsdImg = rot90(ebsdImg,-1);

elseif e==-vector3d.X && n==-vector3d.Y && o==vector3d.Z %7
    % rotate
    % ebsdImg = rot90(ebsdImg,);
    % flipud/lr
    ebsdImg = fliplr(ebsdImg);

elseif e==-vector3d.Y && n==-vector3d.X && o==-vector3d.Z %8
    % rotate
    ebsdImg = rot90(ebsdImg,-1);
    % flipud/lr
    % ebsdImg = flip(ebsdImg);

else
    error("Can't determine screen axes for EBSD map image");
end
