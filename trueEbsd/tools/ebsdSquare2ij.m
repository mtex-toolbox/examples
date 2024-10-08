function [ebsdImg] = ebsdSquare2ij(ebsd,ebsdMap,varargin)
% index a gridified ebsd map in the same order as image ij coordinates
% optional input - disImg.ebsdPlottingConvention


%TODO - move this out of the function - keep xy rot/flip tool separate
ebsd = argin_check(ebsd,"EBSDsquare");
if isa(ebsdMap,'char')||isa(ebsdMap,'string')
    ebsdImg = ebsd.(ebsdMap);
else
    ebsdImg = ebsdMap;
end


% by default, images are always plotted in matlab with +x = east, +y = down
% in a right handed set, this means that +z must point into the screeen

e = ebsd.plottingConvention.east;
n = ebsd.plottingConvention.north;
o = ebsd.plottingConvention.outOfScreen;

% handle each case of axis transformations explicitly
% the screen nwse must be aligned to xyz vectors.
if e==vector3d.X && n==-vector3d.Y && o==-vector3d.Z %1
    % do nothing

elseif e==vector3d.Y && n==-vector3d.X && o==vector3d.Z %2
    % rotate
    ebsdImg = rot90(ebsdImg,1);
    % % flipud/lr
    ebsdImg = flipud(ebsdImg);

elseif e==vector3d.X && n==vector3d.Y && o==vector3d.Z %3
    % rotate
    % ebsdImg = rot90(ebsdImg,);
    % flipud/lr
    ebsdImg = flipud(ebsdImg);

elseif e==vector3d.Y && n==vector3d.X && o==-vector3d.Z %4
    % rotate
    ebsdImg = rot90(ebsdImg,3);
    % flipud/lr
    % ebsdImg = flip(ebsdImg);

elseif e==-vector3d.X && n==vector3d.Y && o==-vector3d.Z %5
    % rotate
    ebsdImg = rot90(ebsdImg,2);
    % flipud/lr
    % ebsdImg = flip(ebsdImg);

elseif e==-vector3d.Y && n==vector3d.X && o==vector3d.Z %6
    % rotate
    ebsdImg = rot90(ebsdImg,1);
    % flipud/lr
    ebsdImg = fliplr(ebsdImg);

elseif e==-vector3d.X && n==-vector3d.Y && o==vector3d.Z %7
    % rotate5
    % ebsdImg = rot90(ebsdImg,);
    % flipud/lr
    ebsdImg = fliplr(ebsdImg);

elseif e==-vector3d.Y && n==-vector3d.X && o==-vector3d.Z %8
    % rotate
    ebsdImg = rot90(ebsdImg,1);
    % flipud/lr
    % ebsdImg = flip(ebsdImg);

else
    error("Can't determine screen axes for EBSD map image");
end

