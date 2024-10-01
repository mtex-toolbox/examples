function [ebsd] = ebsdMapOffset(ebsd,offsetPosE, offsetPosS,varargin)
% translate EBSD object by some map direction East and South
% wrapper for EBSD/plus and EBSD/minus


% by default, images are always plotted in matlab with +x = east, +y = down
% in a right handed set, this means that +z must point into the screeen
e = ebsd.plottingConvention.east;
n = ebsd.plottingConvention.north;
o = ebsd.plottingConvention.outOfScreen;

% handle each case of axis transformations explicitly
% the screen nwse must be aligned to xyz vectors.
if e==vector3d.X && n==-vector3d.Y && o==-vector3d.Z %1
    ebsd = ebsd + offsetPosE*vector3d.X;
    ebsd = ebsd + offsetPosS*vector3d.Y;

elseif e==vector3d.Y && n==-vector3d.X && o==vector3d.Z %2
    ebsd = ebsd + offsetPosS*vector3d.X;
    ebsd = ebsd + offsetPosE*vector3d.Y;

elseif e==vector3d.X && n==vector3d.Y && o==vector3d.Z %3
    ebsd = ebsd + offsetPosE*vector3d.X;
    ebsd = ebsd - offsetPosS*vector3d.Y;

elseif e==vector3d.Y && n==vector3d.X && o==-vector3d.Z %4
    ebsd = ebsd - offsetPosS*vector3d.X;
    ebsd = ebsd + offsetPosE*vector3d.Y;

elseif e==-vector3d.X && n==vector3d.Y && o==-vector3d.Z %5
    ebsd = ebsd - offsetPosE*vector3d.X;
    ebsd = ebsd - offsetPosS*vector3d.Y;

elseif e==-vector3d.Y && n==vector3d.X && o==vector3d.Z %6
    ebsd = ebsd - offsetPosS*vector3d.X;
    ebsd = ebsd - offsetPosE*vector3d.Y;

elseif e==-vector3d.X && n==-vector3d.Y && o==vector3d.Z %7
    ebsd = ebsd - offsetPosE*vector3d.X;
    ebsd = ebsd + offsetPosS*vector3d.Y;

elseif e==-vector3d.Y && n==-vector3d.X && o==-vector3d.Z %8
    ebsd = ebsd + offsetPosS*vector3d.X;
    ebsd = ebsd - offsetPosE*vector3d.Y;

else
    error("Can't determine screen axes for EBSD map image");
end
