% Define dislocation types for olivine
% Sept. 23: add [010] slip systems
% See list in Wallis et al., Ultramicroscopy (2016), Table 1
function[dSo]= DisTypesFo(CS)
%CS = ebsd('Forsterite').CS;
sSo1 = slipSystem(Miller(1,0,0,CS,'uvw'),Miller(0,1,0,CS,'hkl'));
sSo2 = slipSystem(Miller(1,0,0,CS,'uvw'),Miller(0,0,1,CS,'hkl'));
sSo3 = slipSystem(Miller(0,0,1,CS,'uvw'),Miller(1,0,0,CS,'hkl'));
sSo4 = slipSystem(Miller(0,0,1,CS,'uvw'),Miller(0,1,0,CS,'hkl'));
% [010]
sSo5 = slipSystem(Miller(0,1,0,CS,'uvw'),Miller(1,0,0,CS,'hkl')); 
sSo6 = slipSystem(Miller(0,1,0,CS,'uvw'),Miller(0,0,1,CS,'hkl')); 
% Screw dislocations are automatically calculated by command dislocationSystem
%assemble
sSo = cat(5,sSo1,sSo2,sSo3,sSo4,sSo5,sSo6);
% to create dislocation systems need to symmetrise slip system
sSs=sSo.symmetrise('antipodal');
% calculate dislocation system
dSo=dislocationSystem(sSs);
% from Heinisch et al. Table 2
dSo.u = [0.73; 0.74; 1.0; 0.88; 2.36; 2.64; 0.65; 1.91; 0.46]; % S: Heinisch
%clear sSo6 sSo5 sSo4 sSo3 sSo2 sSo1 sSo sSs
% for calculation of dislocation density need to modify 
% mtexVersion/geometry/@vector3d/rotate.m