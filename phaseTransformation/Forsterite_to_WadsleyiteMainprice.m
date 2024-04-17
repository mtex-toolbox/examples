%
%
% Forsterite - Wadsleyite
%
%
% David Mainprice 23/01/2018
%

%% Specify Crystal Symmetries

% crystal symmetry
Fo_CS  = crystalSymmetry('222', [4.756 10.207 5.98], 'mineral', 'Forsterite', 'color', 'light green')
Wad_CS = crystalSymmetry('222', [5.6978 11.462 8.2571], 'mineral', 'Wadsleyite', 'color', 'light blue')


%% Define Burgers orientation relation between parent and child phases
%
% Relation (100)ol || (01-1)Wad & [001]ol || [100]Wad
% 
%
%      Parent           Child
% (100)Forsterite  || (01-1)Wadsleyite
% (010)Forsterite  ||  (012)Wadsleyite
% [001]Forsterite  || [100]Wadsleyite

Fo2Wa = orientation.map(...
  Miller(1,0,0,Fo_CS,'hkl'),Miller(0,1,-1,Wad_CS,'hkl'),...
  Miller(0,0,1,Fo_CS,'uvw'),Miller(1,0,0,Wad_CS,'uvw'))

%% axis/angle pair for Forsterite to Wadsleyite misorientation

mis_axis_wrt_Fo = round(axis(Fo2Wa))
mis_axis_wrt_Wad = round(Miller(mis_axis_wrt_Fo,Wad_CS,'hkl'))
mis_angle_Fo_Wad = angle(Fo2Wa)/degree


%% generate Child Wadsleyite orientations using misorientation 'Forsterite2Wadsleyite'

% Forsterite Parent orientation
ori_Fo_Parent = orientation.id(Fo_CS)
% compute a Wadsleyite child orientation related to the Parent Forsterite orientation
Wad_Child = ori_Fo_Parent * inv(Fo2Wa)

%%
% compute all symmetrically possible child orientations
% - strange angles like 144.231 144.231       0       0 etc
% - maybe I made a mistake here ?
% RH: this is incorrect
Wad_Child_orientations_Symm = unique(Wad_Child.symmetrise * inv(Fo2Wa))

% RH: this is correct - you have to symmetrise the parent and then apply
% the orientation relation ship to get all possibe variants
Wad_Child_orientations = unique(ori_Fo_Parent.symmetrise * inv(Fo2Wa))

% Childs using the function variants
% RH: this is correct as well and gives the same result as above
Wad_Child_orientations = ori_Fo_Parent * inv(Fo2Wa.variants)


% all possible child to child misorientations- I get the same result here !
% Wad_Child_mori = (100)/71.5371
% RH: this is correct
Wad_Child_mori = unique(Fo2Wa.variants * inv(Fo2Wa))
%
%
% Wad_Wad misorientation peak - N.B. This not an variant or child rather a twin
% Data analysis from an experiment at T = 1600C & P = 16 GPa
gB_Child_Child = orientation('axis',Miller(0,0,1,Wad_CS,'uvw'),'angle',90.0*degree,Wad_CS,Wad_CS)
%**************************************************************************
%% Parent (Forsterite) to 4 Daughter (Wadsleyite) variants 
% Mis = inv(O1) * O2

% all crystallographically equivalent orientations using symmetry of Parent
% Parent (Forsterite) reference orientation Euler (0,0,0)
ori_Fo_Parent = orientation('Euler',0,0,0,Fo_CS)
% all crystallographically equivalent Child orientations
ori_Wad_Childs = symmetrise(ori_Fo_Parent) * inv(Fo2Wa)

%% Find all combinations of ori_Wad_Childs to define misorientations for child


n=0;
nsymm = length(ori_Wad_Childs);
for i=1:nsymm
  ori_Wad_Childs(i);
  for j=1:nsymm
    if(ne(i,j))
      n=n+1;
      %fprintf('%i %i %i \n',i,j,n);
      Mis_deltaO(n) = inv(ori_Wad_Childs(i))*ori_Wad_Childs(j);
    end
  end
end
Mis_deltaO
%%
%Mis_Axes = round(axis(Mis_unique))
%Mis_Angles = angle(Mis_unique)./degree
%%
%Mis_unique = unique(Mis_deltaO)
%axis(Mis_unique(1))
%angle(Mis_unique(1))/degree
%axis(Mis_unique(2))
%angle(Mis_unique(2))/degree