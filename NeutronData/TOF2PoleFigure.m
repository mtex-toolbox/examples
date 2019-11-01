%% Convert Time of Flight Neutron Data into Pole Figures
%
% This example demonstrates how to convert SKAT time of flight neutron
% intensities into pole figure intensities
%

%% Step 1: Import and Analyze the Spectra

% import all spectra
spec = loadallspectra(fullfile(mtexExamplePath,'NeutronData','Knie_10','Knie_10'));

% compute background
bg = calc_background(spec);

% plot
plot_spectra(spec,bg)

%% Step 2: Peak Detection and Integration 

%Range fits measure period in '2005'
peakpos = [2520 1995 1480 1376 1354 1280 1206 1110]; 
peakrange = [[2509      2550];...
	[1981      2009];...
	[1470      1500];...
	[1366      1390];...
	[1344      1366];...
	[1280      1295];...
	[1196      1211];...
	[1103      1117];...
  [1015      1030]];

% specify scrystal and specimen symmetry
CS = crystalSymmetry('-3m',[4.9,4.9,5.4]);

% file with structure factors for quartz
sf = txt2mat('quartz.txt');

%
sfpos  = {1,[2 3],4,[5 6],7,8,[9,10],11,[13,14]};

for k = 1:length(sfpos)
  i = sfpos{k};
  h{k} = Miller(sf(i,2),sf(i,3),sf(i,4),sf(i,5),CS);
  c{k} = sf(i,end).';
end

 
% extract peaks and calculate spectra sums
[sumdetectr,sumphi,sumspectr,peaks,peaksbg] = proceed_spectra(spec,bg,300:1200,peakrange);


%% Step 3: Set up the Pole Figures

r = DubnaGrid(19);

pf = PoleFigure(h{1},r,fliplr(squeeze(peaks(1,:,:))),CS);

for k = 2:length(h)
  
  pf({k}) = PoleFigure(h{k}.',r,fliplr(squeeze(peaks(k,:,:))),CS);
  
end

pf.c = c;

%plot(rotate(pf,rotation('axis',xvector,'angle',90*degree)))
plot(pf)
