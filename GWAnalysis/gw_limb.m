function [OutData,PW] = gwanalyse_limb(Data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%general-purpose function to measure GWs using vertical profile data
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/10/29
%
%updates:
%  2023/11/03: improved calculation of background temperature; output variable 'Temp' removed and replaced with 'BG'.
%  2023/11/03: added calculation of GW potential energy - new output 'Ep'
%  2024/05/01: added capability to detect multiple waves per spectrum (to Analysis=1 only - hope to add to 2 later)
%
%===================================
%OUTPUTS:
%
%   OutData: struct containing the output fields, on a profile x height (x npeaks, if requested) grid
%     FailReason - reason why no data present (0: data present; 1: MaxdX or Maxdt exceeded; 2: MinFracInProf not met ; 3: MindPhi exceeded; 4. input data was NaN; 5. no peak found )
%     A   - amplitude (input units)
%     Lz  - vertical wavelength (km)
%     Lh  - horizontal wavelength (km)
%     Lat - latitude (deg)
%     Lon - longitude (deg)
%     Alt - altitude (input units)
%     Tp  - temperature perturbation (input units)
%     MF  - momentum flux (mPa)
%     BG  - background temperature (input units)
%     Ep  - potential energy (J.kg^-1)
%     MF and Lh will be present but all NaN if 'filter' is set to anything other than 2
%
%   PWs:     array containing PWs. Format varies depending on filter type used.
%
%===================================
%
%
%INPUTS:
%
%  required:
%    InstrumentData [struct] - data to use for the analysis.  This must contain at least the following:
%         Lat:          latitude of each point, in degrees
%         Lon:          longitude of each point, in degrees
%         Alt:          altitude of each point, units arbitrary but must be internally consistent
%         Time:         time, in Matlab datenum() units
%
%   For all filters except 'Hindley23', we also need:
%         Temp:         temperature, units arbitrary but must be internally consistent
%
%   If using the 'Hindley23' planetary wave filter, which is computed externally, then we instead need:
%         Temp_Residual: residual temperatures after the PW filtering, also arbitrary-but-consistent
%         Temp_PW: planetary wave temperatures for each mode at each location, same units
%      You probably don't want to do this, but we can also use these inputs for the other filters instead of 'Temp'- the
%      routine will merge them to produce 'Temp' if so and then proceed normally.
%
%
%  optional:
%
%    VarName         (type,           default)  description
%    -----------------------------------------------------------------------------
%    Analysis        (integer,              2)  type of analysis to use, from (1) 1DST, (2) Alexander et al 2008 
%    Filter          (char,          'SGolay')  type of detrending filter to use (see below)
%    STScales        (vector,   1:1:NLevels/2)  number of scales to use in 1D ST.
%    STc             (positive real,        1)  value of 'c' to use in ST
%    MinLz           (real,                 0)  minimum vertical wavelength returned by ST
%    MaxLz           (real,             99e99)  maximum vertical wavelength returned by ST
%    STPadSize       (real,                20)  levels of zero-padding to put at each end of the data before S-Transforming
%    RegulariseZ     (logical            true)  interpolate the data to a regular height grid
%    Verbose         (logical,           true)  report to the user what's happening
%    N               (real,              0.02)  assumed Brunt-Vaisala frequency
%    g               (real,              9.81)  assumed acceleration due to gravity
%    FullST          (logical,          false)  return the 2D complex ST for each profile
%    NPeaks          (positive integer,     1)  find this many multiple peaks in spectra
%    NoDistDiscard   (logical,          false)  don't remove profile pairs which are too far apart in space or time for pairing
%    PadandTaper     (logical,          false)  pad the data externally and apply a taper within this range
%
%-----------------------------------
%if 'PadandTaper' is set to true, then the following options can be used:
%    TaperLimits     (positive real,  [0,100])  vertical range of data AFTER padding . If data is longer than this at either end, the *default* will switch to the true value
%    TaperLength     (positive real,        5)  length of taper in km. Will be linear from 0 to 1 over this range, and apply outside the original data range
%-----------------------------------
%if 'Analysis' is set to 2, then the following options can be used:
%    MindPhi         (positive real,    0.025)  minimum fractional phase change permitted for Alexander et al 2008 analysis 
%    MaxdX           (positive real,      400)  maximum interprofile distance [km] permitted for Alexander et al 2008 analysis 
%    Maxdt           (positive real,      900)  maximum interprofile time [s] permitted for Alexander et al 2008 analysis 
%    MinFracInProf   (positive real,      0.5)  minimum fraction useful levels remaining in a profile after above two filters to proceed
%-----------------------------------
%if 'Filter' is set to 'Hindley23', then the following options can be used:
%    H23_OutRem (logical,                true)   remove outliers from H23 fitting, as in get_limbsounders()
%-----------------------------------
%if 'Filter' is set to 'SGolay', then the following options can be used:
%    SGOrder         (integer,              2)  order of Savitzky-Golay filter to use
%    SGLength        (real,                18)  frame length of SG filter, in km
%When using this filter singularities like the tropopause must be removed in advance
%-----------------------------------
%if 'Filter' is set to 'PWgrid', then the following options can be used:
%    NPWs            (integer,              3)  number of PWs to fit (plus zonal mean)
%    PWWindow        (positive real,        1)  number of days to use in each PW fit
%    PWTimeRes       (positive real,        1)  time resolution to export PW fits on
%    PWLonGrid       (vector      -180:20:180)  longitude grid to compute PWs on
%    PWMinPC         (positive real,       66)  fraction of longitude bins which must be filled for each lat band
%    PWLatGrid       (vector         -90:5:90)  longitude grid to compute PWs on
%    PWAltGrid       (vector, input grid mean)  altitude grid to compute PWs on
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%


%data structure
addRequired(p,'Data',@isstruct); %input data must be a struct. Will be hand-parsed later.

%analysis to use
addParameter(p,'Analysis',2,@(x) validateattributes(x,{'numeric'},{'integer','<=',2})) %type of analysis to use (see above)

%filter to use
addParameter(p,'Filter','SGolay',@ischar) %type of filter to use

%ST properties
addParameter(p,'STScales',1:1:size(Data.Alt,2)/2,@isvector)   %scales to compute on ST
addParameter(p,'STc',                          1,@(x) validateattributes(x,{'numeric'},{'>',0})) %'c' parameter for ST
addParameter(p,'STPadSize',                   20,@isnumeric ) %levels of zero-padding to put at each end of the data before S-Transforming
addParameter(p,'MinLz',                        0,@isnumeric)  %minimum vertical wavelength returned
addParameter(p,'MaxLz',                    99e99,@isnumeric)  %maximum vertical wavelength returned
addParameter(p,'FullST',                   false,@islogical)  %return full ST obejct for each prifle


%Alex08 horizontal wavelength properties
addParameter(p,'MaxdX',           300,@(x) validateattributes(x,{'numeric'},{'>',0})) %maximum distance between profiles
addParameter(p,'Maxdt',           900,@(x) validateattributes(x,{'numeric'},{'>',0})) %maximum time between profiles
addParameter(p,'MinFracInProf',   0.5,@(x) validateattributes(x,{'numeric'},{'>',0})) %maximum fraction of profile remaining after above filters
addParameter(p,'MindPhi',       0.025,@(x) validateattributes(x,{'numeric'},{'>',0})) %minimum fractional phase change to compute Lh
addParameter(p,'NoDistDiscard', false,@islogical) %mdisable all these checks!

%pad and add linear tapering at top and bottom?
Z = nanmean(Data.Alt,1);
addParameter(p,'PadandTaper',   false, @islogical);                                    %activation flag
addParameter(p,'TaperLimits', [min([0,min(Z)]),max([100,max(Z)])], @isreal   );        %new top and bottom of the data, in km altitude. Will be ignored if this is within the data or if either value is NaN
addParameter(p,'TaperLength',       5, @(x) validateattributes(x,{'numeric'},{'>',0})) %length of linear taper at each end, in km
clear Z

%physical constants
addParameter(p,'N',0.02,@isnumeric)  %Brunt-Vaisala frequency
addParameter(p,'g',9.81,@isnumeric)  %acceleration due to gravity

%number of peaks
addParameter(p,'NPeaks',1,@(x) validateattributes(x,{'numeric'},{'integer','>=',1})) %number of peaks to try and find in each spectrum

%verbosity
addParameter(p,'Verbose',  true,@islogical) %print progress updates?

%interpolate to a regular height grid? (you probably want to keep this set to true - it checks if it's already true,
%and only does the interpolation if needed)
addParameter(p,'RegulariseZ',true,@islogical)

%additional variables - planetary wave grid fit filter ('PWGrid')
addParameter(p,      'NPWs',                  3,@isreal) %number of planetary waves to fit
addParameter(p,   'PWMinPC',                 66,@(x) validateattributes(x,{'numeric'},{'>',0})) %fraction of longitude bins that must be filled
addParameter(p,  'PWWindow',                  1,@(x) validateattributes(x,{'numeric'},{'>',0})) %window width of period used, in days
addParameter(p, 'PWTimeRes',                  1,@(x) validateattributes(x,{'numeric'},{'>',0})) %time resolution
addParameter(p, 'PWLonGrid',        -180:20:180,@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180,'vector'})) %longitude bins
addParameter(p, 'PWLatGrid',           -90:5:90,@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90,'vector'})) %latitude bins
addParameter(p, 'PWAltGrid',nanmean(Data.Alt,1),@(x) validateattributes(x,{'numeric'},{'>=',   0,         'vector'})) %altitude bins

%additional parameters - SGolay filter ('SGolay')
addParameter(p,'SGOrder',  2,@(x) validateattributes(x,{'numeric'},{'>',0})) %SGolay filter order
addParameter(p,'SGLength',18,@(x) validateattributes(x,{'numeric'},{'>',0})) %SGolay filter length (km)

%additional parameters - Hindley 2023 filter ('Hindley23')
addParameter(p,'H23_OutRem',  true,@islogical) %remove outliers?

%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,Data,varargin{:})

%pull out the remaining arguments into struct "Settings", used throughout rest of routine
Settings = p.Results;
Data = Settings.Data; Settings = rmfield(Settings,'Data');
Settings.TimeRange = [floor(min(Data.Time,[],'all','omitnan')),ceil(max(Data.Time,[],'all','omitnan'))];
clearvars -except InstInfo Settings Data



%check contents of input data struct:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%if we have loaded one of the Hindley23 data files but are not using the Hindley 23 filter (not sure why you'd do this 
% normally, but might be useful for cross-testing) then merge the temperature and PW fields to produce an integrated product
if ~strcmpi(Settings.Filter,'Hindley23') && ~isfield(Data,'Temp') && isfield(Data,'Temp_PW') && isfield(Data,'Temp_Residual');
Data.Temp = sum(Data.Temp_PW,3) + Data.Temp_Residual;
Data = rmfield(Data,{'Temp_PW','Temp_Residual'});
if isfield(Data,'Note'); Data = rmfield(Data,'Note'); end
end

%%do we have at least Lat, Lon, Alt? These are used by all filters
if ~isfield(Data,'Lat') ||  ~isfield(Data,'Lon') || ~isfield(Data,'Alt')
error('Missing field - struct must contain Lat, Lon, and Alt fields.')
return
end

%for all except Hindley23, we also need a 'Temp' field containing temperature
if ~strcmpi(Settings.Filter,'Hindley23') &  ~isfield(Data,'Temp');
error('Missing field - struct must contain Temp field.')
elseif strcmpi(Settings.Filter,'Hindley23');  
Data.Temp = NaN(size(Data.Lon)); %this will be ignored anyway in the analysis
end 

%%are these fields all the same size?
if ~isequal(size(Data.Lat),size(Data.Lon)) ||  ~isequal(size(Data.Lat),size(Data.Alt)) || ~isequal(size(Data.Lat),size(Data.Temp));
error('Lat, Lon, Alt and Temp fields must all be the same size.')
return
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ensure the data is on a regular height grid
% for Hindley23 removal, this must be done AFTER detrending, as new data are loaded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.RegulariseZ == true && ~strcmpi(Settings.Filter,'Hindley23'); Data = func_regularise_data_z(Data,Settings); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% detrend the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tg

%start by duplicating temperatures. All filter functions will act on this,
%in order to let us combine filters if we want
Data.Tp = Data.Temp;


%now, apply the filter. Currently only one at a time, but simple to rewrite in the future to apply multiple filters if needed
if     strcmpi(Settings.Filter,   'PWgrid'); [Data,PW] = func_filter_pwgrid(   Data,Settings);
elseif strcmpi(Settings.Filter,   'SGolay'); Data      = func_filter_sgolay(   Data,Settings);
elseif strcmpi(Settings.Filter,'Hindley23'); [Data,PW] = func_filter_hindley23(Data,Settings); 
else
error('Filter type not included in programme, stopping')
return
end

%create empty PW output if we didn't generate one. Note that this output will be 
%fairly meaningless anyway if the PW filter wasn't the first-applied
if ~exist('PW','var'); PW.Comment = 'Planetary wave filter not used, no PW data computed'; end


%regularisation if using Hindley23 detrending
if Settings.RegulariseZ == true && strcmpi(Settings.Filter,'Hindley23'); Data = func_regularise_data_z(Data,Settings); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% taper the profiles to cover a specific height range?
%note any zero-padding will be applied AFTER this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.PadandTaper

disp('******tapering is in testing, do not use******')
if Settings.Verbose == 1; disp('Applying padding and tapering'); end

%copy out vars to simplify code. will remove again after.
NewLength = Settings.TaperLimits;
TaperLength = Settings.TaperLength;

%if we're using scales not specified by the user,
%we will need to override them to use the new length
if mean(Settings.STScales ==  1:1:size(Data.Alt,2)/2) == 1; ReplaceScalesFlag = 1; else ReplaceScalesFlag = 0; end

%first, work out how much padding we need to get to the fill range
%we know the data are evenly spaced, but it's possible the top and bottom rows may be empty
%so let's fill any NaNs
Data.Alt = fillmissing(Data.Alt,'linear',2); %this is safe as the scales are, after above preprocessing, linearly spaced, monotonic, and the same in every profile

CurrentLimits     = minmax(Data.Alt(:));
dZ                = mean(diff(Data.Alt,1,2),'all','omitnan');
ExtraLevelsTop    = ceil((NewLength(2)-CurrentLimits(2))./dZ);
ExtraLevelsBottom = ceil((CurrentLimits(1)-NewLength(1))./dZ);

%now, zero-pad the data out to the requested full height range
Vars = fieldnames(Data);
if ExtraLevelsTop    > 0;
for iVar=1:1:numel(Vars);
% if strcmp(Vars{iVar},'Alt'); continue; end %handled below
sz = size(Data.(Vars{iVar})); sz(2) = ExtraLevelsTop;
Data.(Vars{iVar}) = cat(2,Data.(Vars{iVar}),zeros(sz));
end; clear iVar sz
end
if ExtraLevelsBottom > 0;
for iVar=1:1:numel(Vars);
% if strcmp(Vars{iVar},'Alt'); continue; end %handled below
sz = size(Data.(Vars{iVar})); sz(2) = ExtraLevelsBottom;
Data.(Vars{iVar}) = cat(2,zeros(sz),Data.(Vars{iVar}));
end; clear iVar sz
end
clear Vars

%fix altitudes separately as these values actually matter below
% (the rest will be wiped anyway at the end so can be ignored)
Data.Alt(Data.Alt == 0) = NaN;
Data.Alt = fillmissing(Data.Alt,'linear',2); %same logic as before re: safety.
NewZ     = nanmean(Data.Alt,1);

%find the indices of the original top and bottom.
%this is both to apply the taper and also to put the data back later as it was
PreTaperLimits = [find(NewZ == CurrentLimits(1)),find(NewZ == CurrentLimits(2))];

%finally, taper the added regions
nz = make_odd(TaperLength./dZ);
if nz == 0; nz = 1; end %forces there to be at least one level of taper (assuming we have any padding levels)
tapermultiplier = linspace(1,0,nz+2); tapermultiplier = tapermultiplier(2:end-1);
for iLev=1:1:nz;
if PreTaperLimits(1)-iLev > 0;                Data.Tp(:,PreTaperLimits(1)-iLev) = Data.Tp(:,PreTaperLimits(1)).*tapermultiplier(iLev); end %bottom
if PreTaperLimits(2)+iLev < size(Data.Alt,2); Data.Tp(:,PreTaperLimits(2)+iLev) = Data.Tp(:,PreTaperLimits(2)).*tapermultiplier(iLev); end %top
end; clear iLev
clear CurrentLimits dZ ExtraLevelsBottom ExtraLevelsTop NewZ nz tapermultiplier NewLength TaperLength

%now, replace scales if default set
if ReplaceScalesFlag == 1; 
if Settings.Verbose == 1; disp('Data padded and tapered - replacing scales'); end
Settings.STScales = 1:1:size(Data.Alt,2)/2; 
end


%REMEMBER - WE NEED TO REMOVE THE ADDED REGIONS AFTER THE ST!

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S-Transform profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%produce storage arrays
NProfiles = size(Data.Tp,1);
NLevs     = size(Data.Tp,2);
OutData   = spawn_uniform_struct({'A','Lz','Lh','Lat','Lon','Alt','Tp','MF','Time','FailReason','BG','Ep'},[NProfiles,NLevs,Settings.NPeaks]);
Mask      = ones([NProfiles,NLevs]); %this is used to mask out bad data later

%some approaches require two adjacent profiles to be computed. To avoid duplicate computation in this case,
%it is marginally more efficient if we work backwards and store the 'previous' (i.e. next) profile to permit this.

if Settings.Verbose == 1; textprogressbar('--> Computing gravity waves '); end
for iProf=NProfiles:-1:1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S-Transform the profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fill any NaNs with a zero for purposes of the ST
%we'll replace with NaNs at the end of the programme
NoData = find(isnan(Data.Tp(iProf,:)));
Tp = Data.Tp(iProf,:); Tp(NoData) = 0;
Mask(iProf,NoData) = 0;
clear NoData
if nansum(Tp) == 0; continue; end % no data
if numel(find(abs(Tp) > 0)) < 5; continue; end %not enough points

%zero-pad the data to prevent FFT wraparound
Tp = [zeros(1,Settings.STPadSize),Tp,zeros(1,Settings.STPadSize)];

%compute ST
ThisST = nph_ndst(Tp,                                ...
Settings.STScales,                 ...
nanmean(diff(Data.Alt(iProf,:))),  ...
Settings.STc,                      ...
'minwavelengths',Settings.MinLz,   ...
'maxwavelengths',Settings.MaxLz);
clear Tp 

%remove the zero-padding region from all output variables
Fields = {'IN','F1','A','R','HA','HR','allgws','BoostFactor','C'};
for iF=1:1:numel(Fields); 
F = ThisST.(Fields{iF});
F = F(Settings.STPadSize+1:end-Settings.STPadSize);
ThisST.(Fields{iF}) = F;
end; clear Fields iF F
ThisST.ST = ThisST.ST(:,Settings.STPadSize+1:end-Settings.STPadSize);

%store and retain full ST field?
if Settings.FullST == true


%if we don't have one, create a storage array
if ~isfield(OutData,'FullST');
OutData.FullST = NaN([size(Data.Tp),numel(Settings.STScales)]);
OutData.FullST = complex(OutData.FullST,0);
OutData.Freqs = ThisST.freqs;
end

%plug the data in
OutData.FullST(iProf,:,:) = ThisST.ST';

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute and store results depending on analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Analysis == 1 && Settings.NPeaks == 1;

%simple 1DST approach with one peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutData.A(   iProf,:) = ThisST.A;
OutData.Lz(  iProf,:) = 1./ThisST.F1;
OutData.Lat( iProf,:) = Data.Lat( iProf,:); OutData.Lon(iProf,:) = Data.Lon(iProf,:);
OutData.Alt( iProf,:) = Data.Alt( iProf,:); OutData.Tp( iProf,:) = Data.Tp( iProf,:);
OutData.Time(iProf,:) = Data.Time(iProf,:); OutData.FailReason(iProf,:) = 0;
OutData.BG(  iProf,:) = Data.BG(  iProf,:);
OutData.Ep(  iProf,:) = 0.5 .* (Settings.g ./ Settings.N).^2 .* (ThisST.A ./ Data.BG(iProf,:)).^2;

elseif Settings.Analysis == 1 && Settings.NPeaks > 1;     

%simple 1DST approach with multiple peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute absolute spectrum for this profile, and drop modal frequency
AbsSpec = abs(ThisST.ST);
AbsSpec(1,:) = NaN;

%loop over heights and find peak amplitudes and wavelengths
for iLev=1:1:NLevs
[pkval,pkidx] = findpeaks(AbsSpec(:,iLev));
[~,idx] = sort(pkval,'desc'); pkidx = pkidx(idx);%sort by magnitude
for iPeak=1:1:Settings.NPeaks
if iPeak <= numel(pkidx); 
OutData.A( iProf,iLev,iPeak) = AbsSpec(pkidx(iPeak),iLev);
OutData.Lz(iProf,iLev,iPeak) = 1./ThisST.freqs(pkidx(iPeak));
end
end
end
clear iLev pkval pkidx idx iPeak

%compute supporting values
for iPeak=1:1:Settings.NPeaks
OutData.Lat( iProf,:,iPeak) = Data.Lat( iProf,:); OutData.Lon(iProf,:,iPeak) = Data.Lon(iProf,:);
OutData.Alt( iProf,:,iPeak) = Data.Alt( iProf,:); OutData.Tp( iProf,:,iPeak) = Data.Tp( iProf,:);
OutData.Time(iProf,:,iPeak) = Data.Time(iProf,:); OutData.BG( iProf,:,iPeak) = Data.BG( iProf,:);
OutData.Ep(  iProf,:,iPeak) = 0.5 .* (Settings.g ./ Settings.N).^2 .* (OutData.A( iProf,:,iPeak) ./ Data.BG(iProf,:)).^2;
OutData.FailReason(iProf,:,iPeak) = isnan(squeeze(OutData.A(iProf,:,iPeak))).*5;
end; clear iPeak


elseif Settings.Analysis == 2

%Wright and Gille (GRL,2013) approach, but without the noise filter
%apply_wg13() can apply this filter to this function's output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if this is the first profile to be processed, retain the data then move to the next
if iProf == NProfiles || ~exist('NextST','var'); NextST = ThisST; continue; end

%if one of the profiles has a lot of missing data, the STs will be different sizes. Skip if so
if ~isequal(size(ThisST.ST),size(NextST.ST)); continue; end

%assume success until proven otherwise
OutData.FailReason(iProf,:) = 0;

%compute cospectrum
CoSpectrum = ThisST.ST .* conj(NextST.ST);

%find time and space separation of the two profiles, 
dx = nph_haversine([Data.Lat(iProf,  :);Data.Lon(iProf,  :)]', ...
[Data.Lat(iProf+1,:);Data.Lon(iProf+1,:)]')';
dt = abs(Data.Time(iProf,:) - Data.Time(iProf+1,:)).*60.*60.*24;
if   Settings.NoDistDiscard == 0; Bad = find(dt > Settings.Maxdt | dx > Settings.MaxdX);
else;                             Bad = [];
end

%discard the profile completely if we fall below the minimum acceptable fraction of safe data
if Settings.NoDistDiscard ~= 1 & (numel(Bad)/numel(dx) > 1-Settings.MinFracInProf) ; 
OutData.FailReason(iProf,:) = 2; 
clear Bad dt dx CoSpectrum; 
continue; 
end

%if we pass the above, discard any heights where we failed either individually
if numel(Bad) > 0; CoSpectrum(:,Bad) = NaN; OutData.FailReason(iProf,Bad) = 1; end

%drop modal frequency
CoSpectrum(1,:) = NaN;

%take sqrt(abs()) of cospectrum, i.e. the covarying amplitude 
CVA = sqrt(abs(CoSpectrum));

%create storage arrays
A = NaN(NLevs,Settings.NPeaks);
Lz = A;
Lh = A;

%loop over heights and find peak amplitudes and wavelengths
for iLev=1:1:NLevs

if Settings.NPeaks == 1;
%locate maximum at this height
[~,pkidx] = nanmax(CoSpectrum(:,iLev),[],1,'omitnan');
else

%locate all local maxima at this height, including the main maximum
[pkval,pkidx] = findpeaks(CVA(:,iLev));

%sort by magnitude
[~,idx] = sort(pkval,'desc'); pkidx = pkidx(idx);
clear idx pkval
end

for iPeak=1:1:Settings.NPeaks
if iPeak <= numel(pkidx);

%amplitude and Lz are straight copy-overs
A( iLev,iPeak) = CVA(pkidx(iPeak),iLev);
Lz(iLev,iPeak) = 1./ThisST.freqs(pkidx(iPeak));

%find phase change, and discard values where it's too small to be meaningful
dx(dx > Settings.MaxdX) = NaN;
dPhi = angle(CoSpectrum(pkidx(iPeak),iLev))./(2*pi);

%hence, compute Lh
Lh(iLev,iPeak) = abs(dx(iLev)./dPhi);

end
end; clear iPeak 
end; clear iLev pkval pkidx idx 

%we need pressure to calculate the density. If we don't have it, estimate it from altitude
%the density estimate is itself quite approximate, so this pressure estimate isn't a major source of error
if ~isfield(Data,'Pres'); Data.Pres = h2p(Data.Alt); end

%compute MF
BG = repmat(Data.BG(  iProf,:)',1,Settings.NPeaks);
P  = repmat(Data.Pres(iProf,:)',1,Settings.NPeaks);

MF = cjw_airdensity(P,BG)           ...
.* (Lz ./ Lh)                     ...
.* (Settings.g ./ Settings.N).^2  ...
.* A./BG.^2;

%compute potential energy
Ep = 0.5 .* (Settings.g ./ Settings.N).^2 .* (A ./ BG).^2;

%adjust lat/lon to the midpoint of the profile-pair
[latmean,lonmean] = meanm(Data.Lat(iProf+[0,1],:), ...
Data.Lon(iProf+[0,1],:));

%store results
OutData.A(   iProf,:,:) = A;
OutData.Lz(  iProf,:,:) = Lz;
OutData.Lh(  iProf,:,:) = Lh;
OutData.MF(  iProf,:,:) = MF;
OutData.Ep(  iProf,:,:) = Ep;
OutData.Lat( iProf,:,:) = repmat(latmean',           1,Settings.NPeaks);
OutData.Lon( iProf,:,:) = repmat(lonmean',           1,Settings.NPeaks);
OutData.Alt( iProf,:,:) = repmat(Data.Alt( iProf,:)',1,Settings.NPeaks);
OutData.Time(iProf,:,:) = repmat(Data.Time(iProf,:)',1,Settings.NPeaks);
OutData.Tp(  iProf,:,:) = repmat(Data.Tp(  iProf,:)',1,Settings.NPeaks);
OutData.BG(  iProf,:,:) = BG;

%store the new ST for the next pass
NextST = ThisST;

clear CoSpectrum A idx Lz dx dPhi Lh latmean lonmean MF dt Bad BG P CVA Ep 

end



if Settings.Verbose == 1; if mod(iProf,200); textprogressbar(100.*(NProfiles-iProf)./NProfiles); end; end
end; clear iProf NextST ThisST NextST NLevs NProfiles
if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%apply bad-data mask, to put NaNs back where we filled the data to ST it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f= fieldnames(OutData);
for iF=1:1:numel(f); 
F = OutData.(f{iF}); 

%%two special cases, both arising if we want to export a Full ST

if strcmpi(f{iF},'FullST') 
%complex ST output - 3D array of profs x heights x freqs
m = repmat(Mask,[1,1,numel(OutData.Freqs)]);
F(m == 0) = NaN;
clear m
elseif strcmpi(f{iF},'freqs');
%frequency scales for the complex ST output
continue %do nothing
else
%normal case - 2D array of profs x heights
F(Mask == 0) = NaN; 
end
OutData.(f{iF}) = F; end; 
OutData.FailReason(Mask == 0) = 4;
clear iF f F Mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove taper if it was applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.PadandTaper

OutData = reduce_struct(OutData,PreTaperLimits(1):1:PreTaperLimits(2),{'Freqs'},2);
end


%and we're done!
end




 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> func_filter_hindley23.m
%-----> Included version last modified:05-Nov-2023 23:24:50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hindley23 PW filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data,PW] = func_filter_hindley23(Data,Settings)


%copy the T' and T fields generated by Neil's code
Data.Tp   = Data.Temp_Residual;
Data.Temp = nansum(Data.Temp_PW,3);

%store PWs
PW = Data.Temp_PW;

%remove outliers
if Settings.H23_OutRem == true

%latitude and longitude have physical limits
Bad = [];
Bad = [Bad;find(Data.Lat <  -90 | Data.Lat >  90 | Data.Lon < -180 | Data.Lon > 180)];

%temperature should be >100K always, and  <400K at altitudes below the mesopause
Bad = [Bad;find(Data.Temp < 100)];
Bad = [Bad;find(Data.Temp > 400 & Data.Alt < 100)];

%do it
Bad = unique(Bad);
Fields = fieldnames(Data);
for iF=1:1:numel(Fields)
if ~isequal(size(Data.(Fields{iF})),size(Data.Alt)); continue; end
F = Data.(Fields{iF});
F(Bad) = NaN;
Data.(Fields{iF}) = F;
end
end

%compute background
Data.BG = Data.Temp - Data.Tp;

return
end

 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> func_filter_pwgrid.m
%-----> Included version last modified:05-Nov-2023 18:25:16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% planetary wave filter, grid-based approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data,PW] = func_filter_pwgrid(Data,Settings)


%compute a time grid to work on, and also store metadata about the PWs
PW.Time = min(Settings.TimeRange):Settings.PWTimeRes:max(Settings.TimeRange);
PW.Lon = Settings.PWLonGrid; PW.Lat = Settings.PWLatGrid; PW.Alt = Settings.PWAltGrid; PW.PWs = 0:1:Settings.NPWs;
PW.WindowSize = Settings.PWWindow; PW.MinPercent = Settings.PWMinPC;

%generate a store array for the PWs, and for the output T'
PW.PW = NaN(numel(Settings.PWLonGrid),numel(Settings.PWLatGrid),numel(Settings.PWAltGrid),Settings.NPWs+1,numel(PW.Time));
A = Data.Tp.*NaN; %working variable used internally to simplify logic

%fill it, stepping over day-by-day using a time window as specified
if Settings.Verbose == 1; textprogressbar('--> Computing planetary waves '); end

for iStep=1:1:numel(PW.Time)

%select the data we need by finding the indices of the points in the UseWindow (points to compute from) and
%the Output window (points to store, higher resolution)
%***logic assumes idxO is a subset of idxU***
idxU = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(Settings.PWWindow)./2);
idxO = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(mean(diff(PW.Time)))./2);
if numel(idxU) == 0 || numel(idxO) == 0; continue; end

%reduce data down to just what we want to use
PWCalcData = reduce_struct(Data,idxU,{'OriginalFiles','Note'},0);

%compute the PWs in the use window, and store it in placeholder A
%A will overwrite most loops - this is fine as long as idxO is a subset of idxU
[A(idxU),b] = pwfilter(Settings.NPWs,Settings.PWMinPC,                      ...
PWCalcData.Lon,PWCalcData.Lat,PWCalcData.Tp,   ...
Settings.PWLonGrid,Settings.PWLatGrid,               ...
PWCalcData.Alt,Settings.PWAltGrid);
%store the Tp and PW data
Data.Tp(idxO) = A(idxO);
PW.PW(:,:,:,:,iStep) = permute(b,[2,1,3,4]);

if Settings.Verbose == 1; textprogressbar(100.*iStep./numel(PW.Time)); end

end; clear iDay OutWindow UseWindow idxU PWCalcData a b idxO iStep
if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!');end

%compute background temperature
I = griddedInterpolant({PW.Lon,PW.Lat,PW.Alt,PW.Time},squeeze(nansum(PW.PW,4)));
Data.BG = I(Data.Lon,Data.Lat,Data.Alt,Data.Time);

return
end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> func_filter_sgolay.m
%-----> Included version last modified:05-Nov-2023 18:20:22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D Savitzky-Golay filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Data = func_filter_sgolay(Data,Settings)

%create array to store background temperature
%put this in an if, in case we're calling multiple filters
if ~isfield(Data,'BG'); Data.BG = Data.Tp.*NaN; GetBG = 1; else GetBG = 0; end


%is the data all the same resolution?
if nanstd(flatten(diff(Data.Alt,1,2)))./nanmean(flatten(diff(Data.Alt,1,2))) < 0.01;
%if so, we can do this in a single pass
FrameLen = abs(make_odd(round(Settings.SGLength./nanmean(flatten(diff(Data.Alt,1,2))))));
BG = sgolayfilt(Data.Tp',Settings.SGOrder,FrameLen)';
if GetBG == 1; Data.BG = BG; end
Data.Tp  = Data.Temp-BG;
else

%if not, we need a loop
for iProf=1:1:numel(Data.Alt,1)
FrameLen = abs(make_odd(round(Settings.SGLength./nanmean(Data.Alt(iProf,:)))));
BG = sgolayfilt(Data.Tp(iProf,:),Settings.SGOrder,FrameLen);
if GetBG == 1; Data.BG(iProf,:) = BG; end
Data.Tp(iProf,:) = Data.Tp(iProf,:) - BG;
end
end

if Settings.Verbose == 1; disp(['--> Savitzky-Golay filter applied vertically']); end

return
end



 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> func_regularise_data_z.m
%-----> Included version last modified:05-Nov-2023 18:18:46
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ensure data is on a regular Z grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Data = func_regularise_data_z(Data,Settings)


%first, check if the data is ALREADY regular AND ascending in z. Counts if the full distribution is within 5% of the mean.
%if it's fine, we don't need to proceed
dZ_distrib = unique(diff(Data.Alt,1,2));


if min(diff(Data.Alt,1,2),[],'all') > 0 & range(dZ_distrib) < 0.1.* nanmean(dZ_distrib); return; end


%ok, we need to regularise. First, work out a scale
NewZ = nanmin(Data.Alt(:)):nanmean(abs(dZ_distrib)):nanmax(Data.Alt(:));
clear dZ_distrib


%create struct to store new data
NewData = struct();
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
for iF=1:1:numel(Fields); NewData.(Fields{iF}) = NaN([size(Data.Alt,1),numel(NewZ)]); end

%now interpolate EVERYTHING onto the new scale
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
if Settings.Verbose == 1; disp('--> Data not on a regular and common height grid, interpolating to [ min Z: mean dZ : max Z ] '); end

for iF=1:1:numel(Fields)

%get data fields
Fin  =    Data.(Fields{iF});
Fout = NewData.(Fields{iF}); 

%check we're working on a 2D profiles x heights array, and just pass it straight through if not
%do this by chekcing if it's the same size as the 'Alt' array (arbitrary choice)
if ~isequal(size(Fin),size(Data.Alt)); NewData.(Fields{iF}) = Fin; continue; end


%interpolate profile-by-profile
for iProfile=1:1:size(Data.Alt,1);
Good = find(~isnan(Data.Alt(iProfile,:)));
if numel(Good) < 2; continue; end
Fout(iProfile,:) = interp1(Data.Alt(iProfile,Good),Fin(iProfile,Good),NewZ);
end

%store data
NewData.(Fields{iF}) = Fout;
end

%store new altitudes
NewData.Alt = repmat(NewZ,size(Data.Alt,1),1);

%copy over, and return
Data = NewData;





return
end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> expose_dim.m
%-----> Included version last modified:22-Dec-2024 11:40:02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function [Matrix,DimSize,DimOrder] = expose_dim(Array,Dim)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Take an n-dimensional Matlab array and reshape such that the array
%becomes 2D with the chosen dimension as the first and all other 
%dimensions merged in the second
%
%To restore the original dimensions assuming no changes to size, use this syntax:
%    Restored = permute(reshape(Matrix,DimSize),DimOrder);
%
%If changes have been made to the exposed dimension, e.g. interpolation or selection, use this syntax:
%    Restored = permute(reshape(Matrix,[size(Matrix,1),DimSize(2:end)]),DimOrder);
%
%inputs:
%  Array: n-dimensional input array
%  Dim:   dimension to bring to the front
%
%outputs:
%  Matrix:   2D output array, as described above
%  DimSize:  size of dimensions needed to put the array back using reshape()
%  DimOrder: order of dimensions needed to put the array back using permute()
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the reshaping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get size of original array
sz = size(Array);

%permute desired dimension to front
DimOrder = unique([Dim,1:1:numel(sz)],'stable');

%reshape to make all other dimensions lines
Matrix = reshape(permute(Array,DimOrder),[sz(Dim),prod(sz(DimOrder(2:end)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%work out the information we need to put it back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DimSize is used with reshape() to put all the dimensions back to the right size
DimSize  = [size(Matrix,1),sz(DimOrder(2:end))];

%NewOrder is used with permute() to put all the dimensions back in the right order
DimOrder = 1:1:numel(sz);
DimOrder = [DimOrder(DimOrder < Dim)+1,1,DimOrder(DimOrder > Dim)];

end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> index_dim.m
%-----> Included version last modified:21-Dec-2024 11:25:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function Array = index_dim(Array,Indices,Dim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select elements of an array along a selected dimension
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/12 
%
%inputs:
%  Array   - the array we want to operate on
%  Dim     - the dimension we want to operate on
%  Indices - the indices we want to extract along dimension Dim
%
%outputs:
%  Array   - the selected array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input handling
if nargin      < 3; Dim = 1; end     %assume first dimension if not specified, and set varargin to blank if not set
if numel(Dim) == 0; Dim = 1; end     %if set to blank, then set first dimension

%expose the desired dimension
[y,a,b] = expose_dim(Array,Dim);

%select
y = y(Indices,:);

%put the dimensions back in order 
Array = permute(reshape(y,[size(y,1),a(2:end)]),b);

return
end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> list_non_modal_size.m
%-----> Included version last modified:21-Dec-2024 11:25:03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function VarList = list_non_modal_size(Struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Find the modal field size in a struct, and return a list of all fields not this size
%inputs:
%  Struct   - the struct to operate on
%
%outputs:
%  VarList  - cell array lsiting all variables not the modal size

%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/01
%updated 2023/08/12 to let the user choose a dimension to operate along
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VarList = {};

%get largest number of dimensions for any field
f = fieldnames(Struct); MaxDims = 0;
for iF=1:1:numel(f); if ndims(Struct.(f{iF})) > MaxDims; MaxDims = ndims(Struct.(f{iF})); end; end
clear iF

%get size of each field
Sizes = NaN([numel(f),MaxDims]);
for iF=1:1:numel(f)
sz = size(Struct.(f{iF}));
Sizes(iF,1:numel(sz)) = sz;
end
clear iF sz MaxDims

%find most common size
[U,I,J] = unique(Sizes, 'rows');
Score = NaN(numel(U),1);
for iU=1:1:numel(U); Score(iU) = numel(find(J == iU)); end
[~,idx] = max(Score);
ModalSize = U(idx,:);
clear U I J Score iU idx

%finally, get a list of vars to ignore
for iF=1:1:numel(f); if ~isequal(size(Struct.(f{iF})),ModalSize); VarList{end+1} = f{iF}; end; end

clear iF f IgnoreWrongSize ModalSize Sizes

return
end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> reduce_struct.m
%-----> Included version last modified:22-Dec-2024 11:42:23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function Struct = reduce_struct(Struct,SubSetIndices,VarsToExclude,Dim,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%For a series of identical fields in a struct(), select a given set of 
%indices from every field
%
%For backwards-compatability with older versions, VarsToExclude argument needs to be
%ordered before Dim. If we want to operate on ALL fields, just set this to [].
%
%inputs:
% required:
%  Struct        - the struct to operate on
%  SubSetIndices - the list of indices to select from each field
%  VarsToExclude - fields to ignore when subsettings. Can be set to just {}.
%  Dim           - the dimension to operate on for each field. If 0, will apply to whole dataset
% optional:
%  'IgnoreWrongSize' if set to true will skip any variables that are not the MOST COMMON size. Default false.
%
%outputs:
%  Struct        - the  reduced structure

%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/01
%updated 2023/08/12 to let the user choose a dimension to operate along
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%handle inputs
if ~exist('VarsToExclude','var'); VarsToExclude = {' '}; end  %assume applies to all vars
if ~exist(          'Dim','var'); Dim = 1;               end  %assume operation is along first dimension
if ~exist('varargin','var'); varargin = {}; end
for iV=1:2:numel(varargin);
if strcmpi(varargin{iV},'IgnoreWrongSize'); IgnoreWrongSize = varargin{iV+1}; end
end; clear iV
if ~exist('IgnoreWrongSize','var'); IgnoreWrongSize = false; end

%if IgnoreWrongSize is true, add wrong-size variables to the ignore list
if IgnoreWrongSize == true; VarsToExclude = [VarsToExclude,list_non_modal_size(Struct)];end


Fields = fieldnames(Struct);
for iField=1:1:numel(Fields);

%skip specified variables
if any(strcmp(Fields{iField},VarsToExclude)); continue; end

%reduce desired variables
F = Struct.(Fields{iField});
if strcmpi(class(F),'table')
%tables can only be trimmed in the first two dimensions
if Dim == 1;
F = F(SubSetIndices,:);
elseif Dim == 2;
F = F(:,SubSetIndices);
else
warning('reduce_struct is trying to subset a table in a dimension higher than 2, skipping')
end
else
%anything else, treat as normal
if Dim == 0; F = index_dim(F(:),SubSetIndices,1);
else;        F = index_dim(F,SubSetIndices,Dim);
end
end
Struct.(Fields{iField}) = F;

%done!

end

return
end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> spawn_uniform_struct.m
%-----> Included version last modified:06-Sep-2023 07:42:37
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function OutStruct = spawn_uniform_struct(Fields,Size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spawn a struct containing many identically-shaped NaN fields
%Corwin Wright, c.wright@bath.ac.uk, 2023/09/06 
%
%inputs:
%  Fields    -cell array containing list of fields, e.g. {'A','B'}
%  Size      - size of each array
%
%outputs:
%  OutStruct - the output struct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutStruct = struct();
for iField=1:1:numel(Fields)

OutStruct.(Fields{iField}) = NaN(Size);

end


end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> cjw_airdensity.m
%-----> Included version last modified:07-Jan-2025 16:34:09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function AirDensity = cjw_airdensity(Pressure,Temperature)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute air density from pressure and temperature
%matlab reimplementation of:
% http://www.iac.ethz.ch/staff/dominik/idltools/atmos_phys/air_density.pro
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%10/JAN/2014
%
%inputs
%---------
%
%Temperature - temperature  (K)
%Pressure - pressure (hPa)
%
%outputs
%---------
%
%AirDensity - air density, kg/m^3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Rd = 286.9; %Gas constant for dry air

AirDensity = Pressure .* 1e2  ./ (Rd .* Temperature);

end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> h2p.m
%-----> Included version last modified:02-Apr-2025 11:03:06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 

function POut = h2p(alt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert pressures to altitudes
%
%inputs:
%   alt - altitude levels, in km (any size of array is fine)
%outputs:
%   p - atmospheric pressure values, in hPa
%
%Corwin Wright, c.wright@bath.ac.uk, 2025/04/02 (full rewrite of older 2023 version to properly vectorise,
%which was itself based on much older code by Neil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%define pressure and altitude bounds for each level
p = [1013.25, 226.3210, 54.7489,  8.6802,  1.1091,  0.6694,  0.0396];
z = [   0,     11,      20,      32,      47,      51,      71];  

% calculate scale heights within each level
H = -diff(z)./log(p(2:end)./p(1:end-1));  % Scale heights (km)

%create output array
POut = NaN(size(alt)); 

%apply to each altitude **band**
for i = 1:length(z)-1
idx = (alt >= z(i) & alt < z(i+1));                 % find altitudes in range
POut(idx) = p(i) * exp(-(alt(idx) - z(i)) / H(i));  %compute output pressure
end

%deal with high altitudes (above 71km)
POut(alt >= z(end)) = p(end) * exp(-(alt(alt >= z(end)) - z(end)) / H(end));


end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> nph_haversine.m
%-----> Included version last modified:27-Apr-2021 16:26:09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
% 
% loc1 = [63.2 -30 ; 63.2 0 ; 63.2 45];
% loc2 = [13.95 150; 13.95 180 ; 13.95 225];
% 
% 


function [km] = nph_haversine(loc1,loc2)


% If your input is an array of lats/lons, it should be of the form:

% loc1 = [N_POSITIONS x LATLON],
% e.g. for 100 positions loc1 = <100x2 double>.


% % HAVERSINE     Compute distance between locations using Haversine formula
% %   KM = HAVERSINE(LOC1, LOC2) returns the distance KM in km between
% %   locations LOC1 and LOC2 using the Haversine formula.  LOC1 and LOC2 are
% %   latitude and longitude coordinates.
% %   
% %   Inputs
% %       loc1 = [LAT1 LON1]
% %       loc2 = [LAT2 LON2]
% % 
% %   Notes
% %       The Haversine formula is used to calculate the great-circle
% %       distance between two points, which is the shortest distance over
% %       the earth's surface.
% % 
% % For a cosmic occultation at h_cosmic ~ 800km with a GPS satellite, a
% % haversine distance of ~11436.398km is found. You can check with geometry
% % the following:
% % loc1 = [63.2 -30];
% % loc2 = [13.95 150];
% % 
% % EDIT - made array safe. I think.
% % 

%%

loc1 = deg2rad(loc1); loc2 = deg2rad(loc2);

R = 6371;                                 % Earth's radius in km
delta_lat = loc2(:,1) - loc1(:,1);        % difference in latitude
delta_lon = loc2(:,2) - loc1(:,2);        % difference in longitude
a = sin(delta_lat./2).^2 + cos(loc1(:,1)) .* cos(loc2(:,1)) .* sin(delta_lon./2).^2;
c = 2 .* atan2(sqrt(a), sqrt(1-a));
km = R .* c;

% 
% %%
% 
% % Get azimuth too? BROKEN I THINK :(
% 
% lat1 = loc1(:,1); lon1 = loc1(:,2);
% lat2 = loc2(:,1); lon2 = loc2(:,2);
% 
% 
% az = atan2(cos(lat2) .* sin(lon2-lon1),...
%            cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(lon2-lon1));
% 
% % Azimuths are undefined at the poles, so we choose a convention: zero at
% % the north pole and pi at the south pole.
% az(lat1 <= -pi/2) = 0;
% az(lat2 >=  pi/2) = 0;
% az(lat2 <= -pi/2) = pi;
% az(lat1 >=  pi/2) = pi;
% 

end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> pwfilter.m
%-----> Included version last modified:07-Jan-2025 16:34:21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function [VarOut,PWStore] = pwfilter(NPWs,MinPC,Lon,Lat,Var,LonGrid,LatGrid,Alt,AltGrid)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simple planetary wave filter for application to scattered data as f(lon,lat,alt)
%altitude is optional - omit last two arguments to apply to 2D data only
%
%the first mode will always be the zonal mean, then the other modes
%will start from number 2 (i.e. modes 1,2,3 are in positions 2,3,4, etc)
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/09/23
%
%inputs (not in order - check above):
% NPWs    - number of planetary wave modes to fit. Zonal mean will always be included.
% MinPC   - minimum percentage of filled bins in a longitude band to compute PWs - will set to NaN otherwise
% Lon     - list of longitudes  at data points
% Lat     - list of latitude    at data points
% Alt     - list of altitudes   at data points (optional)
% Var     - list of data values at data points
% LonGrid - longitude grid to compute results on (PW outputs will also be on this grid)
% LatGrid - latitude  grid compute results on (PW outputs will also be on this grid)
% AltGrid - altitude  grid compute results on (PW outputs will also be on this grid)
%
%Lon, Lat, Alt and Var must all be the same size in the same order, and will be flattened
%LonGrid,LatGrid and AltGrid should be vectors, and will be meshgridded
%
%outputs:
% VarOut - the total PW contribution (inc zonal mean) at each point in the raw data
% PWStore - planetary wave estimated amplitude at each point on the output grid, for each mode (lon x lat x alt x PW)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check input validity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if we didn't specify altitude values and ranges, create dummy values
if nargin ~= 9; Alt = ones(size(Lon)); AltGrid = 1;end

%check input types, ranges and consistency.
p = inputParser;
addRequired(p,'NPWs',    @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'}));
addRequired(p,'MinPC',   @(x) validateattributes(x,{'numeric'},{'nonnegative','<=',100,'scalar'}));
addRequired(p,'Lon',     @isnumeric);
addRequired(p,'Lat',     @(x) validateattributes(x,{'numeric'},{'size',size(Lon)}));
addRequired(p,'Var',     @(x) validateattributes(x,{'numeric'},{'size',size(Lon)}));
addRequired(p,'LonGrid', @(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',360,'vector'}));
addRequired(p,'LatGrid', @(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90,'vector'}));
addRequired(p,'Alt',     @(x) validateattributes(x,{'numeric'},{'nonnegative','size',size(Lon)}));
addRequired(p,'AltGrid', @(x) validateattributes(x,{'numeric'},{'nonnegative','vector'}));
parse(p,NPWs,MinPC,Lon,Lat,Var,LonGrid,LatGrid,Alt,AltGrid);
clear p



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flatten input fields, retaining the input size
InSize = size(Lon);
Lon = Lon(:); Lat = Lat(:); Alt = Alt(:); Var = Var(:);

%fill gaps in altgrid (these arise from interaction of get_limbsounders() and heights
%where we have no data, the 'levels' produced will have no data to fill anyway so this
%has no effect on the final results
AltGrid = inpaint_nans(AltGrid);

%make sure everything ascends monotonically
LatGrid= sort(LatGrid,'asc'); LonGrid = sort(LonGrid,'asc'); AltGrid = sort(AltGrid,'asc');

%create output meshgrids
[LonGrid,LatGrid,AltGrid] = meshgrid(LonGrid,LatGrid,AltGrid);

%bin the input data onto the grid
VarGrid = bin2matN(3,Lon,Lat,Alt,Var,LonGrid,LatGrid,AltGrid,'@nanmean');

%reshape to remove nested loop
sz = size(VarGrid); if numel(sz) < 3; sz = [sz,1,1]; end
VarGrid = reshape(permute(VarGrid,[2,1,3]),sz(2),sz(1)*sz(3));

%drop rows that don't meet MinFrac
Filled = sum(~isnan(VarGrid),1);
VarGrid(:,Filled./sz(2) <= MinPC./100) = NaN;

clear sz Filled Bad MinPC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PW analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract longitude axis
LonAxis = LonGrid(1,:,1);

%create storage arrays
PWStore = NaN([size(VarGrid),NPWs+1]); %map of PWs
VGrid   = zeros(size(VarGrid));        %sum of planetary waves at each point

%compute planetary wave value at each data point
for iPW=0:1:NPWs
for iLine =1:1:size(VarGrid,2)

if nansum(VarGrid(:,iLine)) == 0; continue; end %not enough data to pass percentage check above

if iPW==0;
%take and repmat zonal mean
yfit = repmat(nanmean(VarGrid(:,iLine)),[size(VarGrid,1),1]);
else
%fit PW
warning off
yfit = nph_sinefit(LonAxis,VarGrid(:,iLine),360./iPW);
warning on
end

%remove from the data so we don't do so again
VarGrid(:,iLine) = VarGrid(:,iLine) - yfit;

%store results
VGrid(:,iLine) = VGrid(:,iLine)+yfit;
PWStore(:,iLine,iPW+1) = yfit;

clear yfit F
end; clear iLine
end; clear iPW
clear VarGrid LonAxis

%reshape back to original grid shape
sz = size(LatGrid); if numel(sz) < 3; sz = [sz,1,1]; end
VGrid   = permute(reshape(  VGrid,sz([2,1,3])),[2,1,3]);
PWStore = permute(reshape(PWStore,[sz([2,1,3]),NPWs+1]),[2,1,3,4]);
clear sz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove from raw profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create an interpolant object of the PW sum
%this requires permuting each input from meshgrid to ndgrid
%it also requires splitting out by the 2D and 3D case. Bah, so close.

nd = ndims(LonGrid);
if nd == 3;
%3D case

I = griddedInterpolant(permute(LonGrid,[2,1,3]),...
permute(LatGrid,[2,1,3]),...
permute(AltGrid,[2,1,3]),...
permute(VGrid,  [2,1,3]));


%interpolate to the profile locations, difference from the raw data, 
%and reshape to the original field shape
VarOut = reshape(Var-I(Lon,Lat,Alt),InSize);

elseif nd == 2;

I = griddedInterpolant(permute(LonGrid,[2,1]),...
permute(LatGrid,[2,1]),...
permute(VGrid,  [2,1]));


%interpolate to the profile locations, difference from the raw data, 
%and reshape to the original field shape
VarOut = reshape(Var-I(Lon,Lat),InSize);

else
disp('Error - invalid number of output dimensions')
end


clear I InSize Lon Lat Alt Var AltGrid LonGrid LatGrid VGrid nd

end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> flatten.m
%-----> Included version last modified:26-Nov-2020 13:08:52
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function Out = flatten(In)


Out = In(:);


end

 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> inrange.m
%-----> Included version last modified:14-Oct-2019 11:07:32
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function InRange = inrange(Array,MinMax,NoEnds)

if nargin < 3; NoEnds = 0; end

if NoEnds ~= 1;
InRange = find(Array >= min(MinMax) & Array <= max(MinMax));
else
InRange = find(Array >  min(MinMax) & Array <  max(MinMax));
end
return

end

 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> minmax.m
%-----> Included version last modified:10-Jan-2021 17:19:38
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function MinMax = minmax(Array)

TheMin = nanmin(Array(:));
TheMax = nanmax(Array(:));

MinMax = [TheMin,TheMax];
return

end

 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> bin2matN.m
%-----> Included version last modified:14-May-2022 19:30:20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function ZG = bin2matN(NDims,varargin)

%calls the appropriate version of bin2mat for 1d, 2d or 3d data
%in future, extend to ND, but this will do for now

switch NDims
case 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%call bin2mat1d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x    = double(varargin{1});
z    = double(varargin{2});
xi   = double(varargin{3});

if numel(varargin) < 4;
ZG = bin2mat1d(x,z,xi);
else
ZG = bin2mat1d(x,z,xi,varargin{4:end});
end

case 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x    = double(varargin{1});
y    = double(varargin{2});
z    = double(varargin{3});
xi   = double(varargin{4});
yi   = double(varargin{5});

if numel(varargin) < 6;
ZG = bin2mat2d(x,y,z,xi,yi);
else
ZG = bin2mat2d(x,y,z,xi,yi,varargin{6:end});
end

case 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat3d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x    = double(varargin{1});
y    = double(varargin{2});
z    = double(varargin{3});
v    = double(varargin{4});
xi   = double(varargin{5});
yi   = double(varargin{6});
zi   = double(varargin{7});

if numel(varargin) < 8;
ZG = bin2mat3d(x,y,z,v,xi,yi,zi);
else
ZG = bin2mat3d(x,y,z,v,xi,yi,zi,varargin{8:end});
end

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat1d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ZG = bin2mat1d(x,z,xi,varargin)
% BIN2MAT - create a matrix from scattered data without interpolation


%check inputs
error(nargchk(3,inf,nargin,'struct'));

%make sure the vectors are column vectors
x = x(:);
z = z(:);

if all(any(diff(cellfun(@length,{x,z}))));
error('Inputs x and z must be the same size');
end

%process optional input
fun=@mean;
test=1;
if ~isempty(varargin)
fun=varargin{1};
if ~isa(fun,'function_handle');
fun=str2func(fun);
end

%test the function for non-scalar output
test = feval(fun,rand(5,1),varargin{2:end});

end

%grid nodes
m=numel(xi);

%limit values to those within the specified grid
xmin=min(xi);
xmax=max(xi);

gind =(x>=xmin & x<=xmax);

%find the indices for each x and y in the grid
[junk,xind] = histc(x(gind),xi);


%break the data into a cell for each grid node
blc_ind=accumarray([xind],z(gind),[m 1],@(x){x},{NaN});

%evaluate the data in each grid using FUN
if numel(test)>1
ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind,'uni',0);
else
ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function ZG = bin2mat2d(x,y,z,XI,YI,varargin)
% BIN2MAT - create a matrix from scattered data without interpolation
%
%   ZG = BIN2MAT(X,Y,Z,XI,YI) - creates a grid from the data
%   in the (usually) nonuniformily-spaced vectors (x,y,z)
%   using grid-cell averaging (no interpolation). The grid
%   dimensions are specified by the uniformily spaced vectors
%   XI and YI (as produced by meshgrid).
%
%   ZG = BIN2MAT(...,@FUN) - evaluates the function FUN for each
%   cell in the specified grid (rather than using the default
%   function, mean). If the function FUN returns non-scalar output,
%   the output ZG will be a cell array.
%
%   ZG = BIN2MAT(...,@FUN,ARG1,ARG2,...) provides aditional
%   arguments which are passed to the function FUN.
%
%   EXAMPLE
%
%   %generate some scattered data
%    [x,y,z]=peaks(150);
%    ind=(rand(size(x))>0.9);
%    xs=x(ind); ys=y(ind); zs=z(ind);
%
%   %create a grid, use lower resolution if
%   %no gaps are desired
%    xi=min(xs):0.25:max(xs);
%    yi=min(ys):0.25:max(ys);
%    [XI,YI]=meshgrid(xi,yi);
%
%   %calculate the mean and standard deviation
%   %for each grid-cell using bin2mat
%    Zm=bin2mat(xs,ys,zs,XI,YI); %mean
%    Zs=bin2mat(xs,ys,zs,XI,YI,@std); %std
%
%   %plot the results
%    figure
%    subplot(1,3,1);
%    scatter(xs,ys,10,zs,'filled')
%    axis image
%    title('Scatter Data')
%
%    subplot(1,3,2);
%    pcolor(XI,YI,Zm)
%    shading flat
%    axis image
%    title('Grid-cell Average')
%
%    subplot(1,3,3);
%    pcolor(XI,YI,Zs)
%    shading flat
%    axis image
%    title('Grid-cell Std. Dev.')
%
% SEE also RESHAPE ACCUMARRAY FEVAL

% A. Stevens 3/10/2009
% astevens@usgs.gov

%check inputs
error(nargchk(5,inf,nargin,'struct'));

%make sure the vectors are column vectors
x = x(:);
y = y(:);
z = z(:);

if all(any(diff(cellfun(@length,{x,y,z}))));
error('Inputs x, y, and z must be the same size');
end

%process optional input
fun=@mean;
test=1;
if ~isempty(varargin)
fun=varargin{1};
if ~isa(fun,'function_handle');
fun=str2func(fun);
end

%test the function for non-scalar output
test = feval(fun,rand(5,1),varargin{2:end});

end

%grid nodes
xi=XI(1,:);
yi=YI(:,1);
[m,n]=size(XI);

%limit values to those within the specified grid
xmin=min(xi);
xmax=max(xi);
ymin=min(yi);
ymax=max(yi);

gind =(x>=xmin & x<=xmax & ...
y>=ymin & y<=ymax);

%find the indices for each x and y in the grid
[junk,xind] = histc(x(gind),xi);
[junk,yind] = histc(y(gind),yi);

%break the data into a cell for each grid node
blc_ind=accumarray([yind xind],z(gind),[m n],@(x){x},{NaN});

%evaluate the data in each grid using FUN
if numel(test)>1
ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind,'uni',0);
else
ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind);
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat3d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ZG = bin2mat3d(x,y,c,z,XI,YI,CI,varargin)

%BIN2MAT3D
% Bins data (Z) in 3D coordinates (x,y,c)  nto 3D bins specified by XI, YI,
% CI

% A. Stevens 3/10/2009
% astevens@usgs.gov

%check inputs
error(nargchk(7,inf,nargin,'struct'));

%make sure the vectors are column vectors
x = x(:);
y = y(:);
c = c(:);
z = z(:);

if all(any(diff(cellfun(@length,{x,y,z,c}))));
error('Inputs x, y,c and z must be the same size');
end

%process optional input
fun=@mean;
test=1;
if ~isempty(varargin)
fun=varargin{1};
if ~isa(fun,'function_handle');
fun=str2func(fun);
end

%test the function for non-scalar output
% test = feval(fun,rand(5,1),varargin{2:end});

end

%grid nodes
xi=squeeze(XI(1,:,1));
yi=squeeze(YI(:,1,1));
ci=squeeze(CI(1,1,:));
[m,n,l]=size(XI);

%limit values to those within the specified grid
xmin=min(xi);
xmax=max(xi);
ymin=min(yi);
ymax=max(yi);
cmin=min(ci);
cmax=max(ci);
gind =(x>=xmin & x<=xmax & ...
y>=ymin & y<=ymax & c>=cmin & c<=cmax );

%find the indices for each x and y in the grid

[junk,xind] = histc(x(gind),xi);
[junk,yind] = histc(y(gind),yi);
[junk,cind] = histc(c(gind),ci);

%break the data into a cell for each grid node
blc_ind=accumarray([yind xind cind],z(gind),[m n l],@(x){x},{NaN});

%evaluate the data in each grid using FUN
if numel(test)>1
ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind,'uni',0);
else
ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind);
end
end


 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> inpaint_nans.m
%-----> Included version last modified:23-Jan-2025 17:58:13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function B=inpaint_nans(A,method)
% INPAINT_NANS: in-paints over nans in an array
% usage: B=INPAINT_NANS(A)          % default method
% usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use
%       for the interpolation.) All methods are capable
%       of extrapolation, some are better than others.
%       There are also speed differences, as well as
%       accuracy differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor.
%       method  3 uses a better plate equation,
%                 but may be much slower and uses
%                 more memory.
%       method  4 uses a spring metaphor.
%       method  5 is an 8 neighbor average, with no
%                 rationale behind it compared to the
%                 other methods. I do not recommend
%                 its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a
%         linear system in the case of only a few
%         NaNs in a large array.
%         Extrapolation behavior is linear.
%         
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts
%         of the array which do not have any contact with
%         NaNs. Uses a least squares approach, but it
%         does not modify known values.
%         In the case of small arrays, this method is
%         quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%         
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%
%         Note: method 2 has problems in 1-d, so this
%         method is disabled for vector inputs.
%         
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result
%         in more accurate interpolations, at some cost
%         in speed.
%         
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero)
%         connect each node with every neighbor
%         (horizontally, vertically and diagonally)
%         Since each node tries to be like its neighbors,
%         extrapolation is as a constant function where
%         this is consistent with the neighboring nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element.
%         This method is NOT recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1);
%  z0 = exp(x+y);
%  znan = z0;
%  znan(20:50,40:70) = NaN;
%  znan(30:90,5:10) = NaN;
%  znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06


% I always need to know which elements are NaN,
% and what size the array is for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
method = 0;
elseif ~ismember(method,0:5)
error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
case 0
% The same as method == 1, except only work on those
% elements which are NaN, or at least touch a NaN.

% is it 1-d or 2-d?
if (m == 1) || (n == 1)
% really a 1-d case
work_list = nan_list(:,1);
work_list = unique([work_list;work_list - 1;work_list + 1]);
work_list(work_list <= 1) = [];
work_list(work_list >= nm) = [];
nw = numel(work_list);

u = (1:nw)';
fda = sparse(repmat(u,1,3),bsxfun(@plus,work_list,-1:1), ...
repmat([1 -2 1],nw,1),nw,nm);
else
% a 2-d case

% horizontal and vertical neighbors only
talks_to = [-1 0;0 -1;1 0;0 1];
neighbors_list=identify_neighbors(n,m,nan_list,talks_to);

% list of all nodes we have identified
all_list=[nan_list;neighbors_list];

% generate sparse array with second partials on row
% variable for each element in either list, but only
% for those nodes which have a row index > 1 or < n
L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
nl=length(L);
if nl>0
fda=sparse(repmat(all_list(L,1),1,3), ...
repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
repmat([1 -2 1],nl,1),nm,nm);
else
fda=spalloc(n*m,n*m,size(all_list,1)*5);
end

% 2nd partials on column index
L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
nl=length(L);
if nl>0
fda=fda+sparse(repmat(all_list(L,1),1,3), ...
repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
repmat([1 -2 1],nl,1),nm,nm);
end
end

% eliminate knowns
rhs=-fda(:,known_list)*A(known_list);
k=find(any(fda(:,nan_list(:,1)),2));

% and solve...
B=A;
B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);

case 1
% least squares approach with del^2. Build system
% for every array element as an unknown, and then
% eliminate those which are knowns.

% Build sparse matrix approximating del^2 for
% every element in A.

% is it 1-d or 2-d?
if (m == 1) || (n == 1)
% a 1-d case
u = (1:(nm-2))';
fda = sparse(repmat(u,1,3),bsxfun(@plus,u,0:2), ...
repmat([1 -2 1],nm-2,1),nm-2,nm);
else
% a 2-d case

% Compute finite difference for second partials
% on row variable first
[i,j]=ndgrid(2:(n-1),1:m);
ind=i(:)+(j(:)-1)*n;
np=(n-2)*m;
fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
repmat([1 -2 1],np,1),n*m,n*m);

% now second partials on column variable
[i,j]=ndgrid(1:n,2:(m-1));
ind=i(:)+(j(:)-1)*n;
np=n*(m-2);
fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
repmat([1 -2 1],np,1),nm,nm);
end

% eliminate knowns
rhs=-fda(:,known_list)*A(known_list);
k=find(any(fda(:,nan_list),2));

% and solve...
B=A;
B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);

case 2
% Direct solve for del^2 BVP across holes

% generate sparse array with second partials on row
% variable for each nan element, only for those nodes
% which have a row index > 1 or < n

% is it 1-d or 2-d?
if (m == 1) || (n == 1)
% really just a 1-d case
error('Method 2 has problems for vector input. Please use another method.')

else
% a 2-d case
L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
nl=length(L);
if nl>0
fda=sparse(repmat(nan_list(L,1),1,3), ...
repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
repmat([1 -2 1],nl,1),n*m,n*m);
else
fda=spalloc(n*m,n*m,size(nan_list,1)*5);
end

% 2nd partials on column index
L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
repmat([1 -2 1],nl,1),n*m,n*m);
end

% fix boundary conditions at extreme corners
% of the array in case there were nans there
if ismember(1,nan_list(:,1))
fda(1,[1 2 n+1])=[-2 1 1];
end
if ismember(n,nan_list(:,1))
fda(n,[n, n-1,n+n])=[-2 1 1];
end
if ismember(nm-n+1,nan_list(:,1))
fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
end
if ismember(nm,nan_list(:,1))
fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
end

% eliminate knowns
rhs=-fda(:,known_list)*A(known_list);

% and solve...
B=A;
k=nan_list(:,1);
B(k)=fda(k,k)\rhs(k);

end

case 3
% The same as method == 0, except uses del^4 as the
% interpolating operator.

% del^4 template of neighbors
talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
0 1;0 2;1 -1;1 0;1 1;2 0];
neighbors_list=identify_neighbors(n,m,nan_list,talks_to);

% list of all nodes we have identified
all_list=[nan_list;neighbors_list];

% generate sparse array with del^4, but only
% for those nodes which have a row & column index
% >= 3 or <= n-2
L = find( (all_list(:,2) >= 3) & ...
(all_list(:,2) <= (n-2)) & ...
(all_list(:,3) >= 3) & ...
(all_list(:,3) <= (m-2)));
nl=length(L);
if nl>0
% do the entire template at once
fda=sparse(repmat(all_list(L,1),1,13), ...
repmat(all_list(L,1),1,13) + ...
repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
else
fda=spalloc(n*m,n*m,size(all_list,1)*5);
end

% on the boundaries, reduce the order around the edges
L = find((((all_list(:,2) == 2) | ...
(all_list(:,2) == (n-1))) & ...
(all_list(:,3) >= 2) & ...
(all_list(:,3) <= (m-1))) | ...
(((all_list(:,3) == 2) | ...
(all_list(:,3) == (m-1))) & ...
(all_list(:,2) >= 2) & ...
(all_list(:,2) <= (n-1))));
nl=length(L);
if nl>0
fda=fda+sparse(repmat(all_list(L,1),1,5), ...
repmat(all_list(L,1),1,5) + ...
repmat([-n,-1,0,+1,n],nl,1), ...
repmat([1 1 -4 1 1],nl,1),nm,nm);
end

L = find( ((all_list(:,2) == 1) | ...
(all_list(:,2) == n)) & ...
(all_list(:,3) >= 2) & ...
(all_list(:,3) <= (m-1)));
nl=length(L);
if nl>0
fda=fda+sparse(repmat(all_list(L,1),1,3), ...
repmat(all_list(L,1),1,3) + ...
repmat([-n,0,n],nl,1), ...
repmat([1 -2 1],nl,1),nm,nm);
end

L = find( ((all_list(:,3) == 1) | ...
(all_list(:,3) == m)) & ...
(all_list(:,2) >= 2) & ...
(all_list(:,2) <= (n-1)));
nl=length(L);
if nl>0
fda=fda+sparse(repmat(all_list(L,1),1,3), ...
repmat(all_list(L,1),1,3) + ...
repmat([-1,0,1],nl,1), ...
repmat([1 -2 1],nl,1),nm,nm);
end

% eliminate knowns
rhs=-fda(:,known_list)*A(known_list);
k=find(any(fda(:,nan_list(:,1)),2));

% and solve...
B=A;
B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);

case 4
% Spring analogy
% interpolating operator.

% list of all springs between a node and a horizontal
% or vertical neighbor
hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
hv_springs=[];
for i=1:4
hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
end

% delete replicate springs
hv_springs=unique(sort(hv_springs,2),'rows');

% build sparse matrix of connections, springs
% connecting diagonal neighbors are weaker than
% the horizontal and vertical springs
nhv=size(hv_springs,1);
springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
repmat([1 -1],nhv,1),nhv,nm);

% eliminate knowns
rhs=-springs(:,known_list)*A(known_list);

% and solve...
B=A;
B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;

case 5
% Average of 8 nearest neighbors

% generate sparse array to average 8 nearest neighbors
% for each nan element, be careful around edges
fda=spalloc(n*m,n*m,size(nan_list,1)*9);

% -1,-1
L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1)); 
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% 0,-1
L = find(nan_list(:,3) > 1);
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% +1,-1
L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% -1,0
L = find(nan_list(:,2) > 1);
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% +1,0
L = find(nan_list(:,2) < n);
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% -1,+1
L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m)); 
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% 0,+1
L = find(nan_list(:,3) < m);
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% +1,+1
L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
nl=length(L);
if nl>0
fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
repmat([1 -1],nl,1),n*m,n*m);
end

% eliminate knowns
rhs=-fda(:,known_list)*A(known_list);

% and solve...
B=A;
k=nan_list(:,1);
B(k)=fda(k,k)\rhs(k);

end

% all done, make sure that B is the same shape as
% A was when we came in.
B=reshape(B,n,m);

end
% ====================================================
%      end of main function
% ====================================================
% ====================================================
%      begin subfunctions
% ====================================================
function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans
%   themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element
%      nan_list(i,2) == row index of i'th nan element
%      nan_list(i,3) == column index of i'th nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%      
%      For example, talks_to = [-1 0;0 -1;1 0;0 1]
%      means that each node talks only to its immediate
%      neighbors horizontally and vertically.
% 
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
% use the definition of a neighbor in talks_to
nan_count=size(nan_list,1);
talk_count=size(talks_to,1);

nn=zeros(nan_count*talk_count,2);
j=[1,nan_count];
for i=1:talk_count
nn(j(1):j(2),:)=nan_list(:,2:3) + ...
repmat(talks_to(i,:),nan_count,1);
j=j+nan_count;
end

% drop those nodes which fall outside the bounds of the
% original array
L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m); 
nn(L,:)=[];

% form the same format 3 column array as nan_list
neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];

% delete replicates in the neighbors list
neighbors_list=unique(neighbors_list,'rows');

% and delete those which are also in the list of NaNs.
neighbors_list=setdiff(neighbors_list,nan_list,'rows');

else
neighbors_list=[];
end

end


































 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> make_odd.m
%-----> Included version last modified:23-Apr-2025 14:05:44
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function Arr = make_odd(Arr,Down)

%which way?
if exist('Down','var'); Shift = -1; else Shift = 1; end

%integer?
for iEl=1:1:numel(Arr);  if ~isinteger(Arr(iEl)); Arr(iEl) = round(Arr(iEl)); end

%odd number?
NotOdd = find(mod(Arr,2) == 0);

Arr(NotOdd) = Arr(NotOdd) + Shift;

end

end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> nph_sinefit.m
%-----> Included version last modified:10-Aug-2022 13:07:24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 

% nph_sinefit.m
%
% Created based on a really simple sine-fitting method I saw online, simply
% using the \ operator. It's witchcraft, but it seems to work.

% ALSO: added the ability to add weightings for the fit (turns out this was
% quite hard and involved lots of matrix maths):
% [yfit,F] = nph_sinefit(x,y,period_estimate,'weights',w);
%
% EDIT: Added the ability to fit multiple periods at once. Simply set
% period_estimate to be a vector of the periods you want to fit.
% Tis results in an output matrix that has all the factors Acos() and Bsin()
% for each period estimate, i.e. Acos(t1) + Bsin(t1) + Ccos(t2) + Dsin(t2) + ... + E.
% This has meant I've had to change to output structure to make it easier to
% reconstruct the signal afterwards. The new structure looks like, for a
% 3-period fit for example:
% F = [ A     B     0
%       C     D     0
%       E     F     G]; where G is the signal mean part.
%
% For a single-period fit, the output still looks like:
% F = [ A     B     C]; where C is the signal mean part.
%
% Have fun!
%
%
%
% INPUTS:
%
% x - the coordinates of the corresponding y values. This allows for
% irregularly-spaced data!
% y - vector of signal data at the points in x to be sine-fitted.
% period_estimate - estimate of period to which to fit a sinusoid, in
% whatever units x was in.
%
% OUTPUTS:
%
% OUT - the fitted sinusoidal wave, equal to the length of IN and evaluated
% at the same points specified in IN.
%
% F - 3-element vector of coefficients [A B C], where the fit is defined
% as:               yfit = C + A*cos(2PIx/t) + B*sin(2PIx/t);
% where C is obv the mean "offset" of the signal (DC) and t is the estimated period.
%
% If you're interested, phase is computed as phi = atan2(A,B), or atan2(F(1),F(2))
%
% [yfit,F] = nph_sinefit(x,y,period_estimate);
%

function varargout = nph_sinefit(x,IN,period_estimate,varargin)

if any(strcmpi(varargin,'weights'))
method = 'weights';
else
%     if any(strcmpi(varargin,'weights2'))
%     method = 'weights2';
%     else
method = 'backslash';
%     end
end

% first, linearise to column vectors
origsize = size(IN);
IN = IN(:); x = x(:); period_estimate = period_estimate(:);

% then cope with NaNs by just ignoring them:
nanlocs = isnan(x) | isnan(IN);
% finiteinds = intersect(find(~isnan(x)),find(~isnan(IN)));
xfinite = x(~nanlocs); % make a list of x values that are non-nan
y = IN(~nanlocs);


switch method
% BACKSLASH METHOD ====================================================
case 'backslash'

% % % % %         % create fit matrix:
% % % % %         X = ones(length(y),3); % three part fit, FIT = A*cos(2*pi*t/period) + B*sin(2*pi*t/period) + C;
% % % % %         X(:,1) = cos((2*pi*xfinite) ./ period_estimate);
% % % % %         X(:,2) = sin((2*pi*xfinite) ./ period_estimate);

% NEW: allow for simultaneous multi-component fit, like so:
% A*cos(2*pi*t/periods(1)) + B*sin(2*pi*t/periods(1)) ...
% C*cos(2*pi*t/periods(2)) + D*sin(2*pi*t/periods(2)) ...
% ... + E;

% why not use the outer product to assemble the sines and
% cosines... (nice!)
cosparts = cos((2*pi*xfinite) * (1./period_estimate)');
sinparts = sin((2*pi*xfinite) * (1./period_estimate)');
oneparts = ones(size(xfinite));

% assemble them:
X = ones(length(y),(2*length(period_estimate))+1);
sincos = 1;
for i = 1:2:(2*length(period_estimate))
X(:,i)   = cosparts(:,sincos);
X(:,i+1) = sinparts(:,sincos);
sincos = sincos+1;
end
X(:,end) = oneparts;

% collapse the fit using matrix inversion (the witchcraft)
F = X \ y;

% % % % % %         THIS IS SLOWWWW USING MATLAB'S INBUILT FITTING
% % % % % %         % WEIGHTS METHOD ======================================================
% % % % % %     case 'weights'
% % % % % %
% % % % % %         w = varargin{find(strcmpi(varargin,'weights'))+1};
% % % % % %         w = w(:); % column vector
% % % % % %
% % % % % %         % Fit: 'untitled fit 1'.
% % % % % %         [xData, yData, weights] = prepareCurveData( x, y, w );
% % % % % %
% % % % % %         % Set up fittype and options.
% % % % % %         ft = fittype( 'a*cos(2*pi*x./360) + b*sin(2*pi*x./360) + c', 'independent', 'x', 'dependent', 'y' );
% % % % % %         opts = fitoptions( 'Method', 'NonLinearLeastSquares' );
% % % % % %         opts.Display = 'Off';
% % % % % %         opts.StartPoint = [sqrt(2) sqrt(2) 0];
% % % % % %         opts.Weights = weights;
% % % % % %
% % % % % %         % Fit model to data.
% % % % % %         [fitresult, ~] = fit( xData, yData, ft, opts );
% % % % % %
% % % % % %         F = nan(3,1);
% % % % % %         F(1) = fitresult.a;
% % % % % %         F(2) = fitresult.b;
% % % % % %         F(3) = fitresult.c;

% WEIGHTS METHOD ======================================================
case 'weights'

w = varargin{find(strcmpi(varargin,'weights'))+1};
w = w(:); % column vector

% remove any duds found earlier
w = w(~nanlocs);

% Very slightly faster in singles:
x = single(x); y = single(y); w = single(w);

% Use the sneaky outer product again :) tiny bit faster than
% repmat I think.
% W = ones(size(w)) * w' .* eye(length(w));
% ^ still rather slow :(
% EDIT: I've removed the above but want to leave it because of its
% simplicity. It's faster not to generate an NxN matrix of
% weightings computationally, but when writing this up analytically I'd use it.

% NEW: allow for simultaneous multi-component fit, like so:
% A*cos(2*pi*t/periods(1)) + B*sin(2*pi*t/periods(1)) ...
% C*cos(2*pi*t/periods(2)) + D*sin(2*pi*t/periods(2)) ...
% ... + E;

% why not use the outer product to assemble the sines and
% cosines... (nice!)
cosparts = cos((2*pi*xfinite) * (1./period_estimate)');
sinparts = sin((2*pi*xfinite) * (1./period_estimate)');
oneparts = ones(size(xfinite));

% assemble them:
X = ones(length(y),(2*length(period_estimate))+1);
sincos = 1;
for i = 1:2:(2*length(period_estimate))
X(:,i)   = cosparts(:,sincos);
X(:,i+1) = sinparts(:,sincos);
sincos = sincos+1;
end
X(:,end) = oneparts;

% Matrix Inversion!
% F = inv(X' * W * X) * X' * W * y;
% F =    (X' * W * X) \ (X' * W * y); % <- this one is about 50% faster
F =    (X' * (repmat(w,1,size(X,2)) .* X)) \ (X' * (w .* y)); % <- slightly faster still!
%         F =    (X * (repmat(w,1,size(X,2)) .* X)') \ (X' * (w' .* y')); % <- slightly faster still!

% Return to double:
F = double(F);

end
% =========================================================================

% Reshape F ino something for manageable for multiple-period fits:
final_form = zeros(length(period_estimate),3);

Fr = reshape(F(1:(2*length(period_estimate))),[2 length(period_estimate)])';
final_form(1:length(period_estimate),1:2) = Fr;
final_form(length(period_estimate),3) = F(end);

F = final_form;

% Evaluate the fitted sinusoid: (use the original x-range, NaNs included)
yfit = zeros(size(x));
for i = 1:length(period_estimate)
yfit = yfit + ...
F(i,1).*cos((2.*pi.*x) ./ period_estimate(i)) + ...
F(i,2).*sin((2.*pi.*x) ./ period_estimate(i)) + ...
F(i,3);
end

% p = 0;
% yfit = F(end);
% for i = 1:2:(size(X,2)-1)
%     p = p + 1;
%     yfit = ...
%         yfit + ...
%         F(i)*cos((2*pi*x) ./ period_estimate(p)) + ...
%         F(i+1)*sin((2*pi*x) ./ period_estimate(p));
% end

% yfit = F(3) + F(1)*cos((2*pi*x) ./ period_estimate) + F(2)*sin((2*pi*x) ./ period_estimate);

yfit = reshape(yfit,origsize);

% if everything was NAN, output all NANs.
if all(isnan(IN)) || all(isnan(x)) || isempty(IN) || isempty(x)
yfit = nan(size(yfit));
F = nan(size(F));
end

% output(s): might add more in future
switch nargout
case {1,0}
varargout{1} = yfit;
case 2
varargout{1} = yfit;
varargout{2} = F;
end

%% PLOTTING
if nargin > 3
if any(strcmpi(varargin,'plot'))
figure; hold all; grid on;
hold on; plot(x,IN,'.k');
hold on; plot(x,yfit,'.b');
end
end
end







%
%
%
% return
% %
%
% %% least squares sinusoid fitting:
%
% p_est = 25;
%
%  t = (1:100)';
%  X = ones(100,3);
%  X(:,2) = cos((2*pi)/p_est*t);
%  X(:,3) = sin((2*pi)/p_est*t);
%  y = 2*cos((2*pi)/p_est*t-pi/4)+randn(size(t));
%  y = y(:);
%  beta = X\y;
%  yhat = beta(1)+beta(2)*cos((2*pi)/p_est*t)+beta(3)*sin((2*pi)/p_est*t);
%
%  figure;
%  plot(t,y,'b');
%  hold on
%  plot(t,yhat,'r','linewidth',2);
%
%
%  %% try sine fitting wind:
%
% u = reshape(squeeze(HWD.Data.u(5,:,100:110)),[1 24*11]);
%
% IN = u;
%
% IN = IN(:); % linearise to column vector
% t = (1:length(IN))'; % must all be column vectors
% % but does it actuslly have to be regularly spaced?! Investigate!
%
% % use only non-nans: (crucial step!)
% finiteinds = ~isnan(IN);
% t = t(finiteinds);
% IN = IN(finiteinds);
%
% p_est = 24; % period estimate, elements, must be regularly spaced
%
%  X = ones(length(IN),3);
%  X(:,2) = cos((2*pi)/p_est*t);
%  X(:,3) = sin((2*pi)/p_est*t);
%
% %  y = 2*cos((2*pi)/p_est*t-pi/4)+randn(size(t));
% %  y = y(:);
%
% B = X \ IN;
%
% F = B(1)+B(2)*cos((2*pi)/p_est*t)+B(3)*sin((2*pi)/p_est*t);
%
% figure; plot(t,IN);
% hold on; plot(t,F,'r','linewi',2)
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> textprogressbar.m
%-----> Included version last modified:30-Dec-2024 15:50:10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
function textprogressbar(c)
% This function creates a text progress bar. It should be called with a 
% STRING argument to initialize and terminate. Otherwise the number correspoding 
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate 
%                       Percentage number to show progress 
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 25;   %   The total number of dots in a progress bar

%% Main 

if isempty(strCR) && ~ischar(c),
% Progress bar must be initialized with a string
error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c),
% Progress bar - initialization
fprintf('%s',c);
strCR = -1;
elseif ~isempty(strCR) && ischar(c),
% Progress bar  - termination
strCR = [];  
fprintf([c '\n']);
elseif isnumeric(c)
% Progress bar - normal progress
c = floor(c);
percentageOut = [num2str(c) '%%'];
percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
nDots = floor(c/100*strDotsMaximum);
dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
strOut = [percentageOut dotOut];

% Print it on the screen
if strCR == -1,
% Don't do carriage return during first run
fprintf(strOut);
else
% Do it during all the other runs
fprintf([strCR strOut]);
end

% Update carriage return
strCR = repmat('\b',1,length(strOut)-1);

else
% Any other unexpected input
error('Unsupported argument type');
end

end
 
 
 
 
 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> AUTOMERGED: 
%-----> nph_ndst.m
%-----> Included version last modified:22-Apr-2023 11:53:13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 

% nph_ndst.m
% A 1, 2, 3, and 4-D application of the Stockwell Transform, developed by
% Neil Hindley, University of Bath, 2017.


%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ST = nph_ndst(IN,scales,point_spacing,c,varargin)
%
% *nomenclature: N = number of dimensions, n = number of elements in a
% particular dimension. M = the desired number of frequency "scales".
%
% IN [REQUIRED]: N-D vector/matrix to be transformed. 1, 2, and 3-D supported.
%
% scales: [1x1 scalar | NxM vector | 1xN cell structure]. Choose from:
%   - [1x1 scalar] If scales is a scalar [M], the dominant [M] scales are found in the input
% data using the FFT, and then ONLY this many scales are computed in the NDST.
% This method is very quick and is the one now I use most. This is not compatible
% with the 'full' output version. You can also specify min and max wavelengths
% to be analysed (in real units given by point_spacing).
%   - [NxM vector] containing the exact scale COMBINATIONS to be analysed.
% I would normally use this by taking the output ST.scales (which is
% usually going to be NxM if you used the scalar or vector option) from a
% previous NDST run and inputting this directly for scales. This is useful
% if you want to analyse lots of different data with the exact same
% frequencies for each.
%   -  [1xN cell structure], where N is dimensions, each containing a 1xM
% length vector of integer scales with which to analyse the relevant
% dimenion. These correspond to integer fractions of the maximum possible
% wavelength, ie the full length of the input series, up to the Nyquist.
% For 1-D ST, just entering a vector is fine, no need for a cell.
% Scales should be anything withing the range [-(n/2-1):-1 1:1:(n/2-1)],
% where n is the number of elements in a particular dimension. Defaults are
% roughly [-(n/3):(n/3)] for each dimension for the 2DST and 3DST, and
% 1:Nyquist (1:n-1) for the 1DST.

% point_spacing [1xN vector]: containing the real physical separation
% between points, be it in time or distance etc. Data input must have this
% regular sampling for measured wavelengths to be accurate.
%
% c [1xN vector]: scaling parameter for the Guassian windows. See Hindley et al., AMT
% (2016) for details on the effect of this. Would generally recommend
% setting c=0.25 for all dimensions in my work, but up to you.
%
% [OPTIONAL ARGUMENTS]
%
% 'quick' (default) | 'full' - choose whether to output the full
% 2N dimensional S-transform complex output. For 1-D, default is 'full', as
% it's only small, but for 2- and 3-D, default is 'quick', which provides
% the dominant spectral component at each location in the data. Remember,
% the 3-D output is 6-D, which can be pretty massive, so use with caution.
%
% 'zeromean' (default) | 'nozeromean' - choose whether to compute the
% NDST on the zero-mean signal or not. Should really be doing this as a
% rule anyway.
%
% 'boost', EDIT: new hilbert boosting is now done by default!
% 'boost' (default for 2D/3D) | 'noboost' (default for 1D) - choose
% whether to apply the Hilbert Boosting method described below to try to
% get a better estimate of wave-packet amplitudes, which are typically
% underestimated by the 2D and 3D ST, but less so by the 1D ST.

% [OUTPUTS]
%
% ST.ST - (m{1,2,3} x n{1,2,3}) Stockwell Tranform complex cospectrum.
%
% ST.C - N-dimensional complex cospectrum of the dominant (largest spectral
% amplitude) frequency at each location. This will either be boosted/not
% boosted depending on whether you want to boost the amplitude to cope with
% the packet-like nature of waves.
%
% ST.A - abs() of ST.C, the instantaneous amplitude of the dominant
% frequency at each location.
%
% ST.R - real() part of ST.C, provides "reconstruction" of the wavefield as
% the NDST saw it.
%
% ST.F{1,2,3} - dominant N-dimensional spatial frequencies at each location
% (inverse of distance or time, no 2*pi)
%
% ST.freqs - 1xN cell object of the frequencies that were analysed,
% computed from the scales that were inputted.
%
% ST.HA - Absolute part of the Hilbert Transform of the psuedo-bandpassed
% data. All the Gaussian windows are combined and applied to the FFT
% spectrum, justa  rough bandpassed filter, then the Hilbert transform is
% applied to this to find the phase-invariant amplitude at each location.
% This is superior because it is the amplitude of only those wavelengths
% that we have considered in the S-transform.
%
% ST.HR - Real part of the above.
%
% 
% REMOVED:
% ST.BoostFactor - factor by which ST.C has been boosted by the Hilbert
% Boosting method.
%
%


function ST = nph_ndst(IN,varargin)

IN = squeeze(IN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine if 1D 2D 3D input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz = size(IN);
sz = sz(sz ~= 1);
type = length(sz);

% fix the annoying 1xN versus Nx1 problem:
if type == 1
IN = reshape(IN,[1 length(IN)]);
end
osz = size(IN); % original corrected size in


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% VARARGIN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Split varargin into inputs and string option flags:
inputs = {};
options = {};
for v = 1:length(varargin)
if ischar(varargin{v}) || isstring(varargin{v})  % options must be a character string
options{v} = varargin{v};
else
inputs{v} = varargin{v};
end
end

% decide whether to allow zeros and 0.5s etc in scales
allownonintegerscales = 0;

% Specify some defaults for options:
switch type
case 1
fullflag = 1;
boostflag = 0;
zeromeanflag = 1;
otherwise
fullflag = 0;
boostflag = 1;
zeromeanflag = 1;
end


% Default SCALES:
default_scales = cell(1,type);
for i = 1:type
if i < type
default_scales{i} = -15:15;
else
default_scales{i} = 1:15;
% just positive freqs for the last dimension, save repeating yourself.
end
end

% Default POINT SPACING:
default_point_spacing = ones(1,type);

% Default scaling parameter C:
default_c = 0.25 * ones(1,type);

scales = varargin{1};
point_spacing = varargin{2};
c = varargin{3};

% Assign defaults if they're missing:
switch length(inputs)
case 0 % no inputs
scales = default_scales;
point_spacing = default_point_spacing;
c = default_c;
case 1 % just scales
scales = inputs{1};
point_spacing = default_point_spacing;
c = default_c;
case 2 % just scales and point_spacing:
scales = inputs{1};
point_spacing = inputs{2};
c = default_c;
otherwise % scales, point_spacing and c:
scales = inputs{1};
point_spacing = inputs{2};
c = inputs{3};
end


% FULL 2-D, 4-D, 6-D, 8-D etc S-TRANSFORM OBJECT
if any(strcmpi(options,'full'))
fullflag = 1;
else
% else only output the full spectrum for the 1DST as default
switch type
case 1
fullflag = 1;
otherwise
fullflag = 0;
end
end

% ZERO MEAN
% Determine if we want to use zero-mean (as a rule we should really get
% in the habit of giving signals zero-mean before the NDST. This means
% of course that any underlying trend will manifest as a signal, but it
% results in fairer recovery of stuff which sits upon that trend.
if any(strcmpi(options,'nozeromean'))
zeromeanflag = 0;
else
zeromeanflag = 1;
end

% HILBERT BOOSTING
% Determine if we want to use the "Hilbert Boosting" method to try to
% more accurately measure the amplitudes of wave packets, which are
% always underestimated by the ST.
if any(strcmpi(options,'noboost'))
boostflag = 0;
else
boostflag = 1;
end


% % GUIDED FOURIER MODE
% % new exciting mode, activated by specifying a scalar number for the number
% % of scales to use for higher dimensional S-transforms:
% 
% if isnumeric(scales) && type ~= 1 && length(scales) == 1
%     guidedfourierflag = 1;
%     fullflag = 0; % haven't yet sorted out a full spectrum with guided fourier mode
% else
%     guidedfourierflag = 0;
% end

% MIN AND MAX WAVELENGTHS FOR GUIDED FOURIER MODE
% so guided fourier mode is great and all, but you lose control over what
% frequencies you want to study. So if you want to only look at the large
% scale stuff, ignoring small scale, you can't. This is why I've added the
% ability to specify min amd max WAVELENGTHS for use in guided fourier mode.
% THESE ARE ABS VALUES, NOT CURRENTLY SUPPORTING NEGATIVES.
maxwavelengthsflag = 0;
if any(strcmpi(options,'maxwavelengths'))
maxwavelengthsflag = 1;
maxwavelengths = abs(varargin{(find(strcmpi(varargin,'maxwavelengths'))+1)});
end
minwavelengthsflag = 0;
if any(strcmpi(options,'minwavelengths'))
minwavelengthsflag = 1;
minwavelengths = abs(varargin{(find(strcmpi(varargin,'minwavelengths'))+1)});
end


%% PARSE INPUT SIZES AND SCALES ===========================================

% We need to work out what input scales is.

% The shape of the scale input determines the method used.

% There are three options:

% 1. 'scalar'
% A scalar value for N triggers guidedfourier mode. In this mode, only the
% dominant N feqs in the whole input data will be analysed.
% If you want to analyse with just one scale, do [X X] of the same
% scale (like when drawing only one contour line in contour.m) to trigger
% option 2 below.

% 2. 'vector'
% 1xN, 2xN, 3xN or 4xN vectors of scales. This is is a list of scale
% combinations for the NDST to analyse. It will ONLY analyse for these
% EXACT combinations, as if they were selected by guided fourier mode.
% This is useful if you want to find some good scales first, then
% analyse several bits of data all using these exact same scales. Also
% triggers guided fourier mode, but we don't compute new scales we just run
% it that way.

% 3. 'cell'
% Cell object of {1xM,1xN,1xP} etc for N-D inputs.
% This will do every other scale for each scale in each dimension,
% so the "full" output would be XxYxMxN for 2-D or XxYxZxMxNxP for 3-D.
% This is potentially quite slow, but it's a more complete approach as you
% get "maps" of the spectral properties at each location in X,Y,Z. This is
% closest to the original implementation.


ssz = size(scales);
scalesformat = 'unknown'; % starts off unknown.

% OPTION 1 - Scalar number of scales
if isnumeric(scales) && numel(scales) == 1
scalesformat = 'scalar';
end

% OPTION 2 - Vector of scale combinations
if isnumeric(scales) && numel(scales) ~= 1 && ssz(1) == type
scalesformat = 'vector';
end

% OPTION 3 - Cell array of range of scales for each dimension
if iscell(scales) && length(scales) == type
scalesformat = 'cell';
end

% FLAGS SPECIFICATION
% Now decide what flags to do based on the scales input:
switch scalesformat
case 'scalar'
guidedfourierflag = 1;
presetscales = 0;
if type == 1
error('Error: Scalar scale input (guided fourier mode) not currently supported for 1-D. Please enter scales manually as 1:N, where N less than the length of the input data.')
end
case 'vector'
guidedfourierflag = 0;
presetscales = 1;
case 'cell'
guidedfourierflag = 0;
presetscales = 0;

end

% ERROR CHECKING
switch scalesformat
case 'scalar'

case 'vector'
if ssz(1) ~= type
error(['Error: Expected scales to be ' num2str(type) '-by-N vector for ' num2str(type) '-D input.'])
end
case 'cell'
if ssz(2) ~= type
error(['Error: Expected scales to be 1-by-' num2str(type) ' cell for ' num2str(type) '-D input.'])
end
otherwise
error(['Error: Expected scales to be a scalar, or a ' num2str(type) '-by-N vector or 1-by-' num2str(type) ' cell.'])
end

% check c:
if numel(c) ~= type
error(['Error: For ' num2str(type) '-D input data, c must be a 1x' num2str(type) ' vector.']) 
end


% WARNING FOR 3DST AND 4DST
% Gonna run out of memory pretty fast if we try and output some full
% 6-D or 8-D structures...
if fullflag && strcmp(scalesformat,'cell')
switch type
case 3
if ~isempty(inputs) % if scales are inputted
s = inputs{1};
w = whos('IN');
Mb = (w.bytes ./ 1024 ./ 1024);
switch class(IN)
case 'double'
Mb = Mb ./ 2;
end
Mb = Mb .* length(s{1})*length(s{2})*length(s{3});
warning(['Outputing a full 6-D S-transform object at single precision. This might require up to ' num2str(Mb) ' Mb of memory.'])
end
case 4
if ~isempty(inputs) % if scales are inputted
s = inputs{1};
w = whos('IN');
Mb = (w.bytes ./ 1024 ./ 1024);
switch class(IN)
case 'double'
Mb = Mb ./ 2;
end
Mb = Mb .* length(s{1})*length(s{2})*length(s{3})*length(s{4});
warning(['Outputing a full 8-D S-transform object at single precision. This might require up to ' num2str(Mb) ' Mb of memory.'])
warning('If you''re sure, press any key to continue.')
pause
end
end
end


% PARSE AND SORT CELL RANGES OF SCALES
switch scalesformat
case 'cell'
% SORT SCALES INTO ASCENDING ORDER TOO!!!!
for t = 1:type
scales{t} = sort(scales{t});
end
% REMOVE DUPLICATES AND ZEROS
for t = 1:type
scales{t} = unique(scales{t},'stable');
end
% Decide if non-integer scales are allowed:
if ~allownonintegerscales
for t = 1:type
scales{t} = unique(fix(scales{t}),'stable');
scales{t} = scales{t}(scales{t} ~= 0);
end
%             for i = 1:length(scales) % integer scales, not equal to zero.
%                 sc = scales{i};
%                 sc = unique(fix(sc),'stable'); % unique sorts it by default :)
%                 scales(i) = {sc(sc ~= 0)};
%             end
end
end


%% GENERATE SCALES FOR GUIDED FOURIER MODE, IF ENABLED ====================
% alright here's where we use the FFT to find the top XXXX frequencies
% present in the input data, then work out what their scales would be.
% Hold my beer...

switch scalesformat

case 'scalar' % guided fourier mode - choose freqs in house!

nfreqs = scales;

if zeromeanflag
IN = IN - mean(IN(:));
end

% take FFT:
F = fftn(IN);
ab = abs(F(:));
im = imag(F(:));

% sort by absolute spectral power
[ab,ib] = sort(ab,'descend');

% also rearrange the imaginary comps by this sorting:
im = im(ib);

% exclude the DC components with imag parts == 0
ab = ab(im ~= 0);
ib = ib(im ~= 0);

% now reshape:
% this should always work - after the zeros are taken out there should
% always be an even number of complex conjugate pairs remaining.
% note: you need the transpose ' here due to the way reshape re-lists things.
abr = reshape(ab',2,length(ab)/2);
ibr = reshape(ib',2,length(ib)/2);

% Make a coord system in fft space
sz = size(F);
v = struct;
for n = 1:type

switch iseven(sz(n))
case 1
N = (sz(n)/2)-1;
v(n).vec = ifftshift([0 -N:N]);
case 0
N = (sz(n)-1)/2;
v(n).vec = ifftshift(-N:N);
end

end

% what were the scales that related to these locations?
ii = struct;
switch type
case 1
ii(1).i = ind2sub(size(F),ibr(1,:));
case 2
[ii(1).i,ii(2).i] = ind2sub(size(F),ibr(1,:));
case 3
[ii(1).i,ii(2).i,ii(3).i] = ind2sub(size(F),ibr(1,:));
case 4
[ii(1).i,ii(2).i,ii(3).i,ii(4).i] = ind2sub(size(F),ibr(1,:));
end

% RESET SCALES:
scales = cell(1,type);
wavelengths = cell(1,type);
physical_dims = nan(1,type);
goodinds = ones(1,length(ii(1).i));

% LIMIT TO MIN/MAX SCALES, NON-ZERO SCALES:
for n = 1:type

scales{n} = v(n).vec(ii(n).i);

% remove zero scales:
goodinds = all(cat(1,goodinds,scales{n} ~= 0));

% covert to wavelengths if needed:
physical_dims(n) = point_spacing(n) * sz(n);
wavelengths{n} = physical_dims(n) ./ scales{n};

% apply MIN wavelength cutoff:
if minwavelengthsflag
goodinds = all(cat(1,goodinds,abs(wavelengths{n}) >= abs(minwavelengths(n))));
end

% apply MAX wavelength cutoff:
if maxwavelengthsflag
goodinds = all(cat(1,goodinds,abs(wavelengths{n}) <= abs(maxwavelengths(n))));
end

end

% Apply this externally to the above loop so that it's the same for all
% dimensions:
for n = 1:type
scales{n} = scales{n}(goodinds);
wavelengths{n} = wavelengths{n}(goodinds);
end

% if the user has asked for more freqs than there are available, limit it:
for n = 1:type
if nfreqs > length(scales{n})
warning(['Too many frequencies requested. Limiting to ' num2str(length(scales{n})) '.'])
nfreqs = length(scales{n});
end
end

% Finally, select the top NFREQS from our formatted scales:
for n = 1:type
scales{n} = scales{n}(1:nfreqs);
wavelengths{n} = wavelengths{n}(1:nfreqs);
end

% % % %     % Try this: for the collapsed spectrum we care about the order
% % % %     % (slightly) that the frequencies are called. So let's sort them in
% % % %     % order of lowest to highest scales. Not sure if it'll make a
% % % %     % difference but lets's see.
% % % %     scalemag = reshape([scales{1:type}],type,nfreqs);
% % % %     scalemag = sqrt(sum(scalemag.^2,1));
% % % %     [~,ord] = sort(scalemag,'ascend');
% % % %     
% % % %     for n = 1:type
% % % %         scales{n} = scales{n}(ord);
% % % %         wavelengths{n} = wavelengths{n}(ord);
% % % %     end

% Convert to vectors:
scales_vec = zeros(type,nfreqs);
wavelengths_vec = zeros(type,nfreqs);
for n = 1:type
scales_vec(n,:) = scales{n};
wavelengths_vec(n,:) = wavelengths{n};
end
scales = scales_vec;
wavelengths = wavelengths_vec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Still guided fourier mode, but we've specified the scale combinations
% that are to be analysed. Now we just need to apply our min/max
% wavelengths cutoffs and other stuff to make sure we're consistent.
case 'vector'

nfreqs = size(scales,2);
szz = size(IN); szz = szz(szz ~= 1);

% Convert to wavelengths:
wl = repmat(szz',1,nfreqs);
ps = repmat(point_spacing',1,nfreqs);
wavelengths = (wl ./ scales) .* ps;

% start with full indeces:
goodinds = ones(1,nfreqs);

for n = 1:type

% remove zero scales:
goodinds(scales(n,:) == 0) = 0;

% apply MIN wavelength cutoff:
if minwavelengthsflag
minwlinds = abs(wavelengths(n,:)) >= abs(minwavelengths(n));
goodinds(~minwlinds) = 0;
end

% apply MAX wavelength cutoff:
if maxwavelengthsflag
maxwlinds = abs(wavelengths(n,:)) <= abs(maxwavelengths(n));
goodinds(~maxwlinds) = 0;
end

end

scales = scales(:,logical(goodinds));
wavelengths = wavelengths(:,logical(goodinds));

%         disp('')


end % end switching between scales input formats


%% ASSEMBLE ALL SCALES INTO A MATRIX ======================================
% so that you can do one for loop below for all the IFFTNs
% man this is a ball ache % wait - USE GRIDS!!!!!!!!

% This approach feels like a bit of a sludge, but it enables you to
% specify individual frequency combinations. Currently, if you specify a
% scale of 10 for the first dimension, all scales in further dimensions
% will be tried with the scale 10 for the first dimension, regardless of
% whether or not they were important, and you couldn't get round this.
% Changing this will allow for guided fourier mode to be enabled, and
% should also produce a slight speed up from not using nested loops.

% new new approach: do all this by the type of scales input.

switch scalesformat
case {'scalar','vector'}
% this will have engaged the guided fourier mode above

%         scale_inds = repmat(1:size(scales,2),type,1);
%         scale_inds = 

scale_inds = scales;
totalnumscales = size(scales,2);
scale_inds = repmat(1:totalnumscales,type,1);

allscales = scales; % that was easy

% the allscales thing is only used for computing the real
% frequencies later. We then use the scale_inds to call these
% frequencies in the ST computations below.

%         
% 
%         totalnumscales = size(scales,2);
%         allscales = nan(type,totalnumscales);
%         
%         for n = 1:type
%             allscales(n,:) = 1:length(scales{n});
%         end
%                 
%     case 'vector'
%         % already defined into a neat line of scale combinations
%         totalnumscales = ssz(2);
%         allscales = scales; % that was easy

case 'cell'
% takes a bit more work. Use grids to pour in all the possible scale
% combinations from the range of scales you've specified for each
% dimension.

numscales = zeros(1,type);
for t = 1:type
numscales(t) = length(scales{t});
end

totalnumscales = prod(numscales);
%         scale_inds = repmat(1:totalnumscales,type,1);

switch type
case 1
scale_inds = 1:totalnumscales; % that was easy
case 2
[S1,S2] = ndgrid(1:numscales(1),1:numscales(2));
scale_inds = cat(2,S1(:),S2(:))';
case 3
[S1,S2,S3] = ndgrid(1:numscales(1),1:numscales(2),1:numscales(3));
scale_inds = cat(2,S1(:),S2(:),S3(:))';
case 4
[S1,S2,S3,S4] = ndgrid(1:numscales(1),1:numscales(2),1:numscales(3),1:numscales(4));
scale_inds = cat(2,S1(:),S2(:),S3(:),S4(:))';
end


% assign an "allscales" structure for consistency with the other
% scales formats during the NDST processing:

allscales = zeros(type,totalnumscales);
for t = 1:type
allscales(t,:) = scales{t}(scale_inds(t,:));
end



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% INITIALISE OUTPUTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise:
ST = struct;

% record inputs:
ST.IN = IN;
ST.scales = scales;
ST.point_spacing = point_spacing;
ST.c = c;

%  Compute frequencies for output:
switch scalesformat
case {'scalar','vector'}
physical_lengths = repmat(sz(:) .* point_spacing(:),1,totalnumscales);        
for t =  1:type
ST.freqs = allscales ./ physical_lengths;
end
case 'cell'
for t = 1:type
ST.freqs{t} = scales{t} ./ (sz(t)*point_spacing(t));
end
end

% RECORD ANY FLAGS USED
ST.AmplitudeBoosting = boostflag;
ST.GuidedFourierMode = guidedfourierflag;
if ~isempty(options)
ST.Options = options;
end

% INITIALISE MAIN OUTPUT FIELDS:
if fullflag
switch scalesformat
case {'scalar','vector'}
ST.ST = zeros([sz totalnumscales],'single');
case 'cell'
scalerangesize = 1:type;
for t = 1:type
scalerangesize(t) = length(scales{t});
end
ST.ST = zeros([sz scalerangesize],'single');
end
end

% now the rest...
for t = 1:type
ST.C = zeros(osz,'single');
ST.(['F' num2str(t)]) = zeros(osz,'single');
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% COMPUTING THE NDST 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ZERO-MEAN?
if zeromeanflag
IN = IN - mean(IN(:));
end


%% STEP 1: FFTN ===========================================================

F = single(fftn(IN));


% %% STEP 2: GATHER FFT WISDOM ==============================================
% % gather wisdom for the IFFTN() in the loop below...
% if type > 1 % don't bother for the 1DST
%     ST.fftwisdom = gatherfftwisdom(F,'ifftn');
%     fftw('swisdom',ST.fftwisdom); % Apply the wisdom for SINGLE precision.
% end


%% STEP 3: HILBERT MASK ===================================================

FM = F .* nph_hilbertmask(F);


%% STEP 4: ASSEMBLE GAUSSIAN WINDOWS ======================================
% now fully ND compatible, pre-make gaussian vectors into a structure
% EDIT: this needs updating for guided fourier mode. At the moment we're
% compiling gaussian for every scale, not every scale combo. This is likely
% using up several tens of times more run time than necessary.

% nope!
% % % % % % EDIT: NOW COMPUTING FOR ZERO SCALES if you allow them, just an array of
% % % % % % zeros with a one at the right place in the centre of where the gaussian
% % % % % % would be.

GW = struct;

for i = 1:type

% First, make a coord system in fftshifted space (easier to work in)
switch iseven(sz(i))
case 1
N = (sz(i)/2)-1;
%                 x = -(N+1):N;
x = [0 -N:N];
case 0
N = (sz(i)-1)/2;
x = -N:N;
end

% and ifftshift it:
x = ifftshift(x);


% Now define Gaussian vector storage structure:
switch scalesformat
case {'scalar','vector'}
len = length(scales(i,:));
GW(i).gvecA = nan(len,sz(i));
GW(i).gvecB = nan(len,sz(i));
case {'cell'}
len = length(scales{i});
GW(i).gvecA = nan(len,sz(i));
GW(i).gvecB = nan(len,sz(i));
end

% % % % % %%% TRIED THIS, HATED IT. MAYBE I DID SOMETHING WRONG?
% % % % %         % IF SCALE == 0
% % % % %         if scales{i}(j) == 0
% % % % %
% % % % %             % for zero scales:
% % % % %             % just a line of zeros with ones at the zeroth freq. In the
% % % % %             % other dimensions it'll be a Gaussian, just one element wide
% % % % %             % in this dimensions.
% % % % %             gwA = zeros(size(x));
% % % % %             gwB = zeros(size(x));
% % % % %
% % % % %             gwA(find(x == 0,1,'first')) = 1;
% % % % %             gwB(find(x == 0,1,'last')) = 1;
% % % % %
% % % % %         else


% Now for each scale combination
for j = 1:len

% FOUND THE PROBLEM - here, these shouldn't be the REAL freqs, but
% the integer scales!

% Evaluate Gaussians in the normal way:
switch scalesformat
case {'scalar','vector'}
% normal window:
gwA = exp( (-2*pi^2) * (c(i)/scales(i,j))^2 * (x - scales(i,j)).^2 );
% a mirror image window:
gwB = exp( (-2*pi^2) * (c(i)/scales(i,j))^2 * (x + scales(i,j)).^2 );
case 'cell'
% normal window:
gwA = exp( (-2*pi^2) * (c(i)/scales{i}(j))^2 * (x - scales{i}(j)).^2 );
% a mirror image window:
gwB = exp( (-2*pi^2) * (c(i)/scales{i}(j))^2 * (x + scales{i}(j)).^2 );
end
% % % % %     end
% % % % %

% and assign:
GW(i).gvecA(j,:) = gwA;
GW(i).gvecB(j,:) = gwB;

end

% % % % %     disp('bugger!')
% % % % %     if any(isnan(GW(i).gvecA(j,:))) || any(isnan(GW(i).gvecB(j,:)))
% % % % %
% % % % %         %         return
% % % % %     end
% % % % %


end


%% MAKE A GAUSSIAN WINDOW STORE to find what freqs we computed:
gws = zeros(osz);


%% STEP 5: APPLY THE GAUSSIAN WINDOWS AND IFFTN ===========================

% Because we've pre-assembled all our gaussian windows, the scales that we
% analyse for aren't really needed here, only their indeces in the gaussian
% window structure.
% as long as these indeces match up to the output indeces for the full
% output, then we're all good.
% the only thing we need the actual scales for is the frequency maps that
% are outputted.

% preallocate indexing using cells - brand new very exciting, didn't know
% matlab could do this!!
range = cell(1,type);
for t = 1:type
range{t} = 1:sz(t);
end

% trying to combine loop into one for all STs
for isc = 1:totalnumscales

% find indeces for scales for each dimension:
inds = scale_inds(:,isc)';
indrange = num2cell(inds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COMBINE GAUSSIAN WINDOWS FOR EACH SCALE COMBINATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble gaussian window for this scale combination.
% EDIT: Now generalised for N dimensions!!!!

% first build a cell structure containing inputs for ndgrid.
gvecsA = cell(1,type);
gvecsB = cell(1,type);
for n = 1:type
gvecsA{n} = GW(n).gvecA(inds(n),:);
gvecsB{n} = GW(n).gvecB(inds(n),:);
end

% now collect outputs from ndgrid into another cell:
gwA = cell(1,type);
gwB = cell(1,type);
[gwA{1:type}] = ndgrid(gvecsA{:});
[gwB{1:type}] = ndgrid(gvecsB{:});

% and combine these together to make the ND Gaussian:
Apart = ones(size(gwA{1}));
Bpart = ones(size(gwB{1}));
for n = 1:type
Apart = Apart .* gwA{n};
Bpart = Bpart .* gwB{n};
end

% the Gaussian function is just the sum of these:
gw = reshape(Apart + Bpart,[osz]);

% gaussian window storage for hilbert amplitude later:
gwloc = gw > gws; % (:) fixes weird matrix multiplication in 1D
gws(gwloc) = gw(gwloc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COMPUTE THE ST FOR THIS SCALE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% or voice, to use Stockwell's terminology
FM_voice = ifftn(FM .* gw);

% do the collapsed "rapide" spectrum
loc = abs(FM_voice) > abs(ST.C);
ST.C(loc) = FM_voice(loc);
switch scalesformat
case {'scalar','vector'}
for t = 1:type
ST.(['F' num2str(t)])(loc) = ST.freqs(t,inds(t));
end
case 'cell'
for t = 1:type
ST.(['F' num2str(t)])(loc) = ST.freqs{t}(inds(t));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% AND DO THE FULL SPECTRUM IF REQUIRED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% and, if required, the full ST object:
if fullflag
switch scalesformat
case {'scalar','vector'}
r = [range {isc}];
%                 switch type
%                     case 1
%                         ST.ST(:,isc) = FM_voice;
%                     case 2
%                         ST.ST(:,:,isc) = FM_voice;
%                     case 3
%                         ST.ST(:,:,:,isc) = FM_voice;
%                     case 4
%                         ST.ST(:,:,:,:,isc) = FM_voice;
%                 end
case 'cell'
r = [range indrange];
%                 switch type
%                     case 1
%                         ST.ST(:,i1) = FM_voice;
%                     case 2
%                         ST.ST(:,:,i1,i2) = FM_voice;
%                     case 3
%                         ST.ST(:,:,:,i1,i2,i3) = FM_voice;
%                     case 4
%                         ST.ST(:,:,:,:,i1,i2,i3,i4) = FM_voice;
%                 end
end
% and assign using cells to specify indexing!! Amazing!!
ST.ST(r{:}) = FM_voice;

end % end if fullflag
%     
% %     if isc == 38
%         hold on; pcolor(gw); shat;
%         title(num2str(isc))
%         drawnow;
% %             pause
% %     end


end % done! next scale!!


% 
% switch type
%     case 1 % 1DST
%         
%             % find indeces for scales for each dimension:
%             i1 = scale_inds(1,isc);
%             % for this gaussian window...
%             gwA = GW(1).gvecA(i1,:);
%             gwB = GW(1).gvecB(i1,:);
%             gw = gwA + gwB;
%             % gaussian window storage for hilbert amplitude later:
%             gwloc = gw > gws;
%             gws(gwloc) = gw(gwloc);
%             % Always compute the full 2D spectrum, it's not much :)
%             FM_voice = ifftn(FM .* gw);
%             % insert into S-transform
%             ST.ST(:,i1) = FM_voice;
%             % in any case, do the collapsed rapide spectrum too:
%             loc = abs(FM_voice) > abs(ST.C);
%             ST.C(loc) = FM_voice(loc);
%             ST.F1(loc) = ST.freqs(i1);
%         end
%     case 2 % 2DST
%         for isc = 1:totalnumscales
%             % find indeces for scales for each dimension:
%             i1 = scale_inds(1,isc);
%             i2 = scale_inds(2,isc);
%             % assemble gaussian window
%             [gw1A,gw2A] = ndgrid(GW(1).gvecA(i1,:),GW(2).gvecA(i2,:));
%             [gw1B,gw2B] = ndgrid(GW(1).gvecB(i1,:),GW(2).gvecB(i2,:));
%             gw = (gw1A .* gw2A) + (gw1B .* gw2B);
%             % gaussian window storage for hilbert amplitude later:
%             gwloc = gw > gws;
%             gws(gwloc) = gw(gwloc);
%             % apply it
%             FM_voice = ifftn(FM .* gw);
%             if fullflag && ~guidedfourierflag % full 4D spectrum, if required
%                 ST.ST(:,:,i1,i2) = FM_voice;
%             end
%             % in any case, do the collapsed rapide spectrum
%             loc = abs(FM_voice) > abs(ST.C);
%             ST.C(loc) = FM_voice(loc);
%             ST.F1(loc) = ST.freqs{1}(i1);
%             ST.F2(loc) = ST.freqs{2}(i2);
%         end
%     case 3 % 3DST
%         for isc = 1:totalnumscales
%             % find indeces for scales for each dimension:
%             i1 = scale_inds(1,isc);
%             i2 = scale_inds(2,isc);
%             i3 = scale_inds(3,isc);
%             % assemble gaussian window
%             [gw1A,gw2A,gw3A] = ndgrid(GW(1).gvecA(i1,:),GW(2).gvecA(i2,:),GW(3).gvecA(i3,:));
%             [gw1B,gw2B,gw3B] = ndgrid(GW(1).gvecB(i1,:),GW(2).gvecB(i2,:),GW(3).gvecB(i3,:));
%             gw = (gw1A .* gw2A .* gw3A) + (gw1B .* gw2B .* gw3B);
%             % gaussian window storage for hilbert amplitude later:
%             gwloc = gw > gws;
%             gws(gwloc) = gw(gwloc);
%             % apply it
%             FM_voice = ifftn(FM .* gw);
%             if fullflag && ~guidedfourierflag % full 6D spectrum, if required
%                 ST.ST(:,:,:,i1,i2,i3) = FM_voice;
%             end
%             % in any case, do the collapsed rapide spectrum
%             loc = abs(FM_voice) > abs(ST.C);
%             ST.C(loc) = FM_voice(loc);
%             ST.F1(loc) = ST.freqs{1}(i1);
%             ST.F2(loc) = ST.freqs{2}(i2);
%             ST.F3(loc) = ST.freqs{3}(i3);
%         end
%         
%     case 4 % 4DST!! :O 
%         for isc = 1:totalnumscales
%             % find indeces for scales for each dimension:
%             i1 = scale_inds(1,isc);
%             i2 = scale_inds(2,isc);
%             i3 = scale_inds(3,isc);
%             i4 = scale_inds(4,isc);
%             % assemble gaussian window
%             [gw1A,gw2A,gw3A,gw4A] = ndgrid(GW(1).gvecA(i1,:),GW(2).gvecA(i2,:),GW(3).gvecA(i3,:),GW(4).gvecA(i4,:));
%             [gw1B,gw2B,gw3B,gw4B] = ndgrid(GW(1).gvecB(i1,:),GW(2).gvecB(i2,:),GW(3).gvecB(i3,:),GW(4).gvecB(i4,:));
%             gw = (gw1A .* gw2A .* gw3A .* gw4A) + (gw1B .* gw2B .* gw3B .* gw4B);
%             % gaussian window storage for hilbert amplitude later:
%             gwloc = gw > gws;
%             gws(gwloc) = gw(gwloc);
%             % apply it
%             FM_voice = ifftn(FM .* gw);
%             if fullflag && ~guidedfourierflag % full *D spectrum, if required
%                 ST.ST(:,:,:,:,i1,i2,i3,i4) = FM_voice;
%             end
%             % in any case, do the collapsed rapide spectrum
%             loc = abs(FM_voice) > abs(ST.C);
%             ST.C(loc) = FM_voice(loc);
%             ST.F1(loc) = ST.freqs{1}(i1);
%             ST.F2(loc) = ST.freqs{2}(i2);
%             ST.F3(loc) = ST.freqs{3}(i3);
%             ST.F4(loc) = ST.freqs{4}(i4);
%         end
% end


%% Flip things for 1D case:
if type == 1
ST.ST = ST.ST';
end


%% GET ABS AND REAL PARTS =================================================
% for the lazy (yes you corwin)
ST.A = abs(ST.C);
ST.R = real(ST.C);


%% (ST-FILTERED) HILBERT AMPLITUDE ========================================
% Replaces the old (well not that old) Hilbert Boosting method.

% take all those Gaussian windows from earlier and normalise to 1:
gws = gws ./ max(gws(:));
gws(gws < 0) = 0; % check for anomalies.
% this gives us a nice blob showing which parts of the fft spectrum we've
% considered in this ST, given the input scales.

% Keep the coefficients in zero freqs the same as the input data:
gws(imag(F) == 0) = 1;

% Now use this like a filter on your input data, and take the Hilbert
% transform of the result:
H = ifftn(FM .* gws);

% Assign the real() and abs() amplitudes:
ST.HA = abs(H);
ST.HR = real(H);

ST.allgws = gws;

% and you're done. Simple as that. I think back to all the years I've been
% thinking of how to do this, and here we are in a few lines. Mad.


%% HILBERT BOOSTING =======================================================
if boostflag
% %% HILBERT AMPLITUDE BOOSTING ALGORITHM:
%%%%% NEW!!!! %%%%%
% 2 - take hilbert transform (apply hilbert mask in fft space)
% 3 - take complex hilbert spectrum and ST spectrum, then boost the ST
% spectrum by the fraction of their abs() means, such that they have the
% same abs() mean.
% 4 - take the sqrt() of the product of the ST spectrum and the conj() of
% the hilbert spectrum to get covarying amplitude.

C_orig = ST.C;

% Get complex "Hilbert" spectrum of instantaneous amplitudes:
%     H_in = ifftn(ifftshift(fftshift(fftn(IN)) .* nph_hilbertmask(size(IN))));
F = fftn(IN);
H_in = ifftn(F .* nph_hilbertmask(F));
% note - not technically hilbert transform, that's only the complex part of
% the complex instantaneous phase object, where the real part is the
% original signal. Look up defintion of HIlbert transform.

% set both to have the same mean, so that the spectral energy is the same
% just redistributed in different frequencies:
hmean = mean(abs(H_in(:)),'omitnan');
stmean = mean(abs(ST.C(:)),'omitnan');
ST.C = ST.C .* (hmean/stmean);

% Covary the two to get covarying amplitude between them:
cv = sqrt(ST.C .* conj(H_in));
ST.C = ST.C .* (abs(cv) ./ abs(ST.C)); % boost by fractional difference at each location

ST.A = abs(ST.C);
ST.R = real(ST.C);

% NOTE this does NOT apply to the full ST object (the big 2-D/4-D/6-D
% one). The reason for this is that, of course, the instantaneous
% amplitude from the Hilbert transform is not defined for when the
% signal is decomposed for each frequency voice as in the Stockwell
% transform otherwise it would be, well, a Stockwell transform.

ST.BoostFactor = abs(ST.C) ./ abs(C_orig);

end



%==========================================================================

end  %  NDST FIN
%==========================================================================




%==========================================================================
% NESTED FUNCTIONS
%==========================================================================

% HILBERT MASK VERSION 2 ==================================================
% A newer, and much more simple, N-D approach to obtaining the analytic
% signal.
% We simply do what we say we do in the paper: find all the
% complex-conjugate pairs, set one to be zero and double the other. All
% coefficients not in a complex-conjugate pair are left unchanged.
% This approach uses linear indeces, so is N-D! Easy!

% Inputs: the fourier spectrum of the input data. Doesn't matter whether
% you've fftshift-ed it or not, the mask will be based on whatever
% arrangement you feed in.

function m = nph_hilbertmask(F)

% do fourier transform if input is real:
if isreal(F)
F = fftn(F);
end

% mask template of NaNs
m = nan(size(F));

% first, put in the zero freqs:
m(imag(F) == 0) = 1;

% now, find pairs:
[a,ib] = sort(abs(imag(F(:))));

% get rid of the zeros before reshaping:
ib = ib(a ~= 0);
a = a(a ~= 0);

% now reshape:
% this should always work - after the zeros are taken out there should
% always be an even number of complex conjugate pairs remaining.
% note: you need the transpose ' here due to the way reshape re-lists things.
ar = reshape(a',2,length(a)/2);
ibr = reshape(ib',2,length(ib)/2);

% now assign 2s and 0s (doesn't matter which order):
m(ibr(1,:)) = 2;
m(ibr(2,:)) = 0;

% and you're done!!

end



% % % % % % OLD HILBERT MASK (now disused, but useful background reading) ===========
% % % % % %
% % % % % % Matlab's hilb.m function only computes the transform across rows for N-D
% % % % % % matirices, need to write our own:
% % % % % %
% % % % % % m = nph_hilbertmask(sz,varargin);
% % % % % %
% % % % % % CHANGE LOG:
% % % % % %
% % % % % % 20171101 - new version of the Hilbert mask to recover "analytic" signal
% % % % % % from a 3D matrix, for use in the Stockwell Transform. NPH.
% % % % % %
% % % % % % 20180311 - Added an input 'neg' flag for the 3DST to select only negative
% % % % % % z freqs. We could have computed this based on inputting the individual
% % % % % % scale, but this would have meant computing the mask every timestep.
% % % % % %
% % % % % %
% % % % % %%% Background
% % % % % % In essence, all we are doing is creating a mask to double
% % % % % % selected frequencies and set their complex conjugate pairs to zero. We
% % % % % % leave and coefficients not in a complex conjugate pair alone.
% % % % % % On the surface, this is fairly simple. We set all +/-ve X,Y and only
% % % % % % +ve Z freqs to be doubled, and then all +/-ve X,Y and -ve Z freqs to be
% % % % % % zero. However, for a fftshifted 3D fourier spectrum, there are
% % % % % % 4 interlocking 2D sheets which need special attention. Where 3 planes
% % % % % % intersect at one location, they contain the "zeroth" frequency points,
% % % % % % ie ones not in a complex conjugate pair. For odd and even length
% % % % % % dimensions, these take funny shapes. This is just a quirk of the
% % % % % % arrangement of the frequencies in the FFT spectrum to result in an
% % % % % % output that is the same size as the input.
% % % % %
% % % % % % For an FFTSHIFTED spectrum, these are the 4 planes that are a little bit
% % % % % % special.
% % % % % %
% % % % % % For a spectrum with [EVEN EVEN EVEN] dimensions:
% % % % % %   .    .
% % % % % %   |\   |\
% % % % % %   | .----.----.           % Front, left side, base
% % % % % %   .-|\-.-|\-. |           % Mid-horz plane (always a feature)
% % % % % %   |\| .----.----.
% % % % % %   | x-|--x-|--. |
% % % % % %   .-|\|-.|\|  |\|
% % % % % %    \| x----x----.
% % % % % %     x-|--x-|--. |
% % % % % %      \|   \|   \|
% % % % % %       x----x----.
% % % % % %
% % % % % % For a spectrum with [ODD ODD ODD] dimensions:
% % % % % %   .    .
% % % % % %   |\   |\
% % % % % %   | .----.----.
% % % % % %   .-|\-.-|\-. |           % Mid-horz plane (always a feature)
% % % % % %   |\| .----.----.
% % % % % %   | .-|--x-|--. |
% % % % % %   .-|\|-.|\|  |\|
% % % % % %    \| .----.----.
% % % % % %     .-|--.-|--. |
% % % % % %      \|   \|   \|
% % % % % %       .----.----.
% % % % % %
% % % % % % "x" marks the possible zeroth frequencies (not in conj pair)
% % % % % %
% % % % % % Fortunately, each of these planes can follow the same pattern of masking
% % % % % % depending on its odd/even dimensions. Interestingly, each plane follows
% % % % % % the same 2D masking pattern as would be used for the 2DST. From what I
% % % % % % can gather, when the fft encounters (for ND > 1) even-numbered
% % % % % % dimensions, it dumps some complex conjugate pairs in a line at the edge
% % % % % % of the spectrum. For only odd dimensions, the each coeff sits happily with
% % % % % % its complex conjugate pair opposite it as a reflection through the very
% % % % % % centre of the spectrum, where the 0th freq component sits. If any
% % % % % % dimension is even, the fft dumps the extra coeffs in a line at the edge,
% % % % % % who sit opposite their complex conjugate pairs via a reflection with a
% % % % % % 0th freq component in the centre of that extra line (see below).
% % % % % %
% % % % % % The number of 2s should always equal the number of 0s.
% % % % % %
% % % % % % For a dimension with ODD number of elements, simply take the mask below
% % % % % % from inside the line, and for EVEN elements include the elements outside
% % % % % % of the lines. You essentially seem to trim the below to suit your needs
% % % % % % for each of these planes.
% % % % % %
% % % % % % MASK EXAMPLE:
% % % % % % EVEN n_elements case | ODD n_elements case
% % % % % %
% % % % % %         -ve X     +ve X
% % % % % %   2 | 0 0 0 0 2 2 2 2 2
% % % % % %   2 | 0 0 0 0 2 2 2 2 2 +ve Y
% % % % % %   2 | 0 0 0 0 2 2 2 2 2
% % % % % %   2 | 0 0 0 0 2 2 2 2 2
% % % % % %   1 | 0 0 0 0 1 2 2 2 2
% % % % % %   0 | 0 0 0 0 0 2 2 2 2
% % % % % %   0 | 0 0 0 0 0 2 2 2 2
% % % % % %   0 | 0 0 0 0 0 2 2 2 2 -ve Y
% % % % % %   0 | 0 0 0 0 0 2 2 2 2
% % % % % %   --+------------------
% % % % % %   1 | 0 0 0 0 1 2 2 2 2
% % % % % %
% % % % % %
% % % % % % Once you've set these 4 planes correctly, the rest is fairly
% % % % % % straightforward - simply set everything else above the middle horizontal
% % % % % % plane to be 2 (for +ve Z freqs) and set everything below (but not the
% % % % % % base) to be 0 (not -ve Zfreqs). You can of course fiddle this method to
% % % % % % consider different freq combinations, but tbh it's probably easier just
% % % % % % to permute you matrix to what you want and put it through this code :)
% % % % % %
% % % % %
% % % % % %%% Now the function itself!
% % % % % %
% % % % % % INPUTS: sz - size() of the desired mask, either 1, 2 or 3D.
% % % % % %
% % % % % % OUTPUTS: m - the mask itself, same size as input dimensions
% % % % % %
% % % % % % Note, we expect an fftshifted vector/matrix in each case.
% % % % % %
% % % % %
% % % % % function m = nph_hilbertmask(sz,varargin)
% % % % %
% % % % % % sz = size of the input matrix to be transformed.
% % % % % % negflag = 'neg' = some flag to denote that we need to create a negative
% % % % % % frequency-including mask on the 3rd dimension.
% % % % %
% % % % % % add scales, if supplied. This allows us to switch between +ve and -ve
% % % % % % masks for the 3rd dimension in the 3DST. Before, we only considered
% % % % % % positive z freqs.
% % % % % negflag = 0;
% % % % % if nargin == 2
% % % % %     if any(strcmpi(varargin{1},{'neg','negflag','-'}))
% % % % %         negflag = 1;
% % % % %     end
% % % % % end
% % % % %
% % % % %
% % % % %
% % % % % m = nan(sz); % ensure output is same size as specified
% % % % %
% % % % % sz(sz == 1) = []; % cope with 1x8, 8x1x5 etc.
% % % % %
% % % % % mid = fix(sz/2)+1; % find midpoint(s).
% % % % %
% % % % % switch length(sz)
% % % % %
% % % % %     %   1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     case 1 % 1D (straight from bob stockwell)
% % % % %
% % % % %         if isodd(sz)
% % % % %             m(1:mid-1)      = 0;
% % % % %             m(mid)          = 1;
% % % % %             m(mid+1:end)    = 2;
% % % % %         else
% % % % %             m(1)            = 1;
% % % % %             m(2:mid-1)      = 0;
% % % % %             m(mid)          = 1;
% % % % %             m(mid+1:end)    = 2;
% % % % %         end
% % % % %         %         m = [ones(1-rem(sz,2),1); 2*ones(fix((sz-1)/2),1); 1; zeros(fix((sz-1)/2),1)];
% % % % %
% % % % %         %   2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     case 2 % 2D (NPH and NDS method)
% % % % %         % make the 2D mask featured in the preamble:
% % % % %         m = mask2D(sz);
% % % % %
% % % % %         %   3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %     case 3 % 3D (NPH method)
% % % % %
% % % % %         % First, select all +/-ve X,Y and +ve Z freqs:
% % % % %         % SWITCH CASE FOR NEGATIVE Z FREQS:
% % % % %         switch negflag
% % % % %             case 0 % +ve z freqs
% % % % %                 m(:,:,mid(3)+1:end) = 2;
% % % % %                 m(:,:,1:mid(3)-1)   = 0;
% % % % %             case 1 % -ve z freqs
% % % % %                 m(:,:,mid(3)+1:end) = 0;
% % % % %                 m(:,:,1:mid(3)-1)   = 2;
% % % % %         end
% % % % %
% % % % %         % Do the middle slice (always done regardless of odd/even
% % % % %         % dimensions)
% % % % %         m(:,:,mid(3)) = mask2D(sz([1 2]));
% % % % %
% % % % %         % Now, determine what extra sides you need depending on the
% % % % %         % odd/even dimensions:
% % % % %
% % % % %         % Base
% % % % %         if iseven(sz(3))
% % % % %             m(:,:,1) = mask2D(sz([1 2]));
% % % % %         end
% % % % %
% % % % %         % left hand side
% % % % %         if iseven(sz(2))
% % % % %             m(:,1,:) = mask2D(sz([1 3]));
% % % % %         end
% % % % %
% % % % %         % front
% % % % %         if iseven(sz(1))
% % % % %             m(1,:,:) = mask2D(sz([2 3]));
% % % % %         end
% % % % %
% % % % %
% % % % %
% % % % % end
% % % % %
% % % % % % Sanity check:
% % % % % if numel(m == 2) ~= numel(m == 0)
% % % % %     disp('Warning: Problem with the Hilbert mask.')
% % % % %     return
% % % % % end
% % % % %
% % % % %
% % % % % end
% % % % %
% % % % %
% % % % % % 2D MASK ==========================================================================
% % % % % function m2 = mask2D(dims)
% % % % %
% % % % % mid = fix(dims/2)+1;
% % % % %
% % % % % % the order matters!
% % % % % m2 = nan(dims);
% % % % %
% % % % % m2(:,mid(2):end)            = 2; % double half
% % % % % m2(1:mid(1),1:mid(2))       = 0; % zero bottom quarter
% % % % % m2(mid(1):end,1:mid(2)-1)   = 0; % zero top quarter
% % % % %
% % % % % m2(mid(1),mid(2))           = 1; % middle 0th freq
% % % % %
% % % % % % If both sides are odd, the 2D mask is finished now! yayyy!
% % % % %
% % % % % % Bottom and left hand side edges?
% % % % % switch double(iseven(dims(1))) / double(iseven(dims(2)))
% % % % %     case Inf % [1 0] odd, even
% % % % %         m2(1,:) = [zeros(1,mid(2)-1) 1 2*ones(1,mid(2)-1)];
% % % % %     case 0   % [0 1] even, odd
% % % % %         m2(:,1) = [zeros(1,mid(1)-1) 1 2*ones(1,mid(1)-1)]';
% % % % %     case 1   % [1 1] even, even
% % % % %         m2(1,:) = [1 zeros(1,mid(2)-2) 1 2*ones(1,mid(2)-2)];
% % % % %         m2(:,1) = [1 zeros(1,mid(1)-2) 1 2*ones(1,mid(1)-2)]';
% % % % % end
% % % % %
% % % % % % Sanity check:
% % % % % if numel(m2 == 2) ~= numel(m2 == 0)
% % % % %     disp('Warning: Problem with the Hilbert mask.')
% % % % %     return
% % % % % end
% % % % %
% % % % % % flipud(m2 - nph_hilb2(dims))
% % % % % % figure; imagesc(1:dims(1),1:dims(2),m2);
% % % % %
% % % % %
% % % % % end
% % % % % %==========================================================================


%==========================================================================
% ISEVEN and ISODD functions
function logikal = iseven(x)
logikal = mod(x,2)==0;
end
function logikal = isodd(x)
logikal = mod(x,2)==1;
end
%==========================================================================

% % % % % %
% % % % % % % CREATE DEFAULTS =========================================================
% % % % % % function Defaults = create_defaults(IN,type)
% % % % % %
% % % % % % switch type
% % % % % %     case 1 % 1DST
% % % % % %         % Default scales = 1:Nyquist for 1DST;
% % % % % %         Defaults.scales = 1:length(IN)-1;
% % % % % %         Defaults.point_spacing = 1;
% % % % % %         Defaults.c = 0.25;
% % % % % %
% % % % % %     case 2 % 2DST
% % % % % %         % Default scales = 1:N/3
% % % % % %         Defaults.scales = {...
% % % % % %             [-fix(size(IN,1)/3):1:-1 1:fix(size(IN,1)/3)], ...
% % % % % %             [-fix(size(IN,2)/3):1:-1 1:fix(size(IN,2)/3)]};
% % % % % %         Defaults.point_spacing = [1 1];
% % % % % %         Defaults.c = [0.25 0.25];
% % % % % %
% % % % % %     case 3 % 3DST
% % % % % %         % Default scales = 1:N/3
% % % % % %         Defaults.scales = {...
% % % % % %             [-fix(size(IN,1)/3):1:-1 1:fix(size(IN,1)/3)], ...
% % % % % %             [-fix(size(IN,2)/3):1:-1 1:fix(size(IN,2)/3)], ...
% % % % % %             [-fix(size(IN,3)/3):1:-1 1:fix(size(IN,3)/3)]};
% % % % % %         Defaults.point_spacing = [1 1 1];
% % % % % %         Defaults.c = [0.25 0.25 0.25];
% % % % % %
% % % % % % end
% % % % % %
% % % % % % end
% % % % % % %==========================================================================




% FFT WISDOM =============================================================
function fftwisdom = gatherfftwisdom(IN,ffttype)
% So we may be able to significantly improve the speed of the ifftn
% computation below by using the fft wisdom library:

% method = 'exhaustive';
method = 'patient';

switch class(IN)
case 'double'
wisdomtype = 'dwisdom';
case 'single'
wisdomtype = 'swisdom';
end

fftw(wisdomtype,''); % clear any existing wisdom
fftw('planner',method); % choose method

switch lower(ffttype)
case {'fft' 'fftn'}
test_fftn = fftn(IN); % test fftn
case {'ifft' 'ifftn'}
test_ifftn = ifftn(IN); % test ifftn
end

fftwisdom = fftw(wisdomtype); % get wisdom

fftw(wisdomtype,fftwisdom); % apply the wisdom


end

































