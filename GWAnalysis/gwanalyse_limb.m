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




