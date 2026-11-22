function Data =  get_limbsounders(TimeRange,Instrument,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Generalised function to load and format limb sounder data.
%Loads a specific set of file formats - see readme.md in this repository.
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%inputs:
%   required:
%     TimeRange [double] - time range, in Matlab units. Details of how this is handled can be specified using 'TimeHandling' flag below.
%     Instrument [string] - instrument name to load, from specified list (see InstInfo struct below)
%
%
%   optional:
%
%     VarName         (type,       default)  description
%     ----------------------------------------------------------------------------- 
%     OriginalZ       (logical,      false)  return data on original vertical grid rather than interpolated to common scale.
%     KeepOutliers    (logical,      false)  don't remove outliers from the data. NOTE THAT BY DEFAULT THEY WILL BE REMOVED.
%     HeightScale     (numeric,  18:0.5:60)  heightscale to interpolate the data onto, in km, if OriginalZ is not set
%     HeightRange     (numeric,  [0,99e99])  height range to clip data to. Overrides and combines with HeightScale, but is more useful with OriginalZ.
%     LatRange        (numeric,   [-90,90])  latitude  range to select. Maximally permissive - allows profiles which enter the box at any height.
%     LonRange        (numeric, [-180,180])  longitude range to select. Also maximally permissive.
%     TimeHandling    (numeric,          3)  see list below
%     Verbose         (logical,       true)  update progress to screen
%
%       TimeHandling can be set to:
%         1. absolutely strictly - (e.g.) datenum(2010,1,[1,2])     will include all of 2010/01/01 and the first second of 2010/01/02
%         2. generously          - (e.g.) datenum(2010,1,[1.5,2.1]) will include all of 2010/01/01 and 2010/01/02, but not the first second of 2010/01/03 
%         3. fuzzy-ended         - (e.g.) datenum(2010,1,[1,2])     will include all of 2010/01/01 and 2010/01/02, but not the first second of 2010/01/03 
%         (3) is the default because it behaves almost as intuitively as (1) but reduces rutime by not loading a whole day of data to grab one second
%       A single TimeRange value t is treated as [t,t] by (1), i.e. only profiles at exactly that instant (a warning is given if none are found);
%       (2) and (3) return the whole day containing t.
%
%
%   advanced optional - only use these if you are confident you understand what will happen:
%
%     VarName         (type,       default)  description
%     -----------------------------------------------------------------------------
%     AdditionalVars  (cell,            {})  list of additional variables to extract, if available. Each is put on the [profiles x levels]
%                                            grid (per-profile and per-level variables are repeated; [levels x profiles] arrays are
%                                            transposed) and then processed like Temp. Variables whose values are all integers (flags,
%                                            counts) are interpolated with 'nearest' so they keep valid values; others use InterpType.
%                                            (An integer-valued physical variable would also be treated as a flag: use OriginalZ if so.)
%                                            Other shapes (e.g. scalars) stop with an error unless IgnoreNonModal is set.
%     DateWarning     (logical,       true)  warn the user that data aren't available for the requested date
%     GetHindleyPWs   (logical,      false)  get Hindley23 PW data for the instrument rather than raw data
%     IgnoreNonModal  (logical,      false)  pass through AdditionalVars that cannot be put on the [profiles x levels] grid (e.g. metadata)
%                                            without postprocessing (no interpolation, time/lat/lon subsetting or outlier removal)
%     FileSource      (logical,      false)  pass out original point locations as file list plus for each point a file and profile number
%     dx              (logical,      false)  return horizontal distance from last point at same level (be careful, in many cases this may not be meaningful)
%     MiscInfo        (cell,            {})  info to load miscellaneous data, containing {'path','identifier'). See notes below.
%     FilePath        (char,            '')  override path to load files from  
%     InterpType      (char,      'linear')  interpolation stategy to use if putting data on linear grid
%     ExtrapData      (logical,      false)  allow extrapolation of data beyond the limits of the input if interpolating
%     SamplePath      (char,            '')  path to sampled data, used when instrument is "Sampler"
%
%
%       MiscInfo can load files from unspecified sources. To use this option, Instrument should be set to 'Misc', and
%       the data files must be formatted the same way as the GNSS and ACE custom format. The MiscInfo array should contain:
%         MiscInfo{1}: the filesystem path containing the data. Data can be in subdirectories of this or at the top level.
%         MiscInfo{2}: an identifying string in the filename any place before the date. This can be the empty string ''.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%outputs:
%    Data: struct containing all variables, on a [profiles x height] grid
%            (exceptions to this size/shape can exist if some advanced-user optional flags are used)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%change log:
%  2023/09/13 added ACE-FTS as a valid instrument
%  2023/09/13 added ability to select profiles by lat/lon
%  2023/09/19 added MIPAS and SOFIE
%  2023/10/30 added filenames and profile numbers for backtracking to raw data 
%  2023/11/05 split off specific instrument cases into module files
%  2023/11/05 added option to load Hindley23 PW-filtered data rather than raw satellite data
%  2023/11/13 adjusted height-interpolation code to correctly interpolate longitudes near dateline
%  2023/11/25 added option to pass through unusually-shaped variables untouched
%  2024/12/29 added verbose option, and fixed bug in non-modal logic for OriginalZ calls
%  2025/02/27 added options to change interpolation strategy and allow extrapolation as an option
%  2025/04/18 modified loading strategy for pwdata to concatenate in chunks, speeding up large loads
%  2025/07/14 added reader for ALOMAR lidar
%  2026/08/17 added reader for data produced by Wright and Hindley (2018) sampler
%  2026/11/22 [Claude] fixed SABER times being offset by +1h01m
%  2026/11/22 [Claude] re-enabled trimming of output to TimeRange (profile level)
%  2026/11/22 [Claude] Instrument name canonicalised (e.g. 'saber' -> 'SABER'); DateWarning now warns rather than errors
%  2026/11/22 [Claude] outlier mask applied to all slices of 3-D variables (was first slice only)
%  2026/11/22 [Claude] HIRDLS AdditionalVars ScanUpFlag no longer crashes; progress bars no longer overshoot
%  2026/11/22 [Claude] column-vector HeightScale no longer crashes; reversed TimeRange now sorted
%  2026/11/22 [Claude] AdditionalVars put on the [profiles x levels] grid in modules; IgnoreNonModal compares with Lat, not modal size; exact name matching
%  2026/11/22 [Claude] lons < -180 wrapped; integer-valued AdditionalVars interpolated with 'nearest'; HIRDLS/MLS fill values -> NaN; ALOMAR file names
%  2026/11/22 [Claude] removed unused StrictTime option (use TimeHandling = 1)
%  2026/11/22 [Claude] longitudes unwrapped along each profile, not across profiles (fixes dateline-crossing profiles)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instruments we can use this on, and their properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ACE-FTS
%only currently loads temperature and pressure, as generated from GLC data
InstInfo.ACE.TimeRange      = [datenum(2004,1,35),datenum(9999,999,999)]; %still running at time of writing
InstInfo.ACE.HeightRange    = [18,125]; %data are available outside this range but are entirely a priori
InstInfo.ACE.Path           = [LocalDataDir,'/ACE/raw/'];

%AIRS-3D
%this is NOT RECOMMENDED for use as AIRS is a nadir sounder, but it will work if absolutely needed
%this setting will load it as if it was a limb sounder, for below-specified granules and level of data thinning
InstInfo.AIRS.TimeRange      = [datenum(2002,8,1),datenum(9999,999,999)]; %still running at time of writing
InstInfo.AIRS.HeightRange    = [0,90]; %a lot of this has awful resolution, but is fine for bulk T. 
InstInfo.AIRS.Path           = [LocalDataDir,'/AIRS/3d_airs/'];
InstInfo.AIRS.ThinFactor     = [1,1]; %only take profiles this often in [xt, at] directions, to reduce data volume
InstInfo.AIRS.Granules       = 1:1:240;

%IASI-3D
%this is also NOT RECOMMENDED for use 
InstInfo.IASI.TimeRange      = [datenum(2007,1,150),datenum(9999,999,999)];  %still running at time of writing
InstInfo.IASI.HeightRange    = [0,90]; %a lot of this has awful resolution, but is fine for bulk T. 
InstInfo.IASI.Path           = [LocalDataDir,'/IASI/'];
InstInfo.IASI.ThinFactor     = [1,1]; %only take profiles this often in [xt, at] directions, to reduce data volume
InstInfo.IASI.Granules       = 1:1:240;

%ALOMAR Lidar
InstInfo.ALOMAR.TimeRange   = [datenum(2003,1,1),datenum(9999,999,999)];
InstInfo.ALOMAR.HeightRange = [0,80];
InstInfo.ALOMAR.Path        = [LocalDataDir,'/lidar/ALOMAR/'];

%GNSS
InstInfo.GNSS.TimeRange     = [datenum(2000,1,1),datenum(9999,999,999)]; %still running at time of writing
InstInfo.GNSS.HeightRange   = [0,40];
InstInfo.GNSS.Path          = [LocalDataDir,'/GNSS/raw/'];

%HIRDLS
InstInfo.HIRDLS.TimeRange   = [datenum(2005,1,29),datenum(2008,1,77)];
InstInfo.HIRDLS.HeightRange = [0,80];
InstInfo.HIRDLS.Path        = [LocalDataDir,'/HIRDLS/raw/'];

%MIPAS
InstInfo.MIPAS.TimeRange   = [datenum(2002,7,2),datenum(2012,4,8)];
InstInfo.MIPAS.HeightRange = [0,80];
InstInfo.MIPAS.Path        = [LocalDataDir,'/MIPAS/raw/'];

%MLS
%only currently loads temperature files, due to the way I have them stored
InstInfo.MLS.TimeRange      = [datenum(2004,1,275),datenum(9999,999,999)]; %still running at time of writing
InstInfo.MLS.HeightRange    = [0,100];
InstInfo.MLS.Path           = [LocalDataDir,'/MLS/raw/'];

%SABER
InstInfo.SABER.TimeRange   = [datenum(2002,1,1),datenum(9999,999,999)]; %still running at time of writing
InstInfo.SABER.HeightRange = [0,120];
InstInfo.SABER.Path        = [LocalDataDir,'/SABER/raw/'];

%SOFIE
InstInfo.SOFIE.TimeRange   = [datenum(2007,1,135),datenum(9999,999,999)]; %still running at time of writing
InstInfo.SOFIE.HeightRange = [10,110]; %this is the range of the mission time/height cross-sections on their website
InstInfo.SOFIE.Path        = [LocalDataDir,'/SOFIE/raw/'];

%output from the Wright and Hindley (2018) sampler - only works for 2D data
InstInfo.Sampler.TimeRange   = [-1,1].*datenum(9999,999,999); %any date
InstInfo.Sampler.HeightRange = [-1,1].*9999;                  %any height
InstInfo.Sampler.Path        = '';                            %set below

%Miscellaneous data
InstInfo.Misc.TimeRange   = [-1,1].*datenum(9999,999,999); %any date
InstInfo.Misc.HeightRange = [-1,1].*9999;                  %any height
InstInfo.Misc.Path        = '';                            %set below

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

%date
CheckDates  = @(x) validateattributes(x,{'numeric'},{'>=',datenum(1979,1,1)}); %must be in the satellite era
addRequired(p,'TimeRange',CheckDates); %time range

%instrument(s)
CheckInst  = @(x) validateStringParameter(x,fieldnames(InstInfo),mfilename,'Instrument');
addRequired(p,'Instrument',CheckInst)

%height range and height scale
CheckHeights = @(x) validateattributes(x,{'numeric'},{'>=',0}); 
addParameter(p,'HeightScale',18:0.5:60,CheckHeights)
addParameter(p,'HeightRange',[0,99e99],CheckHeights)

%latrange and lonrange
addParameter(p,'LatRange',[ -90, 90],@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90}))
addParameter(p,'LonRange',[-180,180],@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180}))

%additional variables
addParameter(p, 'AdditionalVars',      {},@iscell   )
addParameter(p,      'OriginalZ',   false,@islogical)
addParameter(p,   'KeepOutliers',   false,@islogical)
addParameter(p,     'FileSource',   false,@islogical)
addParameter(p,   'TimeHandling',       3,@isnumeric)
addParameter(p,  'GetHindleyPWs',   false,@islogical)
addParameter(p,    'DateWarning',    true,@islogical)
addParameter(p, 'IgnoreNonModal',   false,@islogical)
addParameter(p,        'Verbose',    true,@islogical)
addParameter(p,       'MiscInfo',      {},@iscell   )
addParameter(p,             'dx',   false,@islogical)
addParameter(p,       'FilePath',      '',@ischar   )
addParameter(p,     'InterpType', 'linear',@ischar  )
addParameter(p,     'ExtrapData',   false,@islogical)
addParameter(p,     'SamplePath',      '',@ischar   )




%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,TimeRange,Instrument,varargin{:})

%pull out the contents into struct "Settings", used throughout rest of routine
Settings = p.Results;
Settings.TimeRange = TimeRange;
clearvars -except InstInfo Settings

%replace the instrument name with its canonical form (validatestring accepts any case and unambiguous
%abbreviations, e.g. 'saber' or 'SAB', but struct fields and the module switch below need e.g. 'SABER')
Settings.Instrument = validatestring(Settings.Instrument,fieldnames(InstInfo));

%extract just the metadata for the instrument we want
InstInfo  = InstInfo.(Settings.Instrument);

%force the height scale to a row vector: the interpolation and the Alt output below assume one
%(a column input previously crashed in interp_1d_ndims)
Settings.HeightScale = Settings.HeightScale(:)';

%clip the modified height scale to the specific height range, in case we specified both
idx = inrange(Settings.HeightScale,Settings.HeightRange);
Settings.HeightScale = Settings.HeightScale(idx);
clear idx



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time range handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%additional check on dates - must be either a single value (a day to load all of) or two values (start and end time)
if numel(Settings.TimeRange)  > 2;  
  error('TimeRange must be either one value (a day of interest) or two values (a time range)');
end

%additional check on dates - warn if entirely outside valid range for instrument (loading continues, and will return no profiles)
if Settings.DateWarning == true                          ...
   &&   (min(Settings.TimeRange) > InstInfo.TimeRange(2) ...
   ||   max(Settings.TimeRange) < InstInfo.TimeRange(1))
  warning(['Data for ',Settings.Instrument,' is only available from ',datestr(InstInfo.TimeRange(1)),' to ',datestr(InstInfo.TimeRange(2))]);
end

%duplicate times if single value given
ScalarTime = isscalar(Settings.TimeRange); %remembered for TimeHandling 1 and 3 below
if isscalar(Settings.TimeRange); Settings.TimeRange = [1,1].*Settings.TimeRange; end

%put times in ascending order: TimeHandling below assumes [start,end] (a reversed range previously returned no profiles)
Settings.TimeRange = sort(Settings.TimeRange);

%handle time strictness
switch Settings.TimeHandling
  case 1; %do nothing, output should be literally what we asked for
  case 2;
    %expand out to include a most generous range of days the user could have meant
    Settings.TimeRange(1) = floor(Settings.TimeRange(1));
    if mod(Settings.TimeRange(2),1) == 0; Settings.TimeRange(2) = Settings.TimeRange(2)+1-1e-8;
    else;                                 Settings.TimeRange(2) = ceil(Settings.TimeRange(2))-1e-8;
    end
  case 3;
    %a single non-integer value means the whole day containing it
    %otherwise, if the last entry is an integer, feather it slightly to avoid loading an extra day
    if ScalarTime && mod(Settings.TimeRange(2),1) ~= 0; Settings.TimeRange = [floor(Settings.TimeRange(1)),floor(Settings.TimeRange(1))+1-1e-8];
    elseif mod(Settings.TimeRange(2),1) == 0;          Settings.TimeRange(2) = Settings.TimeRange(2)+1-1e-8;
    end
  otherwise
    error('Invalid time handling option chosen')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this function can also load Hindley23 PW data for each instrument
%this is assumed to be in the same root path as the data, but with
%/raw replaced with /pws. We call this by overriding the InstInfo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.GetHindleyPWs == true

  %first, alter the path to call the PW directory
  InstInfo.Path = strrep(InstInfo.Path,'/raw/','/pws/');

  %now, produce a copy of the instrument's real name
  InstInfo.Inst = Settings.Instrument;

  %and replace the instrument name with a call to the PW loader
  Settings.Instrument = 'HindleyPWs';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% when using miscellaneous data or an override path, set the path here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%misc
if strcmpi(Settings.Instrument,'Misc');
  if numel(Settings.MiscInfo) == 0; error('Miscellaneous instrument selected but MiscInfo cell not set - see instructions in header'); end
  InstInfo.Path       = Settings.MiscInfo{1};
  InstInfo.Identifier = Settings.MiscInfo{2};
end


%override
if numel(Settings.FilePath) > 0;InstInfo.Path = Settings.FilePath; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading - instrument specific, see modules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Verbose == 1; disp('--> Loading data'); end

%at the end of this we want a struct called Data containing the following
%variables on a grid of [profiles x levels]
%
%Lat  - latitude
%Lon  - longitude, using -180 to 180 range
%Time - time, in Matlab units
%Temp - temperature, in K
%Pres - pressure, in hPa
%Alt  - height, in km
% + any additional vars requested if valid for that instrument

%list of variables
if Settings.GetHindleyPWs == true; Vars = [{'Lat','Lon','Time','Temp_Residual','Temp_PW','Pres','Alt','SourceProf','SourceFile'},Settings.AdditionalVars];
else;                              Vars = [{'Lat','Lon','Time','Temp',                   'Pres','Alt','SourceProf','SourceFile'},Settings.AdditionalVars];
end

%get the data for this instrument using the appropriate module
switch Settings.Instrument
 %standard instruments
  case 'ALOMAR';              [Data,FileList] = module_load_ALOMAR(  Settings,InstInfo,Vars);
  case {'ACE','GNSS','Misc'}; [Data,FileList] = module_load_ACE_GNSS(Settings,InstInfo,Vars);
  case {'AIRS','IASI'};       [Data,FileList] = module_load_AIRS(    Settings,InstInfo,Vars);
  case 'HIRDLS';              [Data,FileList] = module_load_HIRDLS(  Settings,InstInfo,Vars);
  case 'MIPAS';               [Data,FileList] = module_load_MIPAS(   Settings,InstInfo,Vars);
  case 'MLS';                 [Data,FileList] = module_load_MLS(     Settings,InstInfo,Vars);
  case 'SABER';               [Data,FileList] = module_load_SABER(   Settings,InstInfo,Vars);
  case 'SOFIE';               [Data,FileList] = module_load_SOFIE(   Settings,InstInfo,Vars);
 %sampler output
  case 'Sampler';             [Data,FileList] = module_load_sampler( Settings,InstInfo,Vars); 
 %pw data loader  
  case 'HindleyPWs';          [Data,FileList] = module_load_pwdata(  Settings,InstInfo,Vars);
 %fail case
  otherwise
    disp(['Instrument ',Settings.Instrument,' not currently handled by this function, terminating'])
    return
end
clear FileCount;

%in the below postprocessing, variables which are not of the 'standard' shape (profs x height) will
%cause issues. The modules put AdditionalVars onto this grid where they can (module_conform_vars), and
%stop with an error otherwise, unless IgnoreNonModal is set. In that case, list the variables whose
%first two dimensions do not match Lat (the reference [profiles x levels] grid) so they are passed
%through untouched. Tables are not listed: they are subset by profile like other data below.
%(Previously the modal size was guessed with list_non_modal_size, which failed if any field was 3-D.)
VarsToIgnore = {};
if Settings.IgnoreNonModal == true;
  f = fieldnames(Data);
  for iF=1:1:numel(f)
    if istable(Data.(f{iF})); continue; end
    if ~isequal(size(Data.(f{iF}),[1,2]),size(Data.Lat)); VarsToIgnore{end+1} = f{iF}; end
  end; clear iF f
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data to chosen height scale
%  and make lons -180 to 180
%  select actual time range specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time range: do first to reduce size for later steps
%do at profile level to avoid breaking profiles
Data = reduce_struct(Data,inrange(mean(Data.Time,2,'omitnan'),Settings.TimeRange),VarsToIgnore,1);

%strict handling of a single value selects only profiles at exactly that instant: warn if that found nothing
if ScalarTime && Settings.TimeHandling == 1 && size(Data.Time,1) == 0
  warning('get_limbsounders:StrictSingleTime',['A single TimeRange value with TimeHandling=1 selects only profiles at exactly that instant, ' ...
          'and none were found. Use TimeHandling 2 or 3 (default) to get the whole day.']);
end

if Settings.OriginalZ == false

  %unwrap longitudes along each profile (dim 2, levels), to make interpolation safe across the dateline
  %(fixed 2026/11/22: unwrap() of a matrix works down columns, i.e. across profiles, so a profile crossing the
  % dateline was interpolated between e.g. 179 and -179 deg. Only finite values are unwrapped, so NaN levels cannot break it)
  Data.Lon = deg2rad(Data.Lon);
  for iProf=1:1:size(Data.Lon,1)
    Fin = isfinite(Data.Lon(iProf,:));
    Data.Lon(iProf,Fin) = unwrap(Data.Lon(iProf,Fin));
  end; clear iProf Fin

  %do the height interpolation over all variables
  Data2 = struct();
  if Settings.Verbose == 1; textprogressbar('----> Interpolating to output z '); end
  for iVar=1:1:numel(Vars);

    if Settings.Verbose == 1; textprogressbar(iVar./numel(Vars).*100); end

    %skip if this variable is on the ignore list
    if any(strcmp(VarsToIgnore,Vars{iVar})); Data2.(Vars{iVar}) = Data.(Vars{iVar}); continue;  end %(exact name match; was contains(), a substring match)

    %copy over altitude - not computationally necessary, but produces some tiny numerical precision issues
    % in a very small number of profiles (about a 100-millionth of the data) later if we don't. This was a
    %hard bug to find!
    if strcmpi(Vars{iVar},'Alt'); Data2.Alt = repmat(Settings.HeightScale,size(Data2.Lat,1),1); continue; end

    %skip if it's a table (this is the case for GNSS metadata)
    if strcmp(class(Data.(Vars{iVar})),'table'); Data2.(Vars{iVar}) = Data.(Vars{iVar}); continue; end

    b = Data.(Vars{iVar});
    sz = size(b); sz(2) = numel(Settings.HeightScale);
    a = NaN(sz);
    if ismatrix(a); sz(3) = 1; end

    %AdditionalVars whose finite values are all integers (flags, counts, indices) are interpolated with
    %'nearest', so the output contains only values that occur in the data; everything else uses InterpType.
    %(Per-profile integers, e.g. orbit numbers, are constant down each profile, so are unaffected.)
    ThisInterp = Settings.InterpType;
    if any(strcmp(Settings.AdditionalVars,Vars{iVar})) && (isnumeric(b) || islogical(b))
      v = b(isfinite(b));
      if ~isempty(v) && all(v == round(v)); ThisInterp = 'nearest'; end
      clear v
    end
    for iProf=1:1:sz(1)

      %some extra handling here to deal with bad data and higher-dimension data, but all
      %we're actually doing is linear interpolation
      Good = find(isfinite(Data.Alt(iProf,:)) ~=0);
      [~,uidx] = unique(Data.Alt(iProf,:));
      Good = intersect(Good,uidx);
      if numel(Good) > 2;
        for iExtraDim=1:1:sz(3);
          if Settings.ExtrapData == 1;
            a(iProf,:,iExtraDim) = interp_1d_ndims(Data.Alt(iProf,Good),b(iProf,Good,iExtraDim),Settings.HeightScale,2,ThisInterp,'extrap');
          else
            a(iProf,:,iExtraDim) = interp_1d_ndims(Data.Alt(iProf,Good),b(iProf,Good,iExtraDim),Settings.HeightScale,2,ThisInterp);
          end
        end
      end

    end;
    Data2.(Vars{iVar}) = a;
  end
  if Settings.Verbose == 1; textprogressbar('!'); end  
  Data = Data2;

  %%wrap the longitudes back into the normal range
  Data.Lon = rad2deg(wrapToPi(Data.Lon));

  clear iProf iVar Data2 a b sz Good uidx
else

  %clip to height range: points outside the range are set to NaN in every [profiles x levels(x ...)] field,
  %not just Alt, then levels (dim 2) with no in-range height left are removed
  BadMask = Data.Alt < min(Settings.HeightRange) | Data.Alt > max(Settings.HeightRange);
  RowsWithData = find(any(isfinite(Data.Alt) & ~BadMask,1)); %(these are levels, despite the name)

  f = fieldnames(Data);
  for iF=1:1:numel(f)
    if any(strcmp(VarsToIgnore,f{iF})); continue;  end
    if strcmp(class(Data.(f{iF})),'table');  continue; end
    g = Data.(f{iF});
    if isequal(size(g,[1,2]),size(BadMask)); g(repmat(BadMask,[1,1,size(g,3:ndims(g))])) = NaN; end %same mask on every trailing slice
    [g,DimSize,DimOrder] = expose_dim(g,2);
    g = g(RowsWithData,:);
    Data.(f{iF}) = permute(reshape(g,[size(g,1),DimSize(2:end)]),DimOrder);
    clear DimSize DimOrder
  end; clear iF f g RowsWithData BadMask

end

%longitude
Data.Lon(Data.Lon > 180) = Data.Lon(Data.Lon > 180)-360;
Data.Lon(Data.Lon < -180) = Data.Lon(Data.Lon < -180)+360; %(added 2026/11/22: these were previously removed as outliers)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postprocessing based on input options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Verbose == 1; disp('--> Postprocessing'); end

%remove outliers?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.KeepOutliers == 0;

  %create one big array of badness, then remove from ALL variables in one pass
  Bad = [];

  %latitude and longitude have physical limits
  Bad = [Bad;find(Data.Lat <  -90 | Data.Lat >  90 | Data.Lon < -180 | Data.Lon > 180)];
  
  %altitude should always be in the specified range
  if Settings.OriginalZ == false;  Bad = [Bad;find(Data.Alt  < min(Settings.HeightScale) | Data.Alt  > max(Settings.HeightScale))]; end

  % %time should always be in the specified range
  % Bad = [Bad;find(Data.Time < min(Settings.TimeRange  ) | Data.Time > max(Settings.TimeRange  ))];  

  %temperature should be >100K always, and  <400K at altitudes below the mesopause 
  %it may not exist if we're loading PW data
  if isfield(Data,'Temp');
    Bad = [Bad;find(Data.Temp < 100)];
    Bad = [Bad;find(Data.Temp > 400 & Data.Alt < 100)];
  end

  %do it
  %Bad holds linear indices into the [profiles x levels] grid. Convert to a logical mask so that it can be
  %applied to every slice of variables with extra trailing dimensions (e.g. Temp_PW, [profs x levels x modes]);
  %applying the linear indices directly would only mask the first slice
  Bad = unique(Bad);
  BadMask = false(size(Data.Lat)); BadMask(Bad) = true;
  Fields = fieldnames(Data);
  for iF=1:1:numel(Fields)

    %skip if this variable is on the ignore list
    if any(strcmp(VarsToIgnore,Fields{iF}));continue;  end
    %skip if it's a table (this is the case for GNSS metadata)
    if strcmp(class(Data.(Fields{iF})),'table'); continue;  end


    F = Data.(Fields{iF});
    if ~ismatrix(F) && isequal(size(F,[1,2]),size(BadMask))
      F(repmat(BadMask,[1,1,size(F,3:ndims(F))])) = NaN; %same mask on every trailing slice
    else  F(Bad) = NaN;   
    end
    Data.(Fields{iF}) = F;
  end

  clear Bad BadMask Fields F iF
  if Settings.Verbose == 1; disp('----> Outliers removed'); end
end


%latitude and longitude filtering?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(Settings.LatRange) ~=  -90 || max(Settings.LatRange) ~=  90 ...
|| min(Settings.LonRange) ~= -180 || max(Settings.LonRange) ~= 180

  %find the min and max of lon and lat for each profile
  Limits = [nanmin(Data.Lon,[],2),nanmax(Data.Lon,[],2), ...
            nanmin(Data.Lat,[],2),nanmax(Data.Lat,[],2)];

  
  %generate a list of profiles where any point falls inside the limits
  Bad = [];
  Bad = [Bad;find(Limits(:,1) > max(Settings.LonRange))];
  Bad = [Bad;find(Limits(:,2) < min(Settings.LonRange))];
  Bad = [Bad;find(Limits(:,3) > max(Settings.LatRange))];
  Bad = [Bad;find(Limits(:,4) < min(Settings.LatRange))];

  Good = 1:1:size(Limits,1); Good(Bad) = [];
  Data = reduce_struct(Data,Good,VarsToIgnore,1);

  if Settings.Verbose == 1; disp('----> Geographically filtered'); end

end


%do we want to know where the data came from?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.FileSource == 1;
  %add the filelist to the structure
  Data.OriginalFiles = FileList;
else
  %remove tracker variables
  Data = rmfield(Data,{'SourceFile','SourceProf'});
end

clear FileList;



%do we want the distance between adjacent points?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.dx == 1
  Data.dx = NaN.*Data.Lat;
  for iProf=2:1:size(Data.dx,1)
    ProfA = [Data.Lat(iProf-1,:);Data.Lon(iProf-1,:)];
    ProfB = [Data.Lat(iProf,  :);Data.Lon(iProf,  :)];
    Data.dx(iProf,:) = nph_haversine(ProfA',ProfB');
  end
  if Settings.Verbose == 1; disp('----> Pair distances computed'); end

end


%end of programme
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to validate list of allowed instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validateStringParameter(varargin)
validatestring(varargin{:});
end



