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
%     HeightRange     (numeric,  [0,99e99])  height range to clip data to. Combines with HeightScale, but is more useful with OriginalZ.
%     LatRange        (numeric,   [-90,90])  latitude  range to select. Maximally permissive - allows profiles which enter the box at any height.
%     LonRange        (numeric, [-180,180])  longitude range to select. Also maximally permissive.
%     TimeHandling    (numeric,          3)  see list below
%
%       TimeHandling can be set to:
%         1. absolutely strictly - (e.g.) datenum(2010,1,[1,2])     will include all of 2010/01/01 and the first second of 2010/01/02
%         2. generously          - (e.g.) datenum(2010,1,[1.5,2.1]) will include all of 2010/01/01 and 2010/01/02, but not the first second of 2010/01/03 
%         3. fuzzy-ended         - (e.g.) datenum(2010,1,[1,2])     will include all of 2010/01/01 and 2010/01/02, but not the first second of 2010/01/03 
%         (3) is the default because it behaves almost as intuitively as (1) but reduces rutime by not loading a whole day of data to grab one second
%
%
%   advanced optional - only use these if you are confident you understand what will happen:
%
%     VarName         (type,       default)  description
%     -----------------------------------------------------------------------------
%     AdditionalVars  (cell,            {})  list of additional variables to extract, if available 
%     DateWarning     (logical,       true)  warn the user that data aren't available for the requested date
%     GetHindleyPWs   (logical,      false)  get Hindley23 PW data for the instrument rather than raw data
%     IgnoreNonModal  (logical,      false)  pass through variables which are not the MODAL size and shape (e.g. metadata) without postprocessing
%     FileSource      (logical,      false)  pass out original point locations as file list plus for each point a file and profile number
%     MiscInfo        (cell,            {})  info to load miscellaneous data, containing {'path','identifier'). See notes below.
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
InstInfo.AIRS.ThinFactor     = [5,5]; %only take profiles this often in [xt, at] directions, to reduce data volume
InstInfo.AIRS.Granules       = 1:1:240;

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
CheckInst  = @(x) validateStringParameter(x,fieldnames(InstInfo),mfilename,Instrument);
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
addParameter(p,     'StrictTime',   false,@islogical)
addParameter(p,   'TimeHandling',       3,@isnumeric)
addParameter(p,  'GetHindleyPWs',   false,@islogical)
addParameter(p,    'DateWarning',    true,@islogical)
addParameter(p,  'IgnoreNonModal',  false,@islogical)
addParameter(p,        'MiscInfo',      {},@iscell   )




%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,TimeRange,Instrument,varargin{:})

%pull out the contents into struct "Settings", used throughout rest of routine
Settings = p.Results;
Settings.TimeRange = TimeRange;
clearvars -except InstInfo Settings

%extract just the metadata for the instrument we want
InstInfo  = InstInfo.(Settings.Instrument);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time range handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%additional check on dates - must be either a single value (a day to load all of) or two values (start and end time)
if numel(Settings.TimeRange)  > 2;  
  error('TimeRange must be either one value (a day of interest) or two values (a time range)');
end

%additional check on dates - must be in valid range for instrument
if Settings.DateWarning == true                         ...
   &&   min(Settings.TimeRange) > InstInfo.TimeRange(2) ...
   |    max(Settings.TimeRange) < InstInfo.TimeRange(1)
  error(['Data for ',Settings.Instrument,' is only available from ',datestr(InstInfo.TimeRange(1)),' to ',datestr(InstInfo.TimeRange(2))]);
end

%duplicate times if single value given
if numel(Settings.TimeRange) == 1; Settings.TimeRange = [1,1].*Settings.TimeRange; end

%handle time strictness
switch Settings.TimeHandling
  case 1; %do nothing, output should be literally what we asked for
  case 2;
    %expand out to include a most generous range of days the user could have meant
    Settings.TimeRange(1) = floor(Settings.TimeRange(1));
    if mod(Settings.TimeRange(2),1) == 0; Settings.TimeRange(2) = Settings.TimeRange(2)+1-1e-8;
    else                                  Settings.TimeRange(2) = ceil(Settings.TimeRange(2))-1e-8;
    end
  case 3;
    %if the last entry is an integer, feather it slightly to avoid loading an extra day
    if mod(Settings.TimeRange(2),1) == 0; Settings.TimeRange(2) = Settings.TimeRange(2)+1-1e-8; end
  otherwise
    error(['Invalid time handling option chosen'])
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
%% when using miscellaneous data, set the path here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(Settings.Instrument,'Misc');
  if numel(Settings.MiscInfo) == 0; error('Miscellaneous instrument selected but MiscInfo cell not set - see instructions in header'); end
  InstInfo.Path       = Settings.MiscInfo{1};
  InstInfo.Identifier = Settings.MiscInfo{2};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading - instrument specific, see modules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
else                               Vars = [{'Lat','Lon','Time','Temp',                   'Pres','Alt','SourceProf','SourceFile'},Settings.AdditionalVars];
end

%get the data for this instrument using the appropriate module
switch Settings.Instrument
 %standard instruments
  case {'ACE','GNSS','Misc'}; [Data,FileList] = module_load_ACE_GNSS(Settings,InstInfo,Vars);
  case 'AIRS';                [Data,FileList] = module_load_AIRS(    Settings,InstInfo,Vars);
  case 'HIRDLS';              [Data,FileList] = module_load_HIRDLS(  Settings,InstInfo,Vars);
  case 'MIPAS';               [Data,FileList] = module_load_MIPAS(   Settings,InstInfo,Vars);
  case 'MLS';                 [Data,FileList] = module_load_MLS(     Settings,InstInfo,Vars);
  case 'SABER';               [Data,FileList] = module_load_SABER(   Settings,InstInfo,Vars);
  case 'SOFIE';               [Data,FileList] = module_load_SOFIE(   Settings,InstInfo,Vars);
 %pw data loader  
  case 'HindleyPWs';          [Data,FileList] = module_load_pwdata(  Settings,InstInfo,Vars);
 %fail case
  otherwise
    disp(['Instrument ',Settings.Instrument,' not currently handled by this function, terminating'])
    return
end
clear FileCount;

%in the below postprocessing, variables which are not of the 'standard' shape (profs x height) will
%cause issues. So, let's create a list of variables to be ignored

%get a list of non-modal-size variables to ignore
if Settings.IgnoreNonModal == true;
  VarsToIgnore = list_non_modal_size(Data);
else
  VarsToIgnore = {};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data to chosen height scale
%  and make lons -180 to 180
%  select actual time range specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time range: do first to reduce size for later steps
%do at profile level to avoid breaking profiles
Data = reduce_struct(Data,inrange(nanmean(Data.Time,2),Settings.TimeRange),VarsToIgnore,1);


if Settings.OriginalZ == false

  %unwrap longitudes, to make interpolation safe
  Data.Lon = unwrap(deg2rad(Data.Lon));

  %do the height interpolation over all variables
  Data2 = struct();
  for iVar=1:1:numel(Vars);

    %skip if this variable is on the ignore list
    if find(contains(VarsToIgnore,Vars{iVar})); Data2.(Vars{iVar}) = Data.(Vars{iVar}); continue;  end

    b = Data.(Vars{iVar});
    sz = size(b); sz(2) = numel(Settings.HeightScale);
    a = NaN(sz);
    if ndims(a) == 2; sz(3) = 1; end
    for iProf=1:1:sz(1)

      %some extra handling here to deal with bad data and higher-dimension data, but all
      %we're actually doing is linear interpolation
      Good = find(isfinite(Data.Alt(iProf,:)) ~=0);
      [~,uidx] = unique(Data.Alt(iProf,:));
      Good = intersect(Good,uidx);
      if numel(Good) > 2;
        for iExtraDim=1:1:sz(3);
          a(iProf,:,iExtraDim) = interp_1d_ndims(Data.Alt(iProf,Good),b(iProf,Good,iExtraDim),Settings.HeightScale,2);
        end
      end

    end;
    Data2.(Vars{iVar}) = a;
  end
  Data = Data2;

  %%wrap the longitudes back into the normal range
  Data.Lon = rad2deg(wrapToPi(Data.Lon));

  clear iProf iVar Data2 a b sz Good uidx
end

%longitude
Data.Lon(Data.Lon > 180) = Data.Lon(Data.Lon > 180)-360;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postprocessing based on input options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove outliers?
%this section does the time filtering as well - if it isn't run you'll get whole days only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.KeepOutliers == 0;

  %create one big array of badness, then remove from ALL variables in one pass
  Bad = [];

  %latitude and longitude have physical limits
  Bad = [Bad;find(Data.Lat <  -90 | Data.Lat >  90 | Data.Lon < -180 | Data.Lon > 180)];
  
  %altitude should always be in the specified range
  if Settings.OriginalZ == false;  Bad = [Bad;find(Data.Alt  < min(Settings.HeightScale) | Data.Alt  > max(Settings.HeightScale))]; end

  %time should always be in the specified range
  Bad = [Bad;find(Data.Time < min(Settings.TimeRange  ) | Data.Time > max(Settings.TimeRange  ))];  

  %temperature should be >100K always, and  <400K at altitudes below the mesopause 
  %it may not exist if we're loading PW data
  if isfield(Data,'Temp');
    Bad = [Bad;find(Data.Temp < 100)];
    Bad = [Bad;find(Data.Temp > 400 & Data.Alt < 100)];
  end

  %do it
  Bad = unique(Bad);
  Fields = fieldnames(Data);
  for iF=1:1:numel(Fields)

    %skip if this variable is on the ignore list
    if find(contains(VarsToIgnore,Fields{iF}));continue;  end

    F = Data.(Fields{iF});
    F(Bad) = NaN;
    Data.(Fields{iF}) = F;
  end

  clear Bad Fields F iF
end


%latitude and longitude filtering?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(Settings.LatRange) ~=  -90 || max(Settings.LatRange) ~=  90 ...
 | min(Settings.LonRange) ~= -180 || max(Settings.LonRange) ~= 180

  %find the min and max of lon and lat for each profile
  Limits = [min(Data.Lon,[],2),max(Data.Lon,[],2), ...
            min(Data.Lat,[],2),max(Data.Lat,[],2)];

  %generate a list of profiles where any point falls inside the limits
  Bad = [];
  Bad = [Bad;find(Limits(:,1) > max(Settings.LonRange))];
  Bad = [Bad;find(Limits(:,2) < min(Settings.LonRange))];
  Bad = [Bad;find(Limits(:,3) > max(Settings.LatRange))];
  Bad = [Bad;find(Limits(:,4) < min(Settings.LatRange))];

  Good = 1:1:size(Limits,1); Good(Bad) = [];
  Data = reduce_struct(Data,Good,VarsToIgnore,1);

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


%end of programme
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to validate list of allowed instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validateStringParameter(varargin)
validatestring(varargin{:});
end



