function Output = get_context(LonPoints,LatPoints,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unified geophysical context loader for datasets stored on eepc-0184
%
%Can currently load:
%
%  1. 'LowResTopo'   - 0.1 degree topography from easyTopo
%                    - requires easy_topo.mat, path can be set in this function
%  2. 'HighResTopo'  - 30m topography from TessaDEM
%                    - requires TESSA data tiles. This function can generate an SCP script to download them.
%  3. 'Wind'         - ERA5 1.5 degree resolution U and V at chosen pressure levels (can also do 
%                    - requires ERA5 netCDF data files as outputted by CDS API.
%                    - output will have an extra dimension corresponding to the pressure levels requested
%  4. 'Indices'      - climate indices: 'QBO','ENSO','JetFuelPrice','NAM','NAO','TSI','SeaIce','AMO'
%                    - output calculated based on time only, lat and lon must be set but will be ignored
%  5. 'Sentinel'     - downloads and imports high-res surface imagery from the Quarterly Cloudless Sentinel-2 Mosaics.
%                    - BE CAREFUL with this option - use is metered on a monthly basis. 
%  6. 'SurfaceImage' - surface imagery from stored global files (coarser than Sentinel, but less resource-intensive)
%                    - by default uses 0.1 degree Natural Earth map.
%                    - Other options: 'GreyScale', 'Modis','NatEarth','HRNatEarth','HRNatEarthBright', 
%                                     'land_ocean_ice', 'pale','land_ocean_ice_cloud','faded'
%  7. 'Pauses'       - computes tropopause and stratopause PRESSURE (i.e. in hPa) from 1.5 degree 3-hour ERA5 data
%
%
%By default the routine will return no useful data. Any chosen outputs
%must be switched on with flags.
%
%Corwin Wright, c.wright@bath.ac.uk, 28/04/APR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIRED INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     VarName            (type)      description
%     ------------------------------------------------------------------------------------------- 
%     LonPoints          (double)    Longitude points to return data for. 
%                                    Must have the same number of points as LatGrid.
%
%     LatPoints          (double)    Latitude  points to return data for. 
%                                    Must have the same number of points as LonGrid.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL GEOSPATIAL INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%  *  TimePoints             (double,         NaN)  Same size as LonPoints.    Required for output options marked with a *, in Matlab units
%  *  SingleTime             (double,         NaN)  Alternative To TimePoints:  a single time which will be duplicated out to the same size as LonPoints.
%  ^  Pressure               (double,         NaN)  1D array of levels in hPa. Required for output options marked with a ^.
%     PressurePoints         (double,         NaN)  Same size as LonPoints, as an alternative to Pressure levels specified as above.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%     Everything             (logical,      false)  try and return all the below. Overrides all individual choices.
%     LowResTopo             (logical,      false)  return easyTopo 0.1 degree topography data (faster)
%     HighResTopo            (logical,      false)  return TessaDEM 30m topography data        (slower)
%  *^ Wind                   (logical,      false)  return winds from ERA5
%  *  Indices                (logical,      false)  return climate indices
%     Sentinel               (logical,      false)  return Quarterly cloudless Sentinel-2 mosaics, ~10m resolution
%     SurfaceImage           (logical,      false)  returns lower-resolution surface imagery
%  *  Pauses                 (logical,      false)  returns tropopause and stratopause pressure computed from ERA5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%     HighResTopo_LRFill     (logical,       true)  fill gaps (poles and oceans) in Tessa data with easytopo data. 
%     HighResTopo_TileScript (logical       false)  generate SCP script to download required Tessa tiles
%     Sentinel_ID            (2-elmt  cell, empty)  Copernicus client ID/password
%     Sentinel_OutFile       (char,     'out.png')  Image file to write Sentinel data to
%     Sentinel_Reuse         (logical,       true)  Load image in Sentinel_OutFile rather than downloading
%     Sentinel_Gain          (numeric,          5)  Increase Sentinel image gain (larger is brighter, but may saturate)
%     SurfaceImage_Image     (char,  'HRNatEarth')  Low-res surface image to use
%     Era5_Clim              (logical,      false)  use ERA5 1991-2020 climatological data rather than the true date. Takes day-of-year from supplied TimePoints and ignores year supplied.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PATHS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%     LowResTopo_Path        (char, see in parser)  path to easytopo data file
%     HighResTopo_Path       (char, see in parser)  path to TessaDEM data files
%     SurfaceImage_Path      (char, see in parser)  Path to surface image file
%     Era5_Path              (char, see in parser)  path to ERA5 data, used for Winds
%     Indices_Path           (char, see in parser)  path to directory containing climate index data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% required functions  - and storage location at time of writing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nph_getnet.m   - https://github.com/corwin365/MatlabFunctions/blob/master/FileHandling/netCDF/nph_getnet.m
%  map_tessa.m    - https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/TessaDEM/map_tessa.m
%  p2h.m          - https://github.com/corwin365/MatlabFunctions/blob/master/GeophysicalAndTime/p2h.m
%  h2p.m          - https://github.com/corwin365/MatlabFunctions/blob/master/GeophysicalAndTime/h2p.m
%  date2doy.m     - https://github.com/corwin365/MatlabFunctions/blob/master/GeophysicalAndTime/date2doy.m
%
%If you have your data stored in the same directory hierarchy as on eepc-0184, you may find it useful to create a
%function "LocalDataDir.m" which takes no inputs and returns a CHAR representing the root directory of the data 
%storage hierarchy [on eepc-0184, this means it should return the string '/data1/Hub/']. This will let this programme
%compute the paths to all the data files used here automatically. If this function does not exist, or your directory 
%hierarchy differs, you need to set all files paths manually using the optional inputs described above.
%For an example see https://github.com/corwin365/MatlabFunctions/blob/master/System/LocalDataDir.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

%lat/lon grids, required
addRequired(p,'LonPoints',@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180}))
addRequired(p,'LatPoints',@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90,'size',size(LonPoints)}))

%variables used for many, but not all, datasets
addParameter(p,'TimePoints',     NaN,@(x) validateattributes(x,{'numeric'},{'size',size(LonPoints)})); %time of each point, in Matlab units
addParameter(p,'SingleTime',     NaN,@isnumeric);                                                      %single time, as an alternative to TimePoints
addParameter(p,'Pressure',       NaN,@(x) validateattributes(x,{'numeric'},{'<=',1200}));              %EITHER pressure levels to find for every lat/lon point OR pressure levels associated with each point, in hPa. Assumed to be the former unless 'PressureAsPoints' is set to true 

%individual datasets
%%%%%%%%%%%%%%%%%%%%%

%top-level output options
addParameter(p,'Everything',  false,@islogical); %try to load all the below
addParameter(p,'LowResTopo',  false,@islogical); %load easyTopo 0.1 degree topography
addParameter(p,'HighResTopo', false,@islogical); %load TessaDEM 30m topography
addParameter(p,'Wind',        false,@islogical); %load winds from 1.5 degree ERA5
addParameter(p,'Indices',     false,@islogical); %load climate indices
addParameter(p,'Sentinel',    false,@islogical); %download and load Sentinel surface imagery. Requires ID and Password, set via Sentinel_ID 
addParameter(p,'SurfaceImage',false,@islogical); %load surface imagery
addParameter(p,'Pauses',      false,@islogical); %compute tropopause and stratopause from ERA5 data

%other options
addParameter(p,'PressureAsPoints',      false,             @islogical); %flag to provide input pressure as corresponding list of points rather than jsut requesting levels - useful for e.g. 3D traces through a wind field
addParameter(p,'HighResTopo_LRFill',    true,              @islogical); %fill high-res topo using using low-res topography if needed
addParameter(p,'HighResTopo_TileScript',false,             @islogical); %return an SCP script to get the tiles needed for the high-res topo option from eepc-0184
addParameter(p,'Sentinel_ID',           {'',''},           @iscell  );  %sentinel API username   and password
addParameter(p,'Sentinel_Reuse',        true,              @islogical); %reuse downloaded Sentinel imagery if it exists
addParameter(p,'Sentinel_OutFile',      'out.png',         @ischar);    %file to write Sentinel image out to
addParameter(p,'Sentinel_Gain',         5,                 @isnumeric); %gain for Sentinel data
addParameter(p,'SurfaceImage_Image',    'HRNatEarth',      @ischar);    %file to write Sentinel image out to
addParameter(p,'Era5_Clim',             false,             @islogical); %use the 1991-2020 ERA5 climatology generated by Time, rather than true ERA5

%paths
%first, check if we have a LocalDataDir file, and set the data root to '/' if not. Then, specify default paths.
try; LD = LocalDataDir; catch; LD = './'; end
addParameter(p,'Era5_Path',         [LD,'/ERA5/'],                                                @ischar); %path to ERA5 data
addParameter(p,'LowResTopo_Path',   [LD,'/topography/easy_tenth_degree_topography/easy_topo.mat'],@ischar); %path to easyTopo data
addParameter(p,'HighResTopo_Path',  [LD,'/topography/tessa'],                                     @ischar); %path to TessaDEM data
addParameter(p,'Indices_Path',      [LD,'/Miscellany/'],                                          @ischar); %path to climate index data
addParameter(p,'SurfaceImage_Path', [LD,'/topography/'],                                          @ischar); %path to surface imagery
clear LD

%done - parse and restructure inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse the array
parse(p,LonPoints,LatPoints,varargin{:})
Settings = p.Results;

%sanity/safety checks that apply at top level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if we requested the SingleTime option, check we didn't also specify TimePoints, then apply
if ~isnan(Settings.SingleTime)
  if nansum(Settings.TimePoints) == 0
    %go ahead and replace
    Settings.TimePoints = ones(size(Settings.LonPoints)).*Settings.SingleTime;
  end
end
%oand now pull timepoints out to be a top-level variable
TimePoints = Settings.TimePoints; Settings = rmfield(Settings,'TimePoints');
clear p varargin

%if we want PressureAsPoints, check we have the right number of points
if Settings.PressureAsPoints == 1
  if numel(Settings.Pressure) ~= numel(LatPoints)
    error('Pressure requested as points: Pressure array should be same size as Lat/Lon')
    return
  end
end

%override actual options if 'EveryThing' is set
if Settings.Everything == true
  Settings.LowResTopo   = true;
  Settings.HighResTopo  = true;
  Settings.Wind         = true;
  Settings.Indices      = true;
  Settings.Sentinel     = true;
  Settings.SurfaceImage = true;
  Settings.Pauses       = true;
  warning('"Everything" option set - all output options will be attempted')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialise output struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Output = struct();
Output.Lon = LonPoints;
Output.Lat = LatPoints;
Output.Time = TimePoints;
if ~isnan(Settings.Pressure); Output.Pressure = Settings.Pressure; end

%also compute a bounding box, needed for some datasets
BBox = [min(LonPoints,[],'all'),min(LatPoints,[],'all'), ...
        max(LonPoints,[],'all'),max(LatPoints,[],'all')];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sentinel API cloudless imagery
%done first as it may require user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.Sentinel == true
  Sentinel = module_sentinel(Settings,LonPoints,LatPoints,TimePoints,BBox);
  if numel(Sentinel) > 0;
    Output.Sentinel = Sentinel;
  end
  clear Sentinel
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.Indices == true

  %check we fed in a time - this is required
  if sum(isnan(TimePoints)) == numel(TimePoints);
    warning('Indices: no TimePoints provided. Skipping.')
  else

    %let's try to get them all
    Indices = {'QBO','ENSO','JetFuelPrice','NAM','NAO','TSI','SeaIce','AMO'};
    Root = Settings.Indices_Path; 

    for iIndex=1:1:numel(Indices)

      try
        switch Indices{iIndex}
          case 'QBO'
            QBO = load([Root,'/QBO.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(QBO.Time,QBO.QBO,TimePoints);
            clear QBO
          case 'ENSO'
            ENSO = load([Root,'/nino34.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(ENSO.Time,ENSO.Nino34,TimePoints);
            clear ENSO
          case 'Fuel'
            Fuel = load([Root,'/jet_fuel_price.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(Fuel.Time,Fuel.Price,TimePoints);
            clear Fuel
          case 'NAM'
            NAM = load([Root,'/daily_nam.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(NAM.Time,NAM.NAM,TimePoints);
            clear NAM
          case 'NAO'
            NAO = load([Root,'/nao.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(NAO.Time,NAO.NAO,TimePoints);
            clear NAO
          case 'SSTs'
            SSTs = load([Root,'/ssts.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(SSTs.Time,SSTs.SSTs,TTimePoints);
            clear SSTs
          case 'TSI'
            TSI = load([Root,'/tsi.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(TSI.Time,TSI.TSI,TimePoints);
            clear TSI
          case 'Time'
            Output.Indices.(Indices{iIndex}) = TimeScale;
          case 'SeaIce'
            SeaIce = readmatrix([Root,'/N_seaice_extent_daily_v3.0.csv']);
            t = datenum(SeaIce(:,1),SeaIce(:,2),SeaIce(:,3));
            Output.Indices.(Indices{iIndex}) = interp1(t,SeaIce(:,4),TimePoints);
            clear SeaIce t
          case 'AMO'
            AMO = load([Root,'/AMO.mat']);
            Output.Indices.(Indices{iIndex}) = interp1(AMO.Time,AMO.AMO,TimePoints);
            clear AMO
        end
      catch; warning(['Indices: error locating input data for ',Indices{iIndex},'; skipping.'])
      end
    end
  end
  clear Root Indices iIndex
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Wind == true

  %check we fed in a time - this is required
  if sum(isnan(TimePoints)) == numel(TimePoints);
    warning('Wind: no TimePoints provided. Skipping.')
  elseif isnan(Settings.Pressure)
    warning('Wind: no pressure levels provided. Skipping')
  else
    %ok, create an interpolant and grab the wind
    I = create_era5_interpolant(TimePoints,Settings,'Wind',BBox);

    if ~isa(I,'double'); 

      %create point arrays that have an extra pressure axis, if needed
      if Settings.PressureAsPoints == 1; 
        Lon  = LonPoints;
        Lat  = LatPoints;
        Time = TimePoints;
        P    = Settings.Pressure;
      else    
        Lon  = repmat(LonPoints, [ones(ndims(LonPoints ),1);numel(Settings.Pressure)]');
        Lat  = repmat(LatPoints, [ones(ndims(LonPoints ),1);numel(Settings.Pressure)]');
        Time = repmat(TimePoints,[ones(ndims(TimePoints),1);numel(Settings.Pressure)]');        
        P    = repmat(permute(Settings.Pressure',[2:ndims(LonPoints)+1,1]),[size(LonPoints),1]);
      end

      %interpolate the data to the points
      Output.U = I.U(Lon,Lat,Time,P);
      Output.V = I.V(Lon,Lat,Time,P);    
    end

    clear I Lon Lat Time P
    
  end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface imagery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.SurfaceImage == true;

  %get path to the file we want
  switch Settings.SurfaceImage_Image
    case 'GreyScale';            Path = '/imagery/greyscale.png';
    case 'Modis';                Path = '/imagery/MODIS_Map.jpg';
    case 'NatEarth';             Path = '/ne/rasterI/NE1_50M_SR_W.tif';
    case 'HRNatEarth';           Path = '/ne/rasterI/HYP_HR_SR_OB_DR.tif';
    case 'HRNatEarthBright';     Path = '/ne/rasterI/ne_bright.tif';
    case 'land_ocean_ice';       Path = '/imagery/land_ocean_ice_8192.png';
    case 'land_ocean_ice_cloud'; Path = '/imagery/land_ocean_ice_cloud_8192.png';
    case 'faded';                Path = '/imagery/faded.jpg';
    case 'pale';                 Path = '/imagery/pale.png';
    otherwise                    Path = '';
  end
  Path = [Settings.SurfaceImage_Path,Path];

  if ~exist(Path,'file')
    warning('SurfaceImage: file not found, skipping')
  elseif ndims(LonPoints) > 2
    warning('SurfaceImage: cannot be requested in >2 dimensions, skipping')
  else
    
    %load the file
    Image.Map = flipud(imread(Path));

    %greyscale needs duplicating out to have three colours, even though they're all the same
    if strcmp(Settings.SurfaceImage_Image,'GreyScale'); Image.Map = repmat(Image.Map,1,1,3); end

    %create corresponding lat and lon arrays (assumes Mercator projection)
    Image.Lon = linspace(-180,180,1+size(Image.Map,2)); Image.Lon = Image.Lon(1:end-1);
    Image.Lat = linspace( -90,90,1+size(Image.Map,1)); Image.Lat = Image.Lat(1:end-1);

    %interpolate the image onto the desired output grid
    Out = NaN([size(LatPoints),3]);
    for iColour=1:1:3;
      I = griddedInterpolant({Image.Lon,Image.Lat},double(Image.Map(:,:,iColour))');
      Out(:,:,iColour) = I(LonPoints',LatPoints')';
    end

    %return and tidy
    Output.SurfaceImage = uint8(Out);
  end

  clear Path Image Out I iColour

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute tropo/stratopause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Pauses == true
  [Error,Tropopause,Stratopause] = module_pauses(Settings,LonPoints,LatPoints,TimePoints,BBox);
  if Error ~= 1;
    Output.Tropopause = Tropopause;
    Output.Stratopause = Stratopause;
  end
  clear Error Tropopause Stratopause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% low-res topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.LowResTopo == true; 
  Output.LowResTopo = module_lowrestopo(Settings,LonPoints,LatPoints);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% high-res topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.HighResTopo == true
  [Output.HighResTopo,TileScript] = module_highrestopo(Settings,LonPoints,LatPoints);
  if Settings.HighResTopo_TileScript == true; Output.TileScript  = TileScript; end
  clear TileScript
end

