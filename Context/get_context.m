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
%
%
%planning to add:
%  A. tropopause height
%  B. stratopause height
%  C. IMERG convection
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL GEOSPATIAL INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%  *  TimePoints             (double,         NaN)  Same size as LonPoints.    Required for output options marked with a *, in Matlab units
%  ^  Pressure               (double,         NaN)  1D array of levels in hPa. Required for output options marked with a ^, in Matlab units
%
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
%
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
%     Sentinel_Reload        (logical,       true)  Load image in Sentinel_OutFile rather than downloading
%     SurfaceImage_Image     (char,  'HRNatEarth')  Low-res surface image to use
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% required functions (and storage location at time of writing):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  nph_getnet.m   - https://github.com/corwin365/MatlabFunctions/blob/master/FileHandling/netCDF/nph_getnet.m
%  map_tessa.m    - https://github.com/corwin365/MatlabFunctions/blob/master/DatasetSpecific/TessaDEM/map_tessa.m
%
%You will also need to create a function LocalDataDir.m which takes no inputs and returns a string representing
%the root directory of our data storage hierarchy. On eepc-0184, this means it should return the string '/data1/Hub/'.
%It can instead set to return an empty string if all files paths used are set manually as options.
% For a more complex multi-system example see https://github.com/corwin365/MatlabFunctions/blob/master/System/LocalDataDir.m
%
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
addParameter(p,'TimePoints', NaN,@(x) validateattributes(x,{'numeric'},{'size',size(LonPoints)})); %time of each point, in Matlab units
addParameter(p,'Pressure',   NaN,@(x) validateattributes(x,{'numeric'},{'<=',1200}));              %pressure levels for output, in hPa

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

%other options
addParameter(p,'HighResTopo_LRFill',    true,         @islogical); %fill high-res topo using using low-res topography if needed
addParameter(p,'HighResTopo_TileScript',false,        @islogical); %return an SCP script to get the tiles needed for the high-res topo option from eepc-0184
addParameter(p,'Sentinel_ID',           {'',''},      @iscell  );  %sentinel API username   and password
addParameter(p,'Sentinel_Reload',       true,         @islogical); %reuse downloaded Sentinel imagery if it exists
addParameter(p,'Sentinel_OutFile',      'out.png',    @ischar);     %file to write Sentinel image out to
addParameter(p,'SurfaceImage_Image',    'HRNatEarth', @ischar);     %file to write Sentinel image out to

%paths
addParameter(p,'Era5_Path',         [LocalDataDir,'/ERA5/'],                                                @ischar); %path to ERA5 data
addParameter(p,'LowResTopo_Path',   [LocalDataDir,'/topography/easy_tenth_degree_topography/easy_topo.mat'],@ischar); %path to data
addParameter(p,'HighResTopo_Path',  [LocalDataDir,'/topography/tessaDEM/raw/'],                             @ischar); %path to data
addParameter(p,'Indices_Path',      [LocalDataDir,'/Miscellany/'],                                          @ischar); %path to climate index data
addParameter(p,'SurfaceImage_Path', [LocalDataDir,'/topography/'],                                          @ischar); %path to surface imagery

%done - parse and restructure inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parse(p,LonPoints,LatPoints,varargin{:})
Settings = p.Results;
TimePoints = Settings.TimePoints; Settings = rmfield(Settings,'TimePoints');
clear p varargin

%override actual options if 'EveryThing' is set
if Settings.Everything == true
  Settings.LowResTopo   = true;
  Settings.HighResTopo  = true;
  Settings.Wind         = true;
  Settings.Indices      = true;
  Settings.Sentinel     = true;
  Settings.SurfaceImage = true;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% low-res topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.LowResTopo == true

  %first check if the input data file exists
  if ~exist(Settings.LowResTopo_Path,'file')
    warning('LowResTopo: easyTopo data not located, skipping.')
  else
    
    %load the data, create an interpolant, and put it on output grid
    EasyTopo = load(Settings.LowResTopo_Path);
    I = griddedInterpolant(EasyTopo.topo.lats,EasyTopo.topo.lons,EasyTopo.topo.elev);
    Output.LowResTopo = I(LatPoints,LonPoints);

    clear I EasyTopo
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% high-res topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.HighResTopo == true
  
  [Alt,~,~,TileScript] = map_tessa(LonPoints,LatPoints, ...
                                   'ETFill',     Settings.HighResTopo_LRFill,     ...
                                   'DataDir',    Settings.HighResTopo_Path,       ...
                                   'TileScript', Settings.HighResTopo_TileScript, ...
                                   'ETPath',     Settings.LowResTopo_Path);
  Output.HighResTopo = Alt;
  Output.TileScript  = TileScript;
  clear Alt TileScript HRTRes
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Wind == true
  %check we fed in a time - this is required
  if sum(isnan(TimePoints)) == numel(TimePoints);
    warning('Wind: no TimePoints provided. Skipping.')
  elseif isnan(Settings.Pressure)
    warning('Wind: no pressure levels provided. Skipping')
  else
    %ok, create an interpolant and grab the surface wind (1000hPa)
    I = create_era5_interpolant(LonPoints,LatPoints,TimePoints,Settings,'Wind');
    if ~isa(I,'double'); 

      %create point arrays that have an extra pressure axis
      Lon  = repmat(LonPoints, [ones(ndims(LonPoints ),1);numel(Settings.Pressure)]');
      Lat  = repmat(LatPoints, [ones(ndims(LonPoints ),1);numel(Settings.Pressure)]');
      Time = repmat(TimePoints,[ones(ndims(TimePoints),1);numel(Settings.Pressure)]');
      P    = repmat(permute(Settings.Pressure',[2:ndims(LonPoints)+1,1]),[size(LonPoints),1]);

      %interpolate the data to the points
      Output.U = I.U(Lon,Lat,Time,P);
      Output.V = I.V(Lon,Lat,Time,P);    
    end
    clear I Lon Lat Time P
    
  end;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sentinel API cloudless imagery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.Sentinel == true

  %several checks to do here, so let's set a flag to keep track
  Fail = 0;

  %do we already have a file?
  if Settings.Sentinel_Reload == true

    %load the file into memory
    if exist(Settings.Sentinel_OutFile,'file')
      Output.Sentinel = flipud(imread(Settings.Sentinel_OutFile));
      %check the size matches
      if ~isequal(size(Output.Sentinel),size(repmat(LonPoints,1,1,3)));
        warning('Sentinel: previously-downloaded image is not the right size, getting new data')
        delete(Settings.Sentinel_OutFile)
      else
        warning('Sentinel: reusing previously-downloaded image')
        Fail = 1; %because we don't need to get the data
      end
    else
        warning('Sentinel: no previously-downloaded image, getting new data')
    end
  end

  %the data MUST be a valid and regular rectangle. Check this first.
  %we're looking for a single unique diff() value in each dimension
  %and exactly two dimensions, but we need to do some extra work to deal
  %with numerical precision issues
  if ndims(LonPoints) ~= 2
    warning('Sentinel: data points must form a regularly-spaced rectangle to get Sentinel data; shape failed. Skipping.')
    Fail = 1;
  end

  a = unique(diff(LonPoints,1,1)); b = unique(diff(LonPoints,1,2));
  c = unique(diff(LatPoints,1,1)); d = unique(diff(LatPoints,1,2));
  ar = range(a)./a; br = range(b)./b;  cr = range(c)./c; dr = range(d)./d;
  e = [ar;br;cr;dr]; e(isnan(e)) = [];
  if max(e) > 1e-5;
    warning('Sentinel: data points must form a regularly-spaced rectangle to get Sentinel data; diff failed. Skipping.')
    Fail = 1;
  end
  clear ar br cr dr e c d

  %check we fed in username and password - required
  if numel(Settings.Sentinel_ID{1}) == 0 | numel(Settings.Sentinel_ID{2}) == 0;
    warning('Sentinel: no Sentinel user credentials supplied. Skipping.')
    Fail = 1;
  end

  %getting better. Still not there, keep going...

  %we can only request 2500x2500 points. Is this true?
  if size(LonPoints,1) > 2500 | size(LonPoints,2) > 2500;
    warning('Sentinel: API is restricted to max 2500 points in an dimension, skipping.')
    Fail = 1;
  end


  %we need a bounding box for the region...
  BBox = [min(LonPoints,[],'all'),min(LatPoints,[],'all'), ...
    max(LonPoints,[],'all'),max(LatPoints,[],'all')];

  %... and a number of points. This requires working out if our data are lat-major or lon-major
  % if lon is the x-axis, then max(b) will be greater than max(a)
  if max(b) > max(a); Resolution = size(LonPoints');
  else                Resolution = size(LonPoints);
  end

  %hence check if we're in the permitted range of resolutions
  boxwidthx = deg2km(distance(BBox(4),BBox(1),BBox(4),BBox(3),'degrees')).*1000;
  boxwidthy = deg2km(distance(BBox(2),BBox(1),BBox(4),BBox(1),'degrees')).*1000;
  resx = boxwidthx./Resolution(1);
  resy = boxwidthy./Resolution(2);
  if resx > 1600; Fail = 1; warning('Sentinel: longitude resolution is too coarse. Skipping'); end
  if resy > 1600; Fail = 1; warning('Sentinel: latitude resolution is too coarse. Skipping'); end
  clear boxwidthx boxwidthy resx resy

  if Fail == 0; %just to stop the warning coming up if we've already failed
    %finally, if the resolution if REALLY high, make sure the user really wants this
    if numel(LonPoints) > 1000*1000;
      Input = input(['Sentinel: this request will query the Sentinel API for ',num2str(numel(LonPoints)),' points. Are you certain? Enter 1 to confirm.']);
      if Input ~= 1; Fail = 1; end
    end
  end

  %ok, let's go
  if Fail == 0;


    
    %downloading uses the Python API for Copernicus. This took me hours to figure out the syntax for. 

    %generate the API calling script
    Script = sentinel_script(Settings.Sentinel_ID,Settings.Sentinel_OutFile,BBox,Resolution)';
    ScriptFile = "get_sentinel_"+strrep(num2str(datenum(now)),'.','')+".py";
    writelines(Script,ScriptFile)

    %now run the script, and delete it
    pyrunfile(ScriptFile);
    delete(ScriptFile)

    %and load the image into memory
    Output.Sentinel = flipud(imread(Settings.Sentinel_OutFile));


  end

  clear Fail BBox Script ScriptFile a b Resolution Input

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface imagery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create ERA5 interpolant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = create_era5_interpolant(LonPoints,LatPoints,TimePoints,Settings,Prefix)

%fallback
I = NaN;

%find all unique days, as we store ECMWF data in daily files
Days = unique(floor(TimePoints));

%load them all up
for iDay=1:1:numel(Days)
  [y,~,~] = datevec(Days(iDay));
  dn = date2doy(Days(iDay));
  FilePath = [Settings.Era5_Path,sprintf('%04d',y),'/era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
  if ~exist(FilePath,'file')
    warning([Prefix,': Cannot find ERA5 file for ',datestr(Days(iDay)),'; skipping.'])
    continue
  end

  %load file and extract surface U and V
  E5 = rCDF(FilePath);

  if ~exist('Store','var');
    Store.U = permute(E5.u,[4,3,1,2]);
    Store.V = permute(E5.v,[4,3,1,2]);;
    Store.t = Days(iDay) + linspace(0,1,size(E5.u,1));
    Store.Lon = E5.longitude;
    Store.Lat = E5.latitude;
  else
    Store.U = cat(3,Store.U,permute(E5.u,[4,3,1,2]));
    Store.V = cat(3,Store.V,permute(E5.v,[4,3,1,2]));
    Store.t = [Store.t,Days(iDay) + linspace(0,1,size(E5.u,1))];
  end

end; clear iDay
clear y dn FilePath E5 Days

if ~exist('Store','var');
  warning([Prefix,': no ERA5 data found for any dates, output variable will not be returned.'])
  return
end

%get pressure axis
Store.P = ecmwf_prs_v3(size(Store.U,4));

%make sure data ascends monotonically
[Store.Lon,idx] = sort(Store.Lon,'ascend'); Store.U = Store.U(idx,:,:,:); Store.V = Store.V(idx,:,:,:);
[Store.Lat,idx] = sort(Store.Lat,'ascend'); Store.U = Store.U(:,idx,:,:); Store.V = Store.V(:,idx,:,:);
[Store.t,  idx] = sort(Store.t,  'ascend'); Store.U = Store.U(:,:,idx,:); Store.V = Store.V(:,:,idx,:);
[Store.P,  idx] = sort(Store.P,  'ascend'); Store.U = Store.U(:,:,:,idx); Store.V = Store.V(:,:,:,idx);

%create interpolant, and interpolate to the requested locations
clear I
I.U = griddedInterpolant({Store.Lon,Store.Lat,Store.t,Store.P},Store.U,'linear','linear');
I.V = griddedInterpolant({Store.Lon,Store.Lat,Store.t,Store.P},Store.V,'linear','linear');


return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get ERA5 pressure levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Pressure = ecmwf_prs_v3(NLevs,LnSP)

%do we have log surface pressure, or do we assume 1000hPa?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2;LnSP   = log(1000.*100); end

%identify A and B coefficients for this number of levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [0,2.000365,3.102241,4.666084,6.827977,9.746966,13.605424,18.608931,24.985718,32.98571,42.879242,54.955463,69.520576,86.895882,107.415741,131.425507,159.279404,191.338562,227.968948,269.539581,316.420746,368.982361,427.592499,492.616028,564.413452,643.339905,729.744141,823.967834,926.34491,1037.201172,1156.853638,1285.610352,1423.770142,1571.622925,1729.448975,1897.519287,2076.095947,2265.431641,2465.770508,2677.348145,2900.391357,3135.119385,3381.743652,3640.468262,3911.490479,4194.930664,4490.817383,4799.149414,5119.89502,5452.990723,5798.344727,6156.074219,6526.946777,6911.870605,7311.869141,7727.412109,8159.354004,8608.525391,9076.400391,9562.682617,10065.97852,10584.63184,11116.66211,11660.06738,12211.54785,12766.87305,13324.66895,13881.33106,14432.13965,14975.61523,15508.25684,16026.11523,16527.32227,17008.78906,17467.61328,17901.62109,18308.43359,18685.71875,19031.28906,19343.51172,19620.04297,19859.39063,20059.93164,20219.66406,20337.86328,20412.30859,20442.07813,20425.71875,20361.81641,20249.51172,20087.08594,19874.02539,19608.57227,19290.22656,18917.46094,18489.70703,18006.92578,17471.83984,16888.6875,16262.04688,15596.69531,14898.45313,14173.32422,13427.76953,12668.25781,11901.33984,11133.30469,10370.17578,9617.515625,8880.453125,8163.375,7470.34375,6804.421875,6168.53125,5564.382813,4993.796875,4457.375,3955.960938,3489.234375,3057.265625,2659.140625,2294.242188,1961.5,1659.476563,1387.546875,1143.25,926.507813,734.992188,568.0625,424.414063,302.476563,202.484375,122.101563,62.78125,22.835938,3.757813,0,0];
B = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000007,0.000024,0.000059,0.000112,0.000199,0.00034,0.000562,0.00089,0.001353,0.001992,0.002857,0.003971,0.005378,0.007133,0.009261,0.011806,0.014816,0.018318,0.022355,0.026964,0.032176,0.038026,0.044548,0.051773,0.059728,0.068448,0.077958,0.088286,0.099462,0.111505,0.124448,0.138313,0.153125,0.16891,0.185689,0.203491,0.222333,0.242244,0.263242,0.285354,0.308598,0.332939,0.358254,0.384363,0.411125,0.438391,0.466003,0.4938,0.521619,0.549301,0.576692,0.603648,0.630036,0.655736,0.680643,0.704669,0.727739,0.749797,0.770798,0.790717,0.809536,0.827256,0.843881,0.859432,0.873929,0.887408,0.8999,0.911448,0.922096,0.931881,0.94086,0.949064,0.95655,0.963352,0.969513,0.975078,0.980072,0.984542,0.9885,0.991984,0.995003,0.99763,1];

%compute half-levels
%%%%%%%%%%%%%%%%%%%%%

A2 = repmat(squeeze(A)',[1,size(LnSP)]);
B2 = repmat(squeeze(B)',[1,size(LnSP)]);
C = permute(repmat(exp(LnSP),[ones(ndims(LnSP),1);NLevs+1]'),[ndims(LnSP)+1,1:ndims(LnSP)]);

PHalf = A2 + B2.*C;
clear A2 B2 C

%compute pressure
%%%%%%%%%%%%%%%%%

%compute
Pressure = PHalf + circshift(PHalf,1,1);
Pressure = 0.5.*Pressure./100;

%drop bad point
Pressure = reshape(Pressure,size(Pressure,1),[]);
Pressure = Pressure(2:end,:);
sz = size(PHalf); sz(1) = size(Pressure,1);
Pressure = reshape(Pressure,sz);

%put pressure dimension at end
Pressure = permute(Pressure,[2:ndims(Pressure),1]);

%if 1d, flatten along primary axis
if ndims(Altitude) == 2;
  if size(Altitude,1) == 1;
    Altitude = squeeze(Altitude)';
    Pressure = squeeze(Pressure)';
  end
end

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% syntactical shortenings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MinMax = minmax(Array)
MinMax = [nanmin(Array(:)),nanmax(Array(:))];
return

function InRange = inrange(Array,MinMax,NoEnds)
InRange = find(Array >  min(MinMax) & Array <  max(MinMax));
return








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write Sentinel python script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Script = sentinel_script(Sentinel_ID,Sentinel_OutFile,BBox,Resolution)

Script        = "from scipy.io import savemat";
Script(end+1) = "import numpy as np";
Script(end+1) = "from oauthlib.oauth2 import BackendApplicationClient";
Script(end+1) = "from requests_oauthlib import OAuth2Session";
Script(end+1) = "";
Script(end+1) = "";
Script(end+1) = "# Your client credentials";
Script(end+1) = "client_id = '"+Sentinel_ID{1}+"'";
Script(end+1) = "client_secret = '"+Sentinel_ID{2}+"'";
Script(end+1) = "";
Script(end+1) = "# Create a session";
Script(end+1) = "client = BackendApplicationClient(client_id=client_id)";
Script(end+1) = "oauth = OAuth2Session(client=client)";
Script(end+1) = "";
Script(end+1) = "# Get token for the session";
Script(end+1) = "token = oauth.fetch_token(token_url='https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token',";
Script(end+1) = "                          client_secret=client_secret, include_client_id=True)";
Script(end+1) = "";
Script(end+1) = "";
Script(end+1) = "response = oauth.get('https://sh.dataspace.copernicus.eu/configuration/v1/wms/instances')";
Script(end+1) = "";
Script(end+1) = "evalscript = '''";
Script(end+1) = "//VERSION=3";
Script(end+1) = "function setup() {";
Script(end+1) = "  return {";
Script(end+1) = "    input: ['B02', 'B03', 'B04'],";
Script(end+1) = "    output: { bands: 3 }";
Script(end+1) = "  };";
Script(end+1) = "}";
Script(end+1) = "";
Script(end+1) = "function evaluatePixel(sample) {";
Script(end+1) = "  return [2.5 * sample.B04/10000, 2.5 * sample.B03/10000, 2.5 * sample.B02/10000];";
Script(end+1) = "}";
Script(end+1) = "'''";
Script(end+1) = "";
Script(end+1) = "request = {";
Script(end+1) = "  'input': {";
Script(end+1) = "    'bounds': {";
Script(end+1) = "      'bbox': [";
Script(end+1) = "        "+BBox(1)+",";
Script(end+1) = "        "+BBox(2)+",";
Script(end+1) = "        "+BBox(3)+",";
Script(end+1) = "        "+BBox(4)+"";
Script(end+1) = "      ]";
Script(end+1) = "    },";
Script(end+1) = "    'data': [";
Script(end+1) = "      {";
Script(end+1) = "        'dataFilter': {";
Script(end+1) = "          'timeRange': {";
Script(end+1) = "            'from': '2023-01-01T00:00:00Z',";
Script(end+1) = "            'to': '2023-01-02T23:59:59Z'";
Script(end+1) = "          }";
Script(end+1) = "        },";
Script(end+1) = "        'type': 'byoc-5460de54-082e-473a-b6ea-d5cbe3c17cca'";
Script(end+1) = "      }";
Script(end+1) = "    ]";
Script(end+1) = "  },";
Script(end+1) = "  'output': {";
Script(end+1) = "    'width': "+Resolution(1)+",";
Script(end+1) = "    'height': "+Resolution(2)+",";
Script(end+1) = "    'responses': [{'format': {'type': 'image/png'}}],";
Script(end+1) = "  },";
Script(end+1) = "  'evalscript': evalscript,";
Script(end+1) = "}";
Script(end+1) = "";
Script(end+1) = "url = 'https://sh.dataspace.copernicus.eu/api/v1/process'";
Script(end+1) = "response = oauth.post(url, json=request)";
Script(end+1) = "";
Script(end+1) = "if response.status_code == 400:";
Script(end+1) = "  print(response.text)";
Script(end+1) = "";
Script(end+1) = "with open('"+Sentinel_OutFile+"','wb') as f:";
Script(end+1) = "  f.write(response.content)";
Script(end+1) = "";


return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read netCDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FileContents = rCDF(FilePath,Method)

%run getnet
Data = nph_getnet(FilePath);

%move data up to top level
FileContents = Data.Data;

%move metadata to subsidiary level
MetaFields = {'Filename','Name','Dimensions','Variables','Attributes','Groups','Format'};
for iField=1:1:numel(MetaFields); FileContents.MetaData.(MetaFields{iField}) = Data.(MetaFields{iField}); end

%reorder dimensions so that they're in netCDF order (i.e. reverse them)
Fields = fieldnames(FileContents);
for iField=1:1:numel(Fields)
  if strcmp(Fields{iField},'MetaData'); continue;
  else;
    sz = size(FileContents.(Fields{iField}));
    if numel(sz) > 2 | (numel(sz) == 2 & (sz) == 1)
      FileContents.(Fields{iField}) = permute(FileContents.(Fields{iField}),numel(sz):-1:1);
    end
  end
end; clear iField

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% date2doy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [doy,fraction] = date2doy(inputDate)
% Author: Anthony Kendall
% Contact: anthony [dot] kendall [at] gmail [dot] com
% Created: 2008-03-11
% Copyright 2008 Michigan State University.
%modified 2024/04/28 to remove unneeded outputs

%Want inputs in rowwise format
[doy,fraction] = deal(zeros(size(inputDate)));
inputDate = inputDate(:);

%Parse the inputDate
[dateVector] = datevec(inputDate);

%Set everything in the date vector to 0 except for the year
dateVector(:,2:end) = 0;
dateYearBegin = datenum(dateVector);

%Calculate the day of the year
doyRow = inputDate - dateYearBegin;

%Fill appropriately-sized output array
doy(:) = doyRow;
if flagFrac
    fraction(:) = fracRow;
end