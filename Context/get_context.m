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
%  ^  Pressure               (double,         NaN)  1D array of levels in hPa. Required for output options marked with a ^, in Matlab units
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
%     Sentinel_Reload        (logical,       true)  Load image in Sentinel_OutFile rather than downloading
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
addParameter(p,'Pauses',      false,@islogical); %compute tropopause and stratopause from ERA5 data

%other options
addParameter(p,'HighResTopo_LRFill',    true,              @islogical); %fill high-res topo using using low-res topography if needed
addParameter(p,'HighResTopo_TileScript',false,             @islogical); %return an SCP script to get the tiles needed for the high-res topo option from eepc-0184
addParameter(p,'Sentinel_ID',           {'',''},           @iscell  );  %sentinel API username   and password
addParameter(p,'Sentinel_Reload',       true,              @islogical); %reuse downloaded Sentinel imagery if it exists
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

%% low-res topography
%%%%%%%%%%%%%%%%%%%%%%

if Settings.LowResTopo == true; 
  Output.LowResTopo = module_lowrestopo(Settings,LonPoints,LatPoints);
end

%% high-res topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.HighResTopo == true
  [Output.HighResTopo,TileScript] = module_highrestopo(Settings,LonPoints,LatPoints);
  if Settings.HighResTopo_TileScript == true; Output.TileScript  = TileScript; end
  clear TileScript
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
    %ok, create an interpolant and grab the wind
    I = create_era5_interpolant(TimePoints,Settings,'Wind',BBox);

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


  %did we request a date at all? If not, use default (1st April 2023 - this is arbitrary)
  SentinelTime = TimePoints;
  if numel(SentinelTime) == 1 && isnan(SentinelTime); SentinelTime(:) = datenum(2023,4,1); end

  %what date do we want the data for? We must only be requesting one  
  if nanvar(SentinelTime,[],'all') ~= 0
    warning('Sentinel: requesting multiple dates for Sentinel data; using earliest-specified in valid range only')
    SentinelTime = nanmin(SentinelTime(SentinelTime > datenum(2023,1,1)));
  else
    SentinelTime = nanmin(SentinelTime,[],'all');
  end

  %is the requested date post-dataset start in 2023? If not, shift to the same DoY in 2023
  if nanmin(SentinelTime,[],'all') < datenum(2023,1,1);
    warning('Sentinel: requested dates before dataset start; shifting pre-2023 date into 2023')
    [y,~,~] = datevec(SentinelTime); dn = date2doy(SentinelTime);  y(y < 2023) = 2023;
    SentinelTime = datenum(y,1,dn);
    clear y dn
  end

  %work out the required resolution
  % This requires working out if our data are lat-major or lon-major
  % if lon is the x-axis, then max(b) will be greater than max(a)
  if ndims(LonPoints) > 2;
    warning('Sentinel: Sentinel data cannot be requested in >2 dimensions. Skipping.')
    Fail = 1;
  else
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
  end
  clear a b Resolution

  if Fail == 0; %so we skip this check if it's not needed, since it queries the user

    %we can only request 2500x2500 points in a single pass - otherwise we need to make multiple API calls
    

    if size(LonPoints,1) > 2500 | size(LonPoints,2) > 2500;  %use multiple passes

      %work out the properties of the multiple passes
      x = ceil(size(LonPoints,1)./2500); y = ceil(size(LonPoints,2)./2500);
      NPasses = x*y;
      BBoxes = NaN(NPasses,4);

      disp(['Sentinel: API allows max 2500 points per dimension, so this will require ',num2str(NPasses),' calls.'])
      Input = input(['Are you certain? Enter 1 to confirm.']);
      if Input ~= 1; Fail = 1; else
        k = 1;
        for iX=1:1:x
          for iY=1:1:y

            %define the bounding box for this pass
            idxX = [1,0]+[iX-1,iX].*2500;
            idxY = [1,0]+[iY-1,iY].*2500;
            idxX(idxX > size(LonPoints,1)) = size(LonPoints,1); idxX = idxX(1):1:idxX(end);
            idxY(idxY > size(LonPoints,2)) = size(LonPoints,2); idxY = idxY(1):1:idxY(end);

            BBoxes(k,:) = [LonPoints(idxX(  1),idxY(  1)), LatPoints(idxX(  1),idxY(  1)), ...
                           LonPoints(idxX(end),idxY(end)), LatPoints(idxX(end),idxY(end))];

            %identify which merged output points these pixels go into
            a = reshape(1:1:numel(LonPoints),size(LonPoints)); %list of ALL points
            a = a(idxX,:); a = a(:,idxY);
            IndexStore.(['store',num2str(k)]) = [idxX(1),idxX(end),idxY(1),idxY(end)];
            k = k+1;
          end
        end
      end
      clear k iX iY x y idxX idxY
      
    else
      %just use the master bounding box and list of output points
      NPasses = 1;
      BBoxes = BBox;
      IndexStore.(['store1']) = [1,size(LonPoints,1),1,size(LonPoints,2)];
    end
  end



  %ok, let's go
  if Fail == 0;
    for iPass = 1:1:NPasses
      %downloading uses the Python API for Copernicus. This took me hours to figure out the syntax for.

      %generate the API calling script
      Resolution = IndexStore.(['store',num2str(iPass)]); Resolution = Resolution([4,2]) - Resolution([3,1]) + [1,1]; 
      WorkingFileName = ['working',num2str(randi(1000,1,1)),'_',num2str(randi(1000,1,1)),'.png'];
      Script = sentinel_script(Settings.Sentinel_ID,WorkingFileName,Settings.Sentinel_Gain,BBoxes(iPass,:),Resolution,SentinelTime)';
      ScriptFile = "get_sentinel_"+strrep(num2str(datenum(now)),'.','')+".py";
      writelines(Script,ScriptFile)

      %now run the script, and delete it
      pyenv(ExecutionMode="OutOfProcess");
      pyrunfile(ScriptFile);
      delete(ScriptFile);

      %tidy up Python
      terminate(pyenv);

      %and load the image into memory
      SentinelData.(['store',num2str(iPass)]) = flipud(imread(WorkingFileName));
      delete(WorkingFileName)
      clear WorkingFileName

    end


  %merge passes then output
  Output.Sentinel = zeros([size(LonPoints),3],'uint8');
  for iPass=1:1:NPasses
    idx = IndexStore.(['store',num2str(iPass)]);
    Output.Sentinel(idx(1):idx(2),idx(3):idx(4),:) = SentinelData.(['store',num2str(iPass)]);
  end; 

  %also write a file image
  imwrite(flipud(Output.Sentinel),Settings.Sentinel_OutFile);
  

  end
  clear Fail Script ScriptFile a b Resolution Input  iPass idx BBoxes IndexStore NPasses SentinelData SentinelTime
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute tropopause and stratopause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.Pauses == true

  %check we fed in a time - this is required
  if sum(isnan(TimePoints)) == numel(TimePoints);
    warning('Pauses: no TimePoints provided. Skipping.')
  else

    %setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %get ERA5 global temperature
    I = create_era5_interpolant(TimePoints,Settings,'Pauses',BBox);

    %compute pressure. We can ignore lnsp as both 'pauses should be above the region it matters.
    Pressure = ecmwf_prs_v3(137);

    %stratopause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %method: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2011JD016893

    %interpolate the data to a regular ~1km height grid between 25 and 80km altitude 
    NewP = h2p(25:1:80)'; NewZ = p2h(NewP);
    lat = repmat( LatPoints,[ones(ndims(LatPoints),1);numel(NewP)]');
    lon = repmat( LonPoints,[ones(ndims(LatPoints),1);numel(NewP)]');
    t   = repmat(TimePoints,[ones(ndims(LatPoints),1);numel(NewP)]');
    p   = permute(repmat(NewP,[1,size(LatPoints)]),[2:ndims(lat),1]);
    T   = I.T(lon,lat,t,p);
    clear lon lat t p

    %reshape to put height first
    T  = permute(T ,[ndims(T ),1:1:ndims(T )-1]);    

    %smooth by 11km. Remember we don't know how many dimensions we have...
    sz = size(T);
    Ts = reshape(T,sz(1),prod(sz(2:end)));
    Ts = smoothdata(Ts,1,'movmean');
    Ts = reshape(Ts,sz);

    %find maximum in each profile
    [~,idx] = max(Ts,[],1);

    %check 5 levels above and below:
      %5 levels above must have -ve lapse rate
      %5 levels below must have +ve lapse rate
    dTdZ = diff(Ts,1,1);
    dTdZ = cat(1,zeros(size(idx)),dTdZ); %add extra level so points line up, rather than half-levels

    Stratopause = NaN(size(idx)); 
    for iProf=1:1:prod(size(idx));

      Above = idx(iProf)+1:1:idx(iProf)+5; Above = Above(Above > 0 & Above < size(NewP,1));
      Below = idx(iProf)-5:1:idx(iProf)-1; Below = Below(Below > 0 & Below < size(NewP,1));

      Above = -dTdZ(Above,iProf); Below = dTdZ(Below,iProf); %note - sign on Above

      if min(Above) > 0 & min(Below) > 0;
        %remove anything outside +/- 15 km from peak, for safety below
        T(NewZ < NewZ(idx(iProf))-15,iProf) = NaN;
        T(NewZ > NewZ(idx(iProf))+15,iProf) = NaN;

        %then find the maximum in the unsmoothed data
        [~, Stratopause(iProf)] = max(T(:,iProf),[],1);
      end
    end; 

    %convert to height, and fill small gaps from the pass condition above (these are usually <1% of the data)
    Stratopause(~isnan(Stratopause)) = NewZ(Stratopause(~isnan(Stratopause)));
    Stratopause = fillmissing(Stratopause,'linear');

    %drop unnecessary dimensions added above, and convert back to pressure
    Output.Stratopause = h2p(permute(Stratopause,[2:1:ndims(Stratopause),1]));

    %tidy up
    clear Above Below dTdZ idx NewP NewZ T Ts sz Stratopause

    %tropopause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %computes tropopause based on WMO definition, approximately following doi:10.1029/2003GL018240 but 
    %with modifications for speed

    %put temperature data onto native ERA5 pressure grid between 700 hPa and 10 hPa
    %we want to work up in height, so flip Pressure
    Pressure = Pressure(end:-1:1);
    lat = repmat( LatPoints,[ones(ndims(LatPoints),1);numel(Pressure)]');
    lon = repmat( LonPoints,[ones(ndims(LatPoints),1);numel(Pressure)]');
    t   = repmat(TimePoints,[ones(ndims(LatPoints),1);numel(Pressure)]');
    p   = permute(repmat(Pressure',[1,size(LatPoints)]),[2:ndims(lat),1]);
    T   = I.T(lon,lat,t,p);
    clear lon lat t p

    %reshape to put height first, then turn into lines
    T = permute(T ,[ndims(T ),1:1:ndims(T )-1]); 
    sz = size(T);
    T = reshape(T,sz(1),prod(sz(2:end)));

    %compute lapse rate. 
    dT = diff(T,1,1);
    dZ = diff([p2h(Pressure)]);
    Gamma = dT .* NaN;
    for iLev=1:1:numel(dZ)-1; Gamma(iLev,:) = dT(iLev,:)./dZ(iLev); end;
    clear dT dZ iLev

    %create an array to store our tropopause levels, then loop over the data to find them
    %we are working UPWARDS
    Tropopause = NaN(size(T,2),1);

    for iLev=1:1:numel(Pressure)

      %if pressure > 700hPa or <10hPa, or if we've already found the t'pause everywhere, skip
      if Pressure(iLev) > 700;           continue; end
      if Pressure(iLev) <  10;           continue; end
      if sum(isnan(Tropopause(:))) == 0; continue; end

      %check if Gamma is less than 2 anywhere at this level
      idx = find(Gamma(iLev,:) > -2);
      if numel(idx) == 0; continue; end %none at this level

      %remove any columns we already found
      Found = find(~isnan(Tropopause));
      [~,Remove] = intersect(idx,Found);
      idx(Remove) = [];
      clear Remove

      %for each element where the above criterion is met, check if the layer
      %2km higher also meets it

      %find which level is 2km above
      Z = p2h(Pressure(iLev));
      jLev = closest(p2h(Pressure),Z+2);
     
      %find all the columns where the criterion remains met ON AVERAGE for these 2km above
      Good = find(nanmean(Gamma(iLev:jLev,idx),1) > -2);
      if numel(Good) < 2 ; continue; end %this needs to be 2 because of ambiguity in the array operations below.
                                         %This leaves a small number of NaNs (<< 1%), which we interpolate over below
      idx= idx(Good);

      %for the remaining columns, find where the gradient crossed above -2 by linear interpolation
      G  = Gamma(iLev:jLev,idx);
      p  = linspace(Pressure(iLev),Pressure(jLev),100); %put onto 100 levels
      Gi = interp1(Pressure(iLev:jLev),Gamma(iLev:jLev,idx),p)+2;
      [~,minpidx] = min(abs(Gi),[],1);

      Tropopause(idx) = p(minpidx);
    end

    %fill small gaps due to put back to the original shape, and return
    Tropopause = fillmissing(Tropopause,'linear');    
    Output.Tropopause = reshape(Tropopause,sz(2:end));

    %tidy up
    clear I Pressure Found G Gamma Gi Good idx iLev iProf jLev minpidx p sz T Z Tropopause




    
  end

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

function I = create_era5_interpolant(TimePoints,Settings,Prefix,BBox)

%fallback
I = NaN;


%find all unique days, as we store ECMWF data in daily files
Days = unique(floor(TimePoints));

%load them all up
for iDay=1:1:numel(Days)
 

  %load daily data or climatological data as appropriate
  if Settings.Era5_Clim ~= true

    %true date data
    [y,~,~] = datevec(Days(iDay));
    dn = date2doy(Days(iDay));
    FilePath = [Settings.Era5_Path,sprintf('%04d',y),'/era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
    if ~exist(FilePath,'file')
      warning([Prefix,': Cannot find ERA5 file for ',datestr(Days(iDay)),'; skipping. Path tried: ',FilePath])
      continue
    end
  else
    %climatological data
    [~,m,d] = datevec(Days(iDay));
    FilePath = [Settings.Era5_Path,'/clim19912020/clim',sprintf('%02d',m),sprintf('%02d',d),'.nc'];
    if ~exist(FilePath,'file')
      warning([Prefix,': Cannot find climatology ERA5 file for ',datestr(Days(iDay),'dd/mmmm'),'; skipping. Path tried: ',FilePath])
      continue
    end
  end

  %load file and extract U,V and T
  E5 = rCDF(FilePath);

  %trim down to slightly larger than the bounding box
  idx.lat = inrange(E5.latitude, BBox([2,4])+[-1,1].*2);
  idx.lon = inrange(E5.longitude,BBox([1,3])+[-1,1].*2);

  E5.longitude = E5.longitude(idx.lon);
  E5.latitude  = E5.latitude( idx.lat);
  Vars = {'t','u','v'};
  for iVar=1:1:numel(Vars)
    Var = E5.(Vars{iVar});
    Var = Var(:,:,idx.lat,:);
    Var = Var(:,:,:,idx.lon);
    E5.(Vars{iVar}) = Var;
  end
  clear idx Vars Var iVar
    

  if ~exist('Store','var');
    Store.U = permute(E5.u,[4,3,1,2]);
    Store.V = permute(E5.v,[4,3,1,2]);
    Store.T = permute(E5.t,[4,3,1,2]);
    Store.t = Days(iDay) + linspace(0,1-(1./size(E5.u,1)),size(E5.u,1));
    Store.Lon = E5.longitude;
    Store.Lat = E5.latitude;
  else
    Store.U = cat(3,Store.U,permute(E5.u,[4,3,1,2]));
    Store.V = cat(3,Store.V,permute(E5.v,[4,3,1,2]));
    Store.T = cat(3,Store.T,permute(E5.t,[4,3,1,2]));
    Store.t = [Store.t,Days(iDay) + linspace(0,1-(1./size(E5.u,1)),size(E5.u,1))];
  end

end; clear iDay
clear y dn m d FilePath E5 Days

if ~exist('Store','var');
  warning([Prefix,': no ERA5 data found for any dates, output variable will not be returned.'])
  return
end

%get pressure axis
Store.P = ecmwf_prs_v3(size(Store.U,4));

%make sure data ascends monotonically
[Store.Lon,idx] = sort(Store.Lon,'ascend'); Store.U = Store.U(idx,:,:,:); Store.V = Store.V(idx,:,:,:); Store.T = Store.T(idx,:,:,:);
[Store.Lat,idx] = sort(Store.Lat,'ascend'); Store.U = Store.U(:,idx,:,:); Store.V = Store.V(:,idx,:,:); Store.T = Store.T(:,idx,:,:);
[Store.t,  idx] = sort(Store.t,  'ascend'); Store.U = Store.U(:,:,idx,:); Store.V = Store.V(:,:,idx,:); Store.T = Store.T(:,:,idx,:);
[Store.P,  idx] = sort(Store.P,  'ascend'); Store.U = Store.U(:,:,:,idx); Store.V = Store.V(:,:,:,idx); Store.T = Store.T(:,:,:,idx);

%create interpolant, and interpolate to the requested locations
clear I
I.U = griddedInterpolant({Store.Lon,Store.Lat,Store.t,Store.P},Store.U,'linear','linear');
I.V = griddedInterpolant({Store.Lon,Store.Lat,Store.t,Store.P},Store.V,'linear','linear');
I.T = griddedInterpolant({Store.Lon,Store.Lat,Store.t,Store.P},Store.T,'linear','linear');


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

function [Indices, Values] = closest(Values,Lookup)
[Values,Indices] = min(abs(Values-Lookup));
return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write Sentinel python script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Script = sentinel_script(Sentinel_ID,Sentinel_OutFile,Gain,BBox,Resolution,SentinelTime)

TimeStringStart = [datestr(SentinelTime-7,'yyyy-mm-dd'),'T00:00:00Z'];
TimeStringEnd   = [datestr(SentinelTime+7,'yyyy-mm-dd'),'T23:59:59Z'];


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
Script(end+1) = "  let gain = "+num2str(Gain);
Script(end+1) = "  return [gain * sample.B04/10000, gain * sample.B03/10000, gain * sample.B02/10000];";
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
Script(end+1) = "            'from': '"+TimeStringStart+"',";
Script(end+1) = "            'to': '"+TimeStringEnd+"'";
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

