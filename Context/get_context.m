function Output = get_context(LonPoints,LatPoints,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unified geophysical context loader.
%
%Can currently load:
%
%  1. 'LowResTopo' - 0.1 degree topography from easyTopo
%                  - requires easy_topo.mat, path can be set in this function
%  2. 'HighResTopo' - 30m topography from TessaDEM
%                   - requires TESSA data tiles. This function can generate an SCP script to download them.
%  3. 'SurfaceWind' - ERA5 1.5 degree resolution surface U and V
%                   - requires ERA5 netCDF data files as outputted by CDS API.
%
%planning to add:
%  A. tropopause height
%  B. stratopause height
%  C. wind speed at heights chosen by user
%  D. surface imagery
%  E. IMERG convection
%  F. sea surface temperature
%  G. Climate indices (ENSO, TSI, QBO, NAM, NAO).
%
%
%inputs:
%
%+++++ REQUIRED:
%
%     VarName            (type)      description
%     ------------------------------------------------------------------------------------------- 
%     LonPoints          (double)    Longitude points to return data for. 
%                                    Must have the same numnber of points as LatGrid.
%
%     LatPoints          (double)    Latitude  points to return data for. 
%                                    Must have the same numnber of points as LonGrid.
%
%
%++++TO REQUEST OUTPUTS:
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%     TimePoints             (double,         NaN)  required for output options marked with a *, in Matlab units
%
%     LowResTopo             (logical,      false)  return easyTopo 0.1 degree topography data (faster)
%     HighResTopo            (logical,      false)  return TessaDEM 30m topography data        (slower)
%  *  SurfaceWind            (logical,      false)  return surface wind from ERA5
%
%
%++++ADDITIONAL OPTIONS (prefix indicates associated output):
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%     LowResTopo_Path        (char, see in parser)  path to easytopo data file
%     HighResTopo_Path       (char, see in parser)  path to TessaDEM data files
%     HighResTopo_LRFill     (logical,      false)  fill gaps (poles and oceans) in Tessa data with easytopo data. Currently assumes some file paths, so turned off by default for safety.
%     HighResTopo_TileScript (logical       false)  generate SCP script to download required Tessa tiles
%     Era5_Path              (char, see in parser)  path to ERA5 data, used for SurfaceWinds
%
%By default the routine will return no useful data. Any chosen outputs
%must be switched on with flags
%
%Corwin Wright, c.wright@bath.ac.uk, 16/APR/2024
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

%optional flags
addParameter(p,'LowResTopo',  false,@islogical); %load easyTopo 0.1 degree topography
addParameter(p,'HighResTopo', false,@islogical); %load TessaDEM 30m topography
addParameter(p,'SurfaceWind', false,@islogical); %load surface winds from 1.5 degree ERA5

%variables used for many, but not all, datasets
addParameter(p,'TimePoints', NaN,@(x) validateattributes(x,{'numeric'},{'size',size(LonPoints)})); %time of each point, in Matlab units

%individual datasets
%%%%%%%%%%%%%%%%%%%%%

%used in several functions
addParameter(p,'Era5_Path', [LocalDataDir,'/ERA5/'],@ischar); %path to ERA5 data

%low-res topo
addParameter(p,'LowResTopo_Path', [LocalDataDir,'/topography/easy_tenth_degree_topography/','easy_topo.mat'],@ischar); %path to data

%high-res topo
addParameter(p,'HighResTopo_LRFill',    false,@islogical); %fill using low-res topography if needed. Currently this makes some assumptions about path, so you probably want it turned off.
addParameter(p,'HighResTopo_TileScript',false,@islogical); %return an SCP script to get the tiles needed for this option to work from eepc-0184
addParameter(p,'HighResTopo_Path',  [LocalDataDir,'/topography/tessaDEM/raw/'],@ischar); %path to data





%done - parse and restructure inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parse(p,LonPoints,LatPoints,varargin{:})
Settings = p.Results;
TimePoints = Settings.TimePoints; Settings = rmfield(Settings,'TimePoints');
clear p varargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialise output struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Output = struct();
Output.Lon = LonPoints;
Output.Lat = LatPoints;

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
                                   'TileScript', Settings.HighResTopo_TileScript);
  Output.HighResTopo = Alt;
  Output.TileScript  = TileScript;
  clear Alt TileScript HRTRes
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.SurfaceWind == true
  %check we fed in a time - this is required
  if sum(isnan(TimePoints)) == numel(TimePoints);
    warning('SurfaceWinds: no TimePoints provided. Skipping.')
  else
    %ok, create an interpolant and grab the surface wind (1000hPa)
    I = create_era5_interpolant(LonPoints,LatPoints,TimePoints,Settings,'SurfaceWinds');
    if ~strcmp(class(I),'double'); 
      Output.SurfaceU = I.U(LonPoints,LatPoints,TimePoints,ones(size(LonPoints)).*1000);
      Output.SurfaceV = I.V(LonPoints,LatPoints,TimePoints,ones(size(LonPoints)).*1000);    
    end
    clear I
    
  end;
end














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create ERA5 interpolant
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
