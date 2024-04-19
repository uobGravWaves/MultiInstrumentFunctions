function Output = get_context(LonPoints,LatPoints,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unified geophysical context loader.
%
%Can currently load:
%
%  1. 'LowResTopo'  - 0.1 degree topography from easyTopo
%                   - requires easy_topo.mat, path can be set in this function
%  2. 'HighResTopo' - 30m topography from TessaDEM
%                   - requires TESSA data tiles. This function can generate an SCP script to download them.
%  3. 'Wind'        - ERA5 1.5 degree resolution U and V at chosen pressure levels (can also do 
%                   - requires ERA5 netCDF data files as outputted by CDS API.
%                   - output will have an extra dimension corresponding to the pressure levels requested
%  4. 'Indices'     - climate indices: 'QBO','ENSO','JetFuelPrice',,'NAM','NAO','TSI','SeaIce','AMO'
%                   - output calculated based on time only, lat and lon must be set but will be ignored
%
%
%planning to add:
%  A. tropopause height
%  B. stratopause height
%  C. surface imagery
%  D. IMERG convection
%  E. sea surface temperature
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
%  *  TimePoints             (double,         NaN)  Same size as LonPoints.    Required for output options marked with a *, in Matlab units
%  ^  Pressure               (double,         NaN)  1D array of levels in hPa. Required for output options marked with a ^, in Matlab units
%
%     Everything             (logical,      false)  try and return all the below. Overrides all individual choices.
%     LowResTopo             (logical,      false)  return easyTopo 0.1 degree topography data (faster)
%     HighResTopo            (logical,      false)  return TessaDEM 30m topography data        (slower)
%  *^ Wind                   (logical,      false)  return winds from ERA5
%  *  Indices                (logical,      false)  return climate indices
%
%
%++++SUPPORT OPTIONS FOR SPECIFIC OUTPUTS (prefix indicates associated output):
%
%     VarName                (type,       default)  description
%     -------------------------------------------------------------------------------------------
%     LowResTopo_Path        (char, see in parser)  path to easytopo data file
%     HighResTopo_Path       (char, see in parser)  path to TessaDEM data files
%     HighResTopo_LRFill     (logical,      false)  fill gaps (poles and oceans) in Tessa data with easytopo data. Currently assumes some file paths, so turned off by default for safety.
%     HighResTopo_TileScript (logical       false)  generate SCP script to download required Tessa tiles
%     Indices_Path           (char, see in parser)  path to directory containing climate index data
%     Era5_Path              (char, see in parser)  path to ERA5 data, used for Winds
%
%By default the routine will return no useful data. Any chosen outputs
%must be switched on with flags.
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
addParameter(p,'Everything',  false,@islogical); %try to load all the below
addParameter(p,'LowResTopo',  false,@islogical); %load easyTopo 0.1 degree topography
addParameter(p,'HighResTopo', false,@islogical); %load TessaDEM 30m topography
addParameter(p,'Wind',        false,@islogical); %load winds from 1.5 degree ERA5
addParameter(p,'Indices',     false,@islogical); %load climate indices

%variables used for many, but not all, datasets
addParameter(p,'TimePoints', NaN,@(x) validateattributes(x,{'numeric'},{'size',size(LonPoints)})); %time of each point, in Matlab units
addParameter(p,'Pressure',   NaN,@(x) validateattributes(x,{'numeric'},{'<=',1200}));              %pressure levels for output, in hPa


%individual datasets
%%%%%%%%%%%%%%%%%%%%%

%paths
addParameter(p,'Era5_Path', [LocalDataDir,'/ERA5/'],@ischar); %path to ERA5 data
addParameter(p,'LowResTopo_Path', [LocalDataDir,'/topography/easy_tenth_degree_topography/','easy_topo.mat'],@ischar); %path to data
addParameter(p,'HighResTopo_Path',  [LocalDataDir,'/topography/tessaDEM/raw/'],@ischar); %path to data
addParameter(p,'Indices_Path', [LocalDataDir,'/Miscellany/'],@ischar); %path to cliamte index data

%other options
addParameter(p,'HighResTopo_LRFill',    false,@islogical); %fill high-res topo using using low-res topography if needed. Currently this makes some assumptions about path, so you probably want it turned off.
addParameter(p,'HighResTopo_TileScript',false,@islogical); %return an SCP script to get the tiles needed for the high-res topo option from eepc-0184



%done - parse and restructure inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parse(p,LonPoints,LatPoints,varargin{:})
Settings = p.Results;
TimePoints = Settings.TimePoints; Settings = rmfield(Settings,'TimePoints');
clear p varargin

%override actual options if 'EveryThing' is set
if Settings.Everything == true
  Settings.LowResTopo  = true;
  Settings.HighResTopo = true;
  Settings.Wind        = true;
  Settings.Indices     = true;
  warning('"Everything" option set - all output options will be attempted')
end

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
    if ~strcmp(class(I),'double'); 

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
            Output.(Indices{iIndex}) = interp1(QBO.Time,QBO.QBO,TimePoints);
            clear QBO
          case 'ENSO'
            ENSO = load([Root,'/nino34.mat']);
            Output.(Indices{iIndex}) = interp1(ENSO.Time,ENSO.Nino34,TimePoints);
            clear ENSO
          case 'Fuel'
            Fuel = load([Root,'/jet_fuel_price.mat']);
            Output.(Indices{iIndex}) = interp1(Fuel.Time,Fuel.Price,TimePoints);
            clear Fuel
          case 'NAM'
            NAM = load([Root,'/daily_nam.mat']);
            Output.(Indices{iIndex}) = interp1(NAM.Time,NAM.NAM,TimePoints);
            clear NAM
          case 'NAO'
            NAO = load([Root,'/nao.mat']);
            Output.(Indices{iIndex}) = interp1(NAO.Time,NAO.NAO,TimePoints);
            clear NAO
          case 'SSTs'
            SSTs = load([Root,'/ssts.mat']);
            Output.(Indices{iIndex}) = interp1(SSTs.Time,SSTs.SSTs,TTimePoints);
            clear SSTs
          case 'TSI'
            TSI = load([Root,'/tsi.mat']);
            Output.(Indices{iIndex}) = interp1(TSI.Time,TSI.TSI,TimePoints);
            clear TSI
          case 'Time'
            Output.(Indices{iIndex}) = TimeScale;
          case 'SeaIce'
            SeaIce = readmatrix([Root,'/N_seaice_extent_daily_v3.0.csv']);
            t = datenum(SeaIce(:,1),SeaIce(:,2),SeaIce(:,3));
            Output.(Indices{iIndex}) = interp1(t,SeaIce(:,4),TimePoints);
            clear SeaIce t
          case 'AMO'
            AMO = load([Root,'/AMO.mat']);
            Output.(Indices{iIndex}) = interp1(AMO.Time,AMO.AMO,TimePoints);
            clear AMO
        end
      catch; warning(['Indices: error locating input data for ',Indices{iIndex},'; skipping.'])
      end
    end
  end
  clear Root Indices iIndex
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
