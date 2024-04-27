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
%  4. 'Indices'     - climate indices: 'QBO','ENSO','JetFuelPrice','NAM','NAO','TSI','SeaIce','AMO'
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