
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

  %if we have less than 137 levels, put it onto 137 lines to preserve logic below
  if numel(E5.level) ~= 137;
    sz = size(E5.u);
    New = spawn_uniform_struct({'u','v','t'},[sz(1),137,sz(3),sz(4)]);
    New.u(:,E5.level,:,:) = E5.u;
    New.v(:,E5.level,:,:) = E5.u;
    New.t(:,E5.level,:,:) = E5.u;    

    E5.u = New.u;
    E5.t = New.t;
    E5.v = New.v;
    clear sz New
  end

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
