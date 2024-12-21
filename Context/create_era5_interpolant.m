
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

  %there are some minor format differences between older and newer files. If we have a newer file,
  %reformat to be like the older so the existing code works
  if isfield(E5,'model_level'); E5.level = E5.model_level; end;
  if isfield(E5,'valid_time');  E5.time  = E5.valid_time;  end;

  %some files don't contain all three basic variables. If so, create dummy variables for the logic.
  if ~isfield(E5,'u'); E5.u = NaN(numel(E5.time),numel(E5.level),numel(E5.latitude),numel(E5.longitude)); end
  if ~isfield(E5,'v'); E5.v = NaN(numel(E5.time),numel(E5.level),numel(E5.latitude),numel(E5.longitude)); end
  if ~isfield(E5,'t'); E5.t = NaN(numel(E5.time),numel(E5.level),numel(E5.latitude),numel(E5.longitude)); end

  %if we have less than 137 levels, put it onto 137 lines to preserve logic below
  if numel(E5.level) ~= 137;
    sz = size(E5.u);
    New = spawn_uniform_struct({'u','v','t'},[sz(1),137,sz(3),sz(4)]);
    New.u(:,E5.level,:,:) = E5.u;
    New.v(:,E5.level,:,:) = E5.v;
    New.t(:,E5.level,:,:) = E5.t;    

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
end


