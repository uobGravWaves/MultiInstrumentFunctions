function [Data,FileList] = module_load_ALOMAR(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from ALOMAR Lidar for get_limbsounders()
%
%Corwin Wright, c.wright@bath.ac.uk, 14/July/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%where is Alomar?
Site.Lat = 69.38;
Site.Lon = 16.08;

%what resolution do we want?
Res = 'T1'; %this is the higher resolution product. Expose a switch for this later
% Res = 'T2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%list of variables
Vars = {[Res,'_lidar'],[Res,'_lidar_unc'],'time_lidar','z_lidar'};

%files are yearly, so we need to identify all distinct years and then pull out what we need
[y,~,~] = datevec(minmax(Settings.TimeRange));
k = 0; FileList = {};

for Y=min(y):1:max(y);

  %identify annual file and check it exists
  FileName = [InstInfo.Path,num2str(Y),'_V2025-07-09_file4juli.h5'];
  if ~exist(FileName,'file'); continue; end

  %get a list of days actually present in the file
  DayList = h5info(FileName);
  DayList = {DayList.Groups.Name};
  
  %file name
  k = k+1;
  FileList{k} = FileName;


  %load each day
  for iDay=1:1:numel(DayList)

    %get the data
    a = struct();
    for iVar=1:1:numel(Vars)
      a.(Vars{iVar}) = h5read(FileName, "/"+DayList{iDay}+"/"+Vars{iVar});
    end; clear iVar

    %convert time units and tranpose axis
    a.time_lidar = datenum(datetime(a.time_lidar, 'ConvertFrom', 'posixtime'))';

    %profile number and file number
    a.PN = 1:1:numel(a.time_lidar);
    a.FN = ones(size(a.time_lidar)).*k;

    %store data
    if ~exist('Store','var'); Store = a;
    else
      Store = cat_struct(Store,a,2,{'z_lidar'});
    end
    clear a

  end; clear iDay DayList

end; clear Y FileName y k


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subset in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxA = find(Store.time_lidar >= Settings.TimeRange(1));
idxB = find(Store.time_lidar <= Settings.TimeRange(2));
idx = intersect(idxA,idxB);
Store = reduce_struct(Store,idx,{'z_lidar'},2);
clear idxA idxB idx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rename primary variables and create meta variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data.Temp       = Store.([Res,'_lidar'])';
Data.Lat        = ones(size(Data.Temp)).*Site.Lat;
Data.Lon        = ones(size(Data.Temp)).*Site.Lon;
Data.Time       = repmat(Store.time_lidar',1,size(Data.Temp,2));
Data.Alt        = repmat(Store.z_lidar',size(Data.Temp,1),1);
Data.Pres       = h2p(Data.Alt);
Data.SourceFile = repmat(Store.FN',1,size(Data.Temp,2));
Data.SourceProf = repmat(Store.PN',1,size(Data.Temp,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return
end
