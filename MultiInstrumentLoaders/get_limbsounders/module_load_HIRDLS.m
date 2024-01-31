function [Data,FileList] = module_load_HIRDLS(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from HIRDLS for get_limbsounders()
%
%Corwin Wright, c.wright@bath.ac.uk, 05/NOV/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%empty structure to store data
Data = struct();
for iVar=1:1:numel(Vars); Data.(Vars{iVar}) = []; end
FileCount = 0;
FileList = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for DayNumber=floor(min(Settings.TimeRange)):1:floor(max(Settings.TimeRange));

  %work out year and day number and hence filepath
  File = wildcardsearch(InstInfo.Path,['*_',yyyyDDD(DayNumber),'*']);
  if numel(File) == 0; clear y dn File; continue; end

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File{1},filesep); FileList{end+1} = f{end}; clear f;


  %load variables we need
  Store = struct();
  for iVar=1:1:numel(Vars)
    switch Vars{iVar}
      case 'Temp';       Store.Temp         = get_HIRDLS(File{1},'Temperature')';
      case 'Lat';        Store.Lat          = get_HIRDLS(File{1},'Latitude');
      case 'Lon';        Store.Lon          = get_HIRDLS(File{1},'Longitude');
      case 'Alt';        Store.Alt          = get_HIRDLS(File{1},'Altitude')'./1000;
      case 'Pres';       Store.Pres         = get_HIRDLS(File{1},'Pressure');
      case 'Time';       Store.Time         = datenum(1993,1,1,0,0,get_HIRDLS(File{1},'Time'));
      case 'SourceProf'; Store.SourceProf   = single((1:1:numel(Store.Lat))');
      case 'SourceFile'; Store.SourceFile   = ones(size(Store.SourceProf)).*FileCount;
      otherwise;
        try;   Store.(Vars{iVar}) = get_HIRDLS(File{1},Vars{iVar});
        catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
        end
    end
  end

  %reshape 1d variables
  Store.Lat        = repmat(Store.Lat,         [1,size(Store.Temp,2)]);
  Store.Lon        = repmat(Store.Lon,         [1,size(Store.Temp,2)]);
  Store.Time       = repmat(Store.Time,        [1,size(Store.Temp,2)]);
  Store.Pres       = repmat(Store.Pres',       [size(Store.Temp,1),1]);
  Store.SourceProf = repmat(Store.SourceProf,  [1,size(Store.Temp,2)]);
  Store.SourceFile = repmat(Store.SourceFile,  [1,size(Store.Temp,2)]);

  %HIRDLS stores geolocation at the 30km level, but travels while scanning up and down
  %see Wright et al (ACP, 2015) for the logic of what we're going to do here to reverse
  %this choice and get 'true' lat and lon values

  %find the 30km level
  [~,zidx] = min(abs(nanmean(Store.Alt,1)-30));

  %for each profile compute the directional azimuth
  Theta = azimuth(Store.Lat(1:end-1,zidx),Store.Lon(1:end-1,zidx), ...
    Store.Lat(2:end,  zidx),Store.Lon(2:end,  zidx),'degrees');
  Theta(end+1) = Theta(end); %close enough for last point

  %hence find the true(ish) location of each point in 2D space
  KmAlongtrackPerKmVertical = 0.6;
  dx = Store.Lat.*NaN;
  for iLev=1:1:size(Store.Lat,2)
    dx(:,iLev) =  KmAlongtrackPerKmVertical.*(iLev-zidx);
    [Store.Lat(:,iLev),Store.Lon(:,iLev)] = reckon(Store.Lat(:,zidx),Store.Lon(:,zidx),km2deg(dx(:,iLev)),Theta);
  end
  clear zidx Theta KmAlongtrackPerKmVertical dx iLev

  %store in main repository
  Data = cat_struct(Data,Store,1);
  clear Store iVar y dn File


end; clear DayNumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return





function Output = get_HIRDLS(FilePath,Variable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract a single variable from a HIRDLS HDF5 data file
%returns NaN if variable does not exist
%
%adapted from a previous version developed by N Hindley, U Bath, May 2014
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%11/MAY/2014
%
%inputs
%---------
%
%FilePath: path to file
%Variable: variable wanted from file
%
%outputs
%---------
%
%Output: variable contents (or NaN if failed)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------
%list possible variables in an HDF5-EOS-HIRDLS file, by type
%------------------------------------------------------------------------------------

GeolocationFields = {'Altitude';'Latitude';'LocalSolarTime';'Longitude'; ...
                     'OrbitAscendingFlag';'OrbitNumber';'Pressure';      ...
                     'ProfileID';'ScanAzimuthAtNominalAltitude';         ...
                     'ScanElevationAtNominalAltitude';'ScanTable';       ...
                     'ScanUpFlag';'ScienceScanMode';'SecondsInDay';      ...
                     'SolarZenithAngle';'SpacecraftAltitude';            ...
                     'SpacecraftLatitude';'SpacecraftLongitude';         ...
                     'TangentHeightAtNominalAltitude';'Time';            ...
                     'ViewDirectionAtNominalAltitude'};

DataFields = {'10.8MicronCloudAerosolFlag';'10.8MicronExtinction';             ...
              '10.8MicronExtinctionNormChiSq';'10.8MicronExtinctionPrecision'; ...
              '10.8MicronExtinctionQuality';'12.1MicronCloudAerosolFlag';      ...
              '12.1MicronExtinction';'12.1MicronExtinctionNormChiSq';          ...
              '12.1MicronExtinctionPrecision';'12.1MicronExtinctionQuality';   ...
              '17.4MicronCloudAerosolFlag';'17.4MicronExtinction';             ...
              '17.4MicronExtinctionNormChiSq';'17.4MicronExtinctionPrecision'; ...
              '17.4MicronExtinctionQuality';'7.1MicronCloudAerosolFlag';       ...
              '7.1MicronExtinction';'7.1MicronExtinctionNormChiSq';            ...
              '7.1MicronExtinctionPrecision';'7.1MicronExtinctionQuality';     ...
              '8.3MicronCloudAerosolFlag';'8.3MicronExtinction';               ...
              '8.3MicronExtinctionNormChiSq';'8.3MicronExtinctionPrecision';   ...
              '8.3MicronExtinctionQuality';'CFC11';'CFC11NormChiSq';           ...
              'CFC11Precision';'CFC11Quality';'CFC12';'CFC12NormChiSq';        ...
              'CFC12Precision';'CFC12Quality';'CH4';'CH4NormChiSq';            ...
              'CH4Precision';'CH4Quality';'ClONO2';'ClONO2NormChiSq';          ...
              'ClONO2Precision';'ClONO2Quality';'CloudTopPressure';'GPH';      ...
              'GPHPrecision';'H2O';'H2ONormChiSq';'H2OPrecision';'H2OQuality'; ...
              'HNO3';'HNO3NormChiSq';'HNO3Precision';'HNO3Quality';'N2O';      ...
              'N2O5';'N2O5NormChiSq';'N2O5Precision';'N2O5Quality';            ...
              'N2ONormChiSq';'N2OPrecision';'N2OQuality';'NO2';'NO2NormChiSq'; ...
              'NO2Precision';'NO2Quality';'O3';'O3NormChiSq';'O3Precision';    ...
              'O3Quality';'RawGPH';'RawGPHPrecision';'Temperature';            ...
              'TemperatureNormChiSq';'TemperaturePrecision';'TemperatureQuality'};

%------------------------------------------------------------------------------------
%identify type of desired variable
%------------------------------------------------------------------------------------
   
IsGeo  = numel(find(ismember(GeolocationFields,Variable)));
IsData = numel(find(ismember(DataFields,Variable       )));

%------------------------------------------------------------------------------------
%load data
%------------------------------------------------------------------------------------

%open hdf5 file
FileID = H5F.open (FilePath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');


if IsGeo == 1; 
  FieldPath = strcat('HDFEOS/SWATHS/HIRDLS/Geolocation Fields/',Variable);
elseif IsData ==1; 
  FieldPath = strcat('HDFEOS/SWATHS/HIRDLS/Data Fields/',Variable);
else
  disp(['Error: ''',Variable,''' is not present in this file format']); 
  Output = NaN;
  return
end;

%extract data and close file
data_id = H5D.open (FileID,FieldPath);
Output = H5D.read(data_id,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT');
H5F.close(FileID)


return,Output









