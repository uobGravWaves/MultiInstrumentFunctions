function [Data,FileList] = module_load_pwdata(Settings,InstInfo,Vars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from prepared planetary-wave filtered data
%
%Corwin Wright, c.wright@bath.ac.uk, 05/NOV/2023
%
%modified 2025/03/08 to cheat a little on load times by concatenating
%to an intermediate struct, as Juli was trying to load whole years of
%data at once and the programme was getting exponential slower...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%empty structure to store data
Data = struct();
for iVar=1:1:numel(Vars); Data.(Vars{iVar}) = []; end
FileCount = 0;
FileList = {};

%counters for intermediate data stucture, used to speed load times
InnerDayCount = 0; InnerLimit = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Verbose == 1;
  textprogressbar('----> Loading Hindley23-style PW data '); 
  k = 0; n = ceil(range(Settings.TimeRange));
end

for DayNumber=floor(min(Settings.TimeRange)):1:floor(max(Settings.TimeRange));
  if Settings.Verbose == 1; 
    k = k+1;
    textprogressbar(k./n.*100);
  end

  %work out year and day number and hence filepath
  [y,~,~] = datevec(DayNumber); dn = date2doy(DayNumber);
  if ~strcmpi(InstInfo.Inst,'Misc');
    File = [InstInfo.Path,'/',sprintf('%04d',y),filesep,InstInfo.Inst,'_PWs_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.mat'];
  else
    File = wildcardsearch(Settings.MiscInfo{1},['*',Settings.MiscInfo{2},'*',sprintf('%04d',y),'d',sprintf('%03d',dn),'*']);
    if numel(File) == 0; clear y dn File; continue; end
    File = File{1};
  end


  if ~exist(File,'file'); clear y dn File; continue; end

  %load the data
  InstData = load(File);

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File,filesep); FileList{end+1} = f{end}; clear f;
  clear y dn File


  %get the variables we want
  Store = struct();
  for iVar=1:1:numel(Vars);
    if strcmpi(Vars{iVar},      'Pres') == 1 && ~isfield(InstData,'Pres'); 
      InstData.Pres = h2p(InstData.Alt);
    end
    if strcmpi(Vars{iVar},'SourceProf') == 1; 
      InstData.SourceProf = repmat(1:1:size(InstData.Lat,1),size(InstData.Lat,2),1)'; 
    end
    if strcmpi(Vars{iVar},'SourceFile') == 1; 
      InstData.SourceFile = repmat(FileCount,size(InstData.Lat)); 
    end

    Store.(Vars{iVar}) = InstData.(Vars{iVar}); 
  end; clear iVar

  %store in main repository
  %%%%%%%%%%%%%%%%%%%%%%%%%

  %to reduce cat() delays for large datasets, we're going to cat onto a small array at first
  % then write to the big one every few days. This means we're not concatening to the really
  %large array every single loop. The program will still get slower as more data gets added,
  %but this at least slows the rate. It's still a bad idea to do very long time periods though.

  %create inner data struct
  if InnerDayCount == 0;
    InnerData = struct();
    for iVar=1:1:numel(Vars); InnerData.(Vars{iVar}) = []; end
  end
  InnerDayCount = InnerDayCount+1;

  %store the data
  if numel(InnerData.Time) == 0; InnerData.Temp_PW = Store.Temp_PW; InnerData = cat_struct(InnerData,Store,1,{'Temp_PW'});
  else                           InnerData = cat_struct(InnerData,Store,1);
  end

  %if we've reached the limit or finished the last day, copy over to the main store
  if InnerDayCount == InnerLimit | DayNumber == floor(max(Settings.TimeRange))

    if numel(Data.Time) == 0; Data.Temp_PW = InnerData.Temp_PW; Data = cat_struct(Data,InnerData,1,{'Temp_PW'});
    else                      Data = cat_struct(Data,InnerData,1);
    end

    clear InnerData
    InnerDayCount = 0;
  end

  clear Store iVar


end; clear DayNumber
if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!'); end


%catch any InnerData contents lost due to premature exit from the primary loop
if exist('InnerData','var')
  if numel(Data.Time) == 0; Data.Temp_PW = InnerData.Temp_PW; Data = cat_struct(Data,InnerData,1,{'Temp_PW'});
  else                      Data = cat_struct(Data,InnerData,1);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

end