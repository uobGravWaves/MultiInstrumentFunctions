function [Data,FileList] = module_load_pwdata(Settings,InstInfo,Vars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from prepared planetary-wave filtered data
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

  %I don't know why this if-loop is needed, but it works
  if numel(Data.Time) == 0
    Data.Temp_PW = Store.Temp_PW;
    Data = cat_struct(Data,Store,1,{'Temp_PW'});
  else
    Data = cat_struct(Data,Store,1);
  end
  

  clear Store iVar

end; clear DayNumber
if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

end
