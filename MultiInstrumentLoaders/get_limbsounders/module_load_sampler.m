function [Data,FileList] = module_load_ACE_GNSS(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data produced by the Wright and Hindley (2018) sampler
%
%Corwin Wright, c.wright@bath.ac.uk, 16/AUG/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%empty structure to store data
Data = struct();
for iVar=1:1:numel(Vars); Data.(Vars{iVar}) = []; end
FileCount = 0;
FileList = {};

%check if the requested path exists
if ~exist(Settings.SamplePath,'dir');
  error([Settings.SamplePath,' does not exist'])
  return
end
InstInfo.Path = Settings.SamplePath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Verbose == 1;
  textprogressbar(['----> Loading ',Settings.Instrument,' data (from ',Settings.SamplePath,') ']); 
  k = 0; n = ceil(range(Settings.TimeRange));
end

for DayNumber=floor(min(Settings.TimeRange)):1:floor(max(Settings.TimeRange));

  if Settings.Verbose == 1; 
    k = k+1;
    textprogressbar(k./n.*100);
  end

  %work out year and day number and hence filepath
  File = [InstInfo.Path,'/sampled_',num2str(DayNumber),'.mat'];
  if ~exist(File,'file'); continue; end
 
  %load the data, and rename to match expected nomenclature
  InstData = load(File); InstData = InstData.Sampled_Data;
  InstData.SourceFile = File;
  InstData.Temp = InstData.T;   InstData = rmfield(InstData,'T');
  InstData.Pres = InstData.Prs; InstData = rmfield(InstData,'Prs');
  
  %store file information
  FileCount = FileCount+1;
  f = strsplit(File,filesep); FileList{end+1} = f{end}; clear f;
  clear y dn File

  %get the variables we want
  Store = struct();
  for iVar=1:1:numel(Vars)

    switch Vars{iVar}
      case 'SourceFile';Store.SourceFile = ones(size(Store.SourceProf)).*FileCount;
      otherwise;   Store.(Vars{iVar}) = InstData.(Vars{iVar})';
    end

  end; clear iVar

  %store in main repository
  Data = cat_struct(Data,Store,1);

  clear Store iVar

end; clear DayNumber

if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end
