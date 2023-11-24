function [Data,FileList] = module_load_sampleddata(Settings,InstInfo,Vars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from psampler output
%
%Corwin Wright, c.wright@bath.ac.uk, 24/NOV/2023
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


  %find file
  Path = [InstInfo.Path,'/',InstInfo.Sampled.InstName,'/',InstInfo.Sampled.ModelName,'/'];
  Name = ['sampled_',num2str(DayNumber),'*',num2str(InstInfo.Sampled.ModelSubSet),'.mat'];
  File = wildcardsearch(Path,Name);
  if numel(File) == 0; clear Path File; continue; end
  File = File{1};

  %load the data
  InstData = load(File); InstData = InstData.Sampled_Data;

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File,filesep); FileList{end+1} = f{end}; clear f;
  clear y dn File


  %get the variables we want
  Store = struct();
  for iVar=1:1:numel(Vars);
    if strcmpi(Vars{iVar},      'Pres') == 1; InstData.Pres = InstData.Prs; end
    if strcmpi(Vars{iVar},      'Temp') == 1; InstData.Temp = InstData.T; end
    if strcmpi(Vars{iVar},'SourceProf') == 1; InstData.SourceProf = repmat(1:1:size(InstData.Lat,1),size(InstData.Lat,2),1)'; end
    if strcmpi(Vars{iVar},'SourceFile') == 1; InstData.SourceFile = repmat(FileCount,size(InstData.Lat)); end
    Store.(Vars{iVar}) = InstData.(Vars{iVar}); 
  end; clear iVar




  %store in main repository
  %%%%%%%%%%%%%%%%%%%%%%%%%

  if numel(Data.Time) == 0; 
    [Data,IgnoredVars] = cat_struct(Data,Store,1,{},'IgnoreWrongSize',true); 
    if numel(IgnoredVars) >0;  for iVar=1:1:numel(IgnoredVars); Data.(IgnoredVars{iVar}) = Store.(IgnoredVars{iVar}); end; end
  else; 
    Data = cat_struct(Data,Store,1,{},'IgnoreWrongSize',true);
  end



  clear Store iVar

end; clear DayNumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return