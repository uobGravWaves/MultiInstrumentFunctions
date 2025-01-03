function Struct = reduce_struct(Struct,SubSetIndices,VarsToExclude,Dim,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%For a series of identical fields in a struct(), select a given set of 
%indices from every field
%
%For backwards-compatability with older versions, VarsToExclude argument needs to be
%ordered before Dim. If we want to operate on ALL fields, just set this to [].
%
%inputs:
% required:
%  Struct        - the struct to operate on
%  SubSetIndices - the list of indices to select from each field
%  VarsToExclude - fields to ignore when subsettings. Can be set to just {}.
%  Dim           - the dimension to operate on for each field. If 0, will apply to whole dataset
% optional:
%  'IgnoreWrongSize' if set to true will skip any variables that are not the MOST COMMON size. Default false.
%
%outputs:
%  Struct        - the  reduced structure

%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/01
%updated 2023/08/12 to let the user choose a dimension to operate along
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%handle inputs
if ~exist('VarsToExclude','var'); VarsToExclude = {' '}; end  %assume applies to all vars
if ~exist(          'Dim','var'); Dim = 1;               end  %assume operation is along first dimension
if ~exist('varargin','var'); varargin = {}; end
for iV=1:2:numel(varargin);
  if strcmpi(varargin{iV},'IgnoreWrongSize'); IgnoreWrongSize = varargin{iV+1}; end
end; clear iV
if ~exist('IgnoreWrongSize','var'); IgnoreWrongSize = false; end

%if IgnoreWrongSize is true, add wrong-size variables to the ignore list
if IgnoreWrongSize == true; VarsToExclude = [VarsToExclude,list_non_modal_size(Struct)];end


Fields = fieldnames(Struct);
for iField=1:1:numel(Fields);
  
  %skip specified variables
  if any(strcmp(Fields{iField},VarsToExclude)); continue; end

  %reduce desired variables
  F = Struct.(Fields{iField});
  if strcmpi(class(F),'table')
    %tables can only be trimmed in the first two dimensions
    if Dim == 1;
      F = F(SubSetIndices,:);
    elseif Dim == 2;
      F = F(:,SubSetIndices);
    else
      warning('reduce_struct is trying to subset a table in a dimension higher than 2, skipping')
    end
  else
    %anything else, treat as normal
    if Dim == 0; F = index_dim(F(:),SubSetIndices,1);
    else;        F = index_dim(F,SubSetIndices,Dim);
    end
  end
  Struct.(Fields{iField}) = F;
  
  %done!
  
end

return
end