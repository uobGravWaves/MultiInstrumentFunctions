function VarList = list_non_modal_size(Struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Find the modal field size in a struct, and return a list of all fields not this size
%inputs:
%  Struct   - the struct to operate on
%
%outputs:
%  VarList  - cell array lsiting all variables not the modal size

%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/01
%updated 2023/08/12 to let the user choose a dimension to operate along
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VarList = {};

%get largest number of dimensions for any field
f = fieldnames(Struct); MaxDims = 0;
for iF=1:1:numel(f); if ndims(Struct.(f{iF})) > MaxDims; MaxDims = ndims(Struct.(f{iF})); end; end
clear iF

%get size of each field
Sizes = NaN([numel(f),MaxDims]);
for iF=1:1:numel(f)
  sz = size(Struct.(f{iF}));
  Sizes(iF,1:numel(sz)) = sz;
end
clear iF sz MaxDims

%find most common size
[U,I,J] = unique(Sizes, 'rows');
Score = NaN(numel(U),1);
for iU=1:1:numel(U); Score(iU) = numel(find(J == iU)); end
[~,idx] = max(Score);
ModalSize = U(idx,:);
clear U I J Score iU idx

%finally, get a list of vars to ignore
for iF=1:1:numel(f); if ~isequal(size(Struct.(f{iF})),ModalSize); VarList{end+1} = f{iF}; end; end

clear iF f IgnoreWrongSize ModalSize Sizes

return
end