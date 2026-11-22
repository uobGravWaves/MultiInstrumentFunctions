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
%updated 2026/11/22 [Claude] fixed wrong result when any field has more dimensions than the others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VarList = {};

%get largest number of dimensions for any field
f = fieldnames(Struct); MaxDims = 0;
for iF=1:1:numel(f); if ndims(Struct.(f{iF})) > MaxDims; MaxDims = ndims(Struct.(f{iF})); end; end
clear iF

%get size of each field
Sizes = ones([numel(f),MaxDims]); %missing trailing dimensions are size 1, as in MATLAB (was NaN, which made unique() treat every such row as distinct)
for iF=1:1:numel(f)
  sz = size(Struct.(f{iF}));
  Sizes(iF,1:numel(sz)) = sz;
end
clear iF sz MaxDims

%find most common size
[U,I,J] = unique(Sizes, 'rows');
Score = NaN(size(U,1),1);
for iU=1:1:size(U,1); Score(iU) = numel(find(J == iU)); end
[~,idx] = max(Score);
ModalSize = U(idx,:);
clear U I J Score iU idx

%finally, get a list of vars to ignore
for iF=1:1:numel(f); if ~isequal(Sizes(iF,:),ModalSize); VarList{end+1} = f{iF}; end; end

clear iF f IgnoreWrongSize ModalSize Sizes

return
end