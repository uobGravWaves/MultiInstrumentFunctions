function Store = module_conform_vars(Store,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%put user-requested AdditionalVars onto the [profiles x levels] grid
%used by the modules of get_limbsounders(), once per file/day, BEFORE the day is
%concatenated onto the main data (after concatenation, orientation can no longer be fixed)
%
%The grid is taken from Store.Lat, which every module provides as [NProf x NLev].
%Each AdditionalVar is handled as follows:
%  [NProf x NLev (x ...)]     per-point                       -> unchanged
%  [NLev x NProf (x ...)]     per-point, file orientation     -> dims 1 and 2 swapped
%  vector of length NProf     per-profile (e.g. orbit number) -> repeated over levels
%  vector of length NLev      per-level                       -> repeated over profiles
%  table with NProf rows      per-profile metadata            -> unchanged (handled as tables downstream)
%  anything else              (e.g. scalars, other sizes)     -> error, unless Settings.IgnoreNonModal is
%                                                                true, in which case it is passed through
%If NProf == NLev, orientation cannot be determined: a matching [N x N] field is accepted with
%a warning, and a vector of length N is an error (per-profile or per-level is ambiguous).
%
%inputs:
%  Store    - struct of one file/day of data, with Store.Lat [NProf x NLev]
%  Settings - get_limbsounders settings (uses AdditionalVars, IgnoreNonModal, Instrument)
%outputs:
%  Store    - as input, with AdditionalVars conformed
%
%Corwin Wright, c.wright@bath.ac.uk, and Claude, 2026/11/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NProf = size(Store.Lat,1); NLev = size(Store.Lat,2);

for iVar=1:1:numel(Settings.AdditionalVars)
  Name = Settings.AdditionalVars{iVar};
  if ~isfield(Store,Name); error('get_limbsounders:VarNotFound','AdditionalVar ''%s'' was not found in the %s data',Name,Settings.Instrument); end
  F = Store.(Name);
  sz = size(F);

  if istable(F)
    if height(F) == NProf; continue; end
  elseif isnumeric(F) || islogical(F)
    if isequal(sz(1:2),[NProf,NLev])                         %per-point
      if NProf == NLev && NProf > 1; warning('get_limbsounders:SquareGrid','AdditionalVar ''%s'': file has as many profiles as levels (%d), so its orientation cannot be checked; assumed [profiles x levels]',Name,NProf); end
      continue
    elseif isequal(sz(1:2),[NLev,NProf]) && NProf ~= NLev    %per-point, file orientation
      Store.(Name) = permute(F,[2,1,3:ndims(F)]); continue
    elseif isvector(F) && numel(F) == NProf && NProf ~= NLev %per-profile
      Store.(Name) = repmat(F(:),1,NLev); continue
    elseif isvector(F) && numel(F) == NLev && NProf ~= NLev  %per-level
      Store.(Name) = repmat(F(:)',NProf,1); continue
    end
  end

  %anything else: metadata or unrecognised shape
  if Settings.IgnoreNonModal == true; continue; end
  error('get_limbsounders:VarShape',['AdditionalVar ''%s'' from %s has size %s (%s), which is not per-point [%d x %d], per-profile (%d) or per-level (%d). ' ...
        'Set ''IgnoreNonModal'',true to pass it through without processing.'],Name,Settings.Instrument,mat2str(sz),class(F),NProf,NLev,NProf,NLev);
end

return
end
