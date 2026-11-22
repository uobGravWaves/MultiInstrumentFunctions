%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D Savitzky-Golay filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Data = func_filter_sgolay(Data,Settings)

%create array to store background temperature
%put this in an if, in case we're calling multiple filters
if ~isfield(Data,'BG'); Data.BG = Data.Tp.*NaN; GetBG = 1; else GetBG = 0; end


%is the data all the same resolution?
%(fixed 2026/11/22: abs() added - for descending grids the ratio was negative, so irregular descending grids passed this test)
if abs(nanstd(flatten(diff(Data.Alt,1,2)))./nanmean(flatten(diff(Data.Alt,1,2)))) < 0.01;
  %if so, we can do this in a single pass
  FrameLen = abs(make_odd(round(Settings.SGLength./nanmean(flatten(diff(Data.Alt,1,2))))));
  BG = sgolayfilt(Data.Tp',Settings.SGOrder,FrameLen)';
  if GetBG == 1; Data.BG = BG; end
  Data.Tp  = Data.Temp-BG;
else

  %if not, we need a loop
  %(fixed 2026/11/22: FrameLen divided by mean altitude rather than mean level spacing, giving a frame of ~1 level,
  % so this branch crashed; the loop limit was also numel(Data.Alt,1), which is 1, so it would only have filtered the first profile)
  for iProf=1:1:size(Data.Alt,1)
    FrameLen = abs(make_odd(round(Settings.SGLength./nanmean(diff(Data.Alt(iProf,:))))));
    BG = sgolayfilt(Data.Tp(iProf,:),Settings.SGOrder,FrameLen);
    if GetBG == 1; Data.BG(iProf,:) = BG; end
    Data.Tp(iProf,:) = Data.Tp(iProf,:) - BG;
  end
end

if Settings.Verbose == 1; disp(['--> Savitzky-Golay filter applied vertically']); end

return
end



