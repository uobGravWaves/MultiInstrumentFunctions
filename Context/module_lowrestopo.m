function LowResTopo = module_lowrestopo(Settings,LonPoints,LatPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get_context() module to load 0.1 degree easyTopo data
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/05/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first check if the input data file exists
if ~exist(Settings.LowResTopo_Path,'file')
  warning('LowResTopo: easyTopo data not located, skipping.')
else

  %load the data, create an interpolant, and put it on output grid
  EasyTopo = load(Settings.LowResTopo_Path);
  I = griddedInterpolant(EasyTopo.topo.lats,EasyTopo.topo.lons,EasyTopo.topo.elev);
  LowResTopo = I(LatPoints,LonPoints);

  clear I EasyTopo
end

return
end
