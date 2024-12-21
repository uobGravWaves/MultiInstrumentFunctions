function [HighResTopo,TileScript] = module_highrestopo(Settings,LonPoints,LatPoints)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get_context() module to load 30m TessaDEM data
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/05/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Alt,~,~,TileScript] = map_tessa(LonPoints,LatPoints, ...
                                 'ETFill',     Settings.HighResTopo_LRFill,     ...
                                 'DataDir',    Settings.HighResTopo_Path,       ...
                                 'TileScript', Settings.HighResTopo_TileScript, ...
                                 'ETPath',     Settings.LowResTopo_Path);
HighResTopo = Alt;

return
end