function Data = apply_wg13(Data,Instrument,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use noise distributions comptued with make_wg2013_noisefloor() to
%filter out below-noise GWs produced using Alex08 analysis
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2025/01/19

%===================================
%OUTPUTS:
%   OutData: same struct as gwanalyse_limb(), with two changes:
%            1. below-noise waves have NaN in all fields
%            2. An additional field 'NoiseFloor' with the noise floor for this wave
%
%===================================
%
%INPUTS:
%
%  required:
%    Data [struct] - output struct from gwanalyse_limb()
%    Instrument    - name of instrument, used to select noise floor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instruments we can use this on, and their properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%HIRDLS
InstInfo.HIRDLS.NoiseFile  = [LocalDataDir,'/HIRDLS/noisefloor_HIRDLS.mat'];

%MLS
InstInfo.MLS.NoiseFile     = [LocalDataDir,'/MLS/noisefloor_MLS.mat'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

%data structure
addRequired(p,'Data',@isstruct); %input data must be a struct. Will be hand-parsed later.

%instrument(s)
CheckInst  = @(x) validateStringParameter(x,fieldnames(InstInfo),mfilename,Instrument);
addRequired(p,'Instrument',CheckInst)

%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,Data,Instrument,varargin{:})

%pull out the contents into struct "Settings", used throughout rest of routine
Settings = p.Results;
clearvars -except InstInfo Settings Instrument Data

%extract just the metadata for the instrument we want
InstInfo  = InstInfo.(Instrument);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the nosie floor and produce an interpolant from it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load noise
Noise = load(InstInfo.NoiseFile);

%put all the arrays in ascending order
Axes = {'DayScale','LonScale','LatScale','AltScale','LzScale'}; % order is important - must be same as definition of CutOff
for iAx=1:1:numel(Axes)
  [Noise.(Axes{iAx}),idx] = sort(Noise.(Axes{iAx}),'ascend');
  Noise.CutOff = index_dim(Noise.CutOff,idx,iAx);
end

%produce interpolant
N = griddedInterpolant({Noise.DayScale,Noise.LonScale,Noise.LatScale,Noise.AltScale,Noise.LzScale},Noise.CutOff);
clear Noise Axes iAx idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the value of the interpolant (i.e. nosie floor) at each point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dd = date2doy(Data.Time);
Floor = N(dd,Data.Lon,Data.Lat,Data.Alt,Data.Lz);
clear dd N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% apply mask, and return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BelowNoise = find(Data.A < Floor);
Fields = fieldnames(Data);
for iF=1:1:numel(Fields)
  f = Data.(Fields{iF});
  f(BelowNoise) = NaN;
  Data.(Fields{iF}) = f;
end; clear Fields iF f

Data.NoiseFloor = Floor;


return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to validate list of allowed instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validateStringParameter(varargin)
validatestring(varargin{:});
end



