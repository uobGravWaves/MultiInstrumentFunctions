clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce noise floor files to apply WG2013 method to
%get_limbsounders output with a Hindley23 filter
%
%Corwin Wright, c.wright@bath.ac.uk, 2025/Jan/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%instrument selection
Settings.Instrument = 'HIRDLS';
Settings.OutFile = [LocalDataDir,'/corwin/noisefloor_',Settings.Instrument,'.mat'];

%processing choices
Settings.HeightScaleIn = 15:1:80;

%statistical choices
Settings.Percentile    = 99;    %what percentile of the distribution is the stored cutoff?
Settings.NProfilePairs = 2500;  %how many profile-pairs should we use in each box?

%output grid to interpolate to
Settings.LzGrid   = logspace(log10(60),log10(1),40);  %vertical wavelength grid
Settings.AltScale = 15:5:80;                          %height scale

%time handling
Settings.DayRange   = [1,30];% [1,365]; %calendar days to range over
Settings.WindowSize = 14;   %calendar days in window
Settings.DayBinStep = 7;    %days stepped between each window
Settings.Years      = 2000:1:2020;  %years of data to use

%lat and lon bin handling
Settings.LatRange   = [-90,90];
Settings.LatBinStep = 5;
Settings.LatBinSize = 5;

Settings.LonRange   = [-180,180];
Settings.LonBinStep = 30;
Settings.LonBinSize = 30;

%progress outputting
Settings.Verbose    = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% work out bin parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Results.LatScale = Settings.LatRange(1):Settings.LatBinStep:Settings.LatRange(2);
Results.LonScale = Settings.LonRange(1):Settings.LonBinStep:Settings.LonRange(2);
Results.AltScale = Settings.AltScale;
Results.LzScale  = Settings.LzGrid;
Results.DayScale = Settings.DayRange(1):Settings.DayBinStep:Settings.DayRange(2);

Results.CutOff = NaN(numel(Results.DayScale),                         ...
                     numel(Results.LonScale),numel(Results.LatScale), ...
                     numel(Results.AltScale),                         ...
                     numel(Settings.LzGrid));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Results.DayScale);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load all data in the window
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clear Store
  if Settings.Verbose; textprogressbar(['Loading ',Settings.Instrument,' data for ',num2str(Settings.WindowSize),' days around d', num2str(Results.DayScale(iDay)),':  ']); end  
  for iYear=1:1:numel(Settings.Years)

    %get data
    DaysInWindow = 0.5.*[-1,1].*Settings.WindowSize + Results.DayScale(iDay);
    ThisChunk = get_limbsounders(datenum(Settings.Years(iYear),1,DaysInWindow),Settings.Instrument, ...
                                 'DateWarning',false,'GetHindleyPWs',true, ...
                                 'HeightScale',Settings.HeightScaleIn, ...
                                 'TimeHandling',1);

    %append to the stored data
    if numel(ThisChunk.Lat) > 0; 
      if ~exist('Store','var'); Store = ThisChunk; else;Store = cat_struct(Store,ThisChunk,1); end
    end

    if Settings.Verbose; textprogressbar(iYear./numel(Settings.Years).*100); end
  end; clear iYear DaysInWindow ThisChunk
  if Settings.Verbose; textprogressbar('!'); end

  %check we have any data
  if ~exist('Store','var');  continue; end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% now, let's ST all the profiles using Alex08, but retaining the full ST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %assign a nominal lat and lon to each profile
  NomLon = median(Store.Lon,2,'omitnan');
  NomLat = median(Store.Lat,2,'omitnan');

  if Settings.Verbose; textprogressbar('Computing noise floor '); k= 0;end
  for iLat=1:1:numel(Results.LatScale)-1
    for iLon=1:1:numel(Results.LonScale)-1

      %get bin bounds
      x = Results.LonScale(iLon)+0.5.*Settings.LonBinStep + [-1,1].*0.5.*(Settings.LonBinSize);
      y = Results.LatScale(iLat)+0.5.*Settings.LatBinStep + [-1,1].*0.5.*(Settings.LatBinSize);

      %find profiles in bin
      if min(x) < -180
        xi = [find(NomLon > min(x)+360);find(NomLon < max(x))];
      elseif max(x) > 180;
        xi = [find(NomLon > min(x));find(NomLon < max(x) - 360)];
      else
        xi = inrange(NomLon,x);
      end
      idx = intersect(xi,inrange(NomLat,y));
      if numel(idx) == 0; continue; end
      RegionIndices = idx(randi(numel(idx),[Settings.NProfilePairs,1]));

      %get randomly-selected profiles
      Data = reduce_struct(Store,RegionIndices,{},1);

      %st them
      GWs = gwanalyse_limb(Data,'Analysis',2,'Filter','Hindley23','FullST',true,'Verbose',false);

      %take adjacent cospectra, then find the covarying amplitude
      %We can pair the last one with the first, as they're random it doesn't matter
      a =      GWs.FullST( 1:1:Settings.NProfilePairs,   :,:);
      b = conj(GWs.FullST([2:1:Settings.NProfilePairs,1],:,:));   
      CovAmp = sqrt(abs(a.*b));

      %then, identify the cutoff percentile
      CutOff = permute(prctile(CovAmp,Settings.Percentile,1),[2,3,1]);

      %interpolate wavelengths onto the storage grid
      %... in Lz...
      CutOff = interp_1d_ndims(1./GWs.Freqs,CutOff,Settings.LzGrid,2);
      %... and in height...
      z = Settings.HeightScaleIn; good = ~isnan(z); z = z(good); CutOff = CutOff(good,:); 
      if numel(good) < 2; continue; end
      CutOff = interp_1d_ndims(z,CutOff,Results.AltScale,1);      

      %and store
      Results.CutOff(iDay,iLon,iLat,:,:) = CutOff;

      clear a b CutOff GWs Data CovAmp z good RegionIndices

      if Settings.Verbose; k = k+1;textprogressbar(k./((numel(Results.LatScale)-1) .* (numel(Results.LonScale)-1)).*100); end
    end
  end; clear k iLat iLon Store NomLon NomLat
  if Settings.Verbose; textprogressbar(100); textprogressbar('!'); end  



end; clear iDay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the results!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%shift axes to be bin centres, not bin edges

Results.LatScale = Results.LatScale+Settings.LatBinStep./2;
Results.LonScale = Results.LonScale+Settings.LonBinStep./2;
Results.DayScale = Results.DayScale+Settings.DayBinStep./2;

%remove bins that lie off the end of the world
Results.CutOff   = Results.CutOff(:,1:end-1,1:end-1,:,:);
Results.LonScale = Results.LonScale(1:end-1);
Results.LatScale = Results.LatScale(1:end-1);


save(Settings.OutFile,'-struct','Results','-v7.3')