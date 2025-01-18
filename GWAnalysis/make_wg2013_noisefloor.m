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
Settings.OutFile = ['noisefloor_',Settings.Instrument,'.mat'];

%statistical choices
Settings.Percentile    = 99;    %what percentile of the distribution is the stored cutoff?
Settings.NProfilePairs = 10000;  %how many profile-pairs should we use in each box?

%output interpolation choices
Settings.LzGrid   = logspace(log10(60),log10(1),40);  %vertical wavelength grid
Settings.AltScale = 15:5:80;                          %height scale

%time handling
Settings.DayRange   = [1,365]; %calendar days to range over
Settings.WindowSize = 15;    %calendar days in window
Settings.DayBinStep = 7;    %days stepped between each window
Settings.Years      = 2000:1:2020;  %years of data to use

%lat and lon bin handling
Settings.LatRange   = [-90,90];
Settings.LatBinStep = 5;
Settings.LatBinSize = 5;

Settings.LonRange   = [-180,180];
Settings.LonBinStep = 20;
Settings.LonBinSize = 40;

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
  if Settings.Verbose; textprogressbar(['Loading data for ',num2str(Settings.WindowSize),' day range centred on calendar day ', num2str(Results.DayScale(iDay)),':  ']); end
  for iYear=1:1:numel(Settings.Years)

    %get data
    DaysInWindow = 0.5.*[-1,1].*Settings.WindowSize + Results.DayScale(iDay);
    ThisChunk = get_limbsounders(datenum(Settings.Years(iYear),1,[DaysInWindow]),Settings.Instrument, ...
                                 'DateWarning',false,'GetHindleyPWs',true);

    %append to the stored data
    if numel(ThisChunk.Lat) > 0; 
      if ~exist('Store'); Store = ThisChunk; else;Store = cat_struct(Store,ThisChunk,1); end
    end

    if Settings.Verbose; textprogressbar(iYear./numel(Settings.Years).*100); end
  end; clear iYear DaysInWindow ThisChunk
  if Settings.Verbose; textprogressbar('!'); end

  %check we have any data
  if ~exist('Store');  continue; end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% identify the points present in each region
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %create an array to store point indices
  RegionIndices = NaN(numel(Results.LatScale) .* numel(Results.LonScale),Settings.NProfilePairs);

  %assign a nominal lat and lon to each profile
  NomLon = nanmedian(Store.Lon,2);
  NomLat = nanmedian(Store.Lat,2);

  %fill the indices
  iBox = 0;
  for iLat=1:1:numel(Results.LatScale)
    for iLon=1:1:numel(Results.LonScale)

        %get bin bounds
        iBox = iBox+1;
        x = Results.LonScale(iLon)+0.5.*Settings.LonBinStep + [-1,1].*0.5.*(Settings.LonBinSize);
        y = Results.LatScale(iLat)+0.5.*Settings.LatBinStep + [-1,1].*0.5.*(Settings.LatBinSize);

        %find data in bin
        if min(x) < -180      
          xi = [find(NomLon > min(x)+360);find(NomLon < max(x))];
        elseif max(x) > 180;
          xi = [find(NomLon > min(x));find(NomLon < max(x) - 360)];
        else
          xi = inrange(NomLon,x);
        end
        idx = intersect(xi,inrange(NomLat,y));
        if numel(idx) > 0;

          %choose a random order of profiles, with repeat sampling permitted
          RegionIndices(iBox,:) = idx(randi(numel(idx),[Settings.NProfilePairs,1]));
        end
        

    end; clear iLon idx x y xi
  end; clear iLat iBox NomLat NomLon

  %reshape to give access to lons and lats in a logical way
  RegionIndices = reshape(RegionIndices,[numel(Results.LonScale),numel(Results.LatScale),Settings.NProfilePairs]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% now, let's ST all the profiles using Alex08, but retaining the full ST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if Settings.Verbose; textprogressbar('Computing noise floor '); k= 0;end
  for iLat=1:1:numel(Results.LatScale)
    for iLon=1:1:numel(Results.LonScale)

      %get randomly-selected profiles
      if nansum(RegionIndices(iLon,iLat,:),"all") == 0; continue; end
      Data = reduce_struct(Store,RegionIndices(iLon,iLat,:),{},1);

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
      CutOff = interp_1d_ndims(nanmean(GWs.Alt,1),CutOff,Results.AltScale,1);      

      %and store
      Results.CutOff(iDay,iLon,iLat,:,:) = CutOff;

      clear a b CutOff GWs Data CovAmp

      if Settings.Verbose; k = k+1;textprogressbar(k./(numel(Results.LatScale) .* numel(Results.LonScale)).*100); end
    end
  end; clear k iLat iLon RegionIndices Store
  if Settings.Verbose; textprogressbar('!'); end  



end; clear iDay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the results!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%shift axes to be bin centres, not bin edges

Results.LatScale = Results.LatScale+Settings.LatBinStep./2;
Results.LonScale = Results.LonScale+Settings.LonBinStep./2;
Results.DayScale = Results.DayScale+Settings.DayBinStep./2;

Results.CutOff = Results.CutOff(1:end-1,1:end-1,1:end-1,:,:);


save(Settings.OutFile,'-struct','Results','-v7.3')