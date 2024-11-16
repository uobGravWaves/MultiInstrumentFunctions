function [Error,Tropopause,Stratopause] = module_pauses(Settings,LonPoints,LatPoints,TimePoints,BBox)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get_context() module to compute tropopause and
%stratopause pressure level
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/05/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%check we fed in a time - this is required
if sum(isnan(TimePoints)) == numel(TimePoints);
  warning('Pauses: no TimePoints provided. Skipping.')
  Error = 1; Tropopause = []; Stratopause = [];
else

  %setup
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %get ERA5 global temperature
  I = create_era5_interpolant(TimePoints,Settings,'Pauses',BBox);

  %compute pressure. We can ignore lnsp as both 'pauses should be above the region it matters.
  Pressure = ecmwf_prs_v3(137);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %stratopause
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %method: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2011JD016893

  %interpolate the data to a regular ~1km height grid between 25 and 80km altitude
  NewP = h2p(25:1:80)'; NewZ = p2h(NewP);
  lat = repmat( LatPoints,[ones(ndims(LatPoints),1);numel(NewP)]');
  lon = repmat( LonPoints,[ones(ndims(LatPoints),1);numel(NewP)]');
  t   = repmat(TimePoints,[ones(ndims(LatPoints),1);numel(NewP)]');
  p   = permute(repmat(NewP,[1,size(LatPoints)]),[2:ndims(lat),1]);
  T   = I.T(lon,lat,t,p);
  clear lon lat t p

  %reshape to put height first
  T  = permute(T ,[ndims(T ),1:1:ndims(T )-1]);

  %smooth by 11km. Remember we don't know how many dimensions we have...
  sz = size(T);
  Ts = reshape(T,sz(1),prod(sz(2:end)));
  Ts = smoothdata(Ts,1,'movmean');
  Ts = reshape(Ts,sz);

  %find maximum in each profile
  [~,idx] = max(Ts,[],1);

  %check 5 levels above and below:
  %5 levels above must have -ve lapse rate
  %5 levels below must have +ve lapse rate
  dTdZ = diff(Ts,1,1);
  dTdZ = cat(1,zeros(size(idx)),dTdZ); %add extra level so points line up, rather than half-levels

  Stratopause = NaN(size(idx));
  for iProf=1:1:prod(size(idx));

    Above = idx(iProf)+1:1:idx(iProf)+5; Above = Above(Above > 0 & Above < size(NewP,1));
    Below = idx(iProf)-5:1:idx(iProf)-1; Below = Below(Below > 0 & Below < size(NewP,1));

    Above = -dTdZ(Above,iProf); Below = dTdZ(Below,iProf); %note - sign on Above

    if min(Above) > 0 & min(Below) > 0;
      %remove anything outside +/- 15 km from peak, for safety below
      T(NewZ < NewZ(idx(iProf))-15,iProf) = NaN;
      T(NewZ > NewZ(idx(iProf))+15,iProf) = NaN;

      %then find the maximum in the unsmoothed data
      [~, Stratopause(iProf)] = max(T(:,iProf),[],1);
    end
  end;

  %convert to height, and fill small gaps from the pass condition above (these are usually <1% of the data)
  Stratopause(~isnan(Stratopause)) = NewZ(Stratopause(~isnan(Stratopause)));
  Stratopause = fillmissing(Stratopause,'linear');

  %drop unnecessary dimensions added above, and convert back to pressure
  Stratopause = h2p(permute(Stratopause,[2:1:ndims(Stratopause),1]));

  %tidy up
  clear Above Below dTdZ idx NewP NewZ T Ts sz 


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %tropopause
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %computes tropopause based on WMO definition, approximately following doi:10.1029/2003GL018240 but
  %with modifications for speed

  %put temperature data onto native ERA5 pressure grid between 700 hPa and 10 hPa
  %we want to work up in height, so flip Pressure
  Pressure = ecmwf_prs_v3(137);
  Pressure = Pressure(end:-1:1);
  lat = repmat( LatPoints,[ones(ndims(LatPoints),1);numel(Pressure)]');
  lon = repmat( LonPoints,[ones(ndims(LatPoints),1);numel(Pressure)]');
  t   = repmat(TimePoints,[ones(ndims(LatPoints),1);numel(Pressure)]');
  p   = permute(repmat(Pressure,[1,size(LatPoints)]),[2:ndims(lat),1]);
  T   = I.T(lon,lat,t,p);
  clear lon lat t p

  %reshape to put height first, then turn into lines
  T = permute(T ,[ndims(T ),1:1:ndims(T )-1]);
  sz = size(T);
  T = reshape(T,sz(1),prod(sz(2:end)));

  %compute lapse rate.
  dT = diff(T,1,1);
  dZ = diff([p2h(Pressure)]);
  Gamma = dT .* NaN;
  for iLev=1:1:numel(dZ)-1; Gamma(iLev,:) = dT(iLev,:)./dZ(iLev); end;
  clear dT dZ iLev

  %create an array to store our tropopause levels, then loop over the data to find them
  %we are working UPWARDS
  Tropopause = NaN(size(T,2),1);

  for iLev=1:1:numel(Pressure)

    %if pressure > 700hPa or <10hPa, or if we've already found the t'pause everywhere, skip
    if Pressure(iLev) > 700;           continue; end
    if Pressure(iLev) <  10;           continue; end
    if sum(isnan(Tropopause(:))) == 0; continue; end

    %check if Gamma is less than 2 anywhere at this level
    idx = find(Gamma(iLev,:) > -2);
    if numel(idx) == 0; continue; end %none at this level

    %remove any columns we already found
    Found = find(~isnan(Tropopause));
    [~,Remove] = intersect(idx,Found);
    idx(Remove) = [];
    clear Remove

    %for each element where the above criterion is met, check if the layer
    %2km higher also meets it

    %find which level is 2km above
    Z = p2h(Pressure(iLev));
    jLev = closest(p2h(Pressure),Z+2);

    %find all the columns where the criterion remains met ON AVERAGE for these 2km above
    Good = find(nanmean(Gamma(iLev:jLev,idx),1) > -2);
    if numel(Good) < 2 ; continue; end %this needs to be 2 because of ambiguity in the array operations below.
    %This leaves a small number of NaNs (<< 1%), which we interpolate over below
    idx= idx(Good);

    %for the remaining columns, find where the gradient crossed above -2 by linear interpolation
    G  = Gamma(iLev:jLev,idx);
    p  = linspace(Pressure(iLev),Pressure(jLev),100); %put onto 100 levels
    Gi = interp1(Pressure(iLev:jLev),Gamma(iLev:jLev,idx),p)+2;
    [~,minpidx] = min(abs(Gi),[],1);

    Tropopause(idx) = p(minpidx);
  end

  %fill small gaps due to put back to the original shape, and return
  Tropopause = fillmissing(Tropopause,'linear');
  if (ndims(Tropopause) == 2 & size(Tropopause,2) == 1) | (ndims(Tropopause) == 2 & size(Tropopause,1) == 1)
  else Tropopause = reshape(Tropopause,sz(2:end));
  end

  %tidy up
  clear I Pressure Found G Gamma Gi Good idx iLev iProf jLev minpidx p sz T Z


  Error = 0;


end




return
