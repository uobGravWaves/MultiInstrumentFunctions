function Sentinel = module_sentinel(Settings,LonPoints,LatPoints,TimePoints,BBox)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get_context() module to high-resolution Sentinel surface imagery
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/05/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%several checks to do here, so let's set a flag to keep track
Fail = 0;
Sentinel = [];




%do we already have a file?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Sentinel_Reuse == true

  %load the file into memory
  if exist(Settings.Sentinel_OutFile,'file')
    Sentinel = flipud(imread(Settings.Sentinel_OutFile));
    %check the size matches
    if ~isequal(size(Sentinel),size(repmat(LonPoints,1,1,3)));
      warning('Sentinel: previously-downloaded image is not the right size, getting new data')
      delete(Settings.Sentinel_OutFile)
    else
      warning('Sentinel: reusing previously-downloaded image')
      Fail = 1; %because we don't need to get the data
    end
  else
    warning('Sentinel: no previously-downloaded image, getting new data')
  end
end



%the data MUST be a valid and regular rectangle. Check this first.
%we're looking for a single unique diff() value in each dimension
%and exactly two dimensions, but we need to do some extra work to deal
%with numerical precision issues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(LonPoints) ~= 2
  warning('Sentinel: data points must form a regularly-spaced rectangle to get Sentinel data; shape failed. Skipping.')
  Fail = 1;
end

a = unique(diff(LonPoints,1,1)); b = unique(diff(LonPoints,1,2));
c = unique(diff(LatPoints,1,1)); d = unique(diff(LatPoints,1,2));
ar = range(a)./a; br = range(b)./b;  cr = range(c)./c; dr = range(d)./d;
e = [ar;br;cr;dr]; e(isnan(e)) = [];
if max(e) > 1e-5;
  warning('Sentinel: data points must form a regularly-spaced rectangle to get Sentinel data; diff failed. Skipping.')
  Fail = 1;
end
clear ar br cr dr e c d

%check we fed in username and password - required
if numel(Settings.Sentinel_ID{1}) == 0 | numel(Settings.Sentinel_ID{2}) == 0;
  warning('Sentinel: no Sentinel user credentials supplied. Skipping.')
  Fail = 1;
end



%did we request a date at all? If not, use default (1st April 2023 - this is arbitrary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SentinelTime = TimePoints;
if numel(SentinelTime) == 1 && isnan(SentinelTime); SentinelTime(:) = datenum(2023,4,1); end





%what date do we want the data for? We must only be requesting one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nanvar(SentinelTime,[],'all') ~= 0
  warning('Sentinel: requesting multiple dates for Sentinel data; using earliest-specified in valid range only')
  SentinelTime = nanmin(SentinelTime(SentinelTime > datenum(2023,1,1)));
else
  SentinelTime = nanmin(SentinelTime,[],'all');
end


%is the requested date post-dataset start in 2023? If not, shift to the same DoY in 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nanmin(SentinelTime,[],'all') < datenum(2023,1,1);
  warning('Sentinel: requested dates before dataset start; shifting pre-2023 date into 2023')
  [y,~,~] = datevec(SentinelTime); dn = date2doy(SentinelTime);  y(y < 2023) = 2023;
  SentinelTime = datenum(y,1,dn);
  clear y dn
end


%work out the required resolution
% This requires working out if our data are lat-major or lon-major
% if lon is the x-axis, then max(b) will be greater than max(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ndims(LonPoints) > 2;
  warning('Sentinel: Sentinel data cannot be requested in >2 dimensions. Skipping.')
  Fail = 1;
else
  if max(b) > max(a); Resolution = size(LonPoints');
  else                Resolution = size(LonPoints);
  end

  %hence check if we're in the permitted range of resolutions
  boxwidthx = deg2km(distance(BBox(4),BBox(1),BBox(4),BBox(3),'degrees')).*1000;
  boxwidthy = deg2km(distance(BBox(2),BBox(1),BBox(4),BBox(1),'degrees')).*1000;
  resx = boxwidthx./Resolution(1);
  resy = boxwidthy./Resolution(2);
  if resx > 1600; Fail = 1; warning('Sentinel: longitude resolution is too coarse. Skipping'); end
  if resy > 1600; Fail = 1; warning('Sentinel: latitude resolution is too coarse. Skipping'); end
  clear boxwidthx boxwidthy resx resy
end
clear a b Resolution

%ok, we're almsot there. Work out how many tiles 
%we need, and how they are laid out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Fail == 0; %so we skip this check if it's not needed, since it queries the user

  %we can only request 2500x2500 points in a single pass - otherwise we need to make multiple API calls


  if size(LonPoints,1) > 2500 | size(LonPoints,2) > 2500;  %use multiple passes

    %work out the properties of the multiple passes
    x = ceil(size(LonPoints,1)./2500); y = ceil(size(LonPoints,2)./2500);
    NPasses = x*y;
    BBoxes = NaN(NPasses,4);

    disp(['Sentinel: API allows max 2500 points per dimension, so this will require ',num2str(NPasses),' calls.'])
    Input = input(['Are you certain? Enter 1 to confirm.']);
    if Input ~= 1; Fail = 1; else
      k = 1;
      for iX=1:1:x
        for iY=1:1:y

          %define the bounding box for this pass
          idxX = [1,0]+[iX-1,iX].*2500;
          idxY = [1,0]+[iY-1,iY].*2500;
          idxX(idxX > size(LonPoints,1)) = size(LonPoints,1); idxX = idxX(1):1:idxX(end);
          idxY(idxY > size(LonPoints,2)) = size(LonPoints,2); idxY = idxY(1):1:idxY(end);

          BBoxes(k,:) = [LonPoints(idxX(  1),idxY(  1)), LatPoints(idxX(  1),idxY(  1)), ...
            LonPoints(idxX(end),idxY(end)), LatPoints(idxX(end),idxY(end))];

          %identify which merged output points these pixels go into
          a = reshape(1:1:numel(LonPoints),size(LonPoints)); %list of ALL points
          a = a(idxX,:); a = a(:,idxY);
          IndexStore.(['store',num2str(k)]) = [idxX(1),idxX(end),idxY(1),idxY(end)];
          k = k+1;
        end
      end
    end
    clear k iX iY x y idxX idxY

  else
    %just use the master bounding box and list of output points
    NPasses = 1;
    BBoxes = BBox;
    IndexStore.(['store1']) = [1,size(LonPoints,1),1,size(LonPoints,2)];
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% we're there! time to go!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ok, let's go
if Fail == 0;
  for iPass = 1:1:NPasses
    %downloading uses the Python API for Copernicus. This took me hours to figure out the syntax for.

    %generate the API calling script
    Resolution = IndexStore.(['store',num2str(iPass)]); Resolution = Resolution([4,2]) - Resolution([3,1]) + [1,1];
    WorkingFileName = ['working',num2str(randi(1000,1,1)),'_',num2str(randi(1000,1,1)),'.png'];
    Script = sentinel_script(Settings.Sentinel_ID,WorkingFileName,Settings.Sentinel_Gain,BBoxes(iPass,:),Resolution,SentinelTime)';
    ScriptFile = "get_sentinel_"+strrep(num2str(datenum(now)),'.','')+".py";
    writelines(Script,ScriptFile)

    %now run the script, and delete it
    pyenv(ExecutionMode="OutOfProcess");
    pyrunfile(ScriptFile);
    delete(ScriptFile);

    %tidy up Python
    terminate(pyenv);

    %and load the image into memory
    SentinelData.(['store',num2str(iPass)]) = flipud(imread(WorkingFileName));
    delete(WorkingFileName)
    clear WorkingFileName

  end


  %merge passes then output
  Sentinel = zeros([size(LonPoints),3],'uint8');
  for iPass=1:1:NPasses
    idx = IndexStore.(['store',num2str(iPass)]);
    Sentinel(idx(1):idx(2),idx(3):idx(4),:) = SentinelData.(['store',num2str(iPass)]);
  end;

  %also write a file image
  imwrite(flipud(Sentinel),Settings.Sentinel_OutFile);


end
clear Fail Script ScriptFile a b Resolution Input  iPass idx BBoxes IndexStore NPasses SentinelData SentinelTime



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done, return :-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write Sentinel python script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Script = sentinel_script(Sentinel_ID,Sentinel_OutFile,Gain,BBox,Resolution,SentinelTime)

TimeStringStart = [datestr(SentinelTime-14,'yyyy-mm-dd'),'T00:00:00Z'];
TimeStringEnd   = [datestr(SentinelTime+14,'yyyy-mm-dd'),'T23:59:59Z'];


Script        = "from scipy.io import savemat";
Script(end+1) = "import numpy as np";
Script(end+1) = "from oauthlib.oauth2 import BackendApplicationClient";
Script(end+1) = "from requests_oauthlib import OAuth2Session";
Script(end+1) = "";
Script(end+1) = "";
Script(end+1) = "# Your client credentials";
Script(end+1) = "client_id = '"+Sentinel_ID{1}+"'";
Script(end+1) = "client_secret = '"+Sentinel_ID{2}+"'";
Script(end+1) = "";
Script(end+1) = "# Create a session";
Script(end+1) = "client = BackendApplicationClient(client_id=client_id)";
Script(end+1) = "oauth = OAuth2Session(client=client)";
Script(end+1) = "";
Script(end+1) = "# Get token for the session";
Script(end+1) = "token = oauth.fetch_token(token_url='https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token',";
Script(end+1) = "                          client_secret=client_secret, include_client_id=True)";
Script(end+1) = "";
Script(end+1) = "";
Script(end+1) = "response = oauth.get('https://sh.dataspace.copernicus.eu/configuration/v1/wms/instances')";
Script(end+1) = "";
Script(end+1) = "evalscript = '''";
Script(end+1) = "//VERSION=3";
Script(end+1) = "function setup() {";
Script(end+1) = "  return {";
Script(end+1) = "    input: ['B02', 'B03', 'B04'],";
Script(end+1) = "    output: { bands: 3 }";
Script(end+1) = "  };";
Script(end+1) = "}";
Script(end+1) = "";
Script(end+1) = "function evaluatePixel(sample) {";
Script(end+1) = "  let gain = "+num2str(Gain);
Script(end+1) = "  return [gain * sample.B04/10000, gain * sample.B03/10000, gain * sample.B02/10000];";
Script(end+1) = "}";
Script(end+1) = "'''";
Script(end+1) = "";
Script(end+1) = "request = {";
Script(end+1) = "  'input': {";
Script(end+1) = "    'bounds': {";
Script(end+1) = "      'bbox': [";
Script(end+1) = "        "+BBox(1)+",";
Script(end+1) = "        "+BBox(2)+",";
Script(end+1) = "        "+BBox(3)+",";
Script(end+1) = "        "+BBox(4)+"";
Script(end+1) = "      ]";
Script(end+1) = "    },";
Script(end+1) = "    'data': [";
Script(end+1) = "      {";
Script(end+1) = "        'dataFilter': {";
Script(end+1) = "          'timeRange': {";
Script(end+1) = "            'from': '"+TimeStringStart+"',";
Script(end+1) = "            'to': '"+TimeStringEnd+"'";
Script(end+1) = "          }";
Script(end+1) = "        },";
Script(end+1) = "        'type': 'byoc-5460de54-082e-473a-b6ea-d5cbe3c17cca'";
Script(end+1) = "      }";
Script(end+1) = "    ]";
Script(end+1) = "  },";
Script(end+1) = "  'output': {";
Script(end+1) = "    'width': "+Resolution(1)+",";
Script(end+1) = "    'height': "+Resolution(2)+",";
Script(end+1) = "    'responses': [{'format': {'type': 'image/png'}}],";
Script(end+1) = "  },";
Script(end+1) = "  'evalscript': evalscript,";
Script(end+1) = "}";
Script(end+1) = "";
Script(end+1) = "url = 'https://sh.dataspace.copernicus.eu/api/v1/process'";
Script(end+1) = "response = oauth.post(url, json=request)";
Script(end+1) = "";
Script(end+1) = "if response.status_code == 400:";
Script(end+1) = "  print(response.text)";
Script(end+1) = "";
Script(end+1) = "with open('"+Sentinel_OutFile+"','wb') as f:";
Script(end+1) = "  f.write(response.content)";
Script(end+1) = "";


return



return
end