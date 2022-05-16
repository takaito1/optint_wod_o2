close all;
clear all;

% load the data
mask=ncread('basin_mask_01.nc','basin_mask');
mask0=mask(:,:,1:47);
yr=1967:2018;
for t=1:length(yr)
   rn=['intermed_files/O2_OI1x1_',num2str(yr(t)),'.mat'];
   load(rn);
   % apply mask
   o2map(mask0==0)=NaN;
   R2map(mask0==0)=NaN;
   an(:,:,:,t)=o2map;
   R2(:,:,:,t)=R2map;
end

% smoothing for final products
yc=-50; % Southern Ocean boundary
J=find(y>(yc-10)&y<(yc+10));
an(:,J,:,:)=movmean(an(:,J,:,:),[1 1],2,'omitnan');
R2(:,J,:,:)=movmean(R2(:,J,:,:),[1 1],2,'omitnan');
J=find(y>(yc-6)&y<(yc+6));
an(:,J,:,:)=movmean(an(:,J,:,:),[1 1],2,'omitnan');
R2(:,J,:,:)=movmean(R2(:,J,:,:),[1 1],2,'omitnan');

% vertical smoothing (3 point running mean)
an=movmean(an,[1 1],3,'omitnan');
R2=movmean(R2,[1 1],3,'omitnan');

% save them in final array
O2.an=an;
O2.R2=R2;

%% Generate NetCDF file
fn=['o2_grid_1x1_',date,'.nc'];
scope = netcdf.create([fn],'netcdf4');
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
fillValue = -99999;

% 3-3. Define dimensions
X=x;
Y=y;
Z=z;
nz=length(z);
dz=diff(z);
ZF(2:nz)=.5*(z(1:end-1)+z(2:end));
ZF(end+1)=ZF(end)+dz(end);

T=(yr-1950)*365.25; % days since 1950-01-01

dimidX = netcdf.defDim(scope,'lon',length(X));
dimidY = netcdf.defDim(scope,'lat',length(Y));
if length(Z) >= 2
   dimidZ = netcdf.defDim(scope,'depth',length(Z)); % turn on if you have depth directions
   dimidZF = netcdf.defDim(scope,'bnds',2); % turn on if you have depth directions
end
dimidT = netcdf.defDim(scope,'time',length(T));

% 3-4. Define coordinates and time axis
varid = netcdf.defVar(scope,'lon','double',[dimidX]);
netcdf.putAtt(scope,varid,'standard_name','lon');
netcdf.putAtt(scope,varid,'long_name','longitude');
netcdf.putAtt(scope,varid,'units','degrees_east');
netcdf.defVarFill(scope,varid,false,fillValue);
netcdf.putVar(scope,varid,X);

varid = netcdf.defVar(scope,'lat','double',[dimidY]);
netcdf.putAtt(scope,varid,'standard_name','lat');
netcdf.putAtt(scope,varid,'long_name','latitude');
netcdf.putAtt(scope,varid,'units','degrees_north');
netcdf.defVarFill(scope,varid,false,fillValue);
netcdf.putVar(scope,varid,Y);

varid = netcdf.defVar(scope,'depth','double',[dimidZ]);
netcdf.putAtt(scope,varid,'standard_name','depth');
netcdf.putAtt(scope,varid,'long_name','depth from the surface ocean');
netcdf.putAtt(scope,varid,'units','m');
netcdf.putAtt(scope,varid,'bounds','depth_bnds');
netcdf.defVarFill(scope,varid,false,fillValue);
netcdf.putVar(scope,varid,Z);

Zbnds(:,1)=ZF(1:end-1);
Zbnds(:,2)=ZF(2:end);
varid = netcdf.defVar(scope,'depth_bnds','double',[dimidZF dimidZ]);
netcdf.putAtt(scope,varid,'standard_name','depth');
netcdf.putAtt(scope,varid,'units','m');
netcdf.defVarFill(scope,varid,false,fillValue);
netcdf.putVar(scope,varid,Zbnds');
   
varid = netcdf.defVar(scope,'time','double',[dimidT]);
netcdf.putAtt(scope,varid,'standard_name','time');
netcdf.putAtt(scope,varid,'long_name','days');
netcdf.putAtt(scope,varid,'units','days since 1950-01-01');
netcdf.defVarFill(scope,varid,false,fillValue);
netcdf.putVar(scope,varid,T);

% 3-5. Define variable attributes
netcdf.putAtt(scope,NC_GLOBAL,'title',['objectively mapped dissolved oxygen based on World Ocean Database 2018'])
netcdf.putAtt(scope,NC_GLOBAL,'Conventions','CF-1.6')
netcdf.putAtt(scope,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'))
netcdf.close(scope)

%% 3-6. Define and store variables
V = O2.an;
unit='micro-molO2/L';
varname = 'o2';
long_name = 'objective map of dissolved oxygen anomaly';
%%%%%%%%%%%%%%%%%%%%%%%%%%
scope = netcdf.open([fn],'WRITE'); % here we use 'WRITE' because the file already exists
varid = netcdf.defVar(scope,varname,'double',[dimidX,dimidY,dimidZ,dimidT]);
netcdf.putAtt(scope,varid,'long_name',long_name);
netcdf.putAtt(scope,varid,'units',unit);
netcdf.defVarFill(scope,varid,false,fillValue);
% enter fill value
V(isnan(V)) = fillValue;
V(V==0) = fillValue;
netcdf.putVar(scope,varid,V);
netcdf.close(scope) % now insert three-dimensional data to NetCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = O2.R2;
unit='ND';
varname = 'o2_R2';
long_name = 'coefficient of determination for objective map';
%%%%%%%%%%%%%%%%%%%%%%%%%%
scope = netcdf.open([fn],'WRITE'); % here we use 'WRITE' because the file already exists
varid = netcdf.defVar(scope,varname,'double',[dimidX,dimidY,dimidZ,dimidT]);
netcdf.putAtt(scope,varid,'long_name',long_name);
netcdf.putAtt(scope,varid,'units',unit);
netcdf.defVarFill(scope,varid,false,fillValue);
% enter fill value
V(isnan(V)) = fillValue;
V(V==0) = fillValue;
netcdf.putVar(scope,varid,V);
netcdf.close(scope) % now insert three-dimensional data to NetCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%

