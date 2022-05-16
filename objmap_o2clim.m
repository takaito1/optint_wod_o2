%%%% objmap_o2clim
%%%% perform objmective map and save each z levels
%%%%

function dummy = objmap_o2clim(m,klevel)

month={'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
fn='o2_clim_osdctd_pass2.nc';
Nlev=47;
! mkdir -p intermed_files

% display info
disp('Optimal interpolation of WOD statistical mean data.');
disp('OI is applied to each depth level (k) and each basin separately. ');

% set parameters
S2N=5;       % signal to noise ratio in obs
L=750.0e3;   % e-folding scale for Gaussian covariance function
disp(['Assumed signal to noise ratio is ',num2str(S2N),' (=1/epsilon)']);
disp(['e-folding scale for autocovariance is ',num2str(L*1e-3),' km']);

% get mask
bn=textread('basin_name.txt','%s','delimiter','\n');
mask0=ncread('basin_mask_01.nc','basin_mask');
for k=1:Nlev
   maskk=mask0(:,:,k);
   bind0=unique(maskk(:));
   bindK{k}.data=bind0(2:end);
end
maskNz=mask0(:,:,1:Nlev);

x=ncread('basin_mask_01.nc','lon');
y=ncread('basin_mask_01.nc','lat');
z0=ncread('basin_mask_01.nc','depth');
z=z0(1:Nlev);
[yy,xx]=meshgrid(y,x);
Nx=length(x);
Ny=length(y);
Nz=length(z);

% load data
mn=squeeze(nanmean(ncread(fn,'o_mn',[1 1 1 m],[360 180 Nz 1]),4));

% loop over z levels
for k=[klevel]
 disp('----------------------------------')
 % loop over all basins (1=atl, 2=pac, 3=ind, etc)
 bind=bindK{k}.data;
 wn=['intermed_files/o2map_mon',num2str(m),'_klev',num2str(k),'.mat']
 o2map=zeros(360,180);
 for b=1:length(bind)
   disp(['k=',num2str(k),' : mapping ',bn{bind(b)}]);

   % extract data
   dd=squeeze(mn(:,:,k));
   bb=mask0(:,:,k);
   b0=bb(:);
   d0=dd(:);
   x0=xx(:);
   y0=yy(:);
   Y0=zeros(size(x0));
   Q0=zeros(size(x0));
   Z0=zeros(size(x0));
   E0=zeros(size(x0));

   % use data only from the adjacent basin
   if b<=3|b==10
      I=find(~isnan(d0)&(b0==bind(b)|y0<-30));
      J=find(b0==bind(b));
   else
      I=find(~isnan(d0)&b0==bind(b));
      J=find(b0==bind(b));
   end

   if ~isempty(I)
   d2=d0(I);
   x2=x0(I);
   y2=y0(I);
   N2=length(d2);
   x0=x0(J);
   y0=y0(J);
   N0=length(x0);

   % distance between two points in obs
   % D(m,n): m,n: obs point ref, square matrix (N2,N2)
   f=pi/180;
   dlon=repmat(x2,[1 N2]) - repmat(x2',[N2 1]);
   ds=acos(sin(f*y2)*sin(f*y2')+cos(f*y2)*cos(f*y2').*cos(f*dlon));
   dl=6.371e6*real(ds);

   % obs-obs gaussian covariance function
   UU=exp(-.5*(dl/L).^2) + 1/S2N*eye(N2);

   % D(m,n): m: grid point, n: obs point ref, non-square (N0,N2)
   f=pi/180;
   dlon=repmat(x0,[1 N2]) - repmat(x2',[N0 1]);
   ds=acos(sin(f*y0)*sin(f*y2')+cos(f*y0)*cos(f*y2').*cos(f*dlon));
   dl=6.371e6*real(ds);
   % obs-grid gaussian covariance function
   VU=exp(-.5*(dl/L).^2);

   % mapping matrix
   invUU=inv(UU);

   % lagrange multiplier part
   v=ones(1,N2);
   lmd1=(1-VU*invUU*v');
   lmd2=(v*invUU*v');
   A = (VU+lmd1*v/lmd2)*invUU;
   disp(['... dimension of OI matrix is (grid x obs) ',num2str(size(A))]);

   %% map it
   Y = A*d2;
   Y0(J)=Y;
   Test = reshape(Y0,[Nx,Ny]);
   o2map=o2map+Test;
   o2map(maskNz(:,:,k)==0)=NaN;

   else
      disp(['... no obs is available in this basin']);
   end % IF ~isempty(I) %
 end % loop over basin 
 o2_mn=mn;
 save('-v7.3',wn,'o2map','o2_mn','x','y','z');
end % loop over depth
