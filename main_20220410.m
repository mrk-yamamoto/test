
function main
make_rcvData_HFR_1D
end

function make_rcvData_HFR_1D

f0  = 7.5e6;%params_Vera.fc;%6.3e6;                  %  Transducer center frequency [Hz]
fs  = 1/3.3333e-8;%1/params_Vera.dt;%100e6;                %  Sampling frequency [Hz]
c   = 1480;%params_Vera.C;%1540;                  %  Speed of sound [m/s]
lambda  = c/f0;             %  Wavelength [m]
width   = 2e-4;%params_Vera.dele;%lambda;            %  Width of element
% element_height=5/1000;   %  Height of element [m]
% kerf    = 2e-4;%params_Vera.xele(2)-params_Vera.xele(1);%0.27/1000-lambda;%0.0/1000;           %  Kerf [m]
z_focus = 30e-3;         % 
N_elements  = 192;          %  Number of physical elements
N_active    = 128;%params_Vera.Nele_R;%128;             %  Number of active elements
xmit_N_active= 128;%params_Vera.Nele_R;%128;       %  Number of active transmit elements for constant F#
rec_N_active = 128;%params_Vera.Nele_R;%128;        %  Number of active receive elements for constant F#
% Nxsub=1;
% Nysub=10;
% d_x     = (width+kerf);
dz      = lambda/2;
Ndcm_z  = ceil(dz/(c/fs));
theta   = [-18,-12,-6,0,6,12,18]/180*pi;

fs      = fs*5; %% 2020/5/1 without this, tail of PSF does not appeal

%% imaging parameter 
%% DAS
beamWx  = 1*lambda; % 0.2 mm
beamWz  = 1*lambda; % 0.2 mm

PPR     = 6/6%[3/6,4/6]%:1/6:1%1/6:1/6:1
xR      = 1/5;%1/10;

% pxl def
dx_pxl  = beamWx*PPR*xR;
dz_pxl  = beamWz*PPR*xR;

x_pxl   = 0:dx_pxl:beamWx*50;
z_pxl   = 0:dz_pxl:beamWz*50;
x_pxl   = x_pxl-mean(x_pxl);
z_pxl   = z_pxl-mean(z_pxl)+z_focus;
DAS     = zeros(length(x_pxl),length(z_pxl));
 
for icmp=4%1:7
 
    %% transmit delay
    x = width*(1:N_active);
    x = x-mean(x);
%   [elex,eley] = meshgrid(x,y);
    elex = x;

    %% 1 point focus
    delay_z = sqrt(elex.^2+ z_focus^2)-z_focus;
    %% plane wave
    delay_z = elex*tan(theta(icmp)); % 20220410
    delay_T = delay_z;
    
    %% dist
    % ptm_pos = [0,0,z_focus];
    r_sct = [0,0,z_focus];
    ptm_pos = r_sct;
    ptm_amp = 1;
    dist = sqrt((elex-r_sct(1)).^2 + r_sct(3)^2)';

    t = 0:1/fs:6/f0;
    t = t-mean(t);
    pulse = exp(sqrt(-1)*2*pi*f0*t).*exp(-(t.^2)*(2*pi*f0/6)^2/2);
    npulse = length(pulse);  

    t_ariv = zeros(N_active,ceil(2*(z_focus+N_active*width)/c*fs));          
    IQ= t_ariv*0; 
    for ix=1:length(delay_z) % Receive
        disp(ix/length(delay_z)*100)
    for jx=1:length(delay_z) % Trans
          d = dist(jx) - delay_T(jx) + dist(ix);
          t_ariv(ix,ceil(d/c*fs)) = t_ariv(ix,ceil(d/c*fs)) + 1/dist(jx)/dist(ix);
    end
    tmp = conv(t_ariv(ix,:),pulse);
    t0=npulse/2;
    te=t0+size(t_ariv,2)-1;
    IQ(ix,:)=tmp(t0:te);
    end
 
    flg_show_rcvData=1;
    if flg_show_rcvData==1
        figure;
        subplot(2,1,1);plot(real(IQ(64,:)));
        subplot(2,1,2);plot(t_ariv(64,:),'.');
        pause();
    end
    
    %% imaging
    %% DAS
    for ix=1:length(x_pxl), 
        [icmp,ix,length(x_pxl)]
        for iz=1:length(z_pxl)
 
        for ielex=1:size(IQ,1) % receive element
            delayT_elei = x_pxl(ix)*tan(theta(icmp)); 
            d_pxl_Rele = sqrt(z_pxl(iz)^2 + (x_pxl(ix)-elex(1,ielex))^2);
            ii=ceil( (d_pxl_Rele + z_pxl(iz) + delayT_elei) *fs/c );
            DAS(ix,iz) = DAS(ix,iz) + IQ(ielex,ii);
    end,end,end
    img=log(abs(DAS)+1e-10);
    img=img-max(max(max(img)));
%   imagesc(z_pxl,x_pxl,img);caxis([-10,0]);
    pause();
end % icmp
log_env=20*log(abs(DAS)+1e-10)/log(10);
log_env=log_env-max(max(max(log_env)));
imagesc(z_pxl,x_pxl,log_env); colormap('gray');caxis([-10,0]);
pause();
% subplot(1,2,1);plot(log_env(:,51)); ylim([-18,0]);
% subplot(1,2,2);plot(log_env(64,:)); ylim([-18,0]);

real_DAS=real(DAS);
imag_DAS=imag(DAS);
close; figure; imagesc(log_env); colormap('gray');caxis([-60,0]);colorbar(); pause();

%str='HFR'
%save(['./',str,'psf_nyx10_PPR', num2str(PPR*6),'_6.mat'], 'y_pxl', 'z_pxl', 'x_pxl', 'log_env', 'DAS','imag_DAS','real_DAS'); 
 
end

 function make_PSF_HFR_2D

 f0=7.5e6;%params_Vera.fc;%6.3e6;                  %  Transducer center frequency [Hz]
 fs=1/3.3333e-8;%1/params_Vera.dt;%100e6;                %  Sampling frequency [Hz]
 c=1480;%params_Vera.C;%1540;                  %  Speed of sound [m/s]
 lambda=c/f0;             %  Wavelength [m]
 width=2e-4;%params_Vera.dele;%lambda;            %  Width of element
% element_height=5/1000;   %  Height of element [m]
 kerf=2e-4;%params_Vera.xele(2)-params_Vera.xele(1);%0.27/1000-lambda;%0.0/1000;           %  Kerf [m]
 % focus=[0 0 30]/1000;     %  Fixed focal point [m]
 z_focus = 30e-3;         % 
% y_focus = 30e-3;        
 N_elements=192;          %  Number of physical elements
 N_active=128;%params_Vera.Nele_R;%128;             %  Number of active elements
 xmit_N_active=128;%params_Vera.Nele_R;%128;       %  Number of active transmit elements for constant F#
 rec_N_active=128;%params_Vera.Nele_R;%128;        %  Number of active receive elements for constant F#
 % Nxsub=1;
 % Nysub=10;
 d_x = (width+kerf);
 dz=lambda/2;
 Ndcm_z=ceil(dz/(c/fs));
 theta=[-18,-12,-6,0,6,12,18]/180*pi;

 fs=fs*5; %% 2020/5/1 without this, tail of PSF does not appeal
 %% DAS
 beamWx=1*lambda; % 0.2 mm
 beamWz=2*lambda/2; % 0.4 mm /2 <
 beamWy=3*lambda; % 0.6 mm
 
PPR=6/6%[3/6,4/6]%:1/6:1%1/6:1/6:1

 % pxl def
 dx_pxl=beamWx*PPR/10;
 dy_pxl=beamWy*PPR/10;
 dz_pxl=beamWz*PPR/10;

x_pxl=-5/1000:dx_pxl:5/1000;
z_pxl=28/1000:dz_pxl:31/1000;
 y_pxl=0:dy_pxl:1/1000;
x_pxl=x_pxl-mean(x_pxl);
 DAS=zeros(length(y_pxl),length(x_pxl),length(z_pxl));
 
 for icmp=1:7
 for iy0=1:length(y_pxl)

 x=width*(1:N_active);
 x=x-mean(x);
 y=element_height/(Nysub-1)*(1:Nysub);
 y=y-mean(y);

 %% transmit delay
 [elex,eley]=meshgrid(x,y);
 %% 1 point focus
 delay_z=sqrt(elex(1,:).^2+ z_focus^2)-z_focus;
 if flg_HFR>0
 %% plane wave
   delay0_z=sqrt(elex(1,:).^2+ z_focus^2)-z_focus;
%   delay_z=elex(1,:)*tan(theta(icmp));
   delay_z=delay_z-max(delay_z-delay0_z);
 end
 delay_y=sqrt(eley(:,1).^2+ y_focus^2)-y_focus;
 delay_T=0*eley';
 for ix=1:N_active, for iy=1:Nysub
         delay_T(ix,iy)=delay_z(ix)+delay_y(iy);
 end,end
  
 %% receive delay < by pixel. set delayY only
  delay_y=sqrt(eley(:,1).^2+ y_focus^2)-y_focus;
  delay_Ry=0*eley';
  for ix=1:N_active, for iy=1:Nysub
         delay_Ry(ix,iy)=delay_y(iy);   % >0 at edge, =0 at center
  end,end
 
 %% dist
 % ptm_pos = [0,0,z_focus];
 r_sct=[0,0,z_focus];
 % r_sct=[0,y_pxl(iy0),z_focus];
 ptm_pos = r_sct;
 ptm_amp = 1;
 dist=sqrt(elex.^2 + (eley-r_sct(2)).^2 + r_sct(3)^2)';

 t=0:1/fs:6/f0;
 t=t-mean(t);
 pulse = exp(sqrt(-1)*2*pi*f0*t).*exp(-(t.^2)*(2*pi*f0/6)^2/2);
 npulse=length(pulse);  

 t_ariv= zeros(N_active,ceil(2*(z_focus+N_active*width)/c*fs));          
 IQ= t_ariv*0;          
 for ix=1:length(delay_z),for iy=1:length(delay_y) % Receive
 for jx=1:length(delay_z),for jy=1:length(delay_y) % Trans
          d = dist(jx,jy) - delay_T(jx,jy) + dist(ix,iy) - delay_Ry(ix,iy);
%% 2020/5/1          t_ariv(ix,ceil(d/c*fs)) = t_ariv(ix,ceil(d/c*fs)) + 1/(dist(jx,jy)+dist(ix,iy));
          t_ariv(ix,ceil(d/c*fs)) = t_ariv(ix,ceil(d/c*fs)) + 1/dist(jx,jy)/dist(ix,iy);
 end,end
 end
 tmp = conv(t_ariv(ix,:),pulse);
 t0=npulse/2;
 te=t0+size(t_ariv,2)-1;
 IQ(ix,:)=tmp(t0:te);
 end
 % figure;
 % subplot(2,1,1);plot(real(IQ(64,:)));
 % subplot(2,1,2);plot(t_ariv(64,:),'.');

 %% DAS
 for ix=1:length(x_pxl), 
         [icmp,iy0,ix,length(y_pxl),length(x_pxl)]
     for iz=1:length(z_pxl)
 
     for ielex=1:size(IQ,1) % receive element
                 delay_pxl=x_pxl*tan(theta(icmp));
%% ??? >-2020/5/1 V 5/1- ii=ceil(( sqrt(z_pxl(iz)^2 + (x_pxl(ix)-elex(1,ielex))^2) + (z_pxl(iz)-delay_pxl(ix)) )*fs/c);
                ii=ceil(( sqrt(z_pxl(iz)^2 + (x_pxl(ix)-elex(1,ielex))^2) + z_pxl(iz) )*fs/c);
             DAS(iy0,ix,iz) = DAS(iy0,ix,iz) + IQ(ielex,ii);
 end,end,end,end
 img=log(abs(DAS)+1e-10);
 img=img-max(max(max(img)));
% imagesc(z_pxl,x_pxl,img);caxis([-10,0]);
 pause(1);
 end % icmp
 log_env=20*log(abs(DAS)+1e-10)/log(10);
 log_env=log_env-max(max(max(log_env)));
%imagesc(z_pxl,x_pxl,log_env);caxis([-10,0]);
 pause(1);
% subplot(1,2,1);plot(log_env(:,51)); ylim([-18,0]);
% subplot(1,2,2);plot(log_env(64,:)); ylim([-18,0]);

[ny,nx,nz]= size(log_env);
tmp=zeros((ny-1)*2+1,nx,nz);
tmp(1:ny,:,:)=log_env(ny:-1:1,:,:);
tmp(ny:2*ny-1,:,:)=log_env;
log_env=tmp;

tmp=zeros((ny-1)*2+1,nx,nz);
tmp(1:ny,:,:)=DAS(ny:-1:1,:,:);
tmp(ny:2*ny-1,:,:)=DAS;
DAS=tmp;

 f_y_pxl=fliplr(y_pxl);
y_pxl=[-f_y_pxl,y_pxl(2:length(y_pxl))];
 

real_DAS=real(DAS);
imag_DAS=imag(DAS);
close; figure; imagesc(log_env); colorbar(); pause();
if flg_HFR==0
   str='TRfocus';
else
   str='HFR'
end
save(['./',str,'psf_nyx10_PPR', num2str(PPR*6),'_6.mat'], 'y_pxl', 'z_pxl', 'x_pxl', 'log_env', 'DAS','imag_DAS','real_DAS'); 

 
 end

 function make_PSF_DAS%(flg_HFR,dx_pxl)

% fname_Set = 'D:\Yamamoto2019\Set';
% [params_Vera, Xpxl Zpxl Pxl_xz dx_pxl dz_pxl] = load_IQData_Setting_VSX(fname_Set);
flg_HFR=1;

 f0=7.5e6;%params_Vera.fc;%6.3e6;                  %  Transducer center frequency [Hz]
 fs=1/3.3333e-8;%1/params_Vera.dt;%100e6;                %  Sampling frequency [Hz]
 c=1480;%params_Vera.C;%1540;                  %  Speed of sound [m/s]
 lambda=c/f0;             %  Wavelength [m]
 width=2e-4;%params_Vera.dele;%lambda;            %  Width of element
 element_height=5/1000;   %  Height of element [m]
 kerf=2e-4;%params_Vera.xele(2)-params_Vera.xele(1);%0.27/1000-lambda;%0.0/1000;           %  Kerf [m]
 theta = [-18,-12,-6,0,6,12,18]/180*pi;
 % focus=[0 0 30]/1000;     %  Fixed focal point [m]
 z_focus = 30e-3;         % 単に画像原点
 y_focus = 30e-3;        
 N_elements=192;          %  Number of physical elements
 N_active=128;%params_Vera.Nele_R;%128;             %  Number of active elements
 xmit_N_active=128;%params_Vera.Nele_R;%128;       %  Number of active transmit elements for constant F#
 rec_N_active=128;%params_Vera.Nele_R;%128;        %  Number of active receive elements for constant F#
 Nxsub=1;
 Nysub=10;
 d_x = (width+kerf);
 dz=lambda/2;
 Ndcm_z=ceil(dz/(c/fs));
 theta=[-18,-12,-6,0,6,12,18]/180*pi;

 fs=fs*5; %% 2020/5/1 without this, tail of PSF does not appeal
 %% DAS
 beamWx=1*lambda; % 0.2 mm hanchi-haba
 beamWz=2*lambda/2; % 0.4 mm /2 <折り返し�?��時間�?SFはλ/2
 beamWy=3*lambda; % 0.6 mm
 
 PPR=2/6;
% for PPR=6/6%[3/6,4/6]%:1/6:1%1/6:1/6:1
 % pxl def
 dx_pxl=beamWx*PPR/10;
%  dy_pxl=beamWy*PPR;
 dy_pxl=beamWy*PPR/10;
 dz_pxl=beamWz*PPR/10;
%x_pxl=-.5/1000:dx_pxl:.5/1000;
z_pxl=29.5/1000:dz_pxl:30.5/1000;
x_pxl=-1/1000:dx_pxl:1/1000;

%x_pxl=-2/1000:dx_pxl:2/1000;
%z_pxl=28/1000:dz_pxl:31/1000;
 y_pxl=0:dy_pxl:1/1000;
x_pxl=x_pxl-mean(x_pxl);
 DAS=zeros(length(x_pxl),length(z_pxl));
 %y_pxl=[f_y_pxl,y_pxl(2:length(y_pxl))];
 
 
 x=width*(1:N_active);
 elex=x-mean(x);
 
 
 %% dist
 % ptm_pos = [0,0,z_focus];
 r_sct=[0,0,z_focus];
 % r_sct=[0,y_pxl(iy0),z_focus];
 ptm_pos = r_sct;
 ptm_amp = 1;
 dist=sqrt(elex.^2 + r_sct(3)^2)';

 t=0:1/fs:6/f0;
 t=t-mean(t);
 pulse = exp(sqrt(-1)*2*pi*f0*t).*exp(-(t.^2)*(2*pi*f0/6)^2/2);
 npulse=length(pulse);  

 for ix_lasta=1:length(x_pxl)
 %% RF
 %% transmit delay
 %% 1 point focus
 delay_T=sqrt((elex-x_pxl(ix_lasta)).^2+ z_focus^2)-z_focus;
 
 t_ariv= zeros(N_active,ceil(2*(z_focus+N_active*width)/c*fs));          
 IQ= t_ariv*0;          
 for ix=1:length(delay_T)%,for iy=1:length(delay_y) % Receive
 for jx=1:length(delay_T)%,for jy=1:length(delay_y) % Trans
          d = dist(jx) - delay_T(jx) + dist(ix);
%% 2020/5/1          t_ariv(ix,ceil(d/c*fs)) = t_ariv(ix,ceil(d/c*fs)) + 1/(dist(jx,jy)+dist(ix,iy));
          t_ariv(ix,ceil(d/c*fs)) = t_ariv(ix,ceil(d/c*fs)) + 1/dist(jx)/dist(ix);
 end%,end
 %end
 tmp = conv(t_ariv(ix,:),pulse);
 t0=npulse/2;
 te=t0+size(t_ariv,2)-1;
 IQ(ix,:)=tmp(ceil(t0):ceil(te));
 end
% close; figure;
%  subplot(2,1,1);imagesc(t_ariv);%plot(real(IQ(64,:)));
%  subplot(2,1,2);imagesc(real(IQ));

 %% DAS
         [ix_lasta,length(x_pxl)]
     for iz=1:length(z_pxl)
 
     for ielex=1:size(IQ,1) % receive element
         ii=ceil(( sqrt(z_pxl(iz)^2 + (x_pxl(ix_lasta)-elex(1,ielex))^2) + z_pxl(iz) )*fs/c);
         DAS(ix_lasta,iz) = DAS(ix_lasta,iz) + IQ(ielex,ii);
 end,end
 end % ix_lasta
 img=log(abs(DAS)+1e-10);
 img=img-max(max(max(img)));
% imagesc(z_pxl,x_pxl,img);caxis([-10,0]);
 %pause(1);
 
log_env=20*log(abs(DAS)+1e-10)/log(10);
 log_env=log_env-max(max(max(log_env)));
%imagesc(z_pxl,x_pxl,log_env);caxis([-10,0]);
 pause(1);
% subplot(1,2,1);plot(log_env(:,51)); ylim([-18,0]);
% subplot(1,2,2);plot(log_env(64,:)); ylim([-18,0]);
 
real_DAS=real(DAS);
imag_DAS=imag(DAS);
close; figure; imagesc(z_pxl,x_pxl,log_env); caxis([-60,0]);colorbar(); pause();
   str='TRfocus';
save(['./',str,'psf_nyx10_PPR', num2str(PPR*6),'_6.mat'], 'y_pxl', 'z_pxl', 'x_pxl', 'log_env', 'DAS','imag_DAS','real_DAS'); 


 end

 
 