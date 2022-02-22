clear; clc; close all; 

f = fopen("~/Documents/Modelling/Modelling Projects Git/FD2D/build/OUTPUT/field_P_00001","r");
u = fread(f,"float64");
fclose(f);

nz = 52;
nx = 232;
npml = 60;

nx = nx + 2*npml;
nz = nz + npml;

%u = reshape(u,400,nz,[]); 
u = reshape(u,[],nz,nx);


maxamp = max(max(max(u)));       % max amplitude calculation for clip
clip = 98;  % clip value
ampclip = (1-clip/100)*maxamp;
data    = u; 
data(data > ampclip)= ampclip;   % clipping positive
data(data < -ampclip)= -ampclip; % clipping negative


for is=1:200
    imagesc(reshape(data(is,:,:),nz,nx));   
    colormap(redblue);
    %caxis([-6e-15 6e-15])
    caxis([-max(data(50,:,:),[],'all')  max(data(50,:,:),[],'all') ])

    colorbar

    pause(.009) 
    
end
