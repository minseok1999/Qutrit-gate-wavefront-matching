
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

clear all;
clc;
close all;
close all hidden;

matfile = 'unitary_trial';	
Hologram=load(fullfile('C:\Users\rmins\Desktop\MPLC\MPLC_system\WaveFrontMatching_Code\trials',matfile));	
Hologram=Hologram.Hologram;
tic;
global lamda k
lamda=808e-3;
k=2*pi/lamda;
%sampling point interval in micrometer=pixel0%
pixel0=8;
%phase modulation%
%%%%%%	
phi1=Hologram(:,:,1);
%%%%%%
toc

dim=size(phi1);
mx0=dim(1);
my0=dim(2);
%%%%% L0(mm)=diameter=8micrometerxdim(1)(mm)%%
L0=(dim(1))*pixel0;
%%%%everything written in units of micrometer%
[xx,yy]=meshgrid(-(mx0-1)/2:(mx0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;
%No aperture%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%

waist=0.9e3;

toc

    % Grid in cylindrical coordinates:
    Rad = sqrt(x0.^2+y0.^2);
    Angle = angle(x0+1i.*y0)+pi; % Matrix with all angles starting left-center

g=(GenModesLG([1 0], waist, Rad, Angle)+GenModesLG([-1 0], waist, Rad, Angle))./sqrt(2);

%%%%%%%%%

g=g.*exp(1i.*phi1);
%Input plot
%Input defining
h=abs(g).^2; h=h./max(h(:)); %%%
%unit in mm in plotting, not micrometer%%
fig(1)=figure('Name','Input Gaussian amplitude');
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g);%
fig(2)=figure('Name','Modulated phase');
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%Input plot done%%%

tic;
d=800e3;                                 
% L=lamda*d/pixel0;                       
L1=L0;
x1start=-L1./2;
x1end=L1./2;
y1start=-L1./2;
y1end=L1./2; 
mx1=625;
my1=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);   
%%%%%%%
g1=g1.*exp(1i.*Hologram(:,:,2));
%%%%%
L2=L0;
x2start=-L2./2;
x2end=L2./2;
y2start=-L2./2;
y2end=L2./2; 
mx2=625;
my2=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g2,pixel2]=Scalar_Bluestein(g1,mx1,my1,pixel1,d,x2start,x2end,y2start,y2end,mx2,my2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=linspace(x2start,x2end,mx2);                
y2=linspace(y2start,y2end,my2);                
[x2,y2]=meshgrid(x2,y2);     

%%%%%%%
g2=g2.*exp(1i.*Hologram(:,:,3));

%%%%%
L3=L0;
x3start=-L3./2;
x3end=L3./2;
y3start=-L3./2;
y3end=L3./2; 
mx3=625;
my3=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g3,pixel3]=Scalar_Bluestein(g2,mx2,my2,pixel2,d,x3start,x3end,y3start,y3end,mx3,my3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3=linspace(x3start,x3end,mx3);                
y3=linspace(y3start,y3end,my3);                
[x3,y3]=meshgrid(x3,y3);     

%%%%%%%
g3=g3.*exp(1i.*Hologram(:,:,4));


%%%%%
L4=L0;
x4start=-L4./2;
x4end=L4./2;
y4start=-L4./2;
y4end=L4./2; 
mx4=625;
my4=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g4,pixel4]=Scalar_Bluestein(g3,mx3,my3,pixel3,d,x4start,x4end,y4start,y4end,mx4,my4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x4=linspace(x4start,x4end,mx4);                
y4=linspace(y4start,y4end,my4);                
[x4,y4]=meshgrid(x4,y4);     

%%%%%%%
g4=g4.*exp(1i.*Hologram(:,:,5));




%%%%%%%

L5=L1;
x5start=-L5./2;
x5end=L5./2;
y5start=-L5./2;
y5end=L5./2; 
mx5=625;
my5=625;
d=d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g5,pixel5]=Scalar_Bluestein(g4,mx4,my4,pixel4,d,x5start,x5end,y5start,y5end,mx5,my5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x5=linspace(x5start,x5end,mx5);                
y5=linspace(y5start,y5end,my5);                
[x5,y5]=meshgrid(x5,y5);    

%%%mode overlap with desired beam

 % Grid in cylindrical coordinates:
    Rad = sqrt(x5.^2+y5.^2);
    Angle = angle(x5+1i.*y5)+pi; % Matrix with all angles starting left-center
    
gdot=((i).*GenModesLG([1 0], waist, Rad, Angle)+(1+2i).*GenModesLG([-1 0], waist, Rad, Angle)+(1+i).*GenModesLG([0 0], waist, Rad, Angle))./sqrt(8);
overlap=(abs(sum(g5.*conj(gdot),'all'))^2)./(sum(abs(g5).^2,'all'))./(sum(abs(gdot).^2,'all'))

%check
h=abs(g5).^2; h=h./max(h(:)); %%%
%%%%%%
fig(3)=figure('Name','far field amplitude');
surf(x5./1e3,y5./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g5);%
fig(4)=figure('Name','far field phase');
surf(x5./1e3,y5./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%



h=abs(gdot).^2; h=h./max(h(:)); %%%
%%%%%%
fig(5)=figure('Name','desired amplitude');
surf(x5./1e3,y5./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(gdot);%
fig(6)=figure('Name','desired phase');
surf(x5./1e3,y5./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%

for k=1:6
  
 saveas(fig(k), fullfile('C:\Users\rmins\Desktop\hologram',sprintf('1-1superposition_5planes_%d',k)), 'jpeg')
  
end