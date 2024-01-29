
clear all
close all
clc

load('BMApprxCensored_sigma0.9.mat')
psi = psi';
D   = labels;

ind1 = find(psi(3,:)<-0.2); L1 = length(ind1);

psi = psi + 1e5;
psi(:,ind1) = 0;
L = length(nonzeros(psi))/6;
psi = reshape(nonzeros(psi),6,L) - 1e5;

D = D + 1e5;
D(ind1) = 0;
D = nonzeros(D) - 1e5;

scatter3(psi(1,:),psi(2,:),psi(3,:),50,D,'fill')
set(gca,'FontSize',27,'FontName','Courier')
xlabel('DC_1'), ylabel('DC_2'), zlabel('DC_3')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
axis tight
grid off