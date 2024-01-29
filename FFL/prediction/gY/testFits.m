
clear all
close all
clc

gene    = 3;   % gene to simulate
emuInd1 = 1;   % emulator 1

Data = importdata('DataBranch_1_3.txt');
Time = importdata('TimesBranch_1_3.txt');
Emul = importdata('Interpolators1.txt');

% - - - - - - - - - - - - - - - - - - - - - - - - -
estPara  = importdata('Risulta_100_20');
[a,ind]  = sort(estPara(:,size(estPara,2)));
paramAll = estPara(ind,1:size(estPara,2)-1);
param    = exp(estPara(ind(end),1:size(estPara,2)-1));
% - - - - - - - - - - - - - - - - - - - - - - - - -

L  = size(Emul,1);
dt = 1;

nRep = 3; % replicates
nG   = 3; % number of genes

x = zeros(nRep,L);

emu1 = Emul(:,emuInd1)';

x(:,1) = Emul(1,gene);

for rep = 1:nRep
    for i = 2:L        
        k1 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i-1))^(param(4))))) - param(2)*x(rep,i-1);
        k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4))))) - param(2)*(x(rep,i-1)+0.5*k1);
        k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4))))) - param(2)*(x(rep,i-1)+0.5*k2);
        k4 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i))^(param(4))))) - param(2)*(x(rep,i-1)+k3);        
        x(rep,i) = x(rep,i-1) + (k1+2*k2+2*k3+k4)/6;
    end
end

for i = 1:nRep
    subplot(2,2,i), plot(x(i,:)), hold on, plot(Time,Data(:,i),'r.','MarkerSize',10,'LineWidth',3)
    xlabel('time','FontSize',14), ylabel('gene expression','FontSize',14)
    axis([0 length(x) -2 1000])
end


trueP = [100/0.25 0.25 400 20];

fprintf('\n\n Average relative error: %.2f %%\n\n',mean(abs((param-trueP)./trueP))*100)

