
clear all
close all
clc

target = 'geneC';
folder = {'ANDpm';'ANDpp';'ORpm';'ORpp';'Ap';'Dm';'Dp'};

gene      = 3; % gene to simulate
emuInd1   = 1; % emulator 1
emuInd2   = 4; % emulator 2

nBranches = 2; % number of branches
nParams   = [6 6 6 6 4 4 4];

for br = 1:2

    cd(target)
    if br==1
        Data = importdata('DataBranch_1_3.txt');
        Time = importdata('TimesBranch_1_3.txt');
        Emul = importdata('Interpolators1.txt');
    end
    if br==2
        Data = importdata('DataBranch_2_3.txt');
        Time = importdata('TimesBranch_2_3.txt');
        Emul = importdata('Interpolators2.txt');
    end
    cd ..

    L  = size(Emul,1);
    dt = 1;

    nRep    = 3;              % replicates
    nG      = 6;              % number of genes
    nModels = length(folder); % number of models

    x = zeros(1,L);

    emu1 = Emul(:,emuInd1)';
    emu2 = Emul(:,emuInd2)';

    x(:,1) = Emul(1,gene);

    for rep = 1:nModels
    
        % - - - - - - - - - - - - - - - - - - - - - - - - -
        cd(target)
        cd(folder{rep})
        estPara  = importdata('Risulta_100_20');
        [a,ind]  = sort(estPara(:,size(estPara,2)));
        paramAll = estPara(ind,1:size(estPara,2)-1);
        param    = exp(estPara(ind(end),1:size(estPara,2)-1));
        cd ..
        cd ..
        % - - - - - - - - - - - - - - - - - - - - - - - - -    
    
        for i = 2:L
        
            switch folder{rep}
                case 'ANDpm'
                    k1 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i-1))^(param(4)))) * (1/(1+(param(5)/emu2(i-1))^(-param(6))))) - param(2)*x(1,i-1);
                    k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) * (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(-param(6))))) - param(2)*(x(1,i-1)+0.5*k1);
                    k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) * (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(-param(6))))) - param(2)*(x(1,i-1)+0.5*k2);
                    k4 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i))^(param(4)))) * (1/(1+(param(5)/emu2(i))^(-param(6))))) - param(2)*(x(1,i-1)+k3);        
                    x(1,i) = x(1,i-1) + (k1+2*k2+2*k3+k4)/6;
                case 'ANDpp'
                    k1 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i-1))^(param(4)))) * (1/(1+(param(5)/emu2(i-1))^(param(6))))) - param(2)*x(1,i-1);
                    k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) * (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(param(6))))) - param(2)*(x(1,i-1)+0.5*k1);
                    k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) * (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(param(6))))) - param(2)*(x(1,i-1)+0.5*k2);
                    k4 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i))^(param(4)))) * (1/(1+(param(5)/emu2(i))^(param(6))))) - param(2)*(x(1,i-1)+k3);        
                    x(1,i) = x(1,i-1) + (k1+2*k2+2*k3+k4)/6;
                case 'ORpm'
                    k1 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i-1))^(param(4)))) + (1/(1+(param(5)/emu2(i-1))^(-param(6))))) - param(2)*x(1,i-1);
                    k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) + (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(-param(6))))) - param(2)*(x(1,i-1)+0.5*k1);
                    k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) + (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(-param(6))))) - param(2)*(x(1,i-1)+0.5*k2);
                    k4 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i))^(param(4)))) + (1/(1+(param(5)/emu2(i))^(-param(6))))) - param(2)*(x(1,i-1)+k3);        
                    x(1,i) = x(1,i-1) + (k1+2*k2+2*k3+k4)/6;
                case 'ORpp'
                    k1 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i-1))^(param(4)))) + (1/(1+(param(5)/emu2(i-1))^(param(6))))) - param(2)*x(1,i-1);
                    k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) + (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(param(6))))) - param(2)*(x(1,i-1)+0.5*k1);
                    k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4)))) + (1/(1+(param(5)/(0.5*(emu2(i-1)+emu2(i))))^(param(6))))) - param(2)*(x(1,i-1)+0.5*k2);
                    k4 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i))^(param(4)))) + (1/(1+(param(5)/emu2(i))^(param(6))))) - param(2)*(x(1,i-1)+k3);        
                    x(1,i) = x(1,i-1) + (k1+2*k2+2*k3+k4)/6;
                case 'Ap'
                    k1 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i-1))^(param(4))))) - param(2)*x(1,i-1);
                    k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4))))) - param(2)*(x(1,i-1)+0.5*k1);
                    k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu1(i-1)+emu1(i))))^(param(4))))) - param(2)*(x(1,i-1)+0.5*k2);
                    k4 = param(1) * param(2) * ((1/(1+(param(3)/emu1(i))^(param(4))))) - param(2)*(x(1,i-1)+k3);        
                    x(1,i) = x(1,i-1) + (k1+2*k2+2*k3+k4)/6;
                case 'Dm'
                    k1 = param(1) * param(2) * ((1/(1+(param(3)/emu2(i-1))^(-param(4))))) - param(2)*x(1,i-1);
                    k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu2(i-1)+emu2(i))))^(-param(4))))) - param(2)*(x(1,i-1)+0.5*k1);
                    k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu2(i-1)+emu2(i))))^(-param(4))))) - param(2)*(x(1,i-1)+0.5*k2);
                    k4 = param(1) * param(2) * ((1/(1+(param(3)/emu2(i))^(-param(4))))) - param(2)*(x(1,i-1)+k3);        
                    x(1,i) = x(1,i-1) + (k1+2*k2+2*k3+k4)/6;
                case 'Dp'
                    k1 = param(1) * param(2) * ((1/(1+(param(3)/emu2(i-1))^(param(4))))) - param(2)*x(1,i-1);
                    k2 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu2(i-1)+emu2(i))))^(param(4))))) - param(2)*(x(1,i-1)+0.5*k1);
                    k3 = param(1) * param(2) * ((1/(1+(param(3)/(0.5*(emu2(i-1)+emu2(i))))^(param(4))))) - param(2)*(x(1,i-1)+0.5*k2);
                    k4 = param(1) * param(2) * ((1/(1+(param(3)/emu2(i))^(param(4))))) - param(2)*(x(1,i-1)+k3);        
                    x(1,i) = x(1,i-1) + (k1+2*k2+2*k3+k4)/6;
            end
        
        end
        
        sigma = 20;
    
        for i = 1:nRep
            LL(rep,i) = loglik(x(1,Time),Data(:,i),sigma);
        end
        
        subplot(2,7,rep+(nModels*(br-1))), plot(x(1,:)), hold on, plot(Time,Data(:,1),'r.','MarkerSize',10,'LineWidth',3)
        xlabel('time','FontSize',14), ylabel('gene expression','FontSize',14)
        axis([0 length(x) -2 1000])
    
    end
    
    AIC(:,br) = mean(LL,2) - nParams';
    BIC(:,br) = mean(LL,2) - 0.5*(nParams')*log(length(Data));
    
end

% Mean of AIC and BIC among branches
AIC = mean(AIC,2);
BIC = mean(BIC,2);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

cd(target)
cd('ANDpm')
estPara  = importdata('Risulta_100_20');
[a,ind]  = sort(estPara(:,size(estPara,2)));
paramAll = estPara(ind,1:size(estPara,2)-1);
param    = exp(estPara(ind(end),1:size(estPara,2)-1));
cd ..
cd ..

trueP = [250/0.25 0.25 400 20 200 2];

fprintf('\n\n Average relative error: %.2f %%\n\n',mean(abs((param-trueP)./trueP))*100)
