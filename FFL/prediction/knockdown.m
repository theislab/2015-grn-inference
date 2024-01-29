
clear all
close all
clc

% mRNA knockdown using siRNA to achieve a 3-fold increase in mRNA degradation
% through parameter KK [we do not change the decay for gene gX, so that input
% gX is still the same]

KK = 3;

dt = 1;
T  = 340; 

sigma = 5;

decay           = 0.25;
synth           = 250;
K_d_inhibition  = 200;
K_d_activation  = 400;
n_inhibition    = 2;
n_activation    = 20;
synthS          = 100;
K_d_inhibitionS = 200;
n_inhibitionS   = 10;

init = [315 285 400 0]';

a(:,1) = init;

nRealiz = 1000;

meanTime = zeros(nRealiz,20);

for j = 1:nRealiz
    
    kk = 1;
    
    for i = 1:T-1

        a(1,i+1) = a(1,i) + dt*(synth*(1/(1+ (a(2,i)/K_d_inhibition)^n_inhibition)) - decay*a(1,i)) + sqrt(sigma)*2*randn;
        a(2,i+1) = a(2,i) + dt*(synth*(1/(1+ (a(1,i)/K_d_inhibition)^n_inhibition)) - decay*a(2,i)) + sqrt(sigma)*2*randn;
        
        a(3,i+1) = a(3,i) + dt*(synthS*(1/(1+(K_d_activation/a(1,i))^n_activation)) + synthS*(1/(1+ (a(4,i)/K_d_inhibitionS)^n_inhibitionS)) - KK*decay*a(3,i)) + sqrt(sigma)*randn;
        a(4,i+1) = a(4,i) + dt*(synthS*(1/(1+(K_d_activation/a(1,i))^n_activation))  - KK*decay*a(4,i)) + sqrt(sigma)*randn;
        
        if (sum(a([1],i+1)>960))
            a(:,i+1) = init;
            meanTime(j,kk) = i;
            kk = kk+1;
        end
    
    end

   S(j,1:T,1:4) = a';
   fprintf('Realisation %d of %d\n',j,nRealiz)

end

mT = round(mean(meanTime(:,1)));

for i = 1:T-1

    a(1,i+1) = a(1,i) + dt*(synth*(1/(1+ (a(2,i)/K_d_inhibition)^n_inhibition)) - decay*a(1,i));
    a(2,i+1) = a(2,i) + dt*(synth*(1/(1+ (a(1,i)/K_d_inhibition)^n_inhibition)) - decay*a(2,i));
        
    a(3,i+1) = a(3,i) + dt*(synthS*(1/(1+(K_d_activation/a(1,i))^n_activation)) + synthS*(1/(1+ (a(4,i)/K_d_inhibitionS)^n_inhibitionS)) - KK*decay*a(3,i));
    a(4,i+1) = a(4,i) + dt*(synthS*(1/(1+(K_d_activation/a(1,i))^n_activation))  - KK*decay*a(4,i));

end

varData = 20;

% Noisy data generated from perturbed model
m(1,:) = a(1,1:mT) + varData*randn(1,mT);
m(2,:) = a(3,1:mT) + varData*randn(1,mT);
m(3,:) = a(4,1:mT) + varData*randn(1,mT);

% -------------------------------------------------------------------------

% Simulation generated from inferred network with estimated parameters,
% using same perturbation. 

% Input to the model: we use GP emulator for gX
cd GPs
[me,st] = GPreg(m(1,:)',1:mT,1,1,mT);
cd ..

paramY = [398.5130    0.2489  420.0862   14.7256];
paramZ = [389.4401    0.2490  417.1350   14.1570  264.2697    8.1481];

synthY           = paramY(1)*paramY(2);
decayY           = paramY(2);
K_d_activationY  = paramY(3);
n_activationY    = paramY(4);

synthZ           = paramZ(1)*paramZ(2);
decayZ           = paramZ(2);
K_d_activationZ  = paramZ(3);
n_activationZ    = paramZ(4);

K_d_inhibitionZ  = paramZ(5);
n_inhibitionZ    = paramZ(6);

clear a
a(:,1) = init;

for i = 1:mT-1

    a(1,i+1) = me(i+1);
        
    a(3,i+1) = a(3,i) + dt*(synthZ*(1/(1+(K_d_activationZ/a(1,i))^n_activationZ)) + synthZ*(1/(1+ (a(4,i)/K_d_inhibitionZ)^n_inhibitionZ)) - KK*decayZ*a(3,i));
    a(4,i+1) = a(4,i) + dt*(synthY*(1/(1+(K_d_activationY/a(1,i))^n_activationY))  - KK*decayY*a(4,i));

end

plot(m(2,:),'.'), hold on, plot(a(3,:),'r')
