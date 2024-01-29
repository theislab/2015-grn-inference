
clear all
close all
clc

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Define length of vectors for GP regression

load(strcat('Branch_1_rep_',int2str(1)))
L = length(replicate);
load(strcat('Branch_2_rep_',int2str(1)))
L = [L length(replicate)];

T = max(L);
if mod(T,2)
    T = T+1;
end
save T T

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

load(strcat('Branch_1_rep_',int2str(1)))
A = replicate;
Br1r1 = A(2:end,A(1,:))';

load(strcat('Branch_1_rep_',int2str(2)))
A = replicate;
Br1r2 = A(2:end,A(1,:))';

load(strcat('Branch_1_rep_',int2str(3)))
A = replicate;
Br1r3 = A(2:end,A(1,:))';

L = length(Br1r1);

for i = 1:6
    KK.dataState1(i,:) = [Br1r1(:,i)' Br1r2(:,i)' Br1r3(:,i)'];    
end

KK.asseState1 = round(([1:L 1:L 1:L]/L)*T);
KK.asseState1 = KK.asseState1(:,end:-1:1);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

load(strcat('Branch_2_rep_',int2str(1)))
A = replicate;
Br2r1 = A(2:end,A(1,:))';

load(strcat('Branch_2_rep_',int2str(2)))
A = replicate;
Br2r2 = A(2:end,A(1,:))';

load(strcat('Branch_2_rep_',int2str(3)))
A = replicate;
Br2r3 = A(2:end,A(1,:))';

L = length(Br2r1);

for i = 1:6
    KK.dataState2(i,:) = [Br2r1(:,i)' Br2r2(:,i)' Br2r3(:,i)'];    
end

KK.asseState2 = round(([1:L 1:L 1:L]/L)*T);
KK.asseState2 = KK.asseState2(:,end:-1:1);

save KK KK

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Create data replicates

L1 = length(KK.dataState1)/3;
L2 = length(KK.dataState2)/3;

State = 1;
for Gene = 1:6
    Data = reshape(KK.dataState1(Gene,:),L1,3);
    filename = strcat('DataBranch_',int2str(State),'_',int2str(Gene),'.txt');
    %dlmwrite(filename,Data,'\t');
    times = KK.asseState1(1:L1);
    filename = strcat('TimesBranch_',int2str(State),'_',int2str(Gene),'.txt');
    %dlmwrite(filename,times,'\t');
end

State = 2;
for Gene = 1:6
    Data = reshape(KK.dataState2(Gene,:),L2,3);
    filename = strcat('DataBranch_',int2str(State),'_',int2str(Gene),'.txt');
    %dlmwrite(filename,Data,'\t');
    times = KK.asseState2(1:L2);
    filename = strcat('TimesBranch_',int2str(State),'_',int2str(Gene),'.txt');
    %dlmwrite(filename,times,'\t');
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

figure
for i = 1:6
    subplot(2,3,i), plot(Br1r1(:,i),'g'), hold on, plot(Br1r2(:,i),'r'), hold on, plot(Br1r3(:,i),'b')
    axis([0 T 0 1000])
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% GP emulators

cd GPs
runGPreg
