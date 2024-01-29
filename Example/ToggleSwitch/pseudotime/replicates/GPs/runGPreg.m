
clear all
clc

cd ..
load('KK')
load('T')
cd GPs

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ngenes = 6;
dt     = 1;
oStep  = 1;

for i = 1:ngenes
    
    xaxis = KK.asseState1-1;
    ydata = KK.dataState1(i,:)';
    
    [me,st] = GPreg(ydata,xaxis,dt,oStep,T);
    y1(i,:) = me;
    
    assignin('base',strcat('Branch_',int2str(1),'_',int2str(i)),[xaxis;ydata'])
    
end
for i = 1:ngenes

    xaxis = KK.asseState2-1;
    ydata = KK.dataState2(i,:)';
    
    [me,st] = GPreg(ydata,xaxis,dt,oStep,T);
    y2(i,:) = me;
    
    assignin('base',strcat('Branch_',int2str(2),'_',int2str(i)),[xaxis;ydata'])
    
end

M1 = y1';
M2 = y2';

M1(find(M1<0.01)) = 0.01;
M2(find(M2<0.01)) = 0.01;

figure
plot(KK.asseState1,KK.dataState1,'b.'); hold on, plot(1:T,y1,'r')
figure
plot(KK.asseState2,KK.dataState2,'b.'); hold on, plot(1:T,y2,'r')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Export emulators

cd ..
filename1 = cat(2,'Interpolators1.txt');
filename2 = cat(2,'Interpolators2.txt');
%dlmwrite(filename1,M1,'\t');
%dlmwrite(filename2,M2,'\t');
cd('GPs')



