
clear all
close all
clc

load diffusionInput

data    = DataDiff';
no_dims = 4;
nn      = size(data,1);
t       = 1;
sigma   = 1000;

[psi,E] = diffusion_maps_nn(data, no_dims, nn, t, sigma);

% save psi psi

figure(1)
for i = 1:6
    subplot(2,3,i), scatter3(psi(:,1),psi(:,2),psi(:,4), 50, data(:,i),'fill')
end