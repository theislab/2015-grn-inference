
clear, home

addpath('RT')

% - - - - - - - - - - - - - - - -  FFL  - - - - - - - - - - - - - - - - - -
load('ORgateFFL.mat')
Data = DataDiff';

VIMA2   = genie3(Data,[1 2 3],'ET');

VIMA2.*(VIMA2>0.3)
% - - - - - - - - - - - - - - - Toggle switch - - - - - - - - - - - - - - -
load('toggle.mat')
Data = DataDiff';

VIMB2   = genie3(Data,[1 2 3 4 5 6],'ET');

VIMB2.*(VIMB2>0.3)
% - - - - - - - - - - - - - - - Hematopoietic Net - - - - - - - - - - - - -
load('genenames_vicki18.mat'); names = genenames18; clear genenames18
Data = importdata('vicki_normalized.txt');

VIMC2   = genie3(Data,[1:18],'ET');

VIMC2.*(VIMC2>0.08)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -