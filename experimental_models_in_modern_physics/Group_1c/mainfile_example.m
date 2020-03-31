clc
close all
clear 

file1 = 'Pd_diskrod_200nm_15pt_pulse';

lampfile = 'CRS_700nm_glas';

nbr_particles = 15;

[lambda,particles,background,lamp] = getData(file1,lampfile,nbr_particles);

data = (particles - background)./lamp;

%%
plotfcn = data(:,1,1);
%%
figure(1)
plot(lambda, plotfcn, 'b')