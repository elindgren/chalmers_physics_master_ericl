
%load(strcat(folder,'p.mat'))
%load(strcat(folder,'b.mat'))
%load(strcat(folder,'lamp.mat'))
%load(strcat(folder,'lambda.mat'))

comma = 1;          % if need to covert from comma to dot decimal point

if comma == 1
     comma2dot('Au_100nm');
end
%% 
close all

clear
nbr_part = 3;

% meas = 'p04';
p01 = importdata('Au_100nm.asc');
% p01 = importdata('test_getdata.asc');
 
p = zeros(size(p01,1),nbr_part);
b = zeros(size(p01,1),1);


lambda = p01(:,1);
data = p01(:,2:end);

for i=1:nbr_part
    p(:,i) = data(:,i+1);     
   
end
b = data(:,1);

p_reshaped=reshape(p,1024,size(p,1)/1024,nbr_part);
b_reshaped = reshape(b,1024,size(p,1)/1024,1);



lamp = importdata('CRS_lamp_spectrum.asc');
lambda = lamp(:,1);
lamp(:,1) = [];

save('lamp.mat','lamp')
save('lambda.mat','lambda')


Spectra = (p_reshaped-b_reshaped) ./ (lamp);

figure
plot(lambda,Spectra(:,[1:30:end],1))

%%
% peak_guess = [720 720 700 700 700];
peak_guess = [720 720 780 780 780];
for i=1:nbr_part
    [peak{i},C{i},ex{i},FWHM{i}] = PeakFitLorentzCentroid(Spectra(:,:,i),lambda',750,50,50,1);
end
%% testing getdata
% close all
nbr_part=15
clear Spectra

% file = '190912_Pd120_H2_pulse';
file = 'Pd_diskrod_200nm_15pt_pulse';
% lampfile = '190911_lamp700_slit50';
lampfile = 'CRS_700nm_glas';
[lambda,p_reshaped,b_reshaped,lamp] = getData(file,lampfile,15);

Spectra = (p_reshaped-b_reshaped) ./ (lamp);

figure(8)
% plot(lambda,Spectra(:,1:5:end,1))
surf(Spectra(:,:,14))
shading flat
view(1)

%%
% peak_guess = [720 720 700 700 700];
peak_guess = [720 720 780 780 780];
for i=1:nbr_part
    [peak{i},C{i},ex{i},FWHM{i}] = PeakFitLorentzCentroid(100*Spectra(:,:,i),lambda',680,50,50,1);
end


%% 
nbr_part = 15
figure
ylabel('Peak position [nm]')
xlabel('Time [min]')
hold all
for i=1:nbr_part
    plot((10:10:10*length(ex{i}))./60,peak{i})
end
legend('P1','P2','P3','P4','P5')
legend boxoff

figure
ylabel('Centriod [nm]')
xlabel('Time [min]')
hold all
for i=1:nbr_part
    plot((10:10:10*length(ex{i}))./60,C{i})
end
legend('P1','P2','P3','P4','P5')
legend boxoff

figure
ylabel('FWHM [nm]')
xlabel('Time [min]')
hold all
for i=1:nbr_part
    plot((10:10:10*length(ex{i}))./60,FWHM{i})
end
legend('P1','P2','P3','P4','P5')
legend boxoff

figure
hold all
ylabel('Scattering intensity [a.u.]')
xlabel('Time [min]')
for i=1:nbr_part
    plot((10:10:10*length(ex{i}))./60,ex{i})
end
legend('P1','P2','P3','P4','P5')
legend boxoff
