function [lambda,particles,background,lamp] = getData(file,lampfile,nbr_particles)
% function [lambda,particle,background] = getdata(file,lampfile,nbr_part)
% Provide the name of the ASCII file (without extention) as variable
% 'file', the lamp reference as 'lampfile', as well as the number of 
% particles measured on. The dark background reference is
% assumed to be below the first particle.
%
% Output: wavelength vector 'lambda', scattering intenisty collected from
% particles at all time steps in 'particles', and the background scattering
% at all time steps as 'background'. 'Lamp' is the lamp reference.

comma = 0;          
if comma == 1
    comma2dot(file);
    comma2dot(lampfile);
end

imported_data = importdata(strcat(file,'.asc'));
nbr_pixels = 1024;

lambda = imported_data(1:nbr_pixels,1);
background_data = imported_data(:,2);
particle_data = imported_data(:,3:nbr_particles+2);

lamp = importdata(strcat(lampfile,'.asc'));
lamp = lamp(1:end,2);

timesteps = (size(particle_data,1)/nbr_pixels);

background = reshape(background_data,nbr_pixels,timesteps); 
particles = reshape(particle_data,nbr_pixels,timesteps,nbr_particles);             

if timesteps == 1
    particles = permute(particles,[1 3 2]);
end
end