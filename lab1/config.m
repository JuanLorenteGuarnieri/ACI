% name=['bunny_d=0.5_c=[256x256]',
%       'bunny_d=0.5_l=[1x1]_s=[256x256]',
%       'bunny_d=0.5_l=[16x16]_s=[16x16]',
%       'bunnybox_d=0.5_l=[16x16]_s=[16x16]',
%       'planes_d=0.5_l=[16x16]_s=[16x16]',
%       'Z_d=0.5_l=[1x1]_s=[256x256]'];
%       'usaf_d=0.5_l=[1x1]_s=[256x256]';
name = ['usaf_d=0.5_l=[1x1]_s=[256x256]'];             % Dataset name

signalDegradation = false;                        % Signal degradation option (true/false)
photon_counts = 1e6;                             % Photon counts for Poisson noise (in photons) [1e4, 1e5, 1e6, 1e7, 1e8]
focal_plane = 0.5;                               % Unique focal plane distance (in meters) for multipleFocalPlanes == false

multipleFocalPlanes = true;                     % Multiple focal planes option (true/false)
num_planes = 20; %200
minPlaneDistance = 0.3; %0.3                     % Minimum distance between planes (in meters)
maxPlaneDistance = 0.7; %2.5                     % Maximum distance between planes (in meters)

isMultiFrequency = false;
isNarrowBand = false;                            % Narrowband option (true/false)

dataset_path = strcat('data/', name, '.mat');           % Path to the dataset