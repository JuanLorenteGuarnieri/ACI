% clear; clc; close all;

% To configure the parameters of the execution, edit the file config.m
run('config.m');

function f=fftfreq(npts,dt,alias_dt)
% returns a vector of the frequencies corresponding to the length
% of the signal and the time step.
% specifying alias_dt > dt returns the frequencies that would
% result from subsampling the raw signal at alias_dt rather than
% dt.
if (nargin < 3)
  alias_dt = dt;
end
fmin = -1/(2*dt);
df = 1/(npts*dt);
f0 = -fmin;
alias_fmin = -1/(2*alias_dt);
f0a = -alias_fmin;

ff = mod(linspace(0, 2*f0-df, npts)+f0,  2*f0)  - f0;
fa = mod(                        ff+f0a, 2*f0a) - f0a;
%  return the aliased frequencies
f = fa;
end

%% Load the dataset
dataset = load(dataset_path).data;

% Determine if the dataset is confocal
if ndims(dataset.data) == 3
  isConfocal = true;
  fprintf('The dataset is confocal.\n');
else
  isConfocal = false;
  fprintf('The dataset is non-confocal.\n');
end

%% Extract measurement parameters and positions
t0 = dataset.t0;          % Temporal offset (in meters, optical distance)
deltaT = dataset.deltaT;  % Temporal resolution (in meters)
H = dataset.data;         % Measurement data H

if ~isConfocal
  % For non-confocal data, use separate grids.
  laserPos = dataset.laserPositions;  % MxNx3 array
  spadPos  = dataset.spadPositions;   % PxQx3 array
  origin_laser = dataset.laserOrigin;
  origin_spad  = dataset.spadOrigin;

  [laserRows, laserCols, ~] = size(laserPos);
  [spadRows, spadCols, ~] = size(spadPos);

end

%% 2.3 Signal degradation

if signalDegradation && ~isConfocal
  H_total = sum(H(:));
  H_norm = H / H_total;
  H_poisson = poissrnd(photon_counts * H_norm);
  H = H_poisson;
end


%% 2.1 RSD propagation
if ~isConfocal
  tic; % Start timer

  spad_x = spadPos(:,:,1); % x-coordinates of SPAD positions
  spad_y = spadPos(:,:,3); % y-coordinates of SPAD positions
  [spadRows, spadCols] = size(spad_x); % Number of SPADs
  n = spadRows;

  laser_x = laserPos(:,:,1); % x-coordinates of laser positions
  laser_y = laserPos(:,:,3); % y-coordinates of laser positions
  [laserRows, laserCols] = size(laser_x); % Number of lasers

  delta_x = abs(mean(diff(spad_x(:,1))));
  x_dif = (-(n-1):n-1) * delta_x;
  y_dif = (-(n-1):n-1) * delta_x;
  [X_dif, Y_dif] = meshgrid(x_dif, y_dif);

  x_dif_laser = (-(n/2):(n-1)/2) * delta_x;
  y_dif_laser = (-(n/2):(n-1)/2) * delta_x;
  [X_dif_laser, Y_dif_laser] = meshgrid(x_dif_laser, y_dif_laser);

  H_omega_prime = fft(H,[],5);

  % Temporal frequency component of H
  nT = size(H, 5);
  t = t0 + (-(nT/2):(nT-1)/2) * deltaT; % Time vector

  omega0 = fftfreq(nT, deltaT);
  lambda_c = 1 ./ omega0; % Lambda vector
  n_omega = length(omega0); % Number of frequencies

  
  f_recon_spatial = zeros(spadRows, spadCols); % Initialize the reconstructed image

  lambdac = (4*delta_x); % λc = 2 * delta_x
  sigma = lambdac; %lambdac / (2 * log(2)); % 2* lambdac;
  % P_xt = exp(pre_exp_factor * (1/lambdac) * t) .* exp(-(t.^2) / (2 * sigma^2)); % P(xl,t)
  P_xt= exp(2*pi*1i*(1/lambdac)*t).*exp(-(t.^2)./(2*sigma.^2)); % size: [n_omega x nT]

  P_xomega = fftshift(P_xt); % Shift the zero frequency component to the center
  P_xomega = fft(P_xomega); % Fourier transform along the temporal dimension

  % plot(abs(P_xomega)); 
  % fprintf('P_xomega: %d\n', P_xomega); % Print the first element of P_xomega

  if multipleFocalPlanes
    volume_recon = zeros(num_planes, spadRows, spadCols);
    focal_planes = linspace(minPlaneDistance, maxPlaneDistance, num_planes);

    for idx = 1:num_planes
      focal_plane = focal_planes(idx);

      r = sqrt((X_dif).^2 + (Y_dif).^2 + (focal_plane).^2);
      r_laser = sqrt((X_dif_laser).^2 + (Y_dif_laser).^2 + (focal_plane).^2);
      % Loop over each frequency to compute the convolution
      for idx2 = 1:n_omega % 1:n_omega
        if ~isMultiFrequency && (abs(P_xomega(idx2)) < 13.1064)
          continue;
        end
        if isNarrowBand && abs(P_xomega(idx2)) < 1e-3
          continue;
        end
        omega = omega0(idx2); % Current frequency
        G_kernel = exp(2*pi*1i*omega.*r) ./ r;
        G_kernel_laser = exp(2*pi*1i*omega.*r_laser) ./ r_laser;
        
        % imagesc(abs(G_kernel_laser));
        % colormap hot; colorbar;

        fprintf('Processing weight %d (%d) of %d...\n', idx2, P_xomega(idx2), n_omega);

        h_flat = reshape(H_omega_prime(:,:,:,:,idx2), spadRows, spadCols);
        
        % imagesc(angle(h_flat.* P_xomega(idx)));
        % colormap hot; colorbar;
        % conv = conv2(h_flat, G_kernel, 'same');
        conv = conv_fft2(h_flat.* P_xomega(idx2), G_kernel, 'same');
        conv = conv .* G_kernel_laser; 
        f_recon_spatial = f_recon_spatial + conv';
      end
      % Save the kernel for each focal plane
      volume_recon(idx,:,:) = abs(f_recon_spatial);
      if focal_plane == minPlaneDistance
        figure;
        imagesc(abs(f_recon_spatial));
        colormap hot; colorbar;
        title('RSD: Spatial Convolution (Minimum Plane Distance)');
      end
    end
    recon_time = toc;
    fprintf('RSD reconstruction time: %.2f seconds\n', recon_time);
    %% Visualization of the reconstructed volume
    vs_h = volshow(volume_recon, RenderingStyle="MaximumIntensityProjection", Colormap=hot);
    vs_h.Parent.BackgroundColor = [0 0 0];
    vs_h.Parent.GradientColor = [0 0 0];
    % vs_h.Transformation = makehgtform('scale', [1, (maxPlaneDistance - minPlaneDistance)*10 / num_planes, 1]);
  else
    r = sqrt((X_dif).^2 + (Y_dif).^2 + (focal_plane).^2);
    r_laser = sqrt((X_dif_laser).^2 + (Y_dif_laser).^2 + (focal_plane).^2);
    
    % Loop over each frequency to compute the convolution
    for idx = 1:n_omega % 1:n_omega
      if ~isMultiFrequency && (abs(P_xomega(idx)) < 13.1064)
        continue;
      end
      if isNarrowBand && abs(P_xomega(idx)) < 1e-3
        continue;
      end
      omega = omega0(idx); % Current frequency
      G_kernel = exp(2*pi*1i*omega.*r) ./ r;
      G_kernel_laser = exp(2*pi*1i*omega.*r_laser) ./ r_laser;
      
      imagesc(abs(G_kernel_laser));
      colormap hot; colorbar;

      fprintf('Processing weight %d (%d) of %d...\n', idx, P_xomega(idx), n_omega);

      h_flat = reshape(H_omega_prime(:,:,:,:,idx), spadRows, spadCols);

      % imagesc(angle(h_flat.* P_xomega(idx)));
      % colormap hot; colorbar;

      % conv = conv2(h_flat, G_kernel, 'same');
      conv = conv_fft2(h_flat.* P_xomega(idx), G_kernel, 'same');
      conv = conv .* G_kernel_laser; 
      f_recon_spatial = f_recon_spatial + conv';
    end
    imagesc(abs(f_recon_spatial));
    colormap hot; colorbar;
    title('RSD: Spatial Convolution');
    recon_time = toc;
    fprintf('RSD reconstruction time: %.2f seconds\n', recon_time);

    figure;
    imagesc(abs(h_flat.* P_xomega(idx)));
    colormap hot; colorbar;
    title('Magnitude of H\_omega');
    figure;
    imagesc(abs(G_kernel));
    colormap hot; colorbar;
    title('Magnitude of G\_kernel');
  end
end