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

%% Define the voxel grid of the hidden scene
vol_center = dataset.volumePosition;   % 3x1 vector [x, y, z]
vol_size   = dataset.volumeSize;         % 3x1 vector [sx, sy, sz]
if isscalar(vol_size)
  vol_size = [vol_size, vol_size, vol_size];
end

num_voxels = [nVoxels,nVoxels,nVoxels];  % Voxel resolution (x, y, z)
voxel_res = vol_size ./ num_voxels;

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

  % photon_counts = [1e5, 1e6, 1e7]; % valores típicos, puedes ajustarlos
  photon_counts = 1e7;
  H_poisson = poissrnd(photon_counts * H_norm);
  H = H_poisson;
end


%% 2.1 RSD propagation
if ~isConfocal && ~isMultiFrequency
  fprintf('Executing single-frequency RSD propagation...\n');
  tic; % Start timer

  spad_x = spadPos(:,:,1); % x-coordinates of SPAD positions
  spad_y = spadPos(:,:,3); % y-coordinates of SPAD positions
  [spadRows, spadCols] = size(spad_x); % Number of SPADs
  n = spadRows;

  delta_x = mean(diff(spad_x(:,1)));
  x_dif = (-(n-1):n-1) * delta_x;
  y_dif = (-(n-1):n-1) * delta_x;
  [X_dif, Y_dif] = meshgrid(x_dif, y_dif);

  lambda_c = 4 * delta_x; % λc = 2 * delta_x
  % lambda_c = 0.02;
  omega0 = 1 / lambda_c; % Ω₀ = 1 / λc

  % Temporal frequency component of H
  nT = size(H, 5);
  t = t0 + (0:nT-1) * deltaT; % Time vector
  exp_factor = exp(2*pi*1i*omega0.*t); % e^(2πiΩ₀t)
  H_omega = sum(H .* reshape(exp_factor, [1,1,1,1,nT]), 5) * deltaT; % HΩ'
  H_omega = reshape(H_omega, spadRows, spadCols)';

  % Define convolution kernel G
  if multipleFocalPlanes
    volume_recon = zeros(num_planes, spadRows, spadCols);
    focal_planes = linspace(minPlaneDistance, maxPlaneDistance, num_planes);

    for idx = 1:num_planes
      focal_plane = focal_planes(idx);

      r = sqrt((X_dif).^2 + (Y_dif).^2 + (focal_plane).^2);
      G_kernel = exp(2*pi*1i*omega0.*r) ./ r;
      f_recon_spatial = conv2(H_omega, G_kernel, 'same');

      % Save the kernel for each focal plane
      volume_recon(idx,:,:) = abs(f_recon_spatial);
    end
    %% Visualization of the reconstructed volume
    vs_h = volshow(volume_recon, RenderingStyle="MaximumIntensityProjection", Colormap=hot);
    vs_h.Parent.BackgroundColor = [0 0 0];
    vs_h.Parent.GradientColor = [0 0 0];
  else
    focal_plane = 0.5; % Focal plane distance
    r = sqrt((X_dif).^2 + (Y_dif).^2 + (focal_plane).^2);
    G_kernel = exp(2*pi*1i*omega0.*r) ./ r;
    f_recon_spatial = conv2(H_omega, G_kernel, 'same');

    figure;
    imagesc(abs(H_omega));
    colormap hot; colorbar;
    title('Magnitude of H\_omega');
    figure;
    imagesc(abs(G_kernel));
    colormap hot; colorbar;
    title('Magnitude of G\_kernel');
    figure;
    imagesc(abs(f_recon_spatial));
    colormap hot; colorbar;
    title('RSD: Spatial Convolution');
    recon_time = toc;
    fprintf('RSD reconstruction time: %.2f seconds\n', recon_time);
  end
elseif ~isConfocal && isMultiFrequency
  
  spad_x = spadPos(:,:,1); % x-coordinates of SPAD positions
  spad_y = spadPos(:,:,3); % y-coordinates of SPAD positions
  [spadRows, spadCols] = size(spad_x); % Number of SPADs
  n = spadRows;

  delta_x = mean(diff(spad_x(:,1)));
  x_dif = (-(n-1):n-1) * delta_x;
  y_dif = (-(n-1):n-1) * delta_x;
  [X_dif, Y_dif] = meshgrid(x_dif, y_dif);

  
  % Temporal frequency component of H
  nT = size(H, 5);
  t = t0 + (0:nT-1) * deltaT; % Time vector

  % lambda_c = 4 * delta_x; % λc = 2 * delta_x
  % lambda_c = 0.02;
  % omega0 = 1 / lambda_c; % Ω₀ = 1 / λc
  
  
  % Temporal frequency vector (double precision)
  % lambda_c = double(fftfreq(nT, (4*deltaT))); 
  % omega0 = 1 ./ lambda_c; % Ω₀ = 1 / λc
  % omega0 = linspace(-1/(2*deltaT), 1/(2*deltaT), nT); % Omega vector
  omega0 = fftfreq(nT, deltaT);
  lambda_c = 1 ./ omega0; % Lambda vector
  fprintf("omega0 max: %f\n", min(omega0));
  n_omega = length(omega0); % Number of frequencies

  focal_plane = 0.5; % Focal plane distance
    
  % Compute the convolution kernel G for the current frequency
  r = sqrt((X_dif).^2 + (Y_dif).^2 + (focal_plane).^2);

  pre_exp_factor = 2*pi*1i;

  f_recon_spatial = zeros(spadRows, spadCols); % Initialize the reconstructed image
  
  % Precalcular todos los factores exponenciales
  exp_factors = exp(pre_exp_factor * omega0(:) .* t(:)');  % size: [n_omega x nT]

  % Loop over each frequency to compute the convolution
  for idx = 1:n_omega % 1:n_omega
    % if abs(omega0(idx)) > abs(1/(3.99*delta_x)) || abs(omega0(idx)) < abs(1/(4.01*delta_x)) || omega0(idx) == 0
    %   continue;
    % end
    % if abs(omega0(idx)) > abs(1/(2*deltaT)) || abs(omega0(idx)) < abs(1/(2.2*deltaT)) || omega0(idx) == 0
    %   continue;
    % end
    omega = omega0(idx); % Current frequency
    fprintf('Processing frequency %d (%d) of %d...\n', idx, lambda_c(idx), n_omega);

    % exp_factor = exp(pre_exp_factor*omega.*t); % e^(2πiΩ₀t)
    exp_factor = reshape(exp_factors(idx, :), [1,1,1,1,nT]);
    H_omega = sum(H .* reshape(exp_factor, [1,1,1,1,nT]), 5) * deltaT; % HΩ'
    H_omega = reshape(H_omega, spadRows, spadCols)';
    G_kernel = exp(pre_exp_factor*omega.*r) ./ r;
      
    % Perform the convolution in the spatial domain
    f_recon_spatial = f_recon_spatial + conv2(H_omega, G_kernel, 'same');
  end
  
  imagesc(abs(f_recon_spatial));
  colormap hot; colorbar;
  title('RSD: Spatial Convolution');
  recon_time = toc;
  fprintf('RSD reconstruction time: %.2f seconds\n', recon_time);
elseif ~isConfocal && isMultiFrequency

  %% Preprocesamiento espacial y parámetros
  spad_x = spadPos(:,:,1);     % Coordenadas x de los SPAD
  spad_y = spadPos(:,:,3);     % Coordenadas y de los SPAD
  [spadRows, spadCols] = size(spad_x); % Número de SPADs (se asume una matriz cuadrada)

  delta_x = mean(diff(spad_x(:,1)));
  x_dif = (-(spadRows-1):spadRows-1) * delta_x;
  y_dif = (-(spadRows-1):spadRows-1) * delta_x;
  [X_dif, Y_dif] = meshgrid(x_dif, y_dif);
  
  % Distancia focal fija y precomputación de la matriz de distancias r
  focal_plane = 0.5;
  r = sqrt(X_dif.^2 + Y_dif.^2 + focal_plane.^2);  % Tamaño: [2*n-1 x 2*n-1]

  fprintf('Empezamos precomputacion\n');

  %% Preprocesamiento temporal: cálculo de la transformada en frecuencia
  % Se asume que la dimensión 5 de H es el tiempo.
  nT = size(H, 5);
  t = t0 + (0:nT-1) * deltaT; % Time vector

  % Calculamos el vector de longitudes de onda y las frecuencias ω
  lambda_c = fftfreq(nT, deltaT);  % Se asume que fftfreq está definido o es una función externa
  omega0 = 1 ./ lambda_c;          % ω0 = 1 / λc
  n_omega = numel(omega0);         % Número de frecuencias

  % Generamos la matriz de factores exponenciales para todas las frecuencias
  % exp_factors: Tamaño [n_omega x nT]
  exp_factors = exp(2*pi*1i * omega0(:) * t);  

  % Reorganizamos H para que la dimensión temporal sea la última
  % Se asume que H es de tamaño [spadRows x spadCols x 1 x 1 x nT] o similar.
  % La reorganizamos a una matriz 2D de tamaño [spadRows*spadCols, nT]
  H_2d = reshape(H, spadRows*spadCols, nT);

  % Calculamos los H_omega para TODAS las frecuencias en bloque:
  % Cada frecuencia corresponde a una combinación lineal de las muestras en tiempo.
  % resultado: [spadRows*spadCols x n_omega]
  H_omega_all = deltaT * (H_2d * exp_factors.');  
  % Lo reorganizamos a [spadRows x spadCols x n_omega]
  H_omega_all = reshape(H_omega_all, spadRows, spadCols, n_omega);

  %% Calcular el kernel G para todas las frecuencias de forma vectorizada
  
  fprintf('Calculamos kernel\n');
  % r tiene tamaño [2*n-1 x 2*n-1] y queremos calcular G_kernel para cada ω:
  % Usamos expansión implícita: se crea un arreglo 3D de tamaño [2*n-1 x 2*n-1 x n_omega]
  G_kernel_all = exp(2*pi*1i * (r .* reshape(omega0, [1, 1, n_omega]))) ./ r;  
  % Ahora, G_kernel_all(:,:,f) corresponde al kernel para la frecuencia f.

  %% Convolución vectorizada: uso de FFT 2D en “batch” sobre la tercera dimensión
  % Para conv2 con 'same' se debe definir un tamaño de FFT que cubra la suma convolutiva.
  % Tamaño de G_kernel: [M x N] con M = 2*n-1, N = 2*n-1.
  
  fprintf('Preparamos convolucion\n');
  [M, N] = size(r);
  % El tamaño del resultado de la convolución (lineal) es:
  P = spadRows + M - 1;  % filas
  Q = spadCols + N - 1;  % columnas

  % Calcula la FFT2 de H_omega_all en cada página (frecuencia)
  % La FFT operará en las dimensiones 1 y 2:
  H_fft = fft2(H_omega_all, P, Q);  % Tamaño: [P x Q x n_omega]

  % Calcula la FFT2 de G_kernel_all en cada frecuencia
  G_fft = fft2(G_kernel_all, P, Q);   % Tamaño: [P x Q x n_omega]

  % Multiplicación elemento a elemento (para cada frecuencia)
  conv_fft = H_fft .* G_fft;  % [P x Q x n_omega]
  % Obtener la convolución real (por iFFT2 para cada "página")
  
  fprintf('Convolucionamos\n');
  conv_all = ifft2(conv_fft, P, Q);  % [P x Q x n_omega]
  
  % Extraer la región central que corresponde a 'same'
  % Calculamos los índices de inicio y fin para recortar (suponiendo que P y Q son impares o pares)
  start_row = floor((M)/2) + 1;
  start_col = floor((N)/2) + 1;
  end_row = start_row + spadRows - 1;
  end_col = start_col + spadCols - 1;
  
  
  fprintf('Sumar contribuciones\n');
  % Extraer la región central para cada frecuencia y sumar todas las contribuciones
  f_recon_spatial = sum(conv_all(start_row:end_row, start_col:end_col, :), 3);

  %% Visualización y tiempos
  imagesc(abs(f_recon_spatial));
  colormap hot; colorbar;
  title('RSD: Spatial Convolution (Vectorizado)');
  recon_time = toc;
  fprintf('RSD reconstruction time (vectorizado): %.2f seconds\n', recon_time);
end
