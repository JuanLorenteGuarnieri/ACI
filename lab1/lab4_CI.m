% clear; clc; close all;

% To configure the parameters of the execution, edit the file config.m

run('config.m');

function save_reconstructed_volume(volume_path, G_filtered, filename, nVoxels, isFiltered, isAttenuated, isPhasedFiltered, isPhasedOption1, recon_time)
% Save the reconstructed volume to a file
if ~exist('results', 'dir')
  mkdir('results');
end
save(volume_path, 'G_filtered', '-v7');

% Open the results.csv file and append a new line with the details
fileID = fopen('results.csv', 'a');
fprintf(fileID, '%s,%d,%d,%d,%d,%d,%.2f\n', filename, nVoxels, double(isFiltered), double(isAttenuated), double(isPhasedFiltered), double(isPhasedOption1), recon_time);
fclose(fileID);
end


%% 1. Load the dataset
dataset = load(dataset_path).data;

% Determine if the dataset is confocal
if ndims(dataset.data) == 3
  isConfocal = true;
  fprintf('The dataset is confocal.\n');
else
  isConfocal = false;
  fprintf('The dataset is non-confocal.\n');
end

%% 2. Define the voxel grid of the hidden scene
vol_center = dataset.volumePosition;   % 3x1 vector [x, y, z]
vol_size   = dataset.volumeSize;         % 3x1 vector [sx, sy, sz]
if isscalar(vol_size)
  vol_size = [vol_size, vol_size, vol_size];
end

num_voxels = [nVoxels,nVoxels,nVoxels];  % Voxel resolution (x, y, z)
voxel_res = vol_size ./ num_voxels;
x = linspace(vol_center(1) - vol_size(1)/2 + voxel_res(1)/2, vol_center(1) + vol_size(1)/2 - voxel_res(1)/2, num_voxels(1));
y = linspace(vol_center(2) - vol_size(2)/2 + voxel_res(2)/2, vol_center(2) + vol_size(2)/2 - voxel_res(2)/2, num_voxels(2));
z = linspace(vol_center(3) - vol_size(3)/2 + voxel_res(3)/2, vol_center(3) + vol_size(3)/2 - voxel_res(3)/2, num_voxels(3));

G = zeros(num_voxels);  % Initialize reconstructed volume

%% 3. Extract measurement parameters and positions
t0 = dataset.t0;          % Temporal offset (in meters, optical distance)
deltaT = dataset.deltaT;  % Temporal resolution (in meters)
H = dataset.data;         % Measurement data H

% Normal del relay wall (constante)
normal_wall = [0,1,0];

if isConfocal
  % For confocal data, laser and SPAD positions coincide.
  positions = dataset.spadPositions;  % NxMx3 array
  origin_laser = dataset.laserOrigin;   % 3x1 vector
  origin_spad  = dataset.spadOrigin;    % 3x1 vector

  % Precompute distances d1 y d4 para cada punto del muro
  % d1: distancia del dispositivo láser al punto del muro
  d1 = sqrt(sum((positions - reshape(origin_laser, [1,1,3])).^2, 3));
  % d4: distancia del punto del muro al dispositivo SPAD
  d4 = sqrt(sum((positions - reshape(origin_spad, [1,1,3])).^2, 3));

  [nRows, nCols, ~] = size(positions);
else
  % For non-confocal data, use separate grids.
  laserPos = dataset.laserPositions;  % MxNx3 array
  spadPos  = dataset.spadPositions;   % PxQx3 array
  origin_laser = dataset.laserOrigin;
  origin_spad  = dataset.spadOrigin;

  % Precompute constant distances for laser and SPAD points:
  d1_laser = sqrt(sum((laserPos - reshape(origin_laser, [1,1,3])).^2, 3)); % size: [M x N]
  d4_spad  = sqrt(sum((spadPos - reshape(origin_spad, [1,1,3])).^2, 3));   % size: [P x Q]

  [laserRows, laserCols, ~] = size(laserPos);
  [spadRows, spadCols, ~] = size(spadPos);

  % Para la corrección en no-confocal se asume que la normal en láser y SPAD es la misma.
  normal_laser = normal_wall;
  normal_spad = normal_wall;
end

if isPhasedFiltered
  %% Phasor-based filtering (Section 5)
  % The goal is to filter the temporal dimension of H using a 1D Morlet wavelet.
  % We define a complex-valued filter:
  %   Km(t) = exp(2j*pi*Omega_c*t) .* exp(-t.^2/(2*sigma^2))
  % where Omega_c = 1/lambda_c and sigma is chosen in [lambda_c/(2*log2), 2*lambda_c].
  %
  % Here, we define t (the temporal axis) in optical distance units.
  if isConfocal
    % For confocal data, time dimension is the 3rd dimension.
    nT = size(H, 3);
    t = t0 + (0:nT-1) * deltaT;  % time vector in meters (optical distance)

    % Parameters
    lambda_c = 0.1;         % chosen wavelength in meters (ensure it's at least twice the relay wall spacing)
    Omega_c = 1 / lambda_c;   % central frequency
    if isPhasedOption1
      sigma = lambda_c / (2*log(2));
    else
      sigma = 2 * lambda_c;
    end
    % Create the complex-valued Morlet wavelet filter
    Km = exp(2j*pi*Omega_c*t) .* exp(-t.^2/(2*sigma^2));

    % Apply convolution as multiplication in Fourier domain along the 3rd dimension
    Km_fft = fft(Km);  % FFT of the filter (1D)
    H_fft = fft(H, [], 3);  % FFT of H along time dimension
    H_filtered = ifft(H_fft .* reshape(Km_fft, [1 1 nT]), [], 3);
  else
    % For non-confocal data, time dimension is the 5th dimension.
    nT = size(H, 5);
    t = t0 + (0:nT-1) * deltaT;  % time vector in meters (optical distance)

    lambda_c = 0.1;         % chosen wavelength in meters
    Omega_c = 1 / lambda_c;   % central frequency
    if isPhasedOption1
      sigma = lambda_c / (2*log(2));
    else
      sigma = 2 * lambda_c;
    end

    Km = exp(2j*pi*Omega_c*t) .* exp(-t.^2/(2*sigma^2));

    Km_fft = fft(Km);
    H_fft = fft(H, [], 5);  % FFT along the 5th (temporal) dimension
    H_filtered = ifft(H_fft .* reshape(Km_fft, [1 1 1 1 nT]), [], 5);
  end

  % Replace H with its filtered counterpart for use in the backprojection.
  H = H_filtered;
end

%% 4. Back-projection reconstruction with attenuation and foreshortening correction
fprintf('Starting back-projection reconstruction...\n');
tic;
if isConfocal
  % --- Confocal: partial vectorization over the wall grid ---
  if isAttenuated
    for ix = 1:num_voxels(1)
      for iy = 1:num_voxels(2)
        for iz = 1:num_voxels(3)
          % Voxel center
          xv = [x(ix), y(iy), z(iz)];
          % Vectorized: calculate d2 for each wall point
          diff = positions - reshape(xv, [1,1,3]);  % [nRows x nCols x 3]
          d2 = sqrt(sum(diff.^2, 3));                % d2 = distance from wall point to voxel
          % In confocal, d3 = d2
          % Calculate total_distance = d1 + d4 + 2*d2
          total_distance = d1 + d4 + 2*d2;
          % Calculate the time bin index
          tv = round((total_distance - t0) / deltaT);
          valid = tv >= 1 & tv <= size(H,3);

          if any(valid(:))
            % Calculate correction factor:
            %  - Cosine of the angle: (xv - pos_wall)/d2 · normal_wall
            %   Avoid division by zero since d2 > 0 (if d2==0 it would be treated separately).
            % Vectorized: for each valid point, calculate the cosine.
            dir_laser_voxel = diff ./ repmat(d2, [1,1,3]);  % [nRows x nCols x 3]
            % Dot product in the third dimension
            cos_wall = zeros(nRows, nCols);
            for k = 1:3
              cos_wall = cos_wall + dir_laser_voxel(:,:,k) * normal_wall(k);
            end
            cos_wall = max(cos_wall, eps);  % avoid zero values
            % For each wall point, the foreshortening factor is cos_wall^2
            foreshortening_factor = cos_wall.^2;
            % Quadratic factor: 1 + (d2^2 * d3^2) = 1 + d2^4
            quadratic_factor = 1 + d2.^4;

            % Extract valid contributions from H, applying the factors
            [r, c] = find(valid);
            tIdx = tv(valid);
            % Linear indices for H (dimensions: [nRows, nCols, nT])
            linInd = sub2ind(size(H), r, c, tIdx);
            % Apply point-wise multiplication with the factors
            factors = quadratic_factor(valid) .* foreshortening_factor(valid);
            % Accumulate in the voxel
            G(ix,iy,iz) = sum(H(linInd) .* factors);
          end
        end
      end
    end
  else
    for ix = 1:num_voxels(1)
      for iy = 1:num_voxels(2)
        for iz = 1:num_voxels(3)
          % Voxel center coordinates
          xv = [x(ix), y(iy), z(iz)];
          % Compute distance d2 from each wall point to voxel (vectorized)
          diff = positions - reshape(xv, [1 1 3]);  % size: [nRows x nCols x 3]
          d2 = sqrt(sum(diff.^2, 3));  % distance from each wall point to voxel

          % Total optical path length = d1 + d4 + 2*d2
          total_distance = d1 + d4 + 2*d2;

          % Compute corresponding time bin indices (vectorized)
          tv = round((total_distance - t0) / deltaT);
          % Create a logical mask for valid indices
          valid = tv >= 1 & tv <= size(H,3);

          % Convert 2D indices of wall grid to linear indices to extract H values
          if any(valid(:))
            % Subscript indices for valid wall points
            [r, c] = find(valid);
            % For each valid wall point, get the corresponding time bin from tv
            timeIdx = tv(valid);
            % Convert subscripts to linear indices in H for the third dimension.
            % For each valid wall point, H(r,c,timeIdx) is the contribution.
            linInd = sub2ind(size(H), r, c, timeIdx);
            % Sum contributions and assign to voxel (accumulate)
            G(ix,iy,iz) = sum(H(linInd));
          end
        end
      end
    end
  end

else
  % --- Non-confocal: partial vectorization over laser and SPAD grids ---
  if isAttenuated
    for ix = 1:num_voxels(1)
      for iy = 1:num_voxels(2)
        for iz = 1:num_voxels(3)
          xv = [x(ix), y(iy), z(iz)]; % Voxel center
          voxel_sum = 0;
          for i = 1:laserRows
            for j = 1:laserCols
              % Laser point
              pos_laser = squeeze(laserPos(i,j,:));  % 3x1
              d2 = norm(xv(:) - pos_laser(:));

              % Laser direction factor (constant for this point)
              dir_laser = (xv(:) - pos_laser(:)) / d2;
              cos_laser = max(dot(dir_laser, normal_laser), eps);
              d1_val = d1_laser(i,j); % d1 already precomputed for this point

              % Now, vectorize over the SPAD grid:
              % For each SPAD point, calculate d3 and obtain its factors
              % spadPos has size [spadRows x spadCols x 3]
              diff_spad = spadPos - reshape(xv, [1, 1, 3]);  % [spadRows x spadCols x 3]
              d3 = sqrt(sum(diff_spad.^2, 3));                % [spadRows x spadCols]

              % Calculate cosine for SPAD:
              dir_spad = diff_spad ./ d3;  % Vectorized operation (d3 is [spadRows x spadCols])
              cos_spad = dir_spad(:,:,1)*normal_spad(1) + ...
                dir_spad(:,:,2)*normal_spad(2) + ...
                dir_spad(:,:,3)*normal_spad(3);
              cos_spad = max(cos_spad, eps);

              % d4 already precomputed for SPAD:
              % d4_spad is of size [spadRows x spadCols]

              % Calculate the total distance for each laser-SPAD combination:
              % total_distance = d1 + d2 + d3 + d4
              total_distance = d1_val + d2 + d3 + d4_spad;  % [spadRows x spadCols]

              % Calculate the time bin index for each SPAD:
              tv = round((total_distance - t0) / deltaT);  % [spadRows x spadCols]
              valid = (tv >= 1) & (tv <= size(H,5));

              if any(valid(:))
                % Calculate quadratic factor: 1 + (d2^2 * d3^2)
                quad_factor = 1 + (d2^2) * (d3.^2);  % [spadRows x spadCols]
                % Foreshortening factor: cos_laser * cos_spad (cos_laser is scalar)
                fs_factor = cos_laser * cos_spad;  % [spadRows x spadCols]
                factors = quad_factor .* fs_factor;  % [spadRows x spadCols]

                % Extract contributions from H for this pair (i,j) and each SPAD
                % H has dimension: [laserRows x laserCols x spadRows x spadCols x nT]
                % Use sub2ind for each valid element of the SPAD grid.
                [Ispad, Jspad] = ndgrid(1:spadRows, 1:spadCols);
                Ispad_valid = Ispad(valid);
                Jspad_valid = Jspad(valid);
                tv_valid = tv(valid);

                % Convert subscripts to linear indices in H
                linInd = sub2ind(size(H), i*ones(size(Ispad_valid)), j*ones(size(Ispad_valid)), Ispad_valid, Jspad_valid, tv_valid);

                % Sum contribution for this laser point
                voxel_sum = voxel_sum + sum(H(linInd) .* factors(valid));
              end
            end
          end
          G(ix,iy,iz) = G(ix,iy,iz) + voxel_sum;
        end
      end
    end
  else
    for ix = 1:num_voxels(1)
      for iy = 1:num_voxels(2)
        for iz = 1:num_voxels(3)
          xv = [x(ix), y(iy), z(iz)];
          % For laser grid: compute d2 for each laser point (MxN)
          diff_laser = laserPos - reshape(xv, [1 1 3]);
          d2 = sqrt(sum(diff_laser.^2, 3));
          % For SPAD grid: compute d3 for each spad point (PxQ)
          diff_spad = spadPos - reshape(xv, [1 1 3]);
          d3 = sqrt(sum(diff_spad.^2, 3));

          % Now form the 4D total distance by combining laser and spad grids:
          % Use bsxfun (or simply '+' in modern MATLAB) to add the matrices
          % total_distance = d1_laser (MxN) + d2 (MxN) + d3 (PxQ) + d4_spad (PxQ)
          % We want each combination of laser and spad to be considered.
          % Expand d1_laser and d2 to 4D: [laserRows x laserCols x spadRows x spadCols]
          total_distance = repmat(d1_laser + d2, [1, 1, spadRows, spadCols]) + ...
            repmat(permute(d3 + d4_spad, [3 4 1 2]), [laserRows, laserCols, 1, 1]);

          % Compute time bin indices
          tv = round((total_distance - t0) / deltaT);
          % Create a mask for valid indices
          valid = tv >= 1 & tv <= size(H,5);

          % Sum contributions from valid measurements.
          if any(valid(:))
            % For each valid combination, obtain the value of H.
            % H has dimension [laserRows x laserCols x spadRows x spadCols x nT]
            % We use linear indices: we need to extract H(i,j,k,l,t) where t is obtained from tv.
            [iLaser, jLaser, iSpad, jSpad] = ndgrid(1:laserRows, 1:laserCols, 1:spadRows, 1:spadCols);
            % Only take the valid indices:
            idx_valid = valid(:);
            tIdx = tv(idx_valid);
            % Convert subscripts to linear indices for H.
            linInd = sub2ind(size(H), iLaser(idx_valid), jLaser(idx_valid), iSpad(idx_valid), jSpad(idx_valid), tIdx);
            G(ix,iy,iz) = sum(H(linInd));
          end
        end
      end
    end
  end
end

recon_time = toc;
fprintf('Reconstruction completed in %.2f seconds.\n', recon_time);

%% Apply Laplacian filter to enhance contrasts
f_lap = fspecial3('lap');
if isFiltered
  G_filtered = imfilter(G, -f_lap, 'symmetric');
else
  G_filtered = G;
end
if isPhasedFiltered
  G_filtered = abs(G_filtered);
end
%% 6. Visualization of the reconstructed volume
vs_h = volshow(G_filtered, RenderingStyle="MaximumIntensityProjection", Colormap=hot);
vs_h.Parent.BackgroundColor = [0 0 0];
vs_h.Parent.GradientColor = [0 0 0];

% Alternatively, you can use volumeViewer:
% volumeViewer(G_filtered, 'Colormap', hot);
save_reconstructed_volume(volume_path, G_filtered, name, nVoxels, isFiltered, isAttenuated, isPhasedFiltered, isPhasedOption1, recon_time);