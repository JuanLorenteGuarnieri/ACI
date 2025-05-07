% lab2.m - Polarization-based Imaging Assignment

% --- PARAMETERS ---
datasetPath = 'data/apple'; % 'remote' or 'apple' or 'groot', or your own dataset folder

% --- FILENAMES (update if needed) ---
% List all .tif files in the dataset folder and sort them by name
fileList = dir(fullfile(datasetPath, '*.tif'));
[~, idx] = sort({fileList.name});
fileList = fileList(idx);

% Check that there are at least 9 images
if numel(fileList) < 9
  error('Not enough .tif images found in %s', datasetPath);
end

% Assign files in order: L, ih0, ih90, ih45, ih135, iv0, iv90, iv45, iv135
files = {fileList(1).name, fileList(2).name, fileList(3).name, fileList(4).name, ...
     fileList(5).name, fileList(6).name, fileList(7).name, fileList(8).name, fileList(9).name};

% --- READ IMAGES ---
imgs = cell(1,9);
for i = 1:9
  img = im2double(imread(fullfile(datasetPath, files{i})));
  if size(img,3) > 1
    % img = rgb2gray(img); 
  end
  imgs{i} = img;
end

% Assign images to variables for clarity
L      = imgs{1};
iv0    = imgs{2}; iv90  = imgs{3}; iv45  = imgs{4}; iv135  = imgs{5};
ih0    = imgs{6}; ih90  = imgs{7}; ih45  = imgs{8}; ih135  = imgs{9};


% Compute a global scale factor (e.g., max value across all images)
allPixels = cell2mat(cellfun(@(x) x(:), imgs, 'UniformOutput', false));
scaleMax = max(allPixels(:));

% Normalize images to the same scale for better visualization
for i = 1:9
  imgs{i} = imgs{i}/scaleMax;
end

% --- DISPLAY IMAGES ---
figure('Name','Polarization Images');
titles = {'L', 'iv0', 'iv90', 'iv45', 'iv135', 'ih0', 'ih90', 'ih45', 'ih135'};
for i = 1:9
  subplot(3,3,i);
  imagesc(imgs{i}); axis image off; 
  title(titles{i}, 'Interpreter','none');
end
sgtitle('Polarization Images (i_{h|v}\theta)');

% --- CALCULATE STOKES COMPONENTS ---
% s0 = L
sh1 = ih0 - ih90;
sh2 = ih45 - ih135;
sv1 = iv0 - iv90;
sv2 = iv45 - iv135;

sdp1 = sh1 + sv1;
sdp2 = sh2 + sv2;

ssp1 = sv1 - sdp1;
ssp2 = sv2 - sdp2;

sdp0 = sdp1 + sdp2;
ssp0 = ssp1 + ssp2;
s00 = L - sdp0 - ssp0;

% --- DISPLAY STOKES COMPONENTS ---
figure('Name','Stokes Components');
subplot(2,4,1); imagesc(sh1); axis image off; colormap gray; title('sh_{1}');
subplot(2,4,2); imagesc(sh2); axis image off; colormap gray; title('sh_{2}');
subplot(2,4,3); imagesc(sv1); axis image off; colormap gray; title('sv_{1}');
subplot(2,4,4); imagesc(sv2); axis image off; colormap gray; title('sv_{2}');
subplot(2,4,5); imagesc(sdp1); axis image off; colormap gray; title('sdp{1}');
subplot(2,4,6); imagesc(sdp2); axis image off; colormap gray; title('sdp_{2}');
subplot(2,4,7); imagesc(ssp1); axis image off; colormap gray; title('ssp_{1}');
subplot(2,4,8); imagesc(ssp2); axis image off; colormap gray; title('ssp_{2}');
sgtitle('Stokes Vector Components');

% --- DISPLAY REFLECTANCE SEPARATION ---
figure('Name','Reflectance Separation');
subplot(2,3,2); imagesc(L); axis image off; colormap gray; title('s_{0} = L');
subplot(2,3,4); imagesc(sdp0); axis image off; title('Diffuse Polarized');
subplot(2,3,5); imagesc(ssp0); axis image off; title('Specular Polarized');
subplot(2,3,6); imagesc(s00); axis image off; title('Diffuse Unpolarized');
sgtitle('Reflectance Separation');