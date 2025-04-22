% To configure the parameters of the execution, edit the file config.m

% Load the reconstructed volume from the file
run('config.m');
loaded_data = load(volume_path);
G_loaded = loaded_data.G_filtered;

% Visualize the loaded volume
vs_h_loaded = volshow(G_loaded, RenderingStyle="MaximumIntensityProjection", Colormap=hot);
vs_h_loaded.Parent.BackgroundColor = [0 0 0];
vs_h_loaded.Parent.GradientColor = [0 0 0];