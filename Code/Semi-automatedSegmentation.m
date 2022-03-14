%% Created by Dhruv Khatri, IISER Pune 
% Get all the tiff files in the directory 
all_files = dir('**/*.tif');
% iterate over the  all the files, they could be multi tiff 
%all_data = {}; 
%count_struct =  1 ; 
%{
for f = 1:length(all_files)
    cur_file = strcat(all_files(f).folder,   '/',  all_files(f).name); 
    num_frames = imfinfo(cur_file); 
    for s = 1:length(num_frames)
        im_data = imread(cur_file, s); 
        seg_contour = get_contour(im_data); 
        all_data(count_struct).Filename = all_files(f).name; 
        all_data(count_struct).Frame = s; 
        all_data(count_struct).Contour = seg_contour;  
        count_struct  = count_struct +  1 ; 
        pause(0.5)
    end 
end 
%} 
%% Run the segmentation 
load segmnetation_data.mat
final_struct = {}; 
count_fig = 1 ; 
for f = 1:length(all_data)
    %figure(1), imshow(all_data(f).Contour)
    file_name = all_data(f).Filename ;  
    if f < 34
        file_path_folder = all_files(strcmp({all_files.name}, file_name)).folder;
        
    elseif f > 34 && f < 40
        index = find(strcmp({all_files.name}, file_name)); 
        file_path_folder = all_files(index(2)).folder;
    else
        file_path_folder = all_files(strcmp({all_files.name}, file_name)).folder;
    end 
    frame_number = all_data(f).Frame; 
    
    full_path  = strcat(file_path_folder, '/', file_name); 
    img_data = imread(full_path, frame_number);
    
    contour_data = all_data(f).Contour;
    conn_cont = bwconncomp(contour_data, 8); 
    label_cont = bwlabel(contour_data); 
    overlay_color = labeloverlay(imadjust(img_data), label_cont, 'Colormap','autumn'); 
    figure(1), imshow(overlay_color);
    identifier = split(file_name, '.'); 
    saveas(gcf, strcat(identifier{1}, '_' ,num2str(count_fig), '.png'));
    count_fig = count_fig + 1 ; 
    %[length_values,width_values] = getLengthsAndWidths(contour_data, contour_data);
    %final_struct(f).filename = file_name; 
    %final_struct(f).length = length_values*0.045;
    %final_struct(f).width = width_values*0.045; 
    %dir_split = split(full_path, '/'); 
    %final_struct(f).category = dir_split{end-1}; 
end 
%% Plot category wise data 
load Length_Widths(micron).mat
get_cep = find(strcmp({final_struct.category}, 'ceph_new')); 
get_mg1 = find(strcmp({final_struct.category}, '20220309-mg1655'));
get_dh5 = find(strcmp({final_struct.category}, 'dh5a_new'));

dh5_data = [cat(1, final_struct(get_dh5).length), cat(1, final_struct(get_dh5).width)];

% remove outlier 
indx = find(dh5_data(:,2) > 0.7);
dh5_data(indx, :) = []; 
ceph_data = [cat(1, final_struct(get_cep).length), cat(1, final_struct(get_cep).width)]; 

indx = find(ceph_data(:,2) > 0.7);
ceph_data(indx, :) = []; 


mg1_data = [cat(1, final_struct(get_mg1).length), cat(1, final_struct(get_mg1).width)];

indx = find(mg1_data(:,2) > 0.7);
mg1_data(indx, :) = []; 


all_lengths = [dh5_data(:,1); ceph_data(:,1); mg1_data(:,1)];  
all_widths  = [dh5_data(:,2); ceph_data(:,2); mg1_data(:,2)]; 

%% Fit to saturation 
min_length  = min(all_lengths); 
min_width = min(all_widths); 

max_length  = max(all_lengths); 
max_width = max(all_widths); 



mg1_data(:,1) = (mg1_data(:,1) - min_length) / (max_length);
ceph_data(:,1) = (ceph_data(:,1) - min_length) / (max_length); 
dh5_data(:,1) = (dh5_data(:,1) - min_length) / (max_length); 


mg1_data(:,2) = (mg1_data(:,2) - min_width) / (max_width);
ceph_data(:,2) = (ceph_data(:,2) - min_width) / (max_width); 
dh5_data(:,2) = (dh5_data(:,2) - min_width) / (max_width); 



figure(1), hold on , scatter(mg1_data(:,1), mg1_data(:,2)),
scatter(ceph_data(:,1), ceph_data(:,2)),
scatter(dh5_data(:,1), dh5_data(:,2)); 
%%
all_lengths = (all_lengths -min_length)/max_length; 
all_widths = (all_widths - min_width)/max_width; 

xlabel('Length (um)')
ylabel('Width (um)')

set(gca, 'FontSize', 18)

[cff, ss] = saturation_funct(all_lengths, all_widths);
figure(1), plot(unique(sort(all_lengths)), cff(unique(sort(all_lengths))), 'r-', ...
    'LineWidth', 3.0); 

legend({'mg16', 'mg16 + ceph', 'dh5α', 'Saturation Fit'})


hold off
%%

%% Required Functions 
function [cf, s1] = saturation_funct(xx, yy)
x = xx; % length 
y = yy; % Width

model = '(x*b)/(x + c)';
s = fitoptions('Method', 'NonLinearLeastSquares');
ft1 = fittype(model, 'coefficients', {'b', 'c'}, 'options', s); 
[cf ,s1] = fit(x, y, ft1); 
end 
















function contour_output = get_contour(img_data)
%% Edge detection 
I_edge = edge(img_data, 'Canny', 0.20); 
gap_fill = bwmorph(I_edge, 'bridge'); 
fill_image = imfill(gap_fill, 'holes'); 
conn_holes = bwconncomp(fill_image, 8); 
prop_holes = regionprops(conn_holes, 'Area', 'Extent', ...
    'Solidity', 'PixelIdxList'); 
%% select solid objects

ref_image = zeros(size(img_data)); 
for obj = 1:length(prop_holes)
    get_solidity = prop_holes(obj).Solidity;
    get_area = prop_holes(obj).Area;
    if get_solidity >= 0.80  && get_area > 90
        get_pixel = prop_holes(obj).PixelIdxList;
        ref_image(get_pixel) = 1;
    end 
end 

ref_overlay = imoverlay(imadjust(img_data), ref_image, 'red'); 
fig = figure(1); imshow(ref_overlay); hold on

%contour_output = ref_overlay; 

q = []; 
while true
  c = round(ginput(1));
  sel = get(fig, 'SelectionType');
  if strcmpi(sel, 'alt'); break; end
  scatter(c(1),c(2), 110, 'filled');
  q=[ceil(q);ceil(c)]; 
end
pause(1.0), hold off

conn_comp = bwconncomp(ref_image, 8); 
region_props_final = regionprops(conn_comp, 'PixelList', 'PixelIdxList'); 
final_contours2 = zeros(size(img_data));
for  p = 1:length(region_props_final)
    pixellist = region_props_final(p).PixelList; 
    if ~any(ismember(pixellist, q, 'rows'))
        get_indices = region_props_final(p).PixelIdxList; 
        final_contours2(get_indices) = 1; 
    end  
end

act_cont  = activecontour(img_data, final_contours2, 30,...
    'edge', 'SmoothFactor', 1.0, 'ContractionBias', -0.1); 

contour_output =act_cont; 

end 