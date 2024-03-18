%% Created by Dhruv Khatri, Modified by Tanvi Kale, IISER Pune 
% September 2022
clear all; close all;
se = strel('square', 2);
filter_size =5; 
sigma = 2; 
%% load in the test data 
imgPath = '/Users/tanvikale/Desktop/CA/Bac_cell_size/ImageAnalysis_code_demo/ROI1_normalised.tif';
numFrame = imfinfo(imgPath); 
disp(['Number of frames found ', num2str(length(numFrame))]); 
%% Show segmentation output of each frame
all_data = table(); 

inputImage = imread(imgPath, 1); 
figure(1), imshow(inputImage); 

%% OUTER EDGES
lap_filter_1 = fspecial('log',filter_size,sigma) ;
lap_im_out = imfilter(inputImage, lap_filter_1);
see = strel('square', 2); 
lap_im_out = bwskel(lap_im_out  > 0);
lap_im_out = imdilate(lap_im_out, see); 
figure(1), imshow(lap_im_out); 
    
edge_info = inputImage;
edge_info(lap_im_out) = 0; 
figure(2), imshow(edge_info); 

%% INNER EDGES 
lap_filter = fspecial('log',7,9);
lap_im = imfilter(inputImage, -lap_filter); 
binary_ = imbinarize(lap_im); 
binary_ = imfill(binary_, 'holes'); 
figure(1), imshow(binary_,[], 'Border', 'tight')

%% REfining boundaries with Active contours

op_contours  = activecontour(edge_info, binary_,100, 'Chan-Vese', ...
   'SmoothFactor', 1.0, 'ContractionBias', 0.0);
op_contours = imclearborder(op_contours); 
op_contours  = imfill(op_contours, 4,'holes'); 

op_contours = op_contours .* (~lap_im_out); 

%% %{Remove small objects} 
conn_comps = bwconncomp(op_contours, 8); 
region_connected = regionprops(conn_comps,'Area', 'PixelList', 'PixelIdxList', ...
    'Centroid');
minimum_area =  90; 
remove_objs = find((cat(1, region_connected.Area) < minimum_area)); 
region_connected(remove_objs) = [];
final_contours = logical(zeros(size(op_contours))); 
for obj  = 1:length(region_connected)
    indices  = region_connected(obj).PixelIdxList;
    final_contours(indices) = 1; 
end 

clean_overlay = show_overlay(inputImage, op_contours); 
clean_overlay = imoverlay(inputImage, bwperim(final_contours), 'red'); 
fig = figure(1); imshow(clean_overlay, 'Border','tight')
%saveas(gcf, ['ROI1', num2str(f), '.tiff'])
%figure(1); subplot(1,2,2), imshow(edge_info, 'InitialMagnification', 500)
%figure(1); subplot(1,3,3), imshow(lap_im,[], 'InitialMagnification', 500)

return_table = getLengthsAndWidths2(final_contours, inputImage); 
%return_table.Frame = repmat(f, height(return_table),1);
return_table.Centroid = cat(1, region_connected.Centroid);
all_data= [all_data; return_table]; 


%% Required functions 
function return_table = getLengthsAndWidths2(final_contours2, I)
return_table = array2table(zeros(0, 3)); 
return_table.Properties.VariableNames = {'ObjCoordinates', 'Length', 'Width'}; 

conn_comps = bwconncomp(final_contours2, 4); 
region_connected = regionprops(conn_comps, 'PixelList', 'PixelIdxList');

for im = 1:length(region_connected)
    %%
    im_coordinates = region_connected(im).PixelList; 
    reducedImage = isolateObject(im_coordinates,I);
    %%
    [~, lengths, widths] = get_mid_array(reducedImage);  
    %plot(mid_array(:,2)+ min(im_coordinates(:,1)) - 1, ...
        %mid_array(:,1) + min(im_coordinates(:,2)) - 1 , 'r-', 'MarkerSize', 1);
    %  Compute lengths 
    curr_data = {region_connected(im).PixelIdxList,lengths, widths}; 
    return_table = [return_table; curr_data]; 
    %length_array = [length_array; lengths]; 
    %width_array = [width_array; widths*2]; 
    if im == 82
        disp([lengths, widths]);
    end 

end 
%% Required functions 
function [mid_points, length_total, width_value] = get_mid_array(reducedImage)
%%
get_angle = regionprops(reducedImage, 'Orientation'); 
rotatetheImage = imrotate(reducedImage, -(get_angle.Orientation-90), 'nearest', 'loose'); 
rotatetheImage(rotatetheImage >= 1) = 1; 
size_image_f = size(rotatetheImage); 

[row, col] =find(rotatetheImage); 
[~, major_axis] = max([length(unique(row)), length(unique(col))]); 
length_dimension = major_axis; % If 1 then go through rows else columns 
mid_array = []; 
for n = 1:max(size_image_f)
  get_all_entries = find(rotatetheImage(n,:)); 
   if ~isempty(get_all_entries)
       mid_point =  get_all_entries(1) + (get_all_entries(end) - get_all_entries(1))/2;
       if isequal(major_axis, 1)
        mid_array = [mid_array;n, mid_point]; 
       else 
           mid_array = [mid_array; mid_point, n];
       end 
   else 
       continue
   end 
end 
%% 
output = fit(mid_array(:,1), mid_array(:,2), 'smoothingspline',...
    'SmoothingParam',0.2); 
mid_array(:,2) = smooth(mid_array(:,2), 0.3); 
mid_array(:,2)  = output(mid_array(:,1)); 
mid_points = mid_array; 
length_segment = hypot(diff(mid_array(:,1)), diff(mid_array(:,2))); 
length_total =  sum(length_segment);  
% get perimeter contours

[p_row, p_col] = find(bwperim(rotatetheImage)); 
D = pdist2(mid_array, [p_row, p_col], 'euclidean'); 

[min_distance, ~] =  min(D, [], 2); 
median_value = median(min_distance) ; 
width_value = mean(min_distance(min_distance >= median_value)); 
width_value = width_value*2;
 %figure(1), histogram(min_distance,6), hold on 
 %figure(1), line([median_value, median_value], ylim, 'LineWidth', 2.0, 'Color', 'k');
 %figure(1), line([width_value,width_value], ylim, 'LineWidth', 2.0, 'Color', 'b')
 %hold on 
 
%{
 figure(4), imshow(rotatetheImage, [],'InitialMagnification', 5000); hold on 
plot(mid_array(:,2), mid_array(:,1)  , 'r-', ...
    'MarkerSize', 1, 'LineWidth', 2.5);
 for w = 1:length(mid_array)
     first_point = mid_array(w,:); 
     second_point = [p_row(indices(w)), p_col(indices(w))];     
     plot([mid_array(w,2) p_col(indices(w))], ...
         [mid_array(w,1) p_row(indices(w))], 'b--', 'LineWidth', 2.5); 
 end 
%}
%%%figure(4), subplot(2,2,2), imshow(rotatetheImage)
end 

function reducedImage = isolateObject(pixel_coordinates, I) 
%test_im  = I; 
rescale_dims = [pixel_coordinates(:,1) - min(pixel_coordinates(:,1)) + 1, ...
    pixel_coordinates(:,2) - min(pixel_coordinates(:,2)) + 1]; 
null_image = zeros(max(rescale_dims(:,2)), max(rescale_dims(:,1))); 
burn_indx = sub2ind(size(null_image), rescale_dims(:,2), rescale_dims(:,1)); 
null_image(burn_indx) = 1;
reducedImage = null_image; 

%width = max(pixel_coordinates(:,1)) - min(pixel_coordinates(:,1));
%height = max(pixel_coordinates(:,2)) - min(pixel_coordinates(:,2));
%reducedImage = contour_original; 
% figure(2) ,imshow(null_image); 
end 

end

%%
function im =  show_overlay(im1, im2)
im2 = bwperim(im2);
im2 = bwlabel(im2); 
im = labeloverlay(im1, im2); 
end 
