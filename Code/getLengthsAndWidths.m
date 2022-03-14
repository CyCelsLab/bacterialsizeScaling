function [length_array,width_array] = getLengthsAndWidths(final_contours2, I)
conn_comps = bwconncomp(final_contours2, 8); 
region_connected = regionprops(conn_comps, 'PixelList', 'PixelIdxList');
L = bwlabel(final_contours2); 
%figure(1), subplot(2,2,4), imshow(crop_image); hold on 
length_array = []; 
width_array = []; 
for im = 1:length(region_connected)
    %%
    im_coordinates = region_connected(im).PixelList; 
    %im_indexes = region_connected(im).PixelIdxList; 
    reducedImage = isolateObject(im_coordinates,I);
    disp(im)
    [~, lengths, widths] = get_mid_array(reducedImage);  
    %plot(mid_array(:,2)+ min(im_coordinates(:,1)) - 1, ...
        %mid_array(:,1) + min(im_coordinates(:,2)) - 1 , 'r-', 'MarkerSize', 1);
    %  Compute lengths 
    length_array = [length_array; lengths]; 
    width_array = [width_array; widths*2]; 
    % compute width

end 
%% Required functions 
function [mid_points, length_total, width_value] = get_mid_array(reducedImage)
%%
get_angle = regionprops(reducedImage, 'Orientation'); 
rotatetheImage = imrotate(reducedImage, -(get_angle.Orientation-90), 'nearest'); 
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
% figure(1), histogram(min_distance,6), hold on 
% figure(1), line([median_value, median_value], ylim, 'LineWidth', 2.0, 'Color', 'k');
% figure(1), line([width_value,width_value], ylim, 'LineWidth', 2.0, 'Color', 'b')
% hold on 
[min_distance, indices] =  min(D, [], 2); 
median_value = median(min_distance) ; 
width_value = mean(min_distance(min_distance >= median_value)); 
%
%figure(4), imshow(rotatetheImage, 'InitialMagnification', 5000); hold on 
%plot(mid_array(:,2), mid_array(:,1)  , 'r-', ...
%   'MarkerSize', 1, 'LineWidth', 2.5);

%for w = 1:length(mid_array)
%    first_point = mid_array(w,:); 
%    second_point = [p_row(indices(w)), p_col(indices(w))];     
%    plot([mid_array(w,2) p_col(indices(w))], ...
%        [mid_array(w,1) p_row(indices(w))], 'b--', 'LineWidth', 2.5); 
    
    
%end 
%%
%hold off 
%figure(4), subplot(2,2,2), imshow(rotatetheImage)
end 

function reducedImage = isolateObject(pixel_coordinates, I) 
    test_im  = I; 
rescale_dims = [pixel_coordinates(:,1) - min(pixel_coordinates(:,1)) + 1, ...
    pixel_coordinates(:,2) - min(pixel_coordinates(:,2)) + 1]; 
null_image = zeros(max(rescale_dims(:,2)), max(rescale_dims(:,1))); 
burn_indx = sub2ind(size(null_image), rescale_dims(:,2), rescale_dims(:,1)); 
null_image(burn_indx) = 1;
reducedImage = null_image; 

width = max(pixel_coordinates(:,1)) - min(pixel_coordinates(:,1));
height = max(pixel_coordinates(:,2)) - min(pixel_coordinates(:,2));
%reducedImage = contour_original; 
% figure(2) ,imshow(null_image); 
end 

end

