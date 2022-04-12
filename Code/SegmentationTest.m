img_data = imread('Test.tif'); 
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
%% Overlay the result
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
perim_image = bwperim(act_cont); 
figure(1), subplot(2,2,[1 ,3]), imshow(imoverlay(imadjust(img_data), perim_image, 'red')); 
title('Output overlay (red)', 'Fontsize', 12)
%% Plot the statistics 
[length_array,width_array] = getLengthsAndWidths(act_cont, img_data); 
figure(1), subplot(2,2,2), histogram(width_array); 
xlabel('Width (px)', 'Fontsize', 12)
figure(1), subplot(2,2,4), histogram(length_array); 
xlabel('Length (px)', 'Fontsize', 12)







