%% Author: Dhruv Khatri

% Email: dhruv.khatri@students.iiserpune.ac.in

% Work Place: Cycels, IISER pune
% Created: 30th September 2022
% Last update: 13th March 20223

% PhaseSeg_v3e: Script runs the segmentation pipeline on image data, 
% Image_4608.tif for demonstration 

% Inputs: 
% Image file: ".tiff" format phase contrast image of E. coli, "Image_4608.tif"

% Outputs: Script generates only visual output as matlab "figures"

% Example: Either run the whole script for the final L & W output or run
% them section wise to get output from individual code blocks

% Required m-files in the directory
% Following folders should be added to path before running the script
% folder: interpclosed
% folder: interparc
% folder: bresenham
% All credits to the original authors

clear all 
close all 
%% Input Parameters 
filename = "Image_4608.tif"; % Input filename

filterSize = 20; % kernel size for "LoG" filter
filterSigma = 4.0; % std for the filter kernel

contIteration = 200;  % no. of iterations in the active contour
contContraction = -1.45; 
contSmooth = 0.0; 

areaThresh = 550; % Parameters to filter our small debris in background (px)
majorAxisLength = 50; 

ignoreTips = 6; % distance to ignore from top == ignore cap widths
widthThresh = 0.80; % threshold  below which the invagination is detected
%% Input Image
I = imread(filename);
adjustI = imadjust(I);
%% Inner edges
inFil = fspecial('log', filterSize, filterSigma);
insideImage = imfilter(I, inFil);
optContour = activecontour(imcomplement(adjustI), insideImage, ...
    contIteration,'edge',...
    'ContractionBias',contContraction, 'SmoothFactor',contSmooth); 
optContour = imfill(optContour, 'holes'); 
optContour = imclearborder(optContour); 
figure(1), imshow(imoverlay(adjustI, bwperim(optContour), 'red')); 
%% Filter out 
props = regionprops(optContour, 'all'); 
filteredContour = zeros(size(I)); 
for p = 1:length(props)
    cur_props = props(p).Area;
    cur_axis = props(p).MajorAxisLength; 
    if cur_props > areaThresh && cur_axis > majorAxisLength
        indice = props(p).PixelIdxList;
        filteredContour(indice) = 1;
    end
end
objConnect = bwconncomp(filteredContour); 
perimProps = regionprops(objConnect, 'PixelIdxList', 'Area'); 
figure(1), imshow(adjustI); hold on 
for obj = 1:length(perimProps)
    currObj = perimProps(obj).PixelIdxList;
    null_image  = zeros(size(I)); 
    null_image(currObj) = 1; 
    perim_image = bwperim(null_image);
    skel_image = bwskel(logical(null_image)); 
    ends = find(bwmorph(skel_image, 'endpoints')); 
    coords = find(skel_image); 
    for e = 1:length(ends)
        start = ends(e); 
        D = bwdistgeodesic(skel_image,start, 'quasi-euclidean');
        skel_image(D <ignoreTips) = 0; 
    end 
    [B, L, N] = bwboundaries(perim_image); 
    pout = B{1}; 
    % get two closest points to the skeleton 
    yy_b = pout(:,1);
    xx_b = pout(:,2); 
    nn = 3*numel(xx_b); 
    xyq = interpclosed(xx_b, yy_b, 0:1/nn:1);
    xyq(1,:) = smooth(xyq(1,:), 0.02); % Smooth factor hardcoded
    xyq(2,:) = smooth(xyq(2,:),0.02);
    xx_b = xyq(1,:)'; 
    yy_b = xyq(2,:)'; 
    % Smooth the skeleton 
    [yy_s, xx_s] = ind2sub(size(I), find(skel_image));
    figure(1), plot(xx_b, yy_b, 'r'); 
    scatter(xx_s, yy_s, 2,'g'); %hold off 
    % column is the skeleton, row is the boundary
    dist_matrix = pdist2([xx_b, yy_b], [xx_s, yy_s]); 
    avg_dist = mean(min(dist_matrix, [] , 1)); 
    mat_size = size(dist_matrix); 
    possible_points = []; 
    dist_points = [];
    indx_val = []; 
    % Radius criterion
    for col = 1:mat_size(2)
        selec_col = dist_matrix(:,col); 
        [val, indxm] = min(selec_col); 
        A = [xx_b(indxm), yy_b(indxm)]; 
        B = [xx_s(col), yy_s(col)]; 
        if val < widthThresh*avg_dist
            possible_points = [possible_points; B]; 
            dist_points = [dist_points; val];
            indx_val = [indx_val; indxm]; 
        end 
    end 
    
    if ~isempty(possible_points)
        numP = size(possible_points); 
        possibleIndices = sub2ind(size(I), possible_points(:,2), possible_points(:,1)); 
        flagPoint = ones(numP(1),1); 
        if numP > 1
            for pp = 1:numP(1)
                currPoint = possible_points(pp,:);
                B = [xx_b(indx_val(pp)), yy_b(indx_val(pp))]; 
                currIndx = sub2ind(size(I), currPoint(2), currPoint(1)); 
                allIndices = assignCloseness(currIndx,currIndx, skel_image, 5); 
                [LIA, LocB] = ismember(possibleIndices, allIndices);
                distances_array = dist_points(LIA); 
                if any(distances_array < dist_points(pp))
                    flagPoint(pp) = 0; 
                else 
                    plot([currPoint(1), B(1)], [currPoint(2), B(2)], 'm-', 'LineWidth',3);
                end 
                
            end 
        end 

        refinedCoords = possibleIndices(logical(flagPoint)); 
        refinedBoundary =  [xx_b(indx_val(logical(flagPoint))), yy_b(indx_val(logical(flagPoint)))]; 
        for r = 1:length(refinedCoords)
            B = refinedBoundary(r, :); 
            [yyf, xxf]  = ind2sub(size(I), refinedCoords(r)); 
            if B(1) <= xxf
                x_line = [B(1), xxf];
                y_line = [B(2), yyf]; 
                ccc = 'g-' ;
                segLine = polyfit(x_line, y_line,1); 
                x_array = min(pout(:,2)):max(pout(:,2)); 
                y_array = polyval(segLine, x_array);
                indx_sear = x_array >= B(1); 
                [xi, yi] = polyxpoly(x_array(indx_sear),y_array(indx_sear), xx_b, yy_b);

                if length(xi) > 1
                    dist_select = pdist2([xi,yi], [B(1), B(2)]); 
                    [~, min_indx]= min(dist_select); 
                    [xl, yl] = bresenham(xi(min_indx),yi(min_indx),B(1), B(2));
                    plot(xl,yl, '-', 'LineWidth', 3); 
                    indice_burn = sub2ind(size(I), yl,xl);
                    filteredContour(indice_burn)= 0; 
                else 
                    [xl, yl] = bresenham(xi,yi,B(1), B(2));
                   plot(xl,yl, '-', 'LineWidth', 3); 
                   indice_burn = sub2ind(size(I), yl,xl);
                    filteredContour(indice_burn)= 0; 
                end 
            else
                x_line = [xxf, B(1)];
                y_line = [yyf, B(2)]; 
                ccc = 'y-';
                segLine = polyfit(x_line, y_line,1); 
                x_array = min(pout(:,2)):max(pout(:,2)); 
                y_array = polyval(segLine, x_array); 
                indx_sear = x_array <= B(1); 
                [xi,yi]= polyxpoly(x_array(indx_sear), y_array(indx_sear), xx_b, yy_b); 
                if length(xi) > 1
                    dist_select = pdist2([xi,yi], [B(1), B(2)]); 
                    [~, min_indx]= min(dist_select); 
                    [xl, yl] = bresenham(xi(min_indx),yi(min_indx),B(1), B(2));
                    plot(xl,yl, '-', 'LineWidth', 3); 
                    indice_burn = sub2ind(size(I), yl,xl);
                    filteredContour(indice_burn)= 0; 
                else 
                    [xl, yl] = bresenham(xi,yi,B(1), B(2));
                    plot(xl,yl, '-', 'LineWidth', 3); 
                    indice_burn = sub2ind(size(I), yl,xl);
                    filteredContour(indice_burn)= 0; 
                end 
            end 
        end 
    end 
end 
%% Visualize the final filteredContour
% Filter Once More to remove for artifacts created after spetum detection
connComp = bwconncomp(filteredContour, 4); 
props = regionprops(connComp, 'Area', 'MajorAxisLength' ,'PixelIdxList'); 
filteredContour3 = zeros(size(I)); 
for p = 1:length(props)
    cur_props = props(p).Area;
    cur_axis = props(p).MajorAxisLength; 
    if cur_props > areaThresh && cur_axis > majorAxisLength
        indice = props(p).PixelIdxList;
        filteredContour3(indice) = 1;
    end
end
%% Lengths and Width 
figure(2), subplot(1,2,1),  imshow(imoverlay(adjustI,bwperim(filteredContour3)),[])
figure(2), subplot(1,2,2),  imshow(label2rgb(bwlabel(filteredContour3,4), 'spring', 'c','shuffle'))
return_table = getLengthsAndWidths2(filteredContour3, I);
figure(3), subplot(1,3,1),scatter(return_table.Length, return_table.Width); 
xlabel('Length (px)')
ylabel('Width (px)')
figure(3),subplot(1,3,2), histogram(return_table.Length);
xlabel('Length (px)')
figure(3), subplot(1,3,3) , histogram(return_table.Width)
xlabel('Width (px)')
%% Required functions
% Get the distance from a point (x3, y3) to
% a line defined by two points (x1, y1) and (x2, y2);
% Reference: http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
% https://in.mathworks.com/matlabcentral/answers/514135-how-to-find-the-distance-of-a-point-to-a-line-knowing-all-the-points
function indxIgnore =  assignCloseness(startpoint1,indxIgnore,perim_image, upto_thresh)
%{
for a startpoint/Could be multiple find the two closest neigthbors
 and then keep extending the search while ignoring already visited indices
%} 
closeimage = double(perim_image); 
% returns two closest neighbours 
coordsaround = getAround(startpoint1,size(perim_image)); 
coordsaround = coordsaround(closeimage(coordsaround) == 1); 
coordsaround = setdiff(coordsaround, indxIgnore); 
if ~isempty(coordsaround)
   if length(indxIgnore) < 2*upto_thresh 
       indxIgnore = [indxIgnore; coordsaround];
       indxIgnore = assignCloseness(coordsaround, indxIgnore, perim_image, upto_thresh);   
   end
end 
%closeimage(coordsaround(closeimage(coordsaround) == 1)) = closeimage(startpoint1(1), startpoint1(2))*2;
%if ~isempty(find(closeimage, 1))
%    closeimage = assignCloseness()

%end 
end 
function returncoords = getAround(coordSearch,  imageSize)
returncoords = []; 
for l = 1:length(coordSearch)
[yy, xx] = ind2sub(imageSize, coordSearch(l));
aroundCoords = [yy-1 xx-1;
    yy xx-1;
    yy+1 xx-1;
    yy-1 xx;
    yy+1 xx;
    yy-1 xx+1;
    yy xx+1
    yy+1 xx+1];
returncoords = [returncoords;sub2ind(imageSize, aroundCoords(:,1), aroundCoords(:,2))]; 
end
end
function return_table = getLengthsAndWidths2(final_contours2, I)
return_table = array2table(zeros(0, 3)); 
return_table.Properties.VariableNames = {'ObjCoordinates', 'Length', 'Width'}; 

conn_comps = bwconncomp(final_contours2, 4); 
region_connected = regionprops(conn_comps, 'PixelList', 'PixelIdxList');
%L = bwlabel(final_contours2); 
%figure(1), subplot(2,2,4), imshow(crop_image); hold on 
%length_array = []; 
%width_array = []; 
for im = 1:length(region_connected)
    %%
    im_coordinates = region_connected(im).PixelList; 
    %im_indexes = region_connected(im).PixelIdxList; 
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

function reducedImage = isolateObject(pixel_coordinates, ~) 
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

%% Unused functions 
function curvatureValue = getCurvature(point1, point2, point3, imageSize)

    [y1,x1] = ind2sub(imageSize,point1); 
    [y2,x2] = ind2sub(imageSize, point2);
    [y3,x3] = ind2sub(imageSize,point3);

    curvatureValue = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
    sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2)); 
end 
function distance = GetPointLineDistance(x3,y3,x1,y1,x2,y2)
try
	
	% Find the numerator for our point-to-line distance formula.
	numerator = abs((x2 - x1) * (y1 - y3) - (x1 - x3) * (y2 - y1));
	
	% Find the denominator for our point-to-line distance formula.
	denominator = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
	
	% Compute the distance.
	distance = numerator ./ denominator;
catch ME
	errorMessage = sprintf('Error in program %s.\nError Message:\n%s',...
		mfilename, ME.message);
	uiwait(errordlg(errorMessage));
end
return; % from GetPointLineDistance()
end
function around = removeIndices(xx_plot, yy_plot, sImage)
    x1 = xx_plot(1);
    x2 = xx_plot(2);
    y1 = yy_plot(1);
    y2 = yy_plot(2); 
    m = (y2-y1) / ( x2 - x1); 
    b = y1 - m*x1;
    if  x1 < x2 
        x1 = floor(x1); 
        x2 = ceil(x2); 
        x_array = x1:x2; 
        y_array = (m*x_array + b);
    elseif x1 > x2
        x1 = ceil(x1);
        x2 = floor(x2); 
        x_array = x2:x1; 
        y_array = (m*x_array + b); 
    elseif x2 == x1 
        if y1 < y2
            y1 = floor(y1);
            y2 = ceil(y2); 
            y_array = y1:y2; 
        else
            y1 = ceil(y1);
            y2 = floor(y2); 
            y_array = y2:y1;
        end
        x_array = repmat(x1, length(y_array),1)';
    end 
    x_array = round(x_array);
    y_array = round(y_array); 
    subOmit = sub2ind(sImage, y_array, x_array); 
    around = getAround(subOmit, sImage); 
end 
function im =  show_overlay(im1, im2)
im2 = bwperim(im2); 
im = imoverlay(im1, im2, 'red'); 
end 
function K= getCurve(yy,xx, distance)
%% Taken from 
% https://in.mathworks.com/matlabcentral/answers/57194-how-to-find-the-sharp-turn-in-a-2d-line-curve#answer_6918

K = {}; 


for i = 1:distance 
    yy_mod = yy(i:distance:end); 
    xx_mod = xx(i:distance:end); 
    numPoints = length(xx_mod); 
    k = zeros(length(xx_mod),1); 
    
    
    for t = 1:numPoints
      if t == 1
        index1 = numPoints;
        index2 = t;
        index3 = t + 1;
      elseif t >= numPoints
        index1 = t-1;
        index2 = t;
        index3 = 1;
      else
        index1 = t-1;
        index2 = t;
        index3 = t + 1;
      end
    
      x1 = xx_mod(index1);
      y1 = yy_mod(index1);
      x2 = xx_mod(index2);
      y2 = yy_mod(index2);
      x3 = xx_mod(index3);
      y3 = yy_mod(index3);
    
    
      k(t) = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
      sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));

    end 
    K(i).Points = k; 
    K(i).xxMod = xx_mod; 
    K(i).yyMod = yy_mod; 
end 

end 

%------------- END OF CODE --------------
