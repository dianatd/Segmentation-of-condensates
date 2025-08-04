function contact_angle_table = ContactAnglesDT(im_ch1, im_ch2, sigma, sname, P)
%CONTACTANGLESDT Segments two-phase condensates and calculates contact angles.
%
%   This function takes two fluorescence channel images of phase-separated
%   condensates, segments the distinct phases, identifies the interfaces
%   between them, and calculates the contact angles by fitting circles to
%   the boundaries. It is designed to handle condensates with one or more
%   interfaces.
%
% SYNTAX
%   contact_angle_table = ContactAnglesDT(im_ch1, im_ch2, sigma, sname, P)
%
% INPUTS
%   im_ch1:      (matrix) 2D image matrix for channel 1.
%   im_ch2:      (matrix) 2D image matrix for channel 2.
%   sigma:       (double) Sigma value for the Gaussian filter used in segmentation.
%   sname:       (string/char) A sample name or identifier for the object.
%   P:           (double) Perimeter ratio threshold. Sub-regions with a
%                perimeter less than (max_perimeter / P) are discarded.
%
% OUTPUTS
%   contact_angle_table: (table) A table containing the calculated results
%                        for each valid interface found. Each row
%                        corresponds to one interface.
%
% DEPENDENCIES
%   - CircleFitByPratt.m (Available on MATLAB File Exchange)

%% 1. Define Parameters and Initialize Results Table
FLATTENING_DISK_RADIUS = 150; % Radius for morphological opening for background flattening.
DILATION_DISK_SIZE = 10;      % Size of structuring element for finding interface neighbors.
contact_angle_table = table(); % Initialize an empty table to store results.

%% 2. Segment Both Channels to Find Phase Boundaries
[boundaries_ch1, labeled_mask_ch1] = segment_channel(im_ch1, sigma, FLATTENING_DISK_RADIUS, P);
[boundaries_ch2, labeled_mask_ch2] = segment_channel(im_ch2, sigma, FLATTENING_DISK_RADIUS, P);

% Combine the masks to get the outer boundary of the entire condensate.
combined_mask = (labeled_mask_ch1 > 0) | (labeled_mask_ch2 > 0);
[boundary_outer, ~] = bwboundaries(combined_mask, 4, 'noholes');
if isempty(boundary_outer)
    return; % Exit if no objects are found
end

%% 3. Identify All Interfaces Between the Two Phases
% Create binary contours from the boundary coordinates.
contour_ch1 = create_contour_mask(boundaries_ch1, size(im_ch1));
contour_ch2 = create_contour_mask(boundaries_ch2, size(im_ch2));
contour_outer = create_contour_mask(boundary_outer(1), size(im_ch1)); % Assume one large object

% The interface exists where the individual contours are NOT on the outer edge.
interface_mask = (contour_ch1 + contour_ch2) - contour_outer;
[labeled_interfaces, num_interfaces] = bwlabel(interface_mask > 0);
if num_interfaces == 0
    return; % Exit if no interfaces are found
end

%% 4. Loop Through Each Interface to Calculate Angles
se = strel('disk', DILATION_DISK_SIZE);
for i = 1:num_interfaces
    % --- Isolate the current interface and find its neighboring phases ---
    current_interface_mask = (labeled_interfaces == i);
    
    % Dilate the interface to find which labeled regions from each channel it touches.
    neighbors_mask = imdilate(current_interface_mask, se) & ~current_interface_mask;
    
    neighbor_label_ch1 = unique(labeled_mask_ch1(neighbors_mask));
    neighbor_label_ch2 = unique(labeled_mask_ch2(neighbors_mask));
    
    % Remove the zero label (background) and ensure we found one neighbor in each channel.
    neighbor_label_ch1 = neighbor_label_ch1(neighbor_label_ch1 > 0);
    neighbor_label_ch2 = neighbor_label_ch2(neighbor_label_ch2 > 0);
    if isempty(neighbor_label_ch1) || isempty(neighbor_label_ch2)
        continue; % Skip if the interface doesn't connect two phases
    end
    
    % --- Get the specific boundaries for the three relevant contours ---
    boundary_arc_ch1 = create_contour_mask(boundaries_ch1(neighbor_label_ch1), size(im_ch1)) & contour_outer;
    boundary_arc_ch2 = create_contour_mask(boundaries_ch2(neighbor_label_ch2), size(im_ch2)) & contour_outer;
    
    % --- Fit circles and calculate geometry ---
    fitted_circles = FitCircles(boundary_arc_ch1, boundary_arc_ch2, current_interface_mask);
    
    % This part requires dot indexing because FitCircles returns a struct.
    center_ch1 = fitted_circles.center_ch1;
    radius_ch1 = fitted_circles.radius_ch1;
    center_ch2 = fitted_circles.center_ch2;
    radius_ch2 = fitted_circles.radius_ch2;
    center_interface = fitted_circles.center_interface;

    % Find intersection points of the two phase circles.
    [x_intersect, y_intersect] = circcirc(center_ch1(1), center_ch1(2), radius_ch1, center_ch2(1), center_ch2(2), radius_ch2);
    if all(isnan(x_intersect))
        continue; % Skip if circles don't intersect
    end
    
    % --- Calculate contact angles from tangent slopes ---
    slopes1 = get_tangent_slopes(center_ch1, center_ch2, center_interface, x_intersect(1), y_intersect(1));
    slopes2 = get_tangent_slopes(center_ch1, center_ch2, center_interface, x_intersect(2), y_intersect(2));
    
    ca_ch1_1 = 180 - abs(atand((slopes1.interface - slopes1.ch1) / (1 + slopes1.interface * slopes1.ch1)));
    ca_ch1_2 = 180 - abs(atand((slopes2.interface - slopes2.ch1) / (1 + slopes2.interface * slopes2.ch1)));
    
    ca_ch2_1 = 180 - abs(atand((slopes1.interface - slopes1.ch2) / (1 + slopes1.interface * slopes1.ch2)));
    ca_ch2_2 = 180 - abs(atand((slopes2.interface - slopes2.ch2) / (1 + slopes2.interface * slopes2.ch2)));
    
    % --- Calculate interfacial tension ratios ---
    avg_angle_ch1_rad = deg2rad(mean([ca_ch1_1, ca_ch1_2]));
    avg_angle_ch2_rad = deg2rad(mean([ca_ch2_1, ca_ch2_2]));
    
    intT_A_AB = -sin(avg_angle_ch2_rad) / sin(avg_angle_ch1_rad + avg_angle_ch2_rad);
    intT_B_AB = -sin(avg_angle_ch1_rad) / sin(avg_angle_ch1_rad + avg_angle_ch2_rad);
    
    % --- Store results for this interface in a temporary table ---
    results_this_interface = table(sname, ca_ch1_1, ca_ch1_2, ca_ch2_1, ca_ch2_2, ...
        fitted_circles.radius_ch1, fitted_circles.radius_ch2, fitted_circles.radius_interface, ...
        center_ch1(1), center_ch1(2), center_ch2(1), center_ch2(2), ...
        center_interface(1), center_interface(2), intT_A_AB, intT_B_AB, ...
        'VariableNames', {'SampleName', 'cAngles_ch1_1', 'cAngles_ch1_2', 'cAngles_ch2_1', 'cAngles_ch2_2', ...
        'R_ch1', 'R_ch2', 'R_int', 'C_ch1_x', 'C_ch1_y', 'C_ch2_x', 'C_ch2_y', ...
        'C_Int_x', 'C_Int_y', 'intT_A_AB', 'intT_B_AB'});
    
    % Append this interface's results to the main table.
    contact_angle_table = [contact_angle_table; results_this_interface];
end
end % End of main function ContactAnglesDT

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function [boundaries, labeled_mask] = segment_channel(img, sigma, disk_radius, P)
% Segments a single channel image to find boundaries of significant regions.
    img_gauss = imgaussfilt(img, sigma);
    se = strel('disk', disk_radius);
    img_flat = imadjust(img_gauss - imopen(img_gauss, se));
    
    bin_mask = imclearborder(imbinarize(img_flat));
    
    [boundaries_all, labeled_mask] = bwboundaries(bin_mask, 4, 'holes');
    
    props = regionprops(labeled_mask, 'Perimeter');
    if isempty(props), boundaries = {}; return; end
    
    perimeters = [props.Perimeter];
    max_perimeter = max(perimeters);
    
    labels_to_keep = find(perimeters >= max_perimeter / P);
    
    boundaries = boundaries_all(labels_to_keep);
    labeled_mask(~ismember(labeled_mask, labels_to_keep)) = 0;
end
%--------------------------------------------------------------------------
function contour_mask = create_contour_mask(boundaries, mask_size)
% Converts boundary coordinates from bwboundaries into a binary mask.
    contour_mask = zeros(mask_size);
    for k = 1:length(boundaries)
        boundary_pts = boundaries{k};
        idx = sub2ind(mask_size, boundary_pts(:, 1), boundary_pts(:, 2));
        contour_mask(idx) = 1;
    end
end
%--------------------------------------------------------------------------
function fitted_circles = FitCircles(boundary_arc_ch1, boundary_arc_ch2, interface_mask)
% Fits circles to the three boundaries using the Pratt method.

    % Fit circle for channel 1 arc
    [r_ch1, c_ch1] = find(boundary_arc_ch1);
    fit_ch1 = CircleFitByPratt([c_ch1, r_ch1]);
    
    % Fit circle for channel 2 arc
    [r_ch2, c_ch2] = find(boundary_arc_ch2);
    fit_ch2 = CircleFitByPratt([c_ch2, r_ch2]);
    
    % Get intersection points of the two phase circles
    [x_intersect, y_intersect] = circcirc(fit_ch1(1), fit_ch1(2), fit_ch1(3), fit_ch2(1), fit_ch2(2), fit_ch2(3));
    
    % Get the points of the interface itself
    [r_int, c_int] = find(interface_mask);
    len_Int = length(r_int);
    
    % Combine interface points with heavily weighted intersection points
    interface_points_c = [c_int; repelem(x_intersect(:), len_Int * 10)];
    interface_points_r = [r_int; repelem(y_intersect(:), len_Int * 10)];
    
    fit_interface = CircleFitByPratt([interface_points_c, interface_points_r]);
    
    % Return results in a STRUCTURE, which the main function expects.
    fitted_circles.center_ch1 = fit_ch1(1:2)';
    fitted_circles.radius_ch1 = fit_ch1(3);
    fitted_circles.center_ch2 = fit_ch2(1:2)';
    fitted_circles.radius_ch2 = fit_ch2(3);
    fitted_circles.center_interface = fit_interface(1:2)';
    fitted_circles.radius_interface = fit_interface(3);
end
%--------------------------------------------------------------------------
function slopes = get_tangent_slopes(c1, c2, c_int, x, y)
% Calculates the slopes of the tangents to three circles at a common point (x,y).
    slopes.ch1 = -(x - c1(1)) / (y - c1(2));
    slopes.ch2 = -(x - c2(1)) / (y - c2(2));
    slopes.interface = -(x - c_int(1)) / (y - c_int(2));
end


% function contact_angle_table = ContactAnglesDT(im_ch1, im_ch2, sigma, sname, P)
% %CONTACTANGLESDT Segments two-phase condensates and calculates contact angles.
% %
% %   This function takes two fluorescence channel images of phase-separated
% %   condensates, segments the distinct phases, identifies the interfaces
% %   between them, and calculates the contact angles by fitting circles to
% %   the boundaries. It is designed to handle condensates with one or more
% %   interfaces.
% %
% % SYNTAX
% %   contact_angle_table = ContactAnglesDT(im_ch1, im_ch2, im_bf, sigma, sname, P)
% %
% % INPUTS
% %   im_ch1:      (matrix) 2D image matrix for channel 1 (e.g., 488nm).
% %   im_ch2:      (matrix) 2D image matrix for channel 2 (e.g., 647nm).
% %   sigma:       (double) Sigma value for the Gaussian filter used in
% %                segmentation.
% %   sname:       (string/char) A sample name or identifier for the object
% %                being analyzed.
% %   P:           (double) Perimeter ratio threshold. Sub-regions with a
% %                perimeter less than (max_perimeter / P) are discarded.
% %
% % OUTPUTS
% %   contact_angle_table: (table) A table containing the calculated results
% %                        for each valid interface found. Each row
% %                        corresponds to one interface.
% %
% % DEPENDENCIES
% %   - CircleFitByPratt.m (Available on MATLAB File Exchange)
% %
% % See also: bwboundaries, regionprops, imgaussfilt
% 
% %% 1. Define Parameters and Initialize Results Table
% % Can be changed depending on data
% FLATTENING_DISK_RADIUS = 150; % Radius for morphological opening for background flattening.
% DILATION_DISK_SIZE = 10;      % Size of structuring element for finding interface neighbors.
% 
% % Initialize an empty table to store results for all interfaces.
% contact_angle_table = table();
% 
% %% 2. Segment Both Channels to Find Phase Boundaries
% % segment_channel is a local helper function
% [boundaries_ch1, labeled_mask_ch1] = segment_channel(im_ch1, sigma, FLATTENING_DISK_RADIUS, P);
% [boundaries_ch2, labeled_mask_ch2] = segment_channel(im_ch2, sigma, FLATTENING_DISK_RADIUS, P);
% 
% % Combine the masks to get the outer boundary of the entire condensate.
% combined_mask = (labeled_mask_ch1 > 0) | (labeled_mask_ch2 > 0);
% [boundary_outer, ~] = bwboundaries(combined_mask, 4, 'noholes');
% if isempty(boundary_outer)
%     return; % Exit if no objects are found
% end
% 
% %% 3. Identify All Interfaces Between the Two Phases
% % Create binary contours from the boundary coordinates.
% contour_ch1 = create_contour_mask(boundaries_ch1, size(im_ch1));
% contour_ch2 = create_contour_mask(boundaries_ch2, size(im_ch2));
% contour_outer = create_contour_mask(boundary_outer(1), size(im_ch1)); % Assume one large object
% 
% % The interface exists where the individual contours are NOT on the outer edge.
% interface_mask = (contour_ch1 + contour_ch2) - contour_outer;
% [labeled_interfaces, num_interfaces] = bwlabel(interface_mask > 0);
% 
% if num_interfaces == 0
%     return; % Exit if no interfaces are found
% end
% 
% %% 4. Loop Through Each Interface to Calculate Angles
% se = strel('disk', DILATION_DISK_SIZE);
% 
% for i = 1:num_interfaces
%     % --- Isolate the current interface and find its neighboring phases ---
%     current_interface_mask = (labeled_interfaces == i);
% 
%     % Dilate the interface to find which labeled regions from each channel it touches.
%     neighbors_mask = imdilate(current_interface_mask, se) & ~current_interface_mask;
% 
%     neighbor_label_ch1 = unique(labeled_mask_ch1(neighbors_mask));
%     neighbor_label_ch2 = unique(labeled_mask_ch2(neighbors_mask));
% 
%     % Remove the zero label (background) and ensure we found one neighbor in each channel.
%     neighbor_label_ch1 = neighbor_label_ch1(neighbor_label_ch1 > 0);
%     neighbor_label_ch2 = neighbor_label_ch2(neighbor_label_ch2 > 0);
%     if isempty(neighbor_label_ch1) || isempty(neighbor_label_ch2)
%         continue; % Skip if the interface doesn't connect two phases
%     end
% 
%     % --- Get the specific boundaries for the three relevant contours ---
%     boundary_arc_ch1 = create_contour_mask(boundaries_ch1(neighbor_label_ch1), size(im_ch1)) & contour_outer;
%     boundary_arc_ch2 = create_contour_mask(boundaries_ch2(neighbor_label_ch2), size(im_ch2)) & contour_outer;
% 
%     % --- Fit circles and calculate geometry ---
%     fitted_circles = FitCircles(boundary_arc_ch1, boundary_arc_ch2, current_interface_mask);
% 
%     center_ch1 = fitted_circles.center_ch1;
%     radius_ch1 = fitted_circles.radius_ch1;
%     center_ch2 = fitted_circles.center_ch2;
%     radius_ch2 = fitted_circles.radius_ch2;
%     center_interface = fitted_circles.center_interface;
% 
%     % Find intersection points of the two phase circles.
%     [x_intersect, y_intersect] = circcirc(center_ch1(1), center_ch1(2), radius_ch1, center_ch2(1), center_ch2(2), radius_ch2);
%     if all(isnan(x_intersect))
%         continue; % Skip if circles don't intersect
%     end
% 
%     % --- Calculate contact angles from tangent slopes ---
%     slopes1 = get_tangent_slopes(center_ch1, center_ch2, center_interface, x_intersect(1), y_intersect(1));
%     slopes2 = get_tangent_slopes(center_ch1, center_ch2, center_interface, x_intersect(2), y_intersect(2));
% 
%     ca_ch1_1 = 180 - abs(atand((slopes1.interface - slopes1.ch1) / (1 + slopes1.interface * slopes1.ch1)));
%     ca_ch1_2 = 180 - abs(atand((slopes2.interface - slopes2.ch1) / (1 + slopes2.interface * slopes2.ch1)));
% 
%     ca_ch2_1 = 180 - abs(atand((slopes1.interface - slopes1.ch2) / (1 + slopes1.interface * slopes1.ch2)));
%     ca_ch2_2 = 180 - abs(atand((slopes2.interface - slopes2.ch2) / (1 + slopes2.interface * slopes2.ch2)));
% 
%     % --- Calculate interfacial tension ratios ---
%     avg_angle_ch1_rad = deg2rad(mean([ca_ch1_1, ca_ch1_2]));
%     avg_angle_ch2_rad = deg2rad(mean([ca_ch2_1, ca_ch2_2]));
% 
%     intT_A_AB = -sin(avg_angle_ch2_rad) / sin(avg_angle_ch1_rad + avg_angle_ch2_rad);
%     intT_B_AB = -sin(avg_angle_ch1_rad) / sin(avg_angle_ch1_rad + avg_angle_ch2_rad);
% 
%     % --- Store results for this interface in a temporary table ---
%     results_this_interface = table(sname, ca_ch1_1, ca_ch1_2, ca_ch2_1, ca_ch2_2, ...
%         fitted_circles.radius_ch1, fitted_circles.radius_ch2, fitted_circles.radius_interface, ...
%         center_ch1(1), center_ch1(2), center_ch2(1), center_ch2(2), ...
%         center_interface(1), center_interface(2), intT_A_AB, intT_B_AB, ...
%         'VariableNames', {'SampleName', 'cAngles_ch1_1', 'cAngles_ch1_2', 'cAngles_ch2_1', 'cAngles_ch2_2', ...
%         'R_ch1', 'R_ch2', 'R_int', 'C_ch1_x', 'C_ch1_y', 'C_ch2_x', 'C_ch2_y', ...
%         'C_Int_x', 'C_Int_y', 'intT_A_AB', 'intT_B_AB'});
% 
%     % Append this interface's results to the main table.
%     contact_angle_table = [contact_angle_table; results_this_interface];
% end
% 
% end % End of main function ContactAnglesDT
% 
% %% ========================================================================
% %  LOCAL HELPER FUNCTIONS
% %  ========================================================================
% 
% function [boundaries, labeled_mask] = segment_channel(img, sigma, disk_radius, P)
% % Segments a single channel image to find boundaries of significant regions.
%     img_gauss = imgaussfilt(img, sigma);
%     se = strel('disk', disk_radius);
%     img_flat = imadjust(img_gauss - imopen(img_gauss, se));
% 
%     bin_mask = imclearborder(imbinarize(img_flat));
% 
%     [boundaries_all, labeled_mask] = bwboundaries(bin_mask, 4, 'holes');
% 
%     % Filter out small, insignificant regions based on perimeter.
%     props = regionprops(labeled_mask, 'Perimeter');
%     if isempty(props), boundaries = {}; return; end
% 
%     perimeters = [props.Perimeter];
%     max_perimeter = max(perimeters);
% 
%     % Find labels of regions to keep.
%     labels_to_keep = find(perimeters >= max_perimeter / P);
% 
%     % Filter the boundaries and labeled mask.
%     boundaries = boundaries_all(labels_to_keep);
%     labeled_mask(~ismember(labeled_mask, labels_to_keep)) = 0;
% end
% 
% %--------------------------------------------------------------------------
% 
% function contour_mask = create_contour_mask(boundaries, mask_size)
% % Converts boundary coordinates from bwboundaries into a binary mask.
%     contour_mask = zeros(mask_size);
%     for k = 1:length(boundaries)
%         boundary_pts = boundaries{k};
%         idx = sub2ind(mask_size, boundary_pts(:, 1), boundary_pts(:, 2));
%         contour_mask(idx) = 1;
%     end
% end
% 
% %--------------------------------------------------------------------------
% 
% function fitcircles = FitCircles(cboundary_ch1, cboundary_ch2, cboundary_interface)
%     % Fits circles to the boundaries of the two phases and their interface.
%     [c_ch1, r_ch1] = find(cboundary_ch1);
%     fit_ch1 = CircleFitByPratt([r_ch1, c_ch1]);
%     Center_ch1 = fit_ch1(1:2);
%     Radius_ch1 = fit_ch1(3);
% 
%     [c_ch2, r_ch2] = find(cboundary_ch2);
%     fit_ch2 = CircleFitByPratt([r_ch2, c_ch2]);
%     Center_ch2 = fit_ch2(1:2);
%     Radius_ch2 = fit_ch2(3);
% 
%     [Xout,Yout] = circcirc(Center_ch1(1),Center_ch1(2), Radius_ch1, Center_ch2(1),Center_ch2(2), Radius_ch2);
% 
%     [y_interface, x_interface] = find(cboundary_interface);
%     len_Int = length(y_interface);
% 
%     % Heavily weight the intersection points to ensure the interface circle passes through them
%     R_interface = [x_interface(:)', repelem(Xout, len_Int*10)]';
%     C_interface = [y_interface(:)', repelem(Yout, len_Int*10)]';
% 
%     fit_Interface = CircleFitByPratt([R_interface, C_interface]);
%     Center_interface = fit_Interface(1:2);
%     Radius_interface = fit_Interface(3);
% 
%     fitcircles = [Center_ch1, Radius_ch1, Center_ch2, Radius_ch2, Center_interface, Radius_interface];
% end
% 
% %--------------------------------------------------------------------------
% 
% function slopes = get_tangent_slopes(c1, c2, c_int, x, y)
% % Calculates the slopes of the tangents to three circles at a common point (x,y).
%     slopes.ch1 = -(x - c1(1)) / (y - c1(2)); % Note: this is -1/m_radius
%     slopes.ch2 = -(x - c2(1)) / (y - c2(2));
%     slopes.interface = -(x - c_int(1)) / (y - c_int(2));
% end