%% IMAGE SEGMENTATION AND ANALYSIS SCRIPT
% This script analyzes two-channel fluorescence microscopy images (.tif) to
% quantify properties of phase-separated condensates.
%
% --- SCRIPT OVERVIEW ---
% 1.  **Setup**: Clears workspace and prompts the user for input parameters.
% 2.  **File Discovery**: Locates image files in the specified folder.
% 3.  **Main Image Loop**: Iterates through each image file.
%       a. Loads the two channels (e.g., ch1 and ch2).
%       b. Identifies all condensates (objects) in the image.
%       c. Displays the identified condensates for a quick visual check.
% 4.  **Object Analysis Loop**: Iterates through each identified condensate.
%       a. Crops the individual condensate for analysis using the specified centered method.
%       b. Pre-processes the cropped image (filtering, background flattening).
%       c. Creates binary masks for each channel using the multithresh method.
%       d. Calculates metrics: correlation, area, and intensity partitioning.
%       e. **For de-mixed condensates**: Attempts to calculate contact angles
%          by fitting circles to the phase boundaries and the interface. This
%          step includes a manual validation prompt.
% 5.  **Data Export**: Saves the calculated data into two separate CSV files.
%
% --- DEPENDENCIES ---
% This script requires the following custom functions to be in the MATLAB path:
% - ContactAnglesDT.m: Function to calculate contact angles.
% - CircleFitByPratt.m: Function for fitting circles to data points (MATLAB File Exchange)
% https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method
%
% Last edited: 07/30/2025

%% 1. SETUP AND PARAMETERS
% ========================================================================
close all;  % Close all open figures
clear all;  % Clear all variables from the workspace

% --- User-Defined Parameters ---
% An input dialog will collect all necessary parameters from the user.
prompt = {'Enter a unique name for this analysis run:',...
          'Enter basename (common string in filenames):',...
          'Enter perimeter threshold (pixels):',...
          'Enter Gaussian filter sigma (σ):',...
          'Enter contact angle perimeter ratio (P):',...
          'Enter pixel size (nm/pixel):'};
dlgtitle = 'Analysis Parameters';
dims = [1, 70];
definput = {'Anneal_55-25C', '2024', '200', '4', '50', '189.3'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

% Exit if user cancels the dialog
if isempty(answer); return; end

% Assign user inputs to variables
analysis_name = answer{1};      % A descriptor for the output filenames (e.g., 'Sel_55C-25C_01Cper1h')
basename = answer{2};           % A string contained in all files of interest (e.g., '20230215')
t_size = str2double(answer{3}); % Perimeter threshold. Objects smaller than this will be ignored.
sigma = str2double(answer{4});  % Sigma for Gaussian filtering to find object masks.
P = str2double(answer{5});      % Perimeter ratio used by the ContactAnglesDT function.
pixel_size = str2double(answer{6}); % Pixel size for area calculations (nm/pixel).

% 2. FILE DISCOVERY
% ========================================================================
% Prompt user to select the folder containing the image files
DataFolder = uigetdir('', 'Select the folder containing your image files');
if DataFolder == 0; return; end % Exit if user cancels

% Find all .tif files in the folder that contain the basename
search_pattern = fullfile(DataFolder, [basename, '*.tif']);
TPDfiles = dir(search_pattern);
ImageNames = {TPDfiles.name};

% Exit with an error if no matching files are found
if isempty(ImageNames)
    error('No image files found matching the basename "%s" in the selected folder.', basename);
end

% 3. IMAGE PROCESSING LOOP
% ========================================================================
% This main loop iterates over each image file found.
for point = 1:length(ImageNames)
    current_filename = ImageNames{point};
    fprintf('\nProcessing image %d of %d: %s\n', point, length(ImageNames), current_filename);
    
    % --- Image Loading & Initial Processing ---
    % Construct a unique filename for saving output files from this image
    [~, base_file, ~] = fileparts(current_filename);
    file_name = sprintf('%s_%s_t%d_sigma%d_P%d', analysis_name, base_file, t_size, sigma, P);
    
    % Read the two fluorescence channels from the multi-page .tif file
    Im_ch1 = imread(fullfile(DataFolder, current_filename), 1);
    Im_ch2 = imread(fullfile(DataFolder, current_filename), 2);
    
    % Create a fused, color image for visualization and initial object detection
    Im_fused = imfuse(Im_ch1, Im_ch2); % ch1 -> Green, ch2 -> Magenta
    Im = rgb2gray(Im_fused);           % Convert to grayscale for robust object finding
    
    % --- Global Object Detection ---
    % Find all condensates in the full image to loop through them later
    ImG = imgaussfilt(Im, sigma);
    bw_all = imbinarize(ImG);
    bw = imclearborder(bw_all); % Get rid of condensates touching the edge
    [B, L, ~] = bwboundaries(bw, 4, 'holes'); % Find boundaries (B) and a label matrix (L)
    %ImM = label2rgb(L, @summer, [0 0 0]); %optional labeling of each ROI found based on boundaries such that each condensate has a numbered label

    % --- Visualization of Detected Objects ---
    % Display the fused image with outlines of all valid objects drawn
    figure('Name', ['Detected Condensates in ' current_filename]);
    imshow(Im_fused);
    hold on;
    title(['Detected Condensates in: ', strrep(current_filename, '_', '\_')]);
    for k = 1:length(B)
        if size(B{k}, 1) > t_size
            plot(B{k}(:,2), B{k}(:,1), 'w', 'LineWidth', 2);
        end
    end
    exportgraphics(gcf, ["rho_Selected_condensates_" + file_name +  ".png"], 'Resolution', 300);
    fprintf('Detected %d objects. Pausing for review. Press any key to continue...\n', length(B));
    pause; % Pause to allow the user to inspect the detected objects
    hold off;
    close(gcf);
    
    % --- Initialize Data Storage For This Image ---
    % These arrays will store results for each valid condensate.
    sample_names = []; r_raw = []; r_ImFlat = []; corr2_Masks = [];
    A_ch1 = []; A_ch2 = []; Mm_ch1 = []; Mm_ch2 = []; NMm_ch1 = [];
    NMm_ch2 = []; NMm_ch1_back = []; NMm_ch2_back = [];
    Table_ContactAngles = table(); % Master table for all accepted contact angles in the image
   
    % 4. OBJECT ANALYSIS LOOP
    % ========================================================================
    % This loop iterates over each object 'k' found in the current image.
    for k = 1:length(B)
        % Reset the temporary table for each new condensate
        contactAngles = table();
        
        % Analyze the object only if its perimeter is larger than the threshold
        if size(B{k}, 1) > t_size
            fprintf(' -> Analyzing object %d... ', k);
            
            % --- Isolate and Crop the Current Object ---
            selected_object = (L == k); % Create a mask for only the current object
            boundary = B{k};
            
            % Create a region from the boundary to get its properties
            region = roipoly(bw, boundary(:,2), boundary(:,1));
            stats = regionprops(region, 'BoundingBox', 'Circularity', 'Centroid');
            
            % Define a cropping rectangle that is 2x the bounding box size,
            % centered on the object's centroid for a wider view.
            x_box = stats.BoundingBox(1);
            y_box = stats.BoundingBox(2);
            x_c = stats.Centroid(1);
            y_c = stats.Centroid(2);
            crop_region = [2 * x_box - x_c, 2 * y_box - y_c, stats.BoundingBox(3:4) * 2];
            
            % Crop each channel to this specific object.
            % The original image is pre-masked with 'selected_object' to zero
            % out pixels from other condensates that might be in the crop box.
            Imc_BF = imcrop(selected_object.*im2double(Im_ch1), crop_region); % Use ch1 as a stand-in for Brightfield
            Imc_ch1 = imcrop(selected_object.*im2double(Im_ch1), crop_region);
            Imc_ch2 = imcrop(selected_object.*im2double(Im_ch2), crop_region);
            
            % --- Image Preprocessing on Cropped Object ---
            ImG_ch1 = imgaussfilt(Imc_ch1, sigma);
            ImG_ch2 = imgaussfilt(Imc_ch2, sigma);
            
            % Background flattening using morphological opening
            se = strel('disk', 150);
            ImFlat_ch1 = imadjust(ImG_ch1 - imopen(ImG_ch1, se));
            ImFlat_ch2 = imadjust(ImG_ch2 - imopen(ImG_ch2, se));
    
            % --- Generate Binary Masks for Each Phase ---
            level_ch1 = multithresh(ImFlat_ch1);
            Mask_ch1 = imclearborder(imquantize(ImFlat_ch1, level_ch1) - 1);
    
            level_ch2 = multithresh(ImFlat_ch2);
            Mask_ch2 = imclearborder(imquantize(ImFlat_ch2, level_ch2) - 1);
    
            % Create masks for non-overlapping regions ('m') and the background
            mch1 = Mask_ch1 - (Mask_ch1.*Mask_ch2); % ch1-only regions
            mch2 = Mask_ch2 - (Mask_ch1.*Mask_ch2); % ch2-only regions
            Mask_total = (Mask_ch1 + Mask_ch2) ~= 0; % Combined mask of the whole object
            Mask_background = imcomplement(Mask_total);
    
            % --- Calculate Intensity and Correlation Metrics ---
            sample_names = [sample_names, k]; % Store the object index
            
            % Pearson's correlation coefficient
            R_raw = corr2(Imc_ch1, Imc_ch2);
            r_raw = [r_raw, R_raw];
            r_ImFlat = [r_ImFlat, corr2(ImFlat_ch1, ImFlat_ch2)];
            corr2_Masks = [corr2_Masks, corr2(Mask_ch1, Mask_ch2)];
    
            % Calculate mean background intensity
            intensityBackground_ch1 = sum(sum(Imc_ch1.*Mask_background))/nnz(Mask_background);
            intensityBackground_ch2 = sum(sum(Imc_ch2.*Mask_background))/nnz(Mask_background);
            
            % Calculate intensity partitioning, using different mask definitions
            % depending on whether the object is considered mixed or demixed.
            if corr2(Mask_ch1, Mask_ch2) <= 0.91 % For demixed condensates
                IntensityMajor_ch1 = sum(sum(Imc_ch1.*mch1))/nnz(mch1); % ch1 intensity in ch1-rich area
                IntensityMinor_ch1 = sum(sum(Imc_ch1.*mch2))/nnz(mch2); % ch1 intensity in ch2-rich area
                IntensityMajor_ch2 = sum(sum(Imc_ch2.*mch2))/nnz(mch2); % ch2 intensity in ch2-rich area
                IntensityMinor_ch2 = sum(sum(Imc_ch2.*mch1))/nnz(mch1); % ch2 intensity in ch1-rich area
            else % For mixed condensates
                 % If masks overlap significantly, non-overlapping regions (mch1, mch2) may be empty or too small.
                 % This block uses the total mask for each channel to get a more stable intensity reading.
                IntensityMajor_ch1 = sum(sum(Imc_ch1.*Mask_ch1))/nnz(Mask_ch1);
                IntensityMinor_ch1 = sum(sum(Imc_ch1.*Mask_ch2))/nnz(Mask_ch2);
                IntensityMajor_ch2 = sum(sum(Imc_ch2.*Mask_ch2))/nnz(Mask_ch2);
                IntensityMinor_ch2 = sum(sum(Imc_ch2.*Mask_ch1))/nnz(Mask_ch1);
            end
            
            % Append partitioning ratios
            IntensityTotal_ch1 = IntensityMinor_ch1 + IntensityMajor_ch1;
            IntensityTotal_ch2 = IntensityMinor_ch2 + IntensityMajor_ch2;
            
            Mm_ch1 = [Mm_ch1, IntensityMajor_ch1 / IntensityMinor_ch1];
            Mm_ch2 = [Mm_ch2, IntensityMajor_ch2 / IntensityMinor_ch2];
    
            % Handle case where raw correlation is NaN
            if isnan(R_raw)
                Mm_ch1(end) = NaN;
                Mm_ch2(end) = NaN;
            end
    
            % Normalized partitioning coefficients
            NMm_ch1 = [NMm_ch1, IntensityMajor_ch1 / IntensityTotal_ch1];
            NMm_ch2 = [NMm_ch2, IntensityMajor_ch2 / IntensityTotal_ch2];
            NMm_ch1_back = [NMm_ch1_back, (IntensityMajor_ch1 - intensityBackground_ch1) / (IntensityTotal_ch1 - intensityBackground_ch1)];
            NMm_ch2_back = [NMm_ch2_back, (IntensityMajor_ch2 - intensityBackground_ch2) / (IntensityTotal_ch2 - intensityBackground_ch2)];
    
            % Calculate area in square micrometers
            A_ch1 = [A_ch1, nnz(Mask_ch1) * pixel_size^2 / 1e6];
            A_ch2 = [A_ch2, nnz(Mask_ch2) * pixel_size^2 / 1e6];
           
            % --- Contact Angle Analysis ---
            % This re-cropping step is preserved from the original script,
            % as it is essential to the workflow.
            Imc_BF = imcrop(selected_object.*im2double(Im_ch1), crop_region);
            Imc_ch1 = imcrop(selected_object.*im2double(Im_ch1), crop_region);
            Imc_ch2 = imcrop(selected_object.*im2double(Im_ch2), crop_region);
            Imc_fused = imfuse(Imc_ch1, Imc_ch2);
            
            % Attempt contact angle analysis only if object is de-mixed
            if corr2(Mask_ch1, Mask_ch2) <= 0.91
                try
                    % Call external function to get contact angle data
                    fittedcircles = ContactAnglesDT(Imc_ch1, Imc_ch2, sigma, k, P);
                    
                    % Rename columns for easier access
                    fittedcircles.Properties.VariableNames = {'SampleName', 'cAngles_ch1_1', 'cAngles_ch1_2', ...
                        'cAngles_ch2_1', 'cAngles_ch2_2', 'R_ch1', 'R_ch2', 'R_int', 'C_ch1_x', 'C_ch1_y', ...
                        'C_ch2_x', 'C_ch2_y', 'C_Int_x', 'C_Int_y', 'intT_A_AB', 'intT_B_AB'};
                    
                    % --- Manual Validation Loop ---
                    % For each interface found, plot fits and ask for user validation
                    for num_interface = 1:height(fittedcircles)
                        % Extract circle parameters for plotting
                        center_ch1 = [fittedcircles.C_ch1_x(num_interface), fittedcircles.C_ch1_y(num_interface)];
                        radius_ch1 = fittedcircles.R_ch1(num_interface);
                        center_ch2 = [fittedcircles.C_ch2_x(num_interface), fittedcircles.C_ch2_y(num_interface)];
                        radius_ch2 = fittedcircles.R_ch2(num_interface);
                        center_interface = [fittedcircles.C_Int_x(num_interface), fittedcircles.C_Int_y(num_interface)];
                        radius_interface = fittedcircles.R_int(num_interface);
                        
                        % Create validation figure with multiple subplots
                        figure('Name', 'Contact Angle Validation', 'WindowState', 'maximized');
                        
                        % Subplot 1: Fused image with fitted circles
                        subplot(2, 3, 1);
                        imshow(Imc_fused); hold on;
                        title(sprintf('Fitted Circles (Object %d)', k));
                        viscircles(center_ch1, radius_ch1, 'EdgeColor','g', 'LineStyle','--');
                        viscircles(center_ch2, radius_ch2, 'EdgeColor','m', 'LineStyle','--');
                        viscircles(center_interface, radius_interface, 'EdgeColor','w', 'LineStyle','--');
                        hold off;
            
                        % Subplots 2-6: Show the various masks used
                        subplot(2,3,2); imshow(mch1); title('ch1-only Mask');
                        subplot(2,3,5); imshow(mch2); title('ch2-only Mask');
                        subplot(2,3,3); imshow(Mask_ch1); title('Total ch1 Mask');
                        subplot(2,3,6); imshow(Mask_ch2); title('Total ch2 Mask');
                        
                        % Calculate and plot tangents at intersection points
                        subplot(2,3,4);
                        imshow(Imc_fused); hold on;
                        title('Fitted Tangents');
                        [xout,yout] = circcirc(center_ch1(1), center_ch1(2), radius_ch1, center_ch2(1),center_ch2(2), radius_ch2);
                        if ~any(isnan(xout))
                            x_range1 = linspace(xout(1)-40, xout(1)+40);
                            x_range2 = linspace(xout(2)-40, xout(2)+40);
                            m1_ch1 = -1/((center_ch1(2)-yout(1))/(center_ch1(1)-xout(1)));
                            m2_ch1 = -1/((center_ch1(2)-yout(2))/(center_ch1(1)-xout(2)));
                            y1_ch1 = m1_ch1 * (x_range1 - xout(1)) + yout(1);
                            y2_ch1 = m2_ch1 * (x_range2 - xout(2)) + yout(2);
                            m1_ch2 = -1/((center_ch2(2)-yout(1))/(center_ch2(1)-xout(1)));
                            m2_ch2 = -1/((center_ch2(2)-yout(2))/(center_ch2(1)-xout(2)));
                            y1_ch2 = m1_ch2 * (x_range1 - xout(1)) + yout(1);
                            y2_ch2 = m2_ch2 * (x_range2 - xout(2)) + yout(2);
                            m1_int = -1/((center_interface(2)-yout(1))/(center_interface(1)-xout(1)));
                            m2_int = -1/((center_interface(2)-yout(2))/(center_interface(1)-xout(2)));
                            y1_int = m1_int * (x_range1 - xout(1)) + yout(1);
                            y2_int = m2_int * (x_range2 - xout(2)) + yout(2);
                            
                            plot(x_range1, y1_ch1, 'g', 'LineWidth', 2);
                            plot(x_range1, y1_ch2, 'm', 'LineWidth', 2);
                            plot(x_range1, y1_int, 'w', 'LineWidth', 2);
                            plot(x_range2, y2_ch1, 'g', 'LineWidth', 2);
                            plot(x_range2, y2_ch2, 'm', 'LineWidth', 2);
                            plot(x_range2, y2_int, 'w', 'LineWidth', 2);
                        end
                        hold off;
                        
                        % Ask user to validate the fit
                        answer = questdlg('Is the fit correct?', 'Validate Fit', 'Yes', 'No', 'Yes');
                        
                        % Handle response
                        if strcmp(answer, 'Yes')
                            contactAngles = vertcat(contactAngles, fittedcircles(num_interface, :));
                            Table_ContactAngles = vertcat(Table_ContactAngles, contactAngles);
                            fprintf('Fit accepted. ');
                        else
                            % Placeholder command from original script
                            bla = 5; 
                            fprintf('Fit rejected by user. ');
                        end
                        close(gcf); % Close the validation figure
                   end
                catch ExceptionThrown
                    % Placeholder command from original script
                    bla = 5;
                    fprintf('Contact angle analysis failed for object %d: %s. ', k, ExceptionThrown.message);
                    close(findobj('type','figure')); % Close any open figures from the failed attempt           
                end
            else
                % Placeholder command from original script
                demixed = 0;
                fprintf('Object appears mixed (mask correlation > 0.91), skipping contact angle analysis. ');
            end
            fprintf('\n');
         end
    end
    
    % 5. DATA EXPORT
    % ========================================================================
    % Create a table with the intensity and correlation results
    intensity_table = table(sample_names', r_raw', r_ImFlat', corr2_Masks', ...
                            A_ch1', A_ch2', Mm_ch1', Mm_ch2', NMm_ch1', NMm_ch2', ...
                            NMm_ch1_back', NMm_ch2_back', ...
                            'VariableNames', {'SampleName', 'r_raw', 'r_ImFlat', 'corr2_Masks', ...
                            'Area_ch1', 'Area_ch2', 'Mm_ch1', 'Mm_ch2', 'NMm_ch1', 'NMm_ch2', ...
                            'NMm_ch1_background_subtracted', 'NMm_ch2_background_subtracted'});
    
    % Remove any rows with missing data and write to a CSV file
    intensity_table_clean = rmmissing(intensity_table);
    if ~isempty(intensity_table_clean)
        writetable(intensity_table_clean, [file_name '_rho_corr.csv']);
    end
    
    % Process and save the contact angle data if any was collected
    if height(Table_ContactAngles) > 0
        % Define column names (original script re-defines them here, but it's
        % already done in the validation loop. This is kept for logical consistency).
        Table_ContactAngles.Properties.VariableNames = {'SampleName', 'cAngles_ch1_1', 'cAngles_ch1_2', ...
            'cAngles_ch2_1', 'cAngles_ch2_2', 'R_ch1', 'R_ch2', 'R_int', 'C_ch1_x', 'C_ch1_y', ...
            'C_ch2_x', 'C_ch2_y', 'C_Int_x', 'C_Int_y', 'intT_A_AB', 'intT_B_AB'};
        
        % Remove rows with missing data
        CA_Table_clean = rmmissing(Table_ContactAngles);
        
        % Calculate average angles and other derived metrics
        CA_Table_clean.avg_ch1 = (CA_Table_clean.cAngles_ch1_1 + CA_Table_clean.cAngles_ch1_2)/2;
        CA_Table_clean.avg_ch2 = (CA_Table_clean.cAngles_ch2_1 + CA_Table_clean.cAngles_ch2_2)/2;
        CA_Table_clean.theta = 360 - CA_Table_clean.avg_ch1 - CA_Table_clean.avg_ch2;
        CA_Table_clean.cos_theta2 = cosd(CA_Table_clean.theta/2);
        
        % Write the unique, clean results to a CSV file
        if ~isempty(CA_Table_clean)
            writetable(unique(CA_Table_clean), [file_name '_CA.csv']);
        end

        % Create and save a final image showing only the condensates that
        % were successfully analyzed for contact angles.
        figure('Name', 'Final Analyzed Condensates');
        imshow(Im_fused);
        hold on;
        title(['Successfully Analyzed Condensates in: ', strrep(current_filename, '_', '\_')]);
        selected_condensates = unique(CA_Table_clean.SampleName)';
        for k_sel = selected_condensates
            plot(B{k_sel}(:,2), B{k_sel}(:,1), 'w', 'LineWidth', 2);
        end
        exportgraphics(gcf, ["CA_Selected_condensates_" + file_name +  ".png"], 'Resolution', 300);
        hold off;
        close(gcf);
    end
end
fprintf('\n--- Analysis Complete ---\n');