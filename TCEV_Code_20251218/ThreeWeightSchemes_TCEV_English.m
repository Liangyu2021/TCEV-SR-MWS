clear; clc; close all

% Script identifier (for unified naming)
scriptName = 'Case_Parameters';

% Create output directories
if ~exist('OutputData', 'dir')
    mkdir('OutputData')
end
if ~exist(fullfile('OutputData', 'WeightComparisons'), 'dir')
    mkdir(fullfile('OutputData', 'WeightComparisons'))
end

% Read CSV file (compatible with null and missing values)
fid = fopen('EightCasesRawData.csv');
if fid == -1
    error('Cannot open file EightCasesRawData.csv, please ensure the file exists');
end
rawText = textscan(fid, '%s', 'Delimiter', '\n'); % Read by line
fclose(fid);
rawLines = rawText{1};

% Parse data (handle missing values)
maxColumns = 8; % Set maximum columns based on data
data = nan(length(rawLines)-1, maxColumns); % Skip header row

for i = 2:length(rawLines) % Start from row 2 (skip header)
    lineData = textscan(rawLines{i}, '%s', 'Delimiter', ',');
    tokens = lineData{1};
    for j = 1:min(length(tokens), maxColumns)
        if ~isempty(tokens{j}) && ~strcmp(tokens{j}, '-')
            data(i-1,j) = str2double(tokens{j});
        end
    end
end

caseIDs = {'1', '2', '3', '4', '5', '6', '7', '8'}; % Define case IDs
maxnCases = size(data, 2); % Total number of columns

% Define basic weight methods (uniform 4-column cell array)
allWeightMethods = {
    'linear', @(y, cdfCol) y(:,cdfCol), 'Linear weight'; 
    'quadratic',  @(y, cdfCol) y(:,cdfCol).^2, 'Quadratic weight';
    'exponential',  @(y, cdfCol) exp(y(:,cdfCol)), 'Exponential weight'
};

% Initialize result storage
allResults = cell(length(allWeightMethods), 1);
allFittedValues = cell(length(allWeightMethods), 1);
allInitialParams = cell(length(allWeightMethods), 1);
allPerformanceMetrics = cell(length(allWeightMethods), 1);
allIntervalMetrics = cell(length(allWeightMethods), 1); 

fprintf('Starting weight method analysis, total %d weight schemes...\n', length(allWeightMethods));

% Define F(x) intervals
cdf_intervals = [0, 0.2; 0.2, 0.4; 0.4, 0.6; 0.6, 0.8; 0.8, 1.0; 0, 1.0];
interval_names = {'0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0', '0-1.0'};

for w = 1:length(allWeightMethods)
    methodName = allWeightMethods{w,1};
    weightFunc = allWeightMethods{w,2};
    methodDescription = allWeightMethods{w,3};
    
    fprintf('Processing weight method %d/%d: %s (%s)\n', w, length(allWeightMethods), methodName, methodDescription);
    
    % Initialize current method's result storage
    results = cell2table(cell(0,11), 'VariableNames', {...
        'CaseID', 'alpha1', 'beta1', 'alpha2', 'beta2', 'WeightMethod', ...
        'alpha1_initial', 'beta1_initial', 'alpha2_initial', 'beta2_initial', 'Breakpoint_X'});
    
    performanceMetrics = cell2table(cell(0,7), 'VariableNames', {...
        'CaseID', 'WeightMethod', 'Total_RMSE', 'Breakpoint_X', 'Breakpoint_Y', ...
        'Initial_Error', 'Final_Error'});
    
    intervalMetricsAll = []; % Store interval performance metrics
    fittedValuesAll = [];
    initialParamsAll = [];
    
    for ncase = 1:maxnCases
        % Extract current case data (skip NaN)
        muestra = data(:, ncase);
        validIdx = ~isnan(muestra);
        muestra = muestra(validIdx);
        
        if isempty(muestra) || length(muestra) < 5
            warning('Case %s has no valid data or insufficient data, skipping', caseIDs{ncase});
            continue;
        end
        
        muestra = sort(muestra);
        
        % Probability calculation
        temp = (1:length(muestra))';
        yobs = [muestra, temp/(length(temp)+1), (temp-0.44)/(length(temp)+0.12)];
        eleccCDF = 2; % Use second probability estimation method
        
        % Initial parameter estimation (with error handling)
        try
            bpmin = 3;
            [parAjusteLineal, xybp] = PiecewiseLinearFitting(bpmin, yobs(:,1), -log(-log(yobs(:,eleccCDF))));
            parini = [parAjusteLineal(1), -parAjusteLineal(3)/parAjusteLineal(1), ...
                     parAjusteLineal(2), -parAjusteLineal(4)/parAjusteLineal(2)];
            
            % Calculate PosRotura
            PosRotura = find(yobs(:,1) > xybp(1,1),1)-1;
            if isempty(PosRotura) || PosRotura < 2 || PosRotura > length(yobs)-2
                PosRotura = max(2, min(length(yobs)-2, floor(length(yobs)*0.7)));
            end
        catch
            warning('Case %s initial parameter estimation failed, using default parameters', caseIDs{ncase});
            parini = [1, 1, 1, 1];
            xybp = [median(yobs(:,1)), 0.5];
            PosRotura = max(2, min(length(yobs)-2, floor(length(yobs)*0.7)));
        end
        
        % Calculate initial error
        ycalc_initial = TCEVcdfv2(yobs(:,1), parini);
        initial_error = sqrt(mean((ycalc_initial - yobs(:,eleccCDF)).^2));
        
        % Weight calculation
        weights = weightFunc(yobs, eleccCDF);
        
        % Normalize weights (maintain constant weight sum)
        weights = weights / mean(weights);
        
        % Weighted fitting
        options = optimset('MaxFunEvals', 2000, 'Display', 'off');
        try
            [par1, final_error] = fminsearch(@(p) F01error5V3_weighted(p, 1:4, [], yobs(:,1), yobs(:,eleccCDF), @TCEVcdfv2, PosRotura, weights), parini, options);
        catch
            warning('Case %s optimization failed, using initial parameters', caseIDs{ncase});
            par1 = parini;
            final_error = initial_error;
        end
        
        % Calculate final fitted values
        ycalc_final = TCEVcdfv2(yobs(:,1), par1);
        
        % Calculate total RMSE
        total_rmse = sqrt(mean((ycalc_final - yobs(:,eleccCDF)).^2));
        
        % Calculate RMSE for each F(x) interval
        interval_rmse = zeros(length(cdf_intervals), 1);
        for i = 1:length(cdf_intervals)
            interval_mask = (yobs(:,eleccCDF) >= cdf_intervals(i,1)) & (yobs(:,eleccCDF) <= cdf_intervals(i,2));
            if sum(interval_mask) > 0
                interval_rmse(i) = sqrt(mean((ycalc_final(interval_mask) - yobs(interval_mask,eleccCDF)).^2));
            else
                interval_rmse(i) = NaN;
            end
        end
        
        % Store performance metrics
        performanceMetrics = [performanceMetrics; {
            caseIDs{ncase}, methodName, total_rmse, ...
            xybp(1,1), xybp(2,1), initial_error, final_error
        }];
        
        % Store interval performance metrics
        for i = 1:length(cdf_intervals)
            intervalRow = table(...
                {caseIDs{ncase}}, {methodName}, interval_rmse(i), ...
                cdf_intervals(i,1), cdf_intervals(i,2), {interval_names{i}}, ...
                'VariableNames', {'CaseID', 'WeightMethod', 'Interval_RMSE', ...
                'CDF_Start', 'CDF_End', 'Interval_Name'});
            
            if isempty(intervalMetricsAll)
                intervalMetricsAll = intervalRow;
            else
                intervalMetricsAll = [intervalMetricsAll; intervalRow];
            end
        end
        
        % Store parameters
        results = [results; {
            caseIDs{ncase}, par1(1), par1(2), par1(3), par1(4), methodName, ...
            parini(1), parini(2), parini(3), parini(4), xybp(1,1)
        }];
        
        % Calculate fitted values
        xsim = linspace(min(muestra), max(muestra), 100)';
        ysim1 = TCEVcdfv2(xsim, par1);
        tempTable = table(...
            repmat(caseIDs(ncase), 100, 1), xsim, ysim1, ...
            repmat({methodName}, 100, 1), ...
            'VariableNames', {'CaseID', 'X', 'FittedCDF', 'WeightMethod'});
        
        if isempty(fittedValuesAll)
            fittedValuesAll = tempTable;
        else
            fittedValuesAll = [fittedValuesAll; tempTable];
        end
    end
    
    % Store current method's results
    allResults{w} = results;
    allFittedValues{w} = fittedValuesAll;
    allInitialParams{w} = initialParamsAll;
    allPerformanceMetrics{w} = performanceMetrics;
    allIntervalMetrics{w} = intervalMetricsAll;
end

%% Combine all results and save
% Combine parameter results
combinedResults = table();
for w = 1:length(allWeightMethods)
    if ~isempty(allResults{w})
        combinedResults = [combinedResults; allResults{w}];
    end
end
if ~isempty(combinedResults)
    writetable(combinedResults, fullfile('OutputData', 'WeightComparisons', sprintf('%s_AllParameters.csv', scriptName)));
end

% Combine performance metrics
combinedMetrics = table();
for w = 1:length(allWeightMethods)
    if ~isempty(allPerformanceMetrics{w})
        combinedMetrics = [combinedMetrics; allPerformanceMetrics{w}];
    end
end
if ~isempty(combinedMetrics)
    writetable(combinedMetrics, fullfile('OutputData', 'WeightComparisons', sprintf('%s_AllPerformanceMetrics.csv', scriptName)));
end

% Combine interval performance metrics
combinedIntervalMetrics = table();
for w = 1:length(allWeightMethods)
    if ~isempty(allIntervalMetrics{w})
        combinedIntervalMetrics = [combinedIntervalMetrics; allIntervalMetrics{w}];
    end
end
if ~isempty(combinedIntervalMetrics)
    writetable(combinedIntervalMetrics, fullfile('OutputData', 'WeightComparisons', sprintf('%s_AllIntervalMetrics.csv', scriptName)));
end

method_stats = [];
method_names = {};

for w = 1:length(allWeightMethods)
    if ~isempty(allPerformanceMetrics{w}) && height(allPerformanceMetrics{w}) > 0
        metrics = allPerformanceMetrics{w};
        avg_total_rmse = mean(metrics.Total_RMSE);
        
        method_names{end+1} = allWeightMethods{w,1};
        method_stats(end+1) = avg_total_rmse;
        
        fprintf('%-30s %-10.4f %-8s\n', ...
            allWeightMethods{w,1}, avg_total_rmse);
    end
end

%% Generate F(x) interval heatmaps for each case
fprintf('\n=== Generating F(x) interval heatmaps for each case (including reference methods) ===\n');

if ~isempty(combinedIntervalMetrics)
    % Read reference data
    refData = [];
    if exist('ReferenceCases.csv', 'file')
        refData = readtable('ReferenceCases.csv');
        fprintf('Successfully read reference data, %d rows %d columns\n', size(refData, 1), size(refData, 2));
    else
        fprintf('ReferenceCases.csv file not found, will use only weight method data\n');
    end
    
    % Get all case IDs
    unique_cases = unique(combinedIntervalMetrics.CaseID);
    
    % Define weight method order
    method_order = {
        'linear', 'quadratic', 'exponential'
    };
    
    for c = 1:length(unique_cases)
        caseID = unique_cases{c};
        fprintf('Generating F(x) interval heatmap for case %s (including references)...\n', caseID);
        
        % Extract current case data
        case_data = combinedIntervalMetrics(strcmp(combinedIntervalMetrics.CaseID, caseID), :);
        
        if height(case_data) < 2
            fprintf('Case %s has insufficient data, skipping heatmap generation\n', caseID);
            continue;
        end
        
        % Process reference data
        ref_methods = {};
        ref_display_names = {};
        ref_interval_data = [];
        
        if ~isempty(refData)
            % Find columns containing current case ID
            case_columns = [];
            for i = 1:width(refData)
                colName = refData.Properties.VariableNames{i};
                if contains(colName, caseID) || contains(colName, 'Case')
                    case_columns = [case_columns, i];
                end
            end
            
            % Extract reference F(x) values
            if ~isempty(case_columns)
                % Get observed data for RMSE calculation
                muestra = data(:, str2double(caseID));
                validIdx = ~isnan(muestra);
                muestra = muestra(validIdx);
                muestra = sort(muestra);
                
                % Probability calculation
                temp = (1:length(muestra))';
                yobs = [muestra, temp/(length(temp)+1), (temp-0.44)/(length(temp)+0.12)];
                eleccCDF = 2;
                
                for col_idx = 1:length(case_columns)
                    col_num = case_columns(col_idx);
                    ref_Fx = refData{:, col_num};
                    
                    % Remove NaN values
                    valid_ref = ~isnan(ref_Fx);
                    if sum(valid_ref) == length(yobs)
                        % Calculate RMSE for references in each interval
                        ref_rmse = zeros(6, 1);
                        for i = 1:length(cdf_intervals)
                            interval_mask = (yobs(:,eleccCDF) >= cdf_intervals(i,1)) & ...
                                           (yobs(:,eleccCDF) <= cdf_intervals(i,2));
                            if sum(interval_mask) > 0
                                ref_rmse(i) = sqrt(mean((ref_Fx(interval_mask) - yobs(interval_mask,eleccCDF)).^2));
                            else
                                ref_rmse(i) = NaN;
                            end
                        end
                        
                        % Store reference data
                        ref_method_name = sprintf('Ref%d', col_idx);
                        ref_display_name = sprintf('Reference%d', col_idx);
                        
                        ref_methods{end+1} = ref_method_name;
                        ref_display_names{end+1} = ref_display_name;
                        
                        % Create reference interval data table
                        for i = 1:length(cdf_intervals)
                            ref_row = table(...
                                {caseID}, {ref_method_name}, ref_rmse(i), ...
                                cdf_intervals(i,1), cdf_intervals(i,2), {interval_names{i}}, ...
                                'VariableNames', {'CaseID', 'WeightMethod', 'Interval_RMSE', ...
                                'CDF_Start', 'CDF_End', 'Interval_Name'});
                            
                            if isempty(ref_interval_data)
                                ref_interval_data = ref_row;
                            else
                                ref_interval_data = [ref_interval_data; ref_row];
                            end
                        end
                    end
                end
            end
        end
        
        % Create heatmap data matrix - references at top, then weight methods
        methods_available = unique(case_data.WeightMethod);
        
        % Filter and sort weight methods according to specified order
        ordered_methods = {};
        ordered_display_names = {};
        for i = 1:length(method_order)
            if any(strcmp(methods_available, method_order{i}))
                ordered_methods{end+1} = method_order{i};
                ordered_display_names{end+1} = method_order{i};
            end
        end
        
        % Define interval order (0-1.0 at the end)
        intervals_ordered = {'0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0', '0-1.0'};
        
        % Initialize data matrix
        all_methods = [ref_methods, ordered_methods];
        all_display_names = [ref_display_names, ordered_display_names];
        data_matrix = zeros(length(all_methods), length(intervals_ordered));
        
        % Fill reference data
        for i = 1:length(ref_methods)
            for j = 1:length(intervals_ordered)
                if ~isempty(ref_interval_data)
                    mask = strcmp(ref_interval_data.WeightMethod, ref_methods{i}) & ...
                           strcmp(ref_interval_data.Interval_Name, intervals_ordered{j});
                    if any(mask)
                        data_matrix(i, j) = ref_interval_data.Interval_RMSE(mask);
                    else
                        data_matrix(i, j) = NaN;
                    end
                end
            end
        end
        
        % Fill weight method data
        for i = 1:length(ordered_methods)
            row_idx = length(ref_methods) + i;
            for j = 1:length(intervals_ordered)
                mask = strcmp(case_data.WeightMethod, ordered_methods{i}) & ...
                       strcmp(case_data.Interval_Name, intervals_ordered{j});
                if any(mask)
                    data_matrix(row_idx, j) = case_data.Interval_RMSE(mask);
                else
                    data_matrix(row_idx, j) = NaN;
                end
            end
        end
        
        % Create heatmap (multiply by 100 for display)
        fig = figure('Position', [100, 100, 1000, 800]);
        
        % Heatmap data multiplied by 100
        data_matrix_display = data_matrix * 100;
        
        imagesc(data_matrix_display);
        
        % Set axis labels
        set(gca, 'XTick', 1:length(intervals_ordered), 'XTickLabel', intervals_ordered);
        set(gca, 'YTick', 1:length(all_display_names), 'YTickLabel', all_display_names);
        set(gca, 'FontName', 'Arial', 'FontSize', 24,'FontWeight', 'bold','FontName' , 'Arial');
        
        % Add colorbar and set labels
        c = colorbar;
        c.Label.String = 'RMSE (¡Á10^{-2})';
        c.Label.FontSize = 24;
        c.Label.FontWeight = 'bold';
        c.Label.FontName = 'Arial';
        % Set colorbar to show only 5 tick labels
        c.Ticks = linspace(c.Limits(1), c.Limits(2), 5);
        c.TickLabels = compose('%d', round(c.Ticks));
              
        title(sprintf('Case %s', caseID), ...
              'FontSize', 24, 'FontWeight', 'bold','FontName' , 'Arial');
        xlabel('F(x)','FontSize', 24, 'FontWeight', 'bold','FontName' , 'Arial');
        
        % Add numerical text on heatmap (display values multiplied by 100)
        for i = 1:length(all_methods)
            for j = 1:length(intervals_ordered)
                if ~isnan(data_matrix_display(i, j))
                    text(j, i, sprintf('%.2f', data_matrix_display(i, j)), ...
                         'HorizontalAlignment', 'center', 'FontSize', 24, ...
                         'FontWeight', 'bold', 'Color', ifelse(data_matrix_display(i, j) > max(data_matrix_display(:))/2, 'white', 'black'));
                end
            end
        end
        
         % Add gray border line that precisely matches canvas edges
        annotation(fig, 'rectangle', ...
                   [0 0 1 1], ...
                   'Color', [0.5 0.5 0.5], ...
                   'LineWidth', 1.5, ...
                   'EdgeColor', [0.5 0.5 0.5]);
    
        % Save heatmap
        heatmap_filename = fullfile('OutputData', 'WeightComparisons', ...
            sprintf('%s_Case_%s_FxInterval_Heatmap_with_References.png', scriptName, caseID));
        print(gcf, heatmap_filename, '-dpng', '-r900');
        close(gcf);
        
        fprintf('Case %s F(x) interval heatmap saved: %s\n', caseID, heatmap_filename);        
    end
end

%% Helper function
function result = ifelse(condition, trueValue, falseValue)
    if condition
        result = trueValue;
    else
        result = falseValue;
    end
end

%% Weighted error function
function error = F01error5V3_weighted(par, indicespar, pesospar, xobs, yobs, func, PosRotura, weights)
% Weighted piecewise error function
    ycalc = func(xobs, par);
    
    if ~isempty(PosRotura) && PosRotura > 0 && PosRotura < length(yobs)
        error = sqrt(mean((ycalc - yobs).^2 .* weights));
    end
end