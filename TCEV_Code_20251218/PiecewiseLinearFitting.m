function [parametrosfin, xybp] = PiecewiseLinearFitting(bpmin, varargin)
% This function works for any type of x-values, integer or non-integer
% Ensures both segments divided by the breakpoint contain at least 4 data points

% INPUTS:
% bpmin, Minimum segment length between breakpoints
% varargin, y Dependent variable

% OUTPUTS:
% parametrosfin, a1,a2,c1,c2 %y1=a1*x+c1; y2=a2*x+c2
% xybp, (x,y) Breakpoint coordinates

%%
% Data input and error checking module
if nargin == 2
    y = varargin{1};
    x = (1:length(y))';
elseif nargin > 3
    fprintf('Number of arguments: %d\n', nargin)
    error('Too much input arguments')
elseif nargin == 3
    x = varargin{1}; 
    y = varargin{2};
    if length(x) - length(y) ~= 0
        error('x and y do not have the same length');
    end
end

% Remove NaN values from the data series
ind = find(isnan(y)); % First handle dependent variable
y(ind) = []; 
x(ind) = [];
ind = find(isnan(x)); % Then handle independent variable
y(ind) = []; 
x(ind) = [];
clear ind

% Check for insufficient data, require at least (2*bpmin+1) points
if length(x) < (2 * bpmin + 1)
    error('Too short data series'); 
end

% Sort by x in ascending order
xy = [x y];
xy = sortrows(xy, 1);

%% Piecewise Linear Fitting
% Get initial values from first and last points
mdl1 = fitlm(xy(1:bpmin, 1), xy(1:bpmin, 2));
mdl2 = fitlm(xy(end-bpmin+1:end, 1), xy(end-bpmin+1:end, 2));

% y1=a1*x+c1; y2=a2*x+c2
parini(1, 1) = mdl1.Coefficients.Estimate(2); % a1 
parini(2, 1) = mdl2.Coefficients.Estimate(2); % a2 
parini(3, 1) = mdl1.Coefficients.Estimate(1); % c1 
parini(4, 1) = quantile(xy(:, 1), 0.53); 


% Optimization using fminsearch
[parfin, opt_error] = fminsearch(@PiecewiseLFv2ErrorRMSE, parini, [], xy);

% Ensure both segments divided by the breakpoint have sufficient data points
min_points = 4; 
x_data = xy(:, 1);
bp_optimized = parfin(4, 1); % Optimized breakpoint

% Calculate number of data points on each side of the optimized breakpoint
left_points_opt = sum(x_data <= bp_optimized);
right_points_opt = sum(x_data > bp_optimized);

% If either side has insufficient points, adjust breakpoint position
bp_adjusted = bp_optimized; % Initialize adjusted breakpoint
adjustment_made = false; % Flag to track if adjustment was made

if left_points_opt < min_points
    % Find the min_points-th smallest x value as the new breakpoint
    sorted_x = sort(x_data);
    bp_adjusted = sorted_x(min_points) + eps; % Slight right offset to ensure enough points
    fprintf('Adjusting breakpoint position: insufficient left points (%d < %d)\n', left_points_opt, min_points);
    adjustment_made = true;
elseif right_points_opt < min_points
    % Find the min_points-th largest x value as the new breakpoint
    sorted_x = sort(x_data);
    bp_adjusted = sorted_x(end - min_points + 1) - eps; % Slight left offset to ensure enough points
    fprintf('Adjusting breakpoint position: insufficient right points (%d < %d)\n', right_points_opt, min_points);
    adjustment_made = true;
end

% Use the adjusted breakpoint
if adjustment_made
    parfin(4, 1) = bp_adjusted;
    fprintf('Breakpoint adjusted from %.4f to %.4f\n', bp_optimized, bp_adjusted);
end

% Recalculate number of points on each side after adjustment
left_points_adj = sum(x_data <= parfin(4, 1));
right_points_adj = sum(x_data > parfin(4, 1));

% Calculate final parameters
parametrosfin = parfin(1:3); % (a1, a2, c1)
parametrosfin(4, 1) = parfin(3, 1) + (parfin(1, 1) - parfin(2, 1)) * parfin(4, 1); % c2

% Calculate breakpoint coordinates (using adjusted breakpoint)
xybp(1, 1) = parfin(4, 1);
xybp(2, 1) = parfin(1, 1) * parfin(4, 1) + parfin(3, 1);

%% Calculate fitted values (using adjusted breakpoint)
yajus = nan(size(xy, 1), 1);
for i = 1:length(xy)
    if xy(i, 1) <= parfin(4, 1)
        yajus(i, 1) = parfin(1, 1) * xy(i, 1) + parfin(3, 1); % First segment fit
    else
        yajus(i, 1) = parfin(1, 1) * parfin(4, 1) + parfin(2, 1) * (xy(i, 1) - parfin(4, 1)) + parfin(3, 1); % Second segment fit
    end
end

%% Plotting (using adjusted breakpoint)
dibuja1 = 0; % Force plotting on
if dibuja1 == 1
    fig = figure('Units', 'inches', 'Position', [2 0 5 3.8], ...
                 'Color', 'white', ...
                 'PaperPositionMode', 'auto');
             
    set(gca, 'FontName', 'Arial'); % Set global font to Arial
    
    ax = gca;
    tightInset = get(ax, 'TightInset');
    newPos = [tightInset(1) + 0.05, tightInset(2) + 0.07, ...
              1 - tightInset(1) - tightInset(3) - 0.1, ...
              1 - tightInset(2) - tightInset(4) - 0.1];
    set(ax, 'Position', newPos);

    % Plot original data points
    scatter(xy(:, 1), xy(:, 2), 30, 'd', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.5 0.5 0.5], ...
         'DisplayName', 'Observed Data');
    hold on;
    
    % Generate theoretical fitted curve (using adjusted breakpoint)
    xteo = linspace(min(xy(:, 1)), max(xy(:, 1)), 200)';
    yteo = nan(size(xteo));
    for i = 1:length(xteo)
        if xteo(i) <= parfin(4, 1)
            yteo(i) = parfin(1, 1) * xteo(i) + parfin(3, 1); % part 1
        else
            yteo(i) = parfin(1, 1) * parfin(4, 1) + parfin(2, 1) * (xteo(i) - parfin(4, 1)) + parfin(3, 1); % part 2
        end
    end
    
    % Plot the two line segments separately (different colors, using adjusted breakpoint)
    % First line segment
    x1 = linspace(min(xy(:, 1)), parfin(4, 1), 100);
    y1 = parfin(1, 1) * x1 + parfin(3, 1);
    plot(x1, y1, 'r--', 'LineWidth', 1, 'DisplayName', 'Segment 1');
    
    % Second line segment
    x2 = linspace(parfin(4, 1), max(xy(:, 1)), 100);
    y2 = parfin(2, 1) * (x2 - parfin(4, 1)) + (parfin(1, 1) * parfin(4, 1) + parfin(3, 1));
    plot(x2, y2, 'b--', 'LineWidth', 1, 'DisplayName', 'Segment 2');
    
    % Mark the breakpoint (using adjusted breakpoint)
    scatter(xybp(1, 1), xybp(2, 1), 150, 'o', ...
        'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceAlpha', 0.1, ...
        'SizeData', 500, 'DisplayName', 'Break Point');

    % Add line equation annotations
    eq1 = sprintf('y = %.5fx %s %.3f', parfin(1, 1), ...
        signChar(parfin(3, 1)), abs(parfin(3, 1)));
    eq2 = sprintf('y = %.5fx %s %.3f', parfin(2, 1), ...
        signChar((parfin(1, 1) - parfin(2, 1)) * parfin(4, 1) + parfin(3, 1)), ...
        abs((parfin(1, 1) - parfin(2, 1)) * parfin(4, 1) + parfin(3, 1)));
    
    % Set base vertical position and line spacing
    base_y = 0.95;
    line_spacing = 0.08;
    
    text(0.05, base_y, eq1, 'Units', 'normalized', ...
        'Color', 'r', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Arial');
    text(0.05, base_y - line_spacing, eq2, 'Units', 'normalized', ...
        'Color', 'b', 'FontSize', 15, 'FontWeight', 'normal', 'FontName', 'Arial');
 
    % Axis style settings
    set(ax, 'FontSize', 12, ...
            'Box', 'off', ...
            'TickDir', 'out', ...
            'XColor', [0.5 0.5 0.5], ...
            'YColor', [0.5 0.5 0.5], ...
            'LineWidth', 1);

    xlim(ax, [0 max(xy(:, 1))]);
    set(ax, 'XTickMode', 'auto');

    % Get ticks and re-add black text
    xticks = get(ax, 'XTick');
    yticks = get(ax, 'YTick');
    xticklabels = get(ax, 'XTickLabel');
    yticklabels = get(ax, 'YTickLabel');

    set(ax, 'XTickLabel', [], 'YTickLabel', []);

    x_offset = 0;
    y_offset = -0.1;
    y_text_offset = -3;

    % X-axis tick labels
    for i = 1:length(xticks)
        text(xticks(i) + x_offset, ax.YLim(1) + y_offset, xticklabels{i}, ...
            'Color', 'k', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'FontSize', 12);
    end

    % Y-axis tick labels
    for i = 1:length(yticks)
        text(ax.XLim(1) + y_text_offset, yticks(i), yticklabels{i}, ...
            'Color', 'k', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'FontSize', 12);
    end

    % Set tick length
    ax.XAxis.TickLength = [0.005 0.005];
    ax.YAxis.TickLength = [0.005 0.005];
    ax.Layer = 'top';

    xlabel('Extreme precipitation', 'FontWeight', 'bold', 'Interpreter', 'tex', ...
           'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'top');
    hx = get(gca, 'XLabel');
    set(hx, 'Units', 'pixels');
    pos_x = get(hx, 'Position');
    set(hx, 'Position', [pos_x(1) pos_x(2) - 10 pos_x(3)], 'Units', 'normalized');

    % Y-axis label
    ylabel('-ln(ln(F(x)))', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', ...
           'HorizontalAlignment', 'center');
    hy = get(gca, 'YLabel');
    set(hy, 'Units', 'pixels');
    pos_y = get(hy, 'Position');
    set(hy, 'Position', [pos_y(1) - 14 pos_y(2) pos_y(3)], 'Units', 'normalized');
    
    % Create legend
    lgd = legend('show');
    lgd.Position = [0.6, 0.25, 0.2, 0.1];
    set(lgd, 'FontSize', 12);
    
    % Add gray border line that precisely matches canvas edges
    annotation(fig, 'rectangle', ...
               [0 0 1 1], ...
               'Color', [0.5 0.5 0.5], ...
               'LineWidth', 1.5, ...
               'EdgeColor', [0.5 0.5 0.5]);

    set(gca, 'FontSize', 12, 'FontName', 'Arial');
    
    % Save figure
    print(fig, 'breakpointPiecewiseLinearFit_Plot.png', '-dpng', '-r600', '-painters');
    
    hold off;
end

end

%% Helper function: sign character handling
function s = signChar(value)
    if value >= 0
        s = '+';
    else
        s = '-';
    end
end

%% Error calculation function
function error = PiecewiseLFv2ErrorRMSE(par, xy)
    % Calculate RMSE error for piecewise linear fitting
    
    a1 = par(1);
    a2 = par(2);
    c1 = par(3);
    bp = par(4);
    
    y_pred = zeros(size(xy, 1), 1);
    
    for i = 1:size(xy, 1)
        if xy(i, 1) <= bp
            y_pred(i) = a1 * xy(i, 1) + c1;
        else
            y_pred(i) = a1 * bp + a2 * (xy(i, 1) - bp) + c1;
        end
    end
    
    % Calculate RMSE
    error = sqrt(mean((xy(:, 2) - y_pred).^2));
end