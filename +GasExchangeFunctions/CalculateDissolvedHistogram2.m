function HistFigure = CalculateDissolvedHistogram2(Data, Edges, Thresholds, CMap, HealthyEdges, HealthyFit, SNR, Mean, Median, SD, Bins, Name)
% Generate histogram figure with labeled statistics and bin percentages

% Create figure
HistFigure = figure('Name', 'Histogram'); 
set(HistFigure, 'WindowState', 'minimized');
set(HistFigure, 'Color', 'white', 'Units', 'inches', 'Position', [0.25 0.25 6 4]);
hold on;

% Plot histogram
histogram(Data, Edges, 'Normalization', 'probability', ...
    'EdgeColor', [0.5 0.5 0.5], 'FaceColor', [1 1 1], 'FaceAlpha', 0.75);

% Plot healthy reference fit
if exist('HealthyEdges', 'var') && exist('HealthyFit', 'var')
    plot(HealthyEdges, abs(HealthyFit) / size(Edges,2) * size(HealthyEdges,2), ...
        'k--', 'LineWidth', 2);
end

% Draw bin threshold backgrounds
ylims = ylim;
xlims = xlim;

% First region
pat = patch([xlims(1) xlims(1) Thresholds(1) Thresholds(1)], ...
            [0 ylims(2) ylims(2) 0], CMap(1,:), 'HandleVisibility', 'off');
pat.LineStyle = 'none'; pat.FaceVertexAlphaData = 0.1; pat.FaceAlpha = 'flat';
uistack(pat, 'bottom');

% Intermediate bins
for bin = 2:length(Thresholds)
    pat = patch([Thresholds(bin-1) Thresholds(bin-1) Thresholds(bin) Thresholds(bin)], ...
                [0 ylims(2) ylims(2) 0], CMap(bin,:), 'HandleVisibility', 'off');
    pat.LineStyle = 'none'; pat.FaceVertexAlphaData = 0.1; pat.FaceAlpha = 'flat';
    uistack(pat, 'bottom');
end

% Last region
pat = patch([Thresholds(end) Thresholds(end) xlims(2) xlims(2)], ...
            [0 ylims(2) ylims(2) 0], CMap(end,:), 'HandleVisibility', 'off');
pat.LineStyle = 'none'; pat.FaceVertexAlphaData = 0.1; pat.FaceAlpha = 'flat';
uistack(pat, 'bottom');

% Set axis limits and style
axis([Edges(2) Edges(end-1) -inf inf]);
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'YTickLabel', []);  % Hide y-axis tick labels

% ---------- Add bin percentages as title ----------
% Format percentages as string
binText = sprintf('%2.1f%% | ', Bins);
binText = binText(1:end-3);  % Remove last delimiter

% Truncate title if too long
if length(binText) > 80
    binText = [binText(1:77) '...'];
end
title(binText, 'FontSize', 18, 'FontWeight', 'bold');

% ---------- Add text labels ----------
% Choose left or right side
% x_text = 0;% 
if strcmp(Name,'Ventilation')
    x_text = 0;
    scaleFactor = 1;
    % Plot SNR, Mean, Median, SD below Name
    stats_text = {
        sprintf('SNR     : %.2f', SNR), ...
        sprintf('Mean    : %.2f', Mean*scaleFactor), ...
        sprintf('Median  : %.2f', Median*scaleFactor), ...
        sprintf('SD      : %.2f', SD*scaleFactor)
    };    
elseif strcmp(Name,'Membrane')
    x_text = Thresholds(5);
    scaleFactor = 1000;
    % Plot SNR, Mean, Median, SD below Name
    stats_text = {
        sprintf('SNR     : %.2f', SNR), ...
        sprintf('Mean    : %.1fx10^-^3', Mean*scaleFactor), ...
        sprintf('Median  : %.1fx10^-^3', Median*scaleFactor), ...
        sprintf('SD      : %.1fx10^-^3', SD*scaleFactor)
    };    
elseif strcmp(Name,'RBC')
    x_text = Thresholds(4);
    scaleFactor = 1000;
    % Plot SNR, Mean, Median, SD below Name
    stats_text = {
        sprintf('SNR     : %.2f', SNR), ...
        sprintf('Mean    : %.1fx10^-^3', Mean*scaleFactor), ...
        sprintf('Median  : %.1fx10^-^3', Median*scaleFactor), ...
        sprintf('SD      : %.1fx10^-^3', SD*scaleFactor)
    };    
elseif strcmp(Name,'RBC:Membrane')
    x_text = 0;% x_text = Thresholds(3);
    scaleFactor = 1;
    % Plot SNR, Mean, Median, SD below Name
    stats_text = {
        sprintf('SNR     : %.2f', SNR), ...
        sprintf('Mean    : %.2f', Mean*scaleFactor), ...
        sprintf('Median  : %.2f', Median*scaleFactor), ...
        sprintf('SD      : %.2f', SD*scaleFactor)
    };    
else
    x_text = 0;% x_text = 0;
    scaleFactor = 1;
    % Plot SNR, Mean, Median, SD below Name
    stats_text = {
        sprintf('SNR     : %.2f', SNR), ...
        sprintf('Mean    : %.4f', Mean*scaleFactor), ...
        sprintf('Median  : %.4f', Median*scaleFactor), ...
        sprintf('SD      : %.4f', SD*scaleFactor)
    };    
end
y_text = ylims(2) - 0.05 * range(ylims);

% Plot Name (large font)
text(x_text, y_text, Name, ...
     'FontSize', 35, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');


for i = 1:length(stats_text)
    text(x_text, y_text - i * 0.12 * range(ylims), stats_text{i}, ...
        'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k');
end
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.01));

hold off;
end
