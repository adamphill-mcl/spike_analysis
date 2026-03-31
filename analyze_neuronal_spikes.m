function results = analyze_neuronal_spikes(csvFile, outputDir, opts)
%ANALYZE_NEURONAL_SPIKES Separate likely neuronal spikes from cardiac spikes.
%
% results = analyze_neuronal_spikes(csvFile, outputDir, opts)
%
% Inputs
%   csvFile   : Path to Maestro/MEA spike-list CSV file.
%   outputDir : Folder for output tables, figures, and summary CSV.
%   opts      : Optional struct with fields:
%       blankWindowMs                - Total blanking window duration
%                                      centered on each cardiac seed spike.
%                                      Default: 25
%       amplitudeThresholdMode       - 'perWellPercentile', 'absolute', or
%                                      'maxOfAbsoluteAndPercentile'.
%                                      Default: 'maxOfAbsoluteAndPercentile'
%       amplitudePercentile          - Percentile used for per-well threshold.
%                                      Default: 95
%       absoluteAmplitudeThresholdMv - Fixed threshold in mV. Default: 0.50
%       minSpikesPerElectrode        - Minimum included spikes for an
%                                      electrode to count as active.
%                                      Default: 1
%       totalElectrodesPerWell       - Total electrodes per well used for
%                                      mean firing rate normalization.
%                                      Default: 16
%       makePlots                    - Save per-well plots. Default: true
%       saveTables                   - Save per-well included/excluded CSVs.
%                                      Default: true
%
% Outputs
%   results : struct with fields:
%       summaryTable
%       perWell
%       options
%
% Notes
%   - This workflow uses large-amplitude spikes as cardiac seed events.
%   - Spikes that fall inside the blanking window around those seed events
%     are excluded from the neuronal count.
%   - "Weighted firing rate" is reported as neuronal spikes per second per
%     active electrode in the well.

if nargin < 1 || isempty(csvFile)
    error('You must provide the spike-list CSV path.');
end

if nargin < 2 || isempty(outputDir)
    [csvFolder, csvName] = fileparts(csvFile);
    outputDir = fullfile(csvFolder, [csvName '_neuronal_analysis']);
end

if nargin < 3
    opts = struct();
end

opts = apply_defaults(opts);

if ~isfolder(outputDir)
    mkdir(outputDir);
end

spikeTable = read_spike_list(csvFile);
if isempty(spikeTable)
    error('No valid spike rows were found in %s.', csvFile);
end

recordingStart = min(spikeTable.Time_s);
recordingEnd = max(spikeTable.Time_s);
recordingDuration_s = recordingEnd - recordingStart;
if recordingDuration_s <= 0
    recordingDuration_s = eps;
end

wellIds = unique(spikeTable.Well, 'stable');
perWell = struct([]);
summaryRows = cell(numel(wellIds), 1);

for wellIdx = 1:numel(wellIds)
    wellId = wellIds{wellIdx};
    wellMask = strcmp(spikeTable.Well, wellId);
    wellTable = sortrows(spikeTable(wellMask, :), 'Time_s');

    [cardiacSeedMask, amplitudeThreshold_mV] = find_cardiac_seed_spikes(wellTable, opts);
    [blankedMask, blankingWindows] = build_blanking_mask(wellTable.Time_s, cardiacSeedMask, opts.blankWindowMs);

    includedTable = wellTable(~blankedMask, :);
    excludedTable = wellTable(blankedMask, :);

    activeElectrodeIds = unique(includedTable.ElectrodeID);
    countsByElectrode = groupcounts(categorical(includedTable.ElectrodeID));
    activeElectrodeCount = sum(countsByElectrode >= opts.minSpikesPerElectrode);
    totalElectrodeCount = numel(unique(wellTable.ElectrodeID));
    normalizationElectrodeCount = opts.totalElectrodesPerWell;

    neuronalSpikeCount = height(includedTable);
    excludedSpikeCount = height(excludedTable);
    meanFiringRate_Hz = neuronalSpikeCount / max(recordingDuration_s * normalizationElectrodeCount, eps);
    weightedFiringRate_Hz = neuronalSpikeCount / max(recordingDuration_s * max(activeElectrodeCount, 1), eps);

    blankedTime_s = sum(blankingWindows(:, 2) - blankingWindows(:, 1));
    adjustedRecordingTime_s = max(recordingDuration_s - blankedTime_s, eps);
    blankedFraction = blankedTime_s / recordingDuration_s;
    adjustedMeanFiringRate_Hz = neuronalSpikeCount / max(adjustedRecordingTime_s * normalizationElectrodeCount, eps);
    adjustedWeightedFiringRate_Hz = neuronalSpikeCount / max(adjustedRecordingTime_s * max(activeElectrodeCount, 1), eps);

    perWell(wellIdx).wellId = wellId;
    perWell(wellIdx).allSpikes = wellTable;
    perWell(wellIdx).includedSpikes = includedTable;
    perWell(wellIdx).excludedSpikes = excludedTable;
    perWell(wellIdx).cardiacSeedMask = cardiacSeedMask;
    perWell(wellIdx).blankedMask = blankedMask;
    perWell(wellIdx).blankingWindows_s = blankingWindows;
    perWell(wellIdx).amplitudeThreshold_mV = amplitudeThreshold_mV;
    perWell(wellIdx).recordingDuration_s = recordingDuration_s;
    perWell(wellIdx).adjustedRecordingTime_s = adjustedRecordingTime_s;
    perWell(wellIdx).meanFiringRate_Hz = meanFiringRate_Hz;
    perWell(wellIdx).weightedFiringRate_Hz = weightedFiringRate_Hz;
    perWell(wellIdx).adjustedMeanFiringRate_Hz = adjustedMeanFiringRate_Hz;
    perWell(wellIdx).adjustedWeightedFiringRate_Hz = adjustedWeightedFiringRate_Hz;
    perWell(wellIdx).activeElectrodeIds = activeElectrodeIds;
    perWell(wellIdx).activeElectrodeCount = activeElectrodeCount;
    perWell(wellIdx).totalElectrodeCount = totalElectrodeCount;
    perWell(wellIdx).normalizationElectrodeCount = normalizationElectrodeCount;

    summaryRows{wellIdx} = table( ...
        string(wellId), ...
        height(wellTable), ...
        neuronalSpikeCount, ...
        excludedSpikeCount, ...
        amplitudeThreshold_mV, ...
        totalElectrodeCount, ...
        activeElectrodeCount, ...
        normalizationElectrodeCount, ...
        recordingDuration_s, ...
        blankedTime_s, ...
        adjustedRecordingTime_s, ...
        blankedFraction, ...
        meanFiringRate_Hz, ...
        weightedFiringRate_Hz, ...
        adjustedMeanFiringRate_Hz, ...
        adjustedWeightedFiringRate_Hz, ...
        'VariableNames', { ...
        'Well', ...
        'TotalSpikes', ...
        'IncludedNeuronalSpikes', ...
        'ExcludedBlankedSpikes', ...
        'AmplitudeThreshold_mV', ...
        'TotalElectrodesWithAnySpike', ...
        'ActiveElectrodesAfterBlanking', ...
        'TotalElectrodesUsedForNormalization', ...
        'RecordingDuration_s', ...
        'BlankedTime_s', ...
        'AdjustedRecordingTime_s', ...
        'BlankedFraction', ...
        'MeanFiringRate_Hz', ...
        'WeightedFiringRate_Hz_perActiveElectrode', ...
        'AdjustedMeanFiringRate_Hz', ...
        'AdjustedWeightedFiringRate_Hz_perActiveElectrode'});

    if opts.saveTables
        writetable(includedTable, fullfile(outputDir, sprintf('%s_included_neuronal_spikes.csv', wellId)));
        writetable(excludedTable, fullfile(outputDir, sprintf('%s_excluded_blanked_spikes.csv', wellId)));
    end

    if opts.makePlots
        save_well_plot(wellTable, includedTable, excludedTable, blankingWindows, wellId, ...
            amplitudeThreshold_mV, meanFiringRate_Hz, weightedFiringRate_Hz, outputDir);
    end
end

summaryTable = vertcat(summaryRows{:});
writetable(summaryTable, fullfile(outputDir, 'well_summary.csv'));
save(fullfile(outputDir, 'neuronal_spike_analysis.mat'), 'summaryTable', 'perWell', 'opts');

results = struct();
results.summaryTable = summaryTable;
results.perWell = perWell;
results.options = opts;
results.outputDir = outputDir;

end

function opts = apply_defaults(opts)
defaults.blankWindowMs = 25;
defaults.amplitudeThresholdMode = 'maxOfAbsoluteAndPercentile';
defaults.amplitudePercentile = 95;
defaults.absoluteAmplitudeThresholdMv = 0.50;
defaults.minSpikesPerElectrode = 1;
defaults.totalElectrodesPerWell = 16;
defaults.makePlots = true;
defaults.saveTables = true;

defaultFields = fieldnames(defaults);
for i = 1:numel(defaultFields)
    fieldName = defaultFields{i};
    if ~isfield(opts, fieldName) || isempty(opts.(fieldName))
        opts.(fieldName) = defaults.(fieldName);
    end
end
end

function spikeTable = read_spike_list(csvFile)
raw = readcell(csvFile, 'FileType', 'text', 'TextType', 'string', 'DatetimeType', 'text');

if size(raw, 2) < 5
    error('Expected at least 5 columns in the spike-list CSV.');
end

timeCol = raw(:, 3);
electrodeCol = raw(:, 4);
amplitudeCol = raw(:, 5);

timeVals = nan(size(timeCol));
ampVals = nan(size(amplitudeCol));

for i = 1:numel(timeCol)
    timeVals(i) = to_numeric_scalar(timeCol{i});
    ampVals(i) = to_numeric_scalar(amplitudeCol{i});
end

validElectrode = cellfun(@(x) ischar(x) || isstring(x), electrodeCol);
electrodeText = strings(size(electrodeCol));
electrodeText(validElectrode) = string(electrodeCol(validElectrode));
electrodeText = strtrim(electrodeText);

isSpikeRow = isfinite(timeVals) & isfinite(ampVals) & contains(electrodeText, "_");
spikeTable = table();
spikeTable.Time_s = timeVals(isSpikeRow);
spikeTable.ElectrodeID = cellstr(electrodeText(isSpikeRow));
spikeTable.Amplitude_mV = ampVals(isSpikeRow);
spikeTable.AbsAmplitude_mV = abs(spikeTable.Amplitude_mV);
spikeTable.Well = regexprep(spikeTable.ElectrodeID, '_.*$', '');

spikeTable = sortrows(spikeTable, {'Well', 'Time_s', 'ElectrodeID'});
end

function value = to_numeric_scalar(entry)
value = NaN;
if isnumeric(entry) && isscalar(entry)
    value = double(entry);
elseif ischar(entry) || isstring(entry)
    value = str2double(entry);
end
end

function [cardiacSeedMask, threshold_mV] = find_cardiac_seed_spikes(wellTable, opts)
absAmp = wellTable.AbsAmplitude_mV;
percentileThreshold = prctile(absAmp, opts.amplitudePercentile);

switch lower(opts.amplitudeThresholdMode)
    case 'perwellpercentile'
        threshold_mV = percentileThreshold;
    case 'absolute'
        threshold_mV = opts.absoluteAmplitudeThresholdMv;
    case 'maxofabsoluteandpercentile'
        threshold_mV = max(opts.absoluteAmplitudeThresholdMv, percentileThreshold);
    otherwise
        error('Unknown amplitudeThresholdMode: %s', opts.amplitudeThresholdMode);
end

cardiacSeedMask = absAmp >= threshold_mV;
end

function [blankedMask, windows] = build_blanking_mask(times_s, cardiacSeedMask, blankWindowMs)
halfWindow_s = (blankWindowMs / 1000) / 2;
seedTimes = times_s(cardiacSeedMask);
blankedMask = false(size(times_s));

if isempty(seedTimes)
    windows = zeros(0, 2);
    return;
end

windows = [seedTimes - halfWindow_s, seedTimes + halfWindow_s];
windows = sortrows(windows, 1);
windows = merge_intervals(windows);

for i = 1:size(windows, 1)
    blankedMask = blankedMask | (times_s >= windows(i, 1) & times_s <= windows(i, 2));
end
end

function merged = merge_intervals(intervals)
if isempty(intervals)
    merged = intervals;
    return;
end

merged = intervals(1, :);
for i = 2:size(intervals, 1)
    current = intervals(i, :);
    if current(1) <= merged(end, 2)
        merged(end, 2) = max(merged(end, 2), current(2));
    else
        merged = [merged; current]; %#ok<AGROW>
    end
end
end

function save_well_plot(allSpikes, includedSpikes, excludedSpikes, blankingWindows, ...
    wellId, amplitudeThreshold_mV, meanFiringRate_Hz, weightedFiringRate_Hz, outputDir)

electrodes = unique(allSpikes.ElectrodeID, 'stable');
electrodeMap = containers.Map(electrodes, 1:numel(electrodes));

includedY = map_electrodes_to_index(includedSpikes.ElectrodeID, electrodeMap);
excludedY = map_electrodes_to_index(excludedSpikes.ElectrodeID, electrodeMap);

f = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1400, 800]);
t = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t, 1);
hold(ax1, 'on');
scatter(ax1, excludedSpikes.Time_s, excludedY, 10, [0.85 0.33 0.10], 'filled', ...
    'MarkerFaceAlpha', 0.65, 'DisplayName', 'Excluded / blanked');
scatter(ax1, includedSpikes.Time_s, includedY, 10, [0.00 0.45 0.74], 'filled', ...
    'MarkerFaceAlpha', 0.65, 'DisplayName', 'Included neuronal');
ylim(ax1, [0, max(numel(electrodes) + 1, 2)]);
shade_windows(ax1, blankingWindows);
hold(ax1, 'off');
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Electrode');
yticks(ax1, 1:numel(electrodes));
yticklabels(ax1, electrodes);
title(ax1, sprintf('%s spike raster with blanked windows', wellId));
legend(ax1, 'Location', 'eastoutside');
grid(ax1, 'on');

ax2 = nexttile(t, 2);
hold(ax2, 'on');
scatter(ax2, excludedSpikes.Time_s, excludedSpikes.AbsAmplitude_mV, 10, [0.85 0.33 0.10], 'filled', ...
    'MarkerFaceAlpha', 0.65, 'DisplayName', 'Excluded / blanked');
scatter(ax2, includedSpikes.Time_s, includedSpikes.AbsAmplitude_mV, 10, [0.00 0.45 0.74], 'filled', ...
    'MarkerFaceAlpha', 0.65, 'DisplayName', 'Included neuronal');
yline(ax2, amplitudeThreshold_mV, '--k', sprintf('Threshold = %.3f mV', amplitudeThreshold_mV), ...
    'LabelHorizontalAlignment', 'left', 'DisplayName', 'Cardiac seed threshold');
shade_windows(ax2, blankingWindows);
hold(ax2, 'off');
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Absolute amplitude (mV)');
title(ax2, sprintf('Mean FR = %.3f Hz | Weighted FR = %.3f Hz/electrode', ...
    meanFiringRate_Hz, weightedFiringRate_Hz));
grid(ax2, 'on');

exportgraphics(f, fullfile(outputDir, sprintf('%s_neuronal_filtering.png', wellId)), 'Resolution', 200);
set(f, 'Visible', 'on');
savefig(f, fullfile(outputDir, sprintf('%s_neuronal_filtering.fig', wellId)));
close(f);
end

function y = map_electrodes_to_index(electrodeIds, electrodeMap)
y = zeros(numel(electrodeIds), 1);
for i = 1:numel(electrodeIds)
    y(i) = electrodeMap(electrodeIds{i});
end
end

function shade_windows(ax, windows)
if isempty(windows)
    return;
end

yl = ylim(ax);
for i = 1:size(windows, 1)
    patch(ax, ...
        [windows(i, 1), windows(i, 2), windows(i, 2), windows(i, 1)], ...
        [yl(1), yl(1), yl(2), yl(2)], ...
        [0.85, 0.85, 0.85], ...
        'FaceAlpha', 0.35, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
uistack(findobj(ax, 'Type', 'Scatter'), 'top');
end
