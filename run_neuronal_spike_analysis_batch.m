function batchResults = run_neuronal_spike_analysis_batch(inputDir, outputRootDir, opts)
%RUN_NEURONAL_SPIKE_ANALYSIS_BATCH Run neuronal spike analysis for each CSV in a folder.
%
% batchResults = run_neuronal_spike_analysis_batch(inputDir, outputRootDir, opts)
%
% Example:
%   inputDir = '/Users/adamphill/Documents/Matlab/Spike_analysis';
%   outputRootDir = fullfile(inputDir, 'batch_neuronal_output');
%   opts = struct();
%   batchResults = run_neuronal_spike_analysis_batch(inputDir, outputRootDir, opts);

if nargin < 1 || isempty(inputDir)
    inputDir = uigetdir(pwd, 'Select folder containing spike-list CSV files');
    if isequal(inputDir, 0)
        error('No input folder selected.');
    end
end

if nargin < 2 || isempty(outputRootDir)
    outputRootDir = fullfile(inputDir, 'batch_neuronal_output');
end

if nargin < 3
    opts = struct();
end

if ~isfolder(inputDir)
    error('Input directory does not exist: %s', inputDir);
end

if ~isfolder(outputRootDir)
    mkdir(outputRootDir);
end

csvFiles = dir(fullfile(inputDir, '*.csv'));
if isempty(csvFiles)
    error('No CSV files were found in %s.', inputDir);
end

batchSummary = table();
batchResults = struct('fileName', {}, 'filePath', {}, 'outputDir', {}, 'status', {}, 'message', {}, 'summaryTable', {});

for fileIdx = 1:numel(csvFiles)
    csvName = csvFiles(fileIdx).name;
    csvPath = fullfile(csvFiles(fileIdx).folder, csvName);
    [~, baseName] = fileparts(csvName);
    safeBaseName = matlab.lang.makeValidName(baseName);
    fileOutputDir = fullfile(outputRootDir, safeBaseName);

    fprintf('Processing %d/%d: %s\n', fileIdx, numel(csvFiles), csvName);

    try
        results = analyze_neuronal_spikes(csvPath, fileOutputDir, opts);

        fileSummary = results.summaryTable;
        fileSummary.SourceFile = repmat(string(csvName), height(fileSummary), 1);
        fileSummary = movevars(fileSummary, 'SourceFile', 'Before', 1);
        batchSummary = [batchSummary; fileSummary]; %#ok<AGROW>

        batchResults(end+1).fileName = csvName; %#ok<AGROW>
        batchResults(end).filePath = csvPath;
        batchResults(end).outputDir = fileOutputDir;
        batchResults(end).status = "success";
        batchResults(end).message = "";
        batchResults(end).summaryTable = results.summaryTable;
    catch ME
        warning('Failed to process %s: %s', csvName, ME.message);

        batchResults(end+1).fileName = csvName; %#ok<AGROW>
        batchResults(end).filePath = csvPath;
        batchResults(end).outputDir = fileOutputDir;
        batchResults(end).status = "failed";
        batchResults(end).message = string(ME.message);
        batchResults(end).summaryTable = table();
    end
end

if ~isempty(batchSummary)
    writetable(batchSummary, fullfile(outputRootDir, 'batch_well_summary.csv'));
end

statusTable = struct2table(rmfield(batchResults, 'summaryTable'));
writetable(statusTable, fullfile(outputRootDir, 'batch_run_status.csv'));
save(fullfile(outputRootDir, 'batch_neuronal_spike_analysis.mat'), 'batchResults', 'batchSummary', 'opts', 'inputDir', 'outputRootDir');

end
