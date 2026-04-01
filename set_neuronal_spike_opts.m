function opts = set_neuronal_spike_opts()
%DEFAULT_NEURONAL_SPIKE_OPTS Return default options for neuronal spike analysis.
%
% Example:
%   opts = default_neuronal_spike_opts();
%   results = analyze_neuronal_spikes(csvFile, outputDir, opts);
% batchResults = run_neuronal_spike_analysis_batch([], [], set_neuronal_spike_opts)

opts = struct();

opts.blankWindowMs = 25;
opts.amplitudeThresholdMode = 'perWellPercentile'; %'maxOfAbsoluteAndPercentile';
opts.amplitudePercentile = 90; %changed from 95
opts.absoluteAmplitudeThresholdMv = 0.50;
opts.minSpikesPerElectrode = 1;
opts.totalElectrodesPerWell = 16;
opts.makePlots = true;
opts.saveTables = true;

end
