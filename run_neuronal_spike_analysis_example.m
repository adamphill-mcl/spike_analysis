csvFile = '/Users/adamphill/Documents/Matlab/Spike_analysis/2. D7_MSCT_Ratioexp_Vehicle(001)(000)_spike_list.csv';
outputDir = '/Users/adamphill/Documents/Matlab/Spike_analysis/neuronal_spike_output';

opts = struct();
opts.blankWindowMs = 25;
opts.amplitudeThresholdMode = 'maxOfAbsoluteAndPercentile';
opts.amplitudePercentile = 95;
opts.absoluteAmplitudeThresholdMv = 0.50;
opts.minSpikesPerElectrode = 1;
opts.makePlots = true;
opts.saveTables = true;

results = analyze_neuronal_spikes(csvFile, outputDir, opts);
disp(results.summaryTable);
