function [td,s] = stimGen_tuning(loFreq,hiFreq,toneDur,toneSOA,numRep)
% stimGen_tuning generates tone bursts of pseudo random frequencies between freq range [loFreq - hiFreq] with log step
% each tone lasts 'toneDur' ms and is followed by an inter-tone interval of 'toneSOA'
% each tone will be played 'numRep' times

% INPUT:
% loFreq  - lowest frequency in Hz
% hiFreq  - highest frequency in Hz
% toneDur - tone duration in ms
% toneSOA - stimulus onset asynchrony = inter-tone interval in ms
% numRep  - number of repetitions

% OUTPUT:
% td - time bins
% s - auditory stream

%% some default variable for function testing
% loFreq = 300; %hz
% hiFreq = 12000; %hz
% toneDur = 100; %ms
% toneSOA = 400; %ms
% numRep = 10;

%% generate sequence of frequencies
% fs = 48828.125; % from earlier version
fs = 348000; % max in MATLAB
nTicksPerBurst = ceil((toneDur+toneSOA)/1000* fs);

freqData = loFreq*10.^(log10(2) * (0:0.33:ceil(log10(hiFreq/loFreq)/log10(2))));
numTonePresentations = length(freqData);

trialDur = numRep * numTonePresentations * (toneDur+toneSOA) + 1000; % in ms

frequencyBuffer = [];
for nn = 1:numRep
    permuteMatrix = randperm(numTonePresentations);
    frequencyBuffer = [frequencyBuffer freqData(permuteMatrix)];
end

numTonePresentations = numTonePresentations * numRep;

%% generate sequence of tones
rampLength = 5;
cosRamp = (rampLength/1000)*fs;

multiplier = fs/1000;
s = zeros(1,trialDur*multiplier);
for ii = 1:numTonePresentations
    ind = (ii-1)*(toneSOA+toneDur)*multiplier + 1;
    tmax = toneDur/1000;                               % Time Duration Of Signal (sec)
    t = linspace(0, tmax, tmax*fs);                          % Time Vector
    f = frequencyBuffer(ii);                                 % Original Frequency
    tone = sin(2*pi*f*t);% Original Signal %TODO: multiply voltage - from calib
    tone = pa_ramp(tone, cosRamp, fs); 
    s(ind:ind+toneDur*multiplier-1) = tone';
end
td = 0:1/fs:(length(s)-1)/fs;

%% for testing - listen
% sound(s, fs) 