function [td,s,frequencyBuffer,isH] = stimGen_static_HL(loFreq,hiFreq,toneDur,toneSOA,trialDur,coh,fs);
% stimGen_static_HL generates a sequence of tone bursts of 2 frequencies (loFreq & hiFreq)
% each tone lasts 'toneDur' ms and is followed by an inter-tone interval of 'toneSOA'
% the ratio of loFreq and hiFreq tones is determined by 'coh'
% the length of the audio stream is 'trialDur' ms 

% INPUT:
% loFreq    - lowest frequency in Hz
% hiFreq    - highest frequency in Hz
% toneDur   - tone duration in ms
% toneSOA   - stimulus onset asynchrony = inter-tone interval in ms
% trialDur  - trial duration in ms
% coh       - coherence 0, 0.25, 0.5, 0.75, 1

% OUTPUT:
% td - time bins
% s - auditory stream

%% some default variable for function testing
% loFreq = 312.5; %hz  625 | 1250 | 2500 | 5000
% hiFreq = 625; %hz   1250 | 2500 | 5000 | 10000
% toneDur = 40; %ms
% toneSOA = 10; %ms
% trialDur = 2000;
% coh = 0.7;

%% generate sequence of frequencies
% calc no. of tone presentations
bufferLen = ceil(trialDur/(toneSOA+toneDur));
numTonePresentations = floor(trialDur/(toneSOA+toneDur));
if coh == 0.5
    numToShuffle = round((1-coh)*numTonePresentations + ((rand(1)-0.5)*0.1));  % add some fluctuations to get low & high trials in 0-coherence trials
else
    numToShuffle = round((1-coh)*numTonePresentations);
end
frequencyBuffer = ones(bufferLen,1)*hiFreq;

%stateList generated using randperm instead of flipping a coin
% per channel, per tone presentation period

stateList = zeros(numTonePresentations,1);

stateIdx = randperm(numTonePresentations);
stateIdx = stateIdx(1,1:numToShuffle);

stateList(stateIdx) = 1;

frequencyBuffer(stateIdx) = loFreq; %low

numLo = sum(frequencyBuffer == loFreq);
numHi = sum(frequencyBuffer == hiFreq);

isH = numHi > numLo;

%% generate sequence of tones
% fs = 348000;
rampLength = 5;
cosRamp = (rampLength/1000)*fs;

multiplier = floor(fs/1000);
s = zeros(1,trialDur*multiplier);
ind = 1;
for ii = 1:numTonePresentations
    ind = (ii-1)*(toneSOA+toneDur)*multiplier + 1;
    tmax = toneDur/1000;                               % Time Duration Of Signal (sec)
    t = linspace(0, tmax, tmax*fs);                          % Time Vector
    f = frequencyBuffer(ii);                                 % Original Frequency
    tone = sin(2*pi*f*t);% Original Signal %TODO: multiply voltage - from calib
    tone = pa_ramp(tone, cosRamp, fs); 
    s(ind:ind+length(tone)-1) = tone';
    if ii < numTonePresentations
        ind = ind+length(tone)+(toneSOA/1000*fs);
    else
        ind = ind+length(tone)-1;
        s = s(1:ind);
    end
end

td = 0:1/fs:(length(s)-1)/fs;

%% for testing - listen
%  sound(s, fs) 