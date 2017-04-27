function [td,s,frequencyBuffer,isH] = stimGen_dynamic_IBI_HL(loFreq,hiFreq,toneDur,IBI_range,trialDur,coh,fs);
% stimGen_static_HL generates a sequence of tone bursts of 2 frequencies (loFreq & hiFreq)
% each tone lasts 'toneDur' ms and is followed by an inter-tone interval of 'toneSOA'
% the ratio of loFreq and hiFreq tones is determined by 'coh'
% the length of the audio stream is 'trialDur' ms 

% INPUT:
% loFreq    - lowest frequency in Hz
% hiFreq    - highest frequency in Hz
% toneDur   - tone duration in ms
% IBI_range - inter-burst interval range in ms
% trialDur  - trial duration in ms
% coh       - coherence 0, 0.25, 0.5, 0.75, 1

% OUTPUT:
% td - time bins
% s - auditory stream

%% some default variable for function testing
loFreq = 500; %hz  625 | 1250 | 2500 | 5000
hiFreq = 2000; %hz   1250 | 2500 | 5000 | 10000
toneDur = 20; %ms
IBI_minmax = [30,150]; %ms
trialDur = 2000;
coh = 1; % 0 - 1
isInc = 1;

if ~isInc
    IBI_minmax = fliplr(IBI_minmax);
end

%% generate sequence of frequencies
% calc no. of tone presentations
toneSOA = mean(IBI_minmax);
bufferLen = ceil(trialDur/(toneSOA+toneDur));
numTonePresentations = floor(trialDur/(toneSOA+toneDur));
IBI_steps = (IBI_minmax(2)-IBI_minmax(1))/numTonePresentations;
IBI_range = IBI_minmax(1):IBI_steps:IBI_minmax(2);
% if coh == 0.5
%     numToShuffle = round((1-coh)*numTonePresentations + ((rand(1)-0.5)*0.1));  % add some fluctuations to get low & high trials in 0-coherence trials
% else
%     numToShuffle = round((1-coh)*numTonePresentations);
% end
frequencyBuffer = ones(bufferLen,1)*hiFreq;
numToShuffle = round((1-coh)*bufferLen);


stateIdx = randperm(bufferLen);
stateIdx = stateIdx(1:numToShuffle);
IBI_range(sort(stateIdx)) = IBI_range(stateIdx);


%% generate sequence of tones
fs = 10000;
rampLength = 5;
cosRamp = (rampLength/1000)*fs;

multiplier = floor(fs/1000);
s = zeros(1,trialDur*multiplier);
f = nan(1,trialDur*multiplier);
ind = 1;
for ii = 1:bufferLen
    tmax = toneDur/1000;                               % Time Duration Of Signal (sec)
    t = linspace(0, tmax, tmax*fs);                          % Time Vector
    freq = frequencyBuffer(ii);                                 % Original Frequency
    tone = sin(2*pi*freq*t);% Original Signal %TODO: multiply voltage - from calib
    tone = pa_ramp(tone, cosRamp, fs); 
    s(ind:ind+length(tone)-1) = tone';
    f(ind:ind+length(tone)-1) = freq;
    ind = ind+length(tone)+round(IBI_range(ii)/1000*fs);
end

td = 0:1/fs:(length(s)-1)/fs;

%% for testing - listen & plot
% sound(s, fs) 
%  
% temp_fig_path = '/dataAnalysis/git_public/Penn_auditoryDecision/stimuli/';
% fig_name = 'sample_stim2_plot';
% h = figure('Name',fig_name,'Position',get(0,'ScreenSize'));
% subplot(2,1,1)
% plot(td,f/1000)
% % text(td(end),(hiFreq+100)/1000,['numHi = ' num2str(numHi)],'HorizontalAlignment','right');
% % text(td(end),(loFreq-100)/1000,['numLo = ' num2str(numLo)],'HorizontalAlignment','right');
% xlabel('Time (secs)')
% ylabel('Frequency(Hz)')
% ylim([0 3])