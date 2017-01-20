%% tuning task
loFreq = 300; %hz
hiFreq = 12000; %hz

fs = 48828.125;

toneDur = 100; %ms
toneSOA = 400; %ms
nTicksPerBurst = ceil((toneDur+toneSOA)/1000* fs);

freqData = loFreq*10.^(log10(2) * (0:0.33:ceil(log10(hiFreq/loFreq)/log10(2))));

permuteMatrix = randperm(length(freqData));

trialDur = length(permuteMatrix)*(toneDur+toneSOA)+1000; % in ms
%% static coherence
loFreq = 200;
hiFreq = 1200;
coh = 0.5;
trialDur = 3000;
stimFlag = false;

rampLength = 5;

toneSOA = 10;
toneDur = 40;

numCh = 1;

% calc no. of tone presentations
bufferLen = ceil(trialDur/(toneSOA+toneDur));
numTonePresentations = floor(trialDur/(toneSOA+toneDur));
numToShuffle = floor((1-coh)*numTonePresentations);
frequencyBuffer = ones(bufferLen,numCh)*hiFreq;

%stateList generated using randperm instead of flipping a coin
% per channel, per tone presentation period

stateList = zeros(numTonePresentations,numCh);
for i = 1:numTonePresentations
    stateList(i,:) = randperm(numCh);
end

for i = 1:numCh
    stateIdx = randperm(numTonePresentations);
    stateIdx = stateIdx(1,1:numToShuffle);
    
    stateList(stateIdx,i) = 1;
    
    frequencyBuffer(stateIdx,i) = loFreq; %low
end

%% dynamic coherence
coh = 0:0.1:1;
coh = [coh,fliplr(coh)];

numInCohBloc = 6;

toneSOA = 10;
toneDur = 40;

loFreq = 500;
hiFreq = 1500;
trialDur = length(coh) * numInCohBloc * (toneSOA+toneDur);
stimFlag = false;

numCh = 1;

% calc no. of tone presentations
bufferLen = ceil(trialDur/(toneSOA+toneDur));
numTonePresentations = floor(trialDur/(toneSOA+toneDur));
frequencyBuffer = ones(bufferLen,numCh)*hiFreq;

stateList = zeros(numTonePresentations,numCh);
for i = 1:numTonePresentations
    stateList(i,:) = randperm(numCh);
end

% numInCohBloc = ceil(numTonePresentations/length(coh));
for i = 1:numCh
    for bb = 1:length(coh)
        numToShuffle = round((1-coh(bb))*numInCohBloc);
        stateIdx = randperm(numInCohBloc);
        stateIdx = stateIdx(1,1:numToShuffle);
        
        ind = (bb-1)*numInCohBloc+1;
        
        stateList(ind+stateIdx,i) = 1;
        
        frequencyBuffer(ind+stateIdx,i) = loFreq; %low
    end
end

%% play sound

Fs = 348000;
cosRamp = (rampLength/1000)*Fs;

multiplier = Fs/1000;
s = zeros(1,trialDur*multiplier);
for ii = 1:numTonePresentations
    ind = (ii-1)*(toneSOA+toneDur)*multiplier + 1;
    tmax = toneDur/1000;                               % Time Duration Of Signal (sec)
    t = linspace(0, tmax, tmax*Fs);                          % Time Vector
    f = frequencyBuffer(ii);                                 % Original Frequency
    s(ind:ind+toneDur*multiplier-1) = sin(2*pi*f*t);% Original Signal
end
sound(s, Fs)                                % Listen To Original Signal
% 
% Fn1 = Fs/2;                                 % Compute & Plot Fourier Series
% Ft1 = fft(s)/length(s);
% Fv1 = linspace(0, 1, fix(length(Ft1)/2)+1)*Fn1;
% Ix1 = 1:length(Fv1);
