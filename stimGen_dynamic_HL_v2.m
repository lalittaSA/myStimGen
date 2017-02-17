function [td,s,frequencyBuffer,isH] = stimGen_dynamic_HL_v2(loFreq,hiFreq,toneDur,toneSOA,trialDur,freqType,cohLevel,changeTime,changeRamp,fs)
% stimGen_static_HL_v2 generates a sequence of tone bursts of 2 frequencies (loFreq & hiFreq)
% each tone lasts 'toneDur' ms and is followed by an inter-tone interval of 'toneSOA'
% the length of the audio stream is 'trialDur' ms 
% the patterns of loFreq and hiFreq tones is determined by 'freqList'
% determined by list of letters e.g. LLL, HLL, HHL, LNL, HSH and all other combinations
% H for High | L for Low | N for Noise | S for Silence
% the coherence level can vary between 0-1: 
% 0 means 100% low-frequency tones
% 1 means 100% high-frequency tones
% if length(cohLevel) == length(freqList), coherence level is set for each frequency fraction
% the time points on which one frequency changes to the next is set using 'changeTime' 
% length(changeTime) has to match length(freqList)-1
% frequency change can be abrupt or ramping 
% 'changRamp' has to be larger than toneDur+toneSOA to set ramping slope,
% otherwise the change will be abrupt


% INPUT:
% loFreq    - lowest frequency in Hz (> 300 Hz)
% hiFreq    - highest frequency in Hz (< 10000 Hz and 2 octaves above loFreq)
% toneDur   - tone duration in ms
% toneSOA   - stimulus onset asynchrony = inter-tone interval in ms
% trialDur  - trial duration in ms
% freqList  - frequency patterns [L | LNL | LSL | HNH | HSH | H]
% cohLevel  - coherence level for high-frequency tone [0.5 - 1] (1-cohLevel for low-frequency tone)
% changeTime - time points when one frequency changes from one to the next
% changeRamp - time span (ms) of change slope (< toneDur+toneSOA - abrupt change)

% OUTPUT:
% td - time bins
% s - auditory stream

%% some default variable for function testing
% loFreq = 625; %hz      312.5 |  625 | 1250 | 2500 |  5000
% hiFreq = 1250; %hz     625   | 1250 | 2500 | 5000 | 10000
% toneDur = 40; %ms
% toneSOA = 10; %ms
% trialDur = 2000; %ms
% freqType = 'LLL';
% cohLevel = [1 0.7 1];
% changeTime = [300 1000]; %ms
% changeRamp = 500; %ms

%% check inputs
if any(changeTime>trialDur), disp('invalid changeTime input: later than trial length'); return; end
if length(changeTime) ~= length(freqType)-1, disp('invalid changeTime or freqType input: lengths do not match'); return; end
if any(diff(changeTime)<changeRamp), disp('invalid changeTime interval: at least one interval is smaller than ramping time'); return; end

%% generate sequence of coherence values
% assign high & low & noise coherence levels
if length(cohLevel) == 1
    hiCoh = cohLevel;
    loCoh = 1-cohLevel;
end
noCoh = 0.5;
% add first time bin to changeTime
if changeTime(1) == 0
    changeTime = changeTime(2:end);
    freqType = freqType(2:end);
end
changeTime = [1, changeTime];

%% gradual change
toneDurTot = toneSOA+toneDur;
numTonePresentations = floor(trialDur/toneDurTot);

% convert times (ms) to coherence bins - depending on tone length
changeBins = unique(ceil((changeTime+toneSOA+toneDur/2)/toneDurTot));

coh = zeros(1,numTonePresentations);
ind_begin = 1;
for ff = 1:length(freqType)  
    if length(cohLevel) > 1
        hiCoh = cohLevel(ff);
        loCoh = 1 - cohLevel(ff);
    end
    switch freqType(ff)
        case 'H', tmpCoh = hiCoh;       % high frequency
        case 'L', tmpCoh = loCoh;       % low frequency
        case 'N', tmpCoh = noCoh;       % noisy signal
        case 'S', tmpCoh = nan;         % absence of signal (silence)
    end

    if ff < length(freqType)
        ind_end = changeBins(ff+1);
        coh(ind_begin:ind_end-1) = tmpCoh * ones(1,ind_end-ind_begin);
        ind_begin = ind_end;
    else
        ind_end = numTonePresentations;
        coh(ind_begin:ind_end) = tmpCoh * ones(1,ind_end-ind_begin+1);
    end
    
    
    % in case of ramping change
    if changeRamp > toneDurTot && ff < length(freqType)
        nSteps = floor(changeRamp/toneDurTot)+1;
        if length(cohLevel) > 1
            next_hiCoh = cohLevel(ff+1);
            next_loCoh = 1 - cohLevel(ff+1);
        else
            next_hiCoh = hiCoh;
            next_loCoh = loCoh;
        end
        switch freqType(ff+1)
            case 'H', rampCoh = tmpCoh:(next_hiCoh-tmpCoh)/nSteps:next_hiCoh;
            case 'L', rampCoh = tmpCoh:(next_loCoh-tmpCoh)/nSteps:next_loCoh;
            case 'N', rampCoh = tmpCoh:(noCoh-tmpCoh)/nSteps:noCoh;
            case 'S', rampCoh = repmat(tmpCoh,1,nSteps+1); % no ramping in case of silent pause
        end
        % take only the transitional steps - remove first (current) & last (next) coherence levels
        rampCoh = rampCoh(2:end-1);
        ind_end = ind_begin + length(rampCoh) - 1;
        coh(ind_begin:ind_end) = rampCoh;
        ind_begin = ind_end + 1;
    end
end

% calc no. of tone presentations
bufferLen = ceil(trialDur/toneDurTot);
frequencyBuffer = ones(1,bufferLen)*hiFreq;

%% change high freq to low freq depending on coherence level by permuting tones in a block of the same coherence
% stateList = zeros(numTonePresentations,numCh);
% for i = 1:numTonePresentations
%     stateList(i,:) = randperm(numCh);
% end
% cohChanges = [1 find(diff(coh))+1 numTonePresentations+1];
% for i = 1:numCh
%     ind = 0;
%     for cc = 2:length(cohChanges)
%         ind_coh = cohChanges(cc-1);
%         numToneinBloc = cohChanges(cc) - ind_coh;
%         if ~isnan(coh(ind_coh))
%             numToShuffle = round((1-coh(ind_coh))*numToneinBloc);
%             stateIdx = randperm(numToneinBloc);
%             stateIdx = stateIdx(1,1:numToShuffle);
% 
% %             stateList(ind+stateIdx,i) = 1;
%             if ~isempty(stateIdx)
%                 frequencyBuffer(ind+stateIdx,i) = loFreq; %low
%             end
%         else
%             frequencyBuffer(ind+1:ind+numToneinBloc,i) = nan;
%         end
%         ind = ind + numToneinBloc;
%     end
%     numLo = sum(frequencyBuffer(:,i) == loFreq);
%     numHi = sum(frequencyBuffer(:,i) == hiFreq);
% end

%% generate random numbers and compare against coherence values
tmp_rand = rand(1,numTonePresentations);
ind_change = tmp_rand > coh;
frequencyBuffer(ind_change) = loFreq;
ind_nan = isnan(coh);
frequencyBuffer(ind_nan) = nan;

numLo = sum(frequencyBuffer == loFreq);
numHi = sum(frequencyBuffer == hiFreq);

isH = numHi > numLo;

%% generate sequence of tones
% fs = 10000;%384000;%12000;
rampLength = 5;
cosRamp = (rampLength/1000)*fs;

multiplier = floor(fs/1000);
s = zeros(1,trialDur*multiplier);
f = nan(1,trialDur*multiplier);
ind = 1;
for ii = 1:numTonePresentations
    if ~isnan(frequencyBuffer(ii))
        tmax = toneDur/1000;                               % Time Duration Of Signal (sec)
        t = linspace(0, tmax, tmax*fs);                          % Time Vector
        freq = frequencyBuffer(ii);                                 % Original Frequency
        tone = sin(2*pi*freq*t);% Original Signal %TODO: multiply voltage - from calib
        tone = pa_ramp(tone, cosRamp, fs);
        s(ind:ind+length(tone)-1) = tone';
        f(ind:ind+length(tone)-1) = freq;
    end
    ind = ind+length(tone);
end

td = 0:1/fs:(length(s)-1)/fs;

%% for testing - listen
% sound(s, fs) 

% record
% filename = 'sample_stim2.wav';
% audiowrite(filename,s,fs)
%% for testing - plot
% temp_fig_path = '/dataAnalysis/git_public/Penn_auditoryDecision/stimuli/';
% fig_name = 'sample_stim2_toneplot';
% h = figure('Name',fig_name,'Position',get(0,'ScreenSize'));
% subplot(2,1,1)
% tone = sin(2*pi*freq*t);
% plot(tone)
% tone = pa_ramp(tone, cosRamp, fs);
% subplot(2,1,2)
% plot(tone)
% saveas(h,[temp_fig_path fig_name '.fig']);
% save2pdf([temp_fig_path fig_name '.pdf'],h);
% clf
% 
% % 
% % plot(td,s)
% % xlabel('Time (ms)')
% % ylabel('Amplitude')
% % title('Signal')
% 
% temp_fig_path = '/dataAnalysis/git_public/Penn_auditoryDecision/stimuli/';
% fig_name = 'sample_stim2_plot';
% h = figure('Name',fig_name,'Position',get(0,'ScreenSize'));
% subplot(2,1,1)
% plot(td,f/1000)
% text(td(end),(hiFreq+100)/1000,['numHi = ' num2str(numHi)],'HorizontalAlignment','right');
% text(td(end),(loFreq-100)/1000,['numLo = ' num2str(numLo)],'HorizontalAlignment','right');
% xlabel('Time (secs)')
% ylabel('Frequency(Hz)')
% ylim([0 3])
% % 
% subplot(2,1,2)
% segmentLength = round(numel(s)/5);
% window = hamming(512);
% noverlap = 256;
% nfft = 512;
% [S,F,T,P] = spectrogram(s,window,noverlap,nfft,fs,'yaxis');
% surf(T,F,10*log10(P),'edgecolor','none'); axis tight; view(0,90);
% colormap(hot)
% set(gca,'clim',[-80 -30])
% ylim([0 3])
% 
% periodogram(s,[],[],fs)
% 
% % spectrogram(s,segmentLength,[],[],fs,'yaxis')
% spectrogram(s,round(segmentLength/5),round(80/100*segmentLength/5),[],fs,'yaxis')
% % spectrogram(s,round(numel(s)/25),[],round(numel(s)/25),fs,'yaxis')
% ylim([0 3])
% saveas(h,[temp_fig_path fig_name '.fig']);
% save2pdf([temp_fig_path fig_name '.pdf'],h);
% caf