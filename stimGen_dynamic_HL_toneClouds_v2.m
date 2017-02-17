function [td,s] = stimGen_dynamic_HL_toneClouds_v2(loFreq,hiFreq,toneDur,toneSOA,trialDur,freqList,freqWidth,nFreq,changeTime,changeRamp)
% stimGen_static_toneClouds_v2 generates a sequence of tone bursts of frequency distributions
% Each tone lasts 'toneDur' ms and is followed by an inter-tone interval of 'toneSOA'.
% Length of the audio stream is 'trialDur' ms. 
% Patterns of loFreq and hiFreq tones is determined by 'freqList'
% determined by list of letters e.g. LLL, HLL, HHL, LNL, HSH and all other combinations
% H for High | L for Low | N for Noise | S for Silence.
% Width of frequency distribution is determined by 'freqWidth' (unit = n octaves)
% and the number of frequencies generated by the distribution is 'nFreq'.
% In this version the coherence level is fixed but we vary frequency levels.
% Time points on which one frequency changes to the next is set using 'changeTime' 
% length(changeTime) has to match length(freqList)-1
% Frequency change can be abrupt or ramping 
% 'changRamp' has to be larger than toneDur+toneSOA to set ramping slope,
% otherwise the change will be abrupt.



% INPUT:
% loFreq    - lowest frequency in Hz (> 300 Hz)
% hiFreq    - highest frequency in Hz (< 10000 Hz and 2 octaves above loFreq)
% toneDur   - tone duration in ms
% toneSOA   - stimulus onset asynchrony = inter-tone interval in ms
% trialDur  - trial duration in ms
% freqList  - frequency patterns [L | LNL | LSL | HNH | HSH | H]
% freqWidth - width of frequency distribution, used to generate the 'tone clouds'
% changeTime - time points when one frequency changes from one to the next
% changeRamp - time span (ms) of change slope (< toneDur+toneSOA - abrupt change)

% OUTPUT:
% td - time bins
% s - auditory stream

%% some default variable for function testing
loFreq = 625; %hz      312.5 |  625 | 1250 | 2500 |  5000
hiFreq = 2500; %hz     625   | 1250 | 2500 | 5000 | 10000
freqWidth = 1; %octave
nFreq = 7; %levels
toneDur = 40; %ms
toneSOA = 10; %ms
trialDur = 2000; %ms
freqType = 'HNH';
changeTime = [300 1000]; %ms
changeRamp = 500; %ms

%% check inputs
if any(changeTime>trialDur), disp('invalid changeTime input: later than trial length'); return; end
if length(changeTime) ~= length(freqType)-1, disp('invalid changeTime or freqType input: lengths do not match'); return; end
if any(diff(changeTime)<changeRamp), disp('invalid changeTime interval: at least one interval is smaller than ramping time'); return; end

%% generate sequence of frequencies
% assign intermediate frequency (noise)
noFreq = loFreq * 10.^(log10(2) * ceil(log10(hiFreq/loFreq)/log10(2))/2);
% add first time bin to changeTime
changeTime = [1, changeTime];

%% gradual change
toneDurTot = toneSOA+toneDur;
numTonePresentations = floor(trialDur/toneDurTot);

% convert times (ms) to coherence bins - depending on tone length
changeBins = unique(ceil((changeTime+toneSOA+toneDur/2)/toneDurTot));

freq = zeros(1,numTonePresentations);
ind_begin = 1;
for ff = 1:length(freqType)  
    switch freqType(ff)
        case 'H', tmpFreq = hiFreq;       % high frequency
        case 'L', tmpFreq = loFreq;       % low frequency
        case 'N', tmpFreq = noFreq;       % noisy signal
        case 'S', tmpFreq = nan;         % absence of signal (silence)
    end

    if ff < length(freqType)
        ind_end = changeBins(ff+1);
        freq(ind_begin:ind_end-1) = tmpFreq * ones(1,ind_end-ind_begin);
        ind_begin = ind_end;
    else
        ind_end = numTonePresentations;
        freq(ind_begin:ind_end) = tmpFreq * ones(1,ind_end-ind_begin+1);
    end
    
    
    % in case of ramping change
    if changeRamp > toneDurTot && ff < length(freqType)
        nSteps = floor(changeRamp/toneDurTot)+1;
        switch freqType(ff+1)
            case 'H'
                freqChange = ceil(log10(hiFreq/tmpFreq)/log10(2));
                rampFreq = tmpFreq * 10.^(log10(2) * (0:freqChange/(nSteps-1):freqChange));
            case 'L'
                freqChange = ceil(log10(loFreq/tmpFreq)/log10(2));
                rampFreq = tmpFreq * 10.^(log10(2) * (0:freqChange/(nSteps-1):freqChange));
            case 'N'
                freqChange = ceil(log10(noFreq/tmpFreq)/log10(2));
                rampFreq = tmpFreq * 10.^(log10(2) * (0:freqChange/(nSteps-1):freqChange));
            case 'S', rampFreq = repmat(tmpFreq,1,nSteps+1); % no ramping in case of silent pause
        end
        % take only the transitional steps - remove first (current) & last (next) coherence levels
        rampFreq = rampFreq(2:end-1);
        ind_end = ind_begin + length(rampFreq) - 1;
        freq(ind_begin:ind_end) = rampFreq;
        ind_begin = ind_end + 1;
    end
end

% calc no. of tone presentations
bufferLen = ceil(trialDur/toneDurTot);
frequencyBuffer = nan(bufferLen,1);

%% generate frequency clouds

freqChanges = [1 find(diff(freq))+1 numTonePresentations+1];
ind = 0;
for ff = 2:length(freqChanges)
    ind_freq = freqChanges(ff-1);
    numToneinBloc = freqChanges(ff) - ind_freq;
    if ~isnan(freq(ind_freq))
        % compute freq distribution
        curFreq = freq(ind_freq);
        minFreq = curFreq * 10.^(log10(2) * -(freqWidth));
        maxFreq = curFreq * 10.^(log10(2) * freqWidth);
        freqRange = minFreq * 10.^(log10(2) * (0:(2*freqWidth)/(nFreq-1):ceil(log10(maxFreq/minFreq)/log10(2))));
        
        frequencyBuffer(ind+1:ind+numToneinBloc,1) = freqRange(randi(nFreq,1,numToneinBloc));
    else
        frequencyBuffer(ind+1:ind+numToneinBloc,1) = nan;
    end
    ind = ind + numToneinBloc;
end

numLo = sum(frequencyBuffer < noFreq);
numHi = sum(frequencyBuffer > noFreq);

%% generate sequence of tones
% fs = 384000;%384000;%12000;
rampLength = 5;
cosRamp = (rampLength/1000)*fs;

multiplier = floor(fs/1000);
s = zeros(1,trialDur*multiplier);
f = nan(1,trialDur*multiplier);
ind = 1;
for ii = 1:numTonePresentations
%     ind = (ii-1)*(toneSOA+toneDur)*multiplier + 1;
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
maxFreq = ceil(max(frequencyBuffer));
%% for testing - listen
sound(s, fs) 

% record
filename = 'sample_stim_toneClouds_movingDist.wav';
audiowrite(filename,s,fs)
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
temp_fig_path = '/dataAnalysis/git_public/Penn_auditoryDecision/stimuli/';
fig_name = 'sample_stim_toneClouds_movingDist_plot';
h = figure('Name',fig_name,'Position',get(0,'ScreenSize'));
subplot(2,1,1)
plot(td,f/1000)
line([td(1) td(end)], [loFreq/1000 loFreq/1000],'LineStyle','--','Color',[0.5 0.5 0.5])
line([td(1) td(end)], [noFreq/1000 noFreq/1000],'LineStyle',':','Color',[0.5 0.5 0.5])
line([td(1) td(end)], [hiFreq/1000 hiFreq/1000],'LineStyle','--','Color',[0.5 0.5 0.5])
text(td(end),(hiFreq+100)/1000,['numHi = ' num2str(numHi)],'HorizontalAlignment','right');
text(td(end),(loFreq-100)/1000,['numLo = ' num2str(numLo)],'HorizontalAlignment','right');
xlabel('Time (secs)')
ylabel('Frequency(Hz)')
ylim([0 ceil(maxFreq/1000)])
% 
subplot(2,1,2)
segmentLength = round(numel(s)/5);
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

% spectrogram(s,segmentLength,[],[],fs,'yaxis')
spectrogram(s,round(segmentLength/5),round(80/100*segmentLength/5),[],fs,'yaxis')
% spectrogram(s,round(numel(s)/25),[],round(numel(s)/25),fs,'yaxis')
ylim([0 ceil(maxFreq/1000)])
saveas(h,[temp_fig_path fig_name '.fig']);
save2pdf([temp_fig_path fig_name '.pdf'],h);
close(gcf)