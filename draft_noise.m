fs = 2000;

fn = 0.5*fs;
fd = 240;
fu = 435;

white_noise = randn(1e4, 1);

[B, A] = fir1(10,[fd fu]/fn);

fvtool(B,A,'Fs',fs)
noise_limited = filter(B, A, white_noise);


[B1,A1] = fir1(48,[fd fu]/fn);
fvtool(B1,A1,'Fs',fs)
noise_limited = filter(B1,A1,white_noise);

sound(noise_limited)


d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', 400,500,40000,45000,60,0.5,60,200000);
Hd = design(d);
y = filter(Hd,randn(2e5,1));
fvtool(Hd)
sound(y,200000)



% design FIR filter to filter noise to half of Nyquist rate
b = fir1(64, 0.5);
% generate Gaussian (normally-distributed) white noise
n = randn(1e4, 1);
% apply to filter to yield bandlimited noise
nb = filter(b,1,n);
fvtool(nb)

var = 3.0;  % just an example  
scale = sqrt(var)/std(nb);
nb = scale*nb;  % nb has variance 'var'

sound(nb)


%%
fs = 2e5;
[B1,A1] = fir1(256,0.4);
fvtool(B1,A1,'Fs',fs)

white_noise = randn(fs, 1);
noise_limited = filter(B1,A1,white_noise);

% sound(white_noise,fs)
% sound(noise_limited,fs)
%%

loFreq = 500; %hz      312.5 |  625 | 1250 | 2500 |  5000
hiFreq = 2000; %hz     625   | 1250 | 2500 | 5000 | 10000 

toneDur = 500; %ms
toneIBI = 100; %ms

rampLength = 10;
cosRamp = (rampLength/1000)*fs;


%%
freqType = 'HHHL';
frequencyBuffer = zeros(length(freqType),1);
lastTone_amp = 0.2;
noise_amp = 1;

trialDur = (toneDur*4) + (toneIBI*3);
multiplier = floor(fs/1000);
s = zeros(1,trialDur*multiplier);
f = nan(1,trialDur*multiplier);
ind = 1;

for ii = 1:length(freqType)
    tmax = toneDur/1000;                               % Time Duration Of Signal (sec)
    t = linspace(0, tmax, tmax*fs);                          % Time Vector
    
    switch freqType(ii)
        case 'H', tmpFreq = hiFreq;       % high frequency
        case 'L', tmpFreq = loFreq;       % low frequency
    end
    
    tone = sin(2*pi*tmpFreq*t);
    tone = pa_ramp(tone, cosRamp, fs);
    
    if ii == length(freqType) % if last tone
        tone = (lastTone_amp*tone)+(noise_amp*noise_limited(1:length(tone)));
    end
    
    s(ind:ind+length(tone)-1) = tone';
    f(ind:ind+length(tone)-1) = tmpFreq;
    ind = ind+length(tone)+(toneIBI/1000*fs);
end

noiser = dotsPlayableWave();
noiser.wave = s;
noiser.sampleFrequency = fs;
% noiser.duration = 1;
noiser.intensity = 0.5;
noiser.prepareToPlay;
noiser.play



