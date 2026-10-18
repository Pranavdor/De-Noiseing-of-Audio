clc; clear; close all; rng(1);

%  1. USER SETTINGS
INPUT_FILE   = 'sp01_car_sn5.wav';      % path of NOISY audio; '' = open file dialog
REF_FILE     = '';      % OPTIONAL: matching CLEAN file (enables true SNR improvement)
TARGET_FS    = 16000;   % resample to this rate in Hz (0 = keep original rate)
MAX_DURATION = 30;      % use at most this many seconds of the file
ADD_SYNTHETIC_NOISE = false;  % true = treat input as CLEAN and add white+hum+rumble
SYN_SNR_DB   = 5;       % SNR (dB) of the synthetic noise
MAINS_HZ     = 50;      % mains frequency: 50 Hz (India/Europe), 60 Hz (USA)

F_C1     = 100;         % lower -6 dB cut-off (Hz)
F_C2     = 4000;        % upper -6 dB cut-off (Hz)
TRANS_BW = 80;          % FIR transition width (Hz): smaller = longer, sharper filter
MAX_TAPS = 6000;        % safety limit on FIR order

% Window-specific parameter(s)

% IIR band-pass settings
IIR_TYPE       = 'ellip';   % 'butter' | 'cheby1' | 'cheby2' | 'ellip'
IIR_RP         = 1;         % pass-band ripple (dB)
IIR_RS         = 40;        % stop-band attenuation (dB)
IIR_TRANS_HI   = 500;       % upper transition width (Hz). Butterworth needs this LARGE
                            % (>= ~1500 Hz) or its order explodes; ellip works with 500.
IIR_ZERO_PHASE = true;      % true: filtfilt (zero phase, |H|^2); false: causal filter

USE_NOTCH     = true;
NUM_HARMONICS = 3;          % 50, 100, 150 Hz
NOTCH_Q       = 35;         % higher Q = narrower notch

PLAY_AUDIO = false;         % true = play noisy then de-noised audio

%  2. LOAD AND PREPARING THE INPUT AUDIO
if isempty(INPUT_FILE)
    [fn_,pth_] = uigetfile({'*.wav;*.flac;*.mp3;*.m4a;*.ogg','Audio files'}, ...
                           'Select NOISY audio file (Cancel = built-in demo)');
    if ~isequal(fn_,0), INPUT_FILE = fullfile(pth_,fn_); end
end

if isempty(INPUT_FILE)                       % ---- demo mode ----
    S_ = load('handel');  x = S_.y(:);  fs = S_.Fs;
    if TARGET_FS > 0 && fs ~= TARGET_FS, x = resampleTo(x,fs,TARGET_FS); fs = TARGET_FS; end
    srcName = 'MATLAB demo clip (handel.mat)';
    ADD_SYNTHETIC_NOISE = true;
else
    [x,fs] = loadMono(INPUT_FILE,TARGET_FS);
    srcName = INPUT_FILE;
end
x = x(1:min(numel(x),round(MAX_DURATION*fs)));

ref = [];
if ~isempty(REF_FILE)
    ref = loadMono(REF_FILE,fs);
    Lm = min(numel(x),numel(ref));  x = x(1:Lm);  ref = ref(1:Lm);
end

sc = max(abs(x)) + eps;   x = x/sc;   if ~isempty(ref), ref = ref/sc; end

if ADD_SYNTHETIC_NOISE                       % build a noisy version of a clean clip
    ref = x;                                 % the clean clip becomes the reference
    tt  = (0:numel(x)-1)'/fs;
    white  = randn(size(x));
    hum    = sin(2*pi*MAINS_HZ*tt) + 0.5*sin(2*pi*2*MAINS_HZ*tt) + 0.25*sin(2*pi*3*MAINS_HZ*tt);
    rumble = filter(1,[1 -0.995],randn(size(x)));          % low-frequency rumble
    noise  = white/std(white) + 0.6*hum/std(hum) + 0.5*rumble/std(rumble);
    noise  = noise/std(noise) * sqrt(mean(ref.^2)/10^(SYN_SNR_DB/10));
    x = ref + noise;
end

F_C2 = min(F_C2, fs/2 - 2*TRANS_BW);
if F_C2 <= F_C1, error('F_C2 must be above F_C1 and below Nyquist.'); end
L = numel(x);   t = (0:L-1)'/fs;
hum_f = MAINS_HZ*(1:NUM_HARMONICS);   hum_f = hum_f(hum_f < 0.45*fs);

fprintf('\n================ AUDIO DE-NOISING : Hanning (Hann) WINDOW ================\n');
fprintf('Input            : %s\n',srcName);
fprintf('Sampling rate    : %d Hz   | Duration: %.2f s | Samples: %d\n',fs,L/fs,L);
fprintf('Band to keep     : %.0f - %.0f Hz\n',F_C1,F_C2);

%  3. ANALYSIS OF THE NOISY INPUT  (Figure 1)
[Sx,~,~] = spectrogram(x,hamming(512),384,1024,fs);
Sdb = 10*log10(abs(Sx).^2+eps);   CL = [max(Sdb(:))-90, max(Sdb(:))];   % common colour scale

newFig('1_Input_Signal');
subplot(2,2,1);  plot(t,x,'Color',[0.85 0.2 0.2]); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Noisy input - time domain');
[Xf,ff] = oneSidedFFT(x,fs);
subplot(2,2,2);  plot(ff,20*log10(Xf+eps)); grid on; hold on;
for k = 1:numel(hum_f), xline(hum_f(k),':k'); end
xlim([0 fs/2]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Input spectrum (dotted = mains hum)');
subplot(2,2,3);  plotSpec(x,fs,CL,'Input spectrogram');
subplot(2,2,4);  [pxx,fpx] = pwelch(x,hamming(1024),512,1024,fs);
plot(fpx,10*log10(pxx+eps)); grid on; xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)'); title('Input power spectral density (Welch)');
sgtitle('Step 1 - Noisy input signal analysis');

%  4. WINDOW FUNCTION: LENGTH, SHAPE AND SPECTRUM  (Figure 2)
C = 3.1;     % empirical constant: transition width ~ C*fs/N  (window dependent)
N = ceil(C*fs/TRANS_BW);
N = N + mod(N,2);            % force EVEN order -> ODD length -> Type-I linear phase
N = min(N,MAX_TAPS);  N = N + mod(N,2);
w = hann(N+1);   w = w(:);

Nw = 2^16;  Wf = fft(w,Nw);  Wm = abs(Wf(1:Nw/2+1));  Wm = Wm/Wm(1);
Wdb = 20*log10(Wm+eps);  fw = (0:Nw/2)'/Nw;                    % cycles/sample
i0 = find(diff(Wm)>0,1);                                        % first spectral null
if isempty(i0), i0 = numel(Wm); mainlobe = NaN; psl = -Inf;
else, mainlobe = fw(i0)*(N+1); psl = max(Wdb(i0:end)); end
i3 = find(Wdb<-3,1);  bw3 = 2*fw(i3)*(N+1);

fprintf('\n--- Window: Hanning (Hann) ---\n');
fprintf('FIR order N = %d  (length %d taps)   Group delay = %d samples = %.2f ms\n',N,N+1,N/2,1000*(N/2)/fs);
fprintf('Main-lobe half-width = %.2f bins | -3 dB width = %.2f bins | Peak side-lobe = %.1f dB\n',mainlobe,bw3,psl);

newFig('2_Window_Function');
subplot(1,3,1); plot(0:N,w,'LineWidth',1.5); grid on; xlabel('n'); ylabel('w[n]'); title('Hanning (Hann) window (time domain)');
subplot(1,3,2); plot(fw,Wdb,'LineWidth',1.2); grid on; ylim([-140 5]);
xlabel('Normalised frequency (cycles/sample)'); ylabel('Magnitude (dB)'); title('Window spectrum');
subplot(1,3,3); plot(fw*(N+1),Wdb,'LineWidth',1.2); grid on; xlim([0 8]); ylim([-140 5]);
xlabel('Frequency (bins)'); ylabel('Magnitude (dB)'); title(sprintf('Main lobe zoom (peak side-lobe %.1f dB)',psl));
sgtitle('Step 2 - Window function and its frequency response');

%  5. FIR BAND-PASS DESIGN BY THE WINDOW METHOD  (Figures 3 and 4)
n  = (0:N)';   m = n - N/2;                      % centre the impulse response
fc1 = F_C1/fs;  fc2 = F_C2/fs;                   % cycles/sample
h_ideal = 2*fc2*mysinc(2*fc2*m) - 2*fc1*mysinc(2*fc1*m);   % ideal band-pass (infinite sinc, truncated)
h_fir   = h_ideal .* w;                          % <-- THE WINDOWING STEP

h_ref = fir1(N,[F_C1 F_C2]/(fs/2),'bandpass',w,'noscale');  % MATLAB's own result
fprintf('Check: max |h_manual - h_fir1| = %.2e  (should be ~1e-16)\n',max(abs(h_ref(:)-h_fir(:))));

newFig('3_FIR_Design_Steps');
subplot(2,2,1); plot(n,h_ideal,'LineWidth',1.1); grid on; xlabel('n'); title('(a) Ideal band-pass impulse response (sinc difference)');
subplot(2,2,2); plot(n,w,'r','LineWidth',1.3);   grid on; xlabel('n'); title('(b) Hanning (Hann) window');
subplot(2,2,3); plot(n,h_fir,'Color',[0 0.5 0],'LineWidth',1.1); grid on; xlabel('n'); title('(c) Windowed impulse response  h[n] = h_{ideal}[n]·w[n]');
subplot(2,2,4); plot(n,h_ideal,'Color',[.6 .6 .6]); hold on; plot(n,h_fir,'Color',[0 0.5 0]); grid on;
legend('Truncated ideal','Windowed','Location','best'); xlabel('n'); title('(d) Effect of windowing on the taps');
sgtitle('Step 3 - FIR window-method design, step by step');

[H_w,fh]  = freqz(h_fir,1,16384,fs);
H_id      = freqz(h_ideal,1,16384,fs);
stopMask  = (fh > F_C2+TRANS_BW) | (fh < F_C1-TRANS_BW);
passMask  = (fh > F_C1+TRANS_BW) & (fh < F_C2-TRANS_BW);
stopAtt   = -20*log10(max(abs(H_w(stopMask)))+eps);
passRip   = max(abs(20*log10(abs(H_w(passMask))+eps)));
fprintf('FIR: stop-band attenuation = %.1f dB | pass-band ripple = +/-%.3f dB\n',stopAtt,passRip);

newFig('4_FIR_Frequency_Response');
subplot(2,2,1); plot(fh,20*log10(abs(H_id)+eps),'Color',[.6 .6 .6]); hold on;
plot(fh,20*log10(abs(H_w)+eps),'b','LineWidth',1.4); ylim([-140 10]); grid on;
xline(F_C1,'--r'); xline(F_C2,'--r'); legend('Rectangular (no window)','Hanning (Hann) window','Location','south');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Magnitude response');
subplot(2,2,2); plot(fh,20*log10(abs(H_w)+eps),'b','LineWidth',1.2); grid on;
xlim([F_C1+TRANS_BW F_C2-TRANS_BW]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Pass-band detail (ripple)');
subplot(2,2,3); plot(fh,unwrap(angle(H_w)),'LineWidth',1.1); grid on;
xlabel('Frequency (Hz)'); ylabel('Phase (rad)'); title('Phase response (linear)');
subplot(2,2,4); gdF = gdFromH(H_w,fh,fs); plot(fh,gdF,'LineWidth',1.3); grid on; ylim([0 N]);
xlabel('Frequency (Hz)'); ylabel('Samples'); title(sprintf('Group delay = N/2 = %d samples (constant)',N/2));
sgtitle('Step 3 (cont.) - FIR frequency-domain characteristics');

%  6. IIR BAND-PASS DESIGN  (Figure 5)
Fp  = [F_C1 F_C2];
Fst = [F_C1/2, min(F_C2+IIR_TRANS_HI,0.98*fs/2)];
[sos,zi,pi_,nI] = designIIR(IIR_TYPE,Fp,Fst,IIR_RP,IIR_RS,fs);
nSec = size(sos,1);
fprintf('\n--- IIR band-pass (%s) ---\n',IIR_TYPE);
fprintf('Prototype order n = %d -> band-pass order %d (%d biquad sections)\n',nI,2*nI,nSec);
fprintf('Stable (all poles inside unit circle): %d\n',all(abs(pi_)<1));
if nI > 16, warning('Very high IIR order. Use IIR_TYPE=''ellip'' or increase IIR_TRANS_HI.'); end

[H_i,fi] = freqz(sos,16384,fs);
newFig('5_IIR_Design');
subplot(2,2,1); plot(fi,20*log10(abs(H_i)+eps),'LineWidth',1.3); hold on;
plot(fi,40*log10(abs(H_i)+eps),'--','LineWidth',1.1); ylim([-120 10]); grid on;
xline(Fp(1),':g'); xline(Fp(2),':g'); xline(Fst(1),':r'); xline(Fst(2),':r'); yline(-IIR_RS,':k');
legend('Single pass |H|','Zero-phase |H|^2','Location','south'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title(sprintf('Magnitude (%s, green = pass edges, red = stop edges)',IIR_TYPE));
subplot(2,2,2); plot(fi,unwrap(angle(H_i)),'LineWidth',1.1); grid on; xlabel('Frequency (Hz)'); ylabel('Phase (rad)'); title('Phase (non-linear for IIR)');
subplot(2,2,3); plot(fi,gdFromH(H_i,fi,fs),'LineWidth',1.1); grid on; xlabel('Frequency (Hz)'); ylabel('Samples'); title('Group delay (varies with frequency)');
subplot(2,2,4); zplane(zi,pi_); grid on; title('Pole-zero plot');
sgtitle('Step 4 - IIR band-pass design');

%  7. IIR NOTCH FILTERS FOR MAINS HUM  (Figure 6)
y1 = x;
newFig('6_IIR_Notch_Filters');  subplot(1,2,1); hold on;
if USE_NOTCH
    for k = 1:numel(hum_f)
        w0 = hum_f(k)/(fs/2);  [bn,an] = iirnotch(w0,w0/NOTCH_Q);
        y1 = filtfilt(bn,an,y1);
        [Hn,fn2] = freqz(bn,an,16384,fs);  plot(fn2,20*log10(abs(Hn)+eps),'LineWidth',1.2);
    end
end
grid on; xlim([0 400]); ylim([-60 5]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Notch responses (hum and harmonics)');
subplot(1,2,2); [Xa,fa] = oneSidedFFT(x,fs); [Xb,fb] = oneSidedFFT(y1,fs);
plot(fa,20*log10(Xa+eps)); hold on; plot(fb,20*log10(Xb+eps)); grid on; xlim([0 400]);
legend('Before notch','After notch'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Spectrum near the hum frequencies');
sgtitle('Step 5 - IIR notch filtering');

%  8. APPLY FIR (delay-compensated) THEN IIR BAND-PASS
D  = N/2;
tmp = filter(h_fir,1,[y1; zeros(D,1)]);   y2 = tmp(D+1:end);    % FIR output, delay removed
y3 = applySOS(sos,y2,IIR_ZERO_PHASE);                            % final output
y_out = y3;   removed = x - y_out;

%  9. FIR vs IIR COMPARISON  (Figure 7)
newFig('7_FIR_vs_IIR');
subplot(1,2,1); plot(fh,20*log10(abs(H_w)+eps),'b','LineWidth',1.3); hold on; plot(fi,20*log10(abs(H_i)+eps),'r','LineWidth',1.3);
ylim([-120 10]); grid on; legend(sprintf('FIR Hanning (Hann) (%d taps)',N+1),sprintf('IIR %s (order %d)',IIR_TYPE,2*nI),'Location','south');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Magnitude');
subplot(1,2,2); plot(fh,gdF,'b','LineWidth',1.3); hold on; plot(fi,gdFromH(H_i,fi,fs),'r','LineWidth',1.3); grid on;
set(gca,'YScale','log'); xlabel('Frequency (Hz)'); ylabel('Samples (log)'); title('Group delay'); legend('FIR','IIR');
sgtitle(sprintf('FIR vs IIR cost: FIR = %d mult/sample (linear phase) | IIR = %d mult/sample (non-linear phase)',N+1,5*nSec));

%  10. STAGE-BY-STAGE SIGNAL PLOTS  (Figures 8 - 15)
stageNames = {'Input (noisy)','After IIR notch','After FIR (Hanning (Hann))','Final (after IIR band-pass)'};
S  = [x, y1, y2, y3];   cols = lines(4);   yl = 1.05*max(abs(S(:)));

newFig('8_Stage_by_Stage_Waveforms');
for k = 1:4
    subplot(4,1,k); plot(t,S(:,k),'Color',cols(k,:)); ylim([-yl yl]); grid on; ylabel('Amplitude'); title(stageNames{k});
end
xlabel('Time (s)');  sgtitle('Step 6 - Waveform after every stage');

[~,i0m] = max(abs(x));  seg = max(1,i0m-round(0.03*fs)) : min(L,i0m+round(0.03*fs));
newFig('9_Zoomed_Waveform_and_Removed_Noise');
subplot(2,1,1); plot(t(seg),x(seg),'Color',[.8 .2 .2]); hold on; plot(t(seg),y_out(seg),'Color',[0 .4 .8],'LineWidth',1.3);
grid on; legend('Noisy input','De-noised output'); xlabel('Time (s)'); ylabel('Amplitude'); title('60 ms zoom: input vs output');
subplot(2,1,2); plot(t(seg),removed(seg),'Color',[.4 .4 .4]); grid on; xlabel('Time (s)'); ylabel('Amplitude'); title('Removed component (input - output)');

newFig('10_Spectrum_After_Each_Stage');
subplot(2,1,1); hold on;
for k = 1:4, [Xk,fk] = oneSidedFFT(S(:,k),fs); plot(fk,20*log10(Xk+eps),'Color',cols(k,:)); end
grid on; xlim([0 fs/2]); legend(stageNames,'Location','northeast'); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Magnitude spectrum');
subplot(2,1,2); hold on;
for k = 1:4, [Xk,fk] = oneSidedFFT(S(:,k),fs); plot(fk,20*log10(Xk+eps),'Color',cols(k,:)); end
grid on; xlim([0 600]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Low-frequency zoom (hum / rumble region)');

newFig('11_Spectrograms');
for k = 1:4, subplot(2,2,k); plotSpec(S(:,k),fs,CL,stageNames{k}); end
sgtitle('Spectrograms - same colour scale for all stages');

newFig('12_PSD_Comparison');
hold on;
for k = 1:4, [pk_,fpk] = pwelch(S(:,k),hamming(1024),512,1024,fs); plot(fpk,10*log10(pk_+eps),'Color',cols(k,:),'LineWidth',1.1); end
[pk_,fpk] = pwelch(removed,hamming(1024),512,1024,fs); plot(fpk,10*log10(pk_+eps),'k--','LineWidth',1.1);
grid on; legend([stageNames,{'Removed noise'}],'Location','northeast'); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)'); title('Welch power spectral density');

newFig('13_Removed_Noise_Analysis');
subplot(1,3,1); plot(t,removed,'Color',[.4 .4 .4]); grid on; xlabel('Time (s)'); ylabel('Amplitude'); title('Removed noise - time');
subplot(1,3,2); [pn_,fpn] = pwelch(removed,hamming(1024),512,1024,fs); plot(fpn,10*log10(pn_+eps)); grid on; xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)'); title('Removed noise - PSD');
subplot(1,3,3); histogram(removed,100,'Normalization','pdf'); grid on; xlabel('Amplitude'); ylabel('PDF'); title('Removed noise - histogram');

fl_ = round(0.025*fs);  hp_ = round(0.0125*fs);
Fin = frameSignal(x,fl_,hp_);  Fout = frameSignal(y_out,fl_,hp_);  tf_ = ((0:size(Fin,2)-1)*hp_ + fl_/2)/fs;
newFig('14_Short_Time_Energy_and_ZCR');
subplot(2,1,1); plot(tf_,10*log10(mean(Fin.^2,1)+eps),'r'); hold on; plot(tf_,10*log10(mean(Fout.^2,1)+eps),'b','LineWidth',1.2);
grid on; legend('Input','Output'); xlabel('Time (s)'); ylabel('Energy (dB)'); title('Short-time energy (25 ms frames)');
subplot(2,1,2); plot(tf_,sum(Fin(1:end-1,:).*Fin(2:end,:)<0,1)/fl_,'r'); hold on; plot(tf_,sum(Fout(1:end-1,:).*Fout(2:end,:)<0,1)/fl_,'b','LineWidth',1.2);
grid on; legend('Input','Output'); xlabel('Time (s)'); ylabel('Crossings / sample'); title('Short-time zero-crossing rate');

newFig('15_Amplitude_Histograms');
histogram(x,120,'Normalization','pdf','DisplayStyle','stairs','EdgeColor','r','LineWidth',1.2); hold on;
histogram(y_out,120,'Normalization','pdf','DisplayStyle','stairs','EdgeColor','b','LineWidth',1.2);
grid on; legend('Input','Output'); xlabel('Amplitude'); ylabel('PDF'); title('Amplitude distribution before and after');

%  11. SIGNAL FEATURES FOR EVERY STAGE  (table + Figure 16)
V = zeros(25,4);
for k = 1:4, [V(:,k),featNames] = featureVector(S(:,k),fs); end
T = array2table(V,'RowNames',featNames,'VariableNames',{'Input','After_Notch','After_FIR','Final_Output'});
T.Delta_Final_minus_Input = T.Final_Output - T.Input;
fprintf('\n================ SIGNAL FEATURES (Hanning (Hann) window) ================\n');
disp(T);

pick = {'Mean','Standard deviation','Average power','RMS value','Estimated SNR (dB)', ...
        'Spectral centroid (Hz)','Spectral flatness','Zero-crossing rate (per s)','Kurtosis'};
newFig('16_Feature_Comparison');
for k = 1:numel(pick)
    subplot(3,3,k); idx = find(strcmp(featNames,pick{k}));
    b = bar(V(idx,:)); b.FaceColor = 'flat'; for q = 1:4, b.CData(q,:) = cols(q,:); end
    xticklabels({'In','Notch','FIR','Final'}); title(pick{k}); grid on;
end
sgtitle('Feature comparison across the processing chain');

%  12. REFERENCE-BASED METRICS (Figure 17)
if ~isempty(ref)
    snrK = zeros(1,4); mseK = zeros(1,4); corK = zeros(1,4);
    for k = 1:4
        e = ref - S(:,k);
        snrK(k) = 10*log10(sum(ref.^2)/sum(e.^2));   mseK(k) = mean(e.^2);
        cc = corrcoef(ref,S(:,k));                    corK(k) = cc(1,2);
    end
    TR = table(snrK',mseK',corK','RowNames',{'Input','After_Notch','After_FIR','Final_Output'}, ...
               'VariableNames',{'SNR_dB','MSE','Correlation'});
    fprintf('\n---------- Metrics against the CLEAN reference ----------\n'); disp(TR);
    fprintf('SNR IMPROVEMENT = %.2f dB\n',snrK(4)-snrK(1));
    newFig('17_Reference_Comparison');
    subplot(2,1,1); plot(t,ref,'Color',[.2 .6 .2]); hold on; plot(t,x,'Color',[.8 .2 .2]); plot(t,y_out,'Color',[0 .3 .8]);
    grid on; legend('Clean reference','Noisy input','De-noised'); xlabel('Time (s)'); title('Waveforms');
    subplot(2,1,2); bar(snrK); xticklabels({'Input','Notch','FIR','Final'}); grid on; ylabel('SNR (dB)');
    title(sprintf('SNR vs clean reference (improvement = %.2f dB)',snrK(4)-snrK(1)));
end

%  13. SAVE AUDIO, TABLES AND FIGURES
outDir = fullfile(pwd,'Results_Hanning');  if ~exist(outDir,'dir'), mkdir(outDir); end
nrm = @(s) 0.95*s/(max(abs(s))+eps);
audiowrite(fullfile(outDir,'noisy_input.wav'),nrm(x),fs);
audiowrite(fullfile(outDir,'denoised_Hanning.wav'),nrm(y_out),fs);
audiowrite(fullfile(outDir,'removed_noise_Hanning.wav'),nrm(removed),fs);
writetable(T,fullfile(outDir,'features_Hanning.csv'),'WriteRowNames',true);
figs = findall(0,'Type','figure');
for k = 1:numel(figs)
    nm = regexprep(figs(k).Name,'[^\w]','_');  if isempty(nm), nm = sprintf('Figure_%d',figs(k).Number); end
    saveas(figs(k),fullfile(outDir,[nm '.png']));
end
fprintf('\nDONE. De-noised audio, feature table and all figures saved in:\n  %s\n',outDir);

if PLAY_AUDIO
    sound(nrm(x),fs);      pause(L/fs+1);
    sound(nrm(y_out),fs);  pause(L/fs+1);
end


%  LOCAL FUNCTIONS
function [x,fs] = loadMono(file,targetFs)
    [x,fs] = audioread(file);
    x = mean(x,2);                                  % stereo -> mono
    if targetFs > 0 && fs ~= targetFs
        x = resampleTo(x,fs,targetFs);  fs = targetFs;
    end
end

function y = resampleTo(x,fs,targetFs)
    [p,q] = rat(targetFs/fs);
    y = resample(x,p,q);
end

function s = mysinc(u)                              % sin(pi u)/(pi u), with sinc(0)=1
    s = ones(size(u));
    k = (u ~= 0);
    s(k) = sin(pi*u(k))./(pi*u(k));
end

function f = newFig(name)
    f = figure('Name',name,'NumberTitle','off','Color','w', ...
               'Units','normalized','Position',[0.05 0.08 0.9 0.8]);
end

function [X,f] = oneSidedFFT(x,fs)
    L = numel(x);  Nfft = 2^nextpow2(L);
    X = abs(fft(x,Nfft))/L;  X = X(1:Nfft/2+1);
    f = (0:Nfft/2)'*fs/Nfft;
end

function plotSpec(sig,fs,CL,ttl)
    [S,f,t] = spectrogram(sig,hamming(512),384,1024,fs);
    imagesc(t,f/1000,10*log10(abs(S).^2+eps));  axis xy;
    set(gca,'CLim',CL);  colorbar;
    xlabel('Time (s)');  ylabel('Frequency (kHz)');  title(ttl);
end

function F = frameSignal(x,fl,hop)                  % columns = frames
    x = x(:);  nF = floor((numel(x)-fl)/hop)+1;
    idx = (1:fl)' + (0:nF-1)*hop;
    F = x(idx);
end

function y = applySOS(sos,x,zeroPhase)              % IIR in second-order sections
    if zeroPhase
        y = filtfilt(sos,1,x);                      % forward-backward: zero phase
    else
        y = x;
        for k = 1:size(sos,1)
            y = filter(sos(k,1:3),sos(k,4:6),y);
        end
    end
end

function gd = gdFromH(H,f,fs)                       % group delay (samples) from H(f)
    ph = unwrap(angle(H));
    gd = -diff(ph)./diff(2*pi*f/fs);
    gd = [gd; gd(end)];
end

function [sos,z,p,n] = designIIR(type,Fp,Fst,Rp,Rs,fs)
    Wp = Fp/(fs/2);  Ws = Fst/(fs/2);               % normalised to Nyquist
    switch lower(type)
        case 'butter'
            [n,Wn] = buttord(Wp,Ws,Rp,Rs);   [z,p,k] = butter(n,Wn,'bandpass');
        case 'cheby1'
            [n,Wn] = cheb1ord(Wp,Ws,Rp,Rs);  [z,p,k] = cheby1(n,Rp,Wn,'bandpass');
        case 'cheby2'
            [n,Wn] = cheb2ord(Wp,Ws,Rp,Rs);  [z,p,k] = cheby2(n,Rs,Wn,'bandpass');
        case 'ellip'
            [n,Wn] = ellipord(Wp,Ws,Rp,Rs);  [z,p,k] = ellip(n,Rp,Rs,Wn,'bandpass');
        otherwise
            error('IIR_TYPE must be butter, cheby1, cheby2 or ellip');
    end
    [sos,g] = zp2sos(z,p,k);                        % stable cascade of biquads
    sos(1,1:3) = sos(1,1:3)*g;                      % fold the gain into section 1
end

function [v,names] = featureVector(x,fs)
    x = x(:);  L = numel(x);
    mu = mean(x);  sd = std(x);
    Pavg = mean(x.^2);  E = sum(x.^2);
    sk = mean((x-mu).^3)/(sd^3+eps);                % skewness
    ku = mean((x-mu).^4)/(sd^4+eps);                % kurtosis
    zc = sum(x(1:end-1).*x(2:end) < 0)/(L/fs);      % zero crossings per second
    [X,f] = oneSidedFFT(x,fs);  P = X.^2;  Pn = P/(sum(P)+eps);
    cen  = sum(f.*Pn);                              % spectral centroid
    bw   = sqrt(sum(((f-cen).^2).*Pn));             % spectral bandwidth
    cs   = cumsum(Pn);  roll = f(find(cs>=0.85,1)); % 85 % roll-off
    flat = exp(mean(log(P+eps)))/(mean(P)+eps);     % spectral flatness (1 = white)
    [~,im] = max(P(2:end));  domf = f(im+1);        % dominant frequency (no DC)
    ent  = -sum(Pn.*log2(Pn+eps))/log2(numel(Pn));  % normalised spectral entropy
    fl = round(0.025*fs);  F = frameSignal(x,fl,fl);% 25 ms frames
    En = sort(mean(F.^2,1));  k = max(1,round(0.1*numel(En)));
    Pnoise = mean(En(1:k));                         % noise floor = quietest 10 % frames
    snr = 10*log10(max(Pavg-Pnoise,eps)/(Pnoise+eps));
    names = {'Mean','Median','Standard deviation','Variance','Mean absolute value', ...
             'RMS value','Average power','Average power (dB)','Total energy', ...
             'Peak amplitude','Minimum','Maximum','Peak-to-peak','Crest factor', ...
             'Skewness','Kurtosis','Zero-crossing rate (per s)','Spectral centroid (Hz)', ...
             'Spectral bandwidth (Hz)','Spectral roll-off 85% (Hz)','Spectral flatness', ...
             'Dominant frequency (Hz)','Spectral entropy (norm.)','Noise-floor power (est.)', ...
             'Estimated SNR (dB)'};
    v = [mu; median(x); sd; sd^2; mean(abs(x)); sqrt(Pavg); Pavg; 10*log10(Pavg+eps); E; ...
         max(abs(x)); min(x); max(x); max(x)-min(x); max(abs(x))/(sqrt(Pavg)+eps); ...
         sk; ku; zc; cen; bw; roll; flat; domf; ent; Pnoise; snr];
end
