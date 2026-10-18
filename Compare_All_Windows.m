clc; clear; close all; rng(1);

% 1. SETTINGS
INPUT_FILE   = 'sp01_car_sn5.wav';   
REF_FILE = '';   
TARGET_FS = 16000;   
MAX_DURATION = 30;
ADD_SYNTHETIC_NOISE = false;   
SYN_SNR_DB = 5;   
MAINS_HZ = 50;
F_C1 = 100;  
F_C2 = 4000;  
TRANS_BW = 80;  
MAX_TAPS = 6000;
IIR_TYPE = 'ellip';  
IIR_RP = 1;  IIR_RS = 40;  IIR_TRANS_HI = 500;  IIR_ZERO_PHASE = true;
USE_NOTCH = true;  NUM_HARMONICS = 3;  NOTCH_Q = 35;
P = struct('kaiserAtt',60,'chebAtt',60,'gaussAlpha',2.5,'tukeyR',0.5);
winNames = {'Rectangular','Bartlett','Hanning','Hamming','Blackman','Kaiser'};

% 2. LOAD INPUT
if isempty(INPUT_FILE)
    [fn_,pth_] = uigetfile({'*.wav;*.flac;*.mp3;*.m4a;*.ogg','Audio files'},'Select NOISY audio (Cancel = demo)');
    if ~isequal(fn_,0), INPUT_FILE = fullfile(pth_,fn_); end
end
if isempty(INPUT_FILE)
    S_ = load('handel');  x = S_.y(:);  fs = S_.Fs;
    if TARGET_FS > 0 && fs ~= TARGET_FS, x = resampleTo(x,fs,TARGET_FS); fs = TARGET_FS; end
    ADD_SYNTHETIC_NOISE = true;
else
    [x,fs] = loadMono(INPUT_FILE,TARGET_FS);
end
x = x(1:min(numel(x),round(MAX_DURATION*fs)));
ref = [];
if ~isempty(REF_FILE)
    ref = loadMono(REF_FILE,fs);  Lm = min(numel(x),numel(ref));  x = x(1:Lm);  ref = ref(1:Lm);
end
sc = max(abs(x))+eps;  x = x/sc;  if ~isempty(ref), ref = ref/sc; end
if ADD_SYNTHETIC_NOISE
    ref = x;  tt = (0:numel(x)-1)'/fs;
    white = randn(size(x));
    hum = sin(2*pi*MAINS_HZ*tt)+0.5*sin(2*pi*2*MAINS_HZ*tt)+0.25*sin(2*pi*3*MAINS_HZ*tt);
    rumble = filter(1,[1 -0.995],randn(size(x)));
    noise = white/std(white)+0.6*hum/std(hum)+0.5*rumble/std(rumble);
    noise = noise/std(noise)*sqrt(mean(ref.^2)/10^(SYN_SNR_DB/10));
    x = ref + noise;
end
F_C2 = min(F_C2,fs/2-2*TRANS_BW);  L = numel(x);  t = (0:L-1)'/fs;
hum_f = MAINS_HZ*(1:NUM_HARMONICS);  hum_f = hum_f(hum_f<0.45*fs);

% 3. COMMON STAGES (notch + IIR design are identical for every window)
y1 = x;
if USE_NOTCH
    for k = 1:numel(hum_f)
        w0 = hum_f(k)/(fs/2);  [bn,an] = iirnotch(w0,w0/NOTCH_Q);  y1 = filtfilt(bn,an,y1);
    end
end
Fp = [F_C1 F_C2];  Fst = [F_C1/2, min(F_C2+IIR_TRANS_HI,0.98*fs/2)];
[sos,~,~,nI] = designIIR(IIR_TYPE,Fp,Fst,IIR_RP,IIR_RS,fs);
fprintf('IIR %s band-pass: order %d\n',IIR_TYPE,2*nI);

% 4. LOOP OVER ALL WINDOWS
nW = numel(winNames);  featAll = zeros(25,nW);  info = zeros(nW,5);  FIRresp = cell(nW,1);  Y = cell(nW,1);
outDir = fullfile(pwd,'Results_Compare_All_Windows');  if ~exist(outDir,'dir'), mkdir(outDir); end
nrm = @(s) 0.95*s/(max(abs(s))+eps);
snrRef = nan(nW,1);
for q = 1:nW
    name = winNames{q};
    N = windowOrder(name,fs,TRANS_BW,P);  N = N+mod(N,2);  N = min(N,MAX_TAPS);  N = N+mod(N,2);
    w = makeWin(name,N,P);
    m = (0:N)'-N/2;  fc1 = F_C1/fs;  fc2 = F_C2/fs;
    h = (2*fc2*mysinc(2*fc2*m)-2*fc1*mysinc(2*fc1*m)).*w;
    [H,fh] = freqz(h,1,16384,fs);  FIRresp{q} = 20*log10(abs(H)+eps);
    stopAtt = -20*log10(max(abs(H((fh>F_C2+TRANS_BW)|(fh<F_C1-TRANS_BW))))+eps);
    passRip = max(abs(20*log10(abs(H((fh>F_C1+TRANS_BW)&(fh<F_C2-TRANS_BW)))+eps)));
    D = N/2;  tmp = filter(h,1,[y1; zeros(D,1)]);  y2 = tmp(D+1:end);
    y3 = applySOS(sos,y2,IIR_ZERO_PHASE);  Y{q} = y3;
    [featAll(:,q),featNames] = featureVector(y3,fs);
    if ~isempty(ref), snrRef(q) = 10*log10(sum(ref.^2)/sum((ref-y3).^2)); end
    info(q,:) = [N+1, stopAtt, passRip, 0, 0];
    audiowrite(fullfile(outDir,['denoised_' name '.wav']),nrm(y3),fs);
end
[featIn,~] = featureVector(x,fs);

% 5. TABLES
fprintf('\n=========== FILTER COMPARISON ===========\n');
TF = table(info(:,1),info(:,2),info(:,3),'RowNames',winNames,'VariableNames',{'Taps','StopbandAtt_dB','PassbandRipple_dB'});
disp(TF);
TA = array2table([featIn featAll],'RowNames',featNames,'VariableNames',[{'Input_noisy'},winNames]);
fprintf('\n=========== SIGNAL FEATURES OF THE OUTPUT FOR EVERY WINDOW ===========\n');
disp(TA);
writetable(TF,fullfile(outDir,'filter_comparison.csv'),'WriteRowNames',true);
writetable(TA,fullfile(outDir,'features_all_windows.csv'),'WriteRowNames',true);
iSNR = find(strcmp(featNames,'Estimated SNR (dB)'));
[~,bestEst] = max(featAll(iSNR,:));
fprintf('Highest ESTIMATED SNR : %s (%.2f dB, input %.2f dB)\n',winNames{bestEst},featAll(iSNR,bestEst),featIn(iSNR));
if ~isempty(ref)
    TS = table(snrRef,'RowNames',winNames,'VariableNames',{'SNR_vs_clean_dB'});  disp(TS);
    [~,bestRef] = max(snrRef);
    fprintf('Highest TRUE SNR vs clean reference : %s (%.2f dB)\n',winNames{bestRef},snrRef(bestRef));
end

% 6. PLOTS
cols = lines(nW);  Lc = 101;  Nc = Lc-1;
newFig('A_All_Windows_Time');   hold on;
for q = 1:nW, plot(linspace(0,1,Lc),makeWin(winNames{q},Nc,P),'Color',cols(q,:),'LineWidth',1.3); end
grid on; legend(winNames,'Location','eastoutside'); xlabel('Normalised n/N'); ylabel('w[n]'); title('Window shapes (all with the same length)');

newFig('B_All_Windows_Spectra');   hold on;
for q = 1:nW
    wc = makeWin(winNames{q},Nc,P);  Wc = abs(fft(wc,2^15));  Wc = Wc(1:2^14+1)/max(Wc);
    plot((0:2^14)'/2^15,20*log10(Wc+eps),'Color',cols(q,:),'LineWidth',1.1);
end
grid on; ylim([-150 5]); legend(winNames,'Location','eastoutside'); xlabel('Frequency (cycles/sample)'); ylabel('dB'); title('Window spectra: main-lobe width vs side-lobe level');

newFig('C_FIR_Responses_Overlay');   hold on;
for q = 1:nW, plot(fh,FIRresp{q},'Color',cols(q,:),'LineWidth',1.1); end
grid on; ylim([-160 10]); xline(F_C1,'--k'); xline(F_C2,'--k'); legend([winNames,{'cut-offs'}],'Location','eastoutside');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('Band-pass FIR magnitude responses for every window');

newFig('D_Summary_Bars');
subplot(1,3,1); bar(info(:,1)); xticks(1:nW); xticklabels(winNames); xtickangle(45); grid on; title('FIR length (taps)');
subplot(1,3,2); bar(info(:,2)); xticks(1:nW); xticklabels(winNames); xtickangle(45); grid on; title('Stop-band attenuation (dB)');
subplot(1,3,3); bar(featAll(iSNR,:)); xticks(1:nW); xticklabels(winNames); xtickangle(45); grid on; title('Estimated output SNR (dB)');

newFig('E_Output_Spectra');   hold on;
for q = 1:nW, [pk_,fpk] = pwelch(Y{q},hamming(1024),512,1024,fs); plot(fpk,10*log10(pk_+eps),'Color',cols(q,:)); end
[pk_,fpk] = pwelch(x,hamming(1024),512,1024,fs); plot(fpk,10*log10(pk_+eps),'k--','LineWidth',1.3);
grid on; legend([winNames,{'Noisy input'}],'Location','eastoutside'); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)'); title('Output PSD for every window');

figs = findall(0,'Type','figure');
for k = 1:numel(figs)
    nm = regexprep(figs(k).Name,'[^\w]','_');  if isempty(nm), nm = sprintf('Figure_%d',figs(k).Number); end
    saveas(figs(k),fullfile(outDir,[nm '.png']));
end
fprintf('\nDONE. Everything saved in %s\n',outDir);

%  WINDOW HELPERS
function w = makeWin(name,N,P)
    switch name
        case 'Rectangular',    w = rectwin(N+1);
        case 'Bartlett',       w = bartlett(N+1);
        case 'Hanning',        w = hann(N+1);
        case 'Hamming',        w = hamming(N+1);
        case 'Blackman',       w = blackman(N+1);
        case 'BlackmanHarris', w = blackmanharris(N+1);
        case 'Nuttall',        w = nuttallwin(N+1);
        case 'FlatTop',        w = flattopwin(N+1);
        case 'Gaussian',       w = gausswin(N+1,P.gaussAlpha);
        case 'Tukey',          w = tukeywin(N+1,P.tukeyR);
        case 'Kaiser',         w = kaiser(N+1,kaiserBeta(P.kaiserAtt));
        case 'Chebyshev',      w = chebwin(N+1,P.chebAtt);
        otherwise, error('Unknown window');
    end
    w = w(:);
end

function b = kaiserBeta(A)
    if A > 50, b = 0.1102*(A-8.7);
    elseif A >= 21, b = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
    else, b = 0; end
end

function N = windowOrder(name,fs,TW,P)
    switch name
        case 'Rectangular',    C = 0.9;
        case 'Bartlett',       C = 3.05;
        case 'Hanning',        C = 3.1;
        case 'Hamming',        C = 3.3;
        case 'Blackman',       C = 5.5;
        case 'BlackmanHarris', C = 7.3;
        case 'Nuttall',        C = 7.3;
        case 'FlatTop',        C = 9.2;
        case 'Gaussian',       C = 4.0;
        case 'Tukey',          C = 2.0;
        case 'Kaiser',         N = ceil((P.kaiserAtt-7.95)/(2.285*2*pi*TW/fs)); return
        case 'Chebyshev',      N = ceil((P.chebAtt-7.95)/(2.285*2*pi*TW/fs));   return
    end
    N = ceil(C*fs/TW);
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
