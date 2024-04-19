close all
clear 
clc

data = readmatrix('BozidarObradovic_vrat.csv');
fs = 250;

t = data(:, 1)/(4 * fs);
ecg = data(:, 2);
n = length(ecg);

%baseline filter
fc = 0.5;
[b, a] = butter(3, fc/(fs/2), 'high');
e1 = filtfilt(b, a, ecg);

%filtering between pulse frequencies (40 - 120 bpm)
f1 = 0.67;
f2 = 2;
[b,a] = butter(3, [f1, f2]/(fs/2), 'bandpass');
e3 = filtfilt(b, a, e1);

%normalisation (between -1 and 1)
e4 = e3;
e5  = 2 * (e4 - min(e4))/(max(e4) - min(e4)) - 1;

%define time and window interval
t1 = 10;
t2 = 80;
T = 10;
N = T*fs;
num = (t2-t1)/T;
avg_HR1 = zeros(1, num);
avg_HR2 = zeros(1, num);
avg_HR3 = zeros(1, num);

for ii = 1:num
    part = (t1*fs + (ii-1)*N : t1*fs + ii*N);
    t = (t1 + (ii-1)*T):1/fs:(t1 + ii*T);
    sig = e5(part);
    
    [~, locs] = findpeaks(sig, t, 'MinPeakHeight', 0.2, ...
        'MinPeakDistance',0.5);  
    val = (locs(end) - locs(1))/(length(locs)-1);
    avg_HR1(ii) = 60/val;
    
    y = autocorr(sig, 'NumLags', 1000);
    [~, lag] = findpeaks(y, 'MinPeakHeight', 0.2,...
        'MinPeakDistance', 100);
    val = (lag(end) - lag(1))/(length(lag)-1);
    avg_HR2(ii) = fs/val * 60;
    
    M = length(sig);
    z = abs(fft(sig));
    f = linspace(0, 1, fix(M/2)+1)*(fs/2); 
    vv = 1:length(f); 
    z = z(vv);    
    [~, ind] = max(z(2:end));
    fHR = f(ind+1);
    avg_HR3(ii) = 60*fHR;
    
end

HR1 = mean(avg_HR1);
s1 = std(avg_HR1);

HR2 = mean(avg_HR2);
s2 = std(avg_HR2);

HR3 = mean(avg_HR3);
s3 = std(avg_HR3);





