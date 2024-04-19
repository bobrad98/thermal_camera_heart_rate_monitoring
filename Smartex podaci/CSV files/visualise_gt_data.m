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

%define time interval
t1 = 10;
t2 = 20;
part = (t1*fs:t2*fs);
ecg = ecg(part);
t11 = t(part);
e5 = e5(part);

figure()
plot(t11, e3(part))
xlabel('Vreme [s]');
ylabel('Amplituda [a.u.]')
grid on
% exportgraphics(gcf,'fig_filtEKG.png','Resolution',300)

%find peak, avgHR and stdHR
%1
[pks, locs] = findpeaks(e5, 'MinPeakHeight', 0.2, 'MinPeakDistance',0.5);

figure()
plot(t11, e5);
hold on
scatter(t11(locs), pks, 'filled');
grid on
xlabel('Vreme [s]');
ylabel('Amplituda [a.u.]');
% exportgraphics(gcf,'fig_1EKG.png','Resolution',300)

val = (locs(end) - locs(1))/(length(locs)-1);
avg_HR1 = 60/val;

%2
y = autocorr(e5, 'NumLags', 1000);
%autocorr(e5, 'NumLags', 1000)
[pks, lag] = findpeaks(y, 'MinPeakHeight', 0.2, 'MinPeakDistance', 100);
figure()
plot(1:1001, y);
hold on
scatter(lag, pks, 'filled');
grid on
xlabel('Kašnjenje');
ylabel('Autokorelacija');
%exportgraphics(gcf,'fig_2EKG.png','Resolution',300)

val = (lag(end) - lag(1))/(length(lag)-1);
avg_HR2 = fs/val * 60;

%3
N = length(e5);

z = abs(fft(e5));
f = linspace(0, 1, fix(N/2)+1)*(fs/2); 
vv = 1:length(f); 
z = z(vv);

[~, ind] = max(z);
 
fHR = f(ind);
avg_HR3 = 60*fHR;

figure()
plot(f, z);
xlabel('Učestanost [Hz]');
ylabel('Magnituda');
grid on
xlim([0 10])
% exportgraphics(gcf,'fig_3EKG.png','Resolution',300)