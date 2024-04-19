close all
clear
clc

file = load('MartaBosnjak_vrat.mat');
data = file.data;
fps = 60; fs = fps;
t1 = 10;
t2 = 80;

frames = (t1*fps:t2*fps); %frames between t1 i t2 seconds
images = data(:, :, frames); 
[m, n, k] = size(images);
temp = zeros(1, k);

%%

for ii = 1:k
    img = data(:, :, ii);
    img = 0.04 * img - 273.15;
    MEAN = mean(img, 'all');
    img = img - MEAN;
    MAX1 = max(img, [], 'all');
    MIN1 = min(img, [], 'all');
    img = (img - MIN1)/(MAX1-MIN1);
    
%     img1 = img > 0.7;
%     border_c = sum(img1, 1);
%     border_r = sum(img1, 2);
%     c1 = 150;
%     [~, c2] = max(border_c);
%     r1 = find(border_r > 140, 1, 'last');
%     r2 = find(border_r > 85, 1, 'last');
    r1 = 180; r2 = 248;
    c1 = 135; c2 = 177;
    ROI = {r1:r2, c1:c2};

    neck = img(ROI{1}, ROI{2});
    J = histeq(neck);
    J1 = J > 0.7;
    SE = strel('disk', 1, 6);
    J1 = imclose(J1, SE);
    J2 = medfilt2(J1, [5 4]);
    J3 = bwskel(J2);
    J4 = imdilate(J3, SE);
    
    size = sum(J4, 'all');
    rr = J4 .* neck;
    
    val = sum(rr, 'all') / size;
    temp(ii) = MIN1 + val*(MAX1 - MIN1) + MEAN;
end

% figure()
% plot(temp)

%filtering between pulse frequencies (40 - 120 bpm)
f1 = 0.67;
f2 = 2;
[b,a] = butter(3, [f1, f2]/(fs/2), 'bandpass');
s1 = filtfilt(b, a, temp);

% figure()
% plot(s1)

%define time and window interval
T = 10;
N = T*fs;
num = (t2-t1)/T;
t = T:1/fs:2*T;
avg_HR1 = zeros(1, num);
avg_HR2 = zeros(1, num);
avg_HR3 = zeros(1, num);

for ii = 1:num
    part = (1+(ii-1)*N:ii*N+1);
    sig = s1(part);
    
    [~, locs] = findpeaks(sig, t, 'MinPeakHeight', 0.005, ...
         'MinPeakDistance', 0.5);  
    val = (locs(end) - locs(1))/(length(locs)-1);
    avg_HR1(ii) = 60/val;
    
    y = autocorr(sig, 'NumLags', 500);
    [~, lag] = findpeaks(y, 'MinPeakHeight', 0.1,...
        'MinPeakDistance', 25); 
    val = (lag(end) - lag(1))/(length(lag)-1);
    avg_HR2(ii) = fs/val * 60;
    
    M = length(sig);
    z = abs(fft(sig));
    f = linspace(0, 1, fix(M/2)+1)*(fs/2); 
    vv = 1:length(f); 
    z = z(vv);
    [~, ind] = max(z);
    fHR = f(ind);
    avg_HR3(ii) = 60*fHR;
    
end

HR1 = mean(avg_HR1);
s1 = std(avg_HR1);

HR2 = mean(avg_HR2);
s2 = std(avg_HR2);

HR3 = mean(avg_HR3);
s3 = std(avg_HR3);




