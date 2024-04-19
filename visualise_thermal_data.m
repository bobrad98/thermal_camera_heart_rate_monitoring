close all
clear
clc

file = load('MartaBosnjak_vrat.mat');
data = file.data;
fps = 60; fs = fps;
t1 = 10; t2 = 20;

frames = (t1*fps:t2*fps); %frames between t1 i t2 seconds
images = data(:, :, frames); 
[m, n, k] = size(images);
temp = zeros(1, k);

%%
ii = 1;

img = images(:, :, ii);
img = 0.04 * img - 273.15;

MEAN = mean(img, 'all');

figure()
imshow(img, []);

img = img - mean(img, 'all');

figure()
imshow(img, [])

MAX1 = max(img, [], 'all');
MIN1 = min(img, [], 'all');
img = (img - MIN1)/(MAX1-MIN1);

figure(5)
imshow(img, [])
% exportgraphics(gcf,'fig_step1.png','Resolution',300)

figure()
imshow(img1, []);
% exportgraphics(gcf,'fig_bin.png','Resolution',300) 
% 
% border_c = sum(img1, 1);
% border_r = sum(img1, 2);
% c1 = 150;
% [~, c2] = max(border_c);
% r1 = find(border_r > 140, 1, 'last');
% r2 = find(border_r > 75, 1, 'last');
r1 = 180; r2 = 248;
c1 = 135; c2 = 177;
ROI = {r1:r2, c1:c2};

figure(7)
imshow(img, [])
hold on
rectangle('Position', [c1, r1, c2-c1, r2-r1])
% exportgraphics(gcf,'fig_rect.png','Resolution',300) 

neck = img(ROI{1}, ROI{2});
figure(8)
imshow(neck, [])
% exportgraphics(gcf,'fig_neck.png','Resolution',300)

figure(9)
histogram(neck, 'BinMethod', 'fd')
c = colorbar;
c.Ticks = [];
c.Location = 'southoutside';
colormap(gray(256))
% exportgraphics(gcf,'fig_h1.png','Resolution',300) 

J = histeq(neck);

figure()
imshow(J, [])
% exportgraphics(gcf,'fig_eqROI.png','Resolution',300)

figure()
histogram(J, 'BinMethod', 'fd')
c = colorbar;
c.Ticks = [];
c.Location = 'southoutside';
colormap(gray(256))
% exportgraphics(gcf,'fig_heq.png','Resolution',300)

J1 = J > 0.8;
SE = strel('disk', 1, 6);
J2 = imclose(J1, SE);
J3 = medfilt2(J2, [5 4]);
J4 = bwskel(J3);
J5 = imdilate(J4, SE);

figure()
subplot(1,5,1)
imshow(J1, [])
title('Binarizacija', 'FontSize', 7);
subplot(1,5,2)
imshow(J2,[])
title('Zatvaranje', 'FontSize', 7);
subplot(1,5,3)
imshow(J3, [])
title('Median filter', 'FontSize', 7)
subplot(1,5,4)
imshow(J4,[])
title('Skeletonizacija', 'FontSize', 7)
subplot(1,5,5)
imshow(J5,[])
title('Dilatacija', 'FontSize', 7)
%exportgraphics(gcf,'fig_op.png','Resolution',300)


%%
% close all
% 
% temp1 = temp;
% t = 10:1/fs:20;
% 
% f1 = 0.67;
% f2 = 2;
% [b,a] = butter(3, [f1, f2]/(fs/2), 'bandpass');
% s1 = filtfilt(b, a, temp1);
% 
% figure()
% plot(t, s1)
% xlabel('Vreme [s]');
% ylabel('Amplituda [a.u.]')
% grid on
% % exportgraphics(gcf,'fig_filtT.png','Resolution',300)
% 
% %1
% [pks, locs] = findpeaks(s1, t, 'MinPeakHeight', 0.005,...
%     'MinPeakDistance', 0.5);
% 
% figure()
% plot(t, s1);
% hold on
% scatter(locs, pks, 'filled');
% grid on
% xlabel('Vreme [s]');
% ylabel('Amplituda [a.u.]');
% % exportgraphics(gcf,'fig_1T.png','Resolution',300)
% 
% % ds = length(pks);
% % avg_HR = ds/10 * 60;
% 
% sed = (locs(end) - locs(1))/(length(locs)-1);
% erer = 60/sed;
% 
% %2
% y = autocorr(s1, 'NumLags', 500);
% [pks, lag] = findpeaks(y, 'MinPeakHeight', 0.1, 'MinPeakDistance', 25);
% 
% figure()
% plot(1:501, y);
% hold on
% scatter(lag, pks, 'filled');
% grid on
% xlabel('Kašnjenje');
% ylabel('Autokorelacija');
% % exportgraphics(gcf,'fig_2T.png','Resolution',300)
% 
% sd = (lag(end) - lag(1))/(length(lag)-1);
% hrhr = fs/sd * 60;
% 
% %3
% N = length(s1);
% 
% z = abs(fft(s1));
% f = linspace(0, 1, fix(N/2)+1)*(fs/2); 
% vv = 1:length(f); 
% z = z(vv);
% [~, ind] = max(z);
% fHR = f(ind);
% avg_HR3 = 60*fHR;
% 
% figure()
% plot(f, z);
% xlabel('Učestanost [Hz]');
% ylabel('Magnituda');
% grid on
% xlim([0 10])
% % exportgraphics(gcf,'fig_3T.png','Resolution',300)