%% Lab 3: ECG
% ELEC3802 _ ECG3 _ Magnify

clc
clear
clf
%% Load data

load('ECG_noisy.mat')
%% A. Examine the frequency spectrum and design filters for ECG
%Examine the DC off-set

% Pick ECG1_noisy

%% Q5 Low order filter
ECG=[ECG1_noisy;ECG2_noisy;ECG3_noisy;ECG4_noisy;ECG5_noisy];
% The variable for which ECG signal
x = 5;
ECG3_noisy = ECG(x,:);
fc = 0.2;
Wn = fc/(Fs/2);
% 1st order 
[b,a] = butter(1,Wn,'high');

ECG3_low_filtered = filter(b,a,ECG3_noisy);

%% Q6 9th order filter

fc_9 = 35;
Wn_9 = fc_9/(Fs/2);
[b_9,a_9] = butter(9,Wn_9,'low');
ECG3_high_filtered = filter(b_9,a_9,ECG3_low_filtered);

%% Q7 Notch filter

noise_freq = 50;
Wo = noise_freq/(Fs/2);  
Bw = Wo/35;
[b_n,a_n] = iirnotch(Wo,Bw);

% Apply the filter to the signal
ECG3_Power_line_filtered = filter(b_n,a_n,ECG3_high_filtered);

%% 
% Q9 Magnify
x0=10;
y0=10;
width=1000;
height=2000;
figure(1)
set(gcf,'position',[x0,y0,width,height]);
N=length(ECG3_Power_line_filtered);
t=[0:N-1]/Fs;
plot(t(500:1000),ECG3_Power_line_filtered(500:1000));   
title('ECG Segment');
xlabel('Time (s)');
ylabel('Intensity');
title(['Filtered ECG',num2str(x),' Magnified Segment']);
saveas(gcf,sprintf('Filtered_ECG%d_Magnify_segment.png',x));