%% Lab 3: ECG
% ELEC3802 _ ECG3 _ Magnify

clc
clear
clf
%% Load data

load('ECG_noisy.mat');
% Set the position and dimesion of the polt
x0=10;
y0=10;
width=800;
height=600;

% Put all ECG signals into a matrix for easy access using a parameter
ECG=[ECG1_noisy;ECG2_noisy;ECG3_noisy;ECG4_noisy;ECG5_noisy];

% Loop for all magnified plots
for i = 1:5
    figure();
    set(gcf,'position',[x0,y0,width,height]);
    ECG_noisy = ECG(i,:);
    fc = 0.2;
    Wn = fc/(Fs/2);
    % High-pass filter processing
    [b,a] = butter(1,Wn,'high');
    ECG_low_filtered = filter(b,a,ECG_noisy);
    % Low-pass filter processing
    fc_9 = 35;
    Wn_9 = fc_9/(Fs/2);
    [b_9,a_9] = butter(9,Wn_9,'low');
    ECG_high_filtered = filter(b_9,a_9,ECG_low_filtered);

    % Notch filter processing
    noise_freq = 50;
    Wo = noise_freq/(Fs/2);  
    Bw = Wo/35;
    [b_n,a_n] = iirnotch(Wo,Bw);
    ECG_Power_line_filtered = filter(b_n,a_n,ECG_high_filtered);

    % Q9 Magnify and plot the filtered signal 
    N=length(ECG_Power_line_filtered);
    t=[0:N-1]/Fs;
    plot(t(500:1000),ECG_Power_line_filtered(500:1000));   
    title('ECG Segment');
    xlabel('Time (s)');
    ylabel('Intensity');
    title(['Filtered ECG',num2str(i),' Magnified Segment']);
    saveas(gcf,sprintf('Filtered_ECG%d_Magnify_segment.png',i));
end
% END of Lab 3: ECG - Magnify each signal for plots %