%% Lab 3: ECG
% ELEC3802
% This document loops through each ECG signal

clc
clear
clf
%% Load data
load('ECG_noisy.mat')

%% Repeat for 5-9
% Combine the matrix into one for looping
ECG=[ECG1_noisy;ECG2_noisy;ECG3_noisy;ECG4_noisy;ECG5_noisy];
[r,c] = size(ECG); % Get the size of the matrix

% Outer loop through each ECG measurement
for i = 1:r
    figure()
    % Inpect the noisy ECG signal
%     periodogram(ECG(i,1:c))
    % Save the file
%     title(['ECG',num2str(i),'_',' noisy']);
%     saveas(gcf,sprintf('ECG%d_noisy.png',i));
    
    % Q5 
    % Low order filter - High Pass Filter
    fc = 0.2;
    Wn = fc/(Fs/2);
    % 1st order 
    [b,a] = butter(1,Wn,'high');

    % Check the frequency response of the filter.
%     figure();
%     freqz(b,a);
%     title(['Frequency Response ',num2str(i)]);
%     saveas(gcf,sprintf('frequencyResponse_%d.png',i));

    % Visualize filter
    ECG_low_filtered = filter(b,a,ECG(i,1:c));
%     figure();
%     periodogram(ECG_low_filtered);
%     title(['ECG High Pass Filter ',num2str(i)]);
%     saveas(gcf,sprintf('ECG_High_Pass_Filter_%d.png',i));
%     freqz(ECG_low_filtered);
%     title(['High-Pass Filter Frequency Response of ECG',num2str(i)]);
%     saveas(gcf,sprintf('High_Pass_Filtered_Frequency_Response_ECG%d.png',i));
    
    % Q6 Low-pass Filter
    fc_9 = 35;
    Wn_9 = fc_9/(Fs/2);
    [b_9,a_9] = butter(9,Wn_9,'low');
    ECG_high_filtered = filter(b_9,a_9,ECG_low_filtered);
%     figure();
%     periodogram(ECG_high_filtered);
%     title(['ECG Low Pass Filter ',num2str(i)]);
%     saveas(gcf,sprintf('ECG_Low_Pass_Filter%d.png',i));
%     
    % Q7 Power Spectrum Density
    % Design a filter with a Q-factor of Q=35 to remove a 50 Hz tone from 
    % system running at 250 Hz.
    noise_freq = 50;
    Wo = noise_freq/(Fs/2);  
    Bw = Wo/35;
    [b_n,a_n] = iirnotch(Wo,Bw);
%     fvtool(b_n,a_n); % Visualize the filter
    
    % Apply the filter to the signal
    ECG_Power_Spectrum_Density = filter(b_n,a_n,ECG_high_filtered);
%     figure();
%     periodogram(ECG_Power_Spectrum_Density);
%     title(['ECG Notch Filter Filter ',num2str(i)]);
%     saveas(gcf,sprintf('ECG_Power_Spectrum_Density_%d.png',i));
    
    % Examine filtered signal
%     figure();
%     plot(t(1:500),ECG1_Power_line_filtered(1:500));   
%     title('ECG Filtered Signal');
%     xlabel('Time (s)');
%     ylabel('Intensity');
%     saveas(gcf,sprintf('ECG_Filtered_Exmained_%d.png',i));
    
    % Q9 Examine the signal in 10 segments

    % Set up the loop
    seg= 10;
    segmentMatix = [];
    N = length(ECG_Power_Spectrum_Density);
    t = N/seg;
    % Set up the size of the plot
    x0=10;
    y0=10;
    width=10000;
    height=2000;
%     figure()
%     set(gcf,'position',[x0,y0,width,height])
    
    % Inner loop to obatin the segments
    for a = 1:seg
        segmentMatix(a,:) = ECG_Power_Spectrum_Density((a-1)*t+1 : t*a );
%         subplot(10,1,a);
%         plot(segmentMatix(a,:));
    end 
%     sgtitle(['Filtered ECG',num2str(i),' Noisy in 10 Equal Segments']);
%     saveas(gcf,sprintf('Filtered_ECG%d_noisy_10-equal_segments.png',i));
    
    %% BPM
    % Calculate the BPM on Segment 2 using T = 3000
    T = 1000;
    bpm = rate(segmentMatix(2,:), Fs, T);
    disp('Caculated BPM from given formula: ');
    disp(round(bpm));

    % Correction on BPM formula
    x = segmentMatix(2,:);
    heart_beat = thresh(x,T);
    % Find the heart beat when the consecutive point is 0
    number = 0;
    for i= 1: length(x)
        if (heart_beat(i) == 1) && (heart_beat(i+1) == 0)
            number = number + 1;
        end
    end
    % from the sampling frequency, compute the length of the
    % segment in minute
    minute_length = length(x)/(Fs*60);
    % compute the beats per minute
    bpmCorrectm = number/minute_length;
    disp('Caculated BPM from corrected formula: ');
    disp(round(bpmCorrectm));

end

% END of Lab 3: ECG - Looping all ECG signals %