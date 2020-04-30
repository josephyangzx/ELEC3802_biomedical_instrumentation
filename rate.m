function bpm = rate(x, Fs, T)
heart_beat = thresh(x, T) ; % find the heart beat (step 10)
% sum to find the number of heart beat occur during the
% segment (since each heart beat is 1)
number = sum(heart_beat);
% from the sampling frequency, compute the length of the
% segment in minute
minute_length = length(x)/(Fs*60);
% compute the beats per minute
bpm = number/minute_length;