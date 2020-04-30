function s = thresh(x, T)
% x: input ECG signal
% T: the threshold
% compute the vector s:
% s(i) = 0 if the x(i) < T or = T, s(i) = 1 if x(i) > T
s = x > T;