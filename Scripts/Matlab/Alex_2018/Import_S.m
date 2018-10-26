del = ' ';
hlines = 5;
A = importdata('File1.s2p',del,hlines);
freq = A.data(:,1);
freqHz = freq.*1e9;