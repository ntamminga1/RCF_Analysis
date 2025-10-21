function [E,dEdx]=ExtractSRIM(SRIMstr)
% This funciton extracts the E and dE/dx values from the SRIM files
% previously generated.
dir = 'C:\Users\ntamm\OneDrive\Desktop\SRIM-2013\SRIM Outputs';% SRIM Output dir location

text = readlines(strcat(dir,'\',SRIMstr,'.txt')); %Read file text
start_line = '  --------------  ---------- ---------- ----------  ----------  ----------';
end_line = '-----------------------------------------------------------';
startindex = find(text==start_line,1)+1; % Start index of interesting data
endindex = find(text==end_line,1); % End index of interesting data +1
% Generate empty arrays
E = zeros([endindex-startindex 1]);
dEdxE = zeros([endindex-startindex 1]);
dEdxN = zeros([endindex-startindex 1]);

for i=1:length(E)
    str = text(i-1+startindex);
    if extractBetween(str,9,11) == 'MeV'
        E(i) = str2double(extractBefore(str,8));
    elseif extractBetween(str,9,11) == 'keV'
        E(i) = str2double(extractBefore(str,8))/1000;
    end
    dEdxE(i) = str2double(extractBetween(str,15,25));
    dEdxN(i) = str2double(extractBetween(str,26,36));
end
dEdx = dEdxN + dEdxE;
end