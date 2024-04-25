function dat = load_cmm_CUDA(base,ending)

if nargin < 2 || isempty(ending)
    ending = ".data";
else
    ending = ending +".data";
end


file=base +"/Monitoring_data/Mesure/Epot"+ending;
fileID = fopen(file,'r');
dat.Epot =fread(fileID,'double');
fclose(fileID);
file=base +"/Monitoring_data/Mesure/Ekin"+ending;
fileID = fopen(file,'r');
dat.Ekin =fread(fileID,'double');
fclose(fileID);
file=base +"/Monitoring_data/Mesure/Mass"+ending;
fileID = fopen(file,'r');
dat.mass =fread(fileID,'double');
fclose(fileID);
file=base +"/Monitoring_data/Mesure/Etot"+ending;
fileID = fopen(file,'r');
dat.Etot =fread(fileID,'double');
fclose(fileID);

file=base+"/Monitoring_data/Mesure/Time_s"+ending;
fileID = fopen(file,'r');
dat.time =fread(fileID,'double');
fclose(fileID);


[dat.Coarse,dat.Fine,dat.caseName,dat.T] = extract_data_from_filename(base)



end


function [Coarse,Fine,caseName,T]=extract_data_from_filename(base)

gridPattern = 'C(\d+)_F(\d+)';
gridSize = regexp(base, gridPattern, 'tokens');
Coarse = str2double(gridSize{1}{1});
Fine = str2double(gridSize{1}{2});

% Extract the case name
casePattern = 'vp_(.*?)_C';
caseName = regexp(base, casePattern, 'tokens');
caseName = caseName{1}{1};

% Extract the total time (T)
timePattern = 'T(\d+)';
totalTime = regexp(base, timePattern, 'tokens');
T = str2double(totalTime{1}{1});
end