% EMG Cortical Silent Period calculations. See Garvey, M.A., Clinical Neurophys 2001
% and Wilke, S., PLOS one 2016

% This uses the average of the 10 MEPs exported from BrainVision Analyzer,
% rectifies the signal, and determines the MEP onset, CSP onset, CSP
% offset, and the CSP-M and CSP durations.
clear all
close all

% Select the subject directory
str = computer;
if strcmp(str,'MACI64') == 1
    directory = uigetdir('/Volumes/mnl/Data/Adaptation/SICI_biman1/');
    % Set Currect Directory to subject dir
    cd([directory '/Peaks'])
    subID = directory(46:end);
    groupID = directory(42:44);
else
    directory = uigetdir('Z:\Data\Adaptation\SICI_biman1\');
    % Set Currect Directory to subject dir
    cd([directory '\Peaks'])
    subID = directory(36:end);
    groupID = directory(32:34);
end

dir_list = dir('*.dat'); % Get only the .dat files
dir_list = {dir_list.name};

exportData = NaN(4,8);

for i = 1:length(dir_list)
data = load(dir_list{i});
filename = dir_list{i};
%Find hemisphere
temp = find(filename == '_');
temp2 = find(filename == '-');
hem = filename((temp+1):(temp+4));
% Find pre/post
pre_post = filename(temp(2)+1:temp(2)+4);
% Find test (1mv, sici, 120, 120_sici)
test = filename(temp(3)+1:temp2-1);
clear temp temp2
if strcmp(hem, 'lhem') == 1
    channel = 1;
    data = data(1,2:end);
    hem = 1; % Set left hemisphere data to "1" (used in R code)
else
    channel = 2;
    data = data(2,2:end);
    hem = 2; % Set right hemisphere data to "2"
end

if strcmp(pre_post, 'pre_') == 1
    pre_post = 1; % Set pre vs. post to "1" (used in R)
else
    pre_post = 2; % Set pre vs. post to "2"
end

sampleRate = 2500; % in Hz
time = (1:length(data))/sampleRate * 1000; % Convert to ms

% Calculate difference vector for Pre-stimulus interval
diff = zeros(1,100); % Pre-stimulus invterval initiation = 100
diff(1) = NaN;
for i = 2:100
    diff(i) = data(i) - data(i-1);
end
diffABS = abs(diff);
mcd = nanmean(diffABS); % This is the "Mean Consecutive Difference"
meanPreStim = nanmean(data(1:100));

% Change for the number of SDs from mean
% 3 SD = 2.66
% 2.5 SD = 2.22
% 2 SD = 1.77
% 1 SD = 0.89
critZ = 0.89;
lvl = meanPreStim - mcd*critZ; % Lower Variation Limit (LVL)

%% Onset Calculation
% Find indicies of EMG data that fall below the lower variation limit
% starting 10ms after the stimulus onset (ie: 110ms = 0.110 seconds)
indSubLVL = find(data((sampleRate * 0.110):end) < lvl) + (sampleRate * 0.110);

%Find 5 consecutive data points under the LVL and determine the index in 'data' where this occurs
for i = 1:(length(indSubLVL) - 4)
   if ~isempty(find(sum((indSubLVL(i):(indSubLVL(i) + 4)) - [indSubLVL(i) indSubLVL(i+1) indSubLVL(i+2) indSubLVL(i+3) indSubLVL(i+4)]) == 0, 1))
       indOnset = indSubLVL(i)-1; % Not sure why it is -1 here, but on visual inspection this is correct!!! I suspect it is because I removed one datum at the beginning? check this again
       break
   end
end % This is so damn simple...best code ever

onsetTime = indOnset/sampleRate * 1000;

%% Offset Calculation
indSuperLVL = find(data(indOnset:end) > lvl) + (indOnset - 1); % Some of the early ones might be right after indOnset. These are clearly not the Offset
% Need to find the 5ms windows following each of these indicies and
% determine the first window to have >=50% data points above the LVL.

% I use 13 because .005(s) * 2500(Hz) = 12.5 samples, so I round up
for i = 1:(length(indSuperLVL) - 13)
    windows(i,:) = [data(indSuperLVL(i)) data(indSuperLVL(i)+1) data(indSuperLVL(i)+2)...
        data(indSuperLVL(i)+3) data(indSuperLVL(i)+4) data(indSuperLVL(i)+5)...
        data(indSuperLVL(i)+6) data(indSuperLVL(i)+7) data(indSuperLVL(i)+8)...
        data(indSuperLVL(i)+9) data(indSuperLVL(i)+10) data(indSuperLVL(i)+11)...
        data(indSuperLVL(i)+12)];
    
    numSuper = length(find(windows(i,:) > lvl));
    if numSuper > 6
        indOffset = indSuperLVL(i);
        break
    end
end % A little more messy, but you can't sequence using indicies and preserve spacing

offsetTime = indOffset/sampleRate * 1000;

cspDuration = offsetTime - onsetTime;

%% MEP Onset Calculation from Wilke
[MEPmax indMEPmax] = max(data);
indMEPonset = indMEPmax - (find(data(indMEPmax:-1:1) < meanPreStim, 1)) + 2;
% The +2 comes from the index starting at indMEPmax - 1 and the search
% ending at "the point we want" - 1. Therefore, the + 2 corrects those.
MEPonsetTime = indMEPonset/sampleRate * 1000;

cspMEPduration = offsetTime - MEPonsetTime;

%% Plotting

flag = 2;
while flag > 1
    figure('Position', [100 100 1920/2 1080/2]); %[bottom left corner coords X and Y, W, H]
    plot(time, data)
    hold on
    line([0 max(time)], [meanPreStim meanPreStim], 'Color','k')
    hold on
    line([0 max(time)], [(meanPreStim - mcd*critZ) (meanPreStim - mcd*critZ)], 'Color','r')
    hold on
    line([0 max(time)], [(meanPreStim + mcd*critZ) (meanPreStim + mcd*critZ)], 'Color','r')
    hold on
    plot(onsetTime, data(indOnset), 'go')
    hold on
    plot(offsetTime, data(indOffset), 'mo')
    hold on
    plot(MEPonsetTime, data(indMEPonset), 'co')
    ylim([-50 200])
    title(['Pre/Post = ' num2str(pre_post), ' ','Hemisphere = ' num2str(hem), ' ', 'Test = ' test])
    
    % Verify Onset/Offset
    button = questdlg('Confirm movement onset/offset?','Onset and offset markers:','Yes','No','Reject','No');
    if strcmp(button,'Yes')
        cspDuration = offsetTime - onsetTime;
        cspMEPduration = offsetTime - MEPonsetTime;
        % Populate the exportData array
        if pre_post == 1 && hem == 1 && strcmp(test,'1mV')==1
            exportData(1,1) = cspDuration;
            exportData(1,5) = cspMEPduration;
        elseif pre_post == 2 && hem == 1 && strcmp(test,'1mV')==1
            exportData(2,1) = cspDuration;
            exportData(2,5) = cspMEPduration;
        elseif pre_post == 1 && hem == 2 && strcmp(test,'1mV')==1
            exportData(3,1) = cspDuration;
            exportData(3,5) = cspMEPduration;
        elseif pre_post == 2 && hem == 2 && strcmp(test,'1mV')==1
            exportData(4,1) = cspDuration;
            exportData(4,5) = cspMEPduration;
            
        elseif pre_post == 1 && hem == 1 && strcmp(test,'sici')==1
            exportData(1,2) = cspDuration;
            exportData(1,6) = cspMEPduration;
        elseif pre_post == 2 && hem == 1 && strcmp(test,'sici')==1
            exportData(2,2) = cspDuration;
            exportData(2,6) = cspMEPduration;
        elseif pre_post == 1 && hem == 2 && strcmp(test,'sici')==1
            exportData(3,2) = cspDuration;
            exportData(3,6) = cspMEPduration;
        elseif pre_post == 2 && hem == 2 && strcmp(test,'sici')==1
            exportData(4,2) = cspDuration;
            exportData(4,6) = cspMEPduration;
            
        elseif pre_post == 1 && hem == 1 && strcmp(test,'120')==1
            exportData(1,3) = cspDuration;
            exportData(1,7) = cspMEPduration;
        elseif pre_post == 2 && hem == 1 && strcmp(test,'120')==1
            exportData(2,3) = cspDuration;
            exportData(2,7) = cspMEPduration;
        elseif pre_post == 1 && hem == 2 && strcmp(test,'120')==1
            exportData(3,3) = cspDuration;
            exportData(3,7) = cspMEPduration;
        elseif pre_post == 2 && hem == 2 && strcmp(test,'120')==1
            exportData(4,3) = cspDuration;
            exportData(4,7) = cspMEPduration;
        
        elseif pre_post == 1 && hem == 1 && strcmp(test,'120_sici')==1
            exportData(1,4) = cspDuration;
            exportData(1,8) = cspMEPduration;
        elseif pre_post == 2 && hem == 1 && strcmp(test,'120_sici')==1
            exportData(2,4) = cspDuration;
            exportData(2,8) = cspMEPduration;
        elseif pre_post == 1 && hem == 2 && strcmp(test,'120_sici')==1
            exportData(3,4) = cspDuration;
            exportData(3,8) = cspMEPduration;
        elseif pre_post == 2 && hem == 2 && strcmp(test,'120_sici')==1
            exportData(4,4) = cspDuration;
            exportData(4,8) = cspMEPduration;
        end
        close(1);
        flag=1;
    elseif strcmp(button,'No') %user can verify if movement onset was computed correctly
        [loc_onset,loc_size]=ginput(3);
        indMEPonset = round(loc_onset(1)*sampleRate/1000); % replace MEPonset with user defined input
        indOnset = round(loc_onset(2)*sampleRate/1000); % replace CSPonset with use defined input
        indOffset = round(loc_onset(3)*sampleRate/1000); % replace CSPoffset with use defined input
        
        MEPonsetTime = indMEPonset/sampleRate*1000;
        onsetTime = indOnset/sampleRate*1000;
        offsetTime = indOffset/sampleRate*1000;
        
        cspDuration = offsetTime - onsetTime;
        cspMEPduration = offsetTime - MEPonsetTime;
        
        close(1);
        flag=2;
    elseif strcmp(button,'Reject')
        indMEPonset = NaN;
        indOnset = NaN;
        indOffset = NaN;
        cspDuration = NaN;
        cspMEPduration = NaN;
        %wrong_trial(i,1)=1;			% Keeps track of rejected trials
        close(1);
        flag=1;
    end
end

end

dlmwrite([subID ,'_csp.txt'],exportData);
cd('Z:\Data\Adaptation\SICI_biman1')