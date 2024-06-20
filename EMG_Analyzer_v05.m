% This reads in the EMG data of every trial (10 per condition) for all
% subjects and compiles the data into a master list. The output is to be
% used with SPSS

clear all
close all

% Select the subject directory
str = computer;
if strcmp(str,'MACI64') == 1
    cd('/Volumes/mnl/Data/Adaptation/SICI_biman1/');
% else
%     cd('Z:\Data\Adaptation\SICI_biman1\');
%     sub = dir('G01\test*');
%     sub = struct2cell(sub);
%     sub = sub(1,:);
%     numSub = length(sub);
%         subNum = zeros(numSub,1);
% 
%         for i = 1:numSub
%             buff1 = sub{i};
%             subNum(i) = str2double(buff1(5:end)); % Selects just the subject number
%         end
% 
%         for j = 1:numSub % For each subject, set current directory to the "Peaks" folder to find .txt data files
%          
%             cd(['G01\',sub{j},'\Peaks\']);
%             files = dir('*.txt');
%             buff2 = struct2cell(files);
%             if size(buff2,2) > 4
%                 filenames(j,:) = buff2(1,:);
%             end
%             cd('Z:\Data\Adaptation\SICI_biman1\');
%         end
% end

else
    directory = uigetdir('Z:\Data\Adaptation\SICI_biman1\');
    % Get Subject and Group IDs
    subjectID = directory(40:end);
    subjectID = str2double(subjectID);
    groupID = directory(33:34);
    groupID = str2double(groupID);
    % Set Currect Directory to subject dir
    cd([directory '\Peaks'])
end

% dialog box will appear asking for user to input coordination (control/iso/mirror), age,
% gender, LQ. These can be found in SICI_coding excel file, and handedness forms
prompt = {'Enter Coordination (c = 1, i = 2, m = 3):','Enter age:','Enter gender (F = 1, M = 2:','Enter laterality quotient:'};
dlg_title = 'Input demographics';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);
answer = cellfun(@(answer)str2double(answer),answer);



% Find all text files
files = dir('*.txt');

% find the subject ID so we can delete existing data output files (this jacks up the Main loop)
sID_end = strfind(files(1).name == '_',1);
sID_end = sID_end(1) - 1;
sID = files(1).name(1:sID_end);

% DELETE peak data and peak data outlier corrected
delete([sID '_peak_data.txt'],[sID '_peak_data_oc.txt'])

% Find the text files again (with previous data files removed)
files = dir('*.txt');

output = zeros(10,38);
% Main loop to generate output file
for i = 1:size(files,1)
    % Read in data
    [Filename,Segment,Min1L,Max1L,peak2peak,area] = importfile(files(i).name);

    condition_start = strfind(Filename{1},'_');
    condition_stop = strfind(Filename{1},'-');
    condition = Filename{1}(condition_start(1)+1:condition_stop-1); % gets the condition name (ie: 'rhem_pre_sici')
    

    % Fill output array with calculated variables
    % This switch statement fills in the rest of output with the desired condition variables
    switch condition
        case 'lhem_pre_1mV'
            output(1:max(Segment),7) = abs(peak2peak);
            output(1:max(Segment),23) = area;
        case 'lhem_pre_sici'
            output(1:max(Segment),8) = abs(peak2peak);
            output(1:max(Segment),24) = area;
        case 'lhem_pre_120'
            output(1:max(Segment),9) = abs(peak2peak);
            output(1:max(Segment),25) = area;
        case 'lhem_pre_120_sici'
            output(1:max(Segment),10) = abs(peak2peak);
            output(1:max(Segment),26) = area;
        case 'lhem_post_1mV'
            output(1:max(Segment),11) = abs(peak2peak);
            output(1:max(Segment),27) = area;
        case 'lhem_post_sici'
            output(1:max(Segment),12) = abs(peak2peak);
            output(1:max(Segment),28) = area;
        case 'lhem_post_120'
            output(1:max(Segment),13) = abs(peak2peak);
            output(1:max(Segment),29) = area;
        case 'lhem_post_120_sici'
            output(1:max(Segment),14) = abs(peak2peak);
            output(1:max(Segment),30) = area;
        case 'rhem_pre_1mV'
            output(1:max(Segment),15) = abs(peak2peak);
            output(1:max(Segment),31) = area;
        case 'rhem_pre_sici'
            output(1:max(Segment),16) = abs(peak2peak);
            output(1:max(Segment),32) = area;
        case 'rhem_pre_120'
            output(1:max(Segment),17) = abs(peak2peak);
            output(1:max(Segment),33) = area;
        case 'rhem_pre_120_sici'
            output(1:max(Segment),18) = abs(peak2peak);
            output(1:max(Segment),34) = area;
        case 'rhem_post_1mV'
            output(1:max(Segment),19) = abs(peak2peak);
            output(1:max(Segment),35) = area;
        case 'rhem_post_sici'
            output(1:max(Segment),20) = abs(peak2peak);
            output(1:max(Segment),36) = area;
        case 'rhem_post_120'
            output(1:max(Segment),21) = abs(peak2peak);
            output(1:max(Segment),37) = area;
        case 'rhem_post_120_sici'
            output(1:max(Segment),22) = abs(peak2peak);
            output(1:max(Segment),38) = area;            
        otherwise
            msgbox('Error: One or more variable was not written to output. Please check .txt file names')
    end

%     % Fill in demographic variables with the user prompted inputs.
%     
%     group = cell(numTrials,1);
%     group(:) = {groupID};
%     coord = cell(numTrials,1);
%     coord(:) = {answer(1)};
%     subID = cell(numTrials,1);
%     subID(:) = {subjectID};
%     age = cell(numTrials,1);
%     age(:) = num2cell(answer(2));
%     gender = cell(numTrials,1);
%     gender(:) = {answer(3)};
%     LQ = cell(numTrials,1);
%     LQ(:) = num2cell(answer(4));
%     
%     % Enter non-calculated values into output
%     output(:,1) = group;
%     output(:,2) = coord;
%     output(:,3) = subID;
%     output(:,4) = age;
%     output(:,5) = gender;
%     output(:,6) = LQ;
    
    
%     Don't forget to uncomment this before you run it. This is because some files have more than 10
%     trials, and those extra rows will not be overwritten if you send in a file with only 10 MEPs.

    clear numTrials
%     clear group
%     clear coord
%     clear subID
%     clear age
%     clear gender
%     clear LQ
    clear condition_start
    clear condition_stop
    clear condition
    clear peak2peak
    clear area
    
    % Find the number of segments for the participant. Need to pass this to the part of the
    % code that fills in the output with the demographic variables
    numSeg(i) = max(Segment);
    clear Segment
    
end

% Fill in demographic variables with the user prompted inputs.
maxSeg = max(numSeg);

group = zeros(maxSeg,1);
group(:) = groupID;
coord = zeros(maxSeg,1);
coord(:) = answer(1);
subID = zeros(maxSeg,1);
subID(:) = subjectID;
age = zeros(maxSeg,1);
age(:) = answer(2);
gender = zeros(maxSeg,1);
gender(:) = answer(3);
LQ = zeros(maxSeg,1);
LQ(:) = answer(4);

% Enter non-calculated values into output
output(:,1) = group;
output(:,2) = coord;
output(:,3) = subID;
output(:,4) = age;
output(:,5) = gender;
output(:,6) = LQ;

output(output == 0) = NaN; % Replace padding 'zeros' with NaNs

%write data to tab del text file (this is the non-outlier corrected data)

dlmwrite([sID,'_peak_data.txt'],output);




% Run an outlier calc

lowerQ = quantile(output,0.25,1); % calc the lower quartile for each column
upperQ = quantile(output,0.75,1); % calc the upper quartile for each column
IQR = upperQ - lowerQ; % calc the interquartile range

length_output = length(output);
[numrow, numcol] = size(output);
%find indicies of the outliers
for i = 1:length_output
    outlierLow{i} = find(output(:,i) < lowerQ(i) - 1.5*IQR(i));
    outlierHi{i} = find(output(:,i) > upperQ(i) + 1.5*IQR(i));
        
end


% Replace outliers with NaNs
output_oc = output;
for i = 1:length_output
    output_oc(outlierLow{1,i},i) = NaN;
    output_oc(outlierHi{1,i},i) = NaN;
end


% Plot the 'raw' data, then mark the outliers
figure(1)
plot(7:23,[output(:,7:22),NaN(numrow,1)],'ob') % I have to pad an extra row because Matlab sucks
hold on
for i = 1:numcol
    if i < 23
        if isempty(outlierLow{1,i}) == 0
            plot(i,output(outlierLow{1,i},i),'xr');
            hold on
        end
        if isempty(outlierHi{1,i}) == 0
            plot(i,output(outlierHi{1,i},i),'xr');
            hold on
        end
    end
end
title('Peak to Peak')

figure(2)
plot(23:39, [output(:,23:end),NaN(numrow,1)],'ob');
hold on
for i = 1:numcol
    if i > 22
        if isempty(outlierLow{1,i}) == 0
            plot(i,output(outlierLow{1,i},i),'xr');
            hold on
        end
        if isempty(outlierHi{1,i}) == 0
            plot(i,output(outlierHi{1,i},i),'xr');
            hold on
        end
    end
end
title('Area')

% Writte text file to subject directory
dlmwrite([sID,'_peak_data_oc.txt'],output_oc);


