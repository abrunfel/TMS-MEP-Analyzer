% This program save data to do MANOVA analysis
clear all
close all

% Select the subject directory
str = computer;
if strcmp(str,'MACI64') == 1
    directory = uigetdir('/Volumes/mnl/Data/Adaptation/SICI_biman1/');
    % Set Currect Directory to subject dir
    cd([directory(1:45),'Post_Step_2/']);
else
    directory = uigetdir('Z:\Data\Adaptation\SICI_biman1\');
    % Set Currect Directory to subject dir
    cd([directory(1:35),'Post_Step_2\']);
end

% pp=['C:\users\fkagerer\Data\Adaptation\Kinsym2\adults\' group '\analysis\post_step2\'];
% 
% eval(['cd ' pp ]);
dir_list = dir('*_lh.mat');    %Store subject *mat data file names in variable (struct array).

dir_list = {dir_list.name}; % filenames
dir_list = sort(dir_list);  % sorts files

A2 = length(dir_list);      % how many files to process?
hand='_lh';


tstamp_start=zeros(166,1); % allocate space
tstamp_end=zeros(166,1);
target_theta = zeros(166,1);
MT_st = zeros(166,1);
rmse_st = zeros(166,1);
ide_st = zeros(166,1);
norm_jerk_st = zeros(166,1);
mov_int_st = zeros(166,1);
EPE_st = zeros(166,1);
EP_X_st = zeros(166,1);
EP_Y_st = zeros(166,1);
velPeak_st = zeros(166,1);


for B = 1:1:A2
   load(char(dir_list(B)));
   subject1 = str2num(fname(5:end));
   group1 = str2num(groupID(2:3));
   
   for trials = 1:1:166
       subject2(trials,:)=subject1;
       group2(trials,:)=group1;
       if wrong_trial(trials) == 0
           tstamp_start(trials,:) = sortData(trials).Right_FS_TimeStamp(onset(trials)); %pulls timestamp of movement onset from time matrix
           tstamp_end(trials,:) = sortData(trials).Right_FS_TimeStamp(offset(trials));  %pulls timestamp of movement offset from time matrix
           end_X_pos(trials) = sortData(trials).Left_HandX(offset(trials));
           end_Y_pos(trials) = sortData(trials).Left_HandY(offset(trials));
       else
           tstamp_start(trials) = NaN; % if trial thrown out, set all to NaN
           tstamp_end(trials) = NaN;
       end
   end
   trial = [1:1:166]';
   
   target_theta(downTrials) = 3*pi/2;
   target_theta(upTrials) = pi/2;
   MT_st(downTrials) = MT_down_st;
   MT_st(upTrials) = MT_up_st;
   rmse_st(downTrials) = rmse_down_st;
   rmse_st(upTrials) = rmse_up_st;
   ide_st(downTrials) = ide_down_st;
   ide_st(upTrials) = ide_up_st;
   norm_jerk_st(downTrials) = norm_jerk_down_st;
   norm_jerk_st(upTrials) = norm_jerk_up_st;
   mov_int_st(downTrials) = mov_int_down_st;
   mov_int_st(upTrials) = mov_int_up_st;
   EPE_st(downTrials) = EPE_down_st;
   EPE_st(upTrials) = EPE_up_st;
   EP_X_st(downTrials) = EP_X_down_st;
   EP_X_st(upTrials) = EP_X_up_st;
   EP_Y_st(downTrials) = EP_Y_down_st;
   EP_Y_st(upTrials) = EP_Y_up_st;
   velPeak_st(downTrials) = velPeak_down_st;
   velPeak_st(upTrials) = velPeak_up_st;

ALL_subjects=[group2 subject2 trial target_theta...
    MT MT_c MT_st...
    rmse rmse_c rmse_st...
    ide ide_c ide_st...
    norm_jerk norm_jerk_c norm_jerk_st...
    mov_int mov_int_c mov_int_st...
    EPE EPE_c EPE_st...
    EP_X EP_X_c EP_X_st...
    EP_Y EP_Y_c EP_Y_st...
    end_X_pos' end_Y_pos'...
    tstamp_start tstamp_end...
    velPeak velPeak_c velPeak_st...
    wrong_trial]; %You store the current matrix

if strcmp(str,'MACI64') == 1
    cd(['/Volumes/mnl/Data/Adaptation/SICI_biman1/',groupID,'/Post_Step_3']);
    dlmwrite([groupID '_lh_raw'], ALL_subjects, '-append', 'delimiter', ',', 'precision','%.6f');
    cd(['/Volumes/mnl/Data/Adaptation/SICI_biman1/',groupID,'/Post_Step_2']);
else
    cd(['Z:\Data\Adaptation\SICI_biman1\',groupID,'\Post_Step_3']);
    dlmwrite([groupID '_lh_raw'], ALL_subjects, '-append', 'delimiter', ',', 'precision','%.6f');
    cd(['Z:\Data\Adaptation\SICI_biman1\',groupID,'\Post_Step_2']);
end
% the notation '%.6f' writes each variable out to six decimal places, should get rid of engineering notation

end

% curdir=pwd;
% eval(['cd ' curdir]);
% %eval(['save C:\Jin_Bo\Visuo_Motor_Statistics\manova_dcdvscontrol_bnvse_20trials' ' age ide_bnvse1 MT_early1 MT_late1 MT_post1 rmse_bnvse1 rmse_early1 rmse_late1 rmse_post1 ide_bnvse1 ide_early1 ide_late1 ide_post1 jerk_bnvse1 jerk_early1 jerk_late1 jerk_post1 prim_bnvse1 prim_early1 prim_late1 prim_post1 total_bnvse1 total_early1 total_late1 total_post1 ide_str1']);
% %eval(['dlmwrite(''C:\Student\Crossmod\Analysis\' group '\' tasks '\ all_m.txt ',[MT_bv1 MT_bnv1 MT_early1 MT_late1 MT_aa1 MT_av1 rmse_bv1 rmse_bnv1 rmse_early1 rmse_late1 rmse_aa1 rmse_av1 ide_bv1 ide_bnv1 ide_early1 ide_late1 ide_aa1 ide_av1],'','']); 
% eval(['dlmwrite(''C:\Student\Adaptation\Crossmod\kines_e2\children\analysis\' group '\' tasks '\' 'rm_phases.txt'' ,[MT_bv1 MT_bnv1 rmse_bv1 rmse_bnv1 ide_bv1 ide_bnv11 ide_bnv21 ide_bnv31 ide_bnv41 norm_jerk_bv1 norm_jerk_bnv1 mov_int_bv1 mov_int_bnv1 mov_int_bnv11 mov_int_bnv21 mov_int_bnv31 mov_int_bnv41],'','');']);
