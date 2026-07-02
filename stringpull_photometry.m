close all
clear variables
clc
%% Import photometry data from TDT files
% Update with current folder path
data = TDTbin2mat(uigetdir);

ts = (0:(length(data.streams.x65G.data)-1))/data.streams.x65G.fs;
Fs=data.streams.x65G.fs;

%Keep original photometry timestamps but cut off parts without video
tssub=ts(round((data.epocs.PC2_.onset(2,1)*Fs):round(data.epocs.PC2_.offset(2,1)*Fs)));

%Shift timestamps to video start and stop s
tssubp=tssub-tssub(1,1);

%Set vid_frames as (total frames - Sync On timestamp) in frame_int CSV timestamp file
frame_int = readtable(uigetfile('*.csv'));
%% 
Events_1 = readtable(uigetfile('*.csv')); 
Events_2 = readtable(uigetfile('*.csv')); 
sync_on_index=Events_2(1,1)-Events_1(1,1)

% Extract the required values from Events_1 and Events_2
sync_on_time = Events_2{1, 1} - Events_1{1, 1}; % Subtract (1,1) of Events_2 and Events_1
sync_off_time = Events_2{2, 1} - Events_1{1, 1}; % Subtract (2,1) of Events_2 and (1,1) of Events_1

% Find the closest timestamps and corresponding frames for sync_on_time
[~, sync_on_index] = min(abs(frame_int.timestamp - sync_on_time));
sync_on_frame = frame_int.frame(sync_on_index);

% Find the closest timestamps and corresponding frames for sync_off_time
[~, sync_off_index] = min(abs(frame_int.timestamp - sync_off_time));
sync_off_frame = frame_int.frame(sync_off_index);

% Calculate the number of frames after "Sync On" and before "Sync Off"
vid_frames = sync_off_frame - sync_on_frame+1;
    
vid_duration=data.epocs.PC2_.offset(2,1)-data.epocs.PC2_.onset(2,1);
vid_fps=vid_frames/vid_duration;

%% Check for dropped frames
% Load the data
timestamps = frame_int.timestamp;
framerate = frame_int.fps;

% Calculate the expected interval for each frame
expected_interval = 1 ./ framerate;

% Calculate the actual interval between consecutive frames
actual_interval = [NaN; diff(timestamps)]; % Insert NaN for first frame, as it has no previous frame

% Set the threshold and tolerance
threshold_factor = 0.9999; % Set to 1 for strict detection
tolerance = 1e-5; % Adjust as needed for precision

% Detect dropped frames considering the tolerance
dropped_frames = actual_interval > (expected_interval * threshold_factor + tolerance);

% Display dropped frames
fprintf('Dropped frames:\n');
for i = 1:length(dropped_frames)
    if dropped_frames(i)
        fprintf('Timestamp: %.3f, Actual Interval: %.3f, Expected Interval: %.3f\n', ...
            timestamps(i), actual_interval(i), expected_interval(i));
    end
end

% Optional: Save dropped frames to a new table
dropped_frame_indices = find(dropped_frames);
dropped_frame_table = table(timestamps(dropped_frame_indices), actual_interval(dropped_frame_indices), expected_interval(dropped_frame_indices), ...
                            'VariableNames', {'timestamp', 'ActualInterval', 'ExpectedInterval'});
%writetable(dropped_frame_table, 'dropped_frames.csv');

%% Extract signal and control data and create plot
% Apply lowpass filter 

green=data.streams.x65G.data(round((data.epocs.PC2_.onset(2,1)*Fs):round(data.epocs.PC2_.offset(2,1)*Fs)));
red=data.streams.x60R.data(round((data.epocs.PC2_.onset(2,1)*Fs):round(data.epocs.PC2_.offset(2,1)*Fs)));
green=double(green);
green=lowpassphotometry(green,Fs,6);
red=lowpassphotometry(red,Fs,6);
red=double(red);

%Trial 1 raw plot
figure;

hold on
%plot(tssubp,ctr_green, 'color','m')
plot(tssubp,green, 'color','[0 0.5 0]')
plot(tssubp,red, 'color','r')
xlabel('Time(s)','FontSize',14)
ylabel('mV at detector','FontSize',14)
hold off
saveas(gcf, fullfile('/save/directory/here/', [data.info.blockname(1,:), '.jpg']));
%close

%% Downsample photometry data to video data

% Generate time vectors for the photometry and video data
fiberTime = linspace(0, vid_duration, length(green));
videoTime = linspace(0, vid_duration, vid_frames);

% Interpolate photometry data to match the video frame time points
green_down = interp1(fiberTime, green, videoTime, 'linear');
red_down = interp1(fiberTime, red, videoTime, 'linear');

%%
%Normalize green and red signals to average of entire trace and detrend

F0_green=mean(green_down);
[normgreen] = detrend((green_down - F0_green) ./ F0_green);
normgreenz=zscore(normgreen);

F0_red=mean(red_down);
[normred] = detrend((red_down - F0_red) ./ F0_red);
normredz=zscore(normred);

%% 
%Plot dFF
figure('units','normalized','outerposition',[0 0 1 4]);
subplot(2,1,1);
plot(videoTime, normgreen, 'color', [0 0.5 0]);
hold on
%plot(tssubp, [normIso1], 'color', 'b')
xlabel('Time(s)','FontSize',14)
ylabel('GluSnFR dF/F (normalized to mean of whole trace)','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
hold off

hold on
subplot(2,1,2);
plot(videoTime, normred, 'color', 'r');
xlabel('Time(s)','FontSize',14)
ylabel('RCaMP dF/F (normalized to mean of whole trace)','FontSize',14)
hold off
fontsize(16,"points")
saveas(gcf, fullfile('/save/directory/here/dFF plot/', [data.info.blockname(1,:), '.jpg']));
%close

%% Import DLC CSV file
% go to the folder where the data is that you want to analyse
% add the filename
file = uigetfile('.csv');
DLC_data = readmatrix(file);
dlcpos = readtable(file);
filename = file(1:21);
pix_mm = 6 / 57.87 ;  %  conversion: width of photometry connector (mm) to pixels from imagej

%Shift DLC CSV to start from Sync On timestamp 
dlcpos = dlcpos(sync_on_index:sync_off_index,:);
%% %% replace low probability points for the left paw
% for each bodypart the name is x position, 1 is y position and 2 is probability
%Position of Left Paw
LeftPawx = DLC_data(:,2);
LeftPawy = DLC_data(:,3);
LeftPawp = DLC_data(:,4);
RealPositionsl = [LeftPawy LeftPawp]; % isolate y position and probability
        for n=1:length(RealPositionsl)   
            if RealPositionsl(n,2) < 0.8
                RealPositionsl(n,1) = NaN;
            end
        end      
figure,
    plot (RealPositionsl(:,1))
    title ('left paw good probability')
 
RealPositionsl(:,1) = fillmissing(RealPositionsl(:,1),'linear');
RealPositionsl = RealPositionsl(sync_on_index:sync_off_index,:);
figure,
    plot (RealPositionsl(:,1))
    title ('left paw linear fill')
%% same for the right side
RightPawx = DLC_data(:,5);
RightPawy = DLC_data(:,6);
RightPawp = DLC_data(:,7);
RealPositionsr = [RightPawy RightPawp];% isolate y position and probability
        for n=1:length(RealPositionsr)   
            if RealPositionsr(n,2) < 0.8
                RealPositionsr(n,1) = NaN;
            end
        end      
figure,
    plot (RealPositionsr(:,1))
    title ('right paw good probability')
 
RealPositionsr(:,1) = fillmissing(RealPositionsr(:,1),'linear');
RealPositionsr = RealPositionsr(sync_on_index:sync_off_index,:);
figure,
    plot (RealPositionsr(:,1))
    title ('right paw linear fill')
%% same for the nose
Nosex = DLC_data(:,8);
Nosey = DLC_data(:,9);
Nosep = DLC_data(:,10);
RealPositionsn = [Nosey Nosep];% isolate y position and probability
        for n=1:length(RealPositionsn)   
            if RealPositionsn(n,2) < 0.7
                RealPositionsn(n,1) = NaN;
            end
        end      
figure,
    plot (RealPositionsn(:,1))
    title ('nose good probability')
 
RealPositionsn(:,1) = fillmissing(RealPositionsn(:,1),'linear');
RealPositionsn = RealPositionsn(sync_on_index:sync_off_index,:);
figure,
    plot (RealPositionsn(:,1))
    title ('nose linear fill')
%%  check relative positions

figure('units','normalized','outerposition',[0 0 1 10]);
subplot(2,1,1);
    hold on
    plot(videoTime,RealPositionsn(:,1)); 
    plot(videoTime,RealPositionsr(:,1)); 
    plot(videoTime,RealPositionsl(:,1));
    hold off
    xlabel('Time (sec)')
    ylabel ('y position (pix)')
    %ylim([0 900])
    title (filename,'interpreter','none')
    legend ('nose','RPaw','LPaw')
    ax = gca;
    ax.YDir = 'reverse';
    fontsize(16,"points")
saveas(gcf,['/save/directory/here/Y position graphs/',filename,'.jpg']);
%% Define the start and stop times for the pull interval
% enter the start and stop time in seconds that includes all pulls
 % has to start before the first top peak on the right and the first bottom
 % on the leftwai
aligned_pulltimes = [4.5 10.2]; % manually adjust this for each trial (in seconds)

pull_start_time = aligned_pulltimes(1,1); 
pull_end_time = aligned_pulltimes(1,2);   

% Convert pull times to frame indices 
pull_start_frame = round(pull_start_time * vid_fps);
pull_end_frame = round(pull_end_time * vid_fps);

% Extract data within the adjusted pulltime range
pull_indices = pull_start_frame:pull_end_frame;
%% normalize position so that values are distance from the nose
% using only high probability points
% position below nose for peaks at the top of the trajectory
rp_xpos = RealPositionsn(pull_indices,2) - RealPositionsr (pull_indices,2); 
rp_n = RealPositionsn(pull_indices,1) - RealPositionsr (pull_indices,1); 
% position below nose for peaks at the bottom of the trajectory
rp_na = abs(RealPositionsn(pull_indices,1) - RealPositionsr(pull_indices,1));
rp_mm = rp_na * pix_mm;

k = kurtosis(rp_mm);  % normal distribution is 3 but over 3 means outliers
% same as above for the left side
lp_n = RealPositionsn(pull_indices,1) - RealPositionsl(pull_indices,1);
lp_na = abs(RealPositionsn(pull_indices,1) - RealPositionsl(pull_indices,1));
lp_mm = lp_na * pix_mm;
%%
figure,
hold on
    plot (pull_indices, rp_mm)
    plot (pull_indices, lp_mm)
hold off
xlabel('pull time (sec)')
ylabel ('y pos below nose (mm)')

title (filename,'interpreter','none')
legend ('Rpaw','Lpaw')
ax = gca;
ax.YDir = 'reverse';
saveas(gcf,['/save/directory/here/Y position graphs norm/',filename,'.jpg']);
save(['/save/directory/here/pulltime/', filename],'aligned_pulltimes', '-ascii'); %% %% find peaks for pull trajectory of the right paw
%% 

% the bottom of the trajectories
[pksbr,locsbr] = (findpeaks(rp_na,'MinPeakWidth', 0.5, 'MinPeakProminence', 20, 'annotate','extents'));
figure;findpeaks(rp_na,'MinPeakWidth', 0.5,'MinPeakProminence', 20,'annotate','extents')
% the top of the trajectories
[pkstr,locstr] = (findpeaks(rp_n,'MinPeakWidth', 0.5,'MinPeakProminence', 20, 'annotate','extents'));
figure;findpeaks(rp_n,'MinPeakWidth', 0.5,'MinPeakProminence', 20,'annotate','extents')

%% %% for the pull down right paw

% if locsb(1) < locst (1)  only if there is a bottom peak first
%    locsb = locsb(2:end); only 
% end
rpull = cell (2,length(pkstr)-1); %make space for pull times and values
rpulldur = zeros (1,length(pkstr)-1);
rpullmag = zeros (1,length(pkstr)-1);
rpullsp = zeros (1,length(pkstr)-1);

for i = 1:length(pkstr)-1 %if the first bottom is before the top in locsb +1
rpull{1,i} = locstr(i) : locsbr (i);% all the pull down frames
rpull{2,i} = rp_n(locstr(i) : locsbr (i));
rpulldur(1,i) = (length (rpull{2,i}))/vid_fps; % duration in s
rpullmag (1,i) = (abs((rp_n(locstr(i))) - (rp_n(locsbr (i)))))*pix_mm; %
rpullsp(1,i) = rpullmag(i)/rpulldur(i);
end
mean_rdur = mean(rpulldur); % average duration
mean_rpull = mean(rpullmag); % average size in mm
mean_rpullsp = mean(rpullsp);
std_rpullsp = std(rpullsp);
std_rpull = std(rpullmag);

%% for the sweep up right paw
rup = cell (2,length(pkstr)-1); %make space for pull times and values
r_up_dur = zeros (1,length(pkstr)-1);
rupmag = zeros (1,length(pkstr)-1);
rupsp = zeros (1,length(pkstr)-1);

for i = 1:(length(pkstr)-1)
rup{1,i} = locsbr(i) : locstr (i+1);% all the up frames
rup{2,i} = rp_n(locsbr(i) : locstr (i+1));
r_up_dur(1,i) = (length (rup{2,i}))/vid_fps; % duration in s
rupmag (1,i) = (abs((rp_n(locsbr(i))) - (rp_n(locstr (i+1)))))*pix_mm; %change in pix
rupsp(1,i) = rupmag(i)/r_up_dur(i);
end
rup_dur_mean = mean(r_up_dur); % average duration
mean_r_up = mean(rupmag);
std_r_up = std(rupmag);% average size in mm
mean_rupsp = mean(rupsp);
std_rupsp = std(rupsp);

stdr=zeros(4,1);
stdr(1,1)=std_rpullsp;
stdr(2,1)=std_rpull;
stdr(3,1)=std_rupsp;
stdr(4,1)=std_r_up;
save(['/save/directory/here/right pull SD/', data.info.blockname(1:21)],'stdr', '-ascii');

durr=zeros(2,1);
durr(1,1)=mean_rdur;
durr(2,1)=rup_dur_mean;
save(['/save/directory/here/right pull duration/', data.info.blockname(1:21)],'durr', '-ascii');

%% find peaks for pull trajectory of the left paw

% the bottom of the trajectories
[pksbl,locsbl] = (findpeaks(lp_na, 'MinPeakWidth', 0.5,'MinPeakProminence', 20, 'annotate','extents'));%,'MinPeakDistance', 5,'MinPeakProminence', 50, 'annotate','extents'));
figure;findpeaks(lp_na,'MinPeakWidth', 0.5,'MinPeakProminence', 20, 'annotate','extents');%,'MinPeakWidth', 5,'MinPeakProminence', 50,'annotate','extents')
% the top of the trajectories
[pkstl,locstl] = (findpeaks(lp_n, 'MinPeakWidth', 0.5,'MinPeakProminence', 20,'annotate','extents'));%,'MinPeakWidth', 5, 'annotate','extents'));
figure;findpeaks(lp_n, 'MinPeakWidth', 0.5,'MinPeakProminence', 20,'annotate','extents');%,'MinPeakWidth', 5,'annotate','extents')

%% %% for the pull down left paw
lpull = cell (2,length(pksbl)-1); %make space for pull times and values
lpulldur = zeros (1,length(pksbl)-1);
lpullmag = zeros (1,length(pksbl)-1);
lpullsp = zeros (1,length(pksbl)-1);

for i = 1:length(pksbl)-1  %if the first bottom is before the top in locsb +1
    lpull{1,i} = locstl(i) : locsbl (i+1);% all the pull down frames
    lpull{2,i} = lp_n(locstl(i) : locsbl (i+1));
    lpulldur(1,i) = (length (lpull{2,i}))/vid_fps; % duration in s
    lpullmag (1,i) = abs ((min(lpull{2,i})-(max(lpull{2,i}))))*pix_mm; %change in pix converted
    lpullsp(1,i) = lpullmag(i)/lpulldur(i); %speed in mm/s
end
mean_ldur = mean(lpulldur); % average duration
mean_lpull= mean(lpullmag); % average size in mm
mean_lpullsp = mean(lpullsp);
std_lpull= std(lpullmag); % average size in mm
std_lpullsp = std(lpullsp);

%% for the sweep up left paw
lup = cell (2,length(pkstl)); %make space for pull times and values
l_up_dur = zeros (1,length(pkstl));
lupmag = zeros (1,length(pkstl));
lup_sp = zeros (1,length(pkstl));
for i = 1:length(pkstl)-1
lup{1,i} = locsbl(i) : locstl (i);% all the pull down frames
lup{2,i} = lp_n(locsbl(i) : locstl (i));
l_up_dur(1,i) = (length (lup{2,i}))/vid_fps; % duration in s
lupmag (1,i) = (abs((lp_n(locsbl(i))) - (lp_n(locstl (i)))))*pix_mm; %change in pix
lup_sp (1,i) = lupmag(i)/l_up_dur(i);  % speed in mm/s
end
lup_dur_mean = mean(l_up_dur); % average duration
mean_l_up = (mean(lupmag)); % average size in mm
mean_lupsp = mean(lup_sp);
std_l_up= std(lupmag);
std_lupsp = std(lup_sp);

stdl=zeros(4,1);
stdl(1,1)=std_lpullsp;
stdl(2,1)=std_lpull;
stdl(3,1)=std_lupsp;
stdl(4,1)=std_l_up;
save(['/save/directory/here/left pull SD/', data.info.blockname(1:21)],'stdl', '-ascii');

durl=zeros(2,1);
durl(1,1)=mean_ldur;
durl(2,1)=lup_dur_mean;
save(['/save/directory/here/left pull duration/', data.info.blockname(1:21)],'durl', '-ascii');
%% 
toplpull = max(lpullmag);
topslpull = max(lpullsp);
toplup = max(lupmag);
topslup = max(lup_sp);

toprpull = max(rpullmag);
topsrpull = max (rpullsp);
toprup = max(rupmag);
topsrup = max(rupsp);

meansize = [mean_rpull ,mean_lpull,mean_r_up,mean_l_up];
writematrix(meansize, fullfile('/save/directory/here/mean size/', strcat(filename, 'meansize')));
meanspeed = [mean_rpullsp ,mean_lpullsp,mean_rupsp,mean_lupsp];
writematrix(meanspeed, fullfile('/save/directory/here/mean speed/', strcat(filename, 'meanspeed')));
%% correlation of left and right paws
rlcor = corr(rp_na,lp_na);
[r,lags] = xcorr(rp_na,lp_na, 'normalized');

% Convert lags from samples to seconds
timeLags = lags / vid_fps;

% Plot the cross-correlation with time lags in seconds
figure;
plot(timeLags, r, 'b', 'LineWidth', 1.5);
xlabel('Lag (seconds)');
ylabel('Cross-Correlation');
title('Cross-Correlation between Left and Right Paw Pulls');
grid on;

[max_r, max_index] = max(r);
max_lag = timeLags(max_index);

directory2 = '/save/directory/here/left_right cross correlation/';

% Save the plot as an image with filename included
plotFilePath = fullfile(directory2, strcat(filename, '_cross_correlation_plot.png'));
saveas(gcf, plotFilePath);

% Save max_r, max_lag and rlcor as a CSV file with filename included
summaryData = table(max_r, max_lag, rlcor, 'VariableNames', {'MaxCorrelation', 'MaxLag', 'Correlation'});
summaryFilePath = fullfile(directory2, strcat(filename, '_cross_correlation_summary.csv'));
writetable(summaryData, summaryFilePath);


%% Plot photometry and DLC data during pull bouts only
% Create a new figure for the plots
figure('units','normalized','outerposition',[0 0 1 1]);

% Plot photometry data in the first subplot
subplot(2, 1, 1);
hold on;
plot(videoTime(pull_indices), normgreenz(pull_indices), 'color', [0 0.5 0], 'DisplayName', 'GluSnFR');
plot(videoTime(pull_indices), normredz(pull_indices), 'color', 'r', 'DisplayName', 'RCaMP');
xlabel('Time (s)', 'FontSize', 14);
ylabel('Z-score', 'FontSize', 14);
title('Photometry Signals During Adjusted Pull Window', 'FontSize', 16);
legend('Location', 'best');
grid on;
hold off;

% Plot RealPositions data in the second subplot
subplot(2, 1, 2);
hold on;
plot(videoTime(pull_indices), rp_mm, 'm', 'DisplayName', 'Right Paw');
plot(videoTime(pull_indices), lp_mm, 'c', 'DisplayName', 'Left Paw');
xlabel('Time (s)', 'FontSize', 14);
ylabel('Y Position (pix)', 'FontSize', 14);
title('Body Part Positions During Adjusted Pull Window', 'FontSize', 16);
legend('Location', 'best');
ax = gca;
ax.YDir = 'reverse'; % Reverse y-axis for RealPositions
grid on;
hold off;

sgtitle('Photometry and RealPositions Data (Aligned to Sync On)', 'FontSize', 18); 
saveas(gcf, fullfile('/save/directory/here/pulls_photometry_plot/', [data.info.blockname(1,:), '.jpg']));

%% %% Compute mean z-score of pull bout
normgreenz_bout_mean = mean(normgreenz(pull_indices));
normredz_bout_mean = mean(normredz(pull_indices));
mean_pull_zscore = [normgreenz_bout_mean normredz_bout_mean];
writematrix(mean_pull_zscore, fullfile('/save/directory/here/mean_pull_zscore/', strcat(filename, 'meanpullzscore')));
%% Detect GluSnFR photometry peaks 
boutduration = pull_end_time - pull_start_time; %in seconds
save(['/save/directory/here/bout duration/', data.info.blockname(1:21)],'boutduration', '-ascii');

%Check peakfinder threshold in whole trace (acts as baseline)
[pksgbLP,locsgbLP,wgbLP,pgbLP] = findpeaks(normgreenz,'MinPeakProminence', 2, 'MinPeakWidth', 4);
figure;
plot(videoTime,normGreen1zLP,videoTime(locsgbLP),pksgbLP,"o")
xlabel('Time(s)','FontSize',14)
ylabel('GluSnFR z-score','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')

mwgbLP = mean(wgbLP);
mpgbLP = mean(pgbLP);
maxgbLP = max(pgbLP);
freqgbLP = numel(locsgbLP)/vid_duration;
%Set values as 0 if no peaks detected
mwgbLP(isnan(mpgbLP))=0;
maxgbLP(isnan(mpgbLP))=0;
mpgbLP(isnan(mpgbLP))=0;

% During pull bout 
[pksgpLP,locsgpLP,wgpLP,pgpLP] = findpeaks(normgreenz(pull_indices),'MinPeakProminence', 2, 'MinPeakWidth', 4);
mwgpLP = mean(wgpLP);
mpgpLP = mean(pgpLP);
maxgpLP = max(pgpLP);
freqgpLP = numel(locsgpLP)/boutduration; %in peaks/s
%Set values as 0 if no peaks detected
mwgpLP(isnan(mpgpLP))=0;
maxgpLP(isnan(mpgpLP))=0;
mpgpLP(isnan(mpgpLP))=0;

figure;
plot(pull_indices,normgreen(pull_indices),pull_indices(locsgpLP),pksgpLP,"o")
xlabel('Time(frames)','FontSize',14)
ylabel('GluSnFR z-score during pull bout','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')

%% %% Detect RCaMP photometry peaks 

%Check peakfinder threshold in whole trace (acts as baseline)
[pksrb,locsrb,wrb,prb] = findpeaks(normredz,'MinPeakProminence', 0.5, 'MinPeakWidth', 1);
figure;
plot(videoTime,normredz,videoTime(locsrb),pksrb,"o")
xlabel('Time(s)','FontSize',14)
ylabel('RCaMP z-score','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')

mwrb = mean(wrb);
mprb = mean(prb);
maxrb = max(prb);
freqrb = numel(locsrb)/vid_duration;
%Set values as 0 if no peaks detected
mwrb(isnan(mprb))=0;
maxrb(isnan(mprb))=0;
mprb(isnan(mprb))=0;

% During pull bout 
[pksrp,locsrp,wrp,prp] = findpeaks(normredz(pull_indices),'MinPeakProminence', 0.5, 'MinPeakWidth', 1);
mwrp = mean(wrp);
mprp = mean(prp);
maxrp = max(prp);
freqrp = numel(locsrp)/boutduration; %in peaks/s
%Set values as 0 if no peaks detected
mwrp(isnan(mprp))=0;
maxrp(isnan(mprp))=0;
mprp(isnan(mprp))=0;
figure;
plot(pull_indices,normredz(pull_indices),pull_indices(locsrp),pksrp,"o")
xlabel('Time(frames)','FontSize',14)
ylabel('RCaMP z-score during pull bout','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')


peakfreq=zeros(4,1);
peakfreq(1,1)=freqgpLP;
peakfreq(2,1)=freqgbLP;
peakfreq(3,1)=freqrp;
peakfreq(4,1)=freqrb;
save(['/save/directory/here/peak frequency/', data.info.blockname(1:21)],'peakfreq', '-ascii');

meanpeakprom=zeros(4,1);
meanpeakprom(1,1)=mpgpLP;
meanpeakprom(2,1)=mpgbLP;
meanpeakprom(3,1)=mprp;
meanpeakprom(4,1)=mprb;
save(['/save/directory/here/mean peak prominence/', data.info.blockname(1:21)],'meanpeakprom', '-ascii');

maxpeakprom=zeros(4,1);
maxpeakprom(1,1)=maxgpLP;
maxpeakprom(2,1)=maxgbLP;
maxpeakprom(3,1)=maxrp;
maxpeakprom(4,1)=maxrb;
save(['/save/directory/here/max peak prominence/', data.info.blockname(1:21)],'maxpeakprom', '-ascii');

meanpeakwidth=zeros(4,1);
meanpeakwidth(1,1)=mwgpLP;
meanpeakwidth(2,1)=mwgbLP;
meanpeakwidth(3,1)=mwrp;
meanpeakwidth(4,1)=mwrb;
save(['/save/directory/here/mean peak width/', data.info.blockname(1:21)],'meanpeakwidth', '-ascii');

%% %% %% Align onset of pull bout (GluSnFR)

pre_event_window = 3; % Seconds before pull onset
post_event_window = 3; % Seconds after pull onset

    % Convert time window to sample points
    pre_event_samples = round(pre_event_window * vid_fps);
    post_event_samples = round(post_event_window * vid_fps);

% Extract the signal for the window around pull onset
signal_window_g = normgreenz(pull_start_frame-pre_event_samples : pull_start_frame+post_event_samples);

% Time axis for plotting (-3 to +3 seconds)
time_axis = linspace(-pre_event_window, post_event_window, length(signal_window_g));     
 
% Plotting the result
figure;
plot(time_axis, signal_window_g, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('GluSnFR zscore');
title(data.info.blockname(1,:));
grid on;
saveas(gcf, fullfile('/save/directory/here/GluSnFR pull onset plot/', [data.info.blockname(1,:), '.jpg']));
writematrix(signal_window_g, fullfile('/save/directory/here/GluSnFR pull onset/', strcat(filename, 'GluSnFR_pull_onset.csv')));

%% %% %% Align onset of pull bout (RCaMP)
% Extract the signal for the window around pull onset
signal_window_r = normredz(pull_start_frame-pre_event_samples : pull_start_frame+post_event_samples);
 
% Plotting the result
figure;
plot(time_axis, signal_window_r, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('RCaMP zscore');
title(data.info.blockname(1,:));
grid on;
saveas(gcf, fullfile('/save/directory/here/RCaMP pull onset plot/', [data.info.blockname(1,:), '.jpg']));
writematrix(signal_window_r, fullfile('/save/directory/here/RCaMP pull onset/', strcat(filename, 'RCaMP_pull_onset.csv')));

%% Compute SD of GluSnFR & RCaMP dF/F during baseline and pull bout 
normgreen_SDpullbout = std(normgreen(pull_indices))
%normgreen_SDbaseline = std(normgreen(1:30*vid_fps))

normred_SDpullbout = std(normred(pull_indices))
%normred_SDbaseline = std(normred(1:30*vid_fps))

dFF_SD=zeros(2,1);
dFF_SD(1,1)=normgreen_SDpullbout;
%dFF_SD(2,1)=normgreen_SDbaseline;
dFF_SD(2,1)=normred_SDpullbout;
%dFF_SD(4,1)=normred_SDbaseline;
save(['/save/directory/here/dFF_SD/', data.info.blockname(1:21)],'dFF_SD', '-ascii');

%% Correlate GluSnFR and RCaMP during pull bout and for whole trace
normgreent=normgreen.';
normredt=normred.';
drcor_pull = corr(normgreent(pull_indices),normredt(pull_indices));
drcor_tot = corr(normgreent,normredt);

drcor = zeros(2,1);
drcor(1,1)=drcor_pull;
drcor(2,1)=drcor_tot;
save(['/save/directory/here/GluSnFR_RCaMP_correlation/', data.info.blockname(1:21)],'drcor', '-ascii');

[rphot,lagsphot] = xcorr(normgreent(pull_indices),normredt(pull_indices), 'normalized');

% Convert lags from samples to seconds
timeLagsphot = lagsphot / vid_fps;

% Plot the cross-correlation with time lags in seconds
figure;
plot(timeLagsphot, rphot, 'b', 'LineWidth', 1.5);
xlabel('Lag (seconds)');
ylabel('Cross-Correlation');
title('Cross-Correlation between GluSnFR and RCaMP during pull bout');
grid on;

[max_r_phot, max_index_phot] = max(rphot);
max_lag_phot = timeLags(max_index_phot);

directory3 = '/save/directory/here/GluSnFR RCaMP cross correlation plot/';

% Save the plot as an image with filename included
plotFilePath1 = fullfile(directory3, strcat(filename, '_GluSnFRRCaMP_crosscorr_plot.png'));
saveas(gcf, plotFilePath1);

% Save max_r_phot and max_lag_phot 
summarydata = table(max_r_phot, max_lag_phot, 'VariableNames', {'MaxCorrelation', 'MaxLag'});
summaryFilePath1 = fullfile(directory3, strcat(filename, '_GluSnFRRCaMP_crosscorr_summary.csv'));
writetable(summarydata, summaryFilePath1);

%% %% Compute correlations between GluSnFR and RCaMP during pull bout in 1s moving window sections
normGreen1tpullboutLP=normgreent(pull_indices);
normredtpullbout=normredt(pull_indices);

% Define window size in samples
window_size = round(vid_fps); % 1-second window
half_window = floor(window_size / 2); % Half window for centering
num_samples = length(normGreen1tpullboutLP); % Total number of samples

% Zero-padding the signals at the edges
%padded_green = [zeros(half_window, 1); normGreen1tpullboutLP; zeros(half_window, 1)];
%padded_red = [zeros(half_window, 1); normredtpullbout; zeros(half_window, 1)];

% Initialize correlation storage
corr_values = nan(num_samples,1); 
time_axis = (1:num_samples) / vid_fps; % Convert indices to time in seconds

% Compute correlation with a centered moving window
for i = 1:num_samples
    % Define window range
    start_idx = max(1, i - half_window);
    end_idx = min(num_samples, i + half_window);
    
    % Extract windowed segments
    window_green = normGreen1tpullboutLP(start_idx:end_idx);
    window_red = normredtpullbout(start_idx:end_idx);

    % Remove NaNs from the window
    valid_idx = ~isnan(window_green) & ~isnan(window_red);
    window_green = window_green(valid_idx);
    window_red = window_red(valid_idx);

    % Compute correlation only if enough valid data points exist
    if length(window_green) > 2 % Avoid empty or single-point correlations
        corr_values(i) = corr(window_green, window_red);
    end
end

% Find the maximum correlation and its corresponding time (ignoring NaNs)
[max_corr, max_idx] = max(corr_values, [], 'omitnan'); % Ignore NaNs
max_time = time_axis(max_idx);

% Plot correlation over time
figure;
plot(time_axis, corr_values, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Correlation');
title('GluSnFR & RCaMP Moving Window Correlation (1s)');
grid on;


% Save the plot and max correlation (1s) during pull bout
plotFilePath3 = fullfile('/save/directory/here/GluSnFR RCaMP 1s correlation plot/', strcat(filename, '_GluSnFRRCaMP_1scorr_plot.png'));
saveas(gcf, plotFilePath3);

drcor_pullLP1s=max(corr_values);
save(['/save/directory/here/GluSnFR RCaMP 1s correlation/', data.info.blockname(1:21)],'drcor_pullLP1s', '-ascii');
save(['/save/directory/here/matlab files/', data.info.blockname(1:21)]);