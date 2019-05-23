clc; clear all;

% This code is applicable to result of all crack growth tests.
disp(datetime('now'));

% Specify all file locations here for convenience. 
% AEfilename_root = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\CTAA04_020419\AE\';
% AEfilename_root = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\CTAA07_020619\AE\';
AEfilename_root = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\CTAA08_020719\AE\';

% Instfilename = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\CTAA04_020419\Instron\Test2\Test2.steps.tracking.csv';
% Instfilename = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\CTAA07_020619\Instron\Test2\Test2.steps.tracking.csv';
Instfilename = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\CTAA08_020719\Instron\Test2\Test2.steps.tracking.csv';

% Specify number of Para data files for each test.
para_num = 8; % CTAA04: 10 -- CTAA07: 2 -- CTAA08: 8 
              
% Lower bound used to find match points in AE and Instron data
% LB = 1; LB2 = 1.9; 
LB = 1; LB2 = 3; % CTAA08

fileID = fopen([AEfilename_root,'Results\Para_Modified.txt'],'w');
fileID2 = fopen([AEfilename_root,'Results\Para_hit.txt'],'w'); % keep hit correction time and cycles data here

% Measure the time for importing data
tic;

for i=1:para_num
    
    if i < 10
        para_filename = ['Para_growth_0',num2str(i),'.TXT'];
    else
        para_filename = ['Para_growth_',num2str(i),'.TXT'];
    end
    
    AEfilename = [AEfilename_root,para_filename];
    AEData = importdata(AEfilename,' ',8);
       
    AEData_t = AEData.data(:,2);
    AEData_l = AEData.data(:,3);
    
    if i==1
        AEData_time = AEData_t;
        AEData_load = AEData_l;
    else
        AEData_time = [AEData_time; AEData_t];
        AEData_load = [AEData_load; AEData_l];
    end
    
end

% Multiply the load values by the multipication factor
AEData_load = AEData_load*10; 

% import csv data from Instron
InstronData = dlmread(Instfilename,',',1,0);

% Measure Instron data import time
t_AEInst = toc;

% Display a message to know when Instron data is imported
disp('AE and Instron data are imported !');
disp(['Elapsed time to import data is ',num2str(t_AEInst,'%.1f'),' seconds'])
disp ('------------------------------------------------------------------')

% explain headers of each column Instron data - Since instead of importdata
% we have to use dlmread for large csv file, we do not read headers and we
% can show them.
% disp (['Headers in Instron data are:' ,char(InstronData1.textdata)]);

% explain headers of each column AE data
disp (['Headers in AE data are:' ,AEData.textdata{5,1},' ',AEData.textdata{5,2},'(s)',' ',AEData.textdata{5,3},'(kN)',' ',AEData.textdata{5,4},'(mm)']);
disp ('------------------------------------------------------------------')

% Measure the overall time it took for the script to run
tic; 

% Merge data.
InstronData_data = InstronData;

% Plot load vs. time from both data sources in one figure
f1 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot(InstronData_data(:,1),InstronData_data(:,9),'b-',AEData_time,AEData_load,'r-')
title('Load vs. Time (Raw data)')
xlabel('Time (S)')
ylabel('Load (kN)')
legend('Instron Data','AE Data','Location','northeast','Orientation','horizontal')
print(f1,'Load vs. Time (Raw data).tiff','-dtiffn','-r300')

% Transfer AE Data to a new variable for time delay remove. It is neccesary
% to make us able for comparing the delay removed and unremoved data.
AEData_time_c = [AEData_time AEData_load];

% First find the intial timedelay in AE data
AEcut_1 = find(AEData_load>0.8,1);
Instcut_1 = find(InstronData_data(:,9)>0.8,1);
timedelay_initial = AEData_time(AEcut_1) - InstronData_data(Instcut_1,1);

% Shift all the AE data with timedelay
for i=1:length(AEData_time_c)
    AEData_time_c(i,1) = AEData_time_c(i,1) - timedelay_initial;
end

% find a point with time greater than zero (or 1 second) to remove negative
% value times
AEcut_2 = find(AEData_time_c(:,1) > 1,1); 
t=0;
for i=1:AEcut_2     % remove AE data that has negative time value
    i=i-t;
    if AEData_time_c(i,1) < 0
        AEData_time_c(i,:)=[];
        t=t+1;
    end
end

%Find the first time correction by matching the first time Instron and AE
%reach first value above LB
AEcut_3 = find(AEData_time_c(:,2)>LB,1);
Instcut_3 = find(InstronData_data(:,9)>LB,1);
timedelay_first = AEData_time_c(AEcut_3,1) - InstronData_data(Instcut_3,1);

% % Shift all the AE data with timedelay
for i=1:length(AEData_time_c)
    AEData_time_c(i,1) = AEData_time_c(i,1) - timedelay_first;
end

%Find index points in Instron data UNDER 50% max load
MaxLoadInst = find(InstronData_data(:,9)>LB); % Use this when AE load vs time data fluctuates

q=1; % Keep initial points at UNDER 50% max loads
for i=2:length(MaxLoadInst)
    if InstronData_data(MaxLoadInst(i),1) - InstronData_data(MaxLoadInst(i-1),1) > 5
        DelayIndexInst(q,:) = MaxLoadInst(i);
        q=q+1;
    end
end

%Find index points in Instron data OVER 50% max load
MaxLoadInst2 = find(InstronData_data(:,9)>LB2); % Use this when AE load vs time data fluctuates

q=1; % Keep initial points at OVER 50% max loads
for i=2:length(MaxLoadInst2)
    if InstronData_data(MaxLoadInst2(i),1) - InstronData_data(MaxLoadInst2(i-1),1) > 5
        DelayIndexInst2(q,:) = MaxLoadInst2(i);
        q=q+1;
    end
end

%Find index points in Instron data OVER 50% max load
MaxLoadInst3 = find(InstronData_data(:,9)>LB2); % Use this when AE load vs time data fluctuates

q=1; % Keep initial points at OVER 50% max loads
for i=2:length(MaxLoadInst3)
    if InstronData_data(MaxLoadInst3(i),1) - InstronData_data(MaxLoadInst3(i-1),1) > 0.5 && InstronData_data(MaxLoadInst3(i),1) - InstronData_data(MaxLoadInst3(i-1),1) < 5
        DelayIndexInst3(q,:) = MaxLoadInst3(i-1);
        q=q+1;
    end
end


n=1;
while 1
    
    MaxLoadAE = find(AEData_time_c(:,2)>LB); %Find index points in AE data
        
    q=1; % Keep initial points at UNDER 50% max loads
    for i=2:length(MaxLoadAE)
        if AEData_time_c(MaxLoadAE(i),1) - AEData_time_c(MaxLoadAE(i-1),1) > 2
            DelayIndexAE(q,:) = MaxLoadAE(i);
            q=q+1;
        end
    end
       
       
    % Check if same number of indexes are taken from AE and Instron
    if length(DelayIndexAE)==length(DelayIndexInst)
        disp([num2str(length(DelayIndexAE)),' data points of AE and Instron data are used to match time indexes.'])
        disp(['Iteration #: ',num2str(n)])
    else
        disp('Error !! number of data points counted from AE and Instron (to match time index) does not match.')
        break
    end
    disp ('-----------------')
    
    % Find time delays in AE by comparing Delay Index times of AE and Instron
    for i=1:length(DelayIndexAE)
      timedelay(i,:) = AEData_time_c(DelayIndexAE(i),1)-InstronData_data(DelayIndexInst(i),1);
    end
    
    timedelay_hist(:,n) = timedelay; % Keep a history of timedelays
    
    % Delete redundant points recorded in AE data
    AEaqrate = 50; % AE data acquisition rate per second

    IndexShift = DelayIndexAE(n)-fix(timedelay(n).*AEaqrate);
%     IndexShift = fix(AEData_time(DelayIndexAE(n),2).*AEaqrate)-fix(timedelay(n).*AEaqrate);
    AEData_time_c(IndexShift:DelayIndexAE(n)-1,:)=[];
    for i=IndexShift:length(AEData_time_c)
        AEData_time_c(i,1) = AEData_time_c(i,1) - timedelay(n);
    end
    
% % % % % % % % % %  Following part of the code removes the redundant time over 50% max load % % % % % % % % % % 
    
    %Find index points OVER 50% max load in AE data
    MaxLoadAE2 = find(AEData_time_c(:,2)>LB2); 
        
    q=1; % Keep initial points at OVER 50% max loads
    for i=2:length(MaxLoadAE2)
        if AEData_time_c(MaxLoadAE2(i),1) - AEData_time_c(MaxLoadAE2(i-1),1) > 4
            DelayIndexAE2(q,:) = MaxLoadAE2(i);
            q=q+1;
        end
    end
        
    if length(DelayIndexAE2) ~= length(DelayIndexInst2)
        disp('Error !! number of data points counted from AE2 and Instron2 (to match time index) does not match.')
        break
    end
    
    % Find time delays in AE by comparing Delay Index times of AE and Instron
    for i=1:length(DelayIndexAE2)
      timedelay2(i,:) = AEData_time_c(DelayIndexAE2(i),1)-InstronData_data(DelayIndexInst2(i),1);
    end
    
    timedelay_hist2(:,n) = timedelay2; % Keep a history of timedelays
    
    IndexShift2_1 = DelayIndexAE2(n)-fix(timedelay2(n).*AEaqrate);
    AEData_time_c(IndexShift2_1:DelayIndexAE2(n)-1,:)=[];
    for i=IndexShift2_1:length(AEData_time_c)
        AEData_time_c(i,1) = AEData_time_c(i,1) - timedelay2(n);
    end
    
    %Find index points OVER 50% max load in AE data
    MaxLoadAE2 = find(AEData_time_c(:,2)>LB2); 
    
    q=1; % Keep initial points at OVER 50% max loads
    for i=2:length(MaxLoadAE2)
        if AEData_time_c(MaxLoadAE2(i),1) - AEData_time_c(MaxLoadAE2(i-1),1) > 4
            DelayIndexAE2(q,:) = MaxLoadAE2(i);
            q=q+1;
        end
    end
    
    %Find index points in AE data with load values less than delay index AE 2 load in AE data
    MaxLoadAE3 = find(AEData_time_c(:,2) - AEData_time_c(DelayIndexAE2(n),2) < 0);
    
    q=1; % Select the ones that happen after 1 second of delay index
    for i=1:length(MaxLoadAE3)
        if AEData_time_c(MaxLoadAE3(i),1) - AEData_time_c(DelayIndexAE2(n),1) > 1
            t(q,:) = MaxLoadAE3(i);
            q=q+1;
        end
    end
    
    % The first one would be Delay Index AE 3
    DelayIndexAE3(n,:) = t(1);
    
    
    for i=1:length(DelayIndexAE3)
      timedelay3(i,:) = AEData_time_c(DelayIndexAE3(i),1)-InstronData_data(DelayIndexInst3(i),1);
    end
    
    timedelay_hist3(:,n) = timedelay3(n); % Keep a history of timedelays
    
    IndexShift3 = DelayIndexAE3(n)-fix(timedelay3(n).*AEaqrate);
    AEData_time_c(IndexShift3:DelayIndexAE3(n)-1,:)=[];
    for i=IndexShift3:length(AEData_time_c)
        AEData_time_c(i,1) = AEData_time_c(i,1) - timedelay3(n);
    end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    n=n+1;
     
    if n > min(min(length(DelayIndexInst),length(DelayIndexInst2)),length(DelayIndexInst3))
        disp('AE and Instron data are matched !')
        break
    end


end

% Delete extra points of AE data after Instron data recording is done
AEcut_4 = max(AEData_time_c(:,1)) - max(InstronData_data(:,1));
AE_endpoint = length(AEData_time_c) - fix(AEcut_4.*AEaqrate);
AEData_time_c(AE_endpoint:length(AEData_time_c),:) = [];

AEcut_5 = find(AEData_time_c(:,1) > 1,1); 
t=0;
for i=1:AEcut_5     % remove AE data that has negative time value
    i=i-t;
    if AEData_time_c(i,1) < 0
        AEData_time_c(i,:)=[];
        t=t+1;
    end
end

% Plot the result for comparing with original AE data
f3 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot(InstronData_data(:,1),InstronData_data(:,9),'b-',AEData_time_c(:,1),AEData_time_c(:,2),'r-')
title('Load vs. Time (Matched data)')
xlabel('Time (S)')
ylabel('Load (kN)')
legend('Instron Data','AE Data','Location','northeast','Orientation','horizontal')
print(f3,'Load vs. Time (Matched data).tiff','-dtiffn','-r300')

% add a column to AEData_time_c to write cycle numbers regarding each point 
% Cycles numbers are read from Instron data. According to AE and Instron
% data aqcuisition rates, 50/sec and 200/sec respectively, each AE point is
% matched with regarding Instron data point.

disp ('------------------------------------------------------------------')

ratio = length(InstronData_data)/length(AEData_time_c);
 
if (fix(ratio) == 20 || ceil(ratio) == 20)
    disp('Instron data points are 20 times AE data points, acccording to data acuisition rates: Instron = 1000 data/second, AE = 50 data/second.')
    for i=1:length(AEData_time_c)
        j=fix(ratio*i);
        AEcycle(i,:) = InstronData_data(j,3);
    end
elseif fix(length(InstronData_data)/length(AEData_time_c)) < 20
        disp('Error ! Too many AE datapoints are removed ! Instron data points must be 20 times AE data points, respecting their data acuisition rates: Instron = 200 data/second, AE = 50 data/second')
else 
        disp('Error ! Not enough AE datapoints are removed ! Instron data points must be 20 times AE data points, respecting their data acuisition rates: Instron = 200 data/second, AE = 50 data/second')
end

% Add cycles to AEData_time as a new column
AEData_time_c = [AEData_time_c AEcycle];






% Now generate a variable with the same size as AEData.data but instead of
% removing noise points, assign zero value to all columns
AEData_time_cz = [AEData_time AEData_load];

% Find Initial point to shift all data
AEcut2_1 = find(AEData_time_cz(:,2)>0.8,1);
IndexShift2_init = fix((AEData_time_cz(AEcut2_1,1) - InstronData_data(Instcut_1,1)).*AEaqrate);

% make the noise point at the beginning equal to zero
for i=1:IndexShift2_init
    AEData_time_cz(i,:) = 0;
end

% make the time shift
for i=(IndexShift2_init + 1):length(AEData_time_cz)
    AEData_time_cz(i,1) = AEData_time_cz(i,1) - timedelay_initial;
end

n=1;
while 1
    
    MaxLoadAE_z1 = find(AEData_time_cz(:,2)>LB); %Find index points in AE data
        
    q=1; % Keep initial points at UNDER 50% max loads
    for i=2:length(MaxLoadAE_z1)
        if AEData_time_cz(MaxLoadAE_z1(i),1) - AEData_time_cz(MaxLoadAE_z1(i-1),1) > 2
            DelayIndexAE_z1(q,:) = MaxLoadAE_z1(i);
            q=q+1;
        end
    end
              
    % Check if same number of indexes are taken from AE and Instron
    if length(DelayIndexAE_z1)==length(DelayIndexInst)
        disp([num2str(length(DelayIndexAE_z1)),' data points of AE and Instron data are used to match time indexes.'])
        disp(['Iteration #: ',num2str(n)])
    else
        disp('Error !! number of data points counted from AE and Instron (to match time index) does not match.')
        break
    end
    disp ('-----------------')
    
    % Find time delays in AE by comparing Delay Index times of AE and Instron
    for i=1:length(DelayIndexAE_z1)
      timedelay_z1(i,:) = AEData_time_cz(DelayIndexAE_z1(i),1)-InstronData_data(DelayIndexInst(i),1);
    end
    
    timedelay_hist_z1(:,n) = timedelay_z1; % Keep a history of timedelays
    
    % Delete redundant points recorded in AE data
    AEaqrate = 50; % AE data acquisition rate per second

    IndexShift_z1 = DelayIndexAE_z1(n)-fix(timedelay_z1(n).*AEaqrate);
    for i=IndexShift_z1:length(AEData_time_cz)
        AEData_time_cz(i,1) = AEData_time_cz(i,1) - timedelay_z1(n);
    end
    AEData_time_cz(IndexShift_z1:DelayIndexAE_z1(n)-1,:) = 0;
    
% % % % % % % % % %  Following part of the code removes the redundant time over 50% max load % % % % % % % % % % 
    
    %Find index points OVER 50% max load in AE data
    MaxLoadAE_z2 = find(AEData_time_cz(:,2)>LB2); 
        
    q=1; % Keep initial points at OVER 50% max loads
    for i=2:length(MaxLoadAE_z2)
        if AEData_time_cz(MaxLoadAE_z2(i),1) - AEData_time_cz(MaxLoadAE_z2(i-1),1) > 4
            DelayIndexAE_z2(q,:) = MaxLoadAE_z2(i);
            q=q+1;
        end
    end
        
    if length(DelayIndexAE_z2) ~= length(DelayIndexInst2)
        disp('Error !! number of data points counted from AE2 and Instron2 (to match time index) does not match.')
        break
    end
    
    % Find time delays in AE by comparing Delay Index times of AE and Instron
    for i=1:length(DelayIndexAE_z2)
      timedelay_z2(i,:) = AEData_time_cz(DelayIndexAE_z2(i),1)-InstronData_data(DelayIndexInst2(i),1);
    end
    
    timedelay_hist_z2(:,n) = timedelay_z2; % Keep a history of timedelays
    
    IndexShift_z2 = DelayIndexAE_z2(n)-fix(timedelay_z2(n).*AEaqrate);
    for i=IndexShift_z2:length(AEData_time_cz)
        AEData_time_cz(i,1) = AEData_time_cz(i,1) - timedelay_z2(n);
    end
    AEData_time_cz(IndexShift_z2:DelayIndexAE_z2(n)-1,:)= 0;
    
    %Find index points in AE data with load values less than delay index AE 2 load in AE data
    MaxLoadAE_z3 = find(AEData_time_cz(:,2) - AEData_time_cz(DelayIndexAE_z2(n),2) < 0);
    
    q=1; % Select the ones that happen after 1 second of delay index
    for i=1:length(MaxLoadAE_z3)
        if AEData_time_cz(MaxLoadAE_z3(i),1) - AEData_time_cz(DelayIndexAE_z2(n),1) > 1
            t(q,:) = MaxLoadAE_z3(i);
            q=q+1;
        end
    end
    
    % The first one would be Delay Index AE 3
    DelayIndexAE_z3(n,:) = t(1);
        
    for i=1:length(DelayIndexAE_z3)
      timedelay3(i,:) = AEData_time_cz(DelayIndexAE_z3(i),1)-InstronData_data(DelayIndexInst3(i),1);
    end
    
    timedelay_hist3(:,n) = timedelay3(n); % Keep a history of timedelays
    
    IndexShift_z3 = DelayIndexAE_z3(n)-fix(timedelay3(n).*AEaqrate);
    for i=IndexShift_z3:length(AEData_time_cz)
        AEData_time_cz(i,1) = AEData_time_cz(i,1) - timedelay3(n);
    end
    AEData_time_cz(IndexShift_z3:DelayIndexAE_z3(n)-1,:)= 0;
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    
    n=n+1;
     
    if n > min(min(length(DelayIndexInst),length(DelayIndexInst2)),length(DelayIndexInst3))
        disp('AE and Instron data are matched !')
        break
    end


end

% Assign zero to extra point at the end
IndexShift2_end = fix((max(AEData_time_cz(:,1))-max(InstronData_data(:,1))).*AEaqrate);
for i=(length(AEData_time_cz)-IndexShift2_end):length(AEData_time_cz)
    AEData_time_cz(i,:) = 0;
end

% For some reason 1 points are not changing to zero. It is done here.
% 1 points from the beginning turned to zero. IN CTAA06.
for i=IndexShift2_init+1:IndexShift2_init+((length(AEData_time) - length(AEData_time_c))-sum(AEData_time_cz(:,1)==0))
    AEData_time_cz(i,:) = 0;
end

% Make sure no negative time is within AEData_time2
for i=1:length(AEData_time_cz)
    if AEData_time_cz(i,1) < 0
        AEData_time_cz(i,:) = 0;
    end
end


% Check if the numbers of rows changed to zero is the same as rows removed
disp ('------------------------------------------------------------------')
if sum(AEData_time_cz(:,1)==0) == length(AEData_time) - length(AEData_time_c)
    disp('Successful ! Number of rows assinged zero value in Para_Hit.txt is the same as number of rows removed from Para_01.txt')
else
    disp('Error ! Number of rows changed to 0 in AEData_time2 variable is not the same as number of points removed from AEData.data')
end

% % Now assign cycle numbers from AEData_time to AEData_time2
j=1;
for i=1:length(AEData_time_cz)
    if AEData_time_cz(i,1) ~= 0
        AEcycle2(i,:) = AEData_time_c(j,3);
        j = j+1;
    elseif AEData_time_cz(i,1) == 0
        AEcycle2(i,:) = 0;
    end
end

% figure(n+2)
% plot(InstronData_data(:,1),InstronData_data(:,9),'b-', AEData_time_cz(:,1),AEData_time_cz(:,2),'ro')
% title('Load vs. Time')
% xlabel('Time (S)')
% ylabel('Load (kN)')
% legend('Instron Data','AE Data','Location','northeast','Orientation','horizontal')

% Add cycles to AEData_time2 as a new column 
AEData_time_cz = [AEData_time_cz AEcycle2];

% Check if length AEData_time2 and original Para file are the same
if length(AEData_time_cz) ~= length(AEData_time)
    disp ('------------------------------------------------------------------')
    disp ('Error ! Para_hit file does not have the same length as original Para file')
end

% Plot the result to reflect cycle match
f4 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot(InstronData_data(:,3),InstronData_data(:,9),'bo',AEData_time_c(:,3),AEData_time_c(:,2),'rx')
title('Load vs. Cycle (Matched data)')
xlabel('Cycle')
ylabel('Load (kN)')
legend('Instron Data','AE Data','Location','northeast','Orientation','horizontal')
print(f4,'Load vs. Cycle (Matched data).tiff','-dtiffn','-r300')

% Plot the result to reflect AE data in Para_hit file match Instron data
f5 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot(InstronData_data(:,1),InstronData_data(:,9),'b-',AEData_time_cz(:,1),AEData_time_cz(:,2),'rx')
title('Load vs. Time (Matched data - Para-hit file)')
xlabel('Time')
ylabel('Load (kN)')
legend('Instron Data','AE Data','Location','northeast','Orientation','horizontal')
print(f5,'Load vs. Time (Matched data - Para_hit file).tiff','-dtiffn','-r300')

% wirte the results into a new file.
% Generate a new Para file that contains cycle number and corrected time.
disp ('------------------------------------------------------------------')
disp ('Exporting results into a text file ...')
AEData_time_colheaders = [AEData.colheaders(2),AEData.colheaders(3),char('CYCLE')];
for i=1:length(AEData_time_colheaders)
    if i==length(AEData_time_colheaders)
        fprintf(fileID,'%s\n',char(AEData_time_colheaders(1,i)));
    else
        fprintf(fileID,'%s\t',char(AEData_time_colheaders(1,i)));
    end
end
for i=1:length(AEData_time_c)
    fprintf(fileID,'%1.4f         %2.4f         %u\r\n',AEData_time_c(i,:));
end
fclose(fileID);

% write another text file regarding AEData_time2 to use for hit data cycle
% assignment and time correction
for i=1:length(AEData_time_colheaders)
    if i==length(AEData_time_colheaders)
        fprintf(fileID2,'%s\n',char(AEData_time_colheaders(1,i)));
    else
        fprintf(fileID2,'%s\t',char(AEData_time_colheaders(1,i)));
    end
end
for i=1:length(AEData_time_cz)
    fprintf(fileID2,'%1.4f         %2.4f         %u\r\n',AEData_time_cz(i,:));
end
fclose(fileID2);

t_script = toc;

% Measure th overall time
t_total = t_AEInst + t_script;

% give a summary of the actions
disp ('------------------------------------------------------------------')
disp('Summary:')
disp([num2str(length(InstronData_data)),' data points from Instron are imported.'])
disp([num2str(length(AEData_time)),' raw data points from AE are imported.'])
disp([num2str(length(AEData_time) - length(AEData_time_c)),' data points are removed from raw AE data.'])
disp([num2str(length(AEData_time_c)),' data points of AE is matched with Instron data.'])
disp('AE data points are matched with Instron data points and results are exported.')
disp(['Total run time is ',num2str(t_total,'%.1f'),' seconds'])

