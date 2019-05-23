% This code is designed to apply load and Delta T filter on hits and record
% valid hit numbers in a file
clc; clear all;

% This code is applicable to result of crack growth tests: CTAA04-CTAA08
disp(datetime('now'));

% Specify all file locations here for convenience. 
AEfilename_root = '~\AE\';

% Specify number of Para data files for each test.
para_num = 8; 
              
% Define maximum load to be used in hit load filtering
maxload = 3.5472;

% Define Delta T value to be used in hit Delta T filtering
deltaT = 10e-6; 

% keep hits that passed load and Delta T filters 
fileID1 = fopen([AEfilename_root,'Results\Ch1_load_deltaT_filtered.txt'],'w');
fileID2 = fopen([AEfilename_root,'Results\Ch2_load_deltaT_filtered.txt'],'w'); 

% Set the limits for figures axises
limit_load = [0 inf 0 maxload*1.1];

% Importing AE data
disp('Importing AE data ...')

% Measure the time for importing data
tic;

% Import Para_growth_## files
para_hit_filename = [AEfilename_root,'Results\Para_hit.txt'];
Parahit = importdata(para_hit_filename,' ',1);
Parahit_data = Parahit.data;

% Import Para_growth_## files
for i=1:para_num
    
    if i < 10
        para_filename = ['Para_growth_0',num2str(i),'.TXT'];
    else
        para_filename = ['Para_growth_',num2str(i),'.TXT'];
    end
    
    AEfilename = [AEfilename_root,para_filename];
    HitData = importdata(AEfilename,' ',8);
       
    Para_t = HitData.data(:,2);
    Para_l = HitData.data(:,3);
    
    if i==1
        ParaData_time = Para_t;
        ParaData_load = Para_l;
    else
        ParaData_time = [ParaData_time; Para_t];
        ParaData_load = [ParaData_load; Para_l];
    end
    
end

% Multiply the load values by the multipication factor
ParaData_load = ParaData_load*10; 

% Generate one variable including all Para data information [time load]
ParaData_data = [ParaData_time ParaData_load];

% Import Test_growth_## files
for i=1:para_num
    
    if i < 10
        hit_filename = ['Test_growth_0',num2str(i),'.TXT'];
    else
        hit_filename = ['Test_growth_',num2str(i),'.TXT'];
    end
    
    AEfilename = [AEfilename_root,hit_filename];
    HitData = importdata(AEfilename,' ',8);
    
    HitData_temp = HitData.data;
    
    if i==1
        HitData_data = HitData_temp;
    else
        HitData_data = [HitData_data; HitData_temp];
    end
    
end

HitData_data(:,3) = HitData_data(:,3)*10;

% Import Sum_growth_## files
for i=1:para_num
    
    if i < 10
        sum_filename = ['Sum_growth_0',num2str(i),'.TXT'];
    else
        sum_filename = ['Sum_growth_',num2str(i),'.TXT'];
    end
    
    AEfilename = [AEfilename_root,sum_filename];
    SumData = importdata(AEfilename,' ',4);
    
    SumData_temp = SumData.data;
    
    if i==1
        SumData_data = SumData_temp;
    else
        SumData_data = [SumData_data(1,1)   SumData_data(1,2)+SumData_temp(1,2); SumData_data(2,1)   SumData_data(2,2)+SumData_temp(2,2)];
    end
    
end

% Measure data import time
t_AE = toc;
tic;

% Display a message to know when data is imported
disp('AE data is imported !');
disp(['Elapsed time to import data is ',num2str(t_AE,'%.1f'),' seconds'])
disp ('------------------------------------------------------------------')

% Display column headers 
for o = 1:length(HitData.colheaders)
    colmn = ['Hit dataset column ',num2str(o), ' header',' is ',char(HitData.colheaders(o)),'.'];
    disp(colmn)
end
disp ('------------------------------------------------------------------')

for k = 1:length(Parahit.textdata)
    colmn = ['Para_hit dataset column ',num2str(k), ' header',' is ',char(Parahit.textdata(k)),'.'];
    disp(colmn)
end
disp ('------------------------------------------------------------------')

%Count number of hits from each channel
ch1cnt = sum(HitData_data(:,5)==1);
ch2cnt = sum(HitData_data(:,5)==2);

% Check if count is done correctly
if ch1cnt == SumData_data(1,2)
    disp('Hits in channel 1 are matching AEWin summary report.')
elseif SumData_data(1,2) - ch1cnt > 0
    disp(['Caution ! ',num2str(SumData_data(1,2) - ch1cnt),' Hits in channel 1 have more than one waveform !'])
elseif SumData_data(1,2) - ch1cnt < 0
    disp('Warning ! Hits counted in channel 1 are NOT matching AEWin summary report !')
    disp([num2str(ch1cnt-SumData_data(1,2)),' hits are missing !'])
end

if ch2cnt == SumData_data(2,2)
    disp('Hits in channel 2 are matching AEWin summary report.')
elseif SumData_data(2,2) - ch2cnt > 0
    disp(['Caution ! ',num2str(SumData_data(2,2) - ch2cnt),' Hits in channel 2 have more than one waveform !'])
elseif SumData_data(2,2) - ch2cnt < 0
    disp('Warning ! Hits counted in channel 2 are NOT matching AEWin summary report !')
    disp([num2str(ch2cnt-SumData_data(2,2)),' hits are missing !'])
end
disp ('------------------------------------------------------------------')

%Prepare matrices to import data regarding each channel
ch1data = zeros(ch1cnt,length(HitData.colheaders));
ch2data = zeros(ch2cnt,length(HitData.colheaders));

%Divide hit data for each channel
i=1;
j=1;
for k = 1:length(HitData_data(:,5))
    if HitData_data(k,5) == 1     
        ch1data(i,:)=HitData_data(k,:);
        i = i+1;
    else
        ch2data(j,:)=HitData_data(k,:);
        j = j+1;
    end
end

% For each channel data, a new column is added that contains the number of
% each data in that channel. This number corresponds the filename number in
% waveform folder. It is useful since removed data will cause miss match in
% order for reading waveforms.
for k = 1:length(ch1data)
    ch1order(k,:) = k;
end

for k = 1:length(ch2data)
    ch2order(k,:) = k;
end

% Add order number to the channel data
ch1data = [ch1data ch1order];
ch2data = [ch2data ch2order];
disp(['Column ',num2str(o+1),' is added to dataset that contains hit orders.']);
disp ('------------------------------------------------------------------')

disp ('Applying load and Delta T filter on hits ...')

%Apply load filter data - Any hit happened in loads less than 20% of the
%peak load is removed - Hit number removed is stored in ch1lfr (or
%ch2lfr) 
% Record load filtered data for channel 1 in ch1lf.
ch1lf = ch1data;
i = 0;
for k = 1:ch1cnt
    k = k - i;
    if le(ch1lf(k,3),0.8*maxload)
        ch1lf(k,:)=[];
        i = i+1;
    end
end

% Record load filtered data for channel 2 in ch2lf.
ch2lf = ch2data;
j = 0;
for k = 1:ch2cnt
    k = k - j;
    if le(ch2lf(k,3),0.8*maxload)
        ch2lf(k,:)=[];
        j = j+1;
    end
end

disp('Load filter applied, now applying Delta T filter ...')

%Find which hit numbers of channel 2 should remain in dataset in variable
%ch2kd. 
for k = 1:length(ch2lf)
    for i = 1:length(ch1lf)
        if le(abs(ch2lf(k,2)-ch1lf(i,2)),deltaT) == 1
            ch2kd(k,:) = ch2lf(k,18);
        end
    end
end

%Clear variable ch2kd by removing zeros
j=0;
for k = 1:length(ch2kd)
    k = k-j;
    if ch2kd(k) == 0
            ch2kd(k,:) = [];
            j = j+1;
    end
end 

%Create variable ch2dtf that only holds data for remained hits based on
%Delta T filter method
i=1;
for k=1:length(ch2lf)
    for m = 1:length(ch2kd)
        if (ch2lf(k,18) - ch2kd(m)) == 0
            ch2dtf(i,:) = ch2lf(k,:);
            i = i+1;
        end
    end        
end

%Do the same steps for data from channel 1. Find which hit numbers of 
%channel 1 should remain in dataset in variable ch1kd. 
for k = 1:length(ch1lf)
    for HitLocCh1 = 1:length(ch2dtf)
        if le(abs(ch1lf(k,2)-ch2dtf(HitLocCh1,2)),deltaT) == 1
            ch1kd(k,:) = ch1lf(k,18);
        end
    end
end

%Clear variable ch1kd by removing zeros
j=0;
for k = 1:length(ch1kd)
    k=k-j;
    if ch1kd(k) == 0
            ch1kd(k) = [];
            j = j+1;
    end
end 

%Create variable ch1dtf that only holds data for remained hits based on
%Delta T filter method
i=1;
for k=1:length(ch1lf)
    for m = 1:length(ch1kd)
        if (ch1lf(k,18) - ch1kd(m)) == 0
            ch1dtf(i,:) = ch1lf(k,:);
            i = i+1;
        end
    end        
end

% Check to see if number length ch1dtf and ch2dtf are the same
if length(ch1dtf) ~= length(ch2dtf)
    disp('Error ! Delta T filtered hits from Channle 1 and Channel two do not have same length !')
    disp ('------------------------------------------------------------------')
else
    disp('Load and Delta T filter is successfully applied on data.')
    disp ('------------------------------------------------------------------')
end

% Lets add cycle and corrected time from Para_hit to the data variables.
% Para files record 50 points per second, therefore have 0.02 seconds
% difference between points. Hit data on the other hand records exact hit
% time that most likely happens in these 0.02 intervals. To add cycle and
% corrected time, first we have to look into original Para file based on
% time, then use the load and extension of original file to look at
% Para_hit file and add corrected time and cycle to ch data files from
% there !
% First, find hit location within original Para file

disp('Finding Channel 1 hits corresponding cycle number and correct time value ...')

% Channel 1
HitLocCh1 = zeros(length(ch1dtf),1);
for i=1:length(ch1dtf)
    HitLocCh1(i,:) = find(abs(ch1dtf(i,2)-ParaData_data(:,1))<0.02,1);
end

HitCycleCh1 = zeros(length(HitLocCh1),1);
for i=1:length(HitLocCh1) % Find cycles
    HitCycleCh1(i,:) = Parahit_data(HitLocCh1(i),3);
    t=1;
    while (HitCycleCh1(i,:) == 0 && HitLocCh1(i)+t<length(Parahit_data))
        HitCycleCh1(i,:) = Parahit_data(HitLocCh1(i)+t,3);
        t=t+1;
    end 
end

HitTimeCh1 = zeros(length(HitLocCh1),1);
for i=1:length(HitLocCh1) % Find correcte time
    HitTimeCh1(i,:) = Parahit_data(HitLocCh1(i),1);
    t=1;
    while  (HitTimeCh1(i,:) == 0 && HitLocCh1(i)+t<length(Parahit_data))
        HitTimeCh1(i,:) = Parahit_data(HitLocCh1(i)+t,1);
        t=t+1;
    end 
end

% Add cycle and corrected time to the ch1data
ch1dtf = [ch1dtf HitCycleCh1 HitTimeCh1];

disp(['Column ',num2str(o+2),' is added to Cahnnel 1 dataset that contains cycle number.']);
disp(['Column ',num2str(o+3),' is added to Cahnnel 1 dataset that contains corrected time values.']);
disp ('------------------------------------------------------------------')

disp('Finding Channel 2 hits corresponding cycle number and correct time value ...')

% Channel 2
HitLocCh2 = zeros(length(ch2dtf),1);
for i=1:length(ch2dtf)
    HitLocCh2(i,:) = find(abs(ch2dtf(i,2)-ParaData_data(:,1))<0.02,1);
end

HitCycleCh2 = zeros(length(HitLocCh2),1);
for i=1:length(HitLocCh2) % Find cycles
    HitCycleCh2(i,:) = Parahit_data(HitLocCh2(i),3);
    t=1;
    while HitCycleCh2(i,:) == 0 && HitLocCh2(i)+t<length(Parahit_data)
        HitCycleCh2(i,:) = Parahit_data(HitLocCh2(i)+t,3);
        t=t+1;
    end        
end

HitTimeCh2 = zeros(length(HitLocCh2),1);
for i=1:length(HitLocCh2) % Find correcte time
    HitTimeCh2(i,:) = Parahit.data(HitLocCh2(i),1);
    t=1;
    while HitTimeCh2(i,:) == 0 && HitLocCh2(i)+t<length(Parahit_data)
        HitTimeCh2(i,:) = Parahit_data(HitLocCh2(i)+t,1);
        t=t+1;
    end
end

% Add cycle and corrected time to the ch1data
ch2dtf = [ch2dtf HitCycleCh2 HitTimeCh2];

disp(['Column ',num2str(o+2),' is added to Cahnnel 2 dataset that contains cycle number.']);
disp(['Column ',num2str(o+3),' is added to Cahnnel 2 dataset that contains corrected time values.']);
disp ('------------------------------------------------------------------')

% %Plot raw hit data for each channel
f1 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot (ch1data(:,2),ch1data(:,3),'rx',ch2data(:,2),ch2data(:,3),'bo')
title ('Raw AE hit data (raw time values)')
xlabel('Time (S)')
ylabel('Load (KN)')
axis(limit_load)
legend('Channel 1','Channel 2','Location','southeast','Orientation','vertical')
movegui(f1,'northwest')
print(f1,'01_Raw AE hit data.tiff','-dtiffn','-r300')

% % Plot load and Delta T filtered hit data
f3 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot (ch1dtf(:,20),ch1dtf(:,3),'rx',ch2dtf(:,20),ch2dtf(:,3),'bo')
title ('Load and Delta T filtered AE hit data (AE Instron matched time)')
xlabel('Time (S)')
ylabel('Load (KN)')
axis(limit_load)
legend('Channel 1','Channel 2','Location','southeast','Orientation','vertical')
movegui(f3,'northeast')
print(f3,'03_Load and Delta T filtered AE hit data.tiff','-dtiffn','-r300')

t_data_filter = toc;

if length(ch1dtf) ~= length(ch2dtf)
    disp('Error ! Number of Delta T filtered hits are not the same in channel 1 and 2 !')
    disp ('------------------------------------------------------------------')
end
disp(['Elapsed time to filter hits and assign cycle number and correct time value is ',num2str(t_data_filter,'%.1f'),' seconds.'])
disp ('------------------------------------------------------------------')

disp('Exporting Load and Delta T filtered data into a text file ...')
ch1dtf_colheader = [HitData.colheaders,'ORDER','CYCLE','TIME'];
for i=1:length(ch1dtf_colheader)
    if i==length(ch1dtf_colheader)
        fprintf(fileID1,'%s\r\n',char(ch1dtf_colheader(1,i)));
    else
        fprintf(fileID1,'%s\t',char(ch1dtf_colheader(1,i)));
    end
end
for i=1:length(ch1dtf)
    fprintf(fileID1,'%u\t %f\t %1.4f\t %1.4f\t %u\t %u\t %u\t %u\t %u\t %u\t %u\t %1.4f\t %u\t %u\t %1.3f\t %u\t %u\t %u\t %u\t %1.4f\r\n',ch1dtf(i,:));
end
fclose(fileID1);

ch2dtf_colheader = [HitData.colheaders,'ORDER','CYCLE','TIME'];
for i=1:length(ch2dtf_colheader)
    if i==length(ch2dtf_colheader)
        fprintf(fileID2,'%s\r\n',char(ch2dtf_colheader(1,i)));
    else
        fprintf(fileID2,'%s\t',char(ch2dtf_colheader(1,i)));
    end
end
for i=1:length(ch2dtf)
    fprintf(fileID2,'%u\t %f\t %1.4f\t %1.4f\t %u\t %u\t %u\t %u\t %u\t %u\t %u\t %1.4f\t %u\t %u\t %1.3f\t %u\t %u\t %u\t %u\t %1.4f\r\n',ch2dtf(i,:));
end
fclose(fileID2);

disp('Done !')
disp('Load and Delta T filtered data from Channel 1 and Channel 2 are exported to text file and ready to be used for AE entropy analysis.')
toc
