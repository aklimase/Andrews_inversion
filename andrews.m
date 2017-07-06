%% "secundo"

%% set paths

readdir='/Users/abaltay/Dropbox/Anza/Alexis_Data/cut_sac_files/';
upperdir='/Users/abaltay/Dropbox/Anza/Alexis_Data';
cd(readdir)
files=dir('*HNE.sac'); %I used these files because the 0.1 files had the wrong month/etc... 

for k=1:length(files)
   [data{k}, time{k}, hdrs{k}]=rsacnewbig(files(k).name);  
end

%% get metadata

for k=1:length(files)
   rsn(k)=k; %record sequence number
   evid(k)=hdrs{k}(79); %event id. In the sac format read by matlab, this info is stored in header variable 79
   sta{k}=char(hdrs{k}(111:115))'; %station name. not a great way to do it, you might read it from the file name. 
   
   % other info that might become useful: 
   mag(k)=hdrs{k}(40); 
   dist(k)=hdrs{k}(51); 
  
   stlat(k)=hdrs{k}(32); stlong(k)=hdrs{k}(33); stelev(k)=hdrs{k}(34); stdepth(k)=hdrs{k}(35); 
   evlat(k)=hdrs{k}(36); evlong(k)=hdrs{k}(37); evdepth(k)=hdrs{k}(39); 
   year(k)=hdrs{k}(71);  doy(k)=hdrs{k}(72); % not this is day-of-year, not day of month. 
   hour(k)=hdrs{k}(73); minute(k)=hdrs{k}(74);
   sec(k)=hdrs{k}(75) + hdrs{k}(76)/1000; 
   delta(k)=hdrs{k}(1); %sampling rate
   channel{k}=char(hdrs{k}(271:273))';
end

%% bookkeeping - create a numeric place holder for each event, station. 

EVID=unique(evid); I=length(EVID);   
STA=unique(sta); J=length(STA);

for i=1:length(evid)
   evid_num(i)=find(evid(i)==EVID); 
end
for j=1:length(sta)
    sta_num(j)=find(strcmp(sta{j}, STA));
end

EVID_NUM=unique(evid_num); %check: should be 1:I
STA_NUM=unique(sta_num); %check: should be 1:J

%% create G matrix
% dimensions: 
I=length(EVID); J=length(STA); K=length(rsn); 

%event matrix
G1=zeros(K,I); 
%site matrix
G2= zeros(K,J); 

for k=1:K
    G1(k,evid_num(k))=1; %put a 1 in the spot that corresponds to that event
    G2(k,sta_num(k))=1; %put a 1 in the spot that corresponds to that station
end

% big G matrix: 
G=[G1 G2]; %check size: should be Kx(I+J)

%% check G
sum(G) % <-- first 27 entries show how many stations record that event. 
% Most are recorded at 10, or 9 stations. What's up with the event only 
% recorded at 2 stations? 
% last 10 entries show the number of events recorded at each station
% they're all between 23 and 27, so that looks about right. 
tester=find(sum(G)<3)
files(evid_num==tester(1)).name
files(evid_num==tester(2)).name
files(evid_num==tester(3)).name

sum(G')
find(sum(G')~=2) %<- empty! 
% these are all 2, which means each record is one event + one site. 
%% set up data matricies 

% data matrix: 
d=zeros(K, F); %but I didn't define f yet... 
for f=1:F 
    d(:,f)=data{k}; %I don't have the spectral data... 
end

%% inversion
for f=1:F
    m(:,f)=inv(G'*G)*G'*d(:,f);  %L-2 solution
    m2(:,f)=pinv(G)*d(:,f); %pseudo inverse of the G matrix
    % just two ideas! 
end


