%%

cd /Users/abaltay/Documents/Global_Eqs/Tohoku/Data/GSN/sod_Swave/allSHorizSpectra
filename='COR.2003_10_31_01_06_32.00.BHT.sac1am'; %a spectral file

[data,freqs,header]=rsacnewlil(filename);
 

%%
fstart=f(2); %second frequency point; 
fend=f(end); %last frequency point. 
% could also define these to be 0.02 to 50 or something. 

num_pts=50; %define number of spectral points. 

c=logspace(log10(fstart), log10(fend), num_pts); 

% at each c point, go there and take a mean out to each side bin
l=length(c); 
newdata(1)=mean(d(1:find(f<=c(2),1,'last'))); %define the first point 
newdata(l)=mean(d(find(f>=c(end-1),1,'first'):end)); %define the last point
for j=2:l-1 % define the 2nd to l-2 points. 
    lower(j)=find(f>=c(j-1),1,'first');
    upper(j)=find(f<=c(j+1),1,'last');
    newdata(j)=mean(d(lower(j):upper(j)));
end

%% make a plot
figure
loglog(f,d)
hold on
plot(c, newdata, 'ok-')

for j=2:l-1 % define the 
    plot(f([lower(j) upper(j)]), [newdata(j) newdata(j)], 'r')
end

legend('raw spectrum', 'smoothed spectrum', 'smoothing bin widths')
xlabel('frequency [Hz]')
ylabel('displacement spectrum')