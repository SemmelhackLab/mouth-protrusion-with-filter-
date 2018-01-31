%this is for all files in the folder(avi), and automatic cropped region
%Output:
% binary image of the first cropped frame for each video
% narea: area ratio of mouth
% avgnarea: filtered area ratio of mouth
% v2: slope of avgnarea

%trouble shooting:
%1. paramecia merging with the mouth may cause error in mouth protrusion
%(area) measurements
%Solution: set a threshold to ignore those data e.g. > 1.2


%template
% obj = VideoReader('D:\Julie\Dec 19th\4th_03_s__12_19_11_58_11_876.avi');
% vid = read(obj);
% frames = obj.NumberOfFrames;
% for x = 1 : frames
% imwrite(vid(:,:,:,x),strcat('frame-',num2str(x),'.tif'));
% end
myFolder = 'E:\DT strike candidates\2018 Jan 11th';%select imput folder
subFolderName = 'D:';%select output folder
%if subFolderName is changed, please also change newChr below
filePattern = fullfile(myFolder, '*.avi'); %select type of files
aviFiles   = dir(filePattern);
w = length(aviFiles);

f = 100; % frames used for normalization (mouth in normal size)
l = 1; % threshold of intensity, change this for 0.7 if the mouth and the line merge together
% r = [281 359 30 13];
bz = 27; %height of roi
nc = 13; %averaging period
%r = [264 309 27 13];% area of investigation, [xmin ymin width height] %[232 315 27 15] [60 206 23 31]
% % c = [307	358
% % 423	
% % 
% % ]; 
% c =  xlsread('D:\Dec 22\annotated dec 22_14 good.xlsx');% target frame with strike behavior observed
% [q u] = size(c);

% n = zeros(f,w);
% narea = zeros(p,w);% narea: normalized area of each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for ij = 5:w
         baseFileName = aviFiles(ij).name;
  fullFileName = fullfile(myFolder, baseFileName);
  %fprintf('Now reading %s\n', fullFileName);
  %imageArray = imread(fullFileName);
  
obj{ij} = VideoReader(fullFileName);
vid{ij} = read(obj{ij});
p = obj{ij}.NumberOfFrames;%total number of frames






%testing if the threshold of intensity is suitable
%test 1: rule out zero matrix
%test 2: rule out extreme ratio of mouth vs background
%test 3: rule out more than one object representing either the mouth or background
%A = imread('C:\Users\admin\Documents\UST 1718 spring\Julie Code\example.tif',1);

 r = roi(vid{ij}(:,:,:,1), bz); %detecting roi 
A2 = imcrop(vid{ij}(:,:,:,1),r);
BW = im2bw(A2, l);

    WB = ~BW;
    
    while nnz(BW) ==0
        l = l-0.05;
    BW = im2bw(A2, l);
    end
    
     if l <0.6
    l = 0.75;
     end   

while nnz(WB)/nnz(BW) >10 | nnz(BW)/nnz(WB) >10 
    l = l-0.05;
    BW = im2bw(A2, l);
    
    WB = ~BW;
    
end

 if l <0.6
    l = 0.75;
     end 

[B,L,N,a] = bwboundaries(BW);
[B2,L2,N2,a2] = bwboundaries(WB);
while N >2 | N2 >2 
    l = l-0.05;
    BW = im2bw(A2, l);
    [B,L,N,a] = bwboundaries(BW);
    WB = ~BW;
    [B2,L2,N2,a2] = bwboundaries(WB);

end

if l <0.8 & nnz(WB)/nnz(BW)< 0.2
    
    l = 0.8;
end

figure;
A2 = imcrop(vid{ij}(:,:,:,1),r);
BW = im2bw(A2, l);
imshow(BW);

% finding normal mouth
for i = 1:p


% A2 = imcrop(vid{ij}(:,:,:,i),r);
% BW = im2bw(A2, l);
% BW2 = bwareaopen(BW,120,26);


BW = im2bw(vid{ij}(:,:,:,i), l);
BW = ~BW;
BW2 = bwareaopen(BW,200,26);
A2 = imcrop(BW2,r);


%for testing reduction of noise
% figure;
% imshow(A2);

n{ij}(i) = bwarea(A2);

end
m = mean(n{ij}(:));
% finding mouth protrusion 
for i = 1:p


narea{ij}(i) = n{ij}(i)/m;

end


% 1. plot area ratio
figure;
hl = 1:p;
hl2 = (hl-1)/300;
plot(hl2,narea{ij}(:));
axis([0 p/300 -inf inf]);

%1b. apply filter
figure;
coeffnc = ones(1, nc)/nc;
avgnarea{ij} = filter(coeffnc, 1, narea{ij});
%plot(hl2,[narea{ij} avgnarea]);
plot(hl2,avgnarea{ij});
axis([0 p/300 0.8 inf]);


% find area ratio >= 1.1

na2{ij} = avgnarea{ij}(find( avgnarea{ij} >= 1.1)); %selected peaks, highere than 1.1
na2x{ij} =hl2(find( avgnarea{ij} >= 1.1));




% %example of moving average filter

% hoursPerDay = 24;
% coeff24hMA = ones(1, hoursPerDay)/hoursPerDay;


% avg24hTempC = filter(coeff24hMA, 1, tempC);
% plot(days,[tempC avg24hTempC])
% legend('Hourly Temp','24 Hour Average (delayed)','location','best')
% ylabel('Temp (\circC)')
% xlabel('Time elapsed from Jan 1, 2011 (days)')
% title('Logan Airport Dry Bulb Temperature (source: NOAA)')


% % % for annotation of peaks, observed by eyes only
% % hold on;
% % cc = c(ij)/300;
% % plot(cc, narea{ij}(c(ij)),'r*');
% hold on;
% rf = 1;
% for u2 = 1:u
%     isnan(c(ij,u2));
%     if isnan(c(ij,u2)) == 1
%     
%     else
% cc = c(ij)/300;
% plot(cc, narea{ij}(c(ij)),'r*');
% hold on;
%     end
% hold off;
% end

% % for annotation of peaks, observed by eyes only
% hold on;
% cc = c(ij)/300;
% plot(cc, narea{ij}(c(ij)),'r*');
hold on;
plot(na2x{ij}, na2{ij},'r*');

hold off;


title(baseFileName);
xlabel('Time(s)');
ylabel('Area ratio');



narea{ij} = narea{ij}';
n{ij} = n{ij}';

figure;
%FX = gradient(F);
v{ij} = gradient(narea{ij});
plot(hl2, v{ij});
title(baseFileName);
xlabel('Time(s)');
ylabel('rate of change of Area ratio');
% find rate of change of area ratio >= 0.015

v2{ij} = v{ij}(find( v{ij} >= 0.015 & v{ij} <= 0.035)); %selected peaks, highere than 1.1
v2x{ij} =hl2(find( v{ij} >= 0.015 & v{ij} <= 0.035));
hold on;
plot(v2x{ij}, v2{ij},'r*');

hold off;

% saveas(figure(2*ij),[pwd '/subFolderName/baseFileName.fig']);
 newChr = strrep(subFolderName,':',':\c12345');
newChr2 = strrep(newChr,'c12345',baseFileName);

newChr3 = strrep(newChr2,'.avi','_binary.fig');
newChr4 = strrep(newChr2,'.avi','_ratio.fig');
newChr5 = strrep(newChr2,'.avi','_ratio filter.fig');
newChr6 = strrep(newChr2,'.avi','_slope filter.fig');

saveas(figure(4*ij-3),newChr3);
saveas(figure(4*ij-2),newChr4);
saveas(figure(4*ij-1),newChr5);
saveas(figure(4*ij),newChr6);


% 
% newChr = strrep(subFolderName,'rame','rame\c');
% newChr = strrep(subFolderName,'c',baseFileName);


% save(gca,'c:\figures\test','jpeg');

end