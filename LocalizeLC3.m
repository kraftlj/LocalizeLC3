%%% Automated analysis of LC3 localization
%%
%     LocalizeLC3 Copyright (C) 2014  Lewis J. Kraft
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Running the script launches the LC3 localization script.  It will
% prompt the user for the location of directory harboring the imaging
% files.  The script assumes there is a CFP channel with Cerulean as a
% control, the YFP channel is Venus-LC3, and the FarRed channel is a DRAQ5
% nuclear label.
%
%
%% This is a script to analyze the cell imaging data
%Inputs are images of cells, and images of nuclei.
%Outputs are:
% Number of nuclei;
% Number of spots per nucleus;
% Number of spots per cytoplasm;
% ratio of nuclear spot to nucleoplasm intensity
% ratio of cytoplasmic spot to cytoplasm intensity
% ratio of nuclear to cytoplasmic intensity
%% Begin
clear all; close all; clc;
tic %start the clock
location = uigetdir;
if exist(fullfile(location,'SaveFolder')) %If the SaveFolder for the analysis exists remove it
    rmdir(fullfile(location,'SaveFolder'),'s');
end
mkdir(fullfile(location,'SaveFolder')); %Make a new SaveFolder
%% Image importing
files = dir(fullfile(location,'*.lsm')); %This catalogs all of the .lsm files in the location directory
wb=waitbar(0,'Overall Time Left'); %This initializes a waitbar
header={'FileNames','CFPNCratio','YFPNCratio','#spots_cyt','AvgSpotSize_cyt','AvgSpotInt_cyt','ISS_cyt','#spots_nuc','AvgSpotSize_nuc','AvgSpotInt_nuc','ISS_nuc','SpotScatter'};
fid = fopen(fullfile(location,'Results.txt'),'w');
fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r\n', header{:});
% for index = 1:20;
    for index = 1:length(files); %This goes through the image files one by one
    wb=waitbar(index./length(files),wb); %This updates the waitbar each time it moves to the next image file
    
    clearvars('-except','location','files','fid','index','wb')
    data=bfopen(fullfile(location,files(index).name));
    CFP=data{1,1}{1};
    YFP=data{1,1}{2};
    FarRED=data{1,1}{3};
    
    %% Find cells and spots
    
    UEB = medfilt2(CFP, [10 10],'symmetric');
    BWYFP = im2bw(UEB, 200/2^16);
    BWYFP = imfill(BWYFP,'holes');
    D = -bwdist(~BWYFP);
    D(~BWYFP) = -Inf;
    D=imhmin(D,20);
    Labels = watershed(D);
%     BWYFP = logical(zeros(size(Labels)));
    BWYFP(Labels==0)=0;
%     BWYFP = imclearborder(BWYFP);
    %     SE = strel('disk',30);
    %BWYFP = imclose(BWYFP,SE);
    
% % % % % % % % % % % % % % % %     I2 = imtophat(YFP, strel('disk', 4));
% % % % % % % % % % % % % % % % BW = im2bw(I2,300/2^16);
% % % % % % % % % % % % % % % % BWspots = imopen(BW,strel('disk',2));
% % % % % % % % % % % % % % % % imagesc(BWspots)
    
 
    BWYFP = bwareaopen(BWYFP,5000);
    %BWYFP = imclearborder(BWYFP);
    UEB = medfilt2(YFP, [10 10],'symmetric');
    spots = YFP-UEB; %This subtracts the uneven cellular background from the spots
    BWspots = im2bw(spots, 150/2^16); %This leaves the spots with value 1 and background with value 0.
    %BWspots = imclearborder(BWspots);
    
    %% Find the nuclei of transfected cells
    UEB = medfilt2(FarRED, [10 10],'symmetric');
    BWnuc = im2bw(UEB, 400/2^16);
    D = -bwdist(~BWnuc);
    D(~BWnuc) = -Inf;
    D=imhmin(D,20);
    Labels = watershed(D);
    BWnuc(Labels==0)=0;
    
    BWnuc(~BWYFP)=0;
    BWnuc=imfill(BWnuc,'holes');
%     BWnuc = imclearborder(BWnuc);
    %     BWnuc = imclose(BWnuc,SE);
    BWnuc = bwareaopen(BWnuc,4000);
   
    
    %% Refine the spots
    
    SE = strel('disk',2);
    BWspots = imopen(BWspots,SE);
    SE = strel('disk',3);
    BWspots = imdilate(BWspots,SE);
    BWspots(~BWYFP)=0;
    BWnucspots = BWspots;
    BWnucspots(~BWnuc) = 0;
    BWcytspots = BWspots;
    BWcytspots(BWnuc) = 0;
    %     D = -bwdist(~BWspots);
    %     D(~BWspots) = -Inf;
    %     D=imhmin(D,0);
    %     SpotLabels = watershed(D);
    %     BWspots = logical(zeros(size(SpotLabels)));
    %     BWspots(SpotLabels>1)=1;
    CCspots = bwconncomp(BWspots,4);
    
    CCnucspots = bwconncomp(BWnucspots,4);
    CCcytspots = bwconncomp(BWcytspots,4);
    
    %% Refine cytoplasmic and nucleoplasmic compartments
    BWYFPcyt = BWYFP;
    BWYFPcyt(BWnuc) = 0;
    BWYFPcyt(BWspots) = 0;
    
    BWYFPnuc = BWnuc;
    BWYFPnuc(BWspots) = 0;
    CCnuc = bwconncomp(BWnuc,4);
    %% Record statistics
    % calculate NC ratios for YFP and CFP images
    YFPNCratio = mean(YFP(BWYFPnuc))./mean(YFP(BWYFPcyt));
    CFPNCratio = mean(CFP(BWYFPnuc))./mean(CFP(BWYFPcyt));
    % calculate stats for nuclear spots
    STATS = regionprops(CCnucspots, YFP, 'MeanIntensity','Area','WeightedCentroid');
    AvgSpotSize_nuc = mean([STATS.Area]);
    AvgSpotInt_nuc = mean([STATS.MeanIntensity]);
    Number_Spots_nuc = CCnucspots.NumObjects./CCnuc.NumObjects;
    ISS_nuc = AvgSpotSize_nuc*AvgSpotInt_nuc;
    % calculate stats for cytoplasmic spots
    STATS = regionprops(CCcytspots, YFP, 'MeanIntensity','Area','WeightedCentroid');
    AvgSpotSize_cyt = mean([STATS.Area]);
    AvgSpotInt_cyt = mean([STATS.MeanIntensity]);
    Number_Spots_cyt = CCcytspots.NumObjects./CCnuc.NumObjects;
    ISS_cyt = AvgSpotSize_cyt*AvgSpotInt_cyt;
    % calculate spot scatter
    STATSspots = regionprops(CCspots, YFP, 'WeightedCentroid');
    if isempty(STATSspots) | length(STATSspots)<=2
        SpotScatter = NaN;
    else  
    centroids = cat(1, STATSspots.WeightedCentroid);
    avgcentroid = mean(centroids,1);
    SpotScatter = std(sqrt((centroids(:,1)-avgcentroid(1)).^2+(centroids(:,2)-avgcentroid(2)).^2));
    end
    %% Save the image analysis results
    rgb(:,:,1)=uint8(BWspots).*255;
    rgb(:,:,2)=uint8(BWYFPcyt).*255;
    rgb(:,:,3)=uint8(BWYFPnuc).*255;
    imwrite(rgb,fullfile(location,'SaveFolder',[files(index).name(1:end-4),'.tif']));
    savedata={files(index).name,CFPNCratio,YFPNCratio,Number_Spots_cyt,AvgSpotSize_cyt,AvgSpotInt_cyt,ISS_cyt,Number_Spots_nuc,AvgSpotSize_nuc,AvgSpotInt_nuc,ISS_nuc,SpotScatter};
    fprintf(fid, '%s\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\r\n', savedata{:});
    end

fclose(fid);
close(wb);

%% End
toc