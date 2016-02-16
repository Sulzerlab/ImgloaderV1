function [ thrIm2 ] = thresholding( thrIm, punctamindiameter, punctamaxdiameter, varargin )

global demo
% Length units are entered in pixels

ws = cell2mat(varargin(find(strcmp([varargin],'ws'))));
demo = cell2mat(varargin(find(strcmp([varargin],'demo')))); % Run program in demo mode

% Demo mode command %
if isempty(demo) == 0

    h = figure();
    imshow(thrIm); title('Pick Area to Threshold');
    demoCoords = fix(getrect());
    thrIm = thrIm(demoCoords(2):demoCoords(2)+demoCoords(4),demoCoords(1):demoCoords(1)+demoCoords(3));
    delete(h);

end
% % % % % % % % % % %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Setting the parameters of thresholding        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin(find(strcmp([varargin],'parameters')))) == 1

    punctaMinArea = round(punctamindiameter^2); % Pixels not length
    punctaMaxArea = round(punctamaxdiameter^2); % Pixels not length

    maskMinArea = 5*punctaMinArea;
    maskMaxArea = 20*punctaMaxArea;
    growth = 1;

    threshDivs = 50;
    
    [dummy,bgEstim] = max(histc(flatmat(thrIm),[0:1/threshDivs:1]));
    
    threshvalues=[0:1/threshDivs:max(max(thrIm))];

    punctaMinProm = 1;
    closeParam = 3;

    parameters = {'threshDivs', threshDivs,...
                  'punctaMinProm', punctaMinProm,...
                  'threshvalues', threshvalues,...
                  'maskMaxSize', maskMaxArea,...
                  'punctaMinArea', punctaMinArea,...
                  'punctaMaxArea', punctaMaxArea,...
                  'closeParam', closeParam};
else
    
    parameters = varargin(find(strcmp([varargin],'parameters')+1));

end


[promIm, keepMasks] = multithreshFn(thrIm, parameters);

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Watershed of the threshold image                        %
%                                                         %
% Inputs: promIm (contour image of puncta)                %
%         keepMasks (masks of puncta reaching criteria)   %
%         punctaMinArea                                   %
%                                                         % 
% Outputs: keepMasks (puncta masks meeting criteria)      %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(ws) == 0

    thrIm2 = bwareaopen(promIm.*keepMasks,punctaMinArea); % To keep puncta
    
    medDist = floor(growth*double(median(median(bwdist(thrIm2)))));
    
    
    bgIm = imdilate(thrIm2,strel('disk',medDist));
    bgIm = imfill(bgIm,'holes');
    
    thrWsImMin = imimposemin(imcomplement(thrIm), ~bgIm | thrIm2); 
    thrWsImMin = imclose(thrWsImMin,strel('square',closeParam));
    
    thrWsImOL = watershed(thrWsImMin);% Actually creates a label matrix
    
    labelAreas = histc(flatmat(thrWsImOL),[1:1:max(max(thrWsImOL))]);
    labelOutliers = find(labelAreas>maskMaxArea | labelAreas<maskMinArea);

    %figure; subplot(2,1,1); bar(labelAreas); hold on; plot(labelOutliers,1.1*labelAreas(labelOutliers),'rp');
    
    labelAreas(labelOutliers) = [];
    %subplot(2,1,2); hist(labelAreas);
    
    thrWsIm = thrWsImOL;
    
    for i = 1:numel(labelOutliers)
        
        thrWsIm(thrWsIm==labelOutliers(i)) = 0;
        
    end
    
    thrWsIm(thrWsIm>0)=1;
    
    [B,L] = bwboundaries(thrWsIm,'noholes');
  
    thrIm2 = zeros(size(thrIm,1),size(thrIm,2));
    
        for j = 1:size(B,1)
            
            if min(B{j}(:,2)) > 1 & min(B{j}(:,1)) > 1 & max(B{j}(:,1)) < size(thrIm,1)-1 & max(B{j}(:,2)) < size(thrIm,2)-1 % Stay away from edge pixels
            
            thrIm2 = thrIm2 + poly2mask(B{j}(:,2),B{j}(:,1),size(thrIm2,1),size(thrIm2,2));
            
            end
        
        end
    
else
    
    thrIm2 = keepMasks;
    
end

        % Demo mode commands %
        if isempty(demo) == 0


            figure; subplot(1,3,1); imagesc(thrWsImMin); title('Minima Imposed Image')
                    subplot(1,3,2); imagesc(thrWsImOL); title('Watershed and Filtered on Size')
                    subplot(1,3,3); imagesc(thrIm2); title('Turned into Masks')

        end
        % % % % % % % % % % %
        
        
end

function [promIm, keepMasks] = multithreshFn(thrIm, parameters)

global demo

promIm = zeros(size(thrIm));

try
    
    threshDivs = cell2mat(parameters(circshift(strcmp('threshDivs',parameters),[0 1])));
    punctaMinProm = cell2mat(parameters(circshift(strcmp('punctaMinProm',parameters),[0 1])));
    threshvalues = cell2mat(parameters(circshift(strcmp('threshvalues',parameters),[0 1])));
    maskMaxSize = cell2mat(parameters(circshift(strcmp('maskMaxSize',parameters),[0 1])));
    punctaMinArea = cell2mat(parameters(circshift(strcmp('punctaMinArea',parameters),[0 1])));
    punctaMaxArea = cell2mat(parameters(circshift(strcmp('punctaMaxArea',parameters),[0 1])));
catch
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looping through multiple thresholds                     %
%                                                         %
% Inputs: thrIm (image to be thresholded),                %
%         threshvalues (threshold values)                 %
%                                                         % 
% Outputs: contour image of puncta intensity (promIm)     %
% (this output captures the prominence of a puncta as     %
% the intensity level of the image)                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Demo mode command %
if isempty(demo) == 0
figure; subplot(1,2,1); imshow(thrIm);
end
% % % % % % % % % % %

for i=1:length(threshvalues)
    
    threshvalues(i);
    
    bwImage = im2bw(thrIm,threshvalues(i));
    
    [B,L] = bwboundaries(bwImage,'noholes');
    
    roistats = regionprops(L,'Area');
    roistats = struct2cell(roistats);
    dummyAreas = cell2mat(roistats);

    
    B(find(dummyAreas>punctaMaxArea | dummyAreas<punctaMinArea)) = [];
    
    if isempty(B) == 0
        
        for j = 1:size(B,1)

            if min(B{j}(:,2)) > 1 & min(B{j}(:,1)) > 1 & max(B{j}(:,1)) < size(thrIm,1)-1 & max(B{j}(:,2)) < size(thrIm,2)-1 % Stay away from edge pixels
            
                promIm = promIm + poly2mask(B{j}(:,2),B{j}(:,1),size(promIm,1),size(promIm,2));
                
            end
        
        end
        
                       
        % Demo mode commands %
        if isempty(demo) == 0


            subplot(1,2,2); imagesc(promIm); title('Accumulated Prominence'); signal = text(4,4,num2str(threshvalues(i)),'Color','w','FontWeight','bold','FontSize',14); %pause(0.5);
            if exist('signal') == 1; delete(signal); end

        end
        % % % % % % % % % % %
        
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering puncta on prominence and area             %
%                                                     %
% Inputs: promIm (contour image of puncta)            %
%                                                     % 
% Outputs: keepMasks (puncta masks meeting criteria)  %
%                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

promImP = promIm;

% Get Canny edges for the prominence image to avoid connected puncta

%promIm = edge(promIm,'canny');
%promIm = imfill(imdilate(promIm,strel('square',1)),'holes');

promIm = promIm.*bwareaopen(promIm,punctaMinArea);

[B,L] = bwboundaries(promIm);

punctaStats = regionprops(L,promIm,'MaxIntensity','Area','PixelValues','PixelList');
punctaStats = struct2cell(punctaStats);

dummyAreas = cell2mat(punctaStats(1,:));
dummypunctaMaxInts = cell2mat(punctaStats(4,:));


% figure;
% 
%     for j = 1:size(punctaStats,2)
%             
%             pixels_image = zeros(size(promIm));
%             pixels_xy = cell2mat(punctaStats(2,j))
%             pixels_values = cell2mat(punctaStats(3,j));
%             
%             %offset_x = floor(20 - (max(pixels_xy(:,2))-min(pixels_xy))/2);
%             %offset_y = floor(20 - (max(pixels_xy(:,1))-min(pixels_xy))/2);
% 
% 
%             for q = 1:numel(pixels_values)
%                 
%             pixels_image(pixels_xy(q,2),pixels_xy(q,1)) = pixels_values(q);
%             
%             end
% 
%             imagesc(pixels_image); pause();
%     end


B(union(find(dummypunctaMaxInts<punctaMinProm),find(dummyAreas<punctaMinArea | dummyAreas>punctaMaxArea)))=[];
keepMasks = zeros(size(thrIm,1),size(thrIm,2));

for j = 1:size(B,1)

    keepMasks = keepMasks + poly2mask(B{j}(:,2),B{j}(:,1),size(keepMasks,1),size(keepMasks,2));

end

keepMasks = bwareaopen(keepMasks, punctaMinArea);

        % Demo mode commands %
        if isempty(demo) == 0


            figure; subplot(1,3,1); imagesc(thrIm); title('Image to threshold') 
                    subplot(1,3,2); imagesc(promImP); title('Pre-Canny Prominence')
                    subplot(1,3,3); imagesc(keepMasks); title('Post-Canny and Size Filtered')
                    
        end
        % % % % % % % % % % %
        

end
