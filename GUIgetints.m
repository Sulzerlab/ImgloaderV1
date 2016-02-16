function [ puncta, groups, bg ] = GUIgetints( images, varargin )

%%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %                                
%                                                          %
% puncta                                                   %
%  .mask: a thresholded image with the masks shown         %
%  .bounds: all of the boundaries of the masks             %
%  .L: a label matrix of the masks                         %
%  .total: the total number of masks                       %
%  .centroids: centroids for each of the masks             %
%  .areas: areas of each of the masks                      %
%  .int: raw intensity values for each of the masks        %
%  .pixels: pixel values of each point in each of the masks%
%  .bgrem_int: background-subtracted intensity values      %
%  .norm_int: normalized intensity values                  %
%  .norm_bgrem_int: normalized and background subtracted   %
%                                                          %
% groups.total: a puncta x 4 matrix with the four columns: %
%   column 1: an index of all of the masks found           %
%   column 2: a logical 1/0 for I(t) > Icut, for any t     %
%   column 3: of those masks with 1 in column 2, PCA group %
%                                                          %
%             Depending on k-groups desired,               %
%             gives 1, 2, ..., n in this column specifying %
%             which cluster the puncta belongs to          %
%                                                          %
%                                                          %
%   column 4: of those masks with 1 in column 2, IQ group  %
%                                                          %
%           1= drop < 1.5 * IQ range                       %
%           2= drop is between 1.5 and 3 of IQ range       %             
%           3= drop is > 3 * IQ range                      %
%                                                          %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



frames = size(images,3);

%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %                                
%                                                          %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargin = varargin{1};

stimperiod = 19;

scale = cell2mat(varargin(1+find(strcmp([varargin],'scale'))));
    if isempty(scale)==1; scale = 1; end;

beginframe = cell2mat(varargin(1+find(strcmp([varargin],'beginframe'))));
    if isempty(beginframe)==1; beginframe = 1; end;

endframe = cell2mat(varargin(1+find(strcmp([varargin],'endframe'))));
    if isempty(endframe)==1; endframe = frames-1; end;

startstim = cell2mat(varargin(1+find(strcmp([varargin],'startstim'))));
    if isempty(startstim)==1; startstim = 20; end;

blineremove = cell2mat(varargin(1+find(strcmp([varargin],'blineremove'))));
    if isempty(blineremove)==1; blineremove = 0; end;
    
bgremove = cell2mat(varargin(1+find(strcmp([varargin],'bgremove'))));
    if isempty(bgremove)==1; bgremove = 0; end;
    
pcabegin = cell2mat(varargin(1+find(strcmp([varargin],'pcabegin'))));
    if isempty(pcabegin)==1; pcabegin = 1; end;
    
pcaend = cell2mat(varargin(1+find(strcmp([varargin],'pcaend'))));
    if isempty(pcaend)==1; pcaend = 1; end;

maskmin = cell2mat(varargin(1+find(strcmp([varargin],'maskmin'))));
    if isempty(maskmin)==1; maskmin = 20; end;
      
maskmax = cell2mat(varargin(1+find(strcmp([varargin],'maskmax'))));
    if isempty(maskmax)==1; maskmax = 50; end;
        
icut = cell2mat(varargin(1+find(strcmp([varargin],'icut'))));
    if isempty(icut)==1; icut = 2; end;

corrcheck = cell2mat(varargin(1+find(strcmp([varargin],'corrcheck'))));
    if isempty(corrcheck)==1; corrcheck = 0; end;
    
pcas = cell2mat(varargin(1+find(strcmp([varargin],'pcas'))));
    if isempty(pcas)==1; pcas = 3; end;
    
kgroups = cell2mat(varargin(1+find(strcmp([varargin],'kgroups'))));
    if isempty(kgroups)==1; kgroups = 2; end;
    
dropcriteria = cell2mat(varargin(1+find(strcmp([varargin],'dropcriteria'))));
    if isempty(dropcriteria)==1; dropcriteria = 35; end;
    
method = cell2mat(varargin(1+find(strcmp([varargin],'method'))));
    if isempty(method)==1; method = 1; end;

dropstart = cell2mat(varargin(1+find(strcmp([varargin],'dropstart'))));
    if isempty(dropstart)==1; dropstart = 21; end;

dropend = cell2mat(varargin(1+find(strcmp([varargin],'dropend'))));
    if isempty(dropend)==1; dropend = 60; end;

xlsfile = (varargin(1+find(strcmp([varargin],'xlsave'))));
    if isempty(xlsfile)==1; xlsave=0; 
    else
        xlsave=1;
    end;
    
maskimage = varargin(1+find(strcmp([varargin],'masks')));

masklife = cell2mat(varargin(1+find(strcmp([varargin],'masklife'))));
    if isempty(masklife)==1; masklife = 1; end;
    
datestr = cell2mat(varargin(1+find(strcmp([varargin],'datestr'))));
if isempty(datestr)==1; datestr = 0; end;
    
baseline = [1:startstim];

% If you have not specified an image to threshold then one will be created
% from images within the interval 

if isempty(maskimage)==0
    
    fprintf('Using preselected masks\n');
    masks = thresholding(images,maskimage);
    
else
    
    if masklife == 1

        puncta.mask = thresholding(mean(images(:,:,[beginframe:startstim]),3),maskmin,maskmax,'ws');
    
    else
    
        masklifeframes(1:2:startstim) = [1:startstim/2]; 
        masklifeframes(2:2:startstim) = [startstim:-1:startstim/2+1];

        
        for i = 1:masklife


            tic;

            puncta.mask(:,:,i) = thresholding(images(:,:,masklifeframes(i)),maskmin,maskmax,'ws');

            endelapsed(i) = toc;
            duration = mean(endelapsed)*masklife;

            if i==1
            h=waitbar(i/masklife,sprintf('Approximate time remaining %2.0f seconds (frame %2.0f)',duration-sum(endelapsed),masklifeframes(i)));
            else
            waitbar(i/masklife,h,sprintf('Approximate time remaining %2.0f seconds (frame %2.0f)',duration-sum(endelapsed),masklifeframes(i)));
            end

            clear duration


        end
        delete(h);

        mastermask = thresholding(mean(images(:,:,[beginframe:startstim]),3),maskmin,maskmax,'ws'); 
        discrimmask = mat2gray(sum(puncta.mask,3),[0 masklife]);
        
        %discrimmask = im2bw(discrimmask,1);
        discrimmask = discrimmask.*mastermask;
        
        [B,L] = bwboundaries(discrimmask,'noholes');

        dummyRegprops = cell2mat(struct2cell(regionprops(L, discrimmask,'MaxIntensity')));

        [dummyKeep] = find(dummyRegprops==1); % If masklife is an absolute requirement (otherwise < 1)

        dummyKeepMat = zeros(size(L)); 

        for i = 1:numel(dummyKeep)

            dummyKeepMat(L==dummyKeep(i))=1;

        end


        % Matching up masks by centroid

        [B1,L1] = bwboundaries(mastermask);

        R1 = regionprops(L1,mastermask,'Centroid');

        R1pixelIdx = struct2cell(regionprops(L1,'PixelIdxList'));

        R1 = circshift(reshape(fix(cell2mat(struct2cell(R1))),2,numel(B1)),[1 0])';
        R1 = sub2ind(size(L1),R1(:,1),R1(:,2));

        toLose = dummyKeepMat(R1);

        R1pixelIdx(toLose==0)=[];
        R1pixelIdx = vertcat(R1pixelIdx{1,:});

        puncta.mask = zeros(size(L1));
        puncta.mask(R1pixelIdx) = 1;

        %puncta.mask = dummyKeepMat.*mastermask;

        % Demo mode
        %figure('Color','w'); imshowpair(puncta.mask,mastermask); title('Purple = Lost, White = Kept');

end
    
end
    

% The background is everything other than the selected masks

bg.mask = ~puncta.mask; 

[puncta.bounds,puncta.L] = bwboundaries(puncta.mask,'noholes');

fprintf('\nThe number of masks is %1.0f\n',numel(puncta.bounds));

puncta.total = numel(puncta.bounds);
roistatcentroids = regionprops(puncta.L,puncta.mask,'Centroid');
puncta.centroids = reshape(cell2mat(struct2cell(roistatcentroids)),puncta.total,2);
roistatareas = regionprops(puncta.L,puncta.mask,'Area');
puncta.areas = reshape(cell2mat(struct2cell(roistatareas)),puncta.total,1);


puncta.int = zeros(frames,puncta.total);

    for i = 1:frames;

        roistats{i} = regionprops(puncta.L,images(:,:,i),'MeanIntensity','PixelList','PixelValues');
        roistats{i} = struct2cell(roistats{i});
        puncta.int(i,:) = cell2mat(roistats{i}(3,:));
        
    end;

    

% %%%%%%%%%%%%%%%%%
% % Pixels        %
% %%%%%%%%%%%%%%%%%
% 
% pixels = zeros(20,20,size(puncta.int,2),frames);
% 
% 
% for i = 1:frames;
% 
% for j = 1:size(puncta.int,2)
%         
% 
%         %pixels(:,:,j,i) = zeros(20,20);
%         
%         pixels_x = roistats{i}{1,j}(:,1) - min(roistats{i}{1,j}(:,1))+1;
%         pixels_y = roistats{i}{1,j}(:,2) - min(roistats{i}{1,j}(:,2))+1;
%         
%         offset_x = floor(20 - (max(pixels_x)-min(pixels_x))/2);
%         offset_y = floor(20 - (max(pixels_y)-min(pixels_y))/2);
%         
% 
%         for q = 1:length(pixels_x)
%         puncta.pixels(pixels_x(q)+offset_x,pixels_y(q)+offset_y,j,i) = roistats{i}{2,j}(q);
%         end
%        
% end
% 
% end
% 
% 
% %%%%%%%%%%%%%%%%%
% % End Pixels    %
% %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%
% Background    %
%%%%%%%%%%%%%%%%%


[bg.bounds,bg.L,Nbgmasks] = bwboundaries(bg.mask,'noholes');

bg.L(bg.L>1)=0; % Set L>1 criteria to get one uniform background
       
    for i = 1:frames;

        roistatsDIL{i} = regionprops(bg.L,images(:,:,i),'PixelValues');
        bg.int(i,:) = mean(roistatsDIL{i}.PixelValues);

        
    end;
 
    
roistatcentroids = regionprops(bg.L,bg.mask,'Centroid');
bg.centroids = [roistatcentroids.Centroid(1),roistatcentroids.Centroid(2)];



%%%%%%%%%%%%%%%%%%%%%%
% Normalization here %
%%%%%%%%%%%%%%%%%%%%%%


puncta.bgrem_int = puncta.int-repmat(bg.int,1,puncta.total); % Get background subtracted intensities and place in bgremINT

puncta.bgrem_int = puncta.bgrem_int+1.01*abs(min(min(puncta.bgrem_int))); % Offset to prevent any negative values

%%%%%%%%%%%%%%%%%%%%%%%
% Add back minimum    %
%%%%%%%%%%%%%%%%%%%%%%%

diagfunc = @(matrix) diag(1./matrix(startstim,:)); % Baseline is the mean of the puncta intensity before stimulation

% Calculate the norms for each series and place on the diagonal %
% The normalized background
bg.norm_int = bg.int * diagfunc(bg.int); 

% The normalized intensities

if bgremove == 1

    puncta.norm_int = puncta.int * diagfunc(puncta.int);
    puncta.norm_bgrem_int = puncta.bgrem_int * diagfunc(puncta.bgrem_int);

else

    puncta.norm_int = puncta.int * diagfunc(puncta.int);
    puncta.norm_bgrem_int = puncta.norm_int;

end

%%%%%%%%%%%%%%%%%%%%%%
% Baseline adjustment%
%%%%%%%%%%%%%%%%%%%%%%

if blineremove == 1

puncta.fit = puncta.norm_bgrem_int;

    for i = 1:puncta.total
        %startstim = frames; % Added this for the unstimulated data
        
        X = [0:startstim-1];
        Y = puncta.norm_bgrem_int([1:startstim],i)'; % Set the current 

        [S,I,R,stats] = logfit(X,Y,'logy');

        S(S>0)=0; % Set any curves with positive slopes to zero -- these will not be fit later

        yApprox = (10^I)*(10^S).^[0:1:frames-1]';

        puncta.yApprox(:,i) = yApprox;
        puncta.Sinit(:,i) = S;

        Check = medfilt1(diff(medfilt1(puncta.norm_bgrem_int([startstim:frames],i) - yApprox([startstim:frames]),10)),5);

        clear endFit

        [dummy,endFit] = max(Check);

        if isempty(endFit) == 0

            endFit = endFit+startstim-1;
            puncta.yApprox([endFit:frames],i)=puncta.yApprox(endFit,i)*ones(frames-endFit+1,1);

        end

        puncta.fit(:,i) = puncta.norm_bgrem_int(:,i)./puncta.yApprox(:,i);
    
    end

puncta.norm_bgrem_intuncorr = puncta.norm_bgrem_int; % The uncorrected intensities
puncta.norm_bgrem_int = puncta.fit; % The corrected intensities

end

%startstim = 20; % Added this for unstimulated
%%%%%%%%%%%%%%%%%%%%%%
% Exp fitting        %
%%%%%%%%%%%%%%%%%%%%%%
    
for i = 1:puncta.total
    
    X = [0:0+stimperiod];
    Y = puncta.norm_bgrem_int([startstim:startstim+stimperiod],i)';

    [puncta.S(i),I,R,stats] = logfit(X,Y,'logy');

end

puncta.S = -(1./(2.303*puncta.S))';

%%%%%%%%%%%%%%%%%%%%%%
% Grouping Hierarchy %
%%%%%%%%%%%%%%%%%%%%%%

groups.total(:,1) = [1:1:puncta.total]'; % Every single puncta indexed into the first column

[dummy,survivors] = find(puncta.norm_bgrem_int([1:end],:)>icut);
fprintf('The number of puncta excluded: %d\n',numel(unique(survivors)));

groups.total(~ismember(groups.total,unique(survivors)),2) = 1; % Include any puncta that don't surpass the Icut by setting this row to 1

clear survivors


analysisgroup = groups.total(groups.total(:,2)==1,1); % Analyze those members that are "1" in the second column

%%%%%%%%%%%%%%%%%%%%%%
% Correlation Groups %
%%%%%%%%%%%%%%%%%%%%%%

% If correlation check is requested, then create a large matrix with the background and the puncta intensities
% and calculate their correlation coefficient  
% The negatively correlated pairs will receive a "2" in the second column
% of groups.total so now groups.total's second column has the key
% 
% 0 = Excluded on Icut
% 1 = All included puncta
% 2 = Puncta anti-correlated with background
%

if corrcheck == 1

    for i = 1:numel(analysisgroup)
        smooth_puncta(:,i) = smooth(puncta.int(:,analysisgroup(i)),11); 
    end

    smooth_bg = repmat(smooth(bg.int,11),1,numel(analysisgroup));

    corrmat = horzcat(smooth_puncta([startstim:startstim+10],:),smooth_bg([startstim:startstim+10],:)); 
    coeffmat = corrcoef(corrmat);
    coeffdiag = diag(coeffmat(1:numel(analysisgroup),numel(analysisgroup)+1:end));

        fprintf('\nAnalyzing only those puncta which anti-correlate with background after frame %1.0f',startstim);
        groups.total(analysisgroup(find(coeffdiag>=0)),2) = 2;
        [analysisgroupindex,blank] = find(groups.total(:,2)>0);
        analysisgroup = groups.total(groups.total(:,2)==1,1);

end

[coeff,score,latent] = princomp(puncta.norm_bgrem_int([pcabegin:pcaend],analysisgroup),'econ');

%%%%%%%%%%%%%%%%%%%%%%
%      K-groups      %
%%%%%%%%%%%%%%%%%%%%%%
   
[kcomb,icomb,sumd,d] = kmeans((coeff(:,[1:pcas])),kgroups,'replicates',100);

groups.total(analysisgroup,3) = kcomb;

clear analysisgroup



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intensity-based methods %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


analysisgroup = groups.total(groups.total(:,2)>0,1); 


switch method
        
% Method 1:

case 1
fprintf('Drop Criteria Method with Requirement of %1.0f\n',dropcriteria);

    for i = 1:numel(analysisgroup)

        dummymean = mean(puncta.norm_bgrem_int([1:startstim],analysisgroup(i)),1);
        dummystd = std(puncta.norm_bgrem_int([1:startstim],analysisgroup(i)),[],1);

        dummydist = (puncta.norm_bgrem_int(:,analysisgroup(i))-dummymean)/dummystd;

        testparam(i) = sum(dummydist([dropstart:dropend])<-2);
        %subplot(2,1,1); plot(smooth(puncta.norm_bgrem_int(:,analysisgroup(i)),5)); subplot(2,1,2); hist(dummydist); title(sum(dummydist([20:50])<-3)); pause();
        if sum(dummydist([dropstart:dropend])<-2) < dropcriteria

            groups.total(analysisgroup(i),4) = 1;

        else

            groups.total(analysisgroup(i),4) = 2;

        end
        %%subplot(2,1,1); plot(puncta.norm_bgrem_int(:,i)); subplot(2,1,2); hist(dummydist); pause();

    end

% Method 2:

case 2
afterstiminterval = [startstim+8:startstim+12];

fprintf('IQ Method\n');

    for i = 1:numel(analysisgroup)


        baseline = prctile(puncta.norm_bgrem_int([1:startstim],analysisgroup(i)),[25 75]);
        avgofts = mean(puncta.norm_bgrem_int(afterstiminterval,analysisgroup(i)),1);

        comparator = (baseline(1) - avgofts)/(baseline(2)-baseline(1)); % 25-prctile compared to 1.5*IQ range

        if comparator < 1.5; groups.total(analysisgroup(i),4) = 1; 
        elseif (1.5 < comparator && comparator < 3); groups.total(analysisgroup(i),4) = 2; 
        elseif (comparator > 3); groups.total(analysisgroup(i),4) = 3; 
        end

    end
% Method 3:

case 3

    for i = 1:numel(analysisgroup)

        dummymean = mean(puncta.norm_bgrem_int([1:startstim],i),1);
        dummystd = std(puncta.norm_bgrem_int([1:startstim],i),[],1);

        dummydist = (puncta.norm_bgrem_int(:,i)-dummymean)/dummystd;

        ksResult(i) = kstest2((dummymean-normrnd(dummymean,dummystd,[1,1000]))/dummystd-3,dummydist([20:end]),'tail','smaller');

        %%subplot(2,1,1); plot(puncta.norm_bgrem_int(:,i)); subplot(2,1,2); hist(dummydist); pause();
        groups.total(analysisgroup(i),4) = 2*ksResult(i);

    end
    
end
    

%if exist('testparam')==1
%figure; cdfplot(testparam)%,coeff(:,1),'.');
%end

%%%%%%%%%%%%%%%%%%%%%%%
% Matching up Groups  %
%%%%%%%%%%%%%%%%%%%%%%%

% Create an ordering for the PCs so that the colors match-up between PC
% groups and IQ groups
%
% This is done by ordering the normalized average after-stimulus fluorescence of
% the puncta. Those traces that have smaller minimum values get placed
% lower

for i =1:kgroups
    
    PCorder(i,:) = -[min(mean(puncta.norm_bgrem_int([startstim:frames],groups.total(:,3)==i),2))]; 
% Negative sign creates a sorting along descent 
% (so that greater means are associated with smaller group index)
end

    [PCorder,PCordernew] = sortrows(PCorder);
    newPCgroups = zeros(max(groups.total(:,1)),1);
    
    for i = 1:kgroups
        newPCgroups(groups.total(:,3)==PCordernew(i))=i;
    end
    
    newPCgroups = horzcat(groups.total(:,3),newPCgroups);
    
    groups.total(:,3) = newPCgroups(:,2);

for i =1:kgroups    
    PCorder(i,:) = [mean(mean(puncta.norm_bgrem_int([frames-10:frames],groups.total(:,3)==i),2))]; 
end


%%%%%%%%%%%%%%%%%%%%%%%
% Group Visualization %
%%%%%%%%%%%%%%%%%%%%%%%

 explained(1) = 0;

for i = 1:length(latent)-1; explained(i) = sum(latent(1:i))/sum(latent); end
fprintf('Five PCs can explain %2.0f percent of variance\n',explained(5)*100);

vis=1; % Set this to zero if you don't want all of the plots to pop-up


kcolors = lines(kgroups+1);
pcacolors = jet(pcas);

if vis==1;

    

    % Of PCA %

    figure('Color','w','position',[0 0 1200 800]);
    pcacolors = hsv(pcas);
    for i = 1:pcas
    subplot(2,3,1); hold on; plot(score(:,i),'Color',pcacolors(i,:),'LineWidth',4);
    set(gca,'LineWidth',4);
    legendstring{i} = sprintf('PCA %1.0f',i);
    end
    

    legend(legendstring);
    clear legendstring;

    % Of K-Groups %

    for i = 1:kgroups

        subplot(2,3,2);plot3(coeff(kcomb==i,1),coeff(kcomb==i,2),coeff(kcomb==i,3),'o','MarkerFaceColor',kcolors(i,:),'MarkerSize',8); hold on;
        legendstring{i} = sprintf('K = %1.0f',i);

    end

    legend(legendstring);
    xlabel('PCA 1'); ylabel('PCA 2'); zlabel('PCA 3');

    subplot(2,3,3); plot(explained(1:end),'bo','MarkerSize',6);
    subplot(2,3,4); boxplot(d);   

    figure('Color','w');

    for i = 1:kgroups

        if kgroups == 2;

            %tauRange = [0:(1/10)*std(puncta.S(groups.total(:,3)==i)):max(puncta.S(groups.total(:,3)==i))];    
            tauRange = [0:10:100];
            subplot(kgroups,2,i); cdfplot(puncta.S(groups.total(:,3)==i)); title(sprintf('PCA group %1.0f',i)); xlim([0 100]);
            subplot(kgroups,2,i+2); cdfplot(puncta.S(groups.total(:,4)==i)); title(sprintf('IQ group %1.0f',i)); xlim([0 100]);

        else

        %tauRange = [0:(1/10)*std(puncta.S(groups.total(:,3)==i)):max(puncta.S(groups.total(:,3)==i))];    
        tauRange = [0:10:100];
        subplot(kgroups,1,i); cdfplot(puncta.S(groups.total(:,3)==i));

        end
    
    end

Figs = fix(max(get(0,'Children')));  
figure(Figs+1); set(Figs+1,'Visible','off','Color','w');
figure(Figs+2); set(Figs+2,'Visible','off','Color','w');
figure(Figs+3); set(Figs+3,'Visible','off','Color','w');
ymin = 0;
ymax = icut;

    for i = 1:3
        figure(Figs+1); subplot(1,3,i); hold on; plot(puncta.norm_bgrem_int(:,groups.total(:,3)==i)); axis([0 frames ymin ymax]); title(sprintf('PCA group %1.0f ALL',i),'FontSize',16)
        figure(Figs+2); subplot(1,3,i); hold on; plot(puncta.norm_bgrem_int(:,groups.total(:,4)==i)); axis([0 frames ymin ymax]); title(sprintf('IQ group %1.0f ALL',i),'FontSize',16)
        figure(Figs+3); subplot(2,1,1); hold on; plot(mean(puncta.norm_bgrem_int(:,groups.total(:,3)==i),2),'LineWidth',4,'Color',kcolors(i,:)); axis([0 frames ymin ymax]); title(sprintf('PCA groups',i),'FontSize',16)
        figure(Figs+3); subplot(2,1,2); hold on; plot(mean(puncta.norm_bgrem_int(:,groups.total(:,4)==i),2),'LineWidth',4,'Color',kcolors(i,:)); axis([0 frames ymin ymax]); title(sprintf('IQ groups',i),'FontSize',16)
        pcalegend{i} = sprintf('Group %1.0f (%1.0f)',i,sum(groups.total(:,3)==i)); 
        iqlegend{i} = sprintf('Group %1.0f (%1.0f)',i,sum(groups.total(:,4)==i));
    end
   
figure(Figs+3); subplot(2,1,1); legend(pcalegend{1},pcalegend{2},pcalegend{3});
figure(Figs+3); subplot(2,1,2); legend(iqlegend{1},iqlegend{2},iqlegend{3});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
% Finding members outside of intersection  %
% b/w categorical variable 1 (PCA group)   %
% and 2 (IQ group)                         %
%
% Notice that if three K-groups are used   %
% then the alignment between PC groups and %
% K-groups is 1:1 so column 3 and column 4 %
% are directly compared                    %
%                                          %
% OTOH, if n K-groups are used then the    %
% first K-group members are identified by  %
% "1" and all others are identified by "2" %
% and the "non-outliers" are identified by %
% "1" and all others are identified by "2" %
%                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch kgroups
    
    case 3
    
    % No re-grouping is necessary because PC k-group 1,2,3 corresponds
    % directly to IQ group 1,2,3
        
    iqcol = 4;
    
    fprintf('Methods agree %3.0f percent\n',100*numel(find(groups.total(:,3)==groups.total(:,iqcol)))/numel(analysisgroup));
    groups.disagree = find(groups.total(:,3)~=groups.total(:,4));
    
    case 2

    % Regrouping is necessary because PC k-group 1,2 corresponds
    % to IQ group 1,{2,3} so create a fifth column (iqcol = 5) that will
    % regroup the puncta giving 1 if IQ group is 1, and 2 if the IQ group
    % is either 2 or 3
    
    iqcol = 5;
    
    groups.total(groups.total(:,4)>1,iqcol)=2; 
    groups.total(groups.total(:,4)==1,iqcol)=1; 
    
    fprintf('Methods agree %3.0f percent',100*numel(find(groups.total(:,3)==groups.total(:,iqcol)))/numel(analysisgroup));
    groups.disagree = find(groups.total(:,3)~=groups.total(:,iqcol));
    
    otherwise
    
    iqcol = 5;
    pcacol = 6;
        
    groups.total(groups.total(:,4)>1,iqcol)=2; 
    groups.total(groups.total(:,4)==1,iqcol)=1;
    
    groups.total(groups.total(:,3)>1,pcacol)=2; 
    groups.total(groups.total(:,3)==1,pcacol)=1;
    
    fprintf('Methods agree %1.0f percent',100*numel(find(groups.total(:,iqcol)==groups.total(:,pcacol)))/puncta.total);
    groups.disagree = find(groups.total(:,pcacol)~=groups.total(:,iqcol));
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creating a verbal key for groups        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1:size(groups.total,1)
        switch groups.total(j,3)
            case 0 
                groups.totalkey{j,1} = 'Exc';
            case 1
                groups.totalkey{j,1} = 'ND';
            case 2
                groups.totalkey{j,1} = 'Slow';
            case 3
                groups.totalkey{j,1} = 'Fast';
        end

        switch groups.total(j,4)
            case 0 
                groups.totalkey{j,2} = 'Exc';
            case 1
                groups.totalkey{j,2} = 'ND';
            case 2
                groups.totalkey{j,2} = 'Slow';
            case 3
                groups.totalkey{j,2} = 'Fast';
        end
            
end

for i = 1:kgroups
    
    groups.agree{i} = find(groups.total(:,3)==i & groups.total(:,iqcol)==i);
           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Separate plots of non-intersection set  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if vis == 1

    figure(Figs+4); set(Figs+4,'Color','w');

    for i = 1:numel(groups.disagree)
        h=subplot(3,(numel(groups.disagree) - mod(numel(groups.disagree),3)+3)/3,i); plot(puncta.norm_bgrem_int(:,groups.disagree(i)),'LineWidth',2,'Color',kcolors(groups.total(groups.disagree(i),3),:)); 
        set(h,'XMinorGrid','on'); 
        axis([0 frames ymin ymax]);
        title(sprintf('Puncta %1.0f IQ %s',groups.disagree(i),groups.totalkey{groups.disagree(i),2}));
    end



    figure(Figs+5); set(Figs+4,'Color','w');

    for i = 1:kgroups

        plot(mean(puncta.norm_bgrem_int(:,groups.agree{i}),2),'Color',kcolors(i,:),'LineWidth',2); hold on;
        plot(mean(puncta.norm_bgrem_int(:,find(groups.total(:,iqcol)==i)),2),'Color',kcolors(i,:),'LineWidth',2,'LineStyle','--'); hold on;
        title('Dashed Line = Before Removal; Solid Line = After Removal');

    end

end

fprintf('\nTotal destaining puncta %2.0f percent\n',(100/numel(find(groups.total(:,2)>0)))*numel(groups.agree{2}));

%%%%%%%%%%%%%%%%%
% Output        %
%%%%%%%%%%%%%%%%%
    if xlsave == 1
           
        datestr(datestr=='-')='';
        
        warning('off','MATLAB:xlswrite:NoCOMServer');

        disagree = zeros(1,size(groups.total,1));
        disagree(groups.agree{2}) = 1;
        
        path = uigetdir();
        
        filestr = strcat(cell2mat(xlsfile{1}),'Intensities',datestr);
        xlswrite([path,'\',filestr],[groups.total,disagree',puncta.norm_bgrem_int']');
        
        fprintf('File saved as %s\n',filestr);
                
        filestr = strcat(cell2mat(xlsfile{1}),'BGintensity',datestr);
        xlswrite([path,'\',filestr],bg.int);
        
        fprintf('File saved as %s\n',filestr);
        
    end    

end
%%%%%%%%%%%%%%%%%
% End Output    %
%%%%%%%%%%%%%%%%%
   

