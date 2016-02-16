function [int,regimages] = getgcampints(time,punctatime)

    % Modifications to Imgloader_build
    %
    % To get x,y,z shift matrices into the base workspace, add lines 916-921
    % to the Imgloader_build code
    %
    % Also import the newest version of GUIdftregister (which outputs
    % xshift and yshift)
    %
    % Have this line to get the x-y registered 
    % assignin('base',sprintf('regimages%s',strcat(num2str(c(4)),num2str(c(5)))),handles.regimages);
    % 
    %
    
    shift = evalin('base',strcat('shift',num2str(time)));
    regindex = evalin('base',strcat('regindex',num2str(time)));
    ch1_images = evalin('base',strcat('regimages',num2str(time)));
    pathdummy = evalin('base', strcat('filename',num2str(time)));
    
    search_dummy = strfind(pathdummy,'\');
    
    path = pathdummy([1:search_dummy(end)]);
    
    images = GUIloadimgs('path',path,'target','Ch3');
    
    regimages = GUIdftregister(images,'index',regindex,'zdepth',size(regindex,2));
        
    puncta = evalin('base',strcat('puncta',num2str(punctatime)));

    puncta_masks = puncta.mask;
    
    % Set the size of regimages to the size of puncta_masks
    
    dw = size(regimages,1) - size(ch1_images,1);
    dh = size(regimages,2) - size(ch1_images,2);
    
    if dw<0; % If ch1_images is larger, decrease size of puncta masks
      puncta_masks = puncta_masks([1:size(regimages,1)],:);
    else
      puncta_masks = cat(1,puncta_masks,zeros(abs(dw),size(puncta_masks,2)));
    end
        
    
    if dh<0; % If ch1_images is larger, decrease size of puncta masks
      puncta_masks = puncta_masks(:,[1:size(regimages,2)]);
    else
      puncta_masks = cat(2,puncta_masks,zeros(size(puncta_masks,1),abs(dh)));
    end 
    
    [B,L] = bwboundaries(puncta_masks);
    
    for i = 1:size(regimages,3);
    
        roistats{i} = regionprops(L,regimages(:,:,i),'MeanIntensity','PixelList','PixelValues');
        roistats{i} = struct2cell(roistats{i});
        int(i,:) = cell2mat(roistats{i}(3,:));
      
    end
    
    
    filestr = strcat('gcampIntensities',num2str(punctatime));
    xlswrite([path,'\',filestr],int);
        
    fprintf('File saved as %s\n',filestr);
    
    assignin('base',strcat('gcampIntensities',num2str(punctatime)),int);

    
    h=figure('Visible','off'); 
    
    for i = 1:size(regimages,3);
        imshow(regimages(:,:,i),'InitialMagnification',500); hold on;
        title(i);
        drawB(B);
        
        try
            saveas(h,[path,'\Analyzed\',sprintf('gcampimages%1.0f',i)],'tiff');
        catch err
            if (strcmp(err.identifier,'MATLAB:saveas:invalidFilename'))
                
                mkdir([path,'\Analyzed']);
                saveas(h,[path,'\Analyzed\',sprintf('gcampimages%1.0f',i)],'tiff');
                
            end
        end
     
    end
     delete(h);
     
end

function [] = drawB(boundaries,varargin)
        
        
    color = cell2mat(varargin(1+find(strcmp([varargin],'color'))));
    if isempty(color)==1; color = 'y'; end;
    
    plabels = cell2mat(varargin(1+find(strcmp([varargin],'labels'))));
    if isempty(plabels)==1; plabels = [1:numel(boundaries)]; end;

        
          for i = 1:numel(boundaries);

            boundaryi = boundaries{i}; %handles.puncta.bounds
            plot(boundaryi(:,2),boundaryi(:,1),color);
            rndRow = ceil(length(boundaryi)\(mod(i,7)+1));
            col = boundaryi(rndRow,2); row = boundaryi(rndRow,1);

            lbl_txt = text(col+1, row-1, num2str(plabels(i)));
            set(lbl_txt,'Color','r','FontSize',12);

          end
            
              
end
    