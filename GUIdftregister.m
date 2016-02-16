function [ regimages, varargout ] = GUIdftregister( raw_images, varargin )

for i = 1:numel(varargin);
    
    if strcmp(varargin{i},'zdepth')==1
        zdepth = varargin{i+1};
    end
    
    if strcmp(varargin{i},'index')==1
        index = varargin{i+1};
    end
    
end


% Scenarios for use:
%
% (1) A MAX projected set of images or images with only one z-slice
%  
% Input: raw_images (your raw image stack having dimensions of rows, columns, and time)
%
% Output: regimages (a registered image stack having dimensions of rows, columns, and time)
%
% Syntax: [regimages] = dftregister(raw_images,'zdepth',zdepth)
%
%
% (2) A stack of images with z-slices and you want to create a MAX project and register it
%  
% Input: raw_images (your raw image stack having dimensions of rows, columns, and time)
%
% Option: interval (how many time points should be averaged -- choose number less than first z-shift in your stack)
%        
% Output: regimages (a MAX projected and registered image stack having dimensions of rows, columns, and time)
%         thrimmax (an image stack having dimensions of rows, columns, and z-slice)
%
% The variable thrimmax can be fed into an analysis if you wish to create
% masks specifically for every z-slice in the stack, then analyze the MAX projection,
% slice by slice 
% 
% Syntax: [regimages, thrimmax] = dftregister(raw_images,'zdepth',zdepth)
%
%
% (3) A stack of images with z-slices and you want to register every z-slice across frames
%
% Input: raw_images (your raw image stack having dimensions of rows, columns, and time)
%        index (the output of measureshiftFFT; an index of which z-slice goes with which after z-shift occurs)
%
% Output: regimages (a MAX projected and registered image CELL stack having dimensions of rows, columns, and time)
%         * The cells of regimages{i} contain all of the frames of z-slice i *
%
% Syntax: [regimages, thrimmax] = dftregister(raw_images, index,'zdepth',zdepth)  

num_images = size(raw_images,3); % This calculates the number of images loaded

if zdepth>1 % This is a z-stack
    
    frames = num_images/zdepth; % The number of frames is the stack divided by the number of slices
    
    if exist('index')~=1; % If you haven't specified an index
        
        fprintf('No index specified, so taking a MAX projection of images...\n');

        for i = 1:frames
        imprime(:,:,i) = (max(raw_images(:,:,[1+(i-1)*zdepth:i*zdepth]),[],3)); % Take the max over the limited third dimension (Z) of the raw images
        end
    
        images = imprime;
        
    else
        
        images = raw_images; % If you have specified an index, the rest of the routine works on the raw images you input
        
    end

else
   
    images = raw_images;  % If you have a stack with zdepth of 1, then the rest of the routine works on the raw images you input
    
end

if exist('index')==1;

    for i = 1:size(index,2);
        regimages(:,:,:,i) = raw_images(:,:,index(:,i)); % Dimensions are [X,Y,t,z]
    end
    
    
    %gfilt = fspecial('gaussian',[3 3],3);
    
    for j = 1:size(regimages,4)
        
        imagesfft = fft2(regimages(:,:,:,j));
        comparefft = fft2(regimages(:,:,1,j));
   
        for i = 1:size(regimages,3)


            %[routput] = dftregistration(imagesfft(:,:,1),imagesfft(:,:,i),100); % Registers the first frame of the ith z-slice with its jth frame
            [routput, Greg(:,:,i,j)] = dftregistration(comparefft,imagesfft(:,:,i),200); % Registers the first frame of the ith z-slice with its jth frame
            yshiftmat(i,j) = routput(3);
            xshiftmat(i,j) = routput(4);
            regimages(:,:,i,j) = real(ifft2(Greg(:,:,i,j)));
            %

            %

        end
    
    end
    
    
    
    maxxshift = ceil(max(xshiftmat,[],2));
    maxyshift = ceil(max(yshiftmat,[],2));
    
    
    [dumb,dHindex] = max(abs(maxyshift)); % This is the maximum amount that any image Y-shifted in the stack, and will be where the images are cut off
    [dumb,dWindex] = max(abs(maxxshift)); % This is the maximum amount that any image X-shifted in the stack, and will be where the images are cut off

    deltaH = maxyshift(dHindex);
    deltaW = maxxshift(dWindex);
    
    regimages = mean(regimages,4);
    regimages = regresize(regimages,findabsmax(deltaH),findabsmax(deltaW));  
    
    
else % Registering the MAX projected images 
    
    
    % Filtering the images
    unfiltim = images;
    
    gfilt = fspecial('gaussian',[20 20],3);
   
    for i = 1:size(images,3)
        images(:,:,i) = imfilter(images(:,:,i),gfilt);
    end
        
    % End filtering
    
    imagesfft = fft2(images);%(imfilter(images,h,'replicate'));

    frames = size(images,3);
    regimages = unfiltim(:,:,1);

    for i = 2:frames
        
        tic;
        
        routput(i,:) = fix(dftregistration(imagesfft(:,:,1),imagesfft(:,:,i),100));
        regimages(:,:,i) = circshift(unfiltim(:,:,i),routput(i,[3:4]));
        
        endelapsed(i) = toc;
        duration = endelapsed(i)*frames;

        %if i==1
        %h=waitbar(i/frames,sprintf('%12.0f seconds remaining',duration-sum(endelapsed)));
        %else
        %waitbar(i/frames,h,sprintf('%12.0f seconds remaining',duration-sum(endelapsed)));
        %end
    
    end
    
    %delete(h);

    [deltaH,dHindex] = max(abs(routput(:,3)));
    [deltaW,dWindex] = max(abs(routput(:,4)));

    deltaH = routput(dHindex,3);
    deltaW = routput(dWindex,4);

    regimages = regresize(regimages,findabsmax(deltaH),findabsmax(deltaW));
    

end

fprintf('Image registration complete');


% %%%%%Individual Z-Slice Thresholding Calculation%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                                                       %
% % This is for finding the images to be thresholded later if you want to %
% % get individual thresholds for each of the z-slices individually to be %
% % later used as masks in the analysis.                                  %
% %                                                                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if nargout>1 & zdepth>1;
%     
%     stacksize = zdepth;
%     
% 
%     frameind = [1:1:frames];
% 
%     shiftindex = repmat(frameind,stacksize,1);
%     shiftindex = [(1:1:num_images)',shiftindex(:)];
% 
%     for i = 1:num_images
% 
%         if i<zdepth
% 
%             thrimages(:,:,i) = regresize(raw_images(:,:,i), deltaH, deltaW);
% 
%         else
% 
%         tmp = circshift(raw_images(:,:,i),routput(shiftindex(i,2),[3:4]));
%         thrimages(:,:,i) = regresize(tmp, deltaH, deltaW);
% 
%         end
% 
%     end
%     
%     interval = zdepth*2;
%     
%     for i = 1:stacksize
% 
%         thrimmax(:,:,i) = mean(thrimages(:,:,[i:stacksize:stacksize*interval]),3);
% 
%     end
% 
% else
%     
% end


varargout{1} = xshiftmat;
varargout{2} = yshiftmat;


end

function value = findabsmax(vector)

vector1 = abs(vector);
[value,loc] = find(vector1==max(vector1),1);
value = vector1(loc)*sign(vector(loc));

end


