function [ index, output ] = GUImeasureshiftFFT( imstack, stacksize )

if exist('stacksize')==0
    
    stacksize = 5;

end

% A comparison is defined as a DFT cross-correlation between an image and a
% set of images one interval away
% 
% If the interval is 10, images [1,...,5] get compared to 11,12,13,14,15 and then 21,22,23,24,25

interval = 2*stacksize; % This is 10, so image 1 gets compared to 11,12,13,14,15 and then 21,22,23,24,25

frames = size(imstack,3); % This is the TOTAL number of images in the image stack

% "Points" contains the number of comparisons the algorithm will have to make
% If the number of frames is not an integer multiple of interval, there
% will be one less comparison than if it were
% 
% This has no effect on the quality of the comparison algorithm

points = (frames-mod(frames,interval))/interval; 

% Two filters are used to smooth the images so that noise does not
% compromise the comparisons

usfilt = fspecial('unsharp',0.5);
gaussfilt = fspecial('gaussian',[10 10],10);



for i = 1:frames

    
    tic;
    
    filtim(:,:,i) = imfilter(imfilter(imstack(:,:,i), gaussfilt, 'replicate'),usfilt, 'replicate');
    imstackfft(:,:,i) = fft2(filtim(:,:,i));

    endelapsed(i) = toc;
    duration = endelapsed(i)*stacksize;

    if i==1
    h=waitbar(i/frames,sprintf('Fourier transforming images...'));
    else
    waitbar(i/frames,h,sprintf('Fourier transforming images...'));
    end
    
    clear endelapsed
    clear duration
    
end

delete(h);

for k = 1:stacksize

    tic;
    
    for i = 1:interval:interval*(points) % Iterate over each set of z slices
    
        for j = 1:stacksize
        
        output = dftregistration(imstackfft(:,:,k),imstackfft(:,:,i+j-1),20);
        
        maxcorrmat(i+j-1,k) = 1-output(1);
        
        end
        
    end
    
    
    endelapsed(k) = toc;
    duration = endelapsed(k)*stacksize;

    if k==1
    h=waitbar(k/stacksize,sprintf('%12.0f seconds remaining',duration-sum(endelapsed)));
    else
    waitbar(k/stacksize,h,sprintf('%12.0f seconds remaining',duration-sum(endelapsed)));
    end
    
end

delete(h);

maxcorrmat = maxcorrmat(maxcorrmat~=0);
maxcorrmat = reshape(maxcorrmat,stacksize,length(maxcorrmat)/stacksize);

%figure; imagesc(maxcorrmat);

[maxcorrmat,index] = max(maxcorrmat);

%figure; plot(reshape(maxcorrmat,points,stacksize));

xlbl = [0:(interval/stacksize):points*interval];


figure('Color','w');
    
for i = 1:stacksize
    
       subplot(stacksize,1,i); scatter([1:points-1],index([points*(i-1)+2:points*(i)]),10.^maxcorrmat([points*(i-1)+2:points*(i)]),'filled'); ylim([1,5]);
       
       set(gca,'XTick',[0:1:points],'XTickLabel',xlbl); title(['\fontsize{24}Z =',num2str(i)]); 

end


%Expand the index matrix since the interval is twice the stacksize
indexadj = zeros(1,points*interval);
indexadj([1:2:length(indexadj)]) = index;
indexadj([2:2:length(indexadj)]) = index;

indexadj = reshape(indexadj,length(indexadj)/stacksize,stacksize);
indexadj1 = repmat([0:stacksize:(points*interval-stacksize)]',1,5);

index = indexadj+indexadj1;
%Then turn it into the correct address in the image set

end


