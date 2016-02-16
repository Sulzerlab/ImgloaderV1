function [ loadedimages, fname ] = GUIloadimgs(varargin) % If loadstyle is set to 1 then the program will MAX project the image sequence

% Last modified 10/21 to take a target to specify filename inputs
% This is to allow for Ch3 or Ch4 (gcamp or ffn)
% Program is no longer compatible with JPEGs (line 49)

global frames;
pathname = cell2mat(varargin(1+find(strcmp([varargin],'path'))));
target = cell2mat(varargin(1+find(strcmp([varargin],'target'))));
if isempty(target)==1; target = 'tif'; end % Default is to find tif files, otherwise a target string may be specified
 

if isempty(pathname) == 1

    [filename, pathname] = uigetfile('*.*');

    fprintf('Loading files from %s',pathname);
    filelist = dir(pathname);
    
else
    
    pathname = sprintf('%s/',pathname);
    filelist = dir(pathname);
    filename = filelist(end).name;

end

for i = 1:numel(varargin)

    if strcmp(varargin{i},'bitdepth')==1
        bitdepth = str2double(varargin{i+1});
    end
    
end

fprintf('Loading files from %s',pathname);

filelist = dir(pathname);

frames = 1;
findex = 0;

for i=3:size(filelist,1)
 
    if filelist(i).isdir ~= true
	        fname = filelist(i).name;
	        % get only if file extension is jpg or tif and check similarity
	        % in file name
	        if isempty(strfind( fname ,target  )) == 0 && strcmp( fname(size(fname,2)-3:size(fname,2)) ,'.tif'  ) == 1 
                
                if exist('im') == 0
                    
                    im(:,:,frames) = imread([pathname fname]);
                    frames = frames + 1;
                    
                else
                    
                    tmp = imread([pathname fname]);
                    
                    if size(tmp,1) == size(im,1) && size(tmp,2) == size(im,2)

                        im(:,:,frames) = tmp;
                        frames = frames + 1;
                    
                    end
                    
                end
                    
                    ffname{frames}= fname;
                    %fprintf('Loaded image %s into image %1.0f.\n',fname,frames-1);

            end
          
    end
    
end


%BitDepth = info.BitDepth
if exist('bitdepth')==0
    
bitdepth = log2(double(max(max(max(im)))));
fprintf('\n[Min Max] of image stack = [%1.0f %1.0f]\n',min(min(min(im))),max(max(max(im))));
fprintf('\nEstimating bit depth of %1.0f bits by rounding up from %1.2f\n',ceil(bitdepth),bitdepth)

else
    
    estbitdepth = log2(double(max(max(max(im)))));
    
    if abs(estbitdepth-bitdepth)/max([bitdepth,estbitdepth]) > 0.25
        
    warndlg(sprintf('Entered bitdepth of %1.0f was replaced with estimate of %1.0f',bitdepth,estbitdepth));
    bitdepth = estbitdepth;
        
    end
    

end

fprintf('Loaded %1.0f images.\n',size(im,3));
fname = strcat(pathname,fname);

loadedimages = mat2gray(im,[0 2^ceil(bitdepth)]);

end