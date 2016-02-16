function simplempeg( images, varargin )
% Usage:
% 
% simplempeg( images, filename )
%
% Will create a movie file with the 
path = uigetdir();

framerate = input('\nEnter frame rate: ');

writerObj = VideoWriter([path,'\',varargin{1},'.avi']);
set(writerObj,'FrameRate',framerate)

open(writerObj);


for i = 1:size(images,3)
   imshow(images(:,:,i));
   frame = getframe;
   writeVideo(writerObj,frame);
end

close(writerObj);


end

