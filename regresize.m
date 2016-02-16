function [ regimages ] = regresize(imagesin, deltaH, deltaW)

H = size(imagesin,1); % Height of the image
W = size(imagesin,2); % Width of the image

if deltaH == 0; deltaH = 1; end
if deltaW == 0; deltaW = 1; end

if deltaH < 0
    
    if deltaW < 0
        
        regimages = imagesin([1:H+deltaH],[1:W+deltaW],:,:);
        
    else
    
        regimages = imagesin([1:H+deltaH],[deltaW:W],:,:);

    end
    
    
end

if deltaH > 0

    if deltaW < 0
      
        regimages = imagesin([deltaH:H],[1:W+deltaW],:,:);
    
    else
        
        regimages = imagesin([deltaH:H],[deltaW:W],:,:);

    end
    
    
end

