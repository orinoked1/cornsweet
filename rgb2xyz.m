% function XYZ = rgb2xyz (RGB)
% 
% Author: Tania Pouli @ University of Bristol - 2009
function XYZ = rgb2xyz(RGB)
load RGB2XYZ;

% If single point or vector
if ndims(RGB) == 2
    
    if (size(RGB,1) == 3)
        XYZ = RGB2XYZ*RGB;
    elseif (size(RGB,2) == 3)
        XYZ = (RGB2XYZ*RGB')';
    else
        disp('ERROR: Run with 3-channel image');
        return;
    end
else
    
    if ndims(RGB) < 3
        disp('ERROR: Run with 3-channel image');
        return;
    end
    sz			= size (RGB);

    RGB			= [reshape(RGB(:,:,1)',1,sz(1)*sz(2));...
                   reshape(RGB(:,:,2)',1,sz(1)*sz(2));...
                   reshape(RGB(:,:,3)',1,sz(1)*sz(2))];
    RGB			= double(RGB);
    XYZtmp		= RGB2XYZ * RGB;


    XYZ(:,:,1)	= reshape (XYZtmp(1,:), sz(2), sz(1))';
    XYZ(:,:,2)	= reshape (XYZtmp(2,:), sz(2), sz(1))';
    XYZ(:,:,3)	= reshape (XYZtmp(3,:), sz(2), sz(1))';

    XYZ(XYZ<0)  = 0;
    if size(XYZ) ~= sz
        disp('WARNING: Oops something went wrong with the reshaping!');
        return;
    end
end