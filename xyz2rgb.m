% function RGB = xyz2rgb (XYZ)
%
% Author: Tania Pouli @ University of Bristol - 2009
function RGB = xyz2rgb (XYZ)
 
load('XYZ2RGB.mat');

if ndims(XYZ) == 2
    
    if (size(XYZ,1) == 3)
        RGB = XYZ2RGB*XYZ;
    elseif (size(XYZ,2) == 3)
        RGB = (XYZ2RGB*XYZ')';
    else
        disp('ERROR: Run with 3-channel image');
        return;
    end
else
    
    if ndims(XYZ) < 3
        disp('ERROR: Run with 3-channel image');
        return;
    end
    sz			= size (XYZ);

    XYZ			= [reshape(XYZ(:,:,1)',1,sz(1)*sz(2));...
                   reshape(XYZ(:,:,2)',1,sz(1)*sz(2));...
                   reshape(XYZ(:,:,3)',1,sz(1)*sz(2))];
    XYZ			= double(XYZ);
    RGBtmp		= XYZ2RGB * XYZ;

    RGB(:,:,1)	= reshape (RGBtmp(1,:), sz(2), sz(1))';
    RGB(:,:,2)	= reshape (RGBtmp(2,:), sz(2), sz(1))';
    RGB(:,:,3)	= reshape (RGBtmp(3,:), sz(2), sz(1))';

    if size(RGB) ~= sz
        disp('WARNING: Oops something went wrong with the reshaping!');
        return;
    end
end

RGB(RGB<0)  = 0;