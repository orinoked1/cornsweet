% function XYZ = yxy2xyz (YXY)
% Author: Tania Pouli @ University of Bristol - 2009
function XYZ = yxy2xyz (YXY)

% If single point
if ndims(YXY) == 2 && any(size(YXY) == 3)
    flip = false;
    if (size(YXY,1) == 3)
       flip = true;
       YXY = YXY';
    end
    
    XYZ(2)          = YXY(1);
    XYZ(1)          = (YXY(1)/YXY(3))*YXY(2);
    XYZ(3)          = (YXY(1)/YXY(3))*(1- YXY(2) - YXY(3));
    
    if flip
        XYZ         = XYZ';
    end
else
    
    if ndims(YXY) < 3
        disp('ERROR: Run with 3-channel image - yxy2xyz');
        return;
    end
    
    sz				= size(YXY);

    YXY				= [reshape(YXY(:,:,1)',1,sz(1)*sz(2));...
                       reshape(YXY(:,:,2)',1,sz(1)*sz(2));...
                       reshape(YXY(:,:,3)',1,sz(1)*sz(2))];
    YXY				= double(YXY);

    W				= YXY(3,:);
    t1				= find(W > 0);
    XYZtmp			= zeros(size(YXY));

    XYZtmp(2,:)		= YXY(1,:);														% Y = Y
    XYZtmp(1,t1)	= (YXY(1,t1) ./ YXY(3,t1)) .* YXY(2,t1);						% X = x*Y/y
    XYZtmp(3,t1)	= (YXY(1,t1) ./ YXY(3,t1)) .* (1.0 - YXY(2,t1) - YXY(3,t1));	% Z = (Y/y)(1 - x - y)

    XYZ(:,:,1)		= reshape (XYZtmp(1,:), sz(2), sz(1))';
    XYZ(:,:,2)		= reshape (XYZtmp(2,:), sz(2), sz(1))';
    XYZ(:,:,3)		= reshape (XYZtmp(3,:), sz(2), sz(1))';

    if size(XYZ) ~= sz
        disp('WARNING: Oops something went wrong with the reshaping!');
        return;
    end
end