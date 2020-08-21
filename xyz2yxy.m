%% function YXY = xyz2yxy (XYZ)
% Author: Tania Pouli @ University of Bristol - 2009
function YXY = xyz2yxy (XYZ)


if ndims(XYZ) == 2 && any(size(XYZ) == 3)
    flip = false;
    if (size(XYZ,1) == 3)
       flip = true;
       XYZ = XYZ';
    end
    
    W		= XYZ(1) + XYZ(2) + XYZ(3);
    
    YXY(1)	= XYZ(2);			% Y = Y
    YXY(2)	= XYZ(1) ./ W;      % x = X/(X+Y+Z)
    YXY(3)	= XYZ(2) ./ W;      % y = Y/(X+Y+Z)
    
    if flip
        YXY         = YXY';
    end
else
    
    if ndims(XYZ) < 3
        disp('ERROR: Run with 3-channel image - xyz2yxy');
        return;
    end
   sz				= size(XYZ);
    XYZ				= [reshape(XYZ(:,:,1)',1,sz(1)*sz(2));...
                       reshape(XYZ(:,:,2)',1,sz(1)*sz(2));...
                       reshape(XYZ(:,:,3)',1,sz(1)*sz(2))];
    XYZ				= double(XYZ);		   
    YXYtmp			= zeros(size(XYZ));

    W				= XYZ(1,:) + XYZ(2,:) + XYZ(3,:);
    t1				= find (W > 0);

    YXYtmp(1,t1)	= XYZ(2,t1);			% Y = Y
    YXYtmp(2,t1)	= XYZ(1,t1) ./ W(t1);	% x = X/(X+Y+Z)
    YXYtmp(3,t1)	= XYZ(2,t1) ./ W(t1);	% y = Y/(X+Y+Z)

    YXY(:,:,1)		= reshape (YXYtmp(1,:), sz(2), sz(1))';
    YXY(:,:,2)		= reshape (YXYtmp(2,:), sz(2), sz(1))';
    YXY(:,:,3)		= reshape (YXYtmp(3,:), sz(2), sz(1))';

    if size(YXY) ~= sz
        disp('WARNING: Oops something went wrong with the reshaping!');
        return;
    end

end