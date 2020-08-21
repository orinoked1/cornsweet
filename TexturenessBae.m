function [ result ] = TexturenessBae( I )

f = hp(9,1);

Habs = abs(imfilter(I,f,'replicate'));

sigmaRange = 0.5;
sigmaSpatial = round(8*sigmaRange);

result = bilateralFilter(I,Habs,0,1,sigmaSpatial,sigmaRange);

end

