function CP_normal = find_normal_directions(L,CP,seg,edges,orination_block_size,n_edges)
% find the normals control point  diraction is to the area with higher L val

%% avarage Luminance depends of segemtation area { L(n+1) is the background }
avgL=zeros(n_edges+1,1);
for i=1:n_edges+1
    avgL(i)=mean(L(seg==i-1));
end 

Orientations = skeletonOrientation(edges,orination_block_size); 
Onormal = Orientations+90; 
CP_normal=cell(n_edges,1);
for  i=1:n_edges
    linearidx = sub2ind(size(Onormal),CP{i,2},CP{i,1});
    normal_direction=Onormal(linearidx);
    % change direction by 180 if from high to low 
    linearidx = sub2ind(size(Onormal),CP{i,2}+ceil(10*sind(-normal_direction)),...
        CP{i,1}+ceil(10*cosd(-normal_direction)));
    label_1 = seg(linearidx)+1;
    linearidx = sub2ind(size(Onormal),CP{i,2}-ceil(10*sind(-normal_direction)),...
        CP{i,1}-ceil(10*cosd(-normal_direction)));
    label_2 = seg(linearidx)+1;
    normal_direction(avgL(label_1)<avgL(label_2))  = normal_direction(avgL(label_1)<avgL(label_2))+180;
    CP_normal{i} = normal_direction;
end