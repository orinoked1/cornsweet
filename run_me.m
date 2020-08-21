close all
clear all

orgim=imread('1.tif');
[img_hight,img_width ] =size(orgim,1:2);
%imageSegmenter(im);
load('BW1')
im = double(orgim);

%% initial param 
max_enhencment_length=50; %50
lambda=1.2; %log scale  1.2
sigma=3;%cornsweet 1
maxCPspace=0; % was 3
percent=0.13; %the percent of the maximum cornsweent effect depends on the max luminance around the CP 0.13
orination_block_size = 21; % 99

%% initail calculations
LAB_orig=rgb2lab(im);
LAB_new=LAB_orig;
L=log10(LAB_new(:,:,1)+1);

seg=bwlabel(BW1);

%% find edges from seg map
edges=edge(BW1,'canny',0.3);
nedges=bwlabel(edges);
n_edges=max(nedges(:));
CP=cell(n_edges,2); % X Y
for  i_edge=1:n_edges
    [cur_y, cur_x]=find(nedges==i_edge);
    [CP{i_edge,1}, CP{i_edge,2}] = findInitialControlPoints( cur_x,cur_y,edges,maxCPspace );
end



%% find the normals control point  diraction is to the area with higher L val
CP_normal = find_normal_directions(L,CP,seg,edges,orination_block_size,n_edges);
%% for each CP find the enhencment extent acording to end conditions
CP_end_pos=cell(n_edges,2);
CP_end_neg=cell(n_edges,2);
CP_end_pos2=cell(n_edges,2);
CP_end_neg2=cell(n_edges,2);

for  i=1:n_edges 
    %% condition (iv) max_enhencment_length
    CP_end_pos{i,1} = CP{i,1}+ceil(max_enhencment_length*cosd(-CP_normal{i}));
    CP_end_pos{i,2} = CP{i,2}+ceil(max_enhencment_length*sind(-CP_normal{i}));
    CP_end_neg{i,1} = CP{i,1}-ceil(max_enhencment_length*cosd(-CP_normal{i}));
    CP_end_neg{i,2} = CP{i,1}-ceil(max_enhencment_length*sind(-CP_normal{i}));
    %% condition (i) image_border
    CP_end_pos{i,1}(CP_end_pos{i,1}>img_width) =img_width;
    CP_end_pos{i,2}(CP_end_pos{i,2}>img_hight) =img_hight;
    CP_end_pos{i,2}(CP_end_pos{i,2}<1) = 1;
    CP_end_neg{i,1}(CP_end_neg{i,1}>img_width) =img_width;
    CP_end_neg{i,2}(CP_end_neg{i,2}>img_hight) =img_hight;
    CP_end_neg{i,2}(CP_end_neg{i,2}<1) = 1;
   % section (v) ri hits the same enhancement B-spline edge Qi lies on.
    for j=1: size(CP{i,2},2) 
        for s=1:size(CP{i,2},2)
            if s~=j
                pos=    linelineinters([CP{i,1}(j),CP{i,2}(j)],[CP_end_pos{i,1}(j),CP_end_pos{i,2}(j)],[CP{i,1}(s),CP{i,2}(s)],[CP_end_pos{i,1}(s),CP_end_pos{i,2}(s)]);
                neg=    linelineinters([CP{i,1}(j),CP{i,2}(j)],[CP_end_neg{i,1}(j),CP_end_neg{i,2}(j)],[CP{i,1}(s),CP{i,2}(s)],[CP_end_neg{i,1}(s),CP_end_neg{i,2}(s)]);
                negpos= linelineinters([CP{i,1}(j),CP{i,2}(j)],[CP_end_pos{i,1}(j),CP_end_pos{i,2}(j)],[CP{i,1}(s),CP{i,2}(s)],[CP_end_neg{i,1}(s),CP_end_neg{i,2}(s)]);
                if pos~=inf
                    CP1dist=(pos(1)-CP{i,1}(j))^2+(pos(2)-CP{i,2}(j))^2;
                    CP2dist=(pos(1)-CP{i,1}(s))^2+(pos(2)-CP{i,2}(s))^2;
                    if ( CP1dist>CP2dist)
                        CP_end_pos{i,1}(j)=pos(1);
                        CP_end_pos{i,2}(j)=pos(2);
                    else 
                        CP_end_pos{i,1}(s)=pos(1);
                        CP_end_pos{i,2}(s)=pos(2);
                    end  
                end
                if neg~=inf
                    CP1dist=(neg(1)-CP{i,1}(j))^2+(neg(2)-CP{i,2}(j))^2;
                    CP2dist=(neg(1)-CP{i,1}(s))^2+(neg(2)-CP{i,2}(s))^2;
                    if ( CP1dist>CP2dist)
                        CP_end_neg{i,1}(j)=neg(1);
                        CP_end_neg{i,2}(j)=neg(2);
                    else 
                        CP_end_neg{i,1}(s)=neg(1);
                        CP_end_neg{i,2}(s)=neg(2);
                    end
                end
                if negpos~=inf
                        CP1dist=(negpos(1)-CP{i,1}(j))^2+(negpos(2)-CP{i,2}(j))^2;
                        CP2dist=(negpos(1)-CP{i,1}(s))^2+(negpos(2)-CP{i,1}(s))^2;
                        if ( CP1dist>CP2dist)
                            CP_end_pos{i,1}(j)=negpos(1);
                            CP_end_pos{i,2}(j)=negpos(2);
                        else 
                            CP_end_neg{i,1}(s)=negpos(1);
                            CP_end_neg{i,2}(s)=negpos(2);
                        end
                end     
            end
        end
    end
    for j=1: size(CP{i,2},2)
        for s=1:size(CP{i,2},2)
            pos=linelineinters([CP{i,1}(j),CP{i,2}(j)],[CP_end_pos{i,1}(j),CP_end_pos{i,2}(j)],[CP{i,1}(s),CP{i,2}(s)],[CP_end_pos{i,1}(s),CP_end_pos{i,2}(s)]);
            neg=linelineinters([CP{i,1}(j),CP{i,2}(j)],[CP_end_neg{i,1}(j),CP_end_neg{i,2}(j)],[CP{i,1}(s),CP{i,2}(s)],[CP_end_neg{i,1}(s),CP_end_neg{i,2}(s)]);
            if pos~=inf
                CP_end_pos{i,1}(j)=pos(1);
                CP_end_pos{i,2}(j)=pos(2);
                CP_end_pos{i,1}(s)=pos(1);
                CP_end_pos{i,2}(s)=pos(2);    
            end
            if neg~=inf
                 CP_end_neg{i,1}(j)=neg(1);
                 CP_end_neg{i,2}(j)=neg(2);
                 CP_end_neg{i,1}(s)=neg(1);
                 CP_end_neg{i,2}(s)=neg(2);
            end
        end
    end       
end 
for i=1:n_edges
    for j=1:n_edges
        if i~=j
            for ii=1: size(CP{i,2},2)
                for jj=1: size(CP{j,2},2)
                    pos=linelineinters([CP{i,1}(ii),CP{i,2}(ii)],[CP_end_pos{i,1}(ii),CP_end_pos{i,2}(ii)],[CP{j,1}(jj),CP{j,2}(jj)],[CP_end_pos{j,1}(jj),CP_end_pos{j,2}(jj)]);
                    neg=linelineinters([CP{i,1}(ii),CP{i,2}(ii)],[CP_end_neg{i,1}(ii),CP_end_neg{i,2}(ii)],[CP{j,1}(jj),CP{j,2}(jj)],[CP_end_neg{j,1}(jj),CP_end_neg{j,2}(jj)]);
                    negpos=linelineinters([CP{i,1}(ii),CP{i,2}(ii)],[CP_end_neg{i,1}(ii),CP_end_neg{i,2}(ii)],[CP{j,1}(jj),CP{j,2}(jj)],[CP_end_pos{j,1}(jj),CP_end_pos{j,2}(jj)]);
                   %posneg=linelineinters([CP{i,1}(ii),CP{i,2}(ii)],[CP_end_pos{i,1}(ii),CP_end_pos{i,2}(ii)],[CP{j,1}(jj),CP{j,2}(jj)],[CP_end_neg{j,1}(jj),CP_end_neg{j,2}(jj)]);
                    if pos~=inf
                        CP1dist=(pos(1)-CP{i,1}(ii))^2+(pos(2)-CP{j,2}(jj))^2;
                        CP2dist=(pos(1)-CP{i,1}(ii))^2+(pos(2)-CP{j,2}(jj))^2;
                        if ( CP1dist>CP2dist)
                            CP_end_pos{i,1}(ii)=pos(1);
                            CP_end_pos{i,2}(ii)=pos(2);
                        else 
                            CP_end_pos{j,1}(jj)=pos(1);
                            CP_end_pos{j,2}(jj)=pos(2);
                        end
                        %{
                        CP_end_pos{i,1}(ii)=pos(1);
                        CP_end_pos{i,2}(ii)=pos(2);
                        CP_end_pos{j,1}(jj)=pos(1);
                        CP_end_pos{j,2}(jj)=pos(2);
                        %}
                    end
                    if neg~=inf
                        CP1dist=(neg(1)-CP{i,1}(ii))^2+(neg(2)-CP{i,2}(ii))^2;
                        CP2dist=(neg(1)-CP{j,1}(jj))^2+(neg(2)-CP{j,2}(jj))^2;
                        if ( CP1dist>CP2dist)
                            CP_end_neg{i,1}(ii)=neg(1);
                            CP_end_neg{i,2}(ii)=neg(2);
                        else 
                            CP_end_neg{j,1}(jj)=neg(1);
                            CP_end_neg{j,2}(jj)=neg(2);
                        end
                       %{
                        CP_end_neg{i,1}(ii)=neg(1);
                        CP_end_neg{i,2}(ii)=neg(2);
                        CP_end_neg{j,1}(jj)=neg(1);
                        CP_end_neg{j,2}(jj)=neg(2);
                        %}
                    end
                    if negpos~=inf
                        CP1dist=(negpos(1)-CP{i,1}(ii))^2+(negpos(2)-CP{i,2}(ii))^2;
                        CP2dist=(negpos(1)-CP{j,1}(jj))^2+(negpos(2)-CP{j,2}(jj))^2;
                        if ( CP1dist>CP2dist)
                            CP_end_neg{i,1}(ii)=negpos(1);
                            CP_end_neg{i,2}(ii)=negpos(2);
                        else 
                            CP_end_pos{j,1}(jj)=negpos(1);
                            CP_end_pos{j,2}(jj)=negpos(2);
                        end
                        %{
                        CP_end_neg{i,1}(ii)=negpos(1);
                        CP_end_neg{i,2}(ii)=negpos(2);
                        CP_end_pos{j,1}(jj)=negpos(1);
                        CP_end_pos{j,2}(jj)=negpos(2);
                        %}
                    end
                    
                end
            end
        end
    end
end%(ii) ri hits the current medial axis Ai
for i=1:n_edges
    for j=1:n_edges
        if i~=j
            for ii=1: size(CP{i,2},2)
                for jj=1: size(CP{j,2},2)
                    pos=linelineinters([CP{i,1}(ii),CP{i,2}(ii)],[CP_end_pos{i,1}(ii),CP_end_pos{i,2}(ii)],[CP{j,1}(jj),CP{j,2}(jj)],[CP_end_pos{j,1}(jj),CP_end_pos{j,2}(jj)]);
                    neg=linelineinters([CP{i,1}(ii),CP{i,2}(ii)],[CP_end_neg{i,1}(ii),CP_end_neg{i,2}(ii)],[CP{j,1}(jj),CP{j,2}(jj)],[CP_end_neg{j,1}(jj),CP_end_neg{j,2}(jj)]);
                    if pos~=inf
                        CP_end_pos{i,1}(ii)=pos(1);
                        CP_end_pos{i,2}(ii)=pos(2);
                        CP_end_pos{j,1}(jj)=pos(1);
                        CP_end_pos{j,2}(jj)=pos(2);
                    end
                    if neg~=inf  
                        CP_end_neg{i,1}(ii)=neg(1);
                        CP_end_neg{i,2}(ii)=neg(2);
                        CP_end_neg{j,1}(jj)=neg(1);
                        CP_end_neg{j,2}(jj)=neg(2);
                    end
                end
            end 
        end
    end 
end

%% Reduce extent distance var
%update distance CP_pos_ext is the distance for the positive normal
%CP_neg_ext for the negative normal

%[CP_end_pos,CP_end_neg]=variancereduce(CP_pos_ext,CP_neg_ext,CP,CP_end_pos,CP_end_neg,CP_normal,n);
%[CP_end_pos,CP_end_neg] = difffromavg(CP_pos_ext,CP_neg_ext,CP,CP_end_pos,CP_end_neg,CP_normal,n);

%% find Z starting Apmlitude 
T=TexturenessBae(LAB_new(:,:,1)/max(max(LAB_new(:,:,1)))); % the values is between 0 and 1 
T=1+(T/5);

CP_Z_pos=cell(n_edges,1);
CP_Z_neg=cell(n_edges,1);
for i=1:n_edges
    for j=1:size(CP{i,2},2)
        Tn=mean(mean(T(CP{i,2}(j)-2:CP{i,2}(j)+2,CP{i,1}(j)-2:CP{i,1}(j)+2)));
        CP_Z_pos{i,1}(j)=(lambda)*Tn;  
        xpos=CP{i,1}(j)+ceil(cosd(-CP_normal{i}(j)));
        ypos=CP{i,2}(j)+ceil(sind(-CP_normal{i}(j)));
        max_neb_L=L(ypos,xpos);  
        %max_neb_L=max(max(L(CP{i,2}(j)-Area:CP{i,2}(j)+Area,CP{i,1}(j)-Area:CP{i,1}(j)+Area)));
        if (10^CP_Z_pos{i,1}(j)>percent*10^max_neb_L)
            CP_Z_pos{i,1}(j)=log10(percent*10^max_neb_L);
        end
        CP_Z_neg{i,1}(j)=lambda*Tn;
        xneg=CP{i,1}(j)-ceil(2*cosd(-CP_normal{i}(j)));
        yneg=CP{i,2}(j)-ceil(2*sind(-CP_normal{i}(j)));
        min_neb_L=L(yneg,xneg);
        %min_neb_L=min(min(L(CP{i,2}(j)-Area:CP{i,2}(j)+Area,CP{i,1}(j)-Area:CP{i,1}(j)+Area)));
        if(10^CP_Z_neg{i,1}(j)>percent*10^min_neb_L)
            CP_Z_neg{i,1}(j)=log10(percent*10^min_neb_L);
        end
    end
end

%% Cornsweet effect lines depends on Segma
CP_line_pos=cell(n_edges,2);
CP_line_neg=cell(n_edges,2);
for i=1:n_edges %find the cordinates of the point that will get effected 
    for j=1:size(CP{i,2},2)
        CP1dist=sqrt((CP{i,1}(j)-CP_end_pos{i,1}(j))^2+(CP{i,2}(j)-CP_end_pos{i,2}(j))^2);
        [CP_line_pos{i,1}{j},CP_line_pos{i,2}{j}]=fillline([CP{i,1}(j),CP{i,2}(j)],[CP_end_pos{i,1}(j),CP_end_pos{i,2}(j)],floor(CP1dist));
        CP_line_pos{i,1}{j}=round(CP_line_pos{i,1}{j});
        CP_line_pos{i,2}{j}=round(CP_line_pos{i,2}{j});
        CP1dist=sqrt((CP{i,1}(j)-CP_end_neg{i,1}(j))^2+(CP{i,2}(j)-CP_end_neg{i,2}(j))^2);
        [CP_line_neg{i,1}{j},CP_line_neg{i,2}{j}]=fillline([CP{i,1}(j),CP{i,2}(j)],[CP_end_neg{i,1}(j),CP_end_neg{i,2}(j)],floor(CP1dist));
        CP_line_neg{i,1}{j}=round(CP_line_neg{i,1}{j});
        CP_line_neg{i,2}{j}=round(CP_line_neg{i,2}{j});
    end
end
 Zplus_effect=zeros(size(seg,2),size(seg,1));
 Zneg_effect=zeros(size(seg,2),size(seg,1));
for i=1:n_edges % fill the values to the effected points 
    for j=1:size(CP{i,2},2)
        if(CP_Z_pos{i,1}(j)~=0)
             D=size(CP_line_pos{i,1}{j},2);
             Dneg=size(CP_line_neg{i,1}{j},2);
             if(D>1)
                x_ind{i,j}=0:D-1;
                y_ind{i,j}=exp(-2*sigma*(x_ind{i,j})/D);
                yy_ind{i,j}=(10^CP_Z_pos{i,1}(j))*(imadjust(y_ind{i,j},[min(y_ind{i,j}),max(y_ind{i,j})],[0 max(y_ind{i,j})]));
                for s=1:D
                  Zplus_effect(CP_line_pos{i,1}{j}(s),CP_line_pos{i,2}{j}(s))=+yy_ind{i,j}(s);
                end
             end
             if(Dneg>1)
                x_ind_neg{i,j}=0:Dneg-1;
                y_ind_neg{i,j}=exp(-2*sigma*(x_ind_neg{i,j})/Dneg);
                yy_ind_neg{i,j}=(10^CP_Z_neg{i,1}(j))*(imadjust(y_ind_neg{i,j},[min(y_ind_neg{i,j}),max(y_ind_neg{i,j})],[0 max(y_ind_neg{i,j})]));
                for s=1:Dneg
                  Zneg_effect(CP_line_neg{i,1}{j}(s),CP_line_neg{i,2}{j}(s))=-yy_ind_neg{i,j}(s);
                end
             end
        end
    end 
end
        
%% blur and conv
%H=3*fspecial('gaussian',3,1);
H=1/25*ones(5,5);
BLplus = imfilter(Zplus_effect,H);
BLneg=imfilter(Zneg_effect,H);
%mask = imdilate(edges,ones(maxextention));
%blur_Zplus(mask==1) = BLplus(mask==1);
%blur_Zneg(mask==1)=BLneg(mask==1);
blur_Z=BLplus+BLneg;
%blur_Z=Zplus_effect+Zneg_effect;
blur_Z=imresize(blur_Z,[size(L,2),size(L,1)]);
Lnew=10.^L+blur_Z';


%BL = imfilter(Lnew,H);
%mask = imdilate(edges,ones(maxextention));
%Lnew(mask==1) = BL(mask==1);
LAB_new(:,:,1)=Lnew;


newim=lab2rgb(LAB_new);
% Resim=uint8(255*(newim./255).^(1/2.2));
Resim = uint8(newim);

%% plotting
%ploting the extend 

figure; plot(x_ind{i,j},y_ind{i,j})
hold on;
plot(-x_ind{i,j},-y_ind{i,j})

figure;
hold on;
for i=1:n_edges
     for ii=1: size(CP{i,2},2)
         plot([CP_end_neg{i,1}(ii),CP{i,1}(ii)],[CP_end_neg{i,2}(ii),CP{i,2}(ii)],'r')
         plot([CP_end_pos{i,1}(ii),CP{i,1}(ii)],[CP_end_pos{i,2}(ii),CP{i,2}(ii)],'b')
         plot(CP_line_pos{i,1}{ii},CP_line_pos{i,2}{ii},'.m')
         plot(CP_line_neg{i,1}{ii},CP_line_neg{i,2}{ii},'.y')
     end 
end
 
% plot the normals 
figure;imshow(edges,[]);
%Overlay normal
hold on;
for i=1:n_edges
    scatter(CP{i,1},CP{i,2},'*b')
    quiver(CP{i,1},CP{i,2},-cosd(CP_normal{i,1}+180),sind(CP_normal{i,1}+180),'r')
    quiver(CP{i,1},CP{i,2},-cosd(CP_normal{i,1}),sind(CP_normal{i,1}),'y')
end
hold off;
%plot the Effect that was added to the luminance
figure;imshow(blur_Z',[])

%plot the results
figure;imshow(Resim,[]),title('after')
figure;imshow(orgim,[]),title('before')

%plot the luminance effect at line between for [x1,x2]
figure;plot(Lnew(350,:,1)), hold on;plot(LAB_orig(350,:,1)),legend('after','before')

