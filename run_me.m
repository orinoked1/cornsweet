 close all
clear all

orgim=imread('results_01.jpg');
[img_hight,img_width ] =size(orgim,1:2);
BW1=imread('BWimg1.jpg')>128;

% load('BW1')
im = double(orgim);

%% initial param
max_enhencment_length=50; %50
lambda=log(1+15); %log scale  1.2
sigma=1;%cornsweet 1
maxCPspace=1; % was 3
percent=0.5; %the percent of the maximum cornsweent effect depends on the max luminance around the CP 0.13
orination_block_size = 21; % 99

%% initail calculations
LAB_orig=rgb2lab(im);
LAB_new=LAB_orig;
log_L=log10(LAB_new(:,:,1)+1);

seg=bwlabel(BW1);

%% find edges from seg map
edges=edge(BW1,'canny',0.3);
nedges=bwlabel(edges);
n_edges=max(nedges(:));
for i_edge=1:n_edges
    if sum(nedges(:)==i_edge)<50
        nedges(nedges==i_edge) = 0;
    end
end
n_edges=max(nedges(:));
CP=cell(n_edges,2); % X Y
for  i_edge=1:n_edges
    [cur_y, cur_x]=find(nedges==i_edge);
    [CP{i_edge,1}, CP{i_edge,2}] = findInitialControlPoints( cur_x,cur_y,edges,maxCPspace );
end



%% find the normals control point  diraction is to the area with higher L val
CP_normal = find_normal_directions(log_L,CP,seg,edges,orination_block_size,n_edges);
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
    CP_end_neg{i,2} = CP{i,2}-ceil(max_enhencment_length*sind(-CP_normal{i}));
    %% condition (i) image_border
    CP_end_pos{i,1}(CP_end_pos{i,1}>img_width) =img_width;
    CP_end_pos{i,2}(CP_end_pos{i,2}>img_hight) =img_hight;
    CP_end_pos{i,2}(CP_end_pos{i,2}<1) = 1;
    CP_end_neg{i,1}(CP_end_neg{i,1}>img_width) =img_width;
    CP_end_neg{i,2}(CP_end_neg{i,2}>img_hight) =img_hight;
    CP_end_neg{i,2}(CP_end_neg{i,2}<1) = 1;
    % section (v) two rayes meet.
    dist_CP = squareform(pdist ([CP{i,1};CP{i,2}]'));
    for j=1: size(CP{i,2},2)
        close_CP_idx = dist_CP(j,:)<max_enhencment_length*2;
        for s=find(close_CP_idx)
            if s~=j
                p2p = findintersection([[CP{i,1}(j),CP{i,2}(j)];[CP_end_pos{i,1}(j),CP_end_pos{i,2}(j)]]...
                    ,[[CP{i,1}(s),CP{i,2}(s)];[CP_end_pos{i,1}(s),CP_end_pos{i,2}(s)]]);
                n2n = findintersection([[CP{i,1}(j),CP{i,2}(j)];[CP_end_neg{i,1}(j),CP_end_neg{i,2}(j)]]...
                    ,[[CP{i,1}(s),CP{i,2}(s)];[CP_end_neg{i,1}(s),CP_end_neg{i,2}(s)]]);
                if p2p~=inf
                    CP_end_pos{i,1}(j)=p2p(1);
                    CP_end_pos{i,2}(j)=p2p(2);
                    CP_end_pos{i,1}(s)=p2p(1);
                    CP_end_pos{i,2}(s)=p2p(2);
                end
                if n2n~=inf
                    CP_end_neg{i,1}(j)=n2n(1);
                    CP_end_neg{i,2}(j)=n2n(2);
                    CP_end_neg{i,1}(s)=n2n(1);
                    CP_end_neg{i,2}(s)=n2n(2);
                end
            end
        end
    end
end

%% find Z starting Apmlitude
% get image texturenes
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
        max_neb_L=log_L(ypos,xpos);
        if (CP_Z_pos{i,1}(j)>percent*10^max_neb_L)
            CP_Z_pos{i,1}(j)=log10(percent*10^max_neb_L);
        end
        CP_Z_neg{i,1}(j)=lambda*Tn;
        xneg=CP{i,1}(j)-ceil(cosd(-CP_normal{i}(j)));
        yneg=CP{i,2}(j)-ceil(sind(-CP_normal{i}(j)));
        min_neb_L=log_L(yneg,xneg);
        if(CP_Z_neg{i,1}(j)>percent*10^min_neb_L)
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
        D=size(CP_line_pos{i,1}{j},2);
        Dneg=size(CP_line_neg{i,1}{j},2);
        if(D>1)
            x_ind{i,j}=0:D-1;
            y_ind{i,j}=exp(-2*sigma*(x_ind{i,j})/D);
            pos_enhencment{i,j}=(10^CP_Z_pos{i,1}(j))*(imadjust(y_ind{i,j},[min(y_ind{i,j}),max(y_ind{i,j})],[0 max(y_ind{i,j})]));
            Zplus_effect(...
                sub2ind(size(Zplus_effect),CP_line_pos{i,1}{j},CP_line_pos{i,2}{j}))=...
                +pos_enhencment{i,j};
            
        end
        if(Dneg>1)
            x_ind_neg{i,j}=0:Dneg-1;
            y_ind_neg{i,j}=exp(-2*sigma*(x_ind_neg{i,j})/Dneg);
            neg_enhencment{i,j}=(10^CP_Z_neg{i,1}(j))*(imadjust(y_ind_neg{i,j},[min(y_ind_neg{i,j}),max(y_ind_neg{i,j})],[0 max(y_ind_neg{i,j})]));
            Zneg_effect(...
                sub2ind(size(Zplus_effect),CP_line_neg{i,1}{j},CP_line_neg{i,2}{j}))...
                = -neg_enhencment{i,j};
            
        end
        
    end
end

%% filtering
H=1/25*ones(5,5);
BLplus = imfilter(Zplus_effect,H);
BLneg=imfilter(Zneg_effect,H);
blur_Z=BLplus+BLneg;
Lnew=10.^log_L+blur_Z';


LAB_new(:,:,1)=Lnew;


newim=lab2rgb(LAB_new);
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

