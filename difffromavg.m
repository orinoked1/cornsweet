function [CP_end_pos,CP_end_neg] = difffromavg(CP_pos_ext,CP_neg_ext,CP,CP_end_pos,CP_end_neg,CP_normal,n)
for i=1:n
    for j=1: size(CP{i,2},2)
        CP_pos_ext{i,1}(j)=sqrt((CP_end_pos{i,1}(i)-CP{i,1}(i))^2+(CP_end_pos{i,2}(i)-CP{i,2}(j))^2);
        CP_neg_ext{i,1}(j)=sqrt((CP_end_neg{i,1}(i)-CP{i,1}(i))^2+(CP_end_neg{i,2}(i)-CP{i,2}(j))^2); 
    end
    maxD = 15; % max difference
    pos_avg=mean(CP_pos_ext{i,1});
    neg_avg=mean(CP_neg_ext{i,1});
    for j=1:size(CP{i,2},2)
        if(CP_pos_ext{i,1}(j)>pos_avg+maxD)
            CP_pos_ext{i,1}(j)=pos_avg+maxD;
        end
        if(CP_pos_ext{i,1}(j)<pos_avg-maxD)
            CP_pos_ext{i,1}(j)=pos_avg-maxD;
        end
        if(CP_neg_ext{i,1}(j)>neg_avg+maxD)
            CP_neg_ext{i,1}(j)=neg_avg+maxD;
        end
        if(CP_neg_ext{i,1}(j)<neg_avg-maxD)
            CP_neg_ext{i,1}(j)=neg_avg-maxD;
        end 
        CP_end_pos{i,1}(j)=CP{i,1}(j)+ceil(CP_pos_ext{i,1}(j)*cosd(-CP_normal{i}(j)));
        CP_end_pos{i,2}(j)=CP{i,2}(j)+ceil(CP_pos_ext{i,1}(j)*sind(-CP_normal{i}(j)));
        CP_end_neg{i,1}(j)=CP{i,1}(j)-ceil(CP_neg_ext{i,1}(j)*cosd(-CP_normal{i}(j)));
        CP_end_neg{i,2}(j)=CP{i,2}(j)-ceil(CP_neg_ext{i,1}(j)*sind(-CP_normal{i}(j)));
    end
end


