function [ CX,CY ] = findInitialControlPoints( X,Y,edges,maxCPspace )

maxSpacing = min(maxCPspace,floor(numel(X)/8)); minSpacing = floor(maxSpacing*.6);

%% Find points with high curvature
% compute local curvature
curvature = zeros(numel(X)-4,1);
for i=11:numel(X)-10
    b1 = (X(i-2)/12)+(X(i-1)/6)-(X(i)/2)+(X(i+1)/6)+(X(i+2)/12);
    b2 = (Y(i-2)/12)+(Y(i-1)/6)-(Y(i)/2)+(Y(i+1)/6)+(Y(i+2)/12);
    c1 = ((X(i+1)-X(i-1))/3)+((X(i+2)-X(i-2))/12);
    c2 = ((Y(i+1)-Y(i-1))/3)+((Y(i+2)-Y(i-2))/12);
    curvature(i-2) = abs(2*((c1*b2-c2*b1)/(c1^2+c2^2)^1.5));
end
 
% displacement
displ = zeros(numel(X)-4,2);
for i=11:numel(X)-10
    displ(i-2,:) = [(X(i-2)/36)+(2*X(i-1)/9)-(X(i)/2)+(2*X(i+1)/9)+(X(i+2)/36)...
        (Y(i-2)/36)+(2*Y(i-1)/9)-(Y(i)/2)+(2*Y(i+1)/9)+(Y(i+2)/36)];
end

t = false(numel(X),1);
t(3:numel(X)-2) = sqrt(displ(:,1).^2+displ(:,2).^2)>=1.5 | curvature>=.8; %| ...

%% for a high-curvature point, find the point with the highest curvature in
% its neighbourhood, fill uniformly if curvature is low
CX = []; CY = [];
trials = [X(t) Y(t)];
curvs = curvature(t(3:numel(X)-2));

tp = [X(1) Y(1)];
if numel(X) > 50
    counter1 = 10; limit = 10;
else
    counter1 = 1; limit = 0;
end
while counter1 <= numel(X)-limit
    point = [X(counter1) Y(counter1)];
    if ~isempty(trials) && ismember(point,trials(1,:),'rows')
        tmp = trials(1,:);
        counter2 = 2;
        while counter2 <= size(trials,1)
            if pdist([trials(1,:);trials(counter2,:)]) <= maxSpacing
                tmp(end+1,:) = trials(counter2,:);
                counter2 = counter2 + 1;
            else
                break;
            end
        end
        [~,M] = max(curvs(1:counter2-1));
        if ~isempty(CX) && pdist([trials(M,:);CX(end) CY(end)]) < minSpacing
            CX(end) = []; CY(end) = [];
        end
        CX(end+1) = trials(M,1); CY(end+1) = trials(M,2);
        [~,counter1] = ismember(trials(M,:),[X Y],'rows');
        tp = trials(M,:);
        trials(1:counter2-1,:) = [];
        curvs(1:counter2-1,:) = [];
    elseif pdist([tp;point]) > maxSpacing
        CX(end+1) = point(1); CY(end+1) = point(2);
        tp = [point(1) point(2)];
    end
    counter1 = counter1 + 1;
end

end

