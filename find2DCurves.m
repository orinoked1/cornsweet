% fit curves to discrete edges. Since we did not find any suitable
% MATLAB code for this, we naively implemented a solution.
% This is extremely slow code and is probably sub-optimal

% C: control points (integers), Cd: double (moved to correct step edge)
% A: angles for each control point (i.e. normal)
% O: orientation in the direction of the given angles (i.e. positive halo
%    or negative halo
function [ C,Cd,A,Ap,o ] = find2DCurves( Px,Py,edge,points,L )
% paramters
l = 0; % length of slope check

C = cell(0);
A = cell(0);
% if closed curve, remove the last (duplicate point)
closed = ismember([Px(1) Py(1)],[Px(end) Py(end)],'rows');
if closed
    points(end,:) = []; 
    Px(end) = []; Py(end) = [];
end

%% Find the 2D curve, and their angles
[C, A] = fitCurveRec([Px' Py'],edge,points,C,A);

%% Merge splits (i.e. same angles for end control points)
for i=1:numel(A)-1
    a1 = A{i}(end);
    a2 = A{i+1}(1);
    angle = (a1+a2)/2;
    [X,Y] = CubicBSplineClamped(C{i});
    pointsSP = [X' Y'];
    normals = CubicBSplineClampedDer(C{i});
    pos = findClosestPoint(C{i}(end,:),pointsSP);
    normal = normals(pos,:);
    angle = checkNormalConsistency([cos(angle) sin(angle)],normal);
    % if angle between the new vector and a1 and a2 is not the same, change
    t1 = acos(dot([cos(angle) sin(angle)]./norm([cos(angle) sin(angle)])...
        ,[cos(a1) sin(a1)]./norm([cos(a1) sin(a1)])));
    t2 = acos(dot([cos(angle) sin(angle)]./norm([cos(angle) sin(angle)])...
        ,[cos(a2) sin(a2)]./norm([cos(a2) sin(a2)])));
    if abs(t1-t2) > .0001
        angle = 2*pi-angle;
    end
    
    A{i+1}(1) = angle;
    A{i}(end) = angle;
end

%% Merge end points for closed curves
if closed
    pointsSP = [fliplr(C{1}(2:3,:)')';fliplr(C{end}(end-1:end,:)')'];
    [X,Y] = CubicBSplineClamped(pointsSP);
    normals = CubicBSplineClampedDer(pointsSP);
    jp = C{1}(1,:);
    D = pdist([jp;[X' Y']]);
    D = D(1:numel(X));
    [~,pos] = min(D);
    normal = normals(pos,:);
    if dot([X(pos) Y(pos)]-jp,normal) > 0, normal = -normal; end
        angle = acos(dot([1 0],normal./norm(normal)));
    if abs(normal(1))<.0001 && abs(normal(2)+1)<.0001, angle = -angle; end
    if abs(cos(angle)-normal(1)) > .0001 || abs(sin(angle)-normal(2)) > .0001
        angle = 2*pi-angle;
    end
    
    C{end}(end+1,:) = C{1}(1,:);
    A{1}(1) = angle;
    A{end}(end+1) = angle;
end

%% move CPs to highest slope, in between pixels
dist = -l:l;
Cd = C;
for i=1:numel(A)
    CP = C{i};
    for j=2:size(CP,1)-1
        diff = zeros(size(dist));
        for k=1:numel(dist)
            Ltarget = getLuminance(CP(j,:),A{i}(j),dist(k),L);
            Lother = getLuminance(CP(j,:),A{i}(j),dist(k)+1,L);
            diff(k) = abs(Ltarget-Lother);
        end
        [~,pos] = max(diff);
        CP(j,:) = [CP(j,1)+(dist(pos)+.5)*cos(A{i}(j)),CP(j,2)+(dist(pos)+.5)*sin(A{i}(j))];
    end
    Cd{i} = CP;
end

%% Interpolate CP vectors to the other points
Ap = cell(numel(C),1);
startp = 2;
for i=1:numel(C)
    CP = C{i};
    Acp = cell(size(CP,1)-1,1);
    for j=1:size(CP,1)-1
        [~,endp] = ismember(CP(j+1,:),points,'rows'); endp = endp-1;
        N1 = [cos(A{i}(j)) sin(A{i}(j))];
        N2 = [cos(A{i}(j+1)) sin(A{i}(j+1))];
        d = endp-startp+2;
        data = zeros(endp-startp+1,3);
        counter = 1;
        for k=startp:endp
            alpha = (d-(k-(startp-1)))/d; % weight for first CP
            beta = (d-((endp+1)-k))/d; % weight for second CP
            N = alpha*N1+beta*N2; N = N./norm(N);
            angle = acos(dot([1 0],N));
            if dot([cos(angle) sin(angle)],N1) < 0, angle = -angle; end
            if abs(cos(angle)-N(1)) > .0001 || abs(sin(angle)-N(2)) > .0001
                angle = 2*pi-angle;
            end
            % move point to the actaul step edge
            diff = zeros(size(dist));
            for b=1:numel(dist)
                Ltarget = getLuminance(points(k,:),angle,dist(b),L);
                Lother = getLuminance(points(k,:),angle,dist(b)+1,L);
                diff(b) = abs(Ltarget-Lother);
            end
            [~,pos] = max(diff);
            tmp = [points(k,1)+(dist(pos)+.5)*cos(angle),points(k,2)+(dist(pos)+.5)*sin(angle)];
            data(counter,:) = [tmp angle];
            counter = counter + 1;
        end
        startp = endp+2;
        Acp{j} = data;
    end
    Ap{i} = Acp;
end

%% Find orientation for the edglet
CP = C{1};
ang = A{1};
nsamples = size(CP,1);
ls = zeros(nsamples-2,2);
imsize = size(edge);
for k=2:nsamples-1 % loop through control points
    normal = 5*[cos(ang(k,1)) sin(ang(k,1))];
    p1 = max(min(round(CP(k,:)+normal),imsize(2)),1);
    p2 = max(min(round(CP(k,:)-normal),imsize(1)),1);
    ls(k-1,1) = L(p1(2),p1(1));
    ls(k-1,2) = L(p2(2),p2(1));
end
L1 = mean(ls(:,1)); L2 = mean(ls(:,2));
if L1 > L2
    o = 1;
else
    o = 0;
end

end

%% Recursively fit B spline curves
% C - control points
% A - corresponding angles
% Ep - full edge (points), Et: trimmed edge
function [ C,A ] = fitCurveRec( Cin,Et,Ep,C,A )

T = 1;
imsize = size(Et);

% find the edge point with the largest curve displacement
% this code is poorly written.
[PX,PY] = CubicBSplineClamped(Cin);
curve = false(imsize); % draw the binary curve
for i=1:numel(PX)-1
    p1 = [round(PX(i)),round(PY(i))];
    p2 = [round(PX(i+1)),round(PY(i+1))];
    [~,~,curve] = bresenham(curve,[p1;p2],0,1);
end
D = zeros(size(Ep,1),1);
[PY,PX] = find(curve);
[py,px] = find(Et);
for i=1:numel(px)
    tmp = sqrt((px(i)-PX).^2+(py(i)-PY).^2);
    D(i) = min(tmp);
end
[M,P] = max(D);
if size(Cin,1)>=4 && M > T % recurse if the displacement is large enough
    [~,pos] = ismember([px(P) py(P)],Ep,'rows');
    [~,ci] = ismember(Cin,Ep,'rows'); 
    t = false(numel(ci),1); tmp = find(ci>0); t(tmp(1):end) = true;
    ci(ci==0&t)=size(Ep,1);
    [isC,loc] = ismember(pos,ci);
    if isC(1) % split is already a control point
        split = loc(1);
    else % insert new control point
        [~,loc] = sort([pos(1) ci']);
        split = find(loc==1);
        tmp = Cin;
        tmp(split,:) = [px(P) py(P)]; 
        tmp(split+1:end+1,:) = Cin(split:end,:);
        Cin = tmp;
    end
    % recurse
    [~,locs1] = ismember(Cin(split-3,:),Ep,'rows'); % 3 or 4?
    [~,locs2] = ismember(Cin(split+3,:),Ep,'rows');
    [~,loc1] = ismember(Cin(4,:),Ep,'rows');
    [~,loc2] = ismember(Cin(end-3,:),Ep,'rows');
    El = Et(:); Er = Et(:); El(:) = false; Er(:) = false;
    El(sub2ind(imsize,Ep(loc1+3:locs1-3,2),Ep(loc1+3:locs1-3,1))) = true; 
    Er(sub2ind(imsize,Ep(locs2+3:loc2-3,2),Ep(locs2+3:loc2-3,1))) = true;
    Er = reshape(Er,imsize); El = reshape(El,imsize);
    [C, A] = fitCurveRec(Cin(1:split,:),El,Ep,C,A);
    [C, A] = fitCurveRec(Cin(split:end,:),Er,Ep,C,A);
else % add control points and angles to cell
    C{end+1} = Cin;
    angles = zeros(size(Cin,1),1);
    [X,Y] = CubicBSplineClamped(Cin);
    points = [X' Y'];
    N = CubicBSplineClampedDer(Cin);
    for i=1:numel(angles) % find angles
        CP = Cin(i,:);
        pos = findClosestPoint(CP,points);
        if pos==numel(X) % find tangent
            T = points(pos,:)-points(pos-1,:);
        else
            T = points(pos+1,:)-points(pos,:);
        end
        tmp = T; % do lots of stupid stuff, checking we have the right angle
        T = [-T(2) T(1)]; T = T./norm(T); % rotate tangent
        normal = N(pos,:)./norm(N(pos,:));
        angle = checkNormalConsistency(T,normal);
        if abs(dot(tmp,[cos(angle) sin(angle)])) > .0001
            angle = 2*pi-angle; 
        end
        angles(i) = angle;
    end
    A{end+1} = angles;
end

end

%% Consistency check
function angle = checkNormalConsistency(NN,normal)
if dot(NN,normal) < 0, NN = -NN; end
angle = acos(dot([1 0],NN));
if abs(NN(1))<.0001 && abs(NN(2)+1)<.0001, angle = -angle; end
end

%% Closest point to cp
function pos = findClosestPoint(cp,points)
D = sqrt((points(:,1)-cp(1)).^2+(points(:,2)-cp(2)).^2);
[~,pos] = min(D);
end

function res = getLuminance(CP,A,d,L)
imsize=size(L);
r = max(min(round(CP(2)+d*sin(A)),imsize(1)),1);
c = max(min(round(CP(1)+d*cos(A)),imsize(2)),1);
res = L(r,c);
end