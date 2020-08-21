function pos = linelineinters(p1,p2,p3,p4)
denom = (p1(1)-p2(1))*(p3(2)-p4(2))-(p1(2)-p2(2))*(p3(1)-p4(1));

if abs(denom) < .0001, pos = inf; return; end
pos = [0 0];
pos(1) = ((p1(1)*p2(2)-p1(2)*p2(1))*(p3(1)-p4(1))-(p1(1)-p2(1))*(p3(1)*p4(2)-p3(2)*p4(1)))/denom;
pos(2) = ((p1(1)*p2(2)-p1(2)*p2(1))*(p3(2)-p4(2))-(p1(2)-p2(2))*(p3(1)*p4(2)-p3(2)*p4(1)))/denom;
if abs(pos(1)-p1(1))<.0001, pos(1)=p1(1); end
if abs(pos(1)-p2(1))<.0001, pos(1)=p2(1); end
if abs(pos(1)-p3(1))<.0001, pos(1)=p3(1); end
if abs(pos(1)-p4(1))<.0001, pos(1)=p4(1); end
if abs(pos(2)-p1(2))<.0001, pos(2)=p1(2); end
if abs(pos(2)-p2(2))<.0001, pos(2)=p2(2); end
if abs(pos(2)-p3(2))<.0001, pos(2)=p3(2); end
if abs(pos(2)-p4(2))<.0001, pos(2)=p4(2); end
if (p1(1)>p2(1) && pos(1)>p1(1))||(p2(1)>p1(1) && pos(1)>p2(1))||(p1(1)<p2(1) && pos(1)<p1(1))||(p2(1)<p1(1) && pos(1)<p2(1))...
        ||(p1(2)>p2(2) && pos(2)>p1(2))||(p2(2)>p1(2) && pos(2)>p2(2))||(p1(2)<p2(2) && pos(2)<p1(2))||(p2(2)<p1(2) && pos(2)<p2(2))...
        ||(p3(1)>p4(1) && pos(1)>p3(1))||(p4(1)>p3(1) && pos(1)>p4(1))||(p3(1)<p4(1) && pos(1)<p3(1))||(p4(1)<p3(1) && pos(1)<p4(1))...
        ||(p3(2)>p4(2) && pos(2)>p3(2))||(p4(2)>p3(2) && pos(2)>p4(2))||(p3(2)<p4(2) && pos(2)<p3(2))||(p4(2)<p3(2) && pos(2)<p4(2))
    pos = inf;
end
end