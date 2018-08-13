function dn = dateasstring(type)
% dn = dateasstring(type)
% 
% Type:
%  0 = date and time (default)
%  1 = date only

if nargin == 0
    type = 0;
end
d=datevec(now);
d1 = num2str(d(1));
d2 = num2str(d(2));
d3 = num2str(d(3));
d4 = num2str(d(4));
d5 = num2str(d(5));
d6 = num2str(floor(d(6)));
if d(2)<10, d2 = ['0' d2]; end
if d(3)<10, d3 = ['0' d3]; end
if d(4)<10, d4 = ['0' d4]; end
if d(5)<10, d5 = ['0' d5]; end
if d(6)<10, d6 = ['0' d6]; end
dn = [d1 d2 d3];
if type == 0
    dn = [dn d4 d5 d6];
end