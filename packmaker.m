function packinfo = packmaker(pack,E)
% This function will create an out put struct file to be
% used in the rcfGUI as the recipe file.
% E is a vector of the "evenly spaced" energy values in MeV
% pack is a cell array of strings indicating the pack
% construction
% Vector coordinates to draw representative pack layers
x=[0,0,0,0;.1,.1,.1,.1;0,.1,.1,0;0,.1,.1,0;0,.1,.1,0;0,.1,.1,0]';
y=[0,5,5,0;0,5,5,0;0,0,5,5;0,0,5,5;0,0,0,0;5,5,5,5]'/3;
z=[0,0,5,5;0,0,5,5;0,0,0,0;5,5,5,5;0,0,5,5;0,0,5,5]'/3;
thick = 0;cscale = 5/3+.01;
% Colormaps to paint pack components
RBG = orderedcolors("gem");
maps = ones(1,3,3);
maps(:,:,1) = RBG(1,:); % HD-V2
maps(:,:,2) = [0.65, 0.65, 0.65]; % Al
maps(:,:,3) = RBG(2,:); % EBT3
ptypes = {'hd','al','eb'};
types = zeros(size(ptypes));
dt = [.2,.2,.3,.3,.3,.15,.4,.3,.15,.15,.15,.15,.15];
st = [1,1,2,2,2,.5,3,2,.5,.5,.5,.5,.5];
for n = 1:length(ptypes)
    types(n) = n*(sum(strncmpi(pack,ptypes{n},2))>0);
end
types = types(types>0);
map = maps(:,:,types(1));
for n = 2:length(types)
    map = [map;maps(:,:,types(n))];
end
clf, set(gcf,'renderer','painters'),hold all
for n = 1:length(pack)
    if (strncmpi(pack{n},'al',2)) && ~isempty(pack{n}(3:end))
        althick = str2double(pack{n}(3:end));
        if althick >= 5000
            x3 = x*8 + thick;
            thick = thick + .9;
        elseif althick >= 4000
            x3 = x*7 + thick;
            thick = thick + .8;
        elseif althick >= 3000
            x3 = x*6 + thick;
            thick = thick + .7;
        elseif althick >= 2000
            x3 = x*5 + thick;
            thick = thick + .6;
        elseif althick >= 1500
            x3 = x*4 + thick;
            thick = thick + .5;
        elseif althick >= 1000
            x3 = x*3 + thick;
            thick = thick + .4;
        elseif althick >= 500
            x3 = x*2 + thick;
            thick = thick + .3;
        elseif althick >=100
            x3 = x*1 + thick;
            thick = thick + .2;
        else
            x3 = x*st(find(strncmpi(ptypes,pack{n},2))) + thick;
            thick = thick + dt(find(strncmpi(ptypes,pack{n},2)));
        end
    else
        x3 = x*st(find(strncmpi(ptypes,pack{n},2))) + thick;
        thick = thick + dt(find(strncmpi(ptypes,pack{n},2)));
    end
    y3 = y;
    z3 = z;
    c3 = z + find(types==find(strncmpi(ptypes,pack{n},2)))*cscale;
    fill3(x3,y3,z3,c3)
end
Lbl = {'HD-V2','Al','EBT3'};
for n=1:length(types)
    fill3([0,.1,.1,0],[5/3,5/3,5/3,5/3],[2,2,2.1,2.1]+.2*(n-1),map(n,:));
    text(.1,5/3-.1,2.1+.2*(n-1),Lbl{types(n)})
end
hold off
colormap(map),axis image, axis off, axis tight,view([-22,22])
% If energy value is given the call filmpack.m to solve
if nargin < 2
    packinfo.Pack = pack;
    packinfo.ERange = [];
    packinfo.Result = [];
    return
else
    pause(1)
    for n = 1:length(E)
        D(n,:) = filmpack(E(n),pack);
    end
    D = [E',D];
    packinfo.Pack = pack;
    packinfo.ERange = [min(E),max(E)];
    packinfo.Result = D;
end