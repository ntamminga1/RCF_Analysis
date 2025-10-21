function [EDose,xrange,EnergyINS]=filmpack(EIN,PackComp)
% Compute the absorbed dose each active layer of film for
% a proton with energy EIN (MeV).
% -------- In microns ---------- DEFAULTS
ALTHICK = 25;
HDGELTHICK = 0.75;
HDACTIVETHICK = 10;
HDBASETHICK = 97;
EBACTIVETHICK = 17;
EBBASETHICK = 97;
EBSURFTHICK = 6;
dEper = .005;
minE = .05;
% ------------------------------
if nargin < 2
    PackComp = {'HD','EB'};
end
Ion = 'H';
EDose=0; xrange=0; EnergyINS=0;
% Get dE/dx (MeV/mm) for all materials used SRIM
[E,dEdx] = ExtractSRIM([Ion,' in Aluminum']);
sdata(1,:) = {E,dEdx};
[E,dEdx] = ExtractSRIM([Ion,' in Active Layer']);
sdata(2,:) = {E,dEdx};
[E,dEdx] = ExtractSRIM([Ion,' in Polyester Base']);
sdata(3,:) = {E,dEdx};
[E,dEdx] = ExtractSRIM([Ion,' in Gelatin']);
sdata(4,:) = {E,dEdx};
[E,dEdx] = ExtractSRIM([Ion,' in EBTACTIVE']);
sdata(5,:) = {E,dEdx};
[E,dEdx] = ExtractSRIM([Ion,' in EBTSURF']);
sdata(6,:) = {E,dEdx};

% Create Thickness/Type Arrays
clear TauArray;
clear TypeArray;
IDArray = [];
counter = 0;
for n = 1:length(PackComp)
    layerstring = PackComp{n};
    if strcmpi(layerstring(1:2),'Al')
        counter = counter + 1;
        if length(layerstring)==2
            TauArray(counter) = ALTHICK;
        else
            TauArray(counter) = ...
                str2num(layerstring(3:end));
        end
        TypeArray(counter) = 1;
        IDArray = [IDArray,n];
    elseif strcmpi(layerstring,'HD')
        counter = counter + 3;
        TauArray(counter-2) = HDGELTHICK;
        TauArray(counter-1) = HDACTIVETHICK;
        TauArray(counter) = HDBASETHICK;
        TypeArray(counter-2) = 4;
        TypeArray(counter-1) = 2;
        TypeArray(counter) = 3;
        IDArray = [IDArray,n,n,n];
    elseif strcmpi(layerstring,'EB')
        counter = counter + 5;
        TauArray(counter-4) = EBBASETHICK;
        TauArray(counter-3) = EBACTIVETHICK;
        TauArray(counter-2) = EBSURFTHICK;
        TauArray(counter-1) = EBACTIVETHICK;
        TauArray(counter) = EBBASETHICK;
        TypeArray(counter-4) = 3;
        TypeArray(counter-3) = 5;
        TypeArray(counter-2) = 6;
        TypeArray(counter-1) = 5;
        TypeArray(counter) = 3;
        IDArray = [IDArray,n,n,n,n,n];
    else
        error('Film Type Not Found')
    end
end
stopflag = 0;
travelX = 0; EDose = zeros(length(PackComp),1);
EnergyIN = EIN(1);
EnergyINS = 0;
oldtype = 0;
for n = 1:length(TypeArray)
    newtype = IDArray(n);
    if oldtype ~= newtype
        EnergyINS(IDArray(n)) = EnergyIN;
    end
    oldtype = newtype;
    nt = TypeArray(n);
    xinlayer = 0;
    while xinlayer < TauArray(n)
        if EnergyIN <= 0, stopflag = 1; break, end
        dEdx = pchip(sdata{nt,1},sdata{nt,2},EnergyIN);
        dx = 1000*dEper*EnergyIN/dEdx;
        dx = min(dx,TauArray(n)-xinlayer);
        xinlayer = xinlayer + dx;
        EnergyOUT = EnergyIN - dx*dEdx/1000;
        if EnergyOUT <= minE, EnergyOUT = 0; end
        travelX = travelX + dx;
        EDose(IDArray(n)) = EDose(IDArray(n)) + ...
            (EnergyIN-EnergyOUT)*(nt==2||nt==6);
        EnergyIN = EnergyOUT;
    end
end
EnergyINS = [EnergyINS,EnergyIN];
if stopflag
    xrange = travelX;
else
    xrange = -1;
end

