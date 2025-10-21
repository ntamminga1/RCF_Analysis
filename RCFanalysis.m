%% Code for analyzing RCF stack data.
% Used for analysis of Titan experimental data in Feb/Mar 2024
% Written by Nathaniel Tamminga 2025, adapted from Chris Willis and
% Dustin Offerman

% Generate the film pack using the packmaker.m and filmpack.m functions
% This generates the response function of each layer of film for protons
% with energy of 1 to 30 MeV. It then calculates the energy range over
% which the film will have substantial response, the mean of this range is
% output in rowmean.
% Edit the pack configuration as necessary. Current code only handles
% aluminum filtering and HD-V2 and EBT-3 films. You must adjust your SRIM
% directory in ExtractSRIM.m
% You can add films and filtering materials to filmpack.m and
% deltaE.m
% Result is the response function of the active layer in MeV/proton
pack = {'al75','hd','al77','hd','al77','hd','al150','hd','al150','hd','al150','eb','al227','eb','al300','eb','al227','eb','al300','eb'};
Energy = 1:0.01:30;
packinfo = packmaker(pack, Energy);
Result = packinfo.Result;
for n=1:length(Energy)
    dE(n,:) = deltaE(Energy(n),pack);
end
rowmean = sum(dE,1) ./ sum(dE~=0,1);

% Function converts pixel value to dose (Gy) from Willis 2016 for HD-V2 film
function [y] = redFunc(x)
y = -276.6*log(x) + 3017.3;
if y<0
    y=0;
end
if y==inf
    y=0;
end
end

% Change this if you want to save the plots you generate
save = false;

%% #########################################################################
% ######################### 50J - 0ps #####################################
% #########################################################################
% The following section of code is for a set of shots with ~50J of total
% pulse engery at 0ps delay.

% ##################### Load the dose (Gy/sr) #############################
% This section is only necessary for loading in the mean dose in each film
% layer

% These shots had variable stack design, thus the need for the estacks variable.
shots = ["S002","S003","S004","S007","S008","S009","S010","S011","S012","S013"]; % Shots of interest
stacks = ["H1", "H2", "H3", "H4", "H5", "E6", "E7", "E8", "E9", "E10"]; % Film configuration
estacks = [4,4,5,5,5,5,5,5,5,5]; % Number of EBT-3 Films in stack, HD-V2 is constant at 5
stacksize = [2,2,2,2,2,2,2,2,4,4]; % Size of stacks (2 in or 4 in)
% Energies of each layer as determined from ray tracing program
TotEnergy = [2.94, 5.48, 7.43, 9.67, 11.61, 13.85, 16.55, 19.24, 21.3, 23.7];


% Loop through all shots and calculate total dose on film
for k = 1:length(shots)
    shot = shots(k);                % Shot Number of interest
    RCFSize = stacksize(k);         % RCF size (2in or 4in)
    HDV2 = 5;                       % Number of HD-V2 layers
    EBT3 = estacks(k);              % Number of EBT-3 layers
    totstacks = HDV2+EBT3;          % Total stack number
    dosePerSR = NaN(1,10);          % Empty dose array
    rightdosePerSR = NaN(1,10);     % Empty dose on right half
    leftdosePerSR = NaN(1,10);      % Empty dose on left half
    h_index = NaN(1,5);             % Empty HD-V2 array

    % Set steradian value for film size
    if RCFSize==2
        SR = 0.222;
    end
    if RCFSize==4
        SR = 0.769;
    end

    % Load the shot dataStruct (Calculated elsewhere)
    shotStruct = load(strcat("..\RCF_Images\dataStruct\",shot,".mat"));

    % Calculate mean dose on each HD-V2 layer
    for i=1:HDV2
        h_index(i) = find([shotStruct.dataStruct.Label]==strcat('H',num2str(i)));
    end

    for i=1:HDV2
        layer = cast(shotStruct.dataStruct(h_index(i)).Data(:,:,1), 'double');
        conv2dose = arrayfun(@(x) redFunc(x), layer);
        midpoint = size(conv2dose, 2)/2;
        totalDose = mean(conv2dose,'all','omitnan');
        leftdose = mean(conv2dose(:,1:midpoint),'all','omitnan');
        rightdose = mean(conv2dose(:,midpoint:end),'all','omitnan');
        dosePerSR(i) = totalDose/SR;
        rightdosePerSR(i) = rightdose/(SR/2);
        leftdosePerSR(i) = leftdose/(SR/2);
    end

    % Calculate mean dose on each EBT-3 layer
    for j = 6:totstacks
        arr = readmatrix(strcat('..\RCF_Images\dataStruct\csv\Dose\',shot,'_',stacks(j),'.csv'));
        midpoint = size(arr, 2)/2;
        totalDose = mean(arr,'all','omitnan');
        leftdose = mean(arr(:,1:midpoint),'all','omitnan');
        rightdose = mean(arr(:,midpoint:end),'all','omitnan');
        dosePerSR(j) = totalDose/SR;
        rightdosePerSR(j) = rightdose/(SR/2);
        leftdosePerSR(j) = leftdose/(SR/2);
    end

    % End with array of the mean dose in each layer for each shot (Gy/sr)
    doseArr(:,k) = dosePerSR;
    leftDose(:,k) = leftdosePerSR;
    rightDose(:,k) = rightdosePerSR;
end

% ######################### Unit Conversions #############################
% Starting with dose in Gy/sr
doseArr = doseArr*6.2415E9; % Gy/sr to MeV/g/sr
% Dose array is indexed (layer, shot) so doseArr(1,:) is all the H1 layer
% while doseArr(:,1) is all shot "S002"

% Multiply l*p (thickness*density) (g/cm^2) for all H layers (1:5)
doseArr(1:5,:) = doseArr(1:5,:)*10E-4*1.25;
% Same for E layers (6:end)
doseArr(6:end,:) = doseArr(6:end,:)*30E-4*1.1;

% Now we have units of MeV/cm^2/s.r., this is an energy flux per steradian
% Now multiply out the area of the RCF
% This is 2in for all but the final 2 sheets
doseArr(:,1:8)=doseArr(:,1:8)*(2*2.54)^2;
doseArr(:,9:end)=doseArr(:,9:end)*(4*2.54)^2;

% This leaves units of MeV/s.r.
% Now account for the RCF bin width (energy range over which the
% active layer works -> rowmean)
% This is a layer changing variable
for n=1:size(doseArr,1)
    doseArr(n,:)=doseArr(n,:)/rowmean(2*n);
end

% Now in units of 1/s.r.
% Now address the response function found in Result
% The response functions are for each layer H1, H2, etc.
for n=1:size(doseArr,1)
    rn = 2*n+1; % Indexing for Result data
    for k=1:size(doseArr,2) %Have to iterate through each shot for the proper math
        dNdE(n,k)=max(doseArr(n,k)/Result(:,rn));
    end
end
% Dividing by the response funciton of MeV/proton, we get the desired units
% of Proton/MeV/s.r.

% ########################## Plotting Figure ##############################
figure(2)
semilogy(TotEnergy, dNdE(:,1),'o', Color="#0072BD", HandleVisibility="off")
hold on
semilogy(TotEnergy, dNdE(:,2),'o', Color="#0072BD", HandleVisibility="off")
semilogy(TotEnergy, dNdE(:,3),'o', Color="#0072BD", HandleVisibility="off")
semilogy(TotEnergy, dNdE(:,6),'o', Color="#0072BD", HandleVisibility="off")
semilogy(TotEnergy, dNdE(:,10),'o', Color="#0072BD", HandleVisibility="off")

semilogy(TotEnergy, dNdE(:,7),'diamond', Color="#D95319", HandleVisibility="off")
semilogy(TotEnergy, dNdE(:,9),'diamond', Color="#D95319", HandleVisibility="off")

semilogy(TotEnergy, dNdE(:,4),'x', Color="#77AC30", HandleVisibility="off")
semilogy(TotEnergy, dNdE(:,5),'x', Color="#77AC30", HandleVisibility="off")
semilogy(TotEnergy, dNdE(:,8),'x', Color="#77AC30", HandleVisibility="off")

% For the legend
empty = NaN(1,10);
semilogy(TotEnergy, empty,'diamond', Color="#D95319", DisplayName="Beam 1")
semilogy(TotEnergy, empty,'x', Color="#77AC30", DisplayName="Beam 2")
semilogy(TotEnergy, empty,'o', Color="#0072BD", DisplayName="Both Beams")

legend

ylim([10^9 2E13])
xlim([0 25])
set(gca,'fontsize',14);
xlabel("Energy (MeV)")
ylabel("dN/dE/d\Omega (Protons/MeV/sr)")
if save
    savefig(gcf, '..\RCF_Images\50JRCF_dNdE.fig')
    saveas(gcf, '..\RCF_Images\50JRCF_dNdE.png')
end
hold off

% ############################ Error Bars #################################
% Calculate the standard error of each layer for each shot configuration
b2cat = cat(1,dNdE(:,4)',dNdE(:,5)',dNdE(:,8)');
b2mean = mean(b2cat,1);
b2max = std(b2cat)/sqrt(3);
b2min = std(b2cat)/sqrt(3);

b1cat = cat(1,dNdE(:,7)',dNdE(:,9)');
b1mean = mean(b1cat,1);
b1max = std(b1cat)/sqrt(2);
b1min = std(b1cat)/sqrt(2);

bbcat = cat(1,dNdE(:,1)',dNdE(:,2)',dNdE(:,3)',dNdE(:,6)',dNdE(:,10)');
nanarr = sum(isnan(bbcat),1);
bbmean = nanmean(bbcat,1);
bbmax = nanstd(bbcat)./sqrt(5-nanarr);
bbmin = nanstd(bbcat)./sqrt(5-nanarr);

% Plot error bar version of previous plot
figure(3)
semilogy(TotEnergy,bbmean,'o', Color="#0072BD",HandleVisibility="off",LineWidth=1.5)
hold on
semilogy(TotEnergy,b1mean,'diamond', Color="#D95319",HandleVisibility="off")
semilogy(TotEnergy,b2mean,'x', Color="#77AC30",HandleVisibility="off")
errorbar(TotEnergy,bbmean,bbmin,bbmax,"LineStyle","none",Color="#0072BD",HandleVisibility="off",LineWidth=1.5)
errorbar(TotEnergy,b1mean,b1min,b1max,"LineStyle","none",Color="#D95319",HandleVisibility="off")
errorbar(TotEnergy,b2mean,b2min,b2max,"LineStyle","none",Color="#77AC30",HandleVisibility="off")

% For the legend
empty = NaN(1,10);
semilogy(TotEnergy, empty,'diamond', Color="#D95319", DisplayName="Beam 1")
semilogy(TotEnergy, empty,'x', Color="#77AC30", DisplayName="Beam 2")
semilogy(TotEnergy, empty,'o', Color="#0072BD", DisplayName="Both Beams",LineWidth=1.5)

legend

ylim([10^9 2E13])
xlabel("Energy (MeV)")
ylabel("dN/dE/d\Omega (Protons/MeV/sr)")
set(gca,'fontsize',14);
if save
    savefig(gcf, '..\RCF_Images\50JRCF_dNdE_ErrorBars.fig')
    saveas(gcf, '..\RCF_Images\50JRCF_dNdE_ErrorBars.png')
end
hold off

% ##################### Right and Left Ratios #############################
% ########################## Unit Conversion ##############################

leftDose = leftDose*6.2415E9; % Gy/s.r.(J/kg/s.r.) to MeV/g/s.r.
% Dose array is indexed (layer, shot) so doseArr(1,:) is all the H1 layer
% while doseArr(:,1) is all shot "S002"
% Multiply l*p for all H layers (1:5)
leftDose(1:5,:) = leftDose(1:5,:)*10E-4*1.25;
% Same for E layers (6:end)
leftDose(6:end,:) = leftDose(6:end,:)*30E-4*1.1;
% Now we have units of MeV/cm^2/s.r., this is an energy flux per steradian
% Now we will multiply out the area of the RCF
% This is 2in for all but the final 2 sheets
leftDose(:,1:8)=leftDose(:,1:8)*(2*2.54)^2/2;
leftDose(:,9:end)=leftDose(:,9:end)*(4*2.54)^2/2;
% This leaves us with units of MeV/s.r.
% Now I will account for the RCF bin width (energy range over which the
% active layer works -> rowmean
% This is a layer changing variable
for n=1:size(leftDose,1)
    leftDose(n,:)=leftDose(n,:)/rowmean(2*n);
end
% Now we are in units of 1/s.r.
% Now I will address the response function found in Result
% The response functions are for each layer H1, H2, etc.
for n=1:size(leftDose,1)
    rn = 2*n+1; % Indexing for Result data
    for k=1:size(leftDose,2) %Have to iterate through each shot for the proper math
        leftdNdE(n,k)=max(leftDose(n,k)/Result(:,rn));
    end
end
% Dividing by the response funciton of MeV/proton, we get the desired units
% of Proton/MeV/s.r.

rightDose = rightDose*6.2415E9; % Gy/s.r.(J/kg/s.r.) to MeV/g/s.r.
% Dose array is indexed (layer, shot) so doseArr(1,:) is all the H1 layer
% while doseArr(:,1) is all shot "S002"
% Multiply l*p for all H layers (1:5)
rightDose(1:5,:) = rightDose(1:5,:)*10E-4*1.25;
% Same for E layers (6:end)
rightDose(6:end,:) = rightDose(6:end,:)*30E-4*1.1;
% Now we have units of MeV/cm^2/s.r., this is an energy flux per steradian
% Now we will multiply out the area of the RCF
% This is 2in for all but the final 2 sheets
rightDose(:,1:8)=rightDose(:,1:8)*(2*2.54)^2/2;
rightDose(:,9:end)=rightDose(:,9:end)*(4*2.54)^2/2;
% This leaves us with units of MeV/s.r.
% Now I will account for the RCF bin width (energy range over which the
% active layer works -> rowmean
% This is a layer changing variable
for n=1:size(rightDose,1)
    rightDose(n,:)=rightDose(n,:)/rowmean(2*n);
end
% Now we are in units of 1/s.r.
% Now I will address the response function found in Result
% The response functions are for each layer H1, H2, etc.
for n=1:size(rightDose,1)
    rn = 2*n+1; % Indexing for Result data
    for k=1:size(rightDose,2) %Have to iterate through each shot for the proper math
        rightdNdE(n,k)=max(rightDose(n,k)/Result(:,rn));
    end
end
% Dividing by the response funciton of MeV/proton, we get the desired units
% of Proton/MeV/s.r.

% ############################## Plot ratio of dose #######################
figure(4)
plot(TotEnergy, log2(rightdNdE(:,1)./leftdNdE(:,1)),'o', Color="#0072BD", HandleVisibility="off")
hold on
plot(TotEnergy, log2(rightdNdE(:,2)./leftdNdE(:,2)),'o', Color="#0072BD", HandleVisibility="off")
plot(TotEnergy, log2(rightdNdE(:,3)./leftdNdE(:,3)),'o', Color="#0072BD", HandleVisibility="off")
plot(TotEnergy, log2(rightdNdE(:,6)./leftdNdE(:,6)),'o', Color="#0072BD", HandleVisibility="off")
plot(TotEnergy, log2(rightdNdE(:,10)./leftdNdE(:,10)),'o', Color="#0072BD", HandleVisibility="off")

plot(TotEnergy, log2(rightdNdE(:,7)./leftdNdE(:,7)),'diamond', Color="#D95319", HandleVisibility="off")
plot(TotEnergy, log2(rightdNdE(:,9)./leftdNdE(:,9)),'diamond', Color="#D95319", HandleVisibility="off")

plot(TotEnergy, log2(rightdNdE(:,4)./leftdNdE(:,4)),'x', Color="#77AC30", HandleVisibility="off")
plot(TotEnergy, log2(rightdNdE(:,5)./leftdNdE(:,5)),'x', Color="#77AC30", HandleVisibility="off")
plot(TotEnergy, log2(rightdNdE(:,8)./leftdNdE(:,8)),'x', Color="#77AC30", HandleVisibility="off")

% For the legend
empty = NaN(1,10);
plot(TotEnergy, empty,'diamond', Color="#D95319", DisplayName="Beam 1")
plot(TotEnergy, empty,'x', Color="#77AC30", DisplayName="Beam 2")
plot(TotEnergy, empty,'o', Color="#0072BD", DisplayName="Both Beams")

legend(Location="northwest")

ylim([-2 2])
set(gca,'fontsize',14);
xlabel("Energy (MeV)")
ylabel("Log2(Ratio)")
if save
    savefig(gcf, '..\RCF_Images\50JRCFratiodNdE.fig')
    saveas(gcf, '..\RCF_Images\50JRCFratiodNdE.png')
end
hold off


%% ########################################################################
% #################### 50J - Variable Timing ##############################
% #########################################################################
shots = ["S014","S015","S016","S017","S018"]; % Shots of interest

% Loop through all shots and calculate total dose on film
for k = 1:length(shots)
    shot = shots(k);      % Shot Number of interest
    totstacks = length(stacks);
    dosePerSR = NaN(1,10);
    h_index = NaN(1,5);
    SR = 0.769;       % Constant RCF size of 4in

    shotStruct = load(strcat("..\RCF_Images\dataStruct\",shot,".mat"));

    for i=1:HDV2
        h_index(i) = find([shotStruct.dataStruct.Label]==strcat('H',num2str(i)));
    end

    for i=1:HDV2
        layer = cast(shotStruct.dataStruct(h_index(i)).Data(:,:,1), 'double');
        conv2dose = arrayfun(@(x) redFunc(x), layer);
        totalDose = mean(conv2dose,'all','omitnan');
        dosePerSR(i) = totalDose/SR;
    end

    for j = 6:totstacks
        arr = readmatrix(strcat('..\RCF_Images\dataStruct\csv\Dose\',shot,'_',stacks(j),'.csv'));
        totalDose = mean(arr,'all','omitnan');
        dosePerSR(j) = totalDose/SR;
    end

    timeDose(:,k) = dosePerSR;
end

% ######################### Unit Conversion ###############################
timeDose = timeDose*6.2415E9; % Gy/s.r.(J/kg/s.r.) to MeV/g/s.r.
% Dose array is indexed (layer, shot) so doseArr(1,:) is all the H1 layer
% while doseArr(:,1) is all shot "S002"
% Multiply l*p for all H layers (1:5)
timeDose(1:5,:) = timeDose(1:5,:)*10E-4*1.25;
% Same for E layers (6:end)
timeDose(6:end,:) = timeDose(6:end,:)*30E-4*1.1;
% Now we have units of MeV/cm^2/s.r., this is an energy flux per steradian
% Now we will multiply out the area of the RCF
timeDose(:,:)=timeDose(:,:)*(4*2.54)^2;
% This leaves us with units of MeV/s.r.
% Now I will account for the RCF bin width (energy range over which the
% active layer works -> rowmean
% This is a layer changing variable
for n=1:size(timeDose,1)
    timeDose(n,:)=timeDose(n,:)/rowmean(2*n);
end
% Now we are in units of 1/s.r.
% Now I will address the response function found in Result
% The response functions are for each layer H1, H2, etc.
for n=1:size(timeDose,1)
    rn = 2*n+1; % Indexing for Result data
    for k=1:size(timeDose,2) %Have to iterate through each shot for the proper math
        timedNdE(n,k)=max(timeDose(n,k)/Result(:,rn));
    end
end
% Dividing by the response funciton of MeV/proton, we get the desired units
% of Proton/MeV/s.r.

% ############################### Plotting ################################
figure(5)
semilogy(TotEnergy, timedNdE(:,5),'diamond', Color="#77AC30",DisplayName="+4ps") % +4ps
hold on
semilogy(TotEnergy, timedNdE(:,4),'+', Color="#7E2F8E", DisplayName="+2ps") % +2ps
semilogy(TotEnergy, timedNdE(:,1),'o', Color="#0072BD", DisplayName="0ps") % 0ps
semilogy(TotEnergy, timedNdE(:,2),'*', Color="#D95319", DisplayName="-2ps") % -2ps
semilogy(TotEnergy, timedNdE(:,3),'x', Color="#EDB120", DisplayName="-4ps") % -4ps

legend

ylim([10^10 2E13])
xlim([0 25])
set(gca,'fontsize',14);
xlabel("Energy (MeV)")
ylabel("dN/dE/d\Omega (Protons/MeV/sr)")
if save
    savefig(gcf, '..\RCF_Images\TimeRCFdNdE.fig')
    saveas(gcf, '..\RCF_Images\TimeRCFdNdE.png')
end
hold off

%% ########################################################################
% ######################### 25J - 0ps #####################################
% #########################################################################
shots = ["S019","S020","S022"]; % Shots of interest

% Loop through all shots and calculate total dose on film
for k = 1:length(shots)
    shot = shots(k);      % Shot Number of interest
    totstacks = length(stacks);
    dosePerSR = NaN(1,10);
    h_index = NaN(1,5);
    SR = 0.769;

    shotStruct = load(strcat("..\RCF_Images\dataStruct\",shot,".mat"));

    for i=1:HDV2
        h_index(i) = find([shotStruct.dataStruct.Label]==strcat('H',num2str(i)));
    end

    for i=1:HDV2
        layer = cast(shotStruct.dataStruct(h_index(i)).Data(:,:,1), 'double');
        conv2dose = arrayfun(@(x) redFunc(x), layer);
        totalDose = mean(conv2dose,'all','omitnan');
        dosePerSR(i) = totalDose/SR;
    end

    for j = 6:totstacks
        arr = readmatrix(strcat('..\RCF_Images\dataStruct\csv\Dose\',shot,'_',stacks(j),'.csv'));
        totalDose = mean(arr,'all','omitnan');
        dosePerSR(j) = totalDose/SR;
    end

    dose25(:,k) = dosePerSR;
end

% ############################# Unit Conversion ###########################
dose25 = dose25*6.2415E9; % Gy/s.r.(J/kg/s.r.) to MeV/g/s.r.
% Dose array is indexed (layer, shot) so doseArr(1,:) is all the H1 layer
% while doseArr(:,1) is all shot "S002"
% Multiply l*p for all H layers (1:5)
dose25(1:5,:) = dose25(1:5,:)*10E-4*1.25;
% Same for E layers (6:end)
dose25(6:end,:) = dose25(6:end,:)*30E-4*1.1;
% Now we have units of MeV/cm^2/s.r., this is an energy flux per steradian
% Now we will multiply out the area of the RCF
dose25(:,:)=dose25(:,:)*(4*2.54)^2;
% This leaves us with units of MeV/s.r.
% Now I will account for the RCF bin width (energy range over which the
% active layer works -> rowmean
% This is a layer changing variable
for n=1:size(dose25,1)
    dose25(n,:)=dose25(n,:)/rowmean(2*n);
end
% Now we are in units of 1/s.r.
% Now I will address the response function found in Result
% The response functions are for each layer H1, H2, etc.
for n=1:size(dose25,1)
    rn = 2*n+1; % Indexing for Result data
    for k=1:size(dose25,2) %Have to iterate through each shot for the proper math
        dNdE25(n,k)=max(dose25(n,k)/Result(:,rn));
    end
end
% Dividing by the response funciton of MeV/proton, we get the desired units
% of Proton/MeV/s.r.

% ############################# Plotting ##################################
figure(6)
semilogy(TotEnergy,bbmean,'o', Color="#0072BD",HandleVisibility="off")
hold on
errorbar(TotEnergy,bbmean,bbmin,bbmax,"LineStyle","none",Color="#0072BD",HandleVisibility="off")
semilogy(TotEnergy, dNdE25(:,1),'diamond', Color="#D95319", HandleVisibility="off",LineWidth=1.5)
semilogy(TotEnergy, dNdE25(:,3),'diamond', Color="#D95319", HandleVisibility="off",LineWidth=1.5)
semilogy(TotEnergy, dNdE25(:,2),'x',Color="#77AC30",HandleVisibility="off",LineWidth=1.5)

% For the legend
empty = NaN(1,10);
semilogy(TotEnergy, empty,'o', Color="#0072BD", DisplayName="Both Beams - 50J")
semilogy(TotEnergy, empty,'diamond', Color="#D95319", DisplayName="Beam 1 - 25J",LineWidth=1.5)
semilogy(TotEnergy, empty,'x',Color="#77AC30",DisplayName="Beam 2 - 25J",LineWidth=1.5)

legend

ylim([10^9 2E13])
xlim([0 25])
set(gca,'fontsize',14);
xlabel("Energy (MeV)")
ylabel("dN/dE/d\Omega (Protons/MeV/sr)")
if save
    savefig(gcf, '..\RCF_Images\50JRCF_25_dnde.fig')
    saveas(gcf, '..\RCF_Images\50JRCF_25_dnde.png')
end
hold off