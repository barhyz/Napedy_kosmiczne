clear; close all;

fuels = ["HTPB"; "GAP"; "PE"; "ABS"];
oxidizers = ["H2O2"; "N2O"; "LOX"];
% fuel_formula = ['C',10,'H',15.4,'O',0.07;
%                 'C',3,'H',5,'N',3,'O',1;
%                 'C',2,'H',4];
% oxid_formula = ['H',2,'O',2;
%                 'N',2,'O',1;
%                 'O',2];
fuel_enthalpy = [-51.88;
    142.26;
    -53.81;
    62.63];
oxid_enthalpy = [-187.4;
    82.05;
    -12.97];
fuel_no=4;
oxid_no=1;

o_f = linspace(1,10,51);
press = 30;%linspace(5,60,51);
%viscosity = zeros(51,51);

x=CEA('problem','hp','equilibrium','o/f',o_f,'case','CEAM-HP1','p,bar',press, ...
    'reactants','fuel','ABS','C',3.85,'H',4.85,'N',0.43,'wt%',100,'h,kJ/mol',fuel_enthalpy(fuel_no),'t(k)',300, ...
    'oxid','H2O2','H',2,'O',2,'wt%',100,'h,kJ/mol',oxid_enthalpy(oxid_no),'t(k)',300,...
    'output','transport','end');

% for i=1:51
%     curr_press = press(i);
%     x=CEA('problem','hp','equilibrium','o/f',o_f,'case','CEAM-HP1','p,bar',curr_press, ...
%         'reactants','fuel','Amoniak','C',3,'H',8,'wt%',100,'h,kJ/mol',enthalpy_formation(fuel_no),'t(k)',300, ...
%         'oxid','N2O','N',2,'O',1,'wt%',100,'h,kJ/mol',82.05,'t(k)',300,...
%         'output','transport','end');
%     viscosity(:,i)=x.output.viscosity;
% end
% 
% x=CEA('problem','hp','equilibrium','o/f',o_f,'case','CEAM-HP1','p,bar',press, ...
%     'reactants','fuel','Amoniak','C',3,'H',8,'wt%',100,'h,kJ/mol',enthalpy_formation(fuel_no),'t(k)',300, ...
%     'oxid','N2O','N',2,'O',1,'wt%',100,'h,kJ/mol',82.05,'t(k)',300,...
%     'output','transport','end');
% 
% x.output.viscosity = viscosity;

temp=x.output.temperature;
kappa=x.output.gamma;
gas_constant=8314./x.output.mw;
c_star=zeros(size(o_f,2),size(press,2));

for i=1:size(press,2)
    for j=1:size(o_f,2)
        c_star(j,i) = sqrt(kappa(j,i) * gas_constant(j,i) * temp(j,i))/...
            (kappa(j,i) * sqrt( (2 / (kappa(j,i)+1) ) ^ ((kappa(j,i)+1) / (kappa(j,i)-1) )));
    end
end

f1=figure;
plot(o_f,c_star(:,1));
hold on
for i=2:size(press,2)
    plot(o_f,c_star(:,i));
end
grid minor;
xlabel('O/F[-]');
ylabel('c*[m/s]');
%legend('p=20[bar]','p=35[bar]','p=50[bar]','p=65[bar]','p=80[bar]');
saveas(f1,strcat('output\',oxidizers(oxid_no),'_',fuels(fuel_no),'_c_star.png'));
saveas(f1,strcat('output\',oxidizers(oxid_no),'_',fuels(fuel_no),'_c_star.fig'));

f2=figure;
plot(o_f,temp(:,1));
hold on
for i=2:size(press,2)
    plot(o_f,temp(:,i));
end
grid minor;
xlabel('O/F[-]');
ylabel('T[K]');
%legend('p=20[bar]','p=35[bar]','p=50[bar]','p=65[bar]','p=80[bar]');
saveas(f2,strcat('output\',oxidizers(oxid_no),'_',fuels(fuel_no),'_temp.png'));
saveas(f2,strcat('output\',oxidizers(oxid_no),'_',fuels(fuel_no),'_temp.fig'));
save(strcat('output\',oxidizers(oxid_no),'_',fuels(fuel_no),'.mat'),'x');

%% c star plot
% Get list of subfolders in the current directory
% folders = dir;
% subfolders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));
% 
% % Initialize a structure to store data
% allData = struct;
% fileIndex = 1;
% 
% % Loop through each subfolder
% for k = 1:length(subfolders)
%     folderPath = fullfile(pwd, subfolders(k).name);  % Full path to subfolder
%     matFiles = dir(fullfile(folderPath, '*.mat'));   % All .mat files in subfolder
% 
%     for i = 1:length(matFiles)
%         matFilePath = fullfile(folderPath, matFiles(i).name);
%         fprintf('Loading %s\n', matFilePath);  % Display which file is being loaded
% 
%         data = load(matFilePath);
% 
%         % Option 1: Store data in a struct array
%         allData(fileIndex).folder = subfolders(k).name;
%         allData(fileIndex).filename = matFiles(i).name;
%         allData(fileIndex).data = data;
% 
%         fileIndex = fileIndex + 1;
%     end
% end

x_c_star = zeros(51,9);
names=["H2O2_HTPB";"H2O2_GAP";"H2O2_PE";"H2O2_ABS";"N2O_HTPB";"N2O_GAP";"N2O_PE";"N2O_ABS";"LOX_HTPB";"LOX_GAP";"LOX_PE";"LOX_ABS"];
figure;

for i = 1:size(names,1)
    load(strcat('output\',names(i),".mat"));
    temp=x.output.temperature;
    kappa=x.output.gamma;
    gas_constant=8314./x.output.mw;
    c_star = sqrt(kappa .* gas_constant .* temp)./...
        (kappa .* sqrt( (2 ./ (kappa+1) ) .^ ((kappa+1) ./ (kappa-1) )));
    x_c_star(:,i)=c_star;
    plot(o_f,c_star);
    hold on;
end
grid minor;
names_label=["H2O2/HTPB";"H2O2/GAP";"H2O2/PE";"H2O2/ABS";"N2O/HTPB";"N2O/GAP";"N2O/PE";"N2O/ABS";"LOX/HTPB";"LOX/GAP";"LOX/PE";"LOX/ABS"];
legend(names_label);
xlabel("O/F [-]");
ylabel("c* [m/s]");

%% Enthalpy of formation
% fuel_enthalpy = 62.63;
% oxid_enthalpy = 82.05;
% o_f = linspace(3,10,51);
% press = 30;%linspace(5,60,51);
% 
% x=CEA('problem','hp','equilibrium','o/f',o_f,'case','CEAM-HP1','p,bar',press, ...
%     'reactants','fuel','ABS','C',3.85,'H',4.85,'N',0.85,'wt%',100,'h,kJ/mol',fuel_enthalpy,'t(k)',300, ...
%     'oxid','N2O','N',2,'O',1,'wt%',100,'h,kJ/mol',oxid_enthalpy,'t(k)',300,...
%     'output','transport');
% 
% f2=figure;
% plot(o_f,x.output.temperature(:));
% hold on
% % for i=2:size(press,2)
% %     plot(o_f,temp(:,i));
% % end
% grid minor;
% xlabel('O/F[-]');
% ylabel('T[K]');





