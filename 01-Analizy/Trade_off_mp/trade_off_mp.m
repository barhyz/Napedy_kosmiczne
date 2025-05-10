clear; close all;

fuels=["Metanol"; "Etanol"; "Izopropanol"; "Amoniak"; "Propan"];
% oxidizer = "HTP";
% 
% fuel_formula = ['C',1,'H',4,'O',1;
%                 'C',2,'H',6,'O',1;
%                 'C',3,'H',8,'O',1;
%                 'H',3,'N',1;
%                 'C',3,'H',8];
% oxid_formula = ['H',2,'O',2;
%                 'N',2,'O',1];
enthalpy_formation= [-238.4;
                     -277.38;
                     -317;
                     -46;
                     -103.5];
% oxid_enthalpy = -238.4;
fuel_no=5;

o_f = linspace(1,12,51);
press = linspace(5,60,51);

x=CEA('problem','hp','equilibrium','o/f',o_f,'case','CEAM-HP1','p,bar',press, ...
    'reactants','fuel','Amoniak','C',3,'H',8,'wt%',100,'h,kJ/mol',enthalpy_formation(fuel_no),'t(k)',300, ...
    'oxid','N2O','N',2,'O',1,'wt%',100,'h,kJ/mol',82.05,'t(k)',300,...
    'output','transport','end');

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
legend('p=20[bar]','p=35[bar]','p=50[bar]','p=65[bar]','p=80[bar]');
saveas(f1,strcat('output_n2o\',fuels(fuel_no),'_c_star.png'));
saveas(f1,strcat('output_n2o\',fuels(fuel_no),'_c_star.fig'));

f2=figure;
plot(o_f,temp(:,1));
hold on
for i=2:size(press,2)
    plot(o_f,temp(:,i));
end
grid minor;
xlabel('O/F[-]');
ylabel('T[K]');
legend('p=20[bar]','p=35[bar]','p=50[bar]','p=65[bar]','p=80[bar]');
saveas(f2,strcat('output_n2o\',fuels(fuel_no),'_temp.png'));
saveas(f2,strcat('output_n2o\',fuels(fuel_no),'_temp.fig'));
save(strcat('output_n2o\',fuels(fuel_no),'.mat'),'x');
