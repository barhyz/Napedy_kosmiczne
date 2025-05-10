function [Tank] = TankModel(Tank, Tank_Oxidizer, downstream_pressure, time, initial_time, Feed_Input, n, warning_flag, Simulation_Settings)

%Stała czasowa
dt = time - initial_time;  %time step

%Sprawdzenie przypadku obliczeniowego - dwufazowy lub jednofazowy


if (Tank.two_phase == true)
    Tank.pressure        = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,2),Tank.temperature,'spline','extrap');
    Tank.density_liquid  = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,3),Tank.temperature,'spline','extrap');
    Tank.density_vapour  = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,4),Tank.temperature,'spline','extrap');
    Tank.enthalpy_liquid = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,5),Tank.temperature,'spline','extrap');
    Tank.enthalpy_vapour = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,6),Tank.temperature,'spline','extrap');
    Tank.energy_liquid   = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,7),Tank.temperature,'spline','extrap');
    Tank.energy_vapour   = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,8),Tank.temperature,'spline','extrap');
    Tank.enthropy_liquid = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,9),Tank.temperature,'spline','extrap');
    Tank.enthropy_vapour = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,10),Tank.temperature,'spline','extrap');
          
    downstream_density_liquid  = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,3),downstream_pressure,'spline','extrap');
    downstream_density_vapour  = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,4),downstream_pressure,'spline','extrap');
    downstream_enthalpy_liquid = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,5),downstream_pressure,'spline','extrap');
    downstream_enthalpy_vapour = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,6),downstream_pressure,'spline','extrap');
    %downstream_energy_liquid = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,7),downstream_pressure,'spline','extrap');
    %downstream_energy_vapour = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,8),downstream_pressure,'spline','extrap');
    downstream_enthropy_liquid = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,9),downstream_pressure,'spline','extrap');
    downstream_enthropy_vapour = interp1(Tank_Oxidizer.properties_interp(:,2),Tank_Oxidizer.properties_interp(:,10),downstream_pressure,'spline','extrap');

    %Tank.quality_fluid = Tank.mass_vapour/(Tank.mass_vapour+Tank.mass_liquid);
    Tank.quality_fluid =0;
    
    Tank.density_fluid  = Tank.quality_fluid * Tank.density_vapour + (1-Tank.quality_fluid) * Tank.density_liquid;
    Tank.enthalpy_fluid = Tank.quality_fluid * Tank.enthalpy_vapour + (1-Tank.quality_fluid) * Tank.enthalpy_liquid;
    Tank.enthropy_fluid = Tank.quality_fluid * Tank.enthropy_vapour + (1-Tank.quality_fluid) * Tank.enthropy_liquid;
    
    downstream_quality_fluid = Tank.quality_fluid +...
    (Tank.enthropy_liquid - downstream_enthropy_liquid)/(downstream_enthropy_vapour-downstream_enthropy_liquid);
    
    downstream_volume_vapour = 1/downstream_density_vapour;
    downstream_volume_liquid = 1/downstream_density_liquid;
    downstream_volume_fluid  = downstream_quality_fluid * downstream_volume_vapour + (1-downstream_quality_fluid) * downstream_volume_liquid;
    downstream_density_fluid  = 1/downstream_volume_fluid;
    downstream_enthalpy_fluid = downstream_quality_fluid * downstream_enthalpy_vapour + (1-downstream_quality_fluid) * downstream_enthalpy_liquid;
   
    Tank.oxidizer_mass_flow_rate_HEM = Feed_Input.discharge_coefficient*Feed_Input.orifices_number*...
    (Feed_Input.orifices_diameter/2)^2*pi*downstream_density_fluid*...
    sqrt(2*(Tank.enthalpy_fluid-downstream_enthalpy_fluid));

    Tank.oxidizer_mass_flow_rate_SPI =  Feed_Input.discharge_coefficient*Feed_Input.orifices_number*...
    (Feed_Input.orifices_diameter/2)^2*pi*...
    sqrt(2*Tank.density_liquid*(Tank.pressure-downstream_pressure));
    
    Tank.oxidizer_mass_flow_rate = (Tank.oxidizer_mass_flow_rate_SPI+...
    Tank.oxidizer_mass_flow_rate_HEM)/2; 

    Tank.oxidizer_total_mass = Tank.oxidizer_total_mass - Tank.oxidizer_mass_flow_rate*dt;
    
    temperature_max=Tank.temperature;
    temperature_min=Tank.temperature - 1;
    Tank.temperature = (temperature_max+temperature_min)/2;
    
    while (true)
            
            new_density_liquid  = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,3),Tank.temperature,'spline','extrap');
            new_density_vapour  = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,4),Tank.temperature,'spline','extrap');
            new_energy_liquid   = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,7),Tank.temperature,'spline','extrap');
            new_energy_vapour   = interp1(Tank_Oxidizer.properties_interp(:,1),Tank_Oxidizer.properties_interp(:,8),Tank.temperature,'spline','extrap');
            
            new_mass_liquid =(Tank_Oxidizer.volume-(Tank.oxidizer_total_mass/new_density_vapour))/...
            ((1/new_density_liquid)-(1/new_density_vapour));
            new_mass_vapour = Tank.oxidizer_total_mass - new_mass_liquid;
            
            RHS = Tank.mass_liquid * (new_energy_liquid - Tank.energy_liquid) + ...
            Tank.energy_liquid * (new_mass_liquid - Tank.mass_liquid) + ...
            Tank.mass_vapour * (new_energy_vapour - Tank.energy_vapour) + ...
            Tank.energy_vapour * (new_mass_vapour - Tank.mass_vapour) + ...
            Tank.oxidizer_mass_flow_rate * Tank.enthalpy_liquid * dt;
            
            if (RHS > 10e-4)
                temperature_max = Tank.temperature;                
                Tank.temperature = (temperature_max + temperature_min)/2;
            elseif (RHS < -10e-4)
                temperature_min = Tank.temperature;
                Tank.temperature = (temperature_max + temperature_min)/2;
            else
                break;                  
            end
            %disp(RHS);
                           
    end
    Tank.mass_liquid = new_mass_liquid;
    Tank.mass_vapour = new_mass_vapour;
    
    if (Tank.mass_liquid < 0.001)
        Tank.two_phase = false;
        Tank.single_phase = true;
        Tank.first_vapour_iteration=true;
        Tank.mass_vapour = Tank.mass_vapour + Tank.mass_liquid;
    end
end





%Stan gazowy
if (Tank.single_phase == true)
    
    
        
    if (Tank.first_vapour_iteration==true)
        

        Tank.initial_density_vapour = Tank.density_vapour;      
        Tank.initial_vapour_mass = Tank.mass_vapour;            
        Tank.initial_vapour_pressure = Tank.pressure;      
        Tank.initial_vapour_temperature = Tank.temperature;     
        
        Tank.initial_compressibility_factor = LinearInterpolation(Tank.initial_vapour_pressure, 0, Tank_Oxidizer.critical_pressure,...
        1, Tank_Oxidizer.critical_compressibility_factor);
        Tank.first_vapour_iteration=false;
    end
    
    
   

    Tank.oxidizer_mass_flow_rate =  Feed_Input.discharge_coefficient*Feed_Input.orifices_number*...
    (Feed_Input.orifices_diameter/2)^2*pi*...
    sqrt(2*Tank.density_vapour*(Tank.pressure-downstream_pressure));
    
    if (isnan(Tank.oxidizer_mass_flow_rate))
    b=0;
    end

 
    Tank.mass_vapour = Tank.mass_vapour - Tank.oxidizer_mass_flow_rate * dt;
    
    Tank.oxidizer_total_mass = Tank.oxidizer_total_mass - Tank.oxidizer_mass_flow_rate * dt;
    
    
    
    Tank.compressibility_factor = LinearInterpolation(Tank.pressure, 0, Tank_Oxidizer.critical_pressure, 1, Tank_Oxidizer.critical_compressibility_factor);
    
    step = 1.0/0.9;
    Aim=0;
    
  
    while (true)
        
    
        Tank.temperature = Tank.initial_vapour_temperature * ((Tank.mass_vapour * Tank.compressibility_factor)/...
        ( Tank.initial_vapour_mass* Tank.initial_compressibility_factor))^(Tank_Oxidizer.gamma - 1.0);
        
      
        Tank.pressure = Tank.initial_vapour_pressure * (Tank.temperature/Tank.initial_vapour_temperature)^...
        (Tank_Oxidizer.gamma/(Tank_Oxidizer.gamma-1.0));
        

        Tank.density_vapour = Tank.initial_density_vapour * (Tank.temperature/Tank.initial_vapour_temperature)^...
        (1/(Tank_Oxidizer.gamma-1.));
        
        Tank.current_compressibility_factor = LinearInterpolation(Tank.pressure, 0, Tank_Oxidizer.critical_pressure, 1, Tank_Oxidizer.critical_compressibility_factor);
        OldAim=Aim;
        
        if (Tank.compressibility_factor <  Tank.current_compressibility_factor)
            Tank.compressibility_factor =  Tank.compressibility_factor * step; %increasing guess
            Aim = 1;
            
        else
            Tank.compressibility_factor = Tank.compressibility_factor / step; %decreasing guess
            Aim = -1;
        end
        
        if (Aim == -OldAim)
            step = sqrt(step);
        end
        
        if( ( abs( Tank.compressibility_factor - Tank.current_compressibility_factor) <0.001) ) %precision of found Z      
            break;
        end
        
        if any(~isreal([Tank.temperature, Tank.pressure, Tank.oxidizer_total_mass, Tank.mass_vapour, Tank.mass_liquid])) && ...
            Simulation_Settings.Booleans.Cd_optimization == false;
            warning('Detected inf/nan/complex numbers in TankModel');
            warning_flag=true;
        end
    end
end



Tank.Ixx=Tank.oxidizer_total_mass*(Tank_Oxidizer.diameter/2)^2*0.5;
Tank.Iyy=(Tank.oxidizer_total_mass*(3*(Tank_Oxidizer.diameter/2)^2+Tank_Oxidizer.length^2))/12;

if any(~isreal([Tank.temperature, Tank.pressure, Tank.oxidizer_total_mass, Tank.mass_vapour, Tank.mass_liquid])) && ...
            Simulation_Settings.Booleans.Cd_optimization == false;
    warning('Detected inf/nan/complex numbers in TankModel');
    warning_flag=true;
end

end


