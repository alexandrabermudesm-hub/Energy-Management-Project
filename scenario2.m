% Load data
%data = readtable('HOURLYDATA1.xlsx', 'Sheet', 'ProduktionsStatistik');
[num, txt, raw] = xlsread('HOURLYDATA1.xlsx');
DH_demand = raw(2:end, 22);
DH_demand = cellfun(@(x) double(x), DH_demand);
DH_production = raw(2:end, 24);
DH_production = cellfun(@(x) double(x), DH_production);
%raw(8748:8766, 22)
%data = readtable('HOURLYDATA1.xlsx');
%DH_production = data.Var24;
%DH_production = str2double(DH_prodstr);
%DH_demand = data.Var22;
%DH_demand = str2double(DH_demandstr);
%data.Properties.VariableNames
%DH_production = data.Var24;          % Nuclear district heating production
%DH_demand = data.Var22; % District heating demand

% We need data for just 8766 hours (1 year)
DH_production = DH_production(1:8766);
DH_demand = DH_demand(1:8766);

% TES SIZE
TES_capacity = 7225.67 * 24 * 1;   
TES_loss_rate = 0.005;         
TES_state = 0;                 
%TES_state = 0.5*TES_capacity;

energy_nuclear = 0;    
energy_TES_used = 0;  
energy_boiler = 0;     

% OPERATING MODES
hours_only_nuclear = 0;        
hours_nuclear_TES_charge = 0;   
hours_nuclear_TES_discharge = 0; 
hours_nuclear_boiler = 0;        

TES_SOC_series = zeros(8766, 1);
el_boiler = zeros(8766,1);
mode = zeros(8766, 1);  

for t = 1:8766
    % Update TES state with losses
    TES_state = TES_state * (1 - TES_loss_rate);
    
    nuclear_supply = DH_production(t);
    demand = DH_demand(t);
    
    if nuclear_supply >= demand
        % Nuclear fully covers demand
        energy_nuclear = energy_nuclear + demand;
        surplus = nuclear_supply - demand;
        
        if TES_state < TES_capacity
            % Charge TES with surplus 
            charge = min(surplus, TES_capacity - TES_state);
            TES_state = TES_state + charge;
            hours_nuclear_TES_charge = hours_nuclear_TES_charge + 1;
            mode(t) = 2; % Nuclear + TES charging
        else
            % TES is already full
            hours_only_nuclear = hours_only_nuclear + 1;
            mode(t) = 1; % Only nuclear
        end
        
    else
        % Nuclear supply is less than demand
        shortage = demand - nuclear_supply;
        energy_nuclear = energy_nuclear + nuclear_supply;
        
        if TES_state >= shortage
            % TES can cover the shortage completely
            TES_state = TES_state - shortage;
            energy_TES_used = energy_TES_used + shortage;
            hours_nuclear_TES_discharge = hours_nuclear_TES_discharge + 1;
            mode(t) = 3; % Nuclear + TES discharging
        else
            % TES cannot cover the full shortage: use full TES and cover remaining with boiler
            discharged = TES_state;
            TES_state = 0;
            energy_TES_used = energy_TES_used + discharged;
            remaining_shortage = shortage - discharged;
            energy_boiler = energy_boiler + remaining_shortage;
            hours_nuclear_boiler = hours_nuclear_boiler + 1;
            mode(t) = 4; % Nuclear + boiler
            el_boiler(t) = remaining_shortage;
        end
    end
    TES_SOC_series(t) = TES_state;
end

energy_nuclear = energy_nuclear + energy_TES_used;



% PIE CHART
figure;
labels = {'Only Nuclear DH Supply', 'Nuclear + TES Charging', 'Nuclear + TES Discharging', 'Nuclear + Boiler'};
sizes = [hours_only_nuclear, hours_nuclear_TES_charge, hours_nuclear_TES_discharge, hours_nuclear_boiler];
pie(sizes, labels);
title('Operating Modes of District Heating System with TES');

% TES SOC
figure;
plot(1:8766, TES_SOC_series);
title('TES State of Charge Over Time');
xlabel('Hour'); 
ylabel('TES Energy [kWh]');
grid on;



% 15-21 January: 
indicesJan = 360:528;
figure;
plot(1:length(indicesJan), DH_demand(indicesJan), 1:length(indicesJan), DH_production(indicesJan), 1:length(indicesJan), TES_SOC_series(indicesJan),1:length(indicesJan),el_boiler(indicesJan));
legend('DH Demand', 'DH Production', 'TES SOC', 'Boiler production');
title('15-21 January');
xlabel('Hour'); ylabel('kW/kWh'); grid on;

% 1-7 April: 
indicesApr = 2160:2328;
figure;
plot(1:length(indicesApr), DH_demand(indicesApr), 1:length(indicesApr), DH_production(indicesApr), 1:length(indicesApr), TES_SOC_series(indicesApr), 1:length(indicesApr),el_boiler(indicesApr));
legend('DH Demand', 'DH Production', 'TES SOC', 'Boiler production');
title('1-7 April');
xlabel('Hour'); ylabel('kW/kWh'); grid on;

% 8-14 July: 
indicesJul = 4512:4680;
figure;
plot(1:length(indicesJul), DH_demand(indicesJul), 1:length(indicesJul), DH_production(indicesJul), 1:length(indicesJul), TES_SOC_series(indicesJul), 1:length(indicesJul),el_boiler(indicesJul));
legend('DH Demand', 'DH Production', 'TES SOC', 'Boiler production');
title('8-14 July');
xlabel('Hour'); ylabel('kW/kWh'); grid on;

% 23-29 October:
indicesOct = 7104:7272;
figure;
plot(1:length(indicesOct), DH_demand(indicesOct), 1:length(indicesOct), DH_production(indicesOct), 1:length(indicesOct), TES_SOC_series(indicesOct), 1:length(indicesOct),el_boiler(indicesOct));
legend('DH Demand', 'DH Production', 'TES SOC', 'Boiler production');
title('23-29 October');
xlabel('Hour');ylabel('kW/kWh');grid on;
% CALCULATIONS

grid_emission_factor = 0.041;  % kg CO2 per kWh
total_emissions = energy_boiler * grid_emission_factor;  % Total emissions in kg CO2


total_DH_covered = energy_nuclear + energy_boiler;  % in kWh
co2_intensity = total_emissions / (total_DH_covered / 1000);  % kg CO2 per MWh


fprintf('\n--- SCENARIO: Nuclear + TES + Boiler ---\n');
fprintf('Total DH Demand Covered: %.2f MWh\n', total_DH_covered / 1000);
fprintf('Covered by Nuclear: %.2f MWh\n', energy_nuclear / 1000);
fprintf('Covered by Boiler: %.2f MWh\n', energy_boiler / 1000);
fprintf('Total CO₂ Emissions: %.2f tons\n', total_emissions / 1000);
fprintf('CO₂ Intensity: %.2f kg CO₂/MWh\n', co2_intensity);

nanIndices = find(isnan(DH_production));