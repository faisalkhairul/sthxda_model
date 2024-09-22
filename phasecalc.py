import math
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import requests
import json
import numpy as np

# Database properti fluida
fluid_database = {
    "water": {"name": "Water", "coolprop_name": "Water"},
    "air": {"name": "Air", "coolprop_name": "Air"},
    "steam": {"name": "Steam", "coolprop_name": "Water"},
    "r134a": {"name": "R134a", "coolprop_name": "R134a"},
    "ethanol": {"name": "Ethanol", "coolprop_name": "Ethanol"},
    "methanol": {"name": "Methanol", "coolprop_name": "Methanol"},
    "propane": {"name": "Propane", "coolprop_name": "Propane"}
}

def get_fluid_properties(fluid, T, P):
    try:
        coolprop_name = fluid_database[fluid.lower()]["coolprop_name"]
        props = {
            'density': PropsSI('D', 'T', T + 273.15, 'P', P, coolprop_name),
            'viscosity': PropsSI('V', 'T', T + 273.15, 'P', P, coolprop_name),
            'conductivity': PropsSI('L', 'T', T + 273.15, 'P', P, coolprop_name),
            'specific_heat': PropsSI('C', 'T', T + 273.15, 'P', P, coolprop_name),
            'saturation_temp': PropsSI('T', 'P', P, 'Q', 0, coolprop_name) - 273.15
        }
        return props
    except:
        return get_nist_data(fluid, T, P)

def get_nist_data(fluid, temperature, pressure):
    base_url = "https://webbook.nist.gov/cgi/fluid.cgi"
    params = {
        "ID": fluid,
        "Action": "Data",
        "T": temperature + 273.15,
        "P": pressure,
        "Format": "JSON"
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        data = json.loads(response.text)
        return {
            'density': data.get('density', 0),
            'viscosity': data.get('viscosity', 0),
            'conductivity': data.get('thermal conductivity', 0),
            'specific_heat': data.get('isobaric heat capacity', 0),
            'saturation_temp': data.get('temperature', 0) - 273.15
        }
    else:
        raise ValueError(f"Unable to fetch data for {fluid} at T={temperature}°C and P={pressure}Pa")

def determine_phase(fluid, T, P):
    try:
        coolprop_name = fluid_database[fluid.lower()]["coolprop_name"]
        density = PropsSI('D', 'T', T + 273.15, 'P', P, coolprop_name)
        density_liquid = PropsSI('D', 'T', T + 273.15, 'Q', 0, coolprop_name)
        density_vapor = PropsSI('D', 'T', T + 273.15, 'Q', 1, coolprop_name)
        
        if density > density_liquid * 0.95:
            return "Liquid"
        elif density < density_vapor * 1.05:
            return "Gas"
        else:
            return "Two-phase"
    except:
        return "Undefined (supercritical or out of range)"

def calculate_lmtd(T_hot_in, T_hot_out, T_cold_in, T_cold_out, flow_type):
    if flow_type == "counter":
        delta_T1 = T_hot_in - T_cold_out
        delta_T2 = T_hot_out - T_cold_in
    else:  # parallel flow
        delta_T1 = T_hot_in - T_cold_in
        delta_T2 = T_hot_out - T_cold_out
    
    if abs(delta_T1 - delta_T2) < 0.01:  # Avoid division by zero
        return delta_T1
    return (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

def calculate_heat_transfer_coefficient(Re, Pr, k, D, L):
    if Re > 10000 and 0.7 <= Pr <= 160:
        Nu = 0.023 * Re**0.8 * Pr**0.4
    else:
        Nu = 3.66 + (0.0668 * (D / L) * Re * Pr) / (1 + 0.04 * ((D / L) * Re * Pr)**(2/3))
    return Nu * k / D

def calculate_pressure_drop(m, rho, mu, D, L, N_passes, N_tubes, D_shell, N_baffles=0, is_shell_side=False):
    if is_shell_side:
        A = math.pi / 4 * (D_shell**2 - N_tubes * D**2)
    else:
        A = math.pi / 4 * D**2 * N_tubes / N_passes
    
    v = m / (rho * A)
    Re = rho * v * D / mu
    
    if Re < 2300:
        f = 64 / Re
    else:
        f = 0.316 * Re**-0.25
    
    dP_straight = f * (L / D) * 0.5 * rho * v**2
    
    if is_shell_side and N_baffles > 0:
        dP_baffles = N_baffles * 0.5 * rho * v**2
    else:
        dP_baffles = 0
    
    if not is_shell_side:
        dP_turns = 4 * (N_passes - 1) * 0.5 * rho * v**2
    else:
        dP_turns = 0
    
    return (dP_straight + dP_baffles + dP_turns) * N_passes

def heat_exchanger_simulation(mode, **params):
    if mode == "design":
        Q_required = params['m_hot'] * params['cp_hot'] * (params['T_hot_in'] - params['T_hot_out'])
        LMTD = calculate_lmtd(params['T_hot_in'], params['T_hot_out'], params['T_cold_in'], params['T_cold_out'], params['flow_type'])
        U_assumed = 500  # Assumed overall heat transfer coefficient, W/m²K
        A_required = Q_required / (U_assumed * LMTD)
        
        print(f"Required heat transfer area: {A_required:.2f} m²")
        print(f"Required heat transfer: {Q_required:.2f} W")
        
    elif mode == "analysis":
        L_eff_tube = params['L'] * params['N_tube_passes']
        L_eff_shell = params['L'] * params['N_shell_passes']

        h_hot = calculate_heat_transfer_coefficient(params['Re_hot'], params['Pr_hot'], params['k_hot'], params['D_tube'], L_eff_tube)
        h_cold = calculate_heat_transfer_coefficient(params['Re_cold'], params['Pr_cold'], params['k_cold'], params['D_shell'] - params['D_tube'], L_eff_shell)
        U = 1 / (1/h_hot + 1/h_cold)

        A = math.pi * params['D_tube'] * L_eff_tube * params['N_tubes']
        C_hot = params['m_hot'] * params['cp_hot']
        C_cold = params['m_cold'] * params['cp_cold']
        C_min = min(C_hot, C_cold)
        C_max = max(C_hot, C_cold)
        
        if C_min <= 0:
            raise ValueError("C_min must be positive")
        
        NTU = U * A / C_min
        Cr = C_min / C_max if C_max != 0 else 0

        if params['flow_type'] == "counter" and params['N_shell_passes'] == 1 and params['N_tube_passes'] == 2:
            try:
                P = (1 - Cr + (1 + Cr**2)**0.5) / 2
                effectiveness = 2 / (1 + Cr + (1 + Cr**2)**0.5 * (1 + math.exp(-NTU * (1 + Cr**2)**0.5)) / (1 - math.exp(-NTU * (1 + Cr**2)**0.5)))
            except (ValueError, ZeroDivisionError, OverflowError):
                effectiveness = 1 - math.exp(NTU**0.22 * (math.exp(-Cr * NTU**0.78) - 1) / Cr)
        elif params['flow_type'] == "counter":
            if Cr < 1:
                try:
                    effectiveness = (1 - math.exp(-NTU * (1 - Cr))) / (1 - Cr * math.exp(-NTU * (1 - Cr)))
                except (ValueError, OverflowError):
                    effectiveness = 1 - Cr if NTU > 20 else NTU / (1 + NTU)
            else:
                effectiveness = NTU / (1 + NTU)
        else:  # parallel flow
            try:
                effectiveness = (1 - math.exp(-NTU * (1 + Cr))) / (1 + Cr)
            except (ValueError, OverflowError):
                effectiveness = 1 / (1 + Cr) if NTU > 20 else (1 - math.exp(-NTU)) / (1 + Cr)
        
        Q_max = C_min * (params['T_hot_in'] - params['T_cold_in'])
        Q_actual = effectiveness * Q_max
        
        T_hot_out = params['T_hot_in'] - Q_actual / C_hot
        T_cold_out = params['T_cold_in'] + Q_actual / C_cold
        
        dP_hot = calculate_pressure_drop(params['m_hot'], params['rho_hot'], params['mu_hot'], params['D_tube'], params['L'], 
                                         params['N_tube_passes'], params['N_tubes'], params['D_shell'], is_shell_side=False)
        dP_cold = calculate_pressure_drop(params['m_cold'], params['rho_cold'], params['mu_cold'], params['D_shell'] - params['D_tube'], 
                                          params['L'], params['N_shell_passes'], params['N_tubes'], params['D_shell'], params['N_baffles'], is_shell_side=True)
        
        print(f"Hot fluid outlet temperature: {T_hot_out:.2f} °C")
        print(f"Cold fluid outlet temperature: {T_cold_out:.2f} °C")
        print(f"Heat transfer rate: {Q_actual:.2f} W")
        print(f"Effectiveness: {effectiveness:.4f}")
        print(f"Overall heat transfer coefficient: {U:.2f} W/m²K")
        print(f"Pressure drop (hot side): {dP_hot:.2f} Pa")
        print(f"Pressure drop (cold side): {dP_cold:.2f} Pa")
        
        hot_phase_in = determine_phase(params['fluid_hot'], params['T_hot_in'], params['P_hot_in'])
        hot_phase_out = determine_phase(params['fluid_hot'], T_hot_out, params['P_hot_in'] - dP_hot)
        cold_phase_in = determine_phase(params['fluid_cold'], params['T_cold_in'], params['P_cold_in'])
        cold_phase_out = determine_phase(params['fluid_cold'], T_cold_out, params['P_cold_in'] - dP_cold)
        
        print(f"Hot fluid phase: In: {hot_phase_in}, Out: {hot_phase_out}")
        print(f"Cold fluid phase: In: {cold_phase_in}, Out: {cold_phase_out}")
        
        plt.figure(figsize=(10, 6))
        if params['flow_type'] == "counter":
            plt.plot([0, 1], [params['T_hot_in'], T_hot_out], 'r-', label='Hot fluid')
            plt.plot([0, 1], [T_cold_out, params['T_cold_in']], 'b-', label='Cold fluid')
        else:
            plt.plot([0, 1], [params['T_hot_in'], T_hot_out], 'r-', label='Hot fluid')
            plt.plot([0, 1], [params['T_cold_in'], T_cold_out], 'b-', label='Cold fluid')
        plt.xlabel('Normalized heat exchanger length')
        plt.ylabel('Temperature (°C)')
        plt.title('Temperature Profile in Heat Exchanger')
        plt.legend()
        plt.grid(True)
        plt.show()
    
    else:
        print("Invalid mode selected.")

if __name__ == "__main__":
    mode = input("Enter mode (design/analysis): ").lower()

    if mode == "design":
        params = {
            'm_hot': float(input("Hot fluid mass flow rate (kg/s): ")),
            'cp_hot': float(input("Hot fluid specific heat capacity (J/kg·K): ")),
            'T_hot_in': float(input("Hot fluid inlet temperature (°C): ")),
            'T_hot_out': float(input("Hot fluid desired outlet temperature (°C): ")),
            'T_cold_in': float(input("Cold fluid inlet temperature (°C): ")),
            'T_cold_out': float(input("Cold fluid desired outlet temperature (°C): ")),
            'flow_type': input("Flow type (counter/parallel): ")
        }
    elif mode == "analysis":
        params = {
            'fluid_hot': input("Hot fluid (e.g., water, air, steam): "),
            'fluid_cold': input("Cold fluid (e.g., water, air, steam): "),
            'm_hot': float(input("Hot fluid mass flow rate (kg/s): ")),
            'm_cold': float(input("Cold fluid mass flow rate (kg/s): ")),
            'T_hot_in': float(input("Hot fluid inlet temperature (°C): ")),
            'T_cold_in': float(input("Cold fluid inlet temperature (°C): ")),
            'P_hot_in': float(input("Hot fluid inlet pressure (Pa): ")),
            'P_cold_in': float(input("Cold fluid inlet pressure (Pa): ")),
            'D_shell': float(input("Shell diameter (m): ")),
            'D_tube': float(input("Tube diameter (m): ")),
            'L': float(input("Heat exchanger length (m): ")),
            'N_tubes': int(input("Number of tubes: ")),
            'flow_type': input("Flow type (counter/parallel): "),
            'N_baffles': int(input("Number of baffles: ")),
            'N_shell_passes': int(input("Number of shell passes: ")),
            'N_tube_passes': int(input("Number of tube passes: "))
        }
        hot_props = get_fluid_properties(params['fluid_hot'], params['T_hot_in'], params['P_hot_in'])
        cold_props = get_fluid_properties(params['fluid_cold'], params['T_cold_in'], params['P_cold_in'])

        v_hot = params['m_hot'] / (hot_props['density'] * (math.pi * params['D_tube']**2 / 4) * params['N_tubes'])
        v_cold = params['m_cold'] / (cold_props['density'] * (math.pi * (params['D_shell']**2 - params['N_tubes'] * params['D_tube']**2) / 4))

        Re_hot = hot_props['density'] * v_hot * params['D_tube'] / hot_props['viscosity']
        Re_cold = cold_props['density'] * v_cold * (params['D_shell'] - params['D_tube']) / cold_props['viscosity']

        Pr_hot = hot_props['specific_heat'] * hot_props['viscosity'] / hot_props['conductivity']
        Pr_cold = cold_props['specific_heat'] * cold_props['viscosity'] / cold_props['conductivity']

        params.update({
            'cp_hot': hot_props['specific_heat'],
            'cp_cold': cold_props['specific_heat'],
            'rho_hot': hot_props['density'],
            'rho_cold': cold_props['density'],
            'mu_hot': hot_props['viscosity'],
            'mu_cold': cold_props['viscosity'],
            'k_hot': hot_props['conductivity'],
            'k_cold': cold_props['conductivity'],
            'Re_hot': Re_hot,
            'Re_cold': Re_cold,
            'Pr_hot': Pr_hot,
            'Pr_cold': Pr_cold
        })
    else:
        print("Invalid mode selected.")
        exit()

    heat_exchanger_simulation(mode, **params)