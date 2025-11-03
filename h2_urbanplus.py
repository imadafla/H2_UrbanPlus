import numpy as np
import pandas as pd

def H2_urbanplus(
    df_demand,
    df_pv,
    df_pv_w,
    energy_col="Total Energy",
    water_col="Water:Facility",
    pv_col="Generated PV Energy",
    pv_w_col="Generated PV Power",
    η_elec=0.6,                   # Electrolysis efficiency (0–1)
    η_fc=0.72,                    # Fuel cell efficiency (0–1)
    HHV_H2=33.33,                 # kWh/kg, Hydrogen higher heating value
    P_elec_max=100,               # Max electrolysis power (kW)
    P_fc_max=50,                  # Max fuel cell power (kW)
    H2_tank_max=20,               # Hydrogen storage capacity (kg)
    W_tank_max=0,                 # Water tank capacity (kg)
    η_recov=0.4,                  # Wastewater recovery efficiency (0–1)
    water_density=1000,           # kg/m³
    electrolysis_heat_factor=0.3, # Heat fraction from electrolysis
    fuel_cell_heat_factor=0.2,    # Heat fraction from fuel cell
    Water_utilization_ratio=0.5,  # Fraction of water drained (0–1)
    c_p_water=4.18,               # kJ/kg·°C (specific heat of water)
    T_in=20,                      # °C, inlet temperature
    T_target=80,                  # °C, target electrolysis water temperature
    heating_eff=0.9,              # Heating system efficiency (0–1)
    *,
    params: dict | None = None
):
    """
    Hydrogen-based urban energy system simulator.
    Developed by Green Energy Park Researchers (Dr. Imad AIT LAASRI & Dr. Niima ES-SAKALI)

    Parameters
    ----------
    df_demand : pandas.DataFrame
        EnergyPlus demand outputs including total energy consumption and water facility use.
    df_pv : pandas.DataFrame
        PV generation (integrated energy output).
    df_pv_w : pandas.DataFrame
        PV generation (instantaneous power output).
    energy_col : str
        Column name for building/district energy demand (Joules).
    water_col : str
        Column name for water demand (m³).
    pv_col : str
        Column name for PV generation energy (Joules).
    pv_w_col : str
        Column name for PV power (Watts).
    params : dict, optional
        Dictionary of overrides for any named parameter above. You can pass only the keys
        you want to change; the rest will keep their default values. Supports ASCII aliases
        for Greek-letters:
            - 'η_elec' or 'eta_elec'
            - 'η_fc'   or 'eta_fc'
            - 'η_recov' or 'eta_recov'

    Returns
    -------
    pandas.DataFrame
        Simulation results with hydrogen production, storage, fuel cell use, water balance, and energy flows.

    Example
    -------
    params = {
        'η_elec': 0.7,
        'η_fc': 0.6,
        'P_elec_max': 100,
        'P_fc_max': 100,
        'H2_tank_max': 1000
    }
    results = H2_urbanplus(df_demand, df_pv, df_pv_w, params=params)
    """

    # ----- Parameter overrides via `params` ---------------------------------
    if params:
        # Allow ASCII aliases for convenience
        alias_to_canonical = {
            'eta_elec': 'η_elec',
            'eta_fc': 'η_fc',
            'eta_recov': 'η_recov',
        }
        # Build the set of valid (canonical) keys
        valid_keys = {
            'energy_col', 'water_col', 'pv_col', 'pv_w_col',
            'η_elec', 'η_fc', 'HHV_H2', 'P_elec_max', 'P_fc_max',
            'H2_tank_max', 'W_tank_max', 'η_recov', 'water_density',
            'electrolysis_heat_factor', 'fuel_cell_heat_factor',
            'Water_utilization_ratio', 'c_p_water', 'T_in', 'T_target',
            'heating_eff'
        }

        # Normalize keys, check for typos
        normalized = {}
        for k, v in params.items():
            k_norm = alias_to_canonical.get(k, k)
            if k_norm not in valid_keys:
                raise ValueError(
                    f"Unknown parameter '{k}'. "
                    f"Allowed keys: {sorted(valid_keys)} "
                    f"(aliases: {list(alias_to_canonical.keys())})"
                )
            normalized[k_norm] = v

        # Apply overrides (only those provided)
        energy_col = normalized.get('energy_col', energy_col)
        water_col = normalized.get('water_col', water_col)
        pv_col = normalized.get('pv_col', pv_col)
        pv_w_col = normalized.get('pv_w_col', pv_w_col)
        η_elec = normalized.get('η_elec', η_elec)
        η_fc = normalized.get('η_fc', η_fc)
        HHV_H2 = normalized.get('HHV_H2', HHV_H2)
        P_elec_max = normalized.get('P_elec_max', P_elec_max)
        P_fc_max = normalized.get('P_fc_max', P_fc_max)
        H2_tank_max = normalized.get('H2_tank_max', H2_tank_max)
        W_tank_max = normalized.get('W_tank_max', W_tank_max)
        η_recov = normalized.get('η_recov', η_recov)
        water_density = normalized.get('water_density', water_density)
        electrolysis_heat_factor = normalized.get('electrolysis_heat_factor', electrolysis_heat_factor)
        fuel_cell_heat_factor = normalized.get('fuel_cell_heat_factor', fuel_cell_heat_factor)
        Water_utilization_ratio = normalized.get('Water_utilization_ratio', Water_utilization_ratio)
        c_p_water = normalized.get('c_p_water', c_p_water)
        T_in = normalized.get('T_in', T_in)
        T_target = normalized.get('T_target', T_target)
        heating_eff = normalized.get('heating_eff', heating_eff)
    # ------------------------------------------------------------------------

    # Convert EnergyPlus outputs
    E_demand = df_demand[energy_col].values / 3.6e6  # J → kWh
    W_demand_m3 = df_demand[water_col].values
    E_pv = df_pv[pv_col].values / 3.6e6  # J → kWh
    P_pv = df_pv_w[pv_w_col].values / 1000  # W → kW

    n_hours = len(E_demand)

    # Initialize arrays
    m_H2_stored = np.zeros(n_hours)
    W_tank = np.zeros(n_hours)
    Tap_Water = np.zeros(n_hours)
    W_lost = np.zeros(n_hours)
    E_fc = np.zeros(n_hours)
    m_H2_prod = np.zeros(n_hours)
    m_H2_used = np.zeros(n_hours)
    E_final = np.zeros(n_hours)
    Q_elec = np.zeros(n_hours)
    Q_fc = np.zeros(n_hours)
    E_demand_ideal = np.zeros(n_hours)
    P_net = np.zeros(n_hours)
    water_needed = np.zeros(n_hours)
    P_elec_fraction = np.zeros(n_hours)
    W_recovery = np.zeros(n_hours)
    P_electrolysis = np.zeros(n_hours)
    P_electrolysis_heating = np.zeros(n_hours)
    water_produced_fc = np.zeros(n_hours)

    for t in range(n_hours):
        # Water recovery
        W_recovery[t] = η_recov * W_demand_m3[t] * water_density * Water_utilization_ratio
        prev_W_tank = W_tank[t-1] if t > 0 else 0
        W_tank[t] = min(prev_W_tank + W_recovery[t], W_tank_max)
        W_lost[t] = max(0, prev_W_tank + W_recovery[t] - W_tank_max)

        # Net balance (per-hour energy; numerically equal to average kW over the hour)
        P_net[t] = E_pv[t] - E_demand[t]

        if P_net[t] > 0:
            # Electrolysis
            P_elec_required = min(P_elec_max, P_net[t], P_pv[t])
            P_elec_fraction[t] = P_net[t] / P_pv[t] if P_pv[t] > 0 else 0

            m_H2 = η_elec * P_elec_required / HHV_H2
            water_needed[t] = 9 * m_H2  # kg of water per kg of H2 (stoichiometry)

            if water_needed[t] > 0:
                Q_heat_water_kWh = (water_needed[t] * c_p_water * (T_target - T_in)) / (3600 * heating_eff)
            else:
                Q_heat_water_kWh = 0

            if W_tank[t] >= water_needed[t]:
                W_used = water_needed[t]
                Tap_Water[t] = 0
                W_tank[t] -= W_used
            else:
                W_used = W_tank[t]
                Tap_Water[t] = water_needed[t] - W_used
                W_tank[t] = 0

            prev_H2 = m_H2_stored[t-1] if t > 0 else 0
            m_H2_stored[t] = min(prev_H2 + m_H2, H2_tank_max)
            m_H2_prod[t] = m_H2
            P_electrolysis[t] = P_elec_required
            P_electrolysis_heating[t] = Q_heat_water_kWh
            Q_elec[t] = electrolysis_heat_factor * P_elec_required

        elif P_net[t] < 0:
            # Fuel Cell
            P_fc_req = min(abs(P_net[t]), P_fc_max)
            H2_needed = P_fc_req / (η_fc * HHV_H2)
            prev_H2 = m_H2_stored[t-1] if t > 0 else 0

            if prev_H2 >= H2_needed:
                m_H2_used[t] = H2_needed
                E_fc[t] = η_fc * H2_needed * HHV_H2
                m_H2_stored[t] = prev_H2 - H2_needed
            else:
                m_H2_used[t] = prev_H2
                E_fc[t] = η_fc * m_H2_used[t] * HHV_H2
                m_H2_stored[t] = 0

            Q_fc[t] = fuel_cell_heat_factor * E_fc[t]
            water_produced_fc[t] = m_H2_used[t] * (18 / 2)  # kg of water from H2 (2 g/mol -> 18 g/mol H2O)

        else:
            P_elec_fraction[t] = 0
            m_H2_stored[t] = m_H2_stored[t-1] if t > 0 else 0

        # Final accounting
        E_final[t] = (E_pv[t] * (1 - P_elec_fraction[t])) + E_fc[t]

    return pd.DataFrame({
        "Hour": np.arange(n_hours),
        "Hydrogen_Produced_kg": m_H2_prod,
        "Hydrogen_Stored_kg": m_H2_stored,
        "Hydrogen_Used_kg": m_H2_used,
        "Electrolysis_Heat_Gain_kWh": Q_elec,
        "FuelCell_Heat_Gain_kWh": Q_fc,
        "Electrolysis Heat Consummed_kW": P_electrolysis_heating,  # kept original label
        "Energy_Consumption_kWh": E_demand,
        "Generated_PV_Power_kWh": E_pv,
        "P_net_kW": P_net,
        "Electrolysis Energy Consummed_kW": P_electrolysis,        # kept original label
        "Fuel Cell Generated Energy_kWh": E_fc,
        "Useful_Energy_After_Storage_kWh": E_final,
        "Water_Waste_filtered_kg": W_recovery,
        "Tap_Water_Used_kg": Tap_Water,
        "Water_Needed_kg": water_needed,
        "Water_Produced_in_FuelCell_kg": water_produced_fc
    })
