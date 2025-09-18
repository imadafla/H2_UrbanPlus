# H2 UrbanPlus

![H2 UrbanPlus Logo](https://raw.githubusercontent.com/imadafla/H2_UrbanPlus/refs/heads/main/logo/H2_UrbanPlus.png)

H2 UrbanPlus is a Python library for simulating hydrogen-based energy and water recovery systems.  
It integrates PV generation, electrolysis, hydrogen storage, fuel cells, and water recovery, making it suitable for building and district energy studies.

Developed by Researchers in [Green Energy Park](https://www.greenenergypark.ma/): [Dr. Imad AIT LAASRI](https://imadcv.vercel.app/) and [Dr. Niima ES-SAKALI](https://www.linkedin.com/in/niimaes-sakali/)

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) after cloning:

```bash
git clone https://github.com/imadafla/H2_UrbanPlus.git
cd H2_UrbanPlus
pip install numpy pandas
```

## Usage `🚀 Quick Start`

### 1. Import the model

```python
from h2_urbanplus import H2_urbanplus
```

### 2. Prepare inputs

The model requires two main time series (hourly resolution recommended):

- PV generation (kWh): Hourly electricity produced by your PV system.

- Electricity demand (kWh): Hourly building/district demand from EnergyPlus.

### 3. Run the simulation

```python
# Example: load EnergyPlus and PV outputs
from h2_urbanplus import H2_urban

df_demand = pd.read_csv("scenario_demand.csv")
df_pv = pd.read_csv("pv_energy.csv")
df_pv_w = pd.read_csv("pv_power.csv")

results = H2_urbanplus(df_demand, df_pv, df_pv_w,
                       energy_col="Total Energy",
                       water_col="Water:Facility",
                       pv_col="Generated PV Energy",
                       pv_w_col="Generated PV Power")
print(results.head())

```

### 4. Outputs

- The function returns a pandas.DataFrame with:

- `Hydrogen_Produced_kg` – hourly hydrogen production

- `Hydrogen_Stored_kg` – hydrogen stored in the tank

- `Hydrogen_Used_kg` – hydrogen consumed by the fuel cell

- `Fuel Cell Generated Energy_kWh` – electricity generated from hydrogen

- `Electrolysis Heat Gain_kWh` & `FuelCell Heat Gain_kWh` – heat gains

- `Tap_Water_Used_kg` & `Water_Produced_in_FuelCell_kg` – water balance

- `Useful_Energy_After_Storage_kWh` – final usable energy after storage


## 🔗 EnergyPlus Coupling Guidelines
### 1- PV Integration

- Use `Generated PV Power (W)` and integrate to obtain `Generated PV Energy (Wh)` 
- `Both are required for the model`

### 2- Energy Demand Extraction

- If no district cooling/heating, you can directly use `Electricity:Facility`.

- If district heating/cooling is used, you must sum electricity, cooling, and heating demand.

- For ideal air load systems, sum heating and cooling demand, then adjust by COP/EER before adding to electricity demand.

### 3- Water Demand

- Extract `Water:Facility` output from EnergyPlus (in m³).

## ⚙️Customization

You can easily change parameters like:

- Electrolyzer/fuel cell efficiency (`η_elec`, `η_fc`)

- Max power capacities (`P_elec_max`, `P_fc_max`)

- Hydrogen tank size (`H2_tank_max`)

- Water tank and recovery efficiency

- Just pass them as keyword arguments in `H2_urbanplus()`.

The `H2_urbanplus` function accepts the following arguments:

```python
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
    heating_eff=0.9               # Heating system efficiency (0–1)
):
```

## License

H2 UrbanPlus is released under a **dual-license model**:

- **Academic / Non-commercial use**:  
  Licensed under [CC-BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) or [GNU AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html).  
  - Free to use for research and personal projects.  
  - Citation is required in academic publications.  
  - AGPLv3 ensures that any modifications shared over a network must also be open.  
