# This file contains the core classes used in the green chemistry simulator.

import numpy as np
from scipy.integrate import solve_ivp

class ChemicalSpecies:
    """A data class representing a specific chemical entity.
    
    Attributes:
        name (str): The unique identifier for the species (e.g., 'NO', 'H2O').
        vdw_a (float): Van der Waals 'a' constant (dm^6 atm mol^-2).
        vdw_b (float): Van der Waals 'b' constant (dm^3 mol^-1)."""
    
    def __init__(self, name, vdw_a=0.0, vdw_b=0.0):
        self.name = name
        self.vdw_a = vdw_a
        self.vdw_b = vdw_b

    def __repr__(self):
        return f"Species('{self.name}')"
    
class Reaction:
    """Represents a single, elementary reaction step (uni-directional).
    To model a reversible equilibrium, instantiate two Reaction objects: one for forward, one for reverse."""
    
    def __init__(self, reactants, products, A, Ea):
        self.reactants = reactants
        self.products = products
        self.A = A
        self.Ea = Ea
        self.R = 8.314 # J mol^-1 K^-1

    def get_rate_constant(self, T):
        """
        Calculates k = A * exp(-Ea / RT).
        Handles both scalar floats and numpy arrays safely.
        """
        # 1. Convert to array (even if scalar) to use numpy features...
        T_vals = np.asarray(T)
        
        # 2. Vectorized calculation with safety checks...
        # np.errstate ignores division by zero warnings if T=0 (we fix it in step 3)
        with np.errstate(divide='ignore', invalid='ignore'):
            k = self.A * np.exp(-self.Ea / (self.R * T_vals))
            
        # 3. Clean up invalid values (T <= 0)...
        # If T is scalar, this returns a 0-d array. If T is array, returns array.
        k = np.where(T_vals <= 0, 0.0, k)
        
        # 4. Return correct type...
        if T_vals.ndim == 0:
            return float(k)
        return k

    def calculate_rate(self, concentrations, T):
        """
        Calculates Rate = k(T) * Product([Conc]^Order)
        """
        k = self.get_rate_constant(T)
        rate = k
        for species_name, order in self.reactants.items():
            conc = concentrations.get(species_name, 0.0)
            rate *= (conc ** order)
        return rate

    def __repr__(self):
        reac_str = " + ".join([f"{v}{k}" for k, v in self.reactants.items()])
        prod_str = " + ".join([f"{v}{k}" for k, v in self.products.items()])
        return f"{reac_str} -> {prod_str}"
    
class GeneralChemicalSystem:
    """The central engine that models the state of a chemical system.
    It holds the species, the reactions, and the current physical conditions (Moles, V, T).
    
    Its primary job is to compute the derivative (dn/dt) for every species by summing the contributions of all reactions in the network."""
    
    def __init__(self, species_list, reaction_list, initial_moles, initial_V, initial_T):
        """Args:
            species_list (list[ChemicalSpecies]): List of all valid species objects.
            reaction_list (list[Reaction]): List of all reaction steps.
            initial_moles (dict): Mapping 'SpeciesName' -> Moles.
            initial_V (float): Initial Volume in dm^3.
            initial_T (float): Initial Temperature in Kelvin."""
        
        self.species_list = species_list
        self.reactions = reaction_list
        self.initial_moles = initial_moles
        self.V = initial_V
        self.T = initial_T
        self.species_names = {s.name for s in species_list}
        self.species_order = sorted(list(self.species_names))

        # validation checks:
        self._validate_reactions()
        for name in initial_moles:
            if name not in self.species_names:
                raise ValueError(f"Unknown species '{name}' found in initial_moles.")
                
        # placeholders for results
        self.solution = None
        self.results = None

    def _validate_reactions(self):
        """Internal method to ensure reaction definitions match the species list."""
        for i, reaction in enumerate(self.reactions):
            # check reactants
            for r_name in reaction.reactants:
                if r_name not in self.species_names:
                    raise ValueError(f"Reaction #{i+1} contains unknown reactant: '{r_name}'")
            # check products
            for p_name in reaction.products:
                if p_name not in self.species_names:
                    raise ValueError(f"Reaction #{i+1} contains unknown product: '{p_name}'")
                    
    def _state_to_array (self, moles_dict):
        """Converts a {Species: Moles} dictionary to a [Moles] numpy array."""
        return np.array([moles_dict.get(s, 0.0) for s in self.species_order])

    def _array_to_state_dict (self, moles_array):
        """Converts a [Moles] numpy array to a {Species: Moles} dictionary."""
        return {s: moles_array[i] for i, s in enumerate(self.species_order)}

    def get_net_rates(self, current_moles, V, T):
        """Calculates the net rate of change (dn/dt) for ALL species in the system.
        
        Logic:
        dn_i/dt = Sum( Stoichiometry_i,j * Rate_j * Volume )
        
        Args:
            current_moles (dict): Mapping 'SpeciesName' -> Moles.
            V (float): Current volume (dm^3).
            T (float): Current temperature (K).
            
        Returns:
            dict: Mapping 'SpeciesName' -> dn/dt (mol s^-1)."""
        
        # 1. Initialize derivatives for all known species to 0.0...
        derivatives = {name: 0.0 for name in self.species_names}
        
        # 2. Calculate concentrations (c = n/V)...
        # needed because the rate aw is defined in terms of concentration
        # we default to 0.0 if a species is missing from the input dict (safety)
        concentrations = {
            name: current_moles.get(name, 0.0) / V 
            for name in self.species_names
        }
        
        # 3. Summation loop...
        for reaction in self.reactions:
            # A. Calculate in mol dm^-3 s^-1.
            # rate = k(T) * [A]^a * [B]^b
            rate_concentration = reaction.calculate_rate(concentrations, T)
            
            # B. Convert to mol s^-1.
            # Reaction happens in the whole volume V.
            rate_moles_per_sec = rate_concentration * V
            
            # C. Update reactants (consumption -> negative change).
            for r_name, coeff in reaction.reactants.items():
                # dn/dt -= coefficient * rate
                derivatives[r_name] -= coeff * rate_moles_per_sec
            
            # D. Update products (formation -> positive change).
            for p_name, coeff in reaction.products.items():
                # dn/dt += coefficient * rate
                derivatives[p_name] += coeff * rate_moles_per_sec
                
        return derivatives

    def _ode_system_adapter(self, t, y):
        # 1. Safety clamp: prevent negative moles creating math errors...
        y = np.maximum(y, 0)
        
        # 2. Convert array -> dict...
        current_moles_dict = self._array_to_state_dict(y)
        
        # 3. Calculate rates...
        # Note: we use the system's current V and T. 
        # in a generic solver, these are constant unless perturbed externally
        rates_dict = self.get_net_rates(current_moles_dict, self.V, self.T)
        
        # 4. Convert dict -> array...
        return self._state_to_array(rates_dict)

    def run_simulation(self, rate_tolerance=1e-7, max_iterations=20):
        """Runs the simulation using an iterative, adaptive timescale. Stops when the maximum net rate falls below rate_tolerance."""
        # 1. Setup initial state...
        y0 = self._state_to_array(self.initial_moles)
        
        solutions_list = []
        time_offset = 0.0
        current_chunk_duration = 0.1 # start small for fast initial kinetics
        equilibrium_found = False
        
        print(f"Starting simulation at T={self.T}K...")

        # 2. The iterative loop...
        for i in range(max_iterations):
            t_span = (time_offset, time_offset + current_chunk_duration)
            
            # solve for this chunk, Radau is robust for stiff systems
            chunk_solution = solve_ivp(
                fun=self._ode_system_adapter, t_span=t_span, y0=y0, dense_output=True, 
                method='Radau', rtol=1e-6, atol=1e-9)
            
            solutions_list.append(chunk_solution)
            
            # update state for next chunk
            y_last = chunk_solution.y[:, -1]
            y0 = y_last
            time_offset = chunk_solution.t[-1]
            
            # 3. Equilibrium check...
            # convert final state to dict to calculate rates
            final_moles = self._array_to_state_dict(y_last)
            final_rates = self.get_net_rates(final_moles, self.V, self.T)
            max_rate = max(abs(v) for v in final_rates.values())
            
            if max_rate < rate_tolerance:
                print(f"Equilibrium reached at t={time_offset:.4f}s. Max Rate: {max_rate:.2e}")
                equilibrium_found = True
                self.equilibrium_time = time_offset
                break
            
            # if not found, extend simulation time geometrically
            # we assume if it didn't finish in 0.1s, it might need 1.0s, then 10.0s...
            current_chunk_duration *= 10
            
        if not equilibrium_found:
            print(f"Warning: Equilibrium not reached after {time_offset:.2f}s (Max Iterations).")

        # 4. Process and store results...
        self._process_final_results(solutions_list)

    def _process_final_results(self, solutions_list):
        """Stitches the chunks together and structures the data."""
        combined_t = np.concatenate([s.t for s in solutions_list])
        combined_y = np.concatenate([s.y for s in solutions_list], axis=1)
        
        # organise data by species name for easy plotting
        species_data = {}
        for idx, name in enumerate(self.species_order):
            species_data[name] = combined_y[idx]
            
        # store as dictionary
        self.results = {
            'time': combined_t,
            'species_data': species_data,
            'final_moles': self._array_to_state_dict(combined_y[:, -1]),
            'info': {'V': self.V, 'T': self.T},
            'equilibrium_time': getattr(self, 'equilibrium_time', None)}

    def calculate_pressure_array(self, moles_data, V_array, T_array):
        """
        Calculates Real Pressure for the entire simulation history.
        moles_data: dict of {species: np.array}
        V_array: np.array of volumes
        T_array: np.array of temperatures
        """
        R_atm = 0.08206
        
        # 1. Calculate n_total for every time step...
        n_total = np.zeros_like(V_array)
        for s in self.species_order:
            n_total += moles_data[s]
            
        # 2. Mixing rules (vectorised)...
        sqrt_a_sum = np.zeros_like(V_array)
        b_sum = np.zeros_like(V_array)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            for s_obj in self.species_list:
                n_i = moles_data[s_obj.name]
                x_i = np.divide(n_i, n_total, where=n_total!=0)
                
                sqrt_a_sum += x_i * np.sqrt(s_obj.vdw_a)
                b_sum += x_i * s_obj.vdw_b
                
        a_mix = sqrt_a_sum ** 2
        b_mix = b_sum
        
        # 3. Van der Waals calculation...
        nb = n_total * b_mix
        V_free = V_array - nb
        
        # safety clamp for V < nb (physical impossibility)
        V_free = np.maximum(V_free, 1e-6) 
        
        P_ideal = (n_total * R_atm * T_array) / V_free
        P_real = P_ideal - ((a_mix * n_total**2) / (V_array**2))
        
        # handle cases where n_total is 0
        P_real = np.where(n_total == 0, 0, P_real)
        
        return P_real

    def _process_final_results(self, solutions_list):
        combined_t = np.concatenate([s.t for s in solutions_list])
        combined_y = np.concatenate([s.y for s in solutions_list], axis=1)
        
        combined_y = np.maximum(combined_y, 0.0) # ensure moles aren't negative
        
        species_data = {}
        for idx, name in enumerate(self.species_order):
            species_data[name] = combined_y[idx]
            
        V_array = np.full_like(combined_t, self.V)
        T_array = np.full_like(combined_t, self.T)
        
        P_real = self.calculate_pressure_array(species_data, V_array, T_array)
        
        # calculate ideal P for reference
        n_total = sum(species_data.values())
        R_atm = 0.08206
        P_ideal = (n_total * R_atm * T_array) / V_array
        
        self.results = {
            'time': combined_t,
            'species_data': species_data,
            'final_moles': self._array_to_state_dict(combined_y[:, -1]),
            'info': {'V': self.V, 'T': self.T},
            'P_real': P_real,
            'P_ideal': P_ideal
        }

    def calculate_thermodynamics(self, fwd_idx, rev_idx, external_results=None):
        """
        Calculates Qc/Kc (concentration) and Qp/Kp (pressure).
        """
        # get reaction objects...
        r_fwd = self.reactions[fwd_idx]
        r_rev = self.reactions[rev_idx]
        
        # 1. Select data source.
        data = external_results if external_results else self.results
        if data is None: return None
        
        T_arr = data.get('temperature', np.full_like(data['time'], self.T))
        V_arr = data.get('volume', np.full_like(data['time'], self.V))
        P_real = data.get('P_real', np.zeros_like(T_arr))
        species_data = data['species_data']
        
        # 2. Calculate Kc.
        k_f = r_fwd.get_rate_constant(T_arr)
        k_r = r_rev.get_rate_constant(T_arr)
        Kc = np.divide(k_f, k_r, out=np.zeros_like(k_f), where=k_r!=0)
        
        # 3. Calculate Qc.
        num = np.ones_like(T_arr)
        for prod, coeff in r_fwd.products.items():
            conc = species_data[prod] / V_arr
            num *= (conc ** coeff)
            
        denom = np.ones_like(T_arr)
        for reac, coeff in r_fwd.reactants.items():
            conc = species_data[reac] / V_arr
            denom *= (conc ** coeff)
            
        Qc = np.divide(num, denom, out=np.zeros_like(num), where=denom > 1e-20)
        
        # 4. Calculate Kp/Qp.
        moles_prod = sum(r_fwd.products.values())
        moles_reac = sum(r_fwd.reactants.values())
        delta_n = moles_prod - moles_reac
        R_atm = 0.08206
        
        Kp = Kc * (R_atm * T_arr)**delta_n
        
        # Qp (calculated based on partial pressures)
        n_total = np.zeros_like(T_arr)
        for name in self.species_names: n_total += species_data.get(name, 0)
            
        num_p = np.ones_like(T_arr)
        for prod, coeff in r_fwd.products.items():
            x_i = np.divide(species_data[prod], n_total, out=np.zeros_like(P_real), where=n_total!=0)
            num_p *= ((x_i * P_real) ** coeff)
            
        denom_p = np.ones_like(T_arr)
        for reac, coeff in r_fwd.reactants.items():
            x_i = np.divide(species_data[reac], n_total, out=np.zeros_like(P_real), where=n_total!=0)
            denom_p *= ((x_i * P_real) ** coeff)
            
        Qp = np.divide(num_p, denom_p, out=np.zeros_like(num_p), where=denom_p > 1e-20)
        
        # safety clamp: prevent 0 values which disappear on log scales
        min_val = 1e-20
        Qc = np.maximum(Qc, min_val)
        Qp = np.maximum(Qp, min_val)
        
        return {'Qc': Qc, 'Kc': Kc, 'Qp': Qp, 'Kp': Kp}