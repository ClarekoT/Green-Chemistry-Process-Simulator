# src/simulation_engine/metrics.py
import numpy as np
class MetricsCalculator:
    """
    Calculates a suite of green chemistry metrics from simulation outputs and metadata.
    This class bridges the gap between the raw physical simulation (moles, time, T)
    and the higher-level green chemistry analysis (PMI, AE, TON, etc.).
    """

    def __init__(self, system_object, initial_metadata):
        """
        Initialises the calculator with data from a completed simulation run.

        Args:
            simulation_results (dict): The complete, final results dictionary produced by a GeneralChemicalSystem instance.
            initial_metadata (dict): A dictionary containing supplementary data not available in the simulation output.
                Expected keys:
                     - 'desired_product_name' (str)
                     - 'product_molecular_weight' (float)
                     - 'catalyst_moles' (float)
                     - 'solvent_scores' (dict): {'SolventName': score}
                     - 'feedstock_scores' (dict): {'FeedstockName': score}
                     - 'q_values' (dict): {'SpeciesName': Q_value}
                     - 'reactant_molecular_weights' (dict): {'ReactantName': MW}
                     - 'initial_masses' (dict): {'SpeciesName': mass_in_grams}
                     - 'overall_reaction_stoichiometry' (dict): {
                        'reactants': {'ReactantName': x, 'ReactantName': y},
                        'products': {'ProductName': z}}
        """
        self.system = system_object # store the whole system
        self.sim_results = system_object.results
        self.metadata = initial_metadata

        # basic validation to ensure essential keys exist
        required_keys = [
            'desired_product_name', 'catalyst_moles', 'solvent_scores',
            'feedstock_scores', 'product_molecular_weight', 'q_values',
            'reactant_molecular_weights', 'initial_masses', 'overall_reaction_stoichiometry'
        ]
        for key in required_keys:
            if key not in self.metadata:
                raise ValueError(f"Missing required key in initial_metadata: '{key}'")

    @property
    def _molecular_weights(self):
        """A helper property to access a unified dictionary of all known molecular weights."""
        if not hasattr(self, '_mw_dict'):
            # Create the dictionary on first access and cache it
            self._mw_dict = self.metadata.get('reactant_molecular_weights', {}).copy()
            prod_name = self.metadata['desired_product_name']
            prod_mw = self.metadata['product_molecular_weight']
            self._mw_dict[prod_name] = prod_mw
        return self._mw_dict

    # --- ATOM ECONOMY ---
    def _calculate_ae(self):
        """
        Calculates the Atom Economy based on the stoichiometry and molecular weights of the reactants and the desired product.
        - Molecular weights is taken from self.metadata.
        - Stoichiometry is contained with the Reaction objects.
        This method is robust for multi-step syntheses because it relies on an explicit definition of the net transformation provided in the metadata,
        rather than inferring it from kinetic elementary steps.
        Returns:
            float: calculate AE as a percentage (0-100).
        """
        stoichiometry = self.metadata.get('overall_reaction_stoichiometry')
        if not stoichiometry:
            raise ValueError("Cannot calculate AE: 'overall_reaction_stoichiometry' is missing from metadata.")

        product_name = self.metadata['desired_product_name']
        product_mw = self.metadata['product_molecular_weight']
        
        # numerator: (MW of desired product * its overall stoichiometric coefficient)
        product_stoich = stoichiometry['products'].get(product_name, 0)
        if product_stoich == 0:
            raise ValueError(f"Cannot calculate AE: Desired product '{product_name}' not found in 'overall_reaction_stoichiometry' products.")
        
        numerator = product_mw * product_stoich

        # denominator: sum of (MW of each reactant * its overall stoichiometric coefficient)
        denominator = 0
        reactant_mws = self.metadata['reactant_molecular_weights']
        for reactant_name, reactant_stoich in stoichiometry['reactants'].items():
            reactant_mw = reactant_mws.get(reactant_name)
            if reactant_mw is None:
                raise ValueError(f"Cannot calculate AE: Molecular weight not found in metadata for reactant '{reactant_name}'.")
            denominator += reactant_mw * reactant_stoich

        # final calculation with division-by-zero safety check.
        if denominator == 0: return 0.0

        atom_economy = (numerator / denominator) * 100
        return atom_economy

    def _score_ae(self, ae_value):
        """Converts a raw AE percentage into a 1-5 score based on project criteria."""
        if ae_value > 95: return 5
        elif 80 <= ae_value <= 95: return 4
        elif 60 <= ae_value < 80: return 3
        elif 40 <= ae_value < 60: return 2
        else: # ae_value < 40
            return 1
        
    # --- PROCESS MASS INTENSITY ---
    def _calculate_pmi(self):
        """
        Calculates the Process Mass Intensity (PMI).

        PMI is the ratio of the total mass of all materials (reactants, solvents, catalysts, etc.) put into a process to the mass of the final desired product.

        Returns:
            float: The calculated PMI value. Returns float('inf') if no product is formed.
        """
        # numerator: sum of all values in the initial_masses dictionary
        initial_masses = self.metadata.get('initial_masses')
        if not initial_masses or not isinstance(initial_masses, dict):
            raise ValueError("Cannot calculate PMI: 'initial_masses' dictionary is missing or invalid in metadata.")
        
        total_input_mass = sum(initial_masses.values())

        # denominator: calculate the final mass of the desired product from simulation results
        product_name = self.metadata['desired_product_name']
        product_mw = self.metadata['product_molecular_weight']
        
        final_moles_dict = self.sim_results.get('final_moles', {})
        final_product_moles = final_moles_dict.get(product_name, 0.0)
        
        final_product_mass = final_product_moles * product_mw

        # final calculation with division-by-zero safety check.
        # (if no product is made, the waste is effectively infinite)
        if final_product_mass <= 0:
            return float('inf')

        pmi = total_input_mass / final_product_mass
        return pmi
    
    def _score_pmi(self, pmi_value):
        """Converts a raw PMI value into a 1-5 score based on project criteria."""
        if pmi_value < 5: return 5 
        elif 5 <= pmi_value < 20: return 4
        elif 20 <= pmi_value < 100: return 3
        elif 100 <= pmi_value < 500: return 2
        else:  # pmi_value >= 500 or pmi_value == float('inf')
            return 1
    
    # --- ENVIRONMENTAL QUOTIENT ---
    def _calculate_eq(self):
        """
        Calculates the Environmental Quotient (EQ).

        EQ is a hazard-weighted measure of waste. It is calculated as the sum of
        the mass of each waste component multiplied by its hazard factor (Q-value),
        all divided by the mass of the desired product.

        Returns:
            float: The calculated EQ value. Returns float('inf') if no product is formed.
        """
        # denominator: calculate the final mass of the desired product
        product_name = self.metadata['desired_product_name']
        final_moles_dict = self.sim_results.get('final_moles', {})
        final_product_moles = final_moles_dict.get(product_name, 0.0)
        final_product_mass = final_product_moles * self._molecular_weights.get(product_name, 0)

        if final_product_mass <= 0:
            return float('inf')

        # numerator: calculate the sum of the weighted waste (Σ mass_i * Q_i)
        total_weighted_waste = 0.0
        q_values = self.metadata.get('q_values', {})

        # iterate through all species present at the end of the reaction
        for species_name, final_moles in final_moles_dict.items():
            # waste includes everything except the desired product
            if species_name == product_name:
                continue

            # get the molecular weight and Q-value for the waste component
            mw = self._molecular_weights.get(species_name)
            if mw is None:
                raise ValueError(f"Cannot calculate EQ: Molecular weight not found for waste component '{species_name}'.")
            
            # default to Q=1 (benign) if a specific Q-value is not provided
            q_value = q_values.get(species_name, 1.0)
            
            mass_of_component = final_moles * mw
            total_weighted_waste += mass_of_component * q_value

        # final calculation
        eq = total_weighted_waste / final_product_mass
        return eq
    
    def _score_eq(self, eq_value):
        """
        Converts a raw EQ value into a 1-5 score.
        """
        # This implementation scores the final aggregated value.
        if eq_value < 4: return 5
        elif 4 <= eq_value < 19: return 4
        elif 19 <= eq_value < 99: return 3 
        elif 99 <= eq_value < 499: return 2
        else:  # eq_value >= 499 or eq_value == float('inf')
            return 1
        
    # --- HAZARD SCORE ---
    def _calculate_hazard_score(self):
        """
        Calculates the Hazard Score based on the "weakest link" principle.

        This metric identifies the highest Q-value among all substances present in the final reaction mixture in any significant amount. It flags the
        single most hazardous material in the process.

        Returns:
            float: The highest Q-value found in the final mixture.
        """
        q_values = self.metadata.get('q_values', {})
        final_moles_dict = self.sim_results.get('final_moles', {})

        # find all species present at the end of the reaction (moles > 0)
        final_species_present = [
            name for name, moles in final_moles_dict.items() if moles > 1e-9 # use a small threshold
        ]

        if not final_species_present:
            return 1.0 # if the reactor is empty, assume benign

        # get the Q-value for each substance present, defaulting to 1 (benign)
        highest_q = 1.0
        for species_name in final_species_present:
            q_value = q_values.get(species_name, 1.0)
            if q_value > highest_q:
                highest_q = q_value
        
        return highest_q
    
    def _score_hazard_score(self, q_value):
        """
        Converts the highest Q-value (the Hazard Score) into a 1-5 score.
        This directly reflects the hazard level of the "weakest link" chemical.
        """
        if q_value <= 1: return 5  # benign
        elif 1 < q_value <= 10: return 4  # low hazard
        elif 10 < q_value <= 100:
            return 3  # moderate hazard
        elif 100 < q_value <= 1000: return 2  # high hazard
        else:  # q_value > 1000
            return 1  # extreme hazard
        
    # --- CATALYST TURNOVER NUMBER ---
    def _calculate_ton(self):
        """
        Calculates the Catalyst Turnover Number (TON).

        TON is the molar amount of product formed per mole of catalyst used. It is a measure of catalyst efficiency and lifetime.

        Returns:
            float: The calculated TON. Returns 0.0 if no catalyst is used.
        """
        # get the initial moles of the catalyst from metadata
        catalyst_moles = self.metadata.get('catalyst_moles')

        # if no catalyst is used, TON is not applicable, return 0
        if catalyst_moles is None or catalyst_moles <= 0:
            return 0.0

        # get the final moles of the desired product from simulation results
        product_name = self.metadata['desired_product_name']
        final_moles_dict = self.sim_results.get('final_moles', {})
        final_product_moles = final_moles_dict.get(product_name, 0.0)

        # final calculation
        ton = final_product_moles / catalyst_moles
        return ton
    
    def _score_ton(self, ton_value):
        """Converts a raw TON value into a 1-5 score."""
        if ton_value > 1_000_000: return 5  # world-class, bulk chemical production
        elif 10_000 <= ton_value <= 1_000_000: return 4  # highly robust and long-lived catalyst
        elif 1_000 <= ton_value < 10_000: return 3  # common target for an "industrially relevant" catalyst
        elif 100 <= ton_value < 1_000: return 2  # indicates a catalyst with a short lifetime or low activity
        else:  # ton_value < 100
            return 1  # barely catalytic, typical of an early-stage catalyst
        
    # --- SOLVENT ASSESSMENT SCORE ---
    def _get_solvent_score(self):
        """
        Determines the overall Solvent Score based on the "weakest link" principle.

        The score is the lowest score among all solvents listed in the metadata.
        A solvent-free process is considered ideal and receives the highest score.

        Returns:
            float: The lowest solvent score found (1-5).
        """
        solvent_scores = self.metadata.get('solvent_scores', {})
        
        # if the dictionary is empty, it's a solvent-free process, which is ideal
        if not solvent_scores:
            return 5.0
        
        # apply the "weakest link" principle: the process is only as green as its worst solvent
        return float(min(solvent_scores.values()))
    
    # --- RENEWABLE FEEDSTOCK SCORE ---
    def _calculate_feedstock_score(self):
        """
        Calculates a mass-weighted average score for feedstock sustainability.

        This method provides a more nuanced score than the "weakest link" by considering the mass contribution of each feedstock to the overall process.
        The score is Σ(mass_i * score_i) / Σ(mass_i).

        Returns:
            float: The mass-weighted average feedstock score.
        """
        initial_masses = self.metadata.get('initial_masses', {})
        feedstock_scores = self.metadata.get('feedstock_scores', {})
        stoichiometry = self.metadata.get('overall_reaction_stoichiometry', {})
        
        if not stoichiometry or 'reactants' not in stoichiometry:
             raise ValueError("Cannot calculate Feedstock Score: 'overall_reaction_stoichiometry' is missing from metadata.")

        feedstocks = stoichiometry['reactants'].keys()
        
        total_weighted_score = 0.0
        total_feedstock_mass = 0.0

        for name in feedstocks:
            mass = initial_masses.get(name)
            if mass is None:
                raise ValueError(f"Cannot calculate Feedstock Score: Initial mass for feedstock '{name}' not found in metadata.")
            
            score = feedstock_scores.get(name)
            if score is None:
                raise ValueError(f"Cannot calculate Feedstock Score: Score for feedstock '{name}' not found in metadata.")

            total_weighted_score += mass * score
            total_feedstock_mass += mass
        
        if total_feedstock_mass == 0:
            # this case happens if a reaction has reactants but all initial masses are zero.
            # return a neutral score of 3, as no judgment can be made.
            return 3.0

        return total_weighted_score / total_feedstock_mass
        
    def _score_feedstock(self, weighted_average_score):
        """
        Converts the weighted-average feedstock score into a final 1-5 integer score.
        
        The final score is determined by rounding the average to the nearest integer,
        ensuring the result is clamped between 1 and 5.
        """
        # round to the nearest integer
        score = int(round(weighted_average_score))
        # clamp the value to be within the 1-5 range
        return max(1, min(5, score))
    
    # --- Public Methods ---
    def calculate_endpoint_metrics(self):
        """
        Calculates all key green chemistry metrics for the final state of the reaction.
        This method provides a "snapshot" of the process's greenness at the point of completion.

        Returns:
            dict: A dictionary where each key is a metric name (e.g., 'PMI', 'AE')
                  and the value is another dictionary containing the 'raw_value' and the normalised 1-to-5 'score'.
        """
        results= {}

        # --- AE ---
        try:
            ae_raw = self._calculate_ae()
            ae_score = self._score_ae(ae_raw)
            results['AE'] = {'raw_value': round(ae_raw, 2), 'score': ae_score}
        except (RuntimeError, ValueError) as e:
            print(f"Could not calculate AE: {e}")
            results['AE'] = {'raw_value': None, 'score': None}

        # --- PMI ---
        try:
            pmi_raw = self._calculate_pmi()
            pmi_score = self._score_pmi(pmi_raw)
            results['PMI'] = {'raw_value': round(pmi_raw, 2), 'score': pmi_score}
        except (RuntimeError, ValueError) as e:
            print(f"Could not calculate PMI: {e}")
            results['PMI'] = {'raw_value': None, 'score': None}

        # --- Solvent Assessment Score ---
        try:
            score_val = self._get_solvent_score()
            # for this metric, the "raw_value" is the score itself
            results['Solvent_Score'] = {'raw_value': score_val, 'score': int(score_val)}
        except (ValueError, RuntimeError) as e:
            print(f"Could not determine Solvent Score: {e}")
            results['Solvent_Score'] = {'raw_value': None, 'score': None}

        # --- EQ ---
        try:
            eq_raw = self._calculate_eq()
            eq_score = self._score_eq(eq_raw)
            results['EQ'] = {'raw_value': round(eq_raw, 2), 'score': eq_score}
        except (RuntimeError, ValueError) as e:
            print(f"Could not calculate EQ: {e}")
            results['EQ'] = {'raw_value': None, 'score': None}

        # --- Hazard Score ---
        try:
            hazard_raw = self._calculate_hazard_score()
            hazard_score = self._score_hazard_score(hazard_raw)
            results['Hazard Score'] = {'raw_value': round(hazard_raw, 2), 'score': hazard_score}
        except (RuntimeError, ValueError) as e:
            print(f"Could not calculate Hazard Score: {e}")
            results['Hazard Score'] = {'raw_value': None, 'score': None}

        # --- Renewable Feedstock Score ---
        try:
            score_val = self._calculate_feedstock_score()
            feedstock_score = self._score_feedstock(score_val)
            # for this metric, the "raw_value" is the score itself too
            results['Feedstock_Score'] = {'raw_value': score_val, 'score': feedstock_score}
        except (ValueError, RuntimeError) as e:
            print(f"Could not determine Feedstock Score: {e}")
            results['Feedstock_Score'] = {'raw_value': None, 'score': None}

        # --- TON ---
        try:
            ton_raw = self._calculate_ton()
            results['TON'] = {'raw_value': round(ton_raw, 2), 'score': self._score_ton(ton_raw)}
        except (ValueError, RuntimeError) as e:
            print(f"Could not calculate TON: {e}")
            results['TON'] = {'raw_value': None, 'score': None}

        return results

    def calculate_time_series_metrics(self):
        """
        Calculates the evolution of metrics over the entire simulation time.
        This is a more complex, step-by-step calculation designed to generate the data needed for an animated, time-resolved dashboard.

        This method is fully vectorised for performance, avoiding explicit Python loops over time steps. It handles singularities at t=0
        by returning np.inf for metrics like PMI and EQ.

        Returns:
            dict: A dictionary where keys are metric names and values are arrays or lists representing the metric's value at each time step.
        """
        # data preparation...
        time_array = self.sim_results['time']
        species_data = self.sim_results['species_data']
        product_name = self.metadata['desired_product_name']
        all_species_names = list(species_data.keys())
        results = {'time': time_array}

        # pre-calculate mass arrays for all lspecies (vectorised)
        mass_arrays = {}
        for name in all_species_names:
            mw = self._molecular_weights.get(name)
            if mw is None:
                raise ValueError(f"Cannot calculate time series: Molecular weight for '{name}' is missing.")
            mass_arrays[name] = species_data[name] * mw

        # get product mass and moles arrays
        product_mass_array = mass_arrays.get(product_name, np.zeros_like(time_array))
        product_moles_array = species_data.get(product_name, np.zeros_like(time_array))
        
        # define a safe division threshold for mass/moles, a product is considered "formed" when its mass is non-negligible
        is_product_formed = product_mass_array > 1e-12

        # calculate dynamic PMI
        total_input_mass = sum(self.metadata.get('initial_masses', {}).values())
        results['PMI'] = np.divide(
            total_input_mass, product_mass_array,
            out=np.full_like(time_array, np.inf), # output is np.inf where...
            where=is_product_formed)              # ...the denominator is zero

        # calculate dynamic TON (cumulative)
        catalyst_moles = self.metadata.get('catalyst_moles', 0.0)
        if catalyst_moles > 0:
            results['TON'] = product_moles_array / catalyst_moles
        else:
            # if no catalyst, TON is not applicable and is zero throughout
            results['TON'] = np.zeros_like(time_array)

        # calculate dynamic EQ
        q_values = self.metadata.get('q_values', {})
        total_weighted_waste_array = np.zeros_like(time_array)
        for name, mass_array in mass_arrays.items():
            if name == product_name:
                continue
            q_value = q_values.get(name, 1.0) # default to benign
            total_weighted_waste_array += mass_array * q_value
        
        results['EQ'] = np.divide(
            total_weighted_waste_array, product_mass_array,
            out=np.full_like(time_array, np.inf), # output is np.inf where...
            where=is_product_formed)              # ...the denominator is zero
        
        # calculate dynamic hazard score
        presence_threshold = 1e-9 # moles
        
        # create an ordered list of species and stack their mole arrays into a 2D matrix
        ordered_species = list(species_data.keys())
        moles_matrix = np.vstack([species_data[s] for s in ordered_species]) # shape: (N_species, N_time)
        
        # create a corresponding column vector of Q-values
        q_values_vector = np.array([q_values.get(s, 1.0) for s in ordered_species]).reshape(-1, 1)
        
        # create a boolean mask where species are present above the threshold
        presence_mask = moles_matrix > presence_threshold
        
        # use the mask to create a matrix of Q-values, with 0 where species are absent
        # broadcasting automatically expands the Q-vector across all time steps
        masked_q_matrix = q_values_vector * presence_mask
        
        # find the maximum Q-value in each column (axis=0) to get the score at each time step
        results['Hazard_Score'] = np.max(masked_q_matrix, axis=0)

        return results
    
    def analyse_waste_composition(self):
        """
        Dissects the waste stream composition over time.
        Returns:
            dict: A dictionary containing the 'time' array and a nested dictionary with the percentage contribution of each waste component over time.
        """
        # data extraction and preparation
        time_array = self.sim_results['time']
        species_data = self.sim_results['species_data']
        product_name = self.metadata['desired_product_name']
        all_species_names = list(species_data.keys())

        # calculate mass array for all species
        mass_arrays = {}
        for name in all_species_names:
            mw = self._molecular_weights.get(name)
            if mw is None:
                raise ValueError(f"Cannot analyse waste: Molecular weight for species '{name}' is missing from metadata.")
            mass_arrays[name] = species_data[name] * mw

        # calculate total mass and waste mass arrays
        total_mass_array = sum(mass_arrays.values())
        product_mass_array = mass_arrays.get(product_name, np.zeros_like(time_array))

        total_waste_array = total_mass_array - product_mass_array

        # calculate percentage contribution of each waste component
        waste_composition_data = {}
        for species_name in all_species_names:
            if species_name == product_name:
                continue
            
            component_mass_array = mass_arrays[species_name]
            
            # use np.divide for safe division, returning 0 where the denominator is 0
            # handles the "perfect reaction" case where waste is zero
            percent_contribution_array = np.divide(
                component_mass_array, total_waste_array, 
                out=np.zeros_like(component_mass_array), where=total_waste_array > 1e-12
            ) * 100
            
            waste_composition_data[species_name] = percent_contribution_array

        return {'time': time_array, 'composition_percent': waste_composition_data}