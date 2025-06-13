Quick Start
===========

To start, we will demonstrate how to calculate the metabolic fluxes at steady state using a `toy model <https://github.com/Chaowu88/freeflux/tree/main/models/toy>`__. This model includes typical reactions in the TCA cycle, with acetyl-CoA and aspartate as the initial substrates.

Loading a Toy Model
-------------------

A model can be defined in either a tab-separated values (.tsv) or Excel spreadsheet (.xlsx) file with the following required format:

.. list-table:: reactions.tsv
   :widths: 25 50 50 15
   :header-rows: 1

   * - #reaction_ID
     - reactant_IDs(atom)
     - product_IDs(atom)
     - reversibility
   * - v1
     - OAA(abcd)+AcCoA(ef)
     - Cit(dcbfea)
     - 0
   * - v2
     - Cit(abcdef)
     - AKG(abcde)+CO2(f)
     - 0
   * - v3
     - AKG(abcde)
     - Glu(abcde)
     - 0
   * - v4
     - AKG(abcde)
     - Suc(bcde)+CO2(a)
     - 0
   * - v5
     - Suc(abcd,dcba)
     - Fum(abcd,dcba)
     - 0
   * - v6
     - Fum(abcd,dcba)
     - OAA(abcd)
     - 1
   * - v7
     - Asp(abcd)
     - OAA(abcd)
     - 0
     
.. Note::
  1. The header must start with a '#' symbol.
  2. The letters in parenthesis indicates the C atom mapping in reactions.
  
Once the model is defined, it can be loaded using the following code:

.. code-block:: python
   
   MODEL_FILE = 'path/to/reactions.tsv'
   MEASURED_MDVS = 'path/to/measured_MDVs.tsv'
   MEASURED_FLUXES = 'path/to/measured_fluxes.tsv'
   
   from freeflux import Model
   
   model = Model('demo')
   model.read_from_file(MODEL_FILE)
   
Setting the Labeling Strategy
-----------------------------

To define the isotopic composition of input substrates (tracers), use the ``set_labeling_strategy`` method. You can specify different tracers for different experiments if you are working with a multi-experiment model.

.. code-block:: python
  
  fit = model.fitter('ss')

  # For an experiment identified as 'exp_glucose_tracer'
  fit.set_labeling_strategy(
      labeled_substrate='Glucose',
      labeling_pattern=['111111'],  # e.g., fully labeled U-13C Glucose
      percentage=[1.0],
      purity=[0.99],
      experiment_id='exp_glucose_tracer'
  )

  # For another experiment 'exp_glutamine_tracer' using a different tracer
  fit.set_labeling_strategy(
      labeled_substrate='Glutamine',
      labeling_pattern=['00110'], # Example: 3,4-13C Glutamine (assuming 5 carbons)
      percentage=[1.0],
      purity=[0.99],
      experiment_id='exp_glutamine_tracer'
  )

  # If experiment_id is omitted, the strategy applies to a default experiment (internally "exp0")
  # This is also suitable for single-experiment models.
  fit.set_labeling_strategy(
      labeled_substrate='AcCoA_from_default_exp', # Example substrate for default experiment
      labeling_pattern=['11'], # e.g., fully labeled AcCoA (2 carbons)
      percentage=[0.5],        # 50% of this tracer
      purity=[0.99]
      # experiment_id='exp0' # This is implicit if parameter is omitted
  )

.. Note::
   1. If the sum of the percentages of the specified isotopomer tracers for a given substrate in an experiment is less than 1, the remainder (1 - sum of percentages) will be considered as the unlabeled form (all atoms in natural abundance) for that substrate in that experiment.
   2. If you have multiple labeled substrates within the same experiment, or different labeling strategies for substrates across different experiments, you should call this method for each substrate and/or each experiment accordingly.
   
Adding Flux Bounds
------------------

We can set the search range for all fluxes or a specific flux using the ``set_flux_bounds`` method. Here, we set the fluxes to the range from -100 to 100:

.. code-block:: python

  fit.set_flux_bounds('all', bounds = [-100, 100])
  
Reading the Measuremnets
------------------------

We need to provide the labeling pattern or mass isotopomer distribution vector (MDV) of measurable metabolites or metabolite fragments for the fitting, as well as the measured exchange reactions, such as substrate uptake, product secretion, and specific growth rate. We can provide them using the following code:

.. code-block:: python
   
   fit.set_measured_MDV(
       'Glu_12345', 
       mean = [0.328,0.276,0.274,0.088,0.03,0.004], 
       sd = [0.01,0.01,0.01,0.01,0.01,0.01]
   )
   fit.set_measured_flux('v1', mean = 10, sd = 1)

.. Note::
   1. Try not to use excessively small values for the unknown measurement standard deviations. The objective function relies on summing weighted residuals based on the inverse variance of measurements. Using small standard deviation values can result in a large objective value and might hinder the optimizer to converge within a limited number of steps to reach a minimum.
   2. If you have a set of measured MDVs and fluxes, it is more convenient to read them from a .tsv or .xlsx file using the ``set_measured_MDVs_from_file`` and ``set_measured_fluxes_from_file`` methods.

.. admonition:: File Format for Measured MDVs (via ``set_measured_MDVs_from_file``)

   When using the ``set_measured_MDVs_from_file`` method, the input ``.tsv`` or ``.xlsx`` file should be structured as follows:

   *   **Comments**: Lines starting with a ``#`` symbol are treated as comments and will be skipped by the parser.
   *   **Header Row**: The first non-comment line is considered the header row, defining the column names.
   *   **Columns**:
       *   ``fragment_ID`` (or similar variation): Identifies the metabolite fragment (e.g., ``Glu_12345``). This column is used as the index for the data.
       *   ``mean``: Comma-separated values for the mean of the mass isotopomer distribution vector (e.g., ``0.6,0.3,0.1``).
       *   ``sd``: Comma-separated values for the standard deviation of the mass isotopomer distribution vector (e.g., ``0.01,0.01,0.01``).
       *   ``experiment_id`` (Optional): Specifies the identifier for the experiment. If this column is omitted, all MDV measurements in the file are assigned to a default experiment (e.g., "exp0").

   The parser attempts to normalize common variations in column names (e.g., ``FragmentID``, ``fragmentid``, ``Mean``, ``StdDev``, ``ExperimentID``, ``exp_id`` are typically understood).

   **Example ``measured_MDVs_single_exp.tsv`` (for a single/default experiment):**

   .. code-block:: tsv

      fragment_ID  mean            sd
      Glu_12345     0.6,0.3,0.1     0.01,0.01,0.01
      Ala_123       0.4,0.5,0.1     0.02,0.02,0.01

   *(This file could also contain comment lines starting with ``#`` before the header or between data rows, which would be ignored.)*

   **Example ``measured_MDVs_multi_exp.tsv`` (with multiple experiments):**

   .. code-block:: tsv

      fragment_ID  mean            sd              experiment_id
      Glu_12345     0.6,0.3,0.1     0.01,0.01,0.01  exp_A
      Ala_123       0.4,0.5,0.1     0.02,0.02,0.01  exp_A
      Glu_12345     0.5,0.4,0.1     0.01,0.01,0.01  exp_B
      Ser_123       0.2,0.6,0.2     0.03,0.02,0.03  exp_B
   
Solving the Fluxes
----------------

Now we can solve the flux distribution in the toy model using the following code:

.. code-block:: python
   
   fit.prepare()
   res = fit.solve()
   
The ``solve`` method returns a FitResults object. You can access the estimated net and total (including both forward and backward fluxes in reversible reactions) fluxes using the ``opt_net_fluxes`` and ``opt_total_fluxes`` attributes.

Working with the "with" Statement
------------------------

The ``fitter`` method returns a context manager, so you can also estimate the fluxes using the with statement, as shown in the following code:

.. code-block:: python
   
   with model.fitter('ss') as fit:
       fit.set_labeling_strategy(
           'AcCoA', 
           labeling_pattern = ['01', '11'], 
           percentage = [0.25, 0.25], 
           purity = [1, 1]
       )
       fit.set_flux_bounds('all', bounds = [-100, 100])
       fit.set_measured_MDV(
           'Glu_12345', 
           mean = [0.328,0.276,0.274,0.088,0.03,0.004], 
           sd = [0.01,0.01,0.01,0.01,0.01,0.01]
       )
       fit.set_measured_flux('v1', mean = 10, sd = 1)
       fit.prepare()
       res = fit.solve()
       
       
