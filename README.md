# MILP_Optimization-High-Temp-Steel
This repository contains a Mixed Integer Linear Programming (MILP) optimisation model for the identification of decarbonisation pathways for high temperature steel processes.
The model was set up as part of the master thesis project, conducted by Samuel Löffler as part of his master degree in Sustainable Energy Engineering at KTH Stockholm. The project was undertaken in collaboration with the Vienna University of Technology.

The optimisation model is set up in an object-oriented structure:
- The file 'BaseCase' serves as a main-file, implementing the constraints of the considered technologies and the objective function
- The classes in 'Optimization_model' define a type of considered technology
- Input_Data contains the necessary input files for running the optimization.

The input data for the electricity price was taken from the [AURES Technical Report on the Modelling of RES auctions](http://aures2project.eu/2022/10/17/modelling-of-res-auctions/). By running the script in 'Repr_period_selection', 4 one-week-periods of the full year were selected that represent the data set with least deviation of the load duration curve.

The data for solar irradiation at the location in Austria was provided by the [PVGIS interface](https://re.jrc.ec.europa.eu/pvg_tools/de/#HR). Again, 4 representative periods were selected here.

The demand data was provided by the company in the use case and can not be published, due to confidentiality obligations. The data was given as hourly heat demand values for a one-week period.

The full documentation on the master thesis project can be found under the following link, as soon as it is published on Diva.  
