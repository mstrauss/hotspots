This repository (find the latest version at https://github.com/mstrauss/hotspots) contains the software used in the article

> Strauss M, Klimek P, Sonneck G, Niederkrotenthaler T.  Suicides on the Austrian Railway Network: Hotspot Analysis and
Effect of Proximity to Psychiatric Institutions. [under review]

The code allows to identify hotspots/clusters (as defined in above paper) on a spatial network.  The original suicide data are not included, but instead we provide random surrogate data.  Still, to run the unmodified code, additional data from the cited sources need to be obtained.

Main Files
----------

The input of scripts in this folder may depend on outputs from other scripts.  The best order to run is:

1. ``script_master_sim``.  Fits the inhomogeneous Poisson model.
  * INPUT FILES (file formats are described below):
    * ``cases_mn.mat``: the railroad network, together with suicide locations (for privacy reasons, our supplied file has a RANDOM set of 1000 cases created on the railroad network)
  * OUTPUT FILES:
    * ``fit.mat``: The mean fitted model, after wiggling case locations.

2. ``script_cases_hotspots``.  Identifies the hotspots and evaluates the fitted model at the hotspots.
  * INPUT FILES:
    * ``fit.mat``.
  * OUTPUT FILES:
    * A CSV file containing the hotspots and associated information.

3. ``script_psych_distribution``.  Tests the null hypothesis (spatial vicinity to psychiatric institutions is not associated with railway suicides).
  * INPUT FILES:
    * ``hospitals.mat``.  Hospital locations.


License
-------

This software is subject to the included license, except where stated otherwise.

If you use this code in your publications, please cite above paper.
