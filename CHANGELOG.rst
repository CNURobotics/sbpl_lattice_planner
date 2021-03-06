^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package sbpl_lattice_planner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CHRISLab Customizations (2019-07-01)
* Use full name space to differentiate planner instances in ROS outputs
* Add Python scripts for generating and visualizing motion primitive files
* New primitives define action cost and use non-uniform angles

0.2.1 (2019-01-16)
------------------
* Reinit on map size, footprint and costmap changes
* Add warning when cost_scaling_factor is too large
  Also see `#33 <https://github.com/ros-planning/navigation_experimental/issues/33>`_.
* Ignore SBPL compile warning (`#31 <https://github.com/ros-planning/navigation_experimental/issues/31>`_)
* Fix example config for TF2 (`#30 <https://github.com/ros-planning/navigation_experimental/issues/30>`_)
* Update to tf2, add dependency
* Contributors: Jonathan Meyer, Martin Günther

0.2.0 (2018-09-03)
------------------
* Initial release into indigo, kinetic, lunar and melodic
* Contributors: Martin Günther, David V. Lu!!, Dave Hershberger, E. Gil Jones, Eitan Marder-Eppstein, Felix Widmaier, Johannes Meyer, Jon Binney, Vincent Rabaud, Austin Hendrix
