# Introduction #

The model is based on the mixture theory framework. The model is described in more detail in `[1]`.

# Tissue model #

The tissue is supposed to consist of:
  * cells (which behave as liquid), volume fraction ![http://cord.googlecode.com/files/eq_phi.png](http://cord.googlecode.com/files/eq_phi.png)
  * extracellular liquid, volume fraction ![http://cord.googlecode.com/files/eq_phi_l.png](http://cord.googlecode.com/files/eq_phi_l.png)
  * extracellular matrix (rigid, non-remodelling), volume fraction ![http://cord.googlecode.com/files/eq_phi_ecm.png](http://cord.googlecode.com/files/eq_phi_ecm.png),
and the saturation condition to hold:

![http://cord.googlecode.com/files/eq_saturation.png](http://cord.googlecode.com/files/eq_saturation.png)

We use the following equation for the cells:

![http://cord.googlecode.com/files/eq_phi_eq.png](http://cord.googlecode.com/files/eq_phi_eq.png)

This equation was derived from mass and stress balance equations under the following assumptions:

  * liquid and cellular phases are incompressible and have the same density ![http://cord.googlecode.com/files/eq_rho.png](http://cord.googlecode.com/files/eq_rho.png)
  * the cells capture the extracellular liquid when growing, and upon death the cell components are returned into the extracellular liquid; ![http://cord.googlecode.com/files/eq_Gamma.png](http://cord.googlecode.com/files/eq_Gamma.png) is growth rate of the cells
  * relaxation of the mechanical stress is much faster than tissue growth
  * the tissue is isotropic; monotonic function ![http://cord.googlecode.com/files/eq_sigma.png](http://cord.googlecode.com/files/eq_sigma.png) describes cell--cell interaction being positive for compressed tissue; ![http://cord.googlecode.com/files/eq_sigma_0.png](http://cord.googlecode.com/files/eq_sigma_0.png) for cell packing density ![http://cord.googlecode.com/files/eq_phi0.png](http://cord.googlecode.com/files/eq_phi0.png)
  * adhesion between the cells and the matrix is described as viscous drag force, ![http://cord.googlecode.com/files/eq_Lambda.png](http://cord.googlecode.com/files/eq_Lambda.png) being a drag constant
  * cells--liquid and liquid--matrix interactions are negligible with respect to cells-matrix interaction

# Nutrient distribution #

_to_ _write_

# Growth and consumption #

_to_ _write_

# References #

  1. S. Astanin, A. Tosin, Mathematical model of tumour cord growth along the source of nutrient, _Mathematical_ _Modelling_ _of_ _Natural_ _Phenomena_, 2007, in press.









