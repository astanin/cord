The code for numerical simulation of avascular tumour growth in the presence of blood vessels. A special case of formation of tumour cords is addressed. The mathematical model is being developed at Politecnico di Torino, Italy. Research was funded within the 5th Marie Curie Research Training Network [Modelling, Mathematical Methods and Computer Simulation of Tumour Growth and Therapy](http://calvino.polito.it/~mcrtn/).

The details on mathematical model published in
  * [S. Astanin, A. Tosin. Mathematical model of tumour cord growth along the source of nutrient // Mathematical Modelling of Natural Phenomena, 2(3), 153–177, 2007](http://dx.doi.org/10.1051/mmnp:2007007)

A new paper which describes bi-population model in detail:
  * [Astanin, S., Preziosi, L., Mathematical modelling of the Warburg effect in tumour cords. J. Theor. Biol. (2009), doi:10.1016/j.jtbi.2009.01.034](http://dx.doi.org/10.1016/j.jtbi.2009.01.034)

General description of the underlying mathematical framework may be found in:
  * [S. Astanin, L. Preziosi. Multiphase Models of Tumour Growth.](http://dx.doi.org/10.1007/978-0-8176-4713-1_9) In "Selected Topics on Cancer Modelling: Genesis — Evolution — Immune Competition — Therapy". 223–253, Birkäuser, 2008.

Some latest features of the model were presented in:
  * S. Astanin. Metabolic switch in bi-population model of tumour cord growth // CancerSim2008: Euroconference on Modelling and Simulation of Cancer Growth and Therapy, 19–21 May 2008, Turin, Italy
  * S. Astanin. Bi-population model of tumour growth and Warburg effect in tumour cords // SIMAI 2008: 9th Congress of the Italian Society for Applied and Industrial Mathematics in cooperation with SIAM, 15–19 September 2008, Rome, Italy.

## Versions ##

  * [current development version (Hg repository)](http://code.google.com/p/cord/source)
  * [cord-0.3.7 (old version)](http://cord.googlecode.com/files/cord-0.3.7.tar.gz)
  * [cord-0.3.6 (old version)](http://cord.googlecode.com/files/cord-0.3.6.tar.gz)
  * [cord-0.3.5 (old version)](http://cord.googlecode.com/files/cord-0.3.5.tar.gz)
  * [cord-0.3.4 (old version)](http://cord.googlecode.com/files/cord-0.3.4.tar.gz)

## See also ##
  * [README](http://code.google.com/p/cord/source/browse/README)
  * [Releases](http://code.google.com/p/cord/source/browse/NEWS)

## Alternative implementations ##

This numerical codes implement only some features of the model. The rely on third-party FVM/FEM tools which make numerical programming much easier. While these codes are slower, they are much shorter and easier to read and modify, and may serve as a good starting point to build upon/extend the model:

  * [FEM implementation in FreeFEM++](FreeFemPlusPlusImplementation.md)
  * [FVM implementation in python using FiPy library](http://code.google.com/p/cord/source/browse/tumour.py?r=d5a1bc29fce42037372ee1277cc8e2c6b6f64504)

## Some examples ##

In fact there are two different models supported by this code. The basic model assumes
that tumour growth is limited only by oxygen availability. The cells stop proliferating and die if there is not enough oxygen. There is only one kind of the tumour cells in this model, and only concentration of oxygen matters. This model can exlain cord expansion and fixed cord width.

A more detailed model assumes that under hypoxic conditions the tumour may spontaneously acquire ability to produce energy (ATP) without oxygen, relying on glycolysis. Thus emerge another “population” of tumour cells, with different metabolism. Also there are at least two chemicals to consider: oxygen and glucose. We refer to this model as binutrient (oxygen+glucose) or bipopulation (aerobic+anaerobic) model.

Please click on the images to watch simulation videos for bipopulation model (Images link to SciVee video site).

### Video 1: emergence of anaerobic cells ###

The first video shows evolution of the region where aerobic cells suffer from hypoxia (ATP deficit) as well as the limit where glycolytic cells start suffering too:

[![](http://cord.googlecode.com/files/zeroatp-clicktoplay.png)](http://www.scivee.tv/node/5153)

### Video 2: tissue movement ###

The second video shows velocity of tissue movement in the same simulation and volume fraction of glycolytic (anaerobic) cells:

[![](http://cord.googlecode.com/files/phi2andv_clicktoplay.png)](http://www.scivee.tv/node/5154)

### Video 3: growth in the vascular network ###

This video shows some (still experimental) modification to the model. There are inner vessels in the simulation domain which act as sources of oxygen. Glycolytic switch is disabled here. The corresponding source code available in _vascular_ branch. Click on the image to view the video at SciVee:

[![](http://cord.googlecode.com/files/vasc_o2_180400000.jpg)](http://www.scivee.tv/node/5889)