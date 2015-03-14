This is an alternative implementation written for FreeFEM++. One programme is the numerical method for the basic model of tumour cord growth (one aerobic population), `cord.edp`. The other programme is implementing the bipopulation model (with additional glycolytic population), `cord_bi.edp`.

I wrote it mostly to learn about FreeFEM++ and to compare the results. FreeFEM++ is a nice tool for quick FEM programming.

Basic model:
  * [cord.edp](http://code.google.com/p/cord/source/browse/cord.edp?r=30ebc7d92463f0b425b26cc617a8ece42986bf66), the programme itslef (you need FreeFEM++ to run it)
  * [plotc](http://code.google.com/p/cord/source/browse/plotc?r=30ebc7d92463f0b425b26cc617a8ece42986bf66), a script for visualization with gnuplot (oxygen concentration)
  * [plotatp0](http://code.google.com/p/cord/source/browse/plotatp0?r=30ebc7d92463f0b425b26cc617a8ece42986bf66), a script for visualization with gnuplot (ATP deficit)
  * [Simulation video at SciVee](http://www.scivee.tv/node/5860)

Bipolation model:
  * [cord\_bi.edp](http://code.google.com/p/cord/source/browse/cord_bi.edp?r=30ebc7d92463f0b425b26cc617a8ece42986bf66), the programme (you need FreeFEM++ to run it)
  * [plotatp0\_bi](http://code.google.com/p/cord/source/browse/plotatp0_bi?r=30ebc7d92463f0b425b26cc617a8ece42986bf66), a script for visualization with gnuplot (ATP deficit)

### See also ###
  * About [FreeFEM++](http://www.freefem.org/ff++/)