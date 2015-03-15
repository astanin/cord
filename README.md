Tumour cord growth model
========================

The model simulates tumour growth in the presence of blood
vessels. A special case of tumour cords is addressed.

## Publications

For general description of the modelling framework based on mixture theory see

 * Astanin, S., & Preziosi, L. (2008). Multiphase models of tumour growth. In Selected topics in cancer modeling (pp. 1-31). Birkhäuser Boston.
   [DOI: 10.1007/978-0-8176-4713-1_9](http://dx.doi.org/10.1007/978-0-8176-4713-1_9)

The basic tumour cord model was explained in

 * Astanin, S., & Tosin, A. (2007). Mathematical model of tumour cord growth along the source of nutrient. Mathematical Modelling of Natural Phenomena, 2(03), 153-177.
   [DOI: 10.1051/mmnp:2007007](http://dx.doi.org/10.1051/mmnp:2007007) (free access)

The full model, which includes transition of the tumour to the anaerobic metabolism, was published in

*  Astanin, S., & Preziosi, L. (2009). Mathematical modelling of the Warburg effect in tumour cords. Journal of theoretical biology, 258(4), 578-590.
  [DOI: 10.1016/j.jtbi.2009.01.034](http://dx.doi.org/10.1016/j.jtbi.2009.01.034)

## Usage

Running simulations is typically something like this:

    $ cord -v -t 100.0 -d 0.5 --dump-period 10.0 2>log ; tail log

Model and method parameters are defined in cord.ini file. By default
the file in the current directory is used. If it is not found, the
programme looks for a file installed with the $prefix/share/cord.
If neither file is found, the default values are used.

Command line options override settings given in the cord.ini file.

For more information about options see

    $ cord --help

Results are saved in HDF5.0 files. OpenDX (.dx) and Gnuplot (.gp)
files may be generated for visualization.

OpenDX (.dx) files may be generated from HDF files like this:

    $ cord_h5todx *.h5

Gnuplot (.gp) files may be generated from HDF files like this:

    $ cord_h5togp *.h5

To produce visualizations from Gnuplot (.gp) files using gnuplot:

    $ cord_gp2eps -gnuplot *.gp

To produce visualizations from Gnuplot (.gp) files using Tioga:

    $ cord_gp2pdf *.gp

To produce visualizations from Gnuplot (.gp) files using GRI:

    $ cord_gp2eps -gri *.gp

HDF5 utils can also convert HDF data to VTK format, suitable
for further processing with VTK or MayaVi. To extract packing density
profile (variable _phi_) from the data file run

    $ h5tovtk -d '/mesh/phi' file.h5

To extract oxygen distribution (variable _c_) run

    $ h5tovtk -d '/mesh/c' file.h5

And to extract level set function (variable _psi_) run

    $ h5tovtk -d '/mesh/psi' file.h5

## Installation

### Requirements

This programme depends on the following libraries:

  * libpopt (command line parsing, required)
  * Blitz++ (internal matrix representation, required)
  * BLAS (either ATLAS, or netlib's BLAS, or GotoBLAS, required)
  * GNU Scientific Library (ODE solver, numerical integration, required)
  * UMFPACK/libufsparse (generic SLE solver, optional, recommended)
  * LAPACK (tridiagonal solver, optional)
  * libhdf5 (data file format, optional)

The programme is written in C++, and was primarily developed with
GNU C++ compiler version 3.4.

Currently the source is distributed together with LSolver package
written by C. Badura (lsolver/ subdirectory) and iniParser by N.
Devillard (iniparser/ subdirectory).

Tested in:

  * Debian GNU/Linux (etch) i386
  * Windows XP SP2 + Cygwin

### How to build

The project uses GNU autoconf/automake build system. The typical
steps to bootstrap the built system are:

```
$ autoheader
$ aclocal
$ automake --add-missing
$ autoconf
```

You may build and install it running `./configure` script,
`make` and `make install`.  See INSTALL file for more details.

Make sure that you have all the dependencies installed (see above).

On my Debian system I usually run `configure' as

    $ ./configure --with-atlas=/usr/bin/sse2

as SSE2-optimized ATLAS libraries are installed in /usr/bin/sse2 on
my system. Alternatively, you may build the package against legacy
BLAS from netlib.org. I do it in my Cygwin installation this way:

    $ ./configure --with-fblas

as libfblas.a is in /usr/local/lib on my system.

If you have some of your libraries installed in non-standard
path, you can set LDFLAGS shell variable before running configure.
For example, if I have popt library in $HOME/usr/lib64, I may
run configure like this:

    $ LDFLAGS="-L$HOME/usr/lib64" ./configure

Another example, if you installed blitz++ library to $HOME/usr,
its headers will be in $HOME/usr/include. The compiler will be
able to find them if you set CPPFLAGS variable, e.g.:

    $ CPPFLAGS="-I$HOME/usr/include" \
         LDFLAGS="-L$HOME/usr/lib" \
         ./configure

### How to build BLAS

If BLAS is not provided in your system (like in Cygwin), you may
build it yourself from source. ATLAS is a better option, but I
BLAS is easier to start with.

 1. Get BLAS (blas.tgz) from http://www.netlib.org/blas/
 2. Unpack the BLAS source package and go to inside:

    ```
    $ tar zxf blas.tgz
    $ cd BLAS
    ```

 3. Compile all the Fortran files:

    ```
    $ g77 -c -O2 *.f
    ```

 4. Assemble all the object files produced into library libblas.a.
    You may later link to this library statically.

    ```
    $ ar rvs libblas.a *.o
    ```

 5. Install the library into the system, e.g. copy it to /usr/local/lib

    ```
    $ cp libblas.a /usr/local/lib/
    ```

You can obtain better performance with ATLAS. You may find instructions
for Cygwin installation in http://www.ifp.uiuc.edu/~nakazato/tips/cygwin_atlas.html.
Most Linux users are likely to be happy with the ready-to-use ATLAS packages
of their distro.

## Credits and license

The model has been developed at Politecnico di Torino, Italy. 
Research and development was funded within
[the 5th Marie Curie Research Training Network](http://calvino.polito.it/~mcrtn/).

The code is distributed under the terms of the GNU General Public
License version 2 or later. See COPYING for details. iniParser and
LSolver distributed with the code are not covered by this license.


## Contacts:

 * Sergey Astanin `<sergey dot astanin at polito dot it>`
 * Luigi Preziosi `<luigi dot preziosi at polito dot it>`
 * Andrea Tosin   `<andrea dot tosin at polito dot it>`
