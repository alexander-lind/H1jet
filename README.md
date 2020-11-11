# H1jet
A fast and easy-to-use program to compute the differential transverse momentum distributions of one of the final-state particles produced in a
  <img src="https://render.githubusercontent.com/render/math?math=2\to 2"> process. Written in Fortran 95 and Python 3. 

The latest version can be obtained with: 
```
git clone https://github.com/alexander-lind/H1jet.git
```

See the User's Manual for more information, available at: 
https://arxiv.org/abs/2011.04694

## Dependencies 
H1jet requires the following external packages: 
 - HOPPET, for QCD DGLAP evolution of PDFs, available [here](https://github.com/gavinsalam/hoppet). 
 - LHAPDF, for PDF sets, available [here](https://lhapdf.hepforge.org).
 - CHAPLIN, for complex harmonic polylogarithms, available [here](https://chaplin.hepforge.org).

## Installation 
To compile H1jet: 
```
./configure
make
```
H1jet requires the CHAPLIN library. It may be necessary to explicitly state the path to the library files with: 
```
./configure LDFLAGS=-L/path/to/chaplin/lib
make
```
To compile with a custom user interface: 
```
./configure USERFILE=/path/to/custom/user_interface.f90 
make
```
See the [User's Manual](https://arxiv.org/abs/2011.04694) for more information on the implementation of a custom user interface. 

To install in a specific location: 
```
./configure --prefix=/path/to/installation
make
make install 
```
By default the prefix is set to the H1jet main directory. 

## Use 
To run H1jet from the main directory: 
```
./bin/h1jet 
```
H1jet will output a brief summary of the settings used along with the Born cross section <img src="https://render.githubusercontent.com/render/math?math=\sigma_0">, followed by the <img src="https://render.githubusercontent.com/render/math?math=\mathrm{d}\sigma/\mathrm{d}p_{T}"> and the integrated cross section <img src="https://render.githubusercontent.com/render/math?math=\sigma(p_{T})"> with a lower bound in <img src="https://render.githubusercontent.com/render/math?math=p_T"> for each <img src="https://render.githubusercontent.com/render/math?math=p_T"> bin. 

To get a complete list of options for H1jet: 
```
./bin/h1jet --help
```
See the [User's Manual](https://arxiv.org/abs/2011.04694) for more details on all options. 

The output from H1jet can be easily plotted with the provided helper script. 
Simply pipe the output from H1jet to the script: 
```
./bin/h1jet | python bin/PlotH1jet.py 
```
This will produce the following plot for default options: 
<br><img src="https://github.com/alexander-lind/H1jet/blob/master/tex/figures/H1jetresult.png?raw=true" alt="Example plot of default H1jet output" width="60%">
