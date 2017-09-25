# Installation
## Prerequisites
#### macOS
- [libgd](https://libgd.github.io/) is required for writing images to disk and can be installed using [brew](https://brew.sh/) or [macports](https://www.macports.org/), e.g. with

```
sudo port install gd2
```

- OpenCL drivers are included in macOS versions later than 10.6 Snow Leopard but we recommend 10.11 El Capitan. Unfortunately, Apples OpenCL drivers tend to be [buggy](http://preta3d.com/os-x-users-unite/) and the only way to change them is by upgrading to a different version of macOS. We cannot say much about the Sierra drivers so far but have made some bad experience with Yosemite.

#### Linux
- [libgd](https://libgd.github.io/) is in the Ubuntu package repository

```
sudo apt-get install libgd-dev
```

- Depending on your hardware you have to install OpenCL drivers according to your vendor. Coming Soon.

## Downloading and Compiling
Either download through the browser or

```
git clone https://github.com/pSpitzner/AtmoCL.git
```

Change into the ```AtmoCL``` folder and

```
cd AtmoCL
make
```

to get the ```executable```. Possibly you have to adjust the path to linked libraries and includes, depending on where you installed OpenCL and libgd. Note that this is platform dependant, edit the right ones.

## Running the Program
If you have a look into the ```run``` file, you will see that this just calls the ```executable``` and does some cleanup work. Start the program using

```
./run
```

and you will see a list of available devices [(adapted from here)](https://github.com/miguelao/clDeviceQuery) and the according driver version, such as below:

```
Available CL Devices (choose via '-dev i')

0 Intel(R) Core(TM) i7-2600K CPU @ 3.40GHz
  Apple - OpenCL 1.2 (Nov  1 2016 21:34:57)

1 AMD Radeon HD Tahiti XT Prototype Compute Engine
  Apple - OpenCL 1.2 (Nov  1 2016 21:34:57)
```

In this example, to choose the GPU one would type

```
./run -dev 1
```

If you run for the first time, a snapshot of the initial profile will be created after equilibrating (and every full hour) and saved in a snapshots folder. Such profiles can be imported using the ```-b``` argument.

```
./run -dev 1 -b ./snapshots/equil/
```

There is also a ```-profiling``` option, which forces the host to wait for the device to finish every step before enqueuing the next. It also creates additional logs for performance measuring and may help debugging. Try this if you are running into problems.

Any generated output is written into the ```output``` folder, which you can change into a symlink to save elsewhere. According subdirectories are created for images ```img```, ```timeseries```, ```verticalprofiles``` and ```logs```.

## Changing Settings

So far, settings are hard coded into the program, so you have to recompile after doing any changes.

```
make clean; make
```

Most settings are located at the beginning of the ```classes/asystem.cpp``` class. They are grouped into a global parameter struct we call ```par```, that is also available in kernels and other classes.

#### System Size
Before increasing the volume to get an actual 3D simulation box, make sure your computer has the power to handle it. Especially on Latops running macOS, your whole OS interface may become unresponsive up to the point where you are unable to kill the process. To be safe, run on a remote box over ssh.

You can change the amount of cells in each dimension

```
par.sx  = 128;
//par.sy  = 1;
par.sy = 128;
par.sz  = 64;
```

and the cell size

```
par.dx  = 50.0;
par.dy  = 50.0;
par.dz  = 20.0;
```

to get the desired system volume.

#### Time Step
```
par.dT = 1.0;
```
is the large time step in seconds from which intermediate step sizes ```par.dtis[s]``` and step amounts ```par.nsi[s]``` are calculated for the different stages.

#### Export Interval
In the ```classes/clexporter.cpp``` class you will find the three parameters setting the interval of exporting images, timeseries and vertical profiles
```
img_every = 60.0;   // write image files every 60 seconds
ts_every  = 60.0;   // make an entry into timeseries
vp_every  = 120.0;  // write a vertical profile every 120s
```

#### Branches and Test Cases
Presently, the repo contains the following branches: ```cell```, ```isdac``` and ```wrf```, the former are test cases and the latter is our first draft to import wrf files.
So far, the [cell](https://doi.org/10.1175/1520-0493(1982)110<0504:TDONSC>2.0.CO;2) works as expected.
[Isdac](http://dx.doi.org/10.1002/2013MS000282) produces qualitatively matching results but needs further work. (```N_i``` isnt fixed to desired value yet, simulation includes ice microphysics at all times).

To run a test case, checkout the branch, make and run, i.e.
```
git checkout cell
make clean; make
./run
```


