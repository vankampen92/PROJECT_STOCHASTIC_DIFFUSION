# README #

This project runs on linux machines out of the box. To run in other platforms makefiles should be modified accordingly. See makefile to see the suitable modifications for Mac users (tested on a Maverick version).

The follwing dependencies are required:

## GNU Scientific Library (GSL)
This is a numerical library for C and C++ programmers. It is free software under the GNU General Public License.
## CPGPLOT Graphix Library.
This is a set of functions written in C relying on cpgplot primitives from pgplot, and plplot. As a consequence, the CPGPLOT Graphix library, in turns, depends on:
### [pgplot](/http://www.astro.caltech.edu/~tjp/pgplot/)
### [plplot](http://plplot.sourceforge.net/)
You can git clone the CPGPLOT library from my repository.

### What is this repository for? ###

* This repository sets up a number of dynamic models of consumer-resource interactions distributed across a network metatapopulation structure. Individual movement is implemented as a random walk betweeen connected patches. Both ODEs and Gillespie simulations are implemented.

* Version: 0.0.0.999
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

First you should install the libraries to meet the dependencies mentioned above. Look for GSL, pgplot and libplplot12 in your usual package handler.
Notice that the linking command from most makefiles contains, at least, the following libraries:

* -lgsl -lgslcblas
* -lWL -lpng -lplplotd -lpgplot -lcpgplot
* -lda_cpgplot_XY_GRID_FILE -lda_cpgplot_BASIC

The first two ones are basic GSL libraries. The following 5 are required to use primitive plotting functions from the graphic libraries cpplot and plplot. The final two are mandatory when using higher-level plotting functions from the CPGPLOT Graphix library. All of them are usually required to produce a graphical output. However, the control variable (see any makefile) 'CPG" can also be set up to 'NON_CPG_REPRESENTATION' and, then, through conditional compilation, the same program is built to just run the numerical computations without graphics. In all cases, the output may be saved in files.

When you git clone the repository on your machine, you should do it from your home directory. As a result, the directory '~/PROJECT_STOCHASTIC_DIFFUSION' will be created on your machine.

If graphic libraries have been correctly installed, this should be enough to make all makefiles work out of the box. Remember though you require to have also git cloned my CPGPLOT repository on your machine. To be clear, you should end up with two directories:

* -$ ./CPGPLOT
* -$ ./PROJECT_STOCHASTIC_DIFFUSION

both in you home directory.

* Summary of set up:
	+ #### 1. Install GSL library
	+ #### 2. Install plplot library
	+ #### 3. Install pgplot library
	+ #### 4. git clone https://github.com/vankampen92/CPGPLOT
	+ #### 5. git clone https://github.com/vankampen92/PROJECT_STOCHASTIC_DIFFUSION.git
	+ #### 6. Tests:
	In order to test if pgplot, plplot and CPGPLOT are correctly installed in your machine, you can expand the tar file PROJECT_CPGPLOT_EXAMPLES.tar, which is in the project root directory on your home directory. Then you will get the directory ~/PROJECT_CPGPLOT_EXAMPLES. In that directory, there is a simple example of how to use the CPGPLOT library. You build it by typing:

		+ ~/PROJECT/CPGPLOT_EXAMPLES/make

		and you will get the exectutable file PLOT. You may run the example with some command arguments (see main.c). You should get a graph with four different subplots. You may also type:

		+ ~/PROJECT/CPGPLOT_EXAMPLES/PLOT -h

		and see other available command line arguments. You may also type:

		+ ~/PROJECT/CPGPLOT_EXAMPLES/PLOT -G29 ?

	and see the different available graphic formats in which plots can be saved. Notice that sometimes the value for these input arguments is overriden by the internal program code. When this happens, it is for a good reason. Please check the code to understand why and make moodgodfications at your own risk. Be creative.   

	+ #### 7. Examples:
	See, for instance, ./MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_STOCHASTICS/main.c, and follow the directions to compile and run the code:

		+ ~$ make MODEL=DIFFUSION

		+ ~$ ./DIFFUSION -y0 0 -y2 1 -HS 1 -HM 10000 -HX 100 -HY 100 -Hu 0.5 -n 1 -v0 10105 -G0 1 -G1 1 -tn 100 -t0 0.0 -t1 30.0 -t4 0 -tR 2 -xn 0 -xN 1000 -HN 1000 -G2 1 -G3 0.0 -G4 30.0 -G5 1 -G6 0.0 -G7 1100.0

	The code depends on some auxiliary libraries in ./Library  and ./Definition_Error_Model subdirectories. You may notice that you need to generate these libraries before, and then execute the command 'make MODEL=DIFFUSION'. In principle, a recursive makefile does this job for you. However, if gcc does not find these libraries, they may have been accidentally deleted and you should build them back up again from sources. Also, the code is linked against R libraries.  You may remove these R links or install R in your system. I recommend this 2n option. This will allow you to create shared libraries that, then, can be called as standard R funcions from RStudio, for example.

	The call on the 2nd line above generates a bunch of stochatic realizations (-tR 10) and presents a single output variable (-n 1), this is, the temporal evolution of the central cell of the 100 times 100 grid. Local populations thrive in the 10000 cells (-HM 10000), organized on a 100 times 100 squared grid (-HX 100 -HY 100). The type of network in controled by the -y2 imput argument value. In this case, grid connections are Von Neumann with periodic boundary conditions (-y2 1).

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###
code
* Drop an email to David Alonso (<dalonso@ceab.csic.es>)
