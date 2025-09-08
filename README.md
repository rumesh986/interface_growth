## Preparing the environment:

- Download and set up [Oomph-Lib](https://github.com/oomph-lib/oomph-lib), the instructions below are adapted from [here](https://oomph-lib.github.io/oomph-lib/doc/the_distribution/html/index.html)
- Clone the repo
	```
	git clone git@github.com:oomph-lib/oomph-lib.git
	```
	From here onwards, `$OOMPH_DIR` will refer to the directory the repo has been cloned into.
- Configure and Run the installation
	```
	cd $OOMPH_DIR
	```
	When your run the actual build command, you will be asked to select a configuration file in `$OOMPH_DIR/config/configure_options/`. You can make your own or copy an existing one. The following options were used for the code in this repo.
	```
	--enable-suppress-doc
	--enable-suppress-pdf-doc
	CXXFLAGS="-O3 -Wall"
	CFLAGS="-O3 -Wall"
	FFLAGS="-fallow-argument-mismatch -O3 -Wall"
	FFLAGS_NO_OPT="-O0"
	```
	
	Note that not all options will be needed. The first two options save some time when compiling the library.

	Once `autogen` has been configured, you can actually build the library with
	
	```
	./autogen.sh
	```

	> [!important]
	> On SCRTP systems you may need to load the appropriate modules first. These can be found in this repo in the file `erf2d/scrtp_modules`. You will need to source this file with `. ./erf2d/scrtp_modules`. As you may not have cloned the repo yet, you may have to download the file seperately or just run the commands within the script yourself.

	You will be asked some questions before the build actually begins. You can use the default options for these (or skip the self-tests to save time).

- Once the library has been built, you can go to `$OOMPH_DIR/user_drivers/` and clone this repo.

```
cd $OOMPH_DIR/user_drivers/
git clone git@github.com:rumesh986/interface_growth.git
```

From here onwards, we will refer to the location of this repo as `$IFACE_DIR`

- go to `$IFACE_DIR/erf2d` to find the latest code for simulating a moving interface.

## Running the code

- The different `erf2d*` files refer to different problems that had been solved, they mostly very with different interface mechanisms and boundary conditions.

- The actual code used in the report is based on `erf2d8` which tracks evolves the interface position based on the heat flux at the interface.
	- The older codes can also be run but their numerical results may not be valid as a lot of bugs had been identified and fixed in `erf2d8` but were not back-ported to the older versions.

- Helper scripts have been written in python to properly build and run the relevant code, as well perform post-processing.
	- To simply run the code, call `./run erf2d8`. This will create a folder in `$IFACE_DIR/erf2d/runs` with the current time. Any graphs created in post-processing will appear in this folder and any sub-folders.
	- To see the available options, run `./run --help` or refer to table below.

| parameter  	| expected input		| default value 			| meaning |
| ------------- | ---------------------	| ------------------------- | ------- |
| --nx --nxs 	| int/array of int 		| 10 						| Number of elements in each phase of simulation domain |
| --dt --dts 	| float/array of float 	| 0.01 						| Timesteps to use for simulations |
| --x --xs	 	| int/array of int 		| 2							| Number of nodes in one dimension of node (for example, this has to be 2 for linear shape functions) |
| --t --ts	 	| int/array of int 		| 1 						| Order of timestepper to use, valid values are 1,2,4 |
| --tsteps	 	| int					| 100 						| Number of steps to run simulation for	|
| --type	 	| 'a', 'x', 't', 's' 	| 'a' 						| Type of analysis to perform, <br> 'a' - standard run and post-process results', <br> 'x' - error analysis for spatial discretization, <br> 't' - error analysis for temporal discretization, <br> 's' - sensitivity analysis |
| --tshift		| float 				| 0.0 						| Time to start simulation from |
| --k1			| float 				| 1.0 						| Thermal Conductivity in phase 1 |
| --k2			| float 				| 0.5 						| Thermal Conductivity in phase 2 |
| --rho1		| float 				| 1.0 						| Density of material in phase 1 |
| --rho2		| float 				| 0.5 						| Density of material in phase 2 |
| --cp1			| float 				| 1.0 						| Specific heat capacity in phase 1 |
| --cp2			| float 				| 0.5 						| Specific heat capacity in phase 2 |
| --L 			| float					| 1.0 						| Latent heat of fusion
| --eps			| float					| 1e-6						| Factor to vary simulation parameters for sensitivity analysis |
| --debug 		|						|							| Run in debug mode, output will be saved to `$IFACE_DIR/runs/debug` folder |
| --nthreads	| int					| Number of available CPUs	| Number of CPUs to use for computation, The actual number of cpus will be limited by the number of tasks that need to be run. |
| --style		| file path				| 'pltstyle.mplstyle'		| Matplotlib style to use for produced images, for report style images, point to `report.mplystyle` instead |

The following commands were used when creating figures from the report
| figure | command |
| ------ | ------- |
| 3, 4	 | `./run erf2d8 --nx 100 --dt 0.0001 --tshift 0.001 --x 2 --t 2 --style report.mplstyle` |
| 5		 | `./run erf2d8 --nx 10  --dt 0.0001 --tshift 0.001 --x 2 --t 2 --style report.mplstyle` |
| 6a	 | `./run erf2d8 --nx 3 4 5 10 20 30 40 50 100 200 400 --dt 0.0001 --xs 2 3 4 --ts 2 --tshift 0.001 --type x --style report.mplstyle` |
| 6b	 | `./run erf2d8 --nx 100 --dt 0.1 0.01 0.05 0.08 0.003 0.002 0.001 0.0002 0.0005 0.0008 0.0001 --xs 2 --ts 1 2 4 --type t --style report.mplstyle` |
| 7		 | `./run erf2d8 --nx 100 --dt 0.0001 --tshift 0.001 --x 2 --t 4 --type s --style report.mplstyle` |