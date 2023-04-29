# Running the program
## How to run the program ?

A **Makefile** is provided in the **build** folder. Every compilations and executions have to
be done inside this folder (otherwise the paths to the files will be broken). So, you
just need to type:

```
$ cd build/
$ make
```

in order to compile and run the program. If this command does not work, please try:

```
$ cd build/
$ cmake ../
$ make 	(or cmake --build .)
```

## How to use the program ?

The program is well commented and every possible operations are already written in the
**main.cpp** file. Each function (both for the 1D Ising and the 2D Ising) takes exactly
the same parameters in the same order, that is:

- `int iterations`: the numbers of Metropolis sweeps to be done (at fixed
temperature T). This parameter can vary to test the program but you can not have 
iterations < 1. 

- `std::mt19937 & gen`: pass the Mersenne Twister as argument. This is always the
same and you can not change this parameter.

- `int latticeSize`: the size of the lattice of spins. In the 2D case, the program
automatically assume that the real size of the lattice is latticeSize power 2 (e.g. if
you enter latticeSize = 4, then you get a lattice of size 4*4). This parameter can vary
to test the program but you can not have latticeSize < 1. 

- `double T_min`: the temperature where the generation of the datas starts. This
parameter can vary to test the program but you can not have T_min = 0 (to avoid
divergence when dividing by T) nor T_min >= T_max. 

- `double T_max`: the temperature where the generation of the datas ends. This
parameter can vary to test the program but you can not have T_max = 0 (to avoid
divergence when dividing by T) nor T_max <= T_min. 

- `double h`: the magnetic field to be applied on the lattice. This parameter can
vary to test the program (even with negative values of h). 

- `double J`: the coupling constant between nearest neighbor. This parameter can
vary to test the program (even with negative values of J).

- `double increment`: the distance between two points to be generated for the plot.
The total number of points produced is `int((T_max-T_min)/increment)`. This parameter can
vary to test the program but its value must be in the interval ]0,1[. 

- `int begin_average`: the point to start the calculation of the average. For
instance, if one see that after 50 000 iterations (`int iterations`) the system has reached
equilibrium, we can let int iterations = 100 000 and begin_average = 50 000. This way,
the average is computed with 50 000 Metropolis sweeps after the system is in equilibrium.
This parameter can vary to test the program but its value must verify begin_average > 0
and begin_average < iterations. 

- `const std::string & filename`: the path to the datas file where the output
should be written. 

## Where will the datas be written ?

By default, the program automatically write the datas into the appropriate subfolder of
the **datas** folder. If, for any reason, you want to change the paths where the datas
are written, you can easily do it in the **main.cpp** file. All the paths are specified at
the beginning in: `const std::string <filename> ("path/to/file")` variables. 

## What happens if a wrong argument is given to a function ?

If a wrong argument is given, that is if:

- latticeSize <= 0 || T_max <= 0 || T_min <= 0
- T_min >= T_max
- iterations < 1
- increment <= 0 || increment >= 1
- begin_average <= 0 || begin_average >= iterations
- fail to open the datas file(s)

an exception is raised and the program brutally stops. With this process, we prevent
any leak of memory because the latter is properly freed even when an unexpected
argument is passed. 

/!\ If you compile the program with g++, you may have a flag from the compiler looking
like: 

```
-warning: narrowing conversion from double to int at lines <number>
```

This conversion is WANTED and if you run the program it will work perfectly as expected.


# Results 

The program can generate the datas of the average magnetization, energy, magnetic susceptibility and specific heat both for the one dimensional and the two dimensional (square lattice) Ising model. In the 2D case, it can also produce the spins configuration at a given temperature as well as the Binder fourth order cumulant (used the determine the critical temperature). Few results (in 2D) are shown hereinafter. 

![spins configurations at the Curie temperature](https://github.com/vguillon/ising-model/blob/main/images/spinsConfigAtTc.png)

![average quantities of the 2D Ising model](https://github.com/vguillon/ising-model/blob/main/images/averages.png)

![Binder cumulant](https://github.com/vguillon/ising-model/blob/main/images/binderCumulant.png) 
