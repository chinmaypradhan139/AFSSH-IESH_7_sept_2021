# AFSSH-IESH_7_sept_2021
Input the desired parameters for IESH in fort.23
input number of trajetories and number of cores in parallel_script.py
to start the simulation run myscript .i.e. ./myscript
On completion of calculation a file name running is created
To average over trajecorties copy avg2.py in running
To average the population calculated using Eq.(9) run python3 avg2.py
To calculate the population calculated using Eq. (18) replace fort.15 by fort.18.
each fort.15 and fort.18 will contain first coloumn as time and second as population.
output.txt will be the output file of averaged population.
To calculate marcus population change input of f(x) and g(x) function in marcus_plot.f90.
compile and run marcus_plot.f90 to generate output file fort.48 containing marcus population and time in atomic units.
