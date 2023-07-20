# WRFsnowdrift

modified code from https://github.com/ChristinaSchmid/snowdrift for snow drift simulations in WRF, changes see READDME.WRFVersions, currently setup for ideal setting -> no variable snicexy, has to be included for real case in module_first_rk_step_part1.F, module_snowdriver.F, module_snower.F, module_snowset.F

Usage:
- copy module_snowdriver.F, module_snower.F, module_snowset.F, module_snowvel.F, module_snowsubl.F into /phys
- add these lines in /phys/Makefile :
   module_snowdriver.o\
   module_snower.o\
   module_snowset.o\
   module_snowsubl.o\
   module_snowvel.o
- search for "cschmid" in Registry.EM_COMMON and add the additional lines
- search for "cschmid" in dyn_em/solve_em.F and module_first_rk_step_part1.F and add the additional lines
