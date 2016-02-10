
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = /gpfs/home/rkn115/work/progs/PICC_Code
MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
PCC		 = mpicc

main.o: main.c

piccode: main.o
	${CLINKER} -o piccode main.o ${PETSC_KSP_LIB}

#----------------------------------------------------------------------------
runex1:
	-@${MPIEXEC} -n 1 ./ex1 -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex1_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex1_1.out ex1_1.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex1_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex1_1.tmp
runex1_2:
	-@${MPIEXEC} -n 1 ./ex1 -pc_type sor -pc_sor_symmetric -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always >\
	   ex1_2.tmp 2>&1;   \
	   if (${DIFF} output/ex1_2.out ex1_2.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex1_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex1_2.tmp
runex1_3:
	-@${MPIEXEC} -n 1 ./ex1 -pc_type eisenstat -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always >\
	   ex1_3.tmp 2>&1;   \
	   if (${DIFF} output/ex1_3.out ex1_3.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex1_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex1_3.tmp
NP = 1
M  = 4
N  = 5
MDOMAINS = 2
NDOMAINS = 1
OVERLAP=1
runex8:
	-@${MPIEXEC} -n ${NP} ./ex8 -m $M -n $N -user_set_subdomains -Mdomains ${MDOMAINS} -Ndomains ${NDOMAINS} -overlap ${OVERLAP} -print_error ${ARGS}

VALGRIND=
runex8g:
	-@${MPIEXEC} -n ${NP} ${VALGRIND} ./ex8g -M $M -N $N -user_set_subdomains -Mdomains ${MDOMAINS} -Ndomains ${NDOMAINS} -overlap ${OVERLAP} -print_error ${ARGS}


BREAKPOINT=
runex8g_debug:
	-@${MPIEXEC} -n ${NP} xterm -e gdb -ex 'set breakpoint pending on' -ex 'b ${BREAKPOINT}' -ex r -ex bt --args ./ex8g -M $M -N $N -user_set_subdomains -Mdomains ${MDOMAINS} -Ndomains ${NDOMAINS} -overlap ${OVERLAP} -print_error ${ARGS}

runex8_1:
	-@${MPIEXEC} -n 1 ./ex8 -print_error > ex8_1.tmp 2>&1; \
	   if (${DIFF} output/ex8_1.out ex8_1.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex8_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex8_1.tmp

runex8g_1:
	-@${MPIEXEC} -n 1 ./ex8g -M 7 -N 9 -user_set_subdomains -Mdomains 1 -Ndomains 3 -overlap 1 -print_error -pc_gasm_print_subdomains > ex8g_1.tmp 2>&1; \
	    if (${DIFF} output/ex8g_1.out ex8g_1.tmp) then true; \
	    else echo ${PWD} ; echo "Possible problem with with ex8g_1, diffs above \n========================================="; fi; \
	    ${RM} -f ex8g_1.tmp

runex8g_2:
	-@${MPIEXEC} -n 2 ./ex8g -M 7 -N 9 -user_set_subdomains -Mdomains 1 -Ndomains 3 -overlap 1 -print_error -pc_gasm_print_subdomains > ex8g_2.tmp 2>&1; \
	   if (${DIFF} output/ex8g_2.out ex8g_2.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex8g_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex8g_2.tmp

runex8g_3:
	-@${MPIEXEC} -n 3 ./ex8g -M 7 -N 9 -user_set_subdomains -Mdomains 1 -Ndomains 3 -overlap 1 -print_error -pc_gasm_print_subdomains > ex8g_3.tmp 2>&1; \
	   if (${DIFF} output/ex8g_3.out ex8g_3.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex8g_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex8g_3.tmp



runex10_basic:
	-@${MPIEXEC} -n 2 ./ex10 -f0 ${PETSC_DIR}/include/datafiles/matrices/spd-${PETSC_SCALAR_TYPE}-${PETSC_INDEX_SIZE}-${PETSC_SCALAR_SIZE} > ex10_1.tmp 2>&1; \
	   if (${DIFF} output/ex10_1.out ex10_1.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex10_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex10_1.tmp


include ${PETSC_DIR}/lib/petsc/conf/test
