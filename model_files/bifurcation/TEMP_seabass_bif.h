//*** TEMP_seabass.h   -- 

#ifndef BIFURCATION
#define BIFURCATION             1
#endif

#if (BIFURCATION == 1)
#define LENGTHCURVES            0
#else
#define LENGTHCURVES            20
#endif

#define POPULATION_NR           1
#define I_STATE_DIM             4
#define I_CONST_DIM             14
#define ENVIRON_DIM             5		
#define OUTPUT_VAR_NR           15	// 1st column is time so (1 + # outputs in .c)
#define PARAMETER_NR            64		
#define TIME_METHOD             RKCK
#define EVENT_NR                0       
#define DYNAMIC_COHORTS         0       
#define ALGEBRAIC_EQS           0
#define CHECK_EXTINCTION        2
