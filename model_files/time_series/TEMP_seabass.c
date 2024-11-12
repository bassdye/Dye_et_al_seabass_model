/***
   NAME
    TEMP_working.c

   NOTES
        Definition file for single cohort and population level model (van Leeuwen et al. 2013),
          adapted with temperature-dependent functions from Ohlberger et al. (2011).
        
        Used to model: cod, sprat, herring, rockling, mullet, seabass 

        Uses power function for simplified ingestion rate (i.e. simple_ingest)

        Differences from Ohlberger et al. (2011) include:
        - No adjustment factors in ingestion and maintenance scaling
          See Appendix - Temperature Dependence for more info
            - i.e. we are scaling from 0 - 1
        - Using parameter values (40 and 400) in temperature dependent equations 
          unlike in Ohlberger et al. (2011) which used 20 and 200 for ingestion
            See Table 2. a (bottom of table) for more info
        - Rm calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)
          - We are using standarized weight (1 + QJ)*X instead; no maintenance costs for gonads (Kooijman, 2010; Soudijn & de Roos 2017)
            - Define either calculation with RM_CALCULATION directive

    The specific definition file for a size structured model including:

      - Constant resource level (Zooplankton resource ~ FEEDING_LEVEL) with FEEDING_LEVEL parameter value (0-1); 1 = max growth. 
      - Allows for constant or changing, seasonal temperature (Pethybridge et al. 2013)
      - Allows for constant or changing, seasonal resource (Pethybridge et al. 2013)
      - Allows for incorporation of marine heat waves (using modified seasonal temperature function)
      - Allows for two resources and an ontogenetic diet shift 
      - Allows for three resources and two ontogenetic diet shifts
      - Allows for constant fish resource (i.e. set to carrying capacity; CONSTANT_FISH definition)
      - Allows for mature fish to experience a minimum winter temperature (T_ADULT_MIN) during the seasonally varying temperature

    Terminology:

        X:      Irreversible mass or bone
        Y:      Reversible mass or fat
        M:      Mass or standard weight, (1+QJ)*X
        W:      Weight or total weight, (X+Y)

    Last revisions: OG AvL - 5 May 2022; TEMP BD - March 27 2024
***/

// General directives
#define IMMEDIATE_HATCHING      1               // 1: Cohort is filled at first day of spawngroup
#define EGGM_EQ_SDM             0               // 0: Following Ohlberger et al. 2011, no egg mortality; 1: egg/yolk-sac mortality follows size-dep. mortality

// #define END_INVASION                 0               // 1: ForcedRunEnd if R0 > 1.5 or R0 < 0.5
// #define R0_OUTPUT               2               // 0: None; 1: output of ReproOutputs[][]; 2: In addition, output of collapsed matrix

#define RM_CALCULATION          1               // 0: Calculate maintenace with total weight (Ohlberger et al. 2011); 1: Calculate with standarized weight (Kooijman, 2010)
#define TEMPERATURE             1               // 0: Constant temperature; 1: Seasonally varying temperature using temp_sin function 
#define HEAT_WAVE               0               // 0: No heat wave; 1: Include heat wave(s) using heat_wave function 
#define TEMP_ADULT_MIN          0               // 0: 1: Mature fish subjected to minimum temperature (T_ADULT_MIN) during seasonally varying temperature using temp_sin function
#define RESOURCE                0               // 0: Constant resource; 1: Seasonally varying resource using resource_sin function 
#define TWO_RESOURCES           0               // 0: Single resource; 1: Include two resources and ontogenetic diet shift
#define THREE_RESOURCES         1               // 0: One (two) resource(s); 1: Include three resources and ontogenetic diet shifts
#define CONSTANT_FISH           1               // 0: Semichemostat growth of the fish resource; 1: Constant resource level equal to carrying capacity, K_FISH

// Species directives      
#define SPECIES_SINGLE          0               // 1: Single cohort, no feedback on environment
#define SPECIES_COHORT          1               // 1: Multiple cohorts, feedback on environment

// #define SPECIES_INVADE          0 

/*#if (BIFURCATION == 0)
#undef  SPECIES_INVADE
#define SPECIES_INVADE          0
#endif*/

#include  "escbox.h"

/*
 *======================================================================
 *
 *              LABELING ENVIRONMENT AND I-STATE VARIABLES
 *
 *======================================================================
 */

#define time           (env[0])

#define plankton       (env[1]) 
#define shrimp         (env[2])
#define fish           (env[3])
#define total_eggs     (env[4])

#define SPECIES           0

#define age             (i_state( 0)) // Second value of second row of .isf. First value is # individuals (refer to EBTmanual.pdf)
#define bone            (i_state( 1))
#define fat             (i_state( 2))
#define gonads          (i_state( 3))

// Constants not set in updateIDcards()
#define IDnbegin        (i_const( 0))
#define IDlifefec       (i_const( 1))
#define IDspawning      (i_const( 2))

// Constant for temporary storage
#define IDscratch       (i_const( 3))

// Size measures set in updateIDcards()
#define IDmass          (i_const( 4))
#define IDweight        (i_const( 5))
#define IDlength        (i_const( 6))
#define IDfatratio      (i_const( 7))

// Energetic rates set in updateIDcards()
#define IDingest        (i_const( 8))
#define IDmaint         (i_const( 9))
#define IDnet_energy    (i_const(10))
#define IDkappa         (i_const(11))

// Mortality & reproduction rates set in updateIDcards()
#define IDfecundity     (i_const(12))
#define IDmort          (i_const(13)) // Initial mortality, change to 0.0 if you want no initial mortality

#define MILLI           0.001
#define MICRO           1.0E-6
#define DIV_1_400       0.0025        // Used in temperature equation

#define MINCONDITION    1.0E-6
#define FOODTINY        1.0E-10

#define MAXCOHORTS      10000         // Used as length for arrays in vdPowx() --> should be int however...
/*#define INVADEDENSITY   1000
#define INVADEPERIOD    30            // Period (in years) with cohort introductions for measuring invasion*/
#define MINSURVIVAL     1.0E-9

/*
 *======================================================================
 *
 *          DEFINING AND LABELLING CONSTANTS AND PARAMETERS
 *
 *======================================================================
 */

// Parameters 
#define YEAR            parameter[ 0]

// Water temperature parameters  
#define TEMP            parameter[ 1]   // User defined constant environmental temperature 
#define T_MEAN          parameter[ 2]   // Temperature sin function mean value
#define T_A             parameter[ 3]   // Temperature sin function amplitude
#define T_OMEGA         parameter[ 4]   // Temperature sin function phase shift

// Species parameters
#define SPAWNSET        parameter[ 5]
#define SPAWNMAX        parameter[ 6]
#define SPAWNPERIOD     parameter[ 7]
#define SPAWNGROUP      parameter[ 8]

#define SPAWNSTART      (SPAWNMAX-0.5*SPAWNPERIOD)
#define SPAWNEND        (SPAWNMAX+0.5*SPAWNPERIOD)

#define EGGPERIOD       parameter[ 9]   // Duration of egg period
#define AGEFEEDING      parameter[10]   // Age at first feeding
#define BIRTHWEIGHT     parameter[11]   // Total weigth at birth
#define MATURELEN       parameter[12]   // Maturation length

#define QS              parameter[13]   // Starvation condition threshold
#define QJ              parameter[14]   // Juvenile condition target
#define QA              parameter[15]   // Adult condition target
#define QSPAWN          parameter[16]   // Threshold condition for spawning

#define MUB             parameter[17]   // Background mortality
#define MUSDC           parameter[18]   // Size-dependent mortality constant
#define MUSDMASS        parameter[19]   // Size-dependent mortality characteristic size
#define MUS             parameter[20]   // Starvation mortality constant
#define EGGMORT         parameter[21]   // Egg mortality
#define YOLKMORT        parameter[22]   // Mortality of yolk-sac larvae

#define CONVEFF         parameter[23]   // Conversion efficiency
#define REPROEFF        parameter[24]   // Reproductive conversion efficiency

#define LWC             parameter[25]   // Weight-length scaling constant
#define LWE             parameter[26]   // Weight-length scaling exponent

#define MAINTC          parameter[27]   // Maintenance scaling constant
#define MAINTE          parameter[28]   // Maintenance scaling exponent

#define GA_M            parameter[29]   // Allometric scalar of Qm; metabolism
#define THETA_M         parameter[30]   // Allometric exponent of Qm
#define YMAX_M          parameter[31]   // Allometric scalar of Tmax
#define VMAX_M          parameter[32]   // Allometric exponent of Tmax
#define YOPT_M          parameter[33]   // Allometric scalar of Topt
#define VOPT_M          parameter[34]   // Allometric exponent of Topt

#define INGEST_C        parameter[35]   // Simplified ingestion scaling constant
#define INGEST_E        parameter[36]   // Simplified ingestion scaling exponent

#define GA_A            parameter[37]   // Allometric scalar of Qa; ingestion
#define THETA_A         parameter[38]   // Allometric exponent of Qa
#define YMAX_A          parameter[39]   // Allometric scalar of Tmax
#define VMAX_A          parameter[40]   // Allometric exponent of Tmax
#define YOPT_A          parameter[41]   // Allometric scalar of Topt
#define VOPT_A          parameter[42]   // Allometric exponent of Topt

// Resource and diet shift parameters 
#define R_PLANKTON      parameter[43]   // Growth rate zooplankon
#define K_PLANKTON      parameter[44]   // Carrying capacity zooplankton 

#define R_SHRIMP        parameter[45]   // Growth rate shrimp
#define K_SHRIMP        parameter[46]   // Carrying capacity shrimp 

#define R_FISH          parameter[47]   // Growth rate fish
#define K_FISH          parameter[48]   // Carrying capacity fish 

#define BEGIN_SHRIMP    parameter[49]   // Length to begin diet shift from zooplankton to shrimp 
#define HALF_SHRIMP     parameter[50]   // Length at 50% benthivory of shrimp 
#define END_SHRIMP      parameter[51]   // Lenth at full diet shift to shrimp (i.e. length feeding only on shrimp) 
#define BEGIN_FISH      parameter[52]   // Length to begin diet shift to piscivory
#define HALF_FISH       parameter[53]   // Length that piscivory is 50% (and shrimp is 50%)

#define FEEDING_LEVEL   parameter[54]   // Constant resource level (for single cohort model)

// Parameters below are not currently used but don't want to remove from script in case we want to use later on 
#define R_MEAN          parameter[55]   // Resource sin function mean value
#define R_A             parameter[56]   // Resource sin function amplitude
#define R_OMEGA         parameter[57]   // Resource sin function phase shift

#define T_MEAN_HEAT     parameter[58]   // Heat wave temperature sin function mean value
#define T_A_HEAT        parameter[59]   // Heat wave temperature sin function amplitude
#define START_HW        parameter[60]   // Start of the heat wave (i.e. day of simulation year)
#define END_HW          parameter[61]   // End of the heat wave (i.e. day of simulation year)
#define T_OMEGA_HEAT    parameter[62]   // Temperature sin heat wave function phase shift

#define T_ADULT_MIN     parameter[63]   // Minimum temperature mature fish experience during winter during seasonally varying temperature using temp_sin function

/*
 *===========================================================================
 *
 *              DEFINING GLOBAL VARIABLES AND USER-IMPLEMENTED FUNCTIONS
 *
 *===========================================================================
 */

extern cohort           abs_tols[POPULATION_NR];

static void             UpdateIDcards(double *env, population *pop);

static int              TimeInYear;

int                     ReproCohorts[POPULATION_NR];
int                     NewCohort = -1;

double                  TotalRepro[POPULATION_NR];
double                  MeanFecundity[POPULATION_NR];

/*#if ((SPECIES_INVADE == 1)) // || (COD_INVADE == 1))
double                  ReproOutputs[INVADEPERIOD][INVADEPERIOD], LastInvasionStart = -1;
#if ((R0_OUTPUT == 1) || (R0_OUTPUT == 2))
static FILE             *reprofile  = NULL;
#endif
#endif*/

#if (LENGTHCURVES > 0)
static double           BirthTimes[LENGTHCURVES];
#endif

/*#if (BIFURCATION == 1)
#define VARIANCES       2                       // Output variances as CVs
//#include  "MeasureBifstats.c"
#endif*/

double          mass[MAXCOHORTS];
double          weight[MAXCOHORTS];
double          length[MAXCOHORTS];
double          sdmmass[MAXCOHORTS];
double          sdmort[MAXCOHORTS];
double          partial_maint[MAXCOHORTS];
double          partial_ingest[MAXCOHORTS];

// Define variables to solve for things like temperature dependent metabolic and ingestion rates
double         rm_weight[MAXCOHORTS];
double         temporary_Tmax_m[MAXCOHORTS], Tmax_m[MAXCOHORTS], temporary_Topt_m[MAXCOHORTS], Topt_m[MAXCOHORTS], temp_Q_m[MAXCOHORTS];
double         Q_m[MAXCOHORTS], log_Q_m[MAXCOHORTS], temp_diff_m[MAXCOHORTS], temp_diff_m_add_2[MAXCOHORTS], Y_m[MAXCOHORTS], W_m[MAXCOHORTS];
double         W_m_pow_2[MAXCOHORTS], div_40_Y_m[MAXCOHORTS], add_1_div_40_Y_m[MAXCOHORTS], pow_point_5_m[MAXCOHORTS];
double         one_pow_point_5_m[MAXCOHORTS], temporary_X_m[MAXCOHORTS], temporary_X_m_2[MAXCOHORTS], X_m[MAXCOHORTS];
double         V_num_m[MAXCOHORTS], V_den_m[MAXCOHORTS], V_m[MAXCOHORTS], V_pow_X_m[MAXCOHORTS], diff_1_V_m[MAXCOHORTS], X_mult_1_V_m[MAXCOHORTS];
double         exp_X_1_V_m[MAXCOHORTS], Rm[MAXCOHORTS];

double         temporary_Tmax_a[MAXCOHORTS], Tmax_a[MAXCOHORTS], temporary_Topt_a[MAXCOHORTS], Topt_a[MAXCOHORTS], temp_Q_a[MAXCOHORTS];
double         Q_a[MAXCOHORTS], log_Q_a[MAXCOHORTS], temp_diff_a[MAXCOHORTS], temp_diff_a_add_2[MAXCOHORTS], Y_a[MAXCOHORTS], W_a[MAXCOHORTS];
double         W_a_pow_2[MAXCOHORTS], div_40_Y_a[MAXCOHORTS], add_1_div_40_Y_a[MAXCOHORTS], pow_point_5_a[MAXCOHORTS];
double         one_pow_point_5_a[MAXCOHORTS], temporary_X_a[MAXCOHORTS], temporary_X_a_2[MAXCOHORTS], X_a[MAXCOHORTS];
double         V_num_a[MAXCOHORTS], V_den_a[MAXCOHORTS], V_a[MAXCOHORTS], V_pow_X_a[MAXCOHORTS], diff_1_V_a[MAXCOHORTS], X_mult_1_V_a[MAXCOHORTS];
double         exp_X_1_V_a[MAXCOHORTS], Ra[MAXCOHORTS];

double         ENV_TEMP, ENV_RESOURCE;

/*
 *===========================================================================
 *
 * UTILITY ROUTINES
 *
 *===========================================================================
 */

double  Sigmoid(double val, double low, double half)
{
  double                result = 0.0, nx;

  // Scale the dependent range between 0 and 1
  nx = 1.5*(val - low)/(half - low);

  if      (nx <= 0.0) result  = 0.0;
  else if (nx <= 1.0) result  = nx*nx*nx/6.0;
  else if (nx <= 2.0) result  = -3.0*nx/2.0 + 3.0*nx*nx/2.0 - nx*nx*nx/3.0 + 0.5;
  else if (nx <= 3.0) result  =  9.0*nx/2.0 - 3.0*nx*nx/2.0 + nx*nx*nx/6.0 - 3.5;
  else                result  = 1.0;

  return result;
}

double Sigmoid_half(double val, double low, double half)
{
  double result = 0.0, nx;

  // Scale the dependent range between 0 and 0.5
  nx = 1.5 * (val - low) / (half - low);

  if      (nx <= 0.0) result = 0.0;
  else if (nx <= 1.0) result = nx * nx * nx / 6.0;
  else if (nx <= 2.0) result = -3.0 * nx / 2.0 + 3.0 * nx * nx / 2.0 - nx * nx * nx / 3.0 + 0.5;
  else if (nx <= 3.0) result = 9.0 * nx / 2.0 - 3.0 * nx * nx / 2.0 + nx * nx * nx / 6.0 - 3.5;
  else                result = 1.0;

  // Scale the result to between 0 and 0.5
  result *= 0.5;

  return result;
}

void    vdPowx(int n, double *x, double b, double *z)
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = pow(x[i],b);

  return;
}

// Functions to solve temperature dependent ingestion and maintenance functions
void    vdPowx_2(int n, double *x, double *b, double *z)
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = pow(x[i],b[i]);

  return;
}

void    vdMul(int n, double x, double *y, double *z) 
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = x*y[i];

  return;
}
   
void    vdMul_2(int n, double *x, double *y, double *z) 
{ 
  int   i;

  for (i=0; i<n; i++)
    z[i] = x[i]*y[i];

  return;
}

void    vdExp(int n, double *x, double *y)
{ 
  int   i;

  for (i=0; i<n; i++)
    y[i] = exp(x[i]);

  return;
}

void    vdLog(int n, double *x, double *y) // removed [i] and *
{
  int   i;

  for(i = 0; i<n;  i++)
    y[i] = log(x[i]); 

  return;
}

void    vdDiff(int n, double *x, double y, double *z) // 
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x[i]-y; 

  return;
}

void    vdDiff_2(int n, double *x, double *y, double *z) // 
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x[i]-y[i]; 

  return;
}

void    vdDiff_3(int n, double x, double *y, double *z) // 
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x-y[i]; 

  return;
}

void    vdAdd(int n, double *x, double y, double *z) // removed * from double y
{
  int   i;

  for(i = 0; i<n; i++)
    z[i] = x[i]+y; 

  return;
}

void    vdDiv(int n, double x, double *y, double *z) //
{
  int   i;

  for(i=0; i<n; i++)
    z[i] = x/y[i];

  return;
}

void    vdDiv_2(int n, double *x, double *y, double *z) //
{
  int   i;

  for(i=0; i<n; i++)
    z[i] = x[i]/y[i];

  return;
}

// Changing, seasonally varying temperature function
double temp_sin(double T)

{
  double temp;
  double PI = 3.141592653589793238463;

  temp = T_MEAN + T_A*sin(2*PI*(T + T_OMEGA) / YEAR); 
  return temp;
}

// Heat wave function
double heat_wave(double T)

{
  double temp;
  double PI = 3.141592653589793238463;

  temp = T_MEAN_HEAT + T_A_HEAT*sin(2*PI*(T + T_OMEGA_HEAT) / YEAR); 
  return temp;
}

// Changing, seasonally varying resource function 
double resource_sin(double T)

{
  double resource;
  double PI = 3.141592653589793238463;

  resource = R_MEAN + R_A*sin(2*PI*(T + R_OMEGA) / YEAR); 
  return resource;

}

/*===========================================================================*/
/*#if ((SPECIES_INVADE == 1)) // || (COD_INVADE == 1))

double DominantEigenvalue(double mat[INVADEPERIOD][INVADEPERIOD], int dim)

// Constructs the next-generation matrix and calculates its dominant eigenvalue

{
  int           ii, jj, Obs, Phase;
  int           n = dim, lda = dim, ldvl = dim, ldvr = dim, info, lwork;
  double        result = -1;
  double        wkopt;
  double        *work;
  double        repout[dim][dim], NextGenMat[dim][dim], wr[dim], wi[dim], vl[dim*dim], vr[dim*dim];
  extern void   dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                        int* lda, double* wr, double* wi, double* vl, int* ldvl,
                        double* vr, int* ldvr, double* work, int* lwork, int* info );

  // Collapse the ReproOutputs matrix to a (dim x dim) matrix (dim=CyclePeriod)
  memset(repout[0], 0, (dim*dim)*sizeof(double));
  for (ii = 0; ii<INVADEPERIOD; ii++)
    {
      // Column index jj refers to reproduction at age jj YEARS
      // With 4 year cycle, first add reproduction at age 4, 8, 12, etc. to reproduction at age 0.
      // Similar for age 1, 2 and 3
      for (jj = dim; jj<INVADEPERIOD; jj++)             // Process columns CyclePeriod-INVADEPERIOD
        mat[ii][jj%dim] += mat[ii][jj];

      // Average rows 0, 4, 8, 12, etc. as they are measurements of the same quantity. Same for other rows.
      // Note that 4 is the 2nd observation (=1+ii/dim) and that ii%dim corresponds to year 0, 1 , 2 or 3
      // in 4-year cycle (the phase)
      Obs   = 1 + ii/dim;
      Phase = ii%dim;
      for (jj = 0; jj<dim; jj++)                        // Process columns 0-CyclePeriod
        UpdateStats(mat[ii][jj], &(repout[Phase][jj]), NULL, NULL, Obs);
    }
  // Matrix repout is condensed ReproOutputs matrix.
  // Rows correspond to phases in cycle, column 0 corresponds to total reproduction at age 0 & 4 & 8 & 12..years,
  // column 1 to total reproduction at age 1 & 5 & 9 & 13 years, etc.
  // NextGenMat should be a phase x phase matrix. Hence, shift the column values in every row to the right
  // row steps (wrap around column dimension)
  for (ii = 0; ii<dim; ii++)
    for (jj = 0; jj<dim; jj++)
      {
        Phase = (ii+jj)%dim;
        NextGenMat[ii][Phase] = repout[ii][jj];
      }

#if (R0_OUTPUT == 2)
  if (reprofile)
    {
      (void)fprintf(reprofile, "Next generation matrix of invaders     :\n");
      for (ii=0; ii<dim; ii++)
        {
          for (jj=0; jj<dim; jj++)
            PrettyPrint(reprofile, NextGenMat[ii][jj]);
          (void)fprintf(reprofile, "\n");
        }
      (void)fprintf(reprofile, "\n");
    }
#endif

  // If period = 1, we should have a single element, corresponding to the life time reproductive effort
  if (dim == 1) return NextGenMat[0][0];

  // Query and allocate the optimal workspace
  lwork = -1;
  dgeev_( "N", "N", &n, &(NextGenMat[0][0]), &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info );
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  // Solve eigenproblem
  dgeev_( "N", "N", &n, &(NextGenMat[0][0]), &lda, wr, wi, vl, &ldvl, vr, &ldvr,   work, &lwork, &info );
  // Check for convergence
  if( info > 0 )
    Warning( "The algorithm failed to compute eigenvalues" );
  else if( info < 0 )
    Warning( "The algorithm to compute eigenvalues was called with an illegal input parameter" );
  else
    {
      result = wr[0];
      for (ii=1; ii<dim; ii++) result = max(result, wr[ii]);
    }

  return result;
}

#endif*/

/*
 *===========================================================================
 *
 * USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
 *
 *===========================================================================
 */

void    UserInit(int argc, char **argv, double *env, population *pop)

{
  register int          i;

  ReportNote("\n    %s%s",
#if (BIFURCATION == 1)
             "Bifurcation of the size-structured consumer-resource model",
#else
             "Dynamics of the size-structured consumer-resource model",
#endif // BIFURCATION
             "\n");

/*#if ((SPECIES_INVADE == 1))
#if (END_INVASION == 1)
  ReportNote("Bifurcation is stopped when R0 > 1.5 or R0 < 0.5");
#else
  ReportNote("Bifurcation is not stopped for any threshold R0");
#endif  // END_INVASION
#endif  // SPECIES_INVADE*/

#if (SPECIES_SINGLE == 1)
  ReportNote("Only a single cohort of species accounted for");

  for (i=0; i<cohort_no[SPECIES]; i++) pop[SPECIES][i][number]=-1.0;

  pop[SPECIES][0][number] = 1.0;
  pop[SPECIES][0][age]    = 0.0;
  pop[SPECIES][0][bone]   = 1/(1+QJ)*BIRTHWEIGHT;
  pop[SPECIES][0][fat]    = (QJ)/(1+QJ)*BIRTHWEIGHT;
  pop[SPECIES][0][gonads] = 0.0;

  for (i=0; i<I_CONST_DIM; i++) popIDcard[SPECIES][0][i] = 0.0;
  
  plankton = FEEDING_LEVEL; // Set plankton to FEEDING_LEVEL; 
#endif

#if (SPECIES_COHORT == 1)
  ReportNote("Multiple cohorts accounted for");
#endif
       
  for (i = 0; i < cohort_no[SPECIES]; i++)
    {
      if (!(iszero(popIDcard[SPECIES][i][IDspawning]) || isequal(popIDcard[SPECIES][i][IDspawning], 1.0)))
        popIDcard[SPECIES][i][IDspawning] = (pop[SPECIES][i][gonads] > 0.0);

      popIDcard[SPECIES][i][IDnbegin] = max(popIDcard[SPECIES][i][IDnbegin], pop[SPECIES][i][number]);
    }

/*#if (BIFURCATION == 1)
//  initMeasureBifstats(argv[1], 0);
#endif */

  // We assume cohort_limit to be 1.0 and that all timing values are in terms of integer number of days
  cohort_limit = 1.0;

  // Fill BirthTimes array with -1.0 (i.e. available entries)
#if (LENGTHCURVES > 0)
  for (i=0; i<LENGTHCURVES; i++) BirthTimes[i] = -1.0;
#endif

  ReportNote("Minimum survival extinction: %E", MINSURVIVAL);
  ReportNote("CHECK_EXT: %d (0 - no action; 1 - message; 2 - message and quit)", CHECK_EXTINCTION);

#if (EGGM_EQ_SDM == 1) 
  ReportNote("egg/yolk-sac mortality follows size-dep. mortality");
  EGGMORT  = MUSDC*exp(-BIRTHWEIGHT/MUSDMASS); // OG mortality function
  YOLKMORT = EGGMORT; // OG mortality function
#endif

#if (EGGM_EQ_SDM == 0)
  ReportNote("No egg/yolk-sac mortality");
#endif

#if (IMMEDIATE_HATCHING == 1)
  ReportNote("Cohort is filled at first day of spawngroup");
#endif

#if (RM_CALCULATION == 0) 
  ReportNote("Rm calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)");
#else     
  ReportNote("Rm calculated with standardized weight (Kooijman, 2010)");  
#endif  // Weight used in Rm calculation

#if (TEMPERATURE == 0)
  ReportNote("Constant environmental temperature");
#else
  ReportNote("Changing environmental temperature: mean temperature = %.1f, amplitude = %.1f, omega = %.1f", T_MEAN, T_A, T_OMEGA); 
#endif

#if (TEMPERATURE == 1 && HEAT_WAVE == 1)
  ReportNote("Marine heat wave: mean temperature = %.1f, amplitude = %.1f, omega = %.1f", T_MEAN_HEAT, T_A_HEAT, T_OMEGA);
#endif

#if (TEMP_ADULT_MIN == 1)
  ReportNote("Minimum winter temp is %.1f C for fish > %.1f cm, during seasonally varying temperature", T_ADULT_MIN, MATURELEN);
#endif

#if (RESOURCE == 0 && SPECIES_SINGLE == 1)
  ReportNote("Constant environmental resource");
#elif (RESOURCE == 1 && SPECIES_SINGLE == 1)
  ReportNote("Changing environmental resource: mean resource = %.1f, amplitude = %.1f, omega = %.1f", R_MEAN, R_A, R_OMEGA);   
#endif

#if (TWO_RESOURCES == 1)
  ReportNote("Two resources with ontogenetic diet shift to shrimp at %.1f cm", HALF_SHRIMP);
#endif

#if (THREE_RESOURCES == 1)
  ReportNote("Three resources with ontogenetic diet shifts to shrimp and fish beginning at %.1f and %.1f cm",
    BEGIN_SHRIMP, BEGIN_FISH); 
#endif

#if (CONSTANT_FISH == 1)
  ReportNote("Constant fish resource = %.1f g/m^3", K_FISH); 
#endif

  SievePop();
  UpdateIDcards(env, pop);

  return;
}

/*
 *===========================================================================
 *
 *      SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
 *
 *===========================================================================
 */

void    SetBpointNo(double *env, population *pop, int *bpoint_no)

{

  bpoint_no[SPECIES] = 0; 

  TimeInYear = floor(fmod(time, YEAR) + MILLI);

#if (SPECIES_SINGLE != 1)
  if (total_eggs > MICRO)
    {
      int       Hatching, CreateCohort, i;
      double    Newborns;

      // Within spawn period? NB: cohort_limit = 1.0!!!
      Hatching     = ((TimeInYear - SPAWNSTART) >= -MILLI) && ((TimeInYear - SPAWNEND) < -(1.0-MILLI));

      if (Hatching)
        {
          // On a spawngroup boundary???
          CreateCohort = ((fmod(TimeInYear - SPAWNSTART, SPAWNGROUP) < MILLI) || 
                          (fmod(TimeInYear - SPAWNSTART, SPAWNGROUP) > SPAWNGROUP - MILLI));

          // Create cohort if both conditions apply
          if (CreateCohort)
            {
              // Compute cumulative newborn produced at the end of the coming spawn group period
              Newborns  = Sigmoid(TimeInYear + SPAWNGROUP, SPAWNSTART, SPAWNMAX);

              // Substract cumulative newborn produced up to now
              Newborns -= Sigmoid(TimeInYear, SPAWNSTART, SPAWNMAX);
              Newborns *= total_eggs;

              NewCohort = AddCohorts(pop, SPECIES, 1);

#if (IMMEDIATE_HATCHING == 0)
              pop[SPECIES][NewCohort][number] = 0.0;
#else // IMMEDIATE_HATCHING == 1
              pop[SPECIES][NewCohort][number] = Newborns;
#endif
              pop[SPECIES][NewCohort][age]    = 0.0;
              pop[SPECIES][NewCohort][bone]   = 1/(1+QJ)*BIRTHWEIGHT;
              pop[SPECIES][NewCohort][fat]    = (QJ)/(1+QJ)*BIRTHWEIGHT;
              pop[SPECIES][NewCohort][gonads] = 0.0;

              for (i=0; i<I_CONST_DIM; i++) popIDcard[SPECIES][NewCohort][i] = 0.0;

              popIDcard[SPECIES][NewCohort][IDnbegin] = Newborns;
              
#if (LENGTHCURVES > 0)
              // Only track the first and last spawn group in a single year
              if (((TimeInYear - SPAWNSTART) < 0.5*SPAWNGROUP) || ((TimeInYear - SPAWNEND) > -1.5*SPAWNGROUP))
                for (i=0; i<LENGTHCURVES; i++)
                  {
                    if (iszero(BirthTimes[i]+1.0))
                      {
                        BirthTimes[i] = time;
                        break;
                      }
                  }
#endif
            }
        }
        NewCohort =  -1; // What is this for? Setting NewCohort back to -1 (i.e. line 247)? Need to read manual
    }
#endif

  return;
} 

/*===========================================================================*/

void    SetBpoints(double *env, population *pop, population *bpoints)


{
    return;
}

/*
 *=====================================================================
 *
 *           SPECIFICATION OF DERIVATIVES
 *
 *=====================================================================
 */

void    Gradient(double *env,     population *pop,     population *ofs,
                 double *envgrad, population *popgrad, population *ofsgrad,
                 population *bpoints)

{
  register int          i;
  double                kappa = 0.0, net_energy = 0.0;
  double                tot_grazing_plankton = 0.0;
  double                tot_grazing_shrimp = 0.0;
  double                tot_grazing_fish= 0.0;
  double                mort = 0.0;

  double                shrimp_fish_frac, zoo_shrimp_frac;


  UpdateIDcards(env, pop);

// #if (BIFURCATION == 1)
//   if (rk_level == 1) measureMinMax(env, pop);
// #endif

  // Species derivatives
  for (i=0; i<cohort_no[SPECIES]; i++)
    {
      popgrad[SPECIES][i][number] =  0.0;
      popgrad[SPECIES][i][age]    =  1.0;
      popgrad[SPECIES][i][bone]   =  0.0;
      popgrad[SPECIES][i][fat]    = -1.0;
      popgrad[SPECIES][i][gonads] =  0.0;
      if (popIDcard[SPECIES][i][IDfatratio] <= MINCONDITION) continue;

      mort       = popIDcard[SPECIES][i][IDmort];
      kappa      = popIDcard[SPECIES][i][IDkappa];
      net_energy = popIDcard[SPECIES][i][IDnet_energy];

      popgrad[SPECIES][i][number] = -mort*pop[SPECIES][i][number];
      if (iszero(popIDcard[SPECIES][i][IDspawning]))
        {
          popgrad[SPECIES][i][bone]       = kappa*net_energy;
          popgrad[SPECIES][i][fat]        = (1-kappa)*net_energy;
        }
      else
        {
          if (net_energy > 0.0)
            {
              popgrad[SPECIES][i][bone]   = kappa*net_energy;
              popgrad[SPECIES][i][gonads] = (1-kappa)*net_energy;
              popgrad[SPECIES][i][fat]    = 0.0;
            }
          else
            {
              if ((pop[SPECIES][i][fat]    > QS*pop[SPECIES][i][bone]) ||
                  (pop[SPECIES][i][gonads] < MICRO))
                {
                  popgrad[SPECIES][i][gonads] = 0.0;
                  popgrad[SPECIES][i][fat]    = net_energy;
                }
              else
                {
                  popgrad[SPECIES][i][gonads] = net_energy;
                  popgrad[SPECIES][i][fat]    = 0.0;
                }
            }
        }

// Two resouces and ontogenetic diet shift at specific length
#if ((SPECIES_SINGLE !=1) && (SPECIES_INVADE !=1) && (TWO_RESOURCES !=0)) // Can we remove the (SPECIES_INVADE !=1) condition?
        if (popIDcard[SPECIES][i][IDlength] < HALF_SHRIMP) 
          tot_grazing_plankton += popIDcard[SPECIES][i][IDingest]*pop[SPECIES][i][number];
        else 
          tot_grazing_shrimp += popIDcard[SPECIES][i][IDingest]*pop[SPECIES][i][number];
#endif

// Three resources, two ontogenentic diet shifts, last diet shift is a percentage of two resources
#if ((SPECIES_SINGLE != 1) && (SPECIES_INVADE !=1) && (THREE_RESOURCES !=0)) // Can we remove the (SPECIES_INVADE !=1) condition?
          if (popIDcard[SPECIES][i][IDlength] < BEGIN_SHRIMP) {
            tot_grazing_plankton += popIDcard[SPECIES][i][IDingest]*pop[SPECIES][i][number];
          }
          else if (popIDcard[SPECIES][i][IDlength] >= BEGIN_SHRIMP && popIDcard[SPECIES][i][IDlength] < END_SHRIMP) {
            zoo_shrimp_frac  = Sigmoid(popIDcard[SPECIES][i][IDlength], BEGIN_SHRIMP, HALF_SHRIMP);
            tot_grazing_plankton += ((popIDcard[SPECIES][i][IDingest])*(1-zoo_shrimp_frac))*pop[SPECIES][i][number];
            tot_grazing_shrimp   += ((popIDcard[SPECIES][i][IDingest])*zoo_shrimp_frac)*pop[SPECIES][i][number];
          }
          else if (popIDcard[SPECIES][i][IDlength] >= END_SHRIMP && popIDcard[SPECIES][i][IDlength] < BEGIN_FISH) {
            tot_grazing_shrimp   += popIDcard[SPECIES][i][IDingest]*pop[SPECIES][i][number];
          }
          else { // Sigmoid transition to 50% fish and 50% shrimp diet
          // tot_grazing_shrimp += (popIDcard[SPECIES][i][IDingest]*PERCENT_SHRIMP)*pop[SPECIES][i][number];
          // tot_grazing_fish   += (popIDcard[SPECIES][i][IDingest]*PERCENT_FISH)*pop[SPECIES][i][number];
          shrimp_fish_frac    = Sigmoid_half(popIDcard[SPECIES][i][IDlength], BEGIN_FISH, HALF_FISH); // 0.0-0.5 not 0-1
          tot_grazing_shrimp += (popIDcard[SPECIES][i][IDingest]*(1-shrimp_fish_frac))*pop[SPECIES][i][number];
          tot_grazing_fish   += (popIDcard[SPECIES][i][IDingest]*shrimp_fish_frac)*pop[SPECIES][i][number];
        }
#endif
      }

#if(IMMEDIATE_HATCHING == 0)
    double                nx, Hatching;

    if (NewCohort >= 0)
    {
      nx = 3.0*(fmod(time, YEAR) - SPAWNSTART)/SPAWNPERIOD;

      if ((nx <= 0.0) || (nx >= 3.0))   Hatching = 0.0;
      else if (nx <= 1.0)               Hatching =                 0.5*nx*nx;
      else if (nx <= 2.0)               Hatching = -1.5 + 3.0*nx -     nx*nx;
      else                              Hatching =  4.5 - 3.0*nx + 0.5*nx*nx;

      Hatching *= total_eggs/(SPAWNPERIOD/3.0);

      mort  = popIDcard[SPECIES][NewCohort][IDmort];

      popgrad[SPECIES][NewCohort][number] = Hatching - mort*pop[SPECIES][NewCohort][number];
      popgrad[SPECIES][NewCohort][age]    = 0.5;
      popgrad[SPECIES][NewCohort][bone]   = 0.0;
      popgrad[SPECIES][NewCohort][fat]    = 0.0;
      popgrad[SPECIES][NewCohort][gonads] = 0.0;
    }
#endif

  // Environment derivatives
  #if (SPECIES_SINGLE == 1) // No feedback of single cohort on resource
    envgrad[0] = 1.0; // time
    envgrad[1] = 0.0; // No change (constant) in PLANKTON = FEEDING_LEVEL
    envgrad[2] = 0.0; // total_eggs
  #elif (SPECIES_COHORT == 1 && TWO_RESOURCES == 1) // Feedback of cohorts on two resources
    envgrad[0] = 1.0; // time
    envgrad[1] = R_PLANKTON * (K_PLANKTON - plankton) - tot_grazing_plankton;
    envgrad[2] = R_SHRIMP * (K_SHRIMP - shrimp) - tot_grazing_shrimp;
    envgrad[3] = 0.0; // total_eggs
  #elif (SPECIES_COHORT == 1 && THREE_RESOURCES == 1) // Feedback of cohorts on three resources
    envgrad[0] = 1.0; // time
    envgrad[1] = R_PLANKTON * (K_PLANKTON - plankton) - tot_grazing_plankton; 
    envgrad[2] = R_SHRIMP * (K_SHRIMP - shrimp) - tot_grazing_shrimp;
    envgrad[3] = R_FISH * (K_FISH - fish) - tot_grazing_fish;
    envgrad[4] = 0.0; // total_eggs
  #endif

  return;
}

/*
 *===========================================================================
 *
 *              SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
 *
 *===========================================================================
 */

void    InstantDynamics(double *env, population *pop, population *ofs)

{
  register int          i, pn;

#if ((LENGTHCURVES > 0))
  int                   j;
#endif

#if (SPECIES_SINGLE == 1) 
  if (FEEDING_LEVEL < FOODTINY) FEEDING_LEVEL = FOODTINY; // Check negative food
#else 
  if (plankton < FOODTINY) plankton = FOODTINY; // Check negative food
#endif

  UpdateIDcards(env, pop);

  for (pn=0; pn<POPULATION_NR; pn++)
    for (i=0; i<cohort_no[pn]; i++)             // remove cohorts with no fat
      {
        if ((popIDcard[pn][i][IDfatratio] > MINCONDITION) &&
            (pop[pn][i][number] > abs_tols[pn][0]) &&
            (pop[pn][i][number] > MINSURVIVAL*popIDcard[pn][i][IDnbegin])) continue;

        pop[pn][i][number] = 0.0;
#if ((LENGTHCURVES))
        switch (pn)
          {
            case 0:
#if (LENGTHCURVES > 0)
              for (j=0; j<LENGTHCURVES; j++)
                {
                  if (fabs(time - BirthTimes[j] - pop[pn][i][age]) < 0.5) // cohort_limit is 1.0
                    {
                      BirthTimes[j] = -2.0;
                      break;
                    }
                }
#endif
              break;
            case 1:
              break;
            }
#endif
      }

#if(IMMEDIATE_HATCHING == 0)
  // Ensure that the newborn cohort is not being discarded
  if (NewCohort >= 0)
    pop[SPECIES][NewCohort][number] = max(pop[SPECIES][NewCohort][number], (1.0 + MILLI)*abs_tols[SPECIES][0]);
#endif

  TimeInYear = floor(fmod(time, YEAR) + MILLI);

  // Species reproduction
  if (isequal(TimeInYear, SPAWNSET) ||
      ((iszero(SPAWNSET) || isequal(SPAWNSET, YEAR)) && isequal(TimeInYear, YEAR)))
    {
      for (i=0; i<cohort_no[SPECIES]; i++)
        {
        if (pop[SPECIES][i][number] <= abs_tols[SPECIES][0]) continue;

        // The following does not lead to cohort splitting due to skipped spawning
          if ((popIDcard[SPECIES][i][IDlength]   >= MATURELEN) &&
              (popIDcard[SPECIES][i][IDfatratio] >= QSPAWN))
            {
              popIDcard[SPECIES][i][IDspawning] = 1.0;
              pop[SPECIES][i][gonads]   = pop[SPECIES][i][fat];
              pop[SPECIES][i][fat]      = QSPAWN*pop[SPECIES][i][bone];
              pop[SPECIES][i][gonads]  -= pop[SPECIES][i][fat];
            }
          else
            popIDcard[SPECIES][i][IDspawning] = 0.0;
        }
    }

  if (isequal(TimeInYear, SPAWNSTART) ||
      ((iszero(SPAWNSTART) || isequal(SPAWNSTART, YEAR)) && isequal(TimeInYear, YEAR)))
    {
      double            fecundity;

      total_eggs            = 0.0;
      ReproCohorts[SPECIES] = 0;
      TotalRepro[SPECIES]   = 0.0;

      /*
       * NB: Note that the initial values for bone and fat of newly spawned
       * eggs is set to the bone and fat values of a first-feeding
       * larvae. Hence, the conversion coefficient REPROEFF should
       * account for the energy content of the yolk sac that is used to
       * maintain the egg and yolk-sac larvae during the period of
       * non-feeding
       */
      for (i=0; i<cohort_no[SPECIES]; i++)
        {
          if (pop[SPECIES][i][number] <= abs_tols[SPECIES][0]) continue;

          if (!iszero(popIDcard[SPECIES][i][IDspawning]))
            {
              fecundity = REPROEFF*max(pop[SPECIES][i][gonads], 0.0)/BIRTHWEIGHT;

              popIDcard[SPECIES][i][IDfecundity] = fecundity;
              popIDcard[SPECIES][i][IDlifefec]  += (fecundity*pop[SPECIES][i][number]/popIDcard[SPECIES][i][IDnbegin]);

/*#if (SPECIES_INVADE == 1)
              int       row, col;

              // row: the year number since start of invasions
              // col: age in years of reproducing cohort
              //      (1 is added because first possible reproduction is at age less than 1 year!)
              row = (int)floor((time - (LastInvasionStart + pop[SPECIES][i][age] - 0.5*cohort_limit))/YEAR);
              col = 1 + (int)floor((pop[SPECIES][i][age]-0.5*cohort_limit)/YEAR);
              ReproOutputs[row][col]          += (fecundity*pop[SPECIES][i][number]/INVADEDENSITY);
#else*/
              total_eggs                        += fecundity*pop[SPECIES][i][number];
// #endif

              ReproCohorts[SPECIES]++;
              TotalRepro[SPECIES]               += pop[SPECIES][i][number];
            }
          else
            popIDcard[SPECIES][i][IDfecundity]   = 0.0;

            popIDcard[SPECIES][i][IDspawning]    = 0.0;
            pop[SPECIES][i][gonads]              = 0.0;
        }

      if (TotalRepro[SPECIES]) MeanFecundity[SPECIES] = total_eggs/TotalRepro[SPECIES];
      else MeanFecundity[SPECIES] = 0.0;

/*#if ((SPECIES_INVADE == 1) && (BIFURCATION == 1))
      double            dummy = fmod(time, BifPeriod);
      if (((dummy + 0.5*cohort_limit) > (BifPeriod - 2*INVADEPERIOD*YEAR)) &&
          ((dummy + 0.5*cohort_limit) < (BifPeriod -   INVADEPERIOD*YEAR)))
        {
          total_eggs = INVADEDENSITY;

          if (LastInvasionStart < 0)
            {
              LastInvasionStart = time;
              memset((void *)ReproOutputs, 0, (INVADEPERIOD*INVADEPERIOD)*sizeof(double));
            }
        }
#endif*/
    }

  return;
} 


/*
 *===========================================================================
 *
 *                      SPECIFICATION OF OUTPUT VARIABLES
 *
 *===========================================================================
 */

void    DefineOutput(double *env, population *pop, double *output)

{
  register int          i, index;
  double                num, biom;

#if (LENGTHCURVES > 0)
  int                   j;
  double                Lengths[LENGTHCURVES];

  for (j=0; j<LENGTHCURVES; j++) Lengths[j] = DBL_MAX;
#endif

  UpdateIDcards(env, pop);

/*#if (BIFURCATION == 1)
  TimeInYear = (int)floor(fmod(time, YEAR) + MILLI);

  LabelState(SPECIES, "par. = %.4f   T = %6.0f years, %3d days", parameter[BifParIndex],
             floor((env[0]+0.1*cohort_limit)/YEAR), TimeInYear);
 #endif // BIFURCATION*/

#define OUTONE 1 // OUTONE 3

#if (SPECIES_SINGLE == 1)
  output[ 0] = FEEDING_LEVEL; // Output column two; first column is set to time by default
  output[OUTONE+ 0] = pop[SPECIES][0][number];
  output[OUTONE+ 1] = pop[SPECIES][0][bone];
  output[OUTONE+ 2] = pop[SPECIES][0][fat];
  output[OUTONE+ 3] = pop[SPECIES][0][gonads];
  output[OUTONE+ 4] = popIDcard[SPECIES][0][IDweight];
  output[OUTONE+ 5] = popIDcard[SPECIES][0][IDlength];
  output[OUTONE+ 6] = popIDcard[SPECIES][0][IDlifefec];
  output[OUTONE+ 7] = popIDcard[SPECIES][0][IDfatratio];
  output[OUTONE+ 8] = popIDcard[SPECIES][0][IDkappa];
  
  // Adjusted by Bass
  output[OUTONE+ 9] = time/YEAR; 
  output[OUTONE+10] = popIDcard[SPECIES][0][IDmaint];
  output[OUTONE+11] = popIDcard[SPECIES][0][IDingest];
  output[OUTONE+12] = popIDcard[SPECIES][0][IDnet_energy];
  output[OUTONE+13] = popIDcard[SPECIES][0][IDmort];
  output[OUTONE+14] = Rm[0];
  output[OUTONE+15] = Ra[0];
  output[OUTONE+16] = ENV_TEMP;
  output[OUTONE+17] = popIDcard[SPECIES][0][IDmass];
  output[OUTONE+18] = ENV_RESOURCE;

#elif(SPECIES_COHORT == 1) // #elif(SPECIES_INVADE !=1) - Bass removed this and put SPECIES_COHORT; can change back later 
  output[ 0] = plankton; // Output column two; first column is set to time by default

  for (i=0; i<cohort_no[SPECIES]; i++)
    {
      num  = pop[SPECIES][i][number];
      biom = (pop[SPECIES][i][bone] + pop[SPECIES][i][fat] + pop[SPECIES][i][gonads])*pop[SPECIES][i][number];

      if (pop[SPECIES][i][age]< (YEAR - MILLI))           // YOY
        {
          output[OUTONE+ 0] += num; 
          output[OUTONE+ 1] += biom; 
        }
      else
        {
          output[OUTONE+ 2] += num;                      // 1+
          output[OUTONE+ 3] += biom;
        }
      if(popIDcard[SPECIES][i][IDlength] < MATURELEN)    // Juveniles
        {
          output[OUTONE+ 4] += num; 
          output[OUTONE+ 5] += biom; 
        }
      else                                               // Adults
        {
          output[OUTONE+ 6] += num; 
          output[OUTONE+ 7] += biom; 
        }

      output[OUTONE+ 8] += num;                          // Total population
      output[OUTONE+ 9] += biom;   

#if (LENGTHCURVES > 0)
      for (j=0; j<LENGTHCURVES; j++) 
        {
          if (fabs(time - BirthTimes[j] - pop[SPECIES][i][age]) < 0.5) // cohort_limit is 1.0
            {
              Lengths[j] = popIDcard[SPECIES][i][IDlength];
              break;
            }
        }
#endif
    }

  output[OUTONE+10] = fish; // popIDcard[SPECIES][0][IDlength]; 
  output[OUTONE+11] = ENV_TEMP; // popIDcard[SPECIES][0][IDlifefec] ~ Same value as total_eggs += fecundity*pop[SPECIES][i][number] on line 931
  output[OUTONE+12] = time/YEAR; 
  output[OUTONE+13] = shrimp; // Column # 16 in output

index = 15;  
#if (LENGTHCURVES > 0)
  for (j=0; j<LENGTHCURVES; j++) 
    {
      output[index++]  = Lengths[j]; // 
      BirthTimes[j]   = max(BirthTimes[j], -1.0);
    }
#endif

#endif 

/*#if (BIFURCATION == 1)
//  outputMeasureBifstats(env, NULL, OUTPUT_VAR_NR, 0);

#if ((SPECIES_INVADE == 1))
 static int             first = 1, direction = 0, fftIndex = 0, outputIndex = OUTONE+12;

#if (SPECIES_INVADE == 1)
 fftIndex    = 2;
 outputIndex = OUTONE+12;
#endif

 double dummy = fmod(env[0], BifPeriod);
 if ((env[0] > 0.5*cohort_limit) && BifPeriodEnd && ((dummy < 0.5*cohort_limit) || (dummy > (BifPeriod-0.5*cohort_limit))))
   {
     int CyclePeriod = (int)(rint(periodFFT[fftIndex]/YEAR) + MILLI);
#if ((R0_OUTPUT == 1) || (R0_OUTPUT == 2))
     if (first)
       {
         strcpy(fn, currun);                            // Open additional file for R0 output
         strcat(fn, ".R0.out");
         reprofile = fopen(fn, "a");
         if (!reprofile)
           fprintf(stderr, "Failed to open file %s\n", fn);
       }
     first = 0;
     if (reprofile)
       {
         int            ii, jj;

         (void)fprintf(reprofile, "Time                                   : ");
         PrettyPrint(reprofile, env[0]);
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "Current value of bifurcation parameter : ");
         PrettyPrint(reprofile, parameter[BifParIndex]);
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "Dominant period of population cycle    : ");
         PrettyPrint(reprofile, CyclePeriod);
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "\n");
         (void)fprintf(reprofile, "Full reproduction matrix of invaders   :\n");
         for (ii=0; ii<INVADEPERIOD; ii++)
           {
             for (jj=0; jj<INVADEPERIOD; jj++)
               PrettyPrint(reprofile, ReproOutputs[ii][jj]);
             (void)fprintf(reprofile, "\n");
           }
         (void)fprintf(reprofile, "\n");
       }
#endif
     if (!CyclePeriod) CyclePeriod = 1;
     output[outputIndex]  = DominantEigenvalue(ReproOutputs, CyclePeriod);
#if ((R0_OUTPUT == 1) || (R0_OUTPUT == 2))
     (void)fprintf(reprofile, "Dominant eigenvalue                    : ");
     PrettyPrint(reprofile, output[outputIndex]);
     (void)fprintf(reprofile, "\n\n\n");
#endif
#if (END_INVASION == 1)
     if (!direction)
       {
          if (output[outputIndex] < 1.0) direction = 1;
          else direction = -1;
       }
     else if (direction == 1)
       {
          if (output[outputIndex] > 1.5) ForcedRunEnd = 1;
       }
     else
       {
          if (output[outputIndex] < 0.5) ForcedRunEnd = 1;
       }
#endif

     // Reset invasion start time to -1
     LastInvasionStart = -1;

#if (DEBUG > 0)
     fprintf(stderr, "R0: %.10f\n", output[outputIndex]);
     ForcedRunEnd = 1;
#endif
   }

#else // i.e. #if ((SPECIES_INVADE != 1)

#endif // #if (SPECIES_INVADE == 1)

#endif // #if (BIFURCATION == 1)

#if (DEBUG > 0)
//  fprintf(stderr, "Time: %6.0f\tCohorts[0]: %5d\tCohorts[1]: %5d\n", env[0], cohort_no[0], cohort_no[1]);
#endif
*/

  return;
}

/*==============================================================================*/

static void     UpdateIDcards(double *env, population *pop)

{
  register int          i;
  double                fatratio;
  double                ingest, simple_ingest;
  double                maint, net_energy, kappa, q;
  double                mort, eggmort, yolkmort;
  double                hatchfrac, shrimp_fish_frac, zoo_shrimp_frac;    
  double                temp_TEMP;             

  // Initialize arrays for use in vdpowx
  for (i=0; i<cohort_no[SPECIES]; i++)
    {
      mass[i]      = pop[SPECIES][i][bone]*(1 + QJ);
      weight[i]    = pop[SPECIES][i][bone] + max(pop[SPECIES][i][fat] + pop[SPECIES][i][gonads], 0.0);

      if (MUSDC > 0) 
        sdmmass[i] = -mass[i]/MUSDMASS;

      popIDcard[SPECIES][i][IDmass]     = mass[i];
      popIDcard[SPECIES][i][IDweight]   = weight[i];
      popIDcard[SPECIES][i][IDfatratio] = (pop[SPECIES][i][fat] + pop[SPECIES][i][gonads])/pop[SPECIES][i][bone];

      // Constant resource level or seasonally varying constant level
      #if (RESOURCE == 0) 
          ENV_RESOURCE = plankton; // plankton is set to FEEDING_LEVEL on line 604
      #else 
          ENV_RESOURCE = resource_sin(time); 
      #endif

      // Set fish resource to constant value i.e. remove semichemostat growth
      #if (CONSTANT_FISH == 1)
        fish = K_FISH;
      #endif

      // Water temperature (constant, seasonally varying and or heat wave)
      #if (TEMPERATURE == 0) // Constant temperature
          ENV_TEMP = TEMP;  
      #elif (TEMPERATURE == 1 && HEAT_WAVE == 1) // Seasonally varying temperature with added heatwave
          if (TimeInYear > START_HW && TimeInYear < END_HW) 
            ENV_TEMP = heat_wave(time);
          else 
            ENV_TEMP = temp_sin(time);
      #elif (TEMPERATURE == 1 && TEMP_ADULT_MIN == 1) // Mature fish subjected to minimum winter temperature 
          temp_TEMP = temp_sin(time); // Seasonally varying temperature
          if (length[i] >= MATURELEN && temp_TEMP <= T_ADULT_MIN)
            ENV_TEMP = T_ADULT_MIN; 
          else
            ENV_TEMP = temp_TEMP; // Seasonally varying temperature
      #else 
            ENV_TEMP = temp_sin(time);; // Seasonally varying temperature
      #endif

      // Rm calculation
      #if (RM_CALCULATION == 0) 
          rm_weight[i] = weight[i];  // Rm calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)
      #else     
          rm_weight[i] = mass[i];    // Rm calculated with standardized weight (Kooijman, 2010)
      #endif  

    }

      // Length and part of simple ingestion calcuation
      vdPowx(cohort_no[SPECIES], mass, LWE, length);  // pow(mass, LWE)
      vdPowx(cohort_no[SPECIES], mass, INGEST_E, partial_ingest); // pow(mass, INGEST_C); solve for part of the simplifed ingestion equation, mass ~ standarized weight 

      // Size-dependent mortality calculation 
      if (MUSDC > 0) 
        vdExp(cohort_no[SPECIES], sdmmass, sdmort);  // exp(-mass[i]/MUSDMASS) 
      
      // Rm calculation with total or standardized weight
      #if (RM_CALCULATION == 0)    // Maintenance calculated with weight = bones + fat + gonads (Ohlberger et al. 2011)  
          vdPowx(cohort_no[SPECIES], weight, MAINTE, partial_maint);  // pow(weight, MAINTE); solve for part of maintenance equation
      #else                        // Maintenance calculated with standardized weight  
          vdPowx(cohort_no[SPECIES], mass, MAINTE, partial_maint);   // pow(mass, MAINTE); solve for part of maintenance equation
      #endif  

      // Steps to calculate Rm - temperature dependent metabolic factor 
      // Solve for Tmax_m which is function x + y
      vdPowx(cohort_no[SPECIES], rm_weight, VMAX_M, temporary_Tmax_m); // (x + y)^VMAX_M... temp_Tmax_m is z[i]
      vdMul(cohort_no[SPECIES], YMAX_M, temporary_Tmax_m, Tmax_m); // YMAX_M * (x + y)^VMAX_M... Tmax_m is z[i]
      
      // Solve for Topt_m which is function x + y
      vdPowx(cohort_no[SPECIES], rm_weight, VOPT_M, temporary_Topt_m); // (x + y)^VOPT_M...temp_Topt_m is z[i]
      vdMul(cohort_no[SPECIES], YOPT_M, temporary_Topt_m, Topt_m); // YOPT_M * (x + y)^VOPT_M... Topt_m is z[i]
      
      // Solve for Q_m = GA_M * (x + y)^THETA_M
      vdPowx(cohort_no[SPECIES], rm_weight, THETA_M, temp_Q_m); // x + y)^THETA_M... temp_Q_m is z[i]
      vdMul(cohort_no[SPECIES], GA_M, temp_Q_m, Q_m); // GA_M * (x + y)^THETA_M .. Q_m is z[i]
      
      // Solve for Y_m = (Tmax_m - Topt_m + 2) * log(Q_m)
      vdLog(cohort_no[SPECIES], Q_m, log_Q_m); // log(Q_m) = log(GA_M * (x + y)^THETA_M) ... log_Q_m is y[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_m, Topt_m, temp_diff_m); // Tmax_m - Topt_m ...temp_diff_m is z[i]
      vdAdd(cohort_no[SPECIES], temp_diff_m, 2, temp_diff_m_add_2); // Tmax_m - Topt_m + 2 ... temp_diff_m_add_2 is z[i]
      vdMul_2(cohort_no[SPECIES], temp_diff_m_add_2, log_Q_m, Y_m); // (Tmax_m - Topt_m + 2) * log(Q_m) ... Y_m is z[i]

      // Solve for W_m = (Tmax_m - Topt_m) * log(Q_m)
      vdMul_2(cohort_no[SPECIES], temp_diff_m, log_Q_m, W_m); // (Tmax_m - Topt_m) * log(Q_m) .. W_m[i]

      // Solve for X_m = W_m^2 * (1/400) * (1 + (1 + (40 / Y_m))^0.5)^2
      vdPowx(cohort_no[SPECIES], W_m, 2, W_m_pow_2); // W_m^2 ...W_m_pow_2 is z[i]
      vdDiv(cohort_no[SPECIES], 40, Y_m, div_40_Y_m); // 40 / Y_m ...Div_40_Y_m is z[i]
      vdAdd(cohort_no[SPECIES], div_40_Y_m, 1, add_1_div_40_Y_m); // 1 + 40/Y_m   ..add_1_div_40_Y_m is z[i]
      vdPowx(cohort_no[SPECIES], add_1_div_40_Y_m, 0.5, pow_point_5_m); // (1 + 40/Y_m)^0.5 ...pow_point_5_m is z[i]
      vdAdd(cohort_no[SPECIES], pow_point_5_m, 1, one_pow_point_5_m); // 1 + (1 + 40/Y_m)^0.5 ...one_pow_point_5_m is z[i]
      vdPowx(cohort_no[SPECIES], one_pow_point_5_m, 2, temporary_X_m); // (1 + (1 + 40/Y_m)^0.5)^2 ..temporary_X_m is z[i]
      vdMul(cohort_no[SPECIES], DIV_1_400, temporary_X_m, temporary_X_m_2); // (1/400) * (1 + (1 + 40/Y_m)^0.5)^2 ...temporary_X_m_2 is z[i]
      vdMul_2(cohort_no[SPECIES], W_m_pow_2, temporary_X_m_2, X_m); // W_m^2 * (1/400) * (1 + (1 + (40 / Y_m))^0.5)^2... X_m is z[i]

      // Solve for V_m = (Tmax_m - TEMP) / (Tmax_m - Topt_m); TEMP is our specific environmental temperature set in parameter file
      vdDiff(cohort_no[SPECIES], Tmax_m, ENV_TEMP, V_num_m); // Tmax_m - TEMP ...V_num_m is z[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_m, Topt_m, V_den_m); // Tmax_m - Topt_m ...V_den_m is z[i]
      vdDiv_2(cohort_no[SPECIES], V_num_m, V_den_m, V_m); // (Tmax_m - temp) / (Tmax_m - Topt_m) ...V_m is z[i]

      // Solve for Rm - temperature dependent metabolic factor 
      // Rm = V_m^(X_m) * (exp(X_m * 1 - V_m))) 
      vdPowx_2(cohort_no[SPECIES], V_m, X_m, V_pow_X_m); // V_m^X_m   ... V_pow_X_m is z[i]
      vdDiff_3(cohort_no[SPECIES], 1, V_m, diff_1_V_m); // 1 - V_m  ... diff_1_V_m is z[i]
      vdMul_2(cohort_no[SPECIES], X_m, diff_1_V_m, X_mult_1_V_m); // X_m * (1 - V_m)  ...X_mult_1_V_m is z[i]
      vdExp(cohort_no[SPECIES], X_mult_1_V_m, exp_X_1_V_m); // exp(X_m * (1 - V_m)) ...exp_X_1_V_m is z[i]
      vdMul_2(cohort_no[SPECIES], V_pow_X_m, exp_X_1_V_m, Rm); // V^(X) * (exp(X * (1 - V))) ... Rm is z[i]


      // Steps to calculate Ra - temperature dependent intake factor 
      // Solve for Tmax_a which is function x(1+qj)
      vdPowx(cohort_no[SPECIES], mass, VMAX_A, temporary_Tmax_a); // (standardized mass)^VMAX_a... temp_Tmax_a is z[i]
      vdMul(cohort_no[SPECIES], YMAX_A, temporary_Tmax_a, Tmax_a); // YMAX_a * (standardized mass)^VMAX_a... Tmax_a is z[i]
      
      // Solve for Topt_a which is function x + y
      vdPowx(cohort_no[SPECIES], mass, VOPT_A, temporary_Topt_a); // (standardized mass)^VOPT_a...temp_Topt_a is z[i]
      vdMul(cohort_no[SPECIES], YOPT_A, temporary_Topt_a, Topt_a); // YOPT_a * (standardized mass)^VOPT_a... Topt_a is z[i]
      
      // Solve for Q_a = GA_a * (x + y)^THETA_a
      vdPowx(cohort_no[SPECIES], mass, THETA_A, temp_Q_a); // standardized mass)^THETA_a... temp_Q_a is z[i]
      vdMul(cohort_no[SPECIES], GA_A, temp_Q_a, Q_a); // GA_a * (standardized mass)^THETA_a .. Q_a is z[i]
      
      // Solve for Y_a = (Tmax_a - Topt_a + 2) * log(Q_a)
      vdLog(cohort_no[SPECIES], Q_a, log_Q_a); // log(Q_a) = log(GA_a * (standardized mass)^THETA_a) ... log_Q_a is y[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_a, Topt_a, temp_diff_a); // Tmax_a - Topt_a ...temp_diff_a is z[i]
      vdAdd(cohort_no[SPECIES], temp_diff_a, 2, temp_diff_a_add_2); // Tmax_a - Topt_a + 2 ... temp_diff_a_add_2 is z[i]
      vdMul_2(cohort_no[SPECIES], temp_diff_a_add_2, log_Q_a, Y_a); // (Tmax_a - Topt_a + 2) * log(Q_a) ... Y_a is z[i]

      // Solve for W_a = (Tmax_a - Topt_a) * log(Q_a)
      vdMul_2(cohort_no[SPECIES], temp_diff_a, log_Q_a, W_a); // (Tmax_a - Topt_a) * log(Q_a) .. W_a[i]

      // Solve for X_a = W_a^2 * (1/400) * (1 + (1 + (40 / Y_a))^0.5)^2
      vdPowx(cohort_no[SPECIES], W_a, 2, W_a_pow_2); // W_a^2 ...W_a_pow_2 is z[i]
      vdDiv(cohort_no[SPECIES], 40, Y_a, div_40_Y_a); // 40 / Y_a ...Div_40_Y_a is z[i]
      vdAdd(cohort_no[SPECIES], div_40_Y_a, 1, add_1_div_40_Y_a); // 1 + 40/Y_a   ..add_1_div_40_Y_a is z[i]
      vdPowx(cohort_no[SPECIES], add_1_div_40_Y_a, 0.5, pow_point_5_a); // (1 + 40/Y_a)^0.5 ...pow_point_5_a is z[i]
      vdAdd(cohort_no[SPECIES], pow_point_5_a, 1, one_pow_point_5_a); // 1 + (1 + 40/Y_a)^0.5 ...one_pow_point_5_a is z[i]
      vdPowx(cohort_no[SPECIES], one_pow_point_5_a, 2, temporary_X_a); // (1 + (1 + 40/Y_a)^0.5)^2 ..temporary_X_a is z[i]
      vdMul(cohort_no[SPECIES], DIV_1_400, temporary_X_a, temporary_X_a_2); // (1/400) * (1 + (1 + 40/Y_a)^0.5)^2 ...temporary_X_a_2 is z[i]
      vdMul_2(cohort_no[SPECIES], W_a_pow_2, temporary_X_a_2, X_a); // W_a^2 * (1/400) * (1 + (1 + (40 / Y_a))^0.5)^2... X_a is z[i]

      // Solve for V_a = (Tmax_a - TEMP) / (Tmax_a - Topt_a); TEMP is our specific environmental temperature set in parameter file
      vdDiff(cohort_no[SPECIES], Tmax_a, ENV_TEMP, V_num_a); // Tmax_a - TEMP ...V_num_a is z[i]
      vdDiff_2(cohort_no[SPECIES], Tmax_a, Topt_a, V_den_a); // Tmax_a - Topt_a ...V_den_a is z[i]
      vdDiv_2(cohort_no[SPECIES], V_num_a, V_den_a, V_a); // (Tmax_a - temp) / (Tmax_a - Topt_a) ...V_a is z[i]

      // Solve for Ra - temperature dependent intake factor 
      // Ra = V_a^(X_a) * (exp(X_a * 1 - V_a)))
      vdPowx_2(cohort_no[SPECIES], V_a, X_a, V_pow_X_a); // V_a^X_a   ... V_pow_X_a is z[i]
      vdDiff_3(cohort_no[SPECIES], 1, V_a, diff_1_V_a); // 1 - V_a  ... diff_1_V_a is z[i]
      vdMul_2(cohort_no[SPECIES], X_a, diff_1_V_a, X_mult_1_V_a); // X_a * (1 - V_a)  ...X_ault_1_V_a is z[i]
      vdExp(cohort_no[SPECIES], X_mult_1_V_a, exp_X_1_V_a); // exp(X_a * (1 - V_a)) ...exp_X_1_V_a is z[i]
      vdMul_2(cohort_no[SPECIES], V_pow_X_a, exp_X_1_V_a, Ra); // V^(X) * (exp(X * (1 - V))) ... Ra is z[i]

// Calculate egg and yolk mortality
#if (EGGM_EQ_SDM == 1)
  eggmort  = MUSDC*exp(-BIRTHWEIGHT/MUSDMASS); // OG mortality
  yolkmort = eggmort; // OG mortality
#else
  eggmort  = EGGMORT; // Parameter egg mortality
  yolkmort = YOLKMORT; // Parameter yolk moartlity 
#endif

  // Update IDcards for species
  for (i=0; i<cohort_no[SPECIES]; i++)
  {
    length[i]    *= LWC;
    fatratio      = popIDcard[SPECIES][i][IDfatratio];
    mort          = MUB; // Background mortality 

    simple_ingest = INGEST_C*partial_ingest[i]*Ra[i]; // simple_ingest = INGEST_C * MASS ^ INGEST_E * Ra 
    maint         = Rm[i]*MAINTC*partial_maint[i]; // maintenance;  Rm is temperature adjustment factor
    net_energy    = -maint;

    ingest        = 0.0;
    kappa         = 0.0;

  // Mortality
    if (pop[SPECIES][i][age] < EGGPERIOD)
    {
      // Egg mortality
      mort       = eggmort;
      maint      = 0.0;
      net_energy = 0.0;
    }
    else if (pop[SPECIES][i][age] < AGEFEEDING)
    {
      // Yolk mortality 
      hatchfrac  = Sigmoid(pop[SPECIES][i][age], EGGPERIOD, EGGPERIOD+0.5*SPAWNGROUP);
      mort       = hatchfrac*yolkmort + (1-hatchfrac)*eggmort;
      maint      = 0.0;
      net_energy = 0.0;
    }
    else if (fatratio > MINCONDITION) 
    {
      // Size-dependent mortality 
      if (MUSDC > 0)  
          mort += MUSDC*sdmort[i]; // MUSDC*exp(-mass[i]/MUSDMASS) 
      // Starvation mortality
      if (fatratio < QS) 
          mort += MUS*(QS/fatratio - 1.0); // s * (qstarv * x / y - 1)

  // Ingestion & resources
    #if (SPECIES_SINGLE == 1) // Single cohort, no feedback on environment
          ingest     = ENV_RESOURCE*simple_ingest; // Ra accounted for in line 1358

    #elif (SPECIES_COHORT == 1 && TWO_RESOURCES == 1) // Multiple cohorts, feedback on environment, two resources with diet shift
      if (length[i] < HALF_SHRIMP) // Ontogenetic diet shift length 
          ingest     = plankton*simple_ingest; 
      else 
          ingest     = shrimp*simple_ingest; 

    #elif (SPECIES_COHORT == 1 && THREE_RESOURCES == 1) // Multiple cohorts, feedback on environment, three resources with two diet shifts
      if (length[i] < BEGIN_SHRIMP) { // Feeding only on plankton 
        ingest           = plankton*simple_ingest; 
      }
      else if (length[i] >= BEGIN_SHRIMP && length[i] < END_SHRIMP) { // Transition between plankton and shrimp
        zoo_shrimp_frac  = Sigmoid(length[i], BEGIN_SHRIMP, HALF_SHRIMP); 
        ingest           = (plankton*simple_ingest*(1-zoo_shrimp_frac)) + ((shrimp*simple_ingest)*zoo_shrimp_frac); 
      } 
      else if (length[i] >= END_SHRIMP && length[i] < BEGIN_FISH) { // Feeding only on shrimp
        ingest           = shrimp*simple_ingest; 
      }
      else { // Sigmoid transition to 50% fish and 50% shrimp diet
        shrimp_fish_frac = Sigmoid_half(length[i], BEGIN_FISH, HALF_FISH); // scaled from 0.0-0.5 not 0-1
        ingest           = ((shrimp*simple_ingest)*(1-shrimp_fish_frac)) + ((fish*simple_ingest)*shrimp_fish_frac); 
      }
      #endif

  // Net energy, kappa
          net_energy = (CONVEFF*ingest - maint); // net_energy = (ke*(1/d1*w^d2)*Ra[i]*FEEDING_LEVEL) - (Rm[i]*MAINTC*maint[i])

          if (length[i] < MATURELEN) 
            q = QJ; // Specifies energy allocation for juveniles or adults
          else 
            q = QA;

          kappa = 0.0;
          if (net_energy > 0.0)
          {
              if (fatratio > q) // fat ratio used for these equations, same as our equations
                kappa = 1/(q+1);
              else
                // Default recovery formula changed to a quadratic dependence on fatratio (below), such that recovery in
                // reversible mass is faster 
                kappa = fatratio*fatratio/((1+q)*q*q);
            }

            hatchfrac   = Sigmoid(pop[SPECIES][i][age], AGEFEEDING, AGEFEEDING+0.5*SPAWNGROUP);
            ingest     *= hatchfrac;
            net_energy *= hatchfrac;
            maint      *= hatchfrac;
            mort       *= hatchfrac;
            mort       += (1.0 - hatchfrac)*yolkmort;
          }

          popIDcard[SPECIES][i][IDlength]     = length[i];
          popIDcard[SPECIES][i][IDmaint]      = maint;
          popIDcard[SPECIES][i][IDingest]   = ingest;
          popIDcard[SPECIES][i][IDnet_energy] = net_energy;
          popIDcard[SPECIES][i][IDkappa]      = kappa;
          popIDcard[SPECIES][i][IDmort]       = mort;
        }

  return;
}  



/*=============================================================================*/
