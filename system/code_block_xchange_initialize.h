/* this is a 'pre-amble' block which must be included before
    the code in 'code_block_xchange...(any of the other routines).h' is called.
    It sets a number of variable and function names based on the user-specified
    master routine name (which should of course itself be unique!). It also sets
    up some definitions for the multi-threading, and pre-defines some of the
    subroutines which will be referenced. */

#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

#define INPUT_STRUCT_NAME MASTER_FUNCTION_NAME##_data_in /* dummy name - must be unique */
#define DATAIN_NAME MASTER_FUNCTION_NAME##_DataIn /* dummy name - must be unique */
#define DATAGET_NAME MASTER_FUNCTION_NAME##_DataGet /* dummy name - must be unique */
#define OUTPUT_STRUCT_NAME MASTER_FUNCTION_NAME##_data_out /* dummy name - must be unique */
#define DATAOUT_NAME MASTER_FUNCTION_NAME##_DataOut /* dummy name - must be unique */
#define DATARESULT_NAME MASTER_FUNCTION_NAME##_DataResult /* dummy name - must be unique */
#define TAG_NAME_A MASTER_FUNCTION_NAME##_TAG_A /* dummy name - must be unique */
#define TAG_NAME_B MASTER_FUNCTION_NAME##_TAG_B /* dummy name - must be unique */
#define PRIMARY_SUBFUN_NAME MASTER_FUNCTION_NAME##_subfun_primary /* dummy name - must be unique */
#define SECONDARY_SUBFUN_NAME MASTER_FUNCTION_NAME##_subfun_secondary /* dummy name - must be unique */
#define DATA_EXCHANGE_FUNCTION MASTER_FUNCTION_NAME##_exchange_data /* dummy name - must be unique */

static inline void *PRIMARY_SUBFUN_NAME(void *p, int loop_iteration);
static inline void *SECONDARY_SUBFUN_NAME(void *p, int loop_iteration);
void DATA_EXCHANGE_FUNCTION(void);

