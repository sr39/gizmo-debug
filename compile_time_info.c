#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        BOX_PERIODIC\n"
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        PMGRID=128\n"
"        MULTIPLEDOMAINS=64\n"
"        PM_PLACEHIGHRESOLUTION=3\n"
"        ADAPTIVE_GRAVSOFT_FORALL=100\n"
"        DM_BARYON_INTERACTION\n"
"        HAVE_HDF5\n"
"        DEBUG\n"
"\n");
}
