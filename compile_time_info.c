#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        USE_FFTW3\n"
"        ALLOW_IMBALANCED_GASPARTICLELOAD\n"
"        GALSF\n"
"        ADAPTIVE_GRAVSOFT_FORGAS\n"
"        ADM\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        GALSF_SFR_VIRIAL_SF_CRITERION=1\n"
"        GALSF_FB_MECHANICAL\n"
"        METALS\n"
"        TURB_DIFF_METALS\n"
"        TURB_DIFF_METALS_LOWORDER\n"
"        COOLING\n"
"        COOL_METAL_LINES_BY_SPECIES\n"
"        COOL_LOW_TEMPERATURES\n"
"        MULTIPLEDOMAINS=64\n"
"\n");
}
