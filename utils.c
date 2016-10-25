/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
  Miscellaneous utils.
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include "utils.h"


/*
  ASCII codes (char int values):
    A = 65
    C = 67
    G = 71
    T = 84
    U = 85
  TODO: No safety checks to make sure nucleotides are valid.
*/
int bpType (char bp1, char bp2)
{
    switch (Abs (ToUppercase (bp1) - ToUppercase (bp2)))
    {
        case 4:   // G-C = 4
        case 19:  // T-A = 19
        case 20:  // U-A = 20
	    return BP_WC;
	    break;
        case 13:  // T-G = 13
        case 14:  // U-G = 14
	    return BP_GU;
	    break;
	default:
	    //Debug ((stderr,
		//    "Non-canon base pair: '%c' and '%c' (abs diff is %d).\n",
		//    bp1, bp2, Abs (ToUppercase (bp1) - ToUppercase (bp2))));
	    return BP_NONE;
    }
}
