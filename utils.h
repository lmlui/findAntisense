/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
  Miscellaneous utils.
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>


/*---------------------------------------------------------------------------*\
  Defines and enums.
\*---------------------------------------------------------------------------*/

// Enable this to disable asserts and debug code.
#define NDEBUG

// Types of base pairs.
enum {BP_NONE, BP_GU, BP_WC};
#define BP_CHAR_NONE '.'
#define BP_CHAR_GU   ':'
#define BP_CHAR_WC   '|'


/*---------------------------------------------------------------------------*\
  Macros for functions.
\*---------------------------------------------------------------------------*/

#define Abs(a) ((a) < 0 ? -(a) : (a))

#define Max(a,b) ((a) >= (b) ? (a) : (b))
#define Min(a,b) ((a) <= (b) ? (a) : (b))

/* Convert lowercase char to uppercase.
   If char is not lowercase, do nothing.  TODO: assert it is uppercase
   TODO: type safety?  how will this behave on non-char types?
*/
#define ToUppercase(c) (((int)(c) >= 97 && (int)(c) <= 122) \
			? (char)((int)(c) - 32) : (c))

/* Dynamically allocate space for a string 'numChars' long, including space
   for terminating character.
   'strPtr' must be already declared.
   Also, null-terminates string, just in case.
   This is a macro instead of a function so that assert() will give us the
   native line number.
*/
#define MakeString(strPtr,numChars)					\
    strPtr = malloc (sizeof (char) * (numChars + 1));			\
    if (strPtr == NULL) {						\
	fprintf (stderr,						\
		 "Cannot allocate memory for string of %u characters.", \
		 (unsigned int)numChars);				\
	assert (false);							\
    }									\
    else								\
	strPtr[numChars] = '\0';

/* A sneaky way to enable/disable calls to debug code.
   VERY IMPORTANT to invoke exactly like this, with TWO parentheses:
       Debug ((stderr, "message %s", msgStr));
   or
       Debug ((stdout, "message %s", msgStr));
*/
#ifdef NDEBUG
#define Debug(args)  /* no debugging messages */
#else
#define Debug(args) fprintf args;
#endif


/*---------------------------------------------------------------------------*\
  Types.
\*---------------------------------------------------------------------------*/

// Use this for any non-negative type: lengths, sizes, array indices, etc.
//
// TODO: This may cause more trouble than it's worth (see Foolin/signedVsUnsigned.c)
// Go through all code using this and ASK CAREFULLLY -- do you REALLY need the
// 2x more space that having an unsigned type buys you?  Would your overflow
// be more/less better on signed versus unsigned?
//
// TODO: How is this different from size_t?  Is this really necessary?
typedef unsigned int len_t;


/*---------------------------------------------------------------------------*\
  Functions.
\*---------------------------------------------------------------------------*/

/* Get the type of the base pair between 'bp1' and 'bp2'.
*/
int bpType (char bp1, char bp2);


#endif
