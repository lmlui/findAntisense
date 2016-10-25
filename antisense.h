/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
  ...
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef ANTISENSE_H
#define ANTISENSE_H


// Enable this to disable asserts and debug code.
// Must go before utils.h to have effect on Debug macro included from it.
//#define NDEBUG

// Need type definitions from here.
#include "utils.h"

// TODO: Validate that it conforms to http://c-faq.com/cpp/index.html
// especially Q10.4.
#define InitDuplexList(name)  \
    DuplexList name;          \
    name.head   = NULL;       \
    name.tail   = NULL;       \
    name.length = 0;


/*---------------------------------------------------------------------------*\
   Types.
\*---------------------------------------------------------------------------*/

// A simple duplex model.
// We assume both sequences in the duplex are of the same length.
typedef struct {
    char  *seqA;     // Pointer to parent sequence nucleotides.  Parent seqs
    char  *seqB;     // MUST NOT BE MODIFIED, since they may be shared among
                     // several duplexes or exist in Python data structures.
    int    seqNumA;  // Numerical ID of parent sequence.
    int    seqNumB;
    int    startA;   // Start of duplex within the respective sequences.
    int    startB;
    int    length;   // Length of duplex.
    char  *annotBp;  // Base pairing annotation (optional).
    double score;    // Score (optional).
} Duplex;

// ...
typedef struct _DuplexListItem DuplexListItem;
struct _DuplexListItem {
    Duplex         *duplex;
    DuplexListItem *next;
};

// ...
typedef struct {
    DuplexListItem *head;
    DuplexListItem *tail;
    int             length;
} DuplexList;

// Config settings that specify how the algorithm(s) should work.
// Made global for efficiency (no need to pass config to repeatedly called
// functions).
//
// TODO: How to set defaults in case setConfig() isn't called?
//
// 'mask' gets applied to seqs from 'seqsB' in findDuplexes() and etc, so
// must be of same length to get a 1-to-1 mapping between mask and seq.
// It must be made up of only the following chars that denote constraints
// at each position:
//     ! = only W/C base pairing allowed
//     ? = any canonical base pairing (W/C or G/U) allowed
//     . = no constraints (anything allowed)
// TODO: Mask chars should be defined in utils.h
//
// The mask is applied REGARDLESS of mode.  We check whether a duplex
// satisfies the mask FIRST, before we even try to score the duplex, for
// efficiency, since we reject duplexes that are invalid according to the
// mask earlier than later.  Think of 'mask' as a filter that overrides
// all others to discard non-compliant duplexes.
//
// TODO: Only one mask can be specified, which means all seqs in 'seqsB'
// must be of same len and share the same mask!  Works for now, but very
// inconvenient for all-vs-all comparison of heterogeneous sets of seqs.
static struct Config {
//    int    mode = 1;  // ARGH why can't I do this?  (TODO)
    int    mode;
    int    minLen;
    int    maxGuBp;
    int    maxMismatches;
    char  *mask;
} config;


/*---------------------------------------------------------------------------*\
   Functions.
\*---------------------------------------------------------------------------*/

/* Set up the global config ...

   You MUST call this function BEFORE using anything else in the module,
   because this function initializes the default config!  (And/or allows you
   to specify your own config.)  If you don't do this, things will CRASH
   when you run the duplex search!

   ...explain optNames, optValues...

  On success, returns a null pointer.  On failure, returns a constant string
  describing the error.

  TODO: This should use variable-length args... somehow.
  Maybe invoked like this, alternating names and values:
      char *err = setConfig ("mode", 0, "minLen", 8);

  TODO: Currenly has iffiness that if you pass an invalid option and an error
  is returned, the correct options before it get set/saved to global config.

  TODO: Errors could be more informative, like include values of botched
  options, but this would be hard.  Can't malloc (will leak mem), but can't
  keep dynamically created strings on the stack, since later calls will change
  those vals and older pointers will be meaningless.
 */
const char *setConfig (char **optNames, char **optValues);

/* ...run the algorithm...
  Call freeDuplexList() on the returned pointer to completely clean up memory.

  seqsA and seqsB are arrays (containing pointers to sequence strings).  They
  must use a NULL pointer as the last item to denote end of the list.
  An array can be empty, in which case the argument is a NULL pointer.
 */
DuplexList findDuplexes (char **seqsA, char **seqsB);

/* ... 
 */
Duplex *makeDuplex (char *seqA, char *seqB, int seqNumA, int seqNumB,
                    int startA, int startB, int length, char *annotBp);

/* ...

  Return modes:
    (1) Return list of all duplexes above threshold (TODO: sort order?).
    (2) Return only maximal duplexes above threshold.
 */
DuplexList findDuplexesBetweenTwoSeqs (char *seqA, char *seqB,
                                       int seqNumA, int seqNumB);

/* Clone a duplex by mallocing a new one and copying all data members.
   String member contents are cloned also, except for seqA and seqB, where
   only the pointers to the strings are copied, since we want to maintain
   shared pointing to a parent sequence among many duplexes.
 */
Duplex *cloneDuplex (Duplex *duplex);

/* ... 

   The scoring function used depends on the value of the global "struct.mode".
 */
DuplexList scoreSubDuplexes (Duplex *duplex);

/* A very simple scoring function that assigns a positive score if the
   duplex satisfies the heuristic criteria in the global "config" or a
   negative score if it doesn't.  The score is written to the "score" field of
   the Duplex struct.
   
   Zero scores are never assigned, so a test by caller for score > 0 versus
   score < 0 is sufficient to determine whether this function accepted or
   rejected this duplex.
   
   Returns the duplex pointer passed in, for convenience.
 */
Duplex *scoreDuplexSimple (Duplex *duplex);

/* ... 
 */
void appendToDuplexList (DuplexList *list, Duplex *duplex);

/* Concatenate two linked lists by making the tail of the first list point
   to the head of the second list.

   Returns DuplexList wrapper struct for the result.
 */
DuplexList concatDuplexLists (DuplexList list1, DuplexList list2);

/* Deep free memory allocated to a duplex linked list.  All DuplexListItem
   structs are freed, and also the Duplex structs that they wrap.  DuplexList
   itself is assumed to be on the stack and is therefore not touched.
 */
void freeDuplexList (DuplexList duplexes);

/* Free memory used by the Duplex struct.  All struct members that are
   pointers are assumed to be malloced and are therefore also freed,
   except for seqA and seqB, which cannot be touched.
 */
void freeDuplex (Duplex *duplex);


#endif

