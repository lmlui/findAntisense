/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
  ...

  TODO:
  - Comment 'Debug' code back in, have better flags controlling its verbosity.
  - Inline some functions for performance.  (But doesn't -O3 do that already?)

\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "antisense.h"


const char *setConfig (char **optNames, char **optValues)
{
    // Set defaults.
    config.mode          = 1;
    config.minLen        = 10;
    config.maxGuBp       = 1;
    config.maxMismatches = 1;
    config.mask          = NULL;

    // Save caller config, overriding defaults where config is specified.
    char **np = optNames;
    char **vp = optValues;
    while (*np != NULL)
    {
        assert (*vp != NULL);  // TODO: ret err instead?
        Debug ((stderr, "Setting config: (%s) = (%s)\n", *np, *vp));

        if (strcmp (*np, "mode") == 0)
        {
            int val = atoi (*vp);
            if (val != 1)
                return "Invalid value for \"mode\" config option.  Must be one of: {1}.";
            config.mode = (int) val;
        }
        else if (strcmp (*np, "minLen") == 0)
        {
            int val = atoi (*vp);
            if (val < 1)
                return "Invalid value for \"minLen\" config option.  Must be >= 1.";
            config.minLen = (int) val;
        }
        else if (strcmp (*np, "maxGuBp") == 0)
        {
            int val = atoi (*vp);
            if (val < 0)
                return "Invalid value for \"maxGuBp\" config option.  Must be >= 0.";
            config.maxGuBp = (int) val;
        }
        else if (strcmp (*np, "maxMismatches") == 0)
        {
            int val = atoi (*vp);
            if (val < 0)
                return "Invalid value for \"maxMismatches\" config option.  Must be >= 0.";
            config.maxMismatches = (int) val;
        }
        else if (strcmp (*np, "mask") == 0)
        {
            // TODO: Validate that mask contains correct symbols.
            config.mask = *vp;
        }
        else {
            return "Invalid option name.";
        }

        np++;
        vp++;
    }

    return NULL;  // success
}


DuplexList findDuplexes (char **seqsA, char **seqsB)
{
    InitDuplexList (duplexes);
    for (int a = 0; seqsA[a] != NULL; a++)
    {
        for (int b = 0; seqsB[b] != NULL; b++)
        {
            Debug ((stderr,
                    "===============\n"
                    "Making all duplexes between\nseqA: (%s)\nand\nseqB: (%s)\nmask: (%s)\n",
                    seqsA[a], seqsB[b], config.mask));
            assert (strlen (seqsB[b]) == strlen (config.mask));
            duplexes = concatDuplexLists (
                duplexes, findDuplexesBetweenTwoSeqs (
                    seqsA[a], seqsB[b], a, b));
        }
    }
    Debug ((stderr, "Done in findDuplexes, returning...\n"));
    return duplexes;
}


Duplex *makeDuplex (char *seqA, char *seqB, int seqNumA, int seqNumB,
                    int startA, int startB, int length, char *annotBp)
{
#ifndef NDEBUG
    //fprintf (stderr, "---------------\nBuilding duplex...\n");  // TODO: Comment back in and use Debug macro
    //fprintf (stderr, "startA: %u\nstartB: %u\n", startA, startB);
    //fprintf (stderr, "len: %u\n", length);
    //char tempA[1000], tempB[1000];
    //strncpy (tempA, seqA + startA, length);
    //strncpy (tempB, seqB + startB, length);
    //tempA[length] = '\0';
    //tempB[length] = '\0';
    //fprintf (stderr, "(%s)\n(%s)\n", tempA, tempB);

    assert (seqA != NULL);    assert (seqB != NULL);
    assert (length <= Min (strlen (seqA), strlen (seqB)));
    assert (startA < strlen (seqA));
    assert (startB < strlen (seqB));
#endif

    Duplex *duplex = malloc (sizeof (Duplex));
    assert (duplex != NULL);
    duplex->seqA    = seqA;     duplex->seqB    = seqB;
    duplex->seqNumA = seqNumA;  duplex->seqNumB = seqNumB;
    duplex->startA  = startA;   duplex->startB  = startB;
    duplex->length = length;
    if (annotBp == NULL)
        duplex->annotBp = NULL;
    else
    {
        MakeString (duplex->annotBp, strlen (annotBp));
        strcpy (duplex->annotBp, annotBp);
    }
    return duplex;
}


DuplexList findDuplexesBetweenTwoSeqs (char *seqA, char *seqB,
                                       int seqNumA, int seqNumB)
{
    //Debug ((stderr, "Got as input: (%s) (%s).\n", seqA, seqB));

    InitDuplexList (duplexes);

    /* "Slide" shorter sequence from left to right past longer sequence
       to generate all possible alignments without gaps...
       TODO: diagram!

    
    curDupNum = the duplex number (1st, 2nd, ...) we are generating.
    numDupls - curDupNum = duplexes left (NOT including current one).
    lenTop - 1 = num duplexes containing top overhang on left OR right.
    ...

    We will have curDupNum in (1 .. numDupls) duplexes total.
    */

    // Precompute.
    int lenA     = strlen (seqA);
    int lenB     = strlen (seqB);
    int numDupls = lenA + lenB - 1;

    // Put shorter sequence on top, for sanity.
    char *seqTop, *seqBottom;
    int   lenTop,  lenBottom;
    bool  swapped;
    if (lenA <= lenB)
    {
        swapped = false;
        seqTop = seqA;    seqBottom = seqB;
        lenTop = lenA;    lenBottom = lenB;
    }
    else
    {
        swapped = true;
        seqTop = seqB;    seqBottom = seqA;
        lenTop = lenB;    lenBottom = lenA;
    }

    int curDupNum = 1;

    // Overhang on left, e.g.:
    //     xxxx
    //       yyyyyy
    for (; curDupNum < lenTop; curDupNum++)
    {
        //Debug ((stderr, "~~~~~\nleft overhang duplex # %u\n", curDupNum));

        int startTop    = lenTop - curDupNum;  // offset from top end
        int startBottom = 0;
        int dupLen      = curDupNum;

        Duplex *duplex;
        if (swapped)
            duplex = makeDuplex (seqBottom, seqTop, seqNumA, seqNumB,
                                 startBottom, startTop, dupLen, NULL);
        else
            duplex = makeDuplex (seqTop, seqBottom, seqNumA, seqNumB,
                                 startTop, startBottom, dupLen, NULL);
        DuplexList goodDupls = scoreSubDuplexes (duplex);
        duplexes = concatDuplexLists (duplexes, goodDupls);
        freeDuplex (duplex);
    }

    // Top completely nested in bottom (no overhang), e.g.:
    //     xxxx      or     xxxx     or      xxxx
    //     yyyyyy          yyyyyy          yyyyyy
    for (; numDupls - (curDupNum - 1) - (lenTop - 1) > 0; curDupNum++)
    {
        //Debug ((stderr, "~~~~~\nno overhang duplex # %u\n", curDupNum));

        int startTop    = 0;
        // lenTop - 1 = # of duplexes already made with left overhang
        // offset = curDupNum - (lenTop - 1) - 1 = curDupNum - lenTop
        int startBottom = curDupNum - lenTop;
        int dupLen      = lenTop;  // shorter seq defines duples length

        Duplex *duplex;
        if (swapped)
            duplex = makeDuplex (seqBottom, seqTop, seqNumA, seqNumB,
                                 startBottom, startTop, dupLen, NULL);
        else
            duplex = makeDuplex (seqTop, seqBottom, seqNumA, seqNumB,
                                 startTop, startBottom, dupLen, NULL);
        DuplexList goodDupls = scoreSubDuplexes (duplex);
        duplexes = concatDuplexLists (duplexes, goodDupls);
        freeDuplex (duplex);
    }

    // Top overhang on right, e.g.:
    //         xxxx
    //     yyyyyy
    for (; curDupNum <= numDupls; curDupNum++)
    {
        //Debug ((stderr, "~~~~~\nright overhang duplex # %u\n", curDupNum));

        int startTop    = 0;
        int dupLen      = numDupls - curDupNum + 1;
        // lenBottom - 1 - (dupLen - 1) = lenBottom - dupLen
        int startBottom = lenBottom - dupLen;

        Duplex *duplex;
        if (swapped)
            duplex = makeDuplex (seqBottom, seqTop, seqNumA, seqNumB,
                                 startBottom, startTop, dupLen, NULL);
        else
            duplex = makeDuplex (seqTop, seqBottom, seqNumA, seqNumB,
                                 startTop, startBottom, dupLen, NULL);
        DuplexList goodDupls = scoreSubDuplexes (duplex);
        duplexes = concatDuplexLists (duplexes, goodDupls);
        freeDuplex (duplex);
    }

    return duplexes;
}


Duplex *cloneDuplex (Duplex *duplex)
{
    Duplex *retDupl = malloc (sizeof (Duplex));
    assert (retDupl != NULL);

    retDupl->seqA    = duplex->seqA;     retDupl->seqB    = duplex->seqB;
    retDupl->seqNumA = duplex->seqNumA;  retDupl->seqNumB = duplex->seqNumB;
    retDupl->startA  = duplex->startA;   retDupl->startB  = duplex->startB;
    retDupl->length  = duplex->length;

    if (duplex->annotBp == NULL)
        retDupl->annotBp = NULL;
    else
    {
        // Copy base pair annotation, if it exists.
        MakeString (retDupl->annotBp, strlen (duplex->annotBp));
        assert (retDupl->annotBp != NULL);
        strcpy (retDupl->annotBp, duplex->annotBp);
        assert (strlen (retDupl->annotBp) == strlen (duplex->annotBp));
    }
    return retDupl;
}


DuplexList scoreSubDuplexes (Duplex *duplex)
{
    Debug ((stderr, "In scoreSubDuplexes()...\n"));
    InitDuplexList (subDupList);

    /* Precompute base pair string for efficiency. */
    // TODO: Making the BP string given two seqs should be a lib func.
    Debug ((stderr, "Precomputing base pair annot...\n"));
    // TODO: It would make more sense if the caller did the pre-computing and
    // we just checked it here, if it wasn't done, THEN do it, but DO NOT
    // write the pre-computed result to "duplex".
    // This function really shouldn't modify data passed in.
    // We can pass the "annotBp" string as a separate arg if need be,
    // but don't mess with the caller's data.
    // "Don't mess with caller data unless absolutely necessary" ought to be
    // a rule, and any exception to it must be CLEARLY documented in the .h
    // file!
    if (duplex->annotBp == NULL)
    {
        MakeString (duplex->annotBp, duplex->length);
    }
    else
    {
        assert (strlen (duplex->annotBp) == duplex->length);
    }
    
    for (int i = 0; i < duplex->length; i++)
    {
        switch (bpType (duplex->seqA[duplex->startA + i],
                        duplex->seqB[duplex->startB + i]))
        {
            case BP_WC:
                duplex->annotBp[i] = BP_CHAR_WC;
                break;
            case BP_GU:
                duplex->annotBp[i] = BP_CHAR_GU;
                break;
            case BP_NONE:
                duplex->annotBp[i] = BP_CHAR_NONE;
                break;
            default:
                assert (false);
        }
    }
    Debug ((stderr, "BP annot: (%s)\n", duplex->annotBp));

    /* Generate and score all possible sub-duplexes of length >= minLen.
       ...describe iteration algorithm, since iteration order is CRUCIAL...
       ...0-based half-open coords...
    */

    assert (config.minLen >= 1);  // code below may depend on this property
   
    // Lower bound (open) on j.
    // Used to restrict the part of the sub-duplex matrix containing
    // sub-duplexes shorter than a good match, so that we return only the
    // maximal match.
    // TODO: Have an option disabling setting 'leftBound', effectively
    // disabling it (it don't do nothing when set to 0), in case user wants
    // ALL duplexes.
    int leftBound = 0;  // init to no-op value

    for (int i = 0; i <= (duplex->length - config.minLen); i++)
    {
        for (int j = duplex->length;
             (j >= i + config.minLen) && (j > leftBound);
             j--)
        {
            //printf ("i = %d  j = %d  leftBound = %d  minLen = %d\n",
            //        i, j, leftBound, config.minLen);  //d!!!

            if (config.mask != NULL)
            {
                // Ensure the sub-duplex is in compliance with the mask.
                // No need to waste time scoring it if it doesn't fit the
                // strict criteria imposed by the mask.
                // This is a check that is independent of the mode.
                // NB: Chars in mask and seqB are in 1-to-1 correspondence.

                // Check mask to make sure there are no required base pair
                // annots outside of this duplex.  If there are, skip
                // immediately, since this sub-duplex cannot contain sites
                // that are required to be base paired -- they are outside
                // of the sub-duplex.
                bool skip = false;

                for (int a = 0; a < i; a++)
                {
                    if (config.mask[duplex->startB + a] == '!' ||
                        config.mask[duplex->startB + a] == '?')
                    {
                        skip = true;
                        break;
                    }
                }
                if (skip) break;  // goto next duplex

                for (int b = j; b < duplex->length; b++)
                {
                    if (config.mask[duplex->startB + b] == '!' ||
                        config.mask[duplex->startB + b] == '?')
                    {
                        skip = true;
                        break;
                    }
                }
                if (skip) break;  // goto next duplex

                // Check to make sure the duplex itself satisfies the mask.
                for (int c = i; c < j; c++)
                {
                    if (config.mask[duplex->startB + c] == '!')
                    {
                        if (duplex->annotBp[c] != BP_CHAR_WC)
                        {
                            // Must have a W/C bp at this position, but we don't.
                            // Instantly reject.
                            skip = true;
                            break;
                        }
                    }
                    else if (config.mask[duplex->startB + c] == '?')
                    {
                        if (duplex->annotBp[c] == BP_CHAR_NONE)
                        {
                            // Must have a canonical bp at this position, but we don't.
                            // Instantly reject.
                            skip = true;
                            break;
                        }
                    }
                    else
                    {
                        assert (config.mask[duplex->startB + c] == '.');
                    }
                }
                if (skip) break;  // goto next duplex
            }

            Duplex *subDupl = cloneDuplex (duplex);
            subDupl->startA += i;
            subDupl->startB += i;
            subDupl->length  = j - i;

            free (subDupl->annotBp);  // wipe cloned full-duplex annot
            MakeString (subDupl->annotBp, subDupl->length);
            strncpy (subDupl->annotBp, duplex->annotBp + i, subDupl->length);

            switch (config.mode)
            {
                case 1:
                    scoreDuplexSimple (subDupl);
                    if (subDupl->score > 0.0)
                    {
                        // Duplex passed threshold.
                        // Save it.
                        appendToDuplexList (&subDupList, subDupl);

                        // Mark off part of the grid that contains
                        // sub-duplexes, so that we only keep the maximal
                        // match.
                        leftBound = j;
                    }
                    else
                    {
                        // Sub-duplex did not pass threshold; get rid of it.
                        free (subDupl);
                    }
                    break;
                default:
                    assert (false);  // unreachable
            }
        }
    }

    Debug ((stderr, "Generated and scored %u sub-duplexes.\n",
            subDupList.length));

    return subDupList;
}


Duplex *scoreDuplexSimple (Duplex *duplex)
{
#ifndef NDEBUG
    Debug ((stderr, "In scoreDuplexSimple(), considering duplex:\n"));
    char *subseqA, *subseqB;
    MakeString (subseqA, duplex->length);
    MakeString (subseqB, duplex->length);
    strncpy (subseqA, duplex->seqA + duplex->startA, duplex->length);
    strncpy (subseqB, duplex->seqB + duplex->startB, duplex->length);
    Debug ((stderr, "  %s\n  %s\n  %s\n", subseqA, duplex->annotBp, subseqB));
    Debug ((stderr, "Seq A start %d\nSeq B start %d\n", duplex->startA, duplex->startB));
#endif

    if (duplex->annotBp[0] == BP_CHAR_NONE
            || duplex->annotBp[duplex->length - 1] == BP_CHAR_NONE)
    {
        // First or last position in duplex is not a base pair;
        // reject immediately, because we only want duplexes that begin and
        // end with a canonical BP.
        duplex->score = -1.0;
    }
    else
    {
        // Determine whether duplex satisfies the base pairing criteria.
        int numGuBp       = 0;
        int numMismatches = 0;
        
        for (int i = 0; i < duplex->length; i++)
        {
            if (duplex->annotBp[i] == BP_CHAR_GU)
                numGuBp++;
            else if (duplex->annotBp[i] == BP_CHAR_NONE)
                numMismatches++;
        }

        if (numGuBp <= config.maxGuBp
                && numMismatches <= config.maxMismatches)
            duplex->score = 1.0;   // accept
        else
            duplex->score = -1.0;  // reject
    }
    Debug ((stderr, "Got score: %f\n", duplex->score));
    return duplex;
}


void appendToDuplexList (DuplexList *list, Duplex *duplex)
{
    assert (list   != NULL);
    assert (duplex != NULL);

    DuplexListItem *item = malloc (sizeof (DuplexListItem));
    assert (item != NULL);
    item->duplex = duplex;
    item->next = NULL;

    if (list->length == 0)
        list->head = item;
    else
    {
        assert (list->tail != NULL);
        assert (list->tail->next == NULL);
        list->tail->next = item;
    }
    list->tail = item;
    list->length++;
}


DuplexList concatDuplexLists (DuplexList list1, DuplexList list2)
{
    DuplexList listRet;

    /* Catch badly-formed lists. */
#ifndef NDEBUG
    if (list1.length == 0)
    {
        assert (list1.head == NULL);
        assert (list1.tail == NULL);
    }
    else
    {
        assert (list1.head != NULL);
        assert (list1.tail != NULL);
        assert ((list1.tail)->next == NULL);
    }
    if (list2.length == 0)
    {
        assert (list2.head == NULL);
        assert (list2.tail == NULL);        
    }
    else
    {
        assert (list2.head != NULL);
        assert (list2.tail != NULL);
        assert ((list2.tail)->next == NULL);
    }
    //fprintf (stderr, "Concatenating duplex lists of lengths %u and %u...\n",
    //         list1.length, list2.length);  // TODO: Comment back in, use Debug macro
#endif

    /* Concatenate. */
    if (list1.length == 0)
    {
        if (list2.length == 0)
        {
            listRet.head   = NULL;
            listRet.tail   = NULL;
            listRet.length = 0;
        }
        else
        {
            listRet.head   = list2.head;
            listRet.tail   = list2.tail;
            listRet.length = list2.length;
        }
    }
    else
    {
        if (list2.length == 0)
        {
            listRet.head   = list1.head;
            listRet.tail   = list1.tail;
            listRet.length = list1.length;
        }
        else
        {
            (list1.tail)->next = list2.head;

            listRet.head   = list1.head;
            listRet.tail   = list2.tail;
            listRet.length = list1.length + list2.length;
        }
    }
    
    return listRet;
}


void freeDuplexList (DuplexList duplexes)
{
    DuplexListItem *iter = duplexes.head;
    int duplsFreed = 0;
    while (iter != NULL)
    {
        freeDuplex (iter->duplex);
        DuplexListItem *next = iter->next;
        free (iter);
        iter = next;
        duplsFreed++;
    }
    assert (duplsFreed == duplexes.length);
}


void freeDuplex (Duplex *duplex)
{
    if (duplex->annotBp != NULL)
        free (duplex->annotBp);
    free (duplex);
}


