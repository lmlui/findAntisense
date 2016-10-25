#file: cantisense.pxd

cdef extern from "antisense.h": 
    #Structs go here
    ctypedef struct Duplex:
        char * seqA
        char * seqB
        int seqNumA
        int seqNumB
        int startA
        int startB
        int length
        char * annotBp
        double score
    
    ctypedef _DuplexListItem DuplexListItem
    cdef struct _DuplexListItem:
        Duplex * duplex
        DuplexListItem * next 
    

    ctypedef struct DuplexList:
        DuplexListItem * head
        DuplexListItem * tail
        int length
    
    cdef struct Config:
        int mode
        int minLen
        int maxGuBp
        int maxMismatches
        char * mask

    ctypedef Config config


    #functions
    const char *setConfig (char **optNames, char **optValues)
    DuplexList findDuplexes (char **seqsA, char **seqsB)
    Duplex *makeDuplex (char *seqA, char *seqB, int seqNumA, int seqNumB, int startA, int startB, int length, char *annotBp)
    DuplexList findDuplexesBetweenTwoSeqs (char *seqA, char *seqB, int seqNumA, int seqNumB)
    Duplex *cloneDuplex (Duplex *duplex)
    DuplexList scoreSubDuplexes (Duplex *duplex)
    Duplex *scoreDuplexSimple (Duplex *duplex)
    void appendToDuplexList (DuplexList *list, Duplex *duplex)
    DuplexList concatDuplexLists (DuplexList list1, DuplexList list2)
    void freeDuplexList (DuplexList duplexes)
    void freeDuplex (Duplex *duplex)



#end of file: cantisense.pxd
