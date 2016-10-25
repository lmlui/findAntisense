#file c_antisense.pyx: Cython source file

cimport cantisense
from libc.stdlib cimport malloc, free

def setConfig(optNames,optValues): #gets passed two lists. Convert to C inpus so can call C function
    cdef char **names_buf = <char**> malloc((len(optNames)+1) * sizeof(char*)) # +1 for terminator 
    cdef char **val_buf = <char**> malloc((len(optValues) + 1)*sizeof(char*))  
    for i in xrange(len(optNames)): #The length of optNames and optVales must be equal, as they're in 1-to-1 correspondence
            names_buf[i] = <char*> optNames[i]
            val_buf[i] = <char*> optValues[i]
    names_buf[len(optNames)] = NULL #Add terminators
    val_buf[len(optValues)] = NULL
    result = cantisense.setConfig(names_buf,val_buf)
    free(names_buf) 
    free(val_buf)
    if result == NULL: return None
    else: return <bytes> result

def findDuplexes(seqsA, seqsB): #takes 2 python lists of strings
    cdef char **seqsA_buf = <char**> malloc((len(seqsA) + 1) * sizeof(char*)) # +1 for terminator 
    cdef char **seqsB_buf = <char**> malloc((len(seqsB) + 1)*sizeof(char*))  
    for i in xrange(len(seqsA)): 
            seqsA_buf[i] = <char*> seqsA[i]
    for i in xrange(len(seqsB)):
            seqsB_buf[i] = <char*> seqsB[i]
    seqsA_buf[len(seqsA)] = NULL #Add terminators
    seqsB_buf[len(seqsB)] = NULL
    duplex_list = cantisense.findDuplexes(seqsA_buf,seqsB_buf)

    #Go from DuplexList to a list of python tuples 
    #Alternatively make Python Duple Objects directly to fix AUV's code
    py_duplexes = []
    it = duplex_list.head   #Cython automatically dereferences 
    while it != NULL:
        c_duplex = it.duplex
        duplex_t = (c_duplex.seqNumA,c_duplex.seqNumB,c_duplex.startA,c_duplex.startB,c_duplex.length,c_duplex.annotBp,c_duplex.score)
        py_duplexes.append(duplex_t)
        it = it.next

    cantisense.freeDuplexList(duplex_list)
    return py_duplexes
    

#end of file
