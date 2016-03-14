/* NINJA OTU Picker: NINJA Is Not Just Another OTU Picker -- preparation program
   Knights Lab (www.knightslab.org/ninja)
   This program generates the databases necessary for OTU mapping with bowtie2.
   
   Compilation information (GCC):
   Ascribes to std=gnu and doesn't require C99 support or UNIX intrinsics.
   Flags: -m64 -Ofast -flto ninja_prep.c -fprofile-generate [-fprofile-use]
*/
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

inline int xcmp(str1, str2) register const char *str1, *str2; {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return (*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}

int cmp(const void *v1, const void *v2) {
	const char *i1 = **(const char ***)v1;
	const char *i2 = **(const char ***)v2;
	return xcmp(i1,i2);
}

int main( int argc, char *argv[] )
{
	const int V_97[20] = {0,0, 0,0, 0,0, 1901,2231, 2263,4051, 0,0, 0,0, 0,0, 0,0, 0,0};
						//null, V1,  V2,  V3,         V4,      V5,   V6,  V7,  V8,  V9
	FILE *in_seqs, *out_seqs, *out_DB, *out_map;
	char *inmode = "rb", *outmode = "wb"; 
	int doSort = 0, V = 0, argS = 1; 
	if ( argc < 4 || argc > 6 )
    {
		printf( "\nNINJA Is Not Just Another OTU Picker: preparation program. Usage:\n");
		printf( "ninja_prep [V4] [sort] in_refs.fa out_DB.fa out_map.db\n" );
		printf( "\nCONTROL PARAMETERS (optional):\n");
		printf( "V[X]: if specified, produce V[X] database only (X = 3 or 4)\n");
		printf( "  Note: V3 and V4 region extraction assumes gg_13_8 97_otus input\n");
		printf( "sort: if specified, database is rearranged lexicographically (faster)\n");
		printf( "\nINPUT PARAMETERS:\n");
		printf( "in_refs.fa: the (aligned) references in fasta format\n");
		printf( "\n" "OUTPUT PARAMETERS:\n");
		printf( "out_DBFile.fa: database file for this program to generate\n");
		printf( "out_mapFile.db: mapping file for this program to generate\n");
		return 1;
    }
	if (argc > 4) {
		if (!strcmp(argv[1],"V4") || !strcmp(argv[2],"V4")) { V = 4; ++argS; }
		if (!strcmp(argv[1],"V3") || !strcmp(argv[2],"V3")) { V = 3; ++argS; }
		if (!strcmp(argv[1],"sort") || !strcmp(argv[2],"sort")) { doSort = 1; ++argS; }
	}
	
	in_seqs = fopen(argv[argS++], inmode);
	out_DB = fopen(argv[argS++], outmode);
	out_map = fopen(argv[argS++], outmode);
	
	if (!in_seqs || !out_DB || !out_map) { //|| !out_seqs
	  fprintf(stderr, "Can't open input/output file(s)! Check parameter order.\n");
	  return 1;
	}
	
	clock_t start; double cpu_time_used; start = clock();
	
	// Read in first OTU's sequence (V4 region can be used here), get max length
	unsigned int otu;
	char line[250000]; // assumption: no one sequence will be over 250k in length
	
	//fscanf(in_seqs, ">%d\n%s\n", &otu, line); // grab first line 
	//int line_len = strlen(line);
	//printf("Line length is: %d (otu here: %d)\n", line_len, otu);
	int lines = 0, dashes = 0;
	
	/// Find file length, then number of lines in it
	fseek(in_seqs, 0, SEEK_END);
	long fsize = ftell(in_seqs);
	fseek(in_seqs, 0, SEEK_SET); 

	char *string = malloc(fsize + 1), *stringP = string - sizeof(char); // allocate fsize direct bytes
	if (!string) { printf("Insufficient memory for read-in.\n"); return 1; }
	fread(string, fsize, 1, in_seqs); //read into string: elems of fsize bytes, 1 elem, using ifp pointer
	fclose(in_seqs); // close the file
	
	// Sequence pre-pass
	int inSeq = 0, seqChars = 0;
	long sz = fsize; do {
		if (*++stringP == '\n') { ++lines; inSeq ^= 1;} 
		else if (inSeq) { if (*stringP == '-') ++dashes; else ++seqChars; }
	} while (--sz);
	int numEntries = lines/2, nE;
	printf("Number of refs: %d\n", numEntries);
	
	char *seqBuf = malloc(10000), *seqBufP = seqBuf,
		**strArr = malloc(numEntries * sizeof(char *)), **strArrP = strArr;
	unsigned *DBid = malloc(numEntries * sizeof (unsigned)), *idPtr = DBid,
		*DBix = malloc(numEntries * sizeof (unsigned)), *ixPtr = DBix;
	if (!seqBuf || !strArr || !DBid || !DBix)
		{ printf("Cannot allocate pre-run memory.\n"); return 1; }
	unsigned curCharNum = 0, bufCharNum = 0;
	// loop thru the file again and put things in their place
	inSeq = 0; seqChars = 0;
	char atoiBuf[64]; 
	char *atoiBufPtr = atoiBuf;
	stringP = string - sizeof(char);
	int V_count, V_beg = V_97[V*2], rng;
	int V_rng = V_97[V*2+1] - V_beg;
	//printf("V_beg %d, V_end %d, range is %d\n", V_beg, V_97[V*2 + 1], V_rng);
	while (*++stringP) { 
		if (*stringP == '\n') { 
			if (!inSeq) { // we're finishing the int bit and prepping for character line
				memset(atoiBufPtr, '\0', 1); // add inset null prior to integer conversion
				*idPtr++ = atoi(atoiBuf); 
				atoiBufPtr = atoiBuf;
				// reset the string buffer pointer
				seqBufP = seqBuf;
				bufCharNum = 0;
				rng = 0;
				// fast-forward to the correct character in the stream
				if (V) { V_count = V_beg; do ++stringP; while (--V_count); }
				//printf("\n");
			}
			else {
				*ixPtr++ = bufCharNum;
				memset(seqBufP,'\0',1);
				*strArrP = malloc(bufCharNum+1); 
#ifdef DEBUG
				if (!*strArrP) { printf("pre-sequence memory depleted!\n"); return 1; }
#endif
				
				strcpy(*strArrP, seqBuf);
				//*(*strArrP+++bufCharNum)=0; // LR >*+++ notation
				memset(*strArrP+++bufCharNum,'\0',1); 
			}
			inSeq ^= 1;  // swap write target
		}
		else if (*stringP == '>') continue;
		else if (*stringP == '-') { if (V) ++rng; }
		else if (inSeq) { 
			*seqBufP++ = *stringP;
			//printf("%c,%u ",*stringP,rng+bufCharNum);
			//printf("%c",*stringP);
			++curCharNum; ++bufCharNum;
			// fast-forward to end if V_bound reached
			if (V && (rng + bufCharNum >= V_rng)) 
				{ while (*++stringP != '\n'); --stringP; continue; }
		}
		else *atoiBufPtr++ = *stringP;
	}
	free(string); free(seqBuf);
	printf("Database generated.\n");

	// Write database to output file
	fprintf(out_DB, ">NINJAdb\n");
#ifdef PROFILE
	printf("->Database generation: %f\n\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif

#ifdef PRINTFASTA
	FILE *out_fa = fopen("out_db_fa.fa", outmode);
	strArrP = strArr - 1; idPtr = DBid - 1; ixPtr = DBix - 1;
	nE = numEntries; do {
		fprintf(out_fa, ">%u\n%s\n",*++idPtr,*++strArrP);
	} while (--nE);
#endif

	if (!doSort) {
		int total = 0-*DBix;
		strArrP = strArr - 1; idPtr = DBid - 1; ixPtr = DBix - 1;
		nE = numEntries; do {
			fprintf(out_DB, "%s", *++strArrP);
			fprintf(out_map, "%u,%u\n",total +=*++ixPtr, *++idPtr);
		} while (--nE);
	}
	else {
		// sort the sequences from an array of their pointers
		printf("Reference DB will be sorted.\n");
		char ***parray = malloc(numEntries * sizeof(char **)), ***pp = parray;
		if (!parray) { printf("Out of post-memory.\n"); return 1; }
		strArrP = strArr; 
		nE = numEntries; do *pp++ = strArrP++; while (--nE); 
		qsort(parray,numEntries, sizeof *parray,cmp);
		
		// add a manual de-duping step here? 13% of the reference V4 DB are duplicates
		
		// rearrange sequences themselves to match their new order in parray
		unsigned *OtuS = malloc(numEntries * sizeof(unsigned int)), *OtuSp = OtuS,
			*newIX = malloc(numEntries * sizeof (unsigned int)), *nip = newIX;
		if (!OtuS || !newIX) { printf("Out of post-memory.\n"); return 1; }
		unsigned totalIx = 0;
		pp = parray; nE = numEntries; do {
			*OtuSp++ = *(DBid + (*pp - strArr)); 
			*nip++ = totalIx;
			totalIx += *(DBix + (*pp - strArr));
		} while (--nE); 
		
		pp = parray; nip = newIX; OtuSp = OtuS; 
		nE = numEntries; do {
			fprintf(out_DB, "%s", **pp++);
			fprintf(out_map, "%d,%d\n",*nip++, *OtuSp++);
		} while (--nE);
		
		free(OtuS); free(newIX); free(parray); // clean-up
	}
	
	// Free currently used memory
	fclose(out_DB); fclose(out_map);
	strArrP = strArr - 1; nE = numEntries; do free(*++strArrP); while (--nE);
	free(strArr); free(DBid); free(DBix); 
	
#ifdef PROFILE
	printf("->All ref processing: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif

	return 0;
}
