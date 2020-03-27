#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioplib/seq.h"
#include "bioplib/sequtil.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

#define MAXBUFF           256
#define HUGEBUFF          1024
#define CHAINTYPE_UNKNOWN 0
#define CHAINTYPE_LIGHT   1
#define CHAINTYPE_HEAVY   2
#define AGLDATADIR        "mydata"

#define CHAINTYPE(x) (                            \
   (x)==CHAINTYPE_LIGHT ? "Light" :               \
    ((x)==CHAINTYPE_HEAVY ? "Heavy" : "Unknown"))

int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, BOOL *verbose, int *chainType);
void ProcessSeq(FILE *out, char *seq, BOOL verbose);
REAL ScanAgainstDB(char *type, char *seq, BOOL verbose, char *match);
REAL CompareSeqs(char *theSeq, char *seq);





int main(int argc, char **argv)
{
   char infile[MAXBUFF+1],
      outfile[MAXBUFF+1];
   BOOL verbose = FALSE;
   int chainType = CHAINTYPE_UNKNOWN;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &verbose, &chainType))
   {
      FILE *in= stdin,
         *out = stdout;
   
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         char header[MAXBUFF+1];
         char *seq = NULL;

         if((seq = blReadFASTA(in, header, MAXBUFF))!=NULL)
         {
            ProcessSeq(out, seq, verbose);
            free(seq);
         }
      }
   }
   else
   {
      Usage();

   }
   
   
   return(0);
   
}


    
/************************************************************************/
void ProcessSeq(FILE *out, char *seq, BOOL verbose)
{
   int chainType = CHAINTYPE_UNKNOWN;
   char lvMatch[MAXBUFF+1],
       hvMatch[MAXBUFF+1],
       lcMatch[MAXBUFF+1],
      hcMatch[MAXBUFF+1];
   REAL lvScore, hvScore, lcScore, hcScore;
   
   lvScore = ScanAgainstDB("light_v", seq, verbose, lvMatch);
   hvScore = ScanAgainstDB("heavy_v", seq, verbose, hvMatch);
   lcScore = ScanAgainstDB("light_c", seq, verbose, lcMatch);
   hcScore = ScanAgainstDB("heavy_c", seq, verbose, hcMatch);

   if((lvScore > hvScore) || (lcScore > hcScore))
   {
      chainType = CHAINTYPE_LIGHT;
   }
   else if((hvScore > lvScore) || (hcScore > lcScore))
   {
      chainType = CHAINTYPE_HEAVY;
   }
   
   fprintf(out, "Chain type: %s\n", CHAINTYPE(chainType));
}

/************************************************************************/
REAL ScanAgainstDB(char *type, char *theSeq, BOOL verbose, char *match)
{
   char filename[MAXBUFF+1];
   BOOL noEnv;
   FILE *dbFp = NULL;
   REAL maxScore = 0.0;
   char bestMatch[MAXBUFF+1];
   
   sprintf(filename, "%s/%s.dat", AGLDATADIR, type);
   
   if(dbFp = blOpenFile(filename, AGLDATADIR, "r", &noEnv))
   {
      char header[MAXBUFF+1];
      char *seq = NULL;

      while((seq = blReadFASTA(dbFp, header, MAXBUFF))!=NULL)
      {
         REAL score;

         if(verbose)
         {
            fprintf(stderr, "Comparing with %s\n", header);
         }
         
         score = CompareSeqs(theSeq, seq);
         if(score > maxScore)
         {
            maxScore = score;
            strncpy(bestMatch, header, MAXBUFF);
         }
         
         free(seq);
      }
      FCLOSE(dbFp);
   }

   strncpy(match, bestMatch, MAXBUFF);
   return(maxScore);
}

/************************************************************************/
REAL CompareSeqs(char *theSeq, char *seq)
{
   int score;
   char align1[HUGEBUFF+1],
      align2[HUGEBUFF+1];
   int alignLen;
   int shortSeqLen = MIN(strlen(theSeq), strlen(seq));
   
   score = blAlign(theSeq, strlen(theSeq),
                   seq, strlen(seq),
                   FALSE, /* verbose */
                   TRUE,  /* identity */
                   1,     /* penalty */
                   align1,
                   align2,
                   &alignLen);
   return ((REAL)score / (REAL)shortSeqLen);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Print usage message.

   17.10.95 Original    By: ACRM
   04.03.15 V1.2
   28.10.15 V1.3
*/
void Usage(void)
{

}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *verbose, int *chainType)
   ---------------------------------------------------------------------
   Input:   int      argc        Argument count
            char     **argv      Argument array
   Output:  char     *infile     Input filename (or blank string)
            char     *outfile    Output filename (or blank string)
            BOOL     *verbose    Verbose
            int      *chainType  Chain type
   Returns: BOOL                 Success

   Parse the command line

   26.03.20 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, int *chainType)
{
   argc--;
   argv++;
   
   infile[0]   = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
            break;
         case 'l':
            if(*chainType != CHAINTYPE_UNKNOWN)
               return(FALSE);
            *chainType = CHAINTYPE_LIGHT;
            break;
         case 'h':
            if(*chainType != CHAINTYPE_UNKNOWN)
               return(FALSE);
            *chainType = CHAINTYPE_HEAVY;
            break;
         default:
            return(FALSE);
            break;
         }
         argc--;
         argv++;
      }
      else
      {
         /* Check that there are only 2-4 arguments left                */
         if(argc > 2)
            return(FALSE);

         if(argc)
         {
            strncpy(infile, argv[0], MAXBUFF);
            argc--; argv++;
            if(argc)
            {
               strncpy(outfile, argv[0], MAXBUFF);
               argc--; argv++;
            }
         }
      }
   }
   
   return(TRUE);
}

