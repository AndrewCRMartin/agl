#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioplib/sequtil.h"
#include "bioplib/general.h"

#define MAXBUFF           256
#define CHAINTYPE_UNKNOWN 0
#define CHAINTYPE_LIGHT   1
#define CHAINTYPE_HEAVY   2

#define CHAINTYPE(x) ( \
   (x)==CHAINTYPE_LIGHT ? "Light" : \
    ((x)==CHAINTYPE_HEAVY ? "Heavy" : "Unknown"))

int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void ProcessSeq(FILE *out, char *seq);
REAL ScanAgainstDB(char *type, char *seq, char *match);




int main(int argc, char **argv)
{
   char infile[MAXBUFF+1],
      outfile[MAXBUFF+1];
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      FILE *in= stdin,
         *out = stdout;
   
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         char header[MAXBUFF+1];
         char *seq = NULL;

         if((seq = blReadFASTA(in, header, MAXBUFF))!=NULL)
         {
            ProcessSeq(out, seq);
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
void ProcessSeq(FILE *out, char *seq)
{
   int chainType = CHAINTYPE_UNKNOWN;
   char lvMatch[MAXBUFF+1],
       hvMatch[MAXBUFF+1],
       lcMatch[MAXBUFF+1],
      hcMatch[MAXBUFF+1];
   REAL lvScore, hvScore, lcScore, hcScore;
   
   lvScore = ScanAgainstDB("light_v", seq, lvMatch);
   hvScore = ScanAgainstDB("heavy_v", seq, hvMatch);
   lcScore = ScanAgainstDB("light_c", seq, lcMatch);
   hcScore = ScanAgainstDB("heavy_c", seq, hcMatch);

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
REAL ScanAgainstDB(char *type, char *seq, char *match)
{
   if(!strncmp(type, "light_v",7))
   {
      return(100.0);
   }
   else if(!strncmp(type, "heavy_v",7))
   {
      return(80.0);
   }
   else if(!strncmp(type, "light_c",7))
   {
      return(20.0);
   }
   else if(!strncmp(type, "heavy_c",7))
   {
      return(30.0);
   }
   return(0.0);
   
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int      argc        Argument count
            char     **argv      Argument array
   Output:  char     *infile     Input filename (or blank string)
            char     *outfile    Output filename (or blank string)
   Returns: BOOL                 Success

   Parse the command line

   26.03.20 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
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

