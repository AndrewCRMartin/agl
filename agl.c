#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioplib/seq.h"
#include "bioplib/sequtil.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "findfields.h"

#define SMALLBUFF         80
#define MAXBUFF           256
#define HUGEBUFF          1024
#define CHAINTYPE_UNKNOWN 0
#define CHAINTYPE_LIGHT   1
#define CHAINTYPE_HEAVY   2
#define AGLDATADIR        "mydata"
#define THRESHOLD_LV      0.5
#define THRESHOLD_LJ      0.5
#define THRESHOLD_LC      0.5
#define THRESHOLD_HV      0.5
#define THRESHOLD_HJ      0.5
#define THRESHOLD_HC      0.5

#define CHAINTYPE(x) (                            \
   (x)==CHAINTYPE_LIGHT ? "Light" :               \
    ((x)==CHAINTYPE_HEAVY ? "Heavy" : "Unknown"))

int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, BOOL *verbose, int *chainType, char *species);
void ProcessSeq(FILE *out, char *seq, BOOL verbose, int chainType, char *species);
REAL ScanAgainstDB(char *type, char *seq, BOOL verbose, char *species, char *match);
REAL CompareSeqs(char *theSeq, char *seq);
BOOL PreferHeader(char *newHeader, char *oldHeader);
void GetDomainID(char *header, char *id, int maxbuff);


int main(int argc, char **argv)
{
   char infile[MAXBUFF+1],
      outfile[MAXBUFF+1],
      species[MAXBUFF+1];
   BOOL verbose = FALSE;
   int chainType = CHAINTYPE_UNKNOWN;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &verbose, &chainType, species))
   {
      FILE *in  = stdin,
           *out = stdout;

      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         char header[MAXBUFF+1];
         char *seq = NULL;

         if((seq = blReadFASTA(in, header, MAXBUFF))!=NULL)
         {
            ProcessSeq(out, seq, verbose, chainType, species);
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
void ProcessSeq(FILE *out, char *seq, BOOL verbose, int chainType, char *species)
{
   char lvMatch[MAXBUFF+1],
        hvMatch[MAXBUFF+1],
        ljMatch[MAXBUFF+1],
        hjMatch[MAXBUFF+1],
        lcMatch[MAXBUFF+1],
        CH1Match[MAXBUFF+1],
        CH2Match[MAXBUFF+1],
        CH3CHSMatch[MAXBUFF+1];
   REAL lvScore     = -1.0,
        hvScore     = -1.0,
        ljScore     = -1.0,
        hjScore     = -1.0,
        lcScore     = -1.0,
        CH1Score    = -1.0,
        CH2Score    = -1.0,
        CH3CHSScore = -1.0;

   /* Find out what the chain type is if it isn't specified             */
   if(chainType == CHAINTYPE_UNKNOWN)
   {
      lvScore  = ScanAgainstDB("light_v", seq, verbose, species, lvMatch);
      hvScore  = ScanAgainstDB("heavy_v", seq, verbose, species, hvMatch);
      lcScore  = ScanAgainstDB("light_c", seq, verbose, species, lcMatch);
      CH1Score = ScanAgainstDB("CH1",     seq, verbose, species, CH1Match);

      if((lvScore > hvScore) || (lcScore > CH1Score))
      {
         chainType = CHAINTYPE_LIGHT;
      }
      else if((hvScore > lvScore) || (CH1Score > lcScore))
      {
         chainType = CHAINTYPE_HEAVY;
      }
      
      fprintf(out, "Chain type: %s\n", CHAINTYPE(chainType));
   }

   switch(chainType)
   {
   case CHAINTYPE_LIGHT:
      if(lvScore < 0.0)
         lvScore = ScanAgainstDB("light_v", seq, verbose, species, lvMatch);
      if(lcScore < 0.0)
         lcScore = ScanAgainstDB("light_c", seq, verbose, species, lcMatch);

      if(lvScore > THRESHOLD_LV)
      {
         fprintf(out, "VL  : %f : %s\n", lvScore, lvMatch);
         ljScore = ScanAgainstDB("light_j", seq, verbose, species, ljMatch);
         if(ljScore > THRESHOLD_LJ)
         {
            fprintf(out, "JL  : %f : %s\n", ljScore, ljMatch);
         }
      }

      if(lcScore > THRESHOLD_LC)
      {
         fprintf(out, "CL  : %f : %s\n", lcScore, lcMatch);
      }
      
      break;
   case CHAINTYPE_HEAVY:
      if(hvScore < 0.0)
         hvScore  = ScanAgainstDB("heavy_v", seq, verbose, species, hvMatch);

      if(CH1Score < 0.0)
         CH1Score = ScanAgainstDB("CH1", seq, verbose, species, CH1Match);

      CH2Score    = ScanAgainstDB("CH2",     seq, verbose, species, CH2Match);
      CH3CHSScore = ScanAgainstDB("CH3-CHS", seq, verbose, species, CH3CHSMatch);
      
      if(hvScore > THRESHOLD_HV)
      {
         fprintf(out, "VH      : %f : %s\n", hvScore, hvMatch);
         hjScore = ScanAgainstDB("heavy_j", seq, verbose, species, hjMatch);
         if(hjScore > THRESHOLD_HJ)
         {
            fprintf(out, "JH      : %f : %s\n", hjScore, hjMatch);
         }
      }

      if(CH1Score > THRESHOLD_HC)
      {
         fprintf(out, "CH1     : %f : %s\n", CH1Score, CH1Match);
      }
      if(CH2Score > THRESHOLD_HC)
      {
         fprintf(out, "CH2     : %f : %s\n", CH2Score, CH2Match);
      }
      if(CH3CHSScore > THRESHOLD_HC)
      {
         fprintf(out, "CH3-CHS : %f : %s\n", CH3CHSScore, CH3CHSMatch);
      }
      
      break;
   default:
      fprintf(stderr, "Can't identify chain type!\n");
      exit (1);
   }
}

/************************************************************************/
REAL ScanAgainstDB(char *type, char *theSeq, BOOL verbose, char *species, char *match)
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
         if((species[0] == '\0') ||
            (strstr(header, species) != NULL))
         {
            REAL score;
         
            if(verbose)
               fprintf(stderr, "Comparing with %s\n", header);
            
            score = CompareSeqs(theSeq, seq);
            if(score > maxScore)
            {
               if(verbose)
                  fprintf(stderr, "Comparing with %s *** %.4f\n", header, score);
            
               maxScore = score;
               strncpy(bestMatch, header, MAXBUFF);
            }
            else if(score == maxScore)
            {
               /* If the scores are the same, choose the one with the better 
                  gene name
               */
               if(PreferHeader(header, bestMatch))
               {
                  if(verbose)
                     fprintf(stderr, "Comparing with %s *** %.4f\n", header, score);
            
                  maxScore = score;
                  strncpy(bestMatch, header, MAXBUFF);
               }
               else if(verbose)
               {
                  fprintf(stderr, "Comparing with %s (rejected %.4f)\n", header, score);
               }
            }
            else if(verbose)
            {
               fprintf(stderr, "Comparing with %s\n", header);
            }
         }
         
         free(seq);
      }
      FCLOSE(dbFp);
   }

   strncpy(match, bestMatch, MAXBUFF);
   return(maxScore);
}


/************************************************************************/
BOOL PreferHeader(char *newHeader, char *oldHeader)
{
   char newID[SMALLBUFF+1],
        oldID[SMALLBUFF+1];
   char newClass[SMALLBUFF],
      newDistal[SMALLBUFF],
      newFamily[SMALLBUFF],
      newSubclass[SMALLBUFF];
   int newAllele;
   char oldClass[SMALLBUFF],
      oldDistal[SMALLBUFF],
      oldFamily[SMALLBUFF],
      oldSubclass[SMALLBUFF];
   int oldAllele;
   int oldNum, newNum;
   
   GetDomainID(newHeader, newID, SMALLBUFF);
   GetDomainID(oldHeader, oldID, SMALLBUFF);

   FindFields(newID, newClass, newSubclass, newFamily, &newAllele, newDistal);
   FindFields(oldID, oldClass, oldSubclass, oldFamily, &oldAllele, oldDistal);

   /* if it's non-distal and the current best is distal, keep the new one */
   if((newDistal[0] == '\0') && (oldDistal[0] != '\0'))
      return(TRUE);
   
   /* If the sub-class is numeric and the old one isn't, keep the new one */
   oldNum = atoi(oldSubclass);
   newNum = atoi(newSubclass);
   
   if((newNum > 0) && (oldNum == 0))
      return(TRUE);

   /* If both are numeric and the new one is lower, keep the new one  */
   if((newNum && oldNum) &&
      (newNum < oldNum))
      return(TRUE);

   /* If the family is numeric and the old one isn't, keep the new one */
   oldNum = atoi(oldFamily);
   newNum = atoi(newFamily);
   
   if((newNum > 0) && (oldNum == 0))
      return(TRUE);

   /* If both are numeric and the new one is lower, keep the new one   */
   if((newNum && oldNum) &&
      (newNum < oldNum))
      return(TRUE);

   /* Keep the lowest allele */
   if(newAllele < oldAllele)
      return(TRUE);

   return(FALSE);
}

/************************************************************************/
void GetDomainID(char *header, char *id, int maxbuff)
{
   int start, i,
      headerLen;
   headerLen = strlen(header);
   
   for(start=0; start<headerLen; start++)
   {
      if(header[start] == '_')
      {
         start++;
         break;
      }
   }
   for(i=0; i<headerLen-start && i<maxbuff; i++)
   {
      id[i] = header[start+i];
      if(id[i] == '_')
      {
         break;
      }
   }
   id[i] = '\0';
}


/************************************************************************/
/* theSeq is the sequence of interest
   seq is the database sequence
*/
REAL CompareSeqs(char *theSeq, char *seq)
{
   int  score;
   char align1[HUGEBUFF+1],
        align2[HUGEBUFF+1];
   int  alignLen;
   int  shortSeqLen = MIN(strlen(theSeq), strlen(seq));
   
   score = blAlign(theSeq, strlen(theSeq),
                   seq, strlen(seq),
                   FALSE, /* verbose */
                   TRUE,  /* identity */
                   5,     /* penalty */
                   align1,
                   align2,
                   &alignLen);
#ifdef DEBUG
   printf("%s\n%s\n", align1, align2);
#endif
   
   
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
                     BOOL *verbose, int *chainType, char *species)
   ---------------------------------------------------------------------
   Input:   int      argc        Argument count
            char     **argv      Argument array
   Output:  char     *infile     Input filename (or blank string)
            char     *outfile    Output filename (or blank string)
            BOOL     *verbose    Verbose
            int      *chainType  Chain type
            char     *species    Species or blank string
   Returns: BOOL                 Success

   Parse the command line

   26.03.20 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, int *chainType, char *species)
{
   argc--;
   argv++;
   
   infile[0]   = outfile[0] = '\0';
   species[0] = '\0';
   
   
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
         case 's':
            argc--; argv++;
            if(!argc)
               return(FALSE);
            strncpy(species, argv[0], MAXBUFF);
            break;
         default:
            return(FALSE);
            break;
         }
         argc--; argv++;
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

