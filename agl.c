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
REAL ScanAgainstDB(char *type, char *seq, BOOL verbose, char *species, char *match, char *bestAlign1, char *bestAlign2);
REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2);
BOOL PreferHeader(char *newHeader, char *oldHeader);
void GetDomainID(char *header, char *id, int maxbuff);
void RemoveSequence(char *seq, char *align1, char *align2, BOOL verbose);
void PrintResult(FILE *out, char *domain, REAL score, char *match);




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
         if(in  != stdin)  fclose(in);
         if(out != stdout) fclose(out);
      }
      else
      {
         fprintf(stderr,"Unable to open input or output file\n");
         return(1);
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
        CH3CHSMatch[MAXBUFF+1],
        bestAlign1[HUGEBUFF+1],
        bestAlign2[HUGEBUFF+1];
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
      lvScore  = ScanAgainstDB("light_v", seq, verbose, species, lvMatch,  bestAlign1, bestAlign2);
      hvScore  = ScanAgainstDB("heavy_v", seq, verbose, species, hvMatch,  bestAlign1, bestAlign2);
      lcScore  = ScanAgainstDB("light_c", seq, verbose, species, lcMatch,  bestAlign1, bestAlign2);
      CH1Score = ScanAgainstDB("CH1",     seq, verbose, species, CH1Match, bestAlign1, bestAlign2);

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
         lvScore = ScanAgainstDB("light_v", seq, verbose, species, lvMatch, bestAlign1, bestAlign2);

      if(lvScore > THRESHOLD_LV)
      {
         RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
         fprintf(out, "VL  : %f : %s\n", lvScore, lvMatch);
         ljScore = ScanAgainstDB("light_j", seq, verbose, species, ljMatch, bestAlign1, bestAlign2);
         if(ljScore > THRESHOLD_LJ)
         {
            RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
            fprintf(out, "JL  : %f : %s\n", ljScore, ljMatch);
         }
      }

      if(lcScore < 0.0)
         lcScore = ScanAgainstDB("light_c", seq, verbose, species, lcMatch, bestAlign1, bestAlign2);
      if(lcScore > THRESHOLD_LC)
      {
         RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
         fprintf(out, "CL  : %f : %s\n", lcScore, lcMatch);
      }
      
      break;
   case CHAINTYPE_HEAVY:
      if(hvScore < 0.0)
         hvScore  = ScanAgainstDB("heavy_v", seq, verbose, species, hvMatch, bestAlign1, bestAlign2);

      if(hvScore > THRESHOLD_HV)
      {
         RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
         PrintResult(out, "VH", hvScore, hvMatch);
         
         hjScore = ScanAgainstDB("heavy_j", seq, verbose, species, hjMatch, bestAlign1, bestAlign2);
         if(hjScore > THRESHOLD_HJ)
         {
            RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
            PrintResult(out, "JH", hjScore, hjMatch);
         }
      }

      if(CH1Score < 0.0)
         CH1Score = ScanAgainstDB("CH1", seq, verbose, species, CH1Match, bestAlign1, bestAlign2);
      if(CH1Score > THRESHOLD_HC)
      {
         RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
         PrintResult(out, "CH1", CH1Score, CH1Match);
      }

      
      CH2Score    = ScanAgainstDB("CH2",     seq, verbose, species, CH2Match, bestAlign1, bestAlign2);
      if(CH2Score > THRESHOLD_HC)
      {
         RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
         PrintResult(out, "CH2", CH2Score, CH2Match);
      }

      CH3CHSScore = ScanAgainstDB("CH3-CHS", seq, verbose, species, CH3CHSMatch, bestAlign1, bestAlign2);
      if(CH3CHSScore > THRESHOLD_HC)
      {
         RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
         PrintResult(out, "CH3-CHS", CH3CHSScore, CH3CHSMatch);
      }
      
      break;
   default:
      fprintf(stderr, "Can't identify chain type!\n");
      exit (1);
   }
}

/************************************************************************/
REAL ScanAgainstDB(char *type, char *theSeq, BOOL verbose, char *species, char *match, char *bestAlign1, char *bestAlign2)
{
   char filename[MAXBUFF+1];
   BOOL noEnv;
   FILE *dbFp = NULL;
   REAL maxScore = 0.0;
   char bestMatch[MAXBUFF+1];
   char align1[HUGEBUFF+1],
        align2[HUGEBUFF+1];
   
   sprintf(filename, "%s/%s.dat", AGLDATADIR, type);

   if(verbose)
   {
      fprintf(stderr,"\n\nChecking %s\n", type);
   }
   
   if((dbFp = blOpenFile(filename, AGLDATADIR, "r", &noEnv))!=NULL)
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
            
            score = CompareSeqs(theSeq, seq, align1, align2);
            if(score > maxScore)
            {
               if(verbose)
                  fprintf(stderr, "Comparing with %s *** %.4f\n", header, score);
            
               maxScore = score;
               strncpy(bestMatch, header, MAXBUFF);
               strncpy(bestAlign1, align1, HUGEBUFF);
               strncpy(bestAlign2, align2, HUGEBUFF);
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
                  strncpy(bestMatch,  header, MAXBUFF);
                  strncpy(bestAlign1, align1, HUGEBUFF);
                  strncpy(bestAlign2, align2, HUGEBUFF);
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
REAL CompareSeqs(char *theSeq, char *seq, char *align1, char *align2)
{
   int  score;
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
   align1[alignLen] = align2[alignLen] = '\0';
   
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

void RemoveSequence(char *seq, char *align1, char *align2, BOOL verbose)
{
   int alnLen = MAX(strlen(align1), strlen(align2));
   int i, j, start, stop;
   char buffer[HUGEBUFF+1];

   if(verbose)
   {
      fprintf(stderr, "Input seq:        %s\n", seq);
      fprintf(stderr, "Alignment:        %s\n", align1);
      fprintf(stderr, "                  %s\n", align2);
      fprintf(stderr, "Alignment length: %d\n", alnLen);
   }
   
   
   /* Find the start of the aligned region */
   for(i=0; i<alnLen; i++)
   {
      if((align1[i] != '-') &&
         (align2[i] != '-'))
      {
         start = i;
         break;
      }
   }
   /* Find the end of the aligned region */
   for(i=alnLen-1; i>start; i--)
   {
      if((align1[i] != '-') &&
         (align2[i] != '-'))
      {
         stop = i;
         break;
      }
   }

   if(verbose)
   {
      fprintf(stderr, "Start: %d Stop: %d\n", start, stop);
   }
   
   /* Copy up to the start */
   j=0;
   for(i=0; i<start; i++)
   {
      if(align1[i] != '-')
      {
         buffer[j++] = align1[i];
      }

   }
   /* Copy after the end */
   for(i=stop+1; i<alnLen; i++)
   {
      if(align1[i] != '-')
      {
         buffer[j++] = align1[i];
      }
   }
   buffer[j] = '\0';
   /* Copy this back into seq */
   i=strlen(seq);
   strncpy(seq, buffer, i);
   seq[i] = '\0';
   
   if(verbose)
   {
      fprintf(stderr, "Output seq:       %s\n", seq);
   }
}

void PrintResult(FILE *out, char *domain, REAL score, char *match)
{
   char id[32],
      frame[8],
      species[32];
   
   match = strchr(match, '_');
   match++;
   strncpy(id, match, 31);
   TERMAT(id, '_');

   match = strchr(match, '_');
   match++;
   strncpy(frame, match, 7);
   TERMAT(frame, '_');
   
   match = strchr(match, '_');
   match++;
   strncpy(species, match, 31);
   
   fprintf(out, "%-7s : %.4f : %-12s : %s : %s\n", domain, score, id, frame, species);
}
