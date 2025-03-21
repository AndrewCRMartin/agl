/************************************************************************/
/**

   Program:    agl (Assign Germ Line)
   \file       agl.c
   
   \version    V1.8
   \date       19.03.25
   \brief      Assigns IMGT germline
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2020-25
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   LICENSE

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file LICENSE.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   31.03.20  Original
   V1.1   14.04.20  Added -a
   V1.2   11.06.21  Now looks in .../share/agl/data below the location of
                    the executable before trying $AGLDATA
   V1.3   14.06.21  Added printing of mismatches with -a
   V1.4   13.06.22  Now uses a window in the alignment
   V1.5   14.11.22  Now supports FASTA input with multiple sequences
   V1.6   26.04.23  Added D-segment option and fixed gap penalty in
                    alignment
   V1.7   22.04.23  Adds info on IGHG1*08 and IGHG1*15
   V1.8   19.03.25  Added hinges

*************************************************************************/
/* Includes
*/
#define USEPATH 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioplib/seq.h"
#include "bioplib/sequtil.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "findfields.h"
#include "agl.h"
#include "whereami/whereami.h"

/************************************************************************/
/* Defines and macros
*/
#define BLOCK_UNDEFINED 0
#define BLOCK_X_V       1
#define BLOCK_D         2

#define CHAINTYPE(x) (                            \
   (x)==CHAINTYPE_LIGHT ? "Light" :               \
    ((x)==CHAINTYPE_HEAVY ? "Heavy" : "Unknown"))

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, BOOL *showAlignment, int *chainType,
                  char *species, char *dataDir, BOOL *doDSegment);
void ProcessSeq(FILE *out, char *seq, BOOL verbose, BOOL showAlignment,
                int chainType, char *species, char *dataDir,
                BOOL doDSegment);
REAL ScanAgainstDB(char *type, char *seq, BOOL verbose, char *species,
                   char *match, char *bestAlign1, char *bestAlign2,
                   char *dataDir);
REAL CompareSeqs(char *theSeq, char *seq, int window,
                 char *align1, char *align2, BOOL noScale);
BOOL PreferHeader(char *newHeader, char *oldHeader);
void GetDomainID(char *header, char *id, int maxbuff);
void RemoveSequence(char *seq, char *align1, char *align2, BOOL verbose);
void PrintResult(FILE *out, char *domain, REAL score, char *match);
int CalculateDbLen(char *seq);
int CalcShortSeqLen(char *align1, char *align2);
void PrintAlignment(FILE *out, char *align1, char *align2);
void DoDSegment(FILE *out, char *seq, char *species, char *dataDir,
                char *hvBestAlign1, char *hvBestAlign2, 
                char *hjBestAlign1, char *hjBestAlign2,
                BOOL verbose, BOOL showAlignment);
void CopyDSegment(char *DSeq, char *seq);
BOOL PrintSpecialMatches(FILE *out, char *CH1, char *CH2, char *CH3);

#ifdef USEPATH
char *FindPath(void);
BOOL DirectoryExists(char *dirName);
#endif

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

   - 31.03.20 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF+1],
        outfile[MAXBUFF+1],
        species[MAXBUFF+1],
        dataDir[MAXBUFF+1];
   BOOL verbose       = FALSE,
        showAlignment = FALSE,
        doDSegment    = FALSE;
   int  chainType     = CHAINTYPE_UNKNOWN;


   if(ParseCmdLine(argc, argv, infile, outfile, &verbose, &showAlignment,
                   &chainType, species, dataDir, &doDSegment))
   {
      FILE *in  = stdin,
           *out = stdout;

      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         char header[MAXBUFF+1];
         char *seq = NULL;
         BOOL ok = FALSE;
         static char fastaBuffer[MAXBUFF+1];

         while((seq = blReadFASTAExtBuffer(in, header, MAXBUFF,
                                           fastaBuffer, MAXBUFF))!=NULL)
         {
            ok = TRUE;
            fprintf(out, "%s\n", header);

            ProcessSeq(out, seq, verbose, showAlignment,
                       chainType, species, dataDir, doDSegment);
            free(seq);
         }

         if(!ok)
         {
            fprintf(stderr, "No sequence found for %s\n", header);
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
/*>void ProcessSeq(FILE *out, char *seq, BOOL verbose, BOOL showAlignment,
                   int chainType, char *species, char *dataDir,
                   BOOL doDSegment)
   ------------------------------------------------------------------
*//**
   \param[in]   out            Output file pointer
   \param[in]   seq            Sequence to analyze
   \param[in]   verbose        Verbose output
   \param[in]   showAlignment  Show the alignment of each region
   \param[in]   chainType      Type of chain (or unknown)
   \param[in]   species        "Homo", "Mus" or blank
   \param[in]   dataDir        Data directory
   \param[in]   doDSegment     Handle D-segment for heavy chains

   - 31.03.20 Original   By: ACRM
   - 14.04.20 Added showAlignment
   - 26.04.23 Added D-segment handling
   - 19.03.25 Added hinge handling
CHECKED - ERROR IS IN ScanAgainstDB
*/
void ProcessSeq(FILE *out, char *seq, BOOL verbose, BOOL showAlignment,
                int chainType, char *species, char *dataDir,
                BOOL doDSegment)

{
   char        lvMatch[MAXBUFF+1],
               hvMatch[MAXBUFF+1],
               ljMatch[MAXBUFF+1],
               hjMatch[MAXBUFF+1],
               lcMatch[MAXBUFF+1],
               CH1Match[MAXBUFF+1],
               CH2Match[MAXBUFF+1],
               hingeMatch[MAXBUFF+1],
               CH3CHSMatch[MAXBUFF+1];
   static char lvBestAlign1[HUGEBUFF+1],
               lvBestAlign2[HUGEBUFF+1],
               hvBestAlign1[HUGEBUFF+1],
               hvBestAlign2[HUGEBUFF+1],
               lcBestAlign1[HUGEBUFF+1],
               lcBestAlign2[HUGEBUFF+1],
               CH1BestAlign1[HUGEBUFF+1],
               CH1BestAlign2[HUGEBUFF+1],
               CH2BestAlign1[HUGEBUFF+1],
               CH2BestAlign2[HUGEBUFF+1],
               CH3BestAlign1[HUGEBUFF+1],
               CH3BestAlign2[HUGEBUFF+1],
               hingeBestAlign1[HUGEBUFF+1],
               hingeBestAlign2[HUGEBUFF+1],
               bestAlign1[HUGEBUFF+1],
               bestAlign2[HUGEBUFF+1];
   REAL        lvScore     = -1.0,
               hvScore     = -1.0,
               ljScore     = -1.0,
               hjScore     = -1.0,
               lcScore     = -1.0,
               CH1Score    = -1.0,
               CH2Score    = -1.0,
               hingeScore  = -1.0,
               CH3CHSScore = -1.0;

   /* Find out what the chain type is if it isn't specified             */
   if(chainType == CHAINTYPE_UNKNOWN)
   {
      lvScore = ScanAgainstDB("light_v", seq, verbose, species,
                              lvMatch,  lvBestAlign1, lvBestAlign2,
                              dataDir);
      hvScore = ScanAgainstDB("heavy_v", seq, verbose, species,
                              hvMatch,  hvBestAlign1, hvBestAlign2,
                              dataDir);
      if((lvScore > hvScore) && (lvScore > THRESHOLD_LV))
      {
         chainType = CHAINTYPE_LIGHT;
      }
      else if((hvScore > lvScore) && (hvScore > THRESHOLD_HV))
      {
         chainType = CHAINTYPE_HEAVY;
      }
      
      if(chainType == CHAINTYPE_UNKNOWN)
      {
         lcScore  = ScanAgainstDB("light_c", seq, verbose, species,
                                  lcMatch,  lcBestAlign1,  lcBestAlign2,
                                  dataDir);
         CH1Score = ScanAgainstDB("CH1",     seq, verbose, species,
                                  CH1Match, CH1BestAlign1, CH1BestAlign2,
                                  dataDir);
      }

      if((lcScore > CH1Score) && (lcScore > THRESHOLD_LC))
      {
         chainType = CHAINTYPE_LIGHT;
      }
      else if((CH1Score > lcScore) && (CH1Score > THRESHOLD_HC))
      {
         chainType = CHAINTYPE_HEAVY;
      }
      fprintf(out, "# Chain type: %s\n", CHAINTYPE(chainType));
   }

   switch(chainType)
   {
   case CHAINTYPE_LIGHT:
      if(lvScore < 0.0)
         lvScore = ScanAgainstDB("light_v", seq, verbose, species,
                                 lvMatch, lvBestAlign1, lvBestAlign2,
                                 dataDir);
      if(lvScore > THRESHOLD_LV)
      {
         PrintResult(out, "VL", lvScore, lvMatch);
         if(showAlignment)
            PrintAlignment(out, lvBestAlign1, lvBestAlign2);

         ljScore = ScanAgainstDB("light_j", seq, verbose, species,
                                 ljMatch, bestAlign1, bestAlign2,
                                 dataDir);
         if(ljScore > THRESHOLD_LJ)
         {
#ifdef REMOVESEQS
            RemoveSequence(seq, lvBestAlign1, lvBestAlign2, verbose);
            RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
#endif
            PrintResult(out, "JL", ljScore, ljMatch);
            if(showAlignment)
               PrintAlignment(out, bestAlign1, bestAlign2);
         }
      }

      if(lcScore < 0.0)
         lcScore = ScanAgainstDB("light_c", seq, verbose, species,
                                 lcMatch, lcBestAlign1, lcBestAlign2,
                                 dataDir);
      if(lcScore > THRESHOLD_LC)
      {
#ifdef REMOVESEQS
         RemoveSequence(seq, lcBestAlign1, lcBestAlign2, verbose);
#endif
         PrintResult(out, "CL", lcScore, lcMatch);
         if(showAlignment)
            PrintAlignment(out, lcBestAlign1, lcBestAlign2);
      }
      
      break;

   case CHAINTYPE_HEAVY:
      if(CH3CHSScore < 0.0)
         CH3CHSScore = ScanAgainstDB("CH3-CHS", seq, verbose, species,
                                     CH3CHSMatch,
                                     CH3BestAlign1, CH3BestAlign2,
                                     dataDir);
#ifdef REMOVESEQS
      if(CH3CHSScore > THRESHOLD_HC)
         RemoveSequence(seq, CH3BestAlign1, CH3BestAlign2, verbose);
#endif      
      if(CH2Score < 0.0)
         CH2Score = ScanAgainstDB("CH2", seq, verbose, species,
                                  CH2Match, CH2BestAlign1, CH2BestAlign2,
                                  dataDir);
#ifdef REMOVESEQS
      if(CH2Score > THRESHOLD_HC)
         RemoveSequence(seq, CH2BestAlign1, CH2BestAlign2, verbose);
#endif
      if(CH1Score < 0.0)
         CH1Score = ScanAgainstDB("CH1", seq, verbose, species,
                                  CH1Match, CH1BestAlign1, CH1BestAlign2,
                                  dataDir);
#ifdef REMOVESEQS
      if(CH1Score > THRESHOLD_HC)
         RemoveSequence(seq, CH1BestAlign1, CH1BestAlign2, verbose);
#endif      
      if(hvScore < 0.0)
         hvScore  = ScanAgainstDB("heavy_v", seq, verbose, species,
                                  hvMatch, hvBestAlign1, hvBestAlign2,
                                  dataDir);

      if(CH1Match[0] && CH2Match[0])
      {
         hingeScore  = ScanAgainstDB("hinges", seq, verbose, species,
                                     hingeMatch, hingeBestAlign1,
                                     hingeBestAlign2, dataDir);
      }
      
      if(hvScore > THRESHOLD_HV)
      {
#ifdef REMOVESEQS
         RemoveSequence(seq, hvBestAlign1, hvBestAlign2, verbose);
#endif
         PrintResult(out, "VH", hvScore, hvMatch);
         if(showAlignment)
            PrintAlignment(out, hvBestAlign1, hvBestAlign2);
         
         hjScore = ScanAgainstDB("heavy_j", seq, verbose, species,
                                 hjMatch, bestAlign1, bestAlign2,
                                 dataDir);
         if(hjScore > THRESHOLD_HJ)
         {
#ifdef REMOVESEQS
            RemoveSequence(seq, bestAlign1, bestAlign2, verbose);
#endif
            PrintResult(out, "JH", hjScore, hjMatch);
            if(showAlignment)
               PrintAlignment(out, bestAlign1, bestAlign2);

            /* If we have found both V and J, then we can try to find
               the D segment (if required)
            */
            if(doDSegment)
            {
               DoDSegment(out, seq, species, dataDir,
                          hvBestAlign1, hvBestAlign2,  /* V             */
                          bestAlign1,   bestAlign2,    /* J             */
                          verbose, showAlignment);
            }
         }
      }

      if(CH1Score > THRESHOLD_HC)
      {
         PrintResult(out, "CH1", CH1Score, CH1Match);
         if(showAlignment)
            PrintAlignment(out, CH1BestAlign1, CH1BestAlign2);
      }
      
      if(hingeScore > THRESHOLD_HINGE)
      {
         PrintResult(out, "HINGE", hingeScore, hingeMatch);
         if(showAlignment)
            PrintAlignment(out, hingeBestAlign1, hingeBestAlign2);
      }
      
      if(CH2Score > THRESHOLD_HC)
      {
         PrintResult(out, "CH2", CH2Score, CH2Match);
         if(showAlignment)
            PrintAlignment(out, CH2BestAlign1, CH2BestAlign2);
      }
      
      if(CH3CHSScore > THRESHOLD_HC)
      {
         PrintResult(out, "CH3-CHS", CH3CHSScore, CH3CHSMatch);
         if(showAlignment)
            PrintAlignment(out, CH3BestAlign1, CH3BestAlign2);
      }

      
      /* 22.02.24 Added to print info on special cases where mixed allelic
         variants actually indicate a higher numbered allele
      */
      PrintSpecialMatches(out, CH1Match, CH2Match, CH3CHSMatch);

      break;

   default:
      fprintf(stderr, "Error: Can't identify chain type!\n");
   }
}


/************************************************************************/
/*>REAL ScanAgainstDB(char *type, char *theSeq, BOOL verbose, 
                      char *species, char *match, char *bestAlign1, 
                      char *bestAlign2, char *dataDir)
   -----------------------------------------------------------------------
*//**
   \param[in]  type       The database type against which we scan
   \param[in]  theSeq     The sequences to test
   \param[in]  verbose    Verbose output
   \param[in]  species    "Homo", "Mus" or blank
   \param[out] match      The best matching entry
   \param[out] bestAlign1 The alignment of our sequence
   \param[out] bestAlign2 The alignment of the database sequence
   \param[in]  dataDir    Data directory
   \return                The score for the match

   Scans a sequence against the specified database.

   - 31.03.20 Original   By: ACRM
   - 11.06.21 Added USEPATH code
   - 13.06.22 Added window size of 2 for C regions and 10 for others
   - 19.03.25 Initialize match
*/
REAL ScanAgainstDB(char *type, char *theSeq, BOOL verbose, char *species,
                   char *match, char *bestAlign1, char *bestAlign2,
                   char *dataDir)
{
   char        filename[MAXBUFF+1],
               bestMatch[MAXBUFF+1];
   BOOL        noEnv;
   FILE        *dbFp = NULL;
   REAL        maxScore = 0.0;
   int         maxDbLen = 0;
   static char align1[HUGEBUFF+1],
               align2[HUGEBUFF+1];

   match[0] = '\0';
   
   if(dataDir[0] != '\0')
   {
      sprintf(filename, "%s/%s.dat", dataDir, type);
   }
   else
   {
#ifdef USEPATH
      char *path;
      /* Get path to executable                                         */
      if((path = FindPath())!=NULL)
      {
         /* Append the location below the binary directory and see if it
            exists. If it does then set the filename to include this
            path. If not, then just save the filename itself
         */
         sprintf(filename, "%s/share/agl/data", path);
         if(DirectoryExists(filename))
         {
            sprintf(filename, "%s/share/agl/data/%s.dat", path, type);
         }
         else
         {
            sprintf(filename, "%s.dat", type);
         }
         
         free(path);
      }
      else
      {
         sprintf(filename, "%s.dat", type);
      }
#else
      sprintf(filename, "%s.dat", type);
#endif
   }

   
   if(verbose)
      fprintf(stderr,"\n\nChecking %s\n", type);

   if((dbFp = blOpenFile(filename, AGLDATADIR, "r", &noEnv))!=NULL)
   {
      char header[MAXBUFF+1];
      char *seq = NULL;
      static char fastaBuffer[MAXBUFF+1];

      while((seq = blReadFASTAExtBuffer(dbFp, header, MAXBUFF,
                                        fastaBuffer, MAXBUFF))!=NULL)
      {
         if((species[0] == '\0') ||
            (strstr(header, species) != NULL))
         {
            REAL score;
            int  dbLen;
            int  window = 10;
            BOOL noScale = FALSE;

            if(header[1] == 'C')
            {
               /* Use a window of 10 by default, or a window of 2 for 
                  constant regions
               */
               window = 2;
            }
            else if(header[1] == 'D')
            {
               /* Don't scale the score if it's a D-segment             */
               noScale = TRUE;
            }
            
            score = CompareSeqs(theSeq, seq, window, align1, align2,
                                noScale);
            dbLen = CalculateDbLen(align2);

#ifdef DEBUG
            if(verbose)
               fprintf(stderr, "Comparing with %s {%f, %f} {%d, %d}\n",
                       header, score, maxScore, dbLen, maxDbLen);
#endif 
            if(/* Length has increased and score not decreased much     */
               ((dbLen > maxDbLen) &&
                (score > (maxScore-THRESHOLD_LEN_SCORE)))   ||
               /* Score has increased and length not shorter            */
               ((score > maxScore) && (dbLen >= maxDbLen))  ||
               /* Score has increased significantly while length is 
                  shorter 
               */
               ((score >  (maxScore + THRESHOLD_SCORE_INC)) &&
                (dbLen >= (maxDbLen-THRESHOLD_SCORE_LEN))))
            {
               if(verbose)
                  fprintf(stderr, "Comparing with %s *** %.4f\n",
                          header, score);
            
               maxScore = score;
               maxDbLen = dbLen;
               strncpy(bestMatch,  header, MAXBUFF);
               strncpy(bestAlign1, align1, HUGEBUFF);
               strncpy(bestAlign2, align2, HUGEBUFF);
            }
            else if(score == maxScore)
            {
               /* If the scores are the same, choose the one with the 
                  better gene name
               */
               if(PreferHeader(header, bestMatch))
               {
                  if(verbose)
                     fprintf(stderr, "Comparing with %s *** %.4f \
(Chosen on name)\n", header, score);
            
                  maxScore = score;
                  strncpy(bestMatch,  header, MAXBUFF);
                  strncpy(bestAlign1, align1, HUGEBUFF);
                  strncpy(bestAlign2, align2, HUGEBUFF);
               }
               else if(verbose)
               {
                  fprintf(stderr, "Comparing with %s (rejected %.4f)\n",
                          header, score);
               }
            }
            else if(verbose)
            {
               fprintf(stderr, "Comparing with %s (%.4f)\n",
                       header, score);
            }
         }
         
         free(seq);
         seq=NULL;
      }
      FCLOSE(dbFp);
   }
   else
   {
      fprintf(stderr, "\nError (agl): Unable to open data file (%s)\n",
              filename);
      if(noEnv)
      {
         fprintf(stderr, "   Set the environment variable '%s' to the \
location of the processed sequence files.\n",
                 AGLDATADIR);
      }
      fprintf(stderr, "\n");
      
      exit(1);
   }

   strncpy(match, bestMatch, MAXBUFF);
   return(maxScore);

}


/************************************************************************/
int CalculateDbLen(char *seq)
{
   int i,
       len=0;

   for(i=0; i<strlen(seq); i++)
   {
      if(seq[i] != '-')
         len++;
   }
   return(len);
}

/************************************************************************/
/*>BOOL PreferHeader(char *newHeader, char *oldHeader)
   ---------------------------------------------------
*//**
   \param[in]  newHeader  The FAA header for the new entry
   \param[in]  oldHeader  The FAA header for the current best entry
   \return                Should be use the new entry?

   Used when two hits score the same. Tests the names of the hits and
   sees if the new one has a preferred name

   - 31.03.20 Original   By: ACRM
CHECKED
*/
BOOL PreferHeader(char *newHeader, char *oldHeader)
{
   char newID[SMALLBUFF+1],
        oldID[SMALLBUFF+1],
        newClass[SMALLBUFF],
        newDistal[SMALLBUFF],
        newFamily[SMALLBUFF],
        newSubclass[SMALLBUFF],
        oldClass[SMALLBUFF],
        oldDistal[SMALLBUFF],
        oldFamily[SMALLBUFF],
        oldSubclass[SMALLBUFF];
   int  newAllele,
        oldAllele,
        oldNum, newNum;

   GetDomainID(newHeader, newID, SMALLBUFF);
   GetDomainID(oldHeader, oldID, SMALLBUFF);

   FindFields(newID, newClass, newSubclass, newFamily, &newAllele,
              newDistal);
   FindFields(oldID, oldClass, oldSubclass, oldFamily, &oldAllele,
              oldDistal);

   /* If it's non-distal and the current best is distal, 
      keep the new one 
   */
   if((newDistal[0] == '\0') && (oldDistal[0] != '\0'))
      return(TRUE);
   
   /* If the sub-class is numeric and the old one isn't, 
      keep the new one 
   */
   oldNum = atoi(oldSubclass);
   newNum = atoi(newSubclass);
   if((newNum > 0) && (oldNum == 0))
      return(TRUE);

   /* If both are numeric and the new one is lower, 
      keep the new one  
   */
   if((newNum && oldNum) && (newNum < oldNum))
      return(TRUE);
   if(newNum > oldNum)
      return(FALSE);

   /* If the family is numeric and the old one isn't,
      keep the new one 
   */
   oldNum = atoi(oldFamily);
   newNum = atoi(newFamily);
   if((newNum > 0) && (oldNum == 0))
      return(TRUE);

   /* If both are numeric and the new one is lower,
      keep the new one   
   */
   if((newNum && oldNum) && (newNum < oldNum))
      return(TRUE);
   if(newNum > oldNum)
      return(FALSE);

   /* Keep the lowest allele                                            
   */
   if(newAllele < oldAllele)
      return(TRUE);

   return(FALSE);
}


/************************************************************************/
/*>void GetDomainID(char *header, char *id, int maxbuff)
   -----------------------------------------------------
*//**
   \param[in]  header   The Fasta header
   \param[out] id       The ID extracted from the header
   \param[out] maxbuff  Maximum buffer size for the ID

   Extracts the IMGT domain identifier from the FASTA header

   - 31.03.20 Original   By: ACRM
*/
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
/*>REAL CompareSeqs(char *theSeq, char *seq, int window, 
                    char *align1, char *align2, BOOL noScale)
   ---------------------------------------------------------------------
*//**
   \param[in]   theSeq   the sequence of interest
   \param[in]   seq      the database sequence
   \param[in]   window   Window size
   \param[out]  align1   Alignment of our sequence
   \param[out]  align2   Alignment of database sequence
   \param[in]   noScale  Don't scale score by length (for D-segment)
   \return               Score for alignment

   - 31.03.20 Original   By: ACRM
   - 13.06.22 Changed to use a window in the alignment to speed it up
              Added window size as a parameter
   - 26.04.23 Added noScale
   - 27.04.23 Changed extension penalty from 5 to 1
*/
REAL CompareSeqs(char *theSeq, char *seq, int window,
                 char *align1, char *align2, BOOL noScale)
{
   int  score;
   int  alignLen;
   int  shortSeqLen = MIN(strlen(theSeq), strlen(seq));
   
   score = blAffinealignWindow(theSeq, strlen(theSeq),
                               seq, strlen(seq),
                               FALSE,  /* verbose            */
                               TRUE,   /* identity           */
                               5,      /* opening penalty    */
                               1,      /* extension penalty  */
                               window, /* window             */
                               align1,
                               align2,
                               &alignLen);
   align1[alignLen] = align2[alignLen] = '\0';

   shortSeqLen = CalcShortSeqLen(align1, align2);

#ifdef DEBUG   
   fprintf(stderr, "\n>>>%s\n", align1);
   fprintf(stderr, ">>>%s %d\n", align2, shortSeqLen);
#endif

   if(noScale)
      return((REAL)score);
   
   return((REAL)score / (REAL)shortSeqLen);
}


/************************************************************************/
int CalcShortSeqLen(char *align1, char *align2)
{
   int start, stop, i, len1=0, len2=0,
      alnLen = strlen(align1);
   

   /* Find where the alignment starts                                   */
   for(start=0; start<alnLen; start++)
   {
      if((align1[start] != '-') && (align2[start] != '-'))
         break;
   }
   /* Find where the alignment stops                                    */
   for(stop=alnLen-1; stop>start; stop--)
   {
      if((align1[stop] != '-') && (align2[stop] != '-'))
         break;
   }
   /* Step between start and stop and calculate the number of residues 
      in each sequence
   */
   for(i=start; i<=stop; i++)
   {
      if(align1[i] != '-')
         len1++;
      if(align2[i] != '-')
         len2++;
   }

   return(MIN(len1, len2));
}


/************************************************************************/
/*>void Usage(void)
   ----------------
 *//**
   Print usage message.

-  31.03.20 Original    By: ACRM
-  14.03.20 V1.1 Added -a
-  11.06.21 V1.2
-  13.06.22 V1.4
-  14.11.22 V1.5
-  26.04.23 V1.6
-  22.02.24 V1.7
-  19.03.25 V1.8
*/
void Usage(void)
{
   printf("\nagl V1.8 (c) 2020-25 UCL, Prof. Andrew C.R. Martin\n\n");

   printf("Usage: agl [-H|-L] [-D] [-s species] [-d datadir] [-v] [-a] \
[file.faa [out.txt]]\n");
   printf("           -H Heavy chain\n");
   printf("           -L Light chain\n");
   printf("           -D Do the D-segment with heavy chains\n");
   printf("           -s Specify a species (Homo or Mus)\n");
   printf("           -d Specify data directory\n");
   printf("           -v Verbose\n");
   printf("           -a Show alignments and number of mismatches\n");

   printf("\nagl (Assign Germ Line) is a program for assigning IMGT \
germlines to\n");
   printf("antibody light or heavy chains. It does a 3-frame \
translation of the\n");
   printf("IMGT germline sequences and uses sequence identity at the \
amino acid\n");
   printf("level to find the best match. If more than one sequence \
scores the\n");
   printf("same, it takes no-distal over distal and then selects on \
the basis\n");
   printf("of the lowest sub-class, family and allele number.\n");

   printf("\nAs of V1.5, the FASTA file may contain multiple \
sequences.\n");
   
   printf("\nIf a data directory is not specified using -d, it will \
first look in the\n");
   printf("share/agl/data directory below the location of the \
executable, then in the\n");
   printf("current directory for files and then in the directory \
specified using the\n");
   printf("environment variable %s\n\n", AGLDATADIR);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *verbose, BOOL *showAlignment, 
                     int *chainType, char *species, char *dataDir,
                     BOOL *doDSegment)
   ---------------------------------------------------------------------
*//**
   \param[in]   argc           Argument count
   \param[in]   **argv         Argument array
   \param[out]  *infile        Input filename (or blank string)
   \param[out]  *outfile       Output filename (or blank string)
   \param[out]  *verbose       Verbose
   \param[out]  *showAlignment Show the alignments
   \param[out]  *chainType     Chain type
   \param[out]  *species       Species or blank string
   \param[out]  *dataDir       Data directory or blank string
   \param[out]  *doDSegment    Handle the D segment for heavy chains
   \return                     Success

   Parse the command line

-  31.03.20 Original    By: ACRM
-  26.04.23 Added -D/doDSegment
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, BOOL *showAlignment, int *chainType,
                  char *species, char *dataDir,
                  BOOL *doDSegment)
{
   argc--;
   argv++;
   
   infile[0]  = outfile[0] = '\0';
   species[0] = dataDir[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
            break;
         case 'a':
            *showAlignment = TRUE;
            break;
         case 'D':
            *doDSegment = TRUE;
            break;
         case 'L':
         case 'l':
            if(*chainType != CHAINTYPE_UNKNOWN)
               return(FALSE);
            *chainType = CHAINTYPE_LIGHT;
            break;
         case 'H':
            if(*chainType != CHAINTYPE_UNKNOWN)
               return(FALSE);
            *chainType = CHAINTYPE_HEAVY;
            break;
         case 'h':
            return(FALSE);
         case 's':
            argc--; argv++;
            if(!argc)
               return(FALSE);
            strncpy(species, argv[0], MAXBUFF);
            break;
         case 'd':
            argc--; argv++;
            if(!argc)
               return(FALSE);
            strncpy(dataDir, argv[0], MAXBUFF);
            break;
         default:
            return(FALSE);
            break;
         }
         argc--; argv++;
      }
      else
      {
         /* Check that there are no more than 2 arguments left          */
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


/************************************************************************/
/*>void RemoveSequence(char *seq, char *align1, char *align2, 
                       BOOL verbose)
   -----------------------------------------------------------
*//**
   \param[in,out]   seq        The sequence
   \param[in]       align1     Our sequence aligned
   \param[in]       align2     Database sequence aligned
   \param[in]       verbose    Verbose output

   Finds the region of the sequence that has been aligned with a 
   database sequence and removes it.

   - 31.03.20 Original   By: ACRM
*/
void RemoveSequence(char *seq, char *align1, char *align2, BOOL verbose)
{
   int         i, j,
               start, stop,
               alnLen = MAX(strlen(align1), strlen(align2));
   static char buffer[HUGEBUFF+1];

   if(verbose)
   {
      fprintf(stderr, "Input seq:        %s\n", seq);
      fprintf(stderr, "Alignment:        %s\n", align1);
      fprintf(stderr, "                  %s\n", align2);
      fprintf(stderr, "Alignment length: %d\n", alnLen);
   }
   
   /* Find the start of the aligned region                              */
   for(i=0; i<alnLen; i++)
   {
      if((align1[i] != '-') &&
         (align2[i] != '-'))
      {
         start = i;
         break;
      }
   }
   
   /* Find the end of the aligned region                                */
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
      fprintf(stderr, "Start: %d Stop: %d\n", start, stop);

#ifdef OLD
   /* Copy up to the start                                              */
   j=0;
   for(i=0; i<start; i++)
   {
      if(align1[i] != '-')
      {
         buffer[j++] = align1[i];
      }

   }
   
   /* Copy after the end                                                */
   for(i=stop+1; i<alnLen; i++)
   {
      if(align1[i] != '-')
      {
         buffer[j++] = align1[i];
      }
   }
   buffer[j] = '\0';
   
   /* Copy this back into seq                                           */
   i = strlen(seq);
   strncpy(seq, buffer, i);
   seq[i] = '\0';
#else
   for(i=0; i<strlen(seq); i++)
   {
      if((i>=start) && (i<=stop))
      {
         seq[i] = 'X';
      }
   }
#endif
   
   if(verbose)
      fprintf(stderr, "Output seq:       %s\n", seq);
}


/************************************************************************/
/*>void PrintResult(FILE *out, char *domain, REAL score, char *match)
   ------------------------------------------------------------------
*//**
   \param[in]   out        Output file pointer
   \param[in]   domain     Domain label
   \param[in]   score      Score for match
   \param[in]   match      Match FASTA header

   Prints the match information

   - 31.03.20 Original   By: ACRM
*/
void PrintResult(FILE *out, char *domain, REAL score, char *match)
{
   char id[LABELBUFF],
        frame[LABELBUFF],
        species[LABELBUFF];
   
   match = strchr(match, '_');
   match++;
   strncpy(id, match, LABELBUFF-1);
   TERMAT(id, '_');

   match = strchr(match, '_');
   match++;
   strncpy(frame, match, LABELBUFF-1);
   TERMAT(frame, '_');
   
   match = strchr(match, '_');
   match++;
   strncpy(species, match, LABELBUFF-1);
   
   fprintf(out, "%-7s : %6.2f%% : %-12s : %s : %s\n",
           domain, 100.0*score, id, frame, species);
}


/************************************************************************/
/*>void PrintAlignment(FILE *out, char *inAlign1, char *inAlign2)
   --------------------------------------------------------------
*//**
   \param[in]   *out       Output file pointer
   \param[in]   *inAlign1  Aligned sequence 1
   \param[in]   *inAlign2  Aligned sequence 2

   Displays the aligned region of the two sequences

-  14.04.20 Original   By: ACRM
-  14.06.21 Added printing of number of mismatches
*/
void PrintAlignment(FILE *out, char *inAlign1, char *inAlign2)
{
   char align1[HUGEBUFF+1],
        align2[HUGEBUFF+1],
        *aln1,
        *aln2;
   int  i,
        nMismatches = 0;
   
   strncpy(align1, inAlign1, HUGEBUFF);
   strncpy(align2, inAlign2, HUGEBUFF);
   aln1=align1+strlen(align1)-1;
   aln2=align2+strlen(align2)-1;
   while((*aln1 == '-')  || (*aln2 == '-'))
   {
      aln1--;
      aln2--;
   }
   *(aln1+1)='\0';
   *(aln2+1)='\0';

   aln1=align1;
   aln2=align2;
   while((*aln1 == '-')||(*aln2 == '-'))
   {
      aln1++;
      aln2++;
   }

   fprintf(out, "    %s\n", aln1);

   fprintf(out, "    ");
   for(i=0; i<strlen(aln1); i++)
   {
      if(aln1[i] == aln2[i])
      {
         fprintf(out, "|");
      }
      else
      {
         fprintf(out, " ");
         nMismatches++;
      }
   }
   fprintf(out,"\n");
   
   fprintf(out, "    %s\n", aln2);
   fprintf(out, "    Mismatches: %d\n\n", nMismatches);
}


/************************************************************************/
/*>void DoDSegment(FILE *out, char *seq, char *species, char *dataDir,
                   char *hvBestAlign1, char *hvBestAlign2, 
                   char *hjBestAlign1, char *hjBestAlign2,
                   BOOL verbose, BOOL showAlignment)
*//**
   \param[in]     *out           Output file pointer
   \param[in,out] *seq           Sequence to analyze
   \param[in]     *species       "Homo", "Mus" or blank
   \param[in]     *dataDir       Data directory
   \param[in]     *hvBestAlign1  Best VH alignment (in seq)
   \param[in]     *hvBestAlign2  Best VH alignment (in database)
   \param[in]     *hjBestAlign1  Best JH alignment (in seq)
   \param[in]     *hjBestAlign2  Best JH alignment (in database)
   \param[in]     verbose        Verbose output
   \param[in]     showAlignment  Show the alignment of the region

   Handle the identification of a D-segment. If we are not deleting
   regions, then this does so and extracts just the region between
   V and J since this is so short.

-  26.04.23 Original   By: ACRM   
*/
void DoDSegment(FILE *out, char *seq, char *species, char *dataDir,
                char *hvBestAlign1, char *hvBestAlign2, 
                char *hjBestAlign1, char *hjBestAlign2,
                BOOL verbose, BOOL showAlignment)
{
   char        hdMatch[MAXBUFF+1],
               DSeq[MAXBUFF+1];
   static char bestAlign1[HUGEBUFF+1],
               bestAlign2[HUGEBUFF+1];
   REAL        hdScore     = -1.0;
   
#ifndef REMOVESEQS
   RemoveSequence(seq, hvBestAlign1, hvBestAlign2, verbose);
   RemoveSequence(seq, hjBestAlign1, hjBestAlign2, verbose);
#endif

   CopyDSegment(DSeq, seq);
   
   hdScore = ScanAgainstDB("heavy_d", DSeq, verbose, species,
                           hdMatch, bestAlign1, bestAlign2,
                           dataDir);

   if(hdScore > THRESHOLD_HD)
   {
      int alnLength = CalcShortSeqLen(bestAlign1, bestAlign2);
      
      PrintResult(out, "DH", hdScore/(REAL)alnLength, hdMatch);
      if(showAlignment)
         PrintAlignment(out, bestAlign1, bestAlign2);
   }
}

/************************************************************************/
/*>void CopyDSegment(char *DSeq, char *seq)
   ----------------------------------------
*//**
   \param[out]   *DSeq    The region between the V and J segments
   \param[in]    *seq     Sequence to analyze

   Takes the sequence (seq) which has the V and J regions masked out
   with X characters and then extracts the region between them and
   copies it into DSeq.

-  26.04.23 Original   By: ACRM   
*/
void CopyDSegment(char *DSeq, char *seq)
{
   int inPos  = 0,
       outPos = 0,
       block  = BLOCK_UNDEFINED;
   
   for(inPos=0, outPos=0;  inPos<strlen(seq);  inPos++)
   {
      switch(block)
      {
      case BLOCK_UNDEFINED:
         if(seq[inPos] == 'X') block = BLOCK_X_V;
         break;
      case BLOCK_X_V:
         if(seq[inPos] != 'X')
         {
            block = BLOCK_D;
            DSeq[outPos++] = seq[inPos];
         }
         break;
      case BLOCK_D:
         if(seq[inPos] == 'X')
         {
            inPos = strlen(seq+1);
            break;
         }
         DSeq[outPos++] = seq[inPos];
      }
   }
   DSeq[outPos] = '\0';
}


/************************************************************************/
/*>BOOL PrintSpecialMatches(FILE *out, char *CH1, char *CH2, char *CH3)
   --------------------------------------------------------------------
*//**
     \param[in]   out    Output file pointer
     \param[in]   CH1    CH1 label
     \param[in]   CH2    CH2 label
     \param[in]   CH3    CH3-CHS label
     \return             True if this was a special case; false otherwise

     Looks for cases where apparent mixed allelic variants actually match
     a different allelic variant indentifier.

     - 22.02.24 Original   By: ACRM
*/
BOOL PrintSpecialMatches(FILE *out, char *CH1, char *CH2, char *CH3)
{
   CH1 = strchr(CH1, '_');
   CH1++;
   TERMAT(CH1, '_');
   
   CH2 = strchr(CH2, '_');
   CH2++;
   TERMAT(CH2, '_');
   
   CH3 = strchr(CH3, '_');
   CH3++;
   TERMAT(CH3, '_');
   
   if(!strncmp(CH1, "IGHG1*03", 8) &&
      !strncmp(CH2, "IGHG1*01", 8) &&
      !strncmp(CH3, "IGHG1*01", 8))
   {
      fprintf(out, "*CH1/2/3 matches IGHG1*08\n");
      return(TRUE);
   }
   else if (!strncmp(CH1, "IGHG1*01", 8) &&
            !strncmp(CH2, "IGHG1*01", 8) &&
            !strncmp(CH3, "IGHG1*03", 8))
   {
      fprintf(out, "*** CH1/2/3 matches IGHG1*15\n");
      return(TRUE);
   }
   else if (!strncmp(CH1, "IGHG1*03", 8) &&
            !strncmp(CH2, "IGHG1*01", 8) &&
            !strncmp(CH3, "IGHG1*03", 8))
   {
      fprintf(out, "*** CH1/2/3 matches IGHG1*03\n");
      return(TRUE);
   }
   else if (!strncmp(CH1, "IGHG1*01", 8) &&
            !strncmp(CH3, "IGHG1*04", 8))
   {
      fprintf(out, "*** CH1/2/3 matches IGHG1*04\n");
      return(TRUE);
   }
   else if (!strncmp(CH1, "IGHG1*01", 8) &&
            !strncmp(CH3, "IGHG1*07", 8))
   {
      fprintf(out, "*** CH1/2/3 matches IGHG1*07\n");
      return(TRUE);
   }
   
   return(FALSE);
}



#ifdef USEPATH
/************************************************************************/
char *FindPath(void)
{
  char *path = NULL;
  int  length, dirnameLength;

  if((length = wai_getExecutablePath(NULL, 0, &dirnameLength)) > 0)
  {
     if((path = (char*)malloc(length + 1))==NULL)
        return(NULL);
     
    wai_getExecutablePath(path, length, &dirnameLength);
    path[dirnameLength] = '\0';
    return(path);
  }

  return(NULL);
}

/************************************************************************/
#include <dirent.h>
#include <errno.h>
BOOL DirectoryExists(char *dirName)
{

   DIR* dir = opendir(dirName);
   if (dir)
   {
      /* Directory exists. */
      closedir(dir);
      return(TRUE);
   }

   return(FALSE);
}
#endif

