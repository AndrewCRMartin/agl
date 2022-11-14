/************************************************************************/
/**

   Program:    agl (Assign Germ Line)
   \file       findfields.c
   
   \version    V1.5
   \date       14.11.22
   \brief      Parse IMGT identifier into fields
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2022
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
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

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
   V1.5   14.11.22  Frees regex buffer if there were no matches

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <assert.h>
#include "bioplib/macros.h"
#include "findfields.h"
#include "agl.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static void CopyMatch(char *out, int maxbuff, char *in,
                      regmatch_t *pMatches, int matchNum);
static int IntMatch(char *in, regmatch_t *pMatches, int matchNum);

/************************************************************************/
#ifdef TEST
int test_findfields(void);
static void testPrint(char *id, char *class, char *subclass, char *family,
                      int allele, char *distal);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Test main program
   
   31.03.20 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   return(test_findfields());
}

/************************************************************************/
/*>static void testPrint(char *id, char *class, char *subclass, 
                         char *family, int allele, char *distal)
   ----------------------------------------------------------------------
*//**
   \param[in]  id        Domains identifier
   \param[in]  class     Class
   \param[in]  subclass  SubClass
   \param[in]  family    Family
   \param[in]  allele    Allele number
   \param[in]  distal    Distal flag

   Prints the extracted information for testing purposes
   
   31.03.20 Original   By: ACRM
*/
static void testPrint(char *id, char *class, char *subclass, char *family,
                      int allele, char *distal)
{
   printf("%s Class: %s, Subclass: %s, Family: %s, Allele: %d, \
Distal: %s\n",
          id, class, subclass, family, allele, distal);
}

/************************************************************************/
/*>int test_findfields(void)
   -------------------------
*//**
   Test harness for trying different identifier patterns

   - 31.03.20 Original   By: ACRM
*/
int test_findfields(void)
{
   char class[SMALLBUFF],
        distal[SMALLBUFF],
        family[SMALLBUFF],
        subclass[SMALLBUFF];
   int  allele;
   
   FindFields("IGKV1D-16*01",
              class, subclass, family, &allele, distal);
   testPrint("IGKV1D-16*01",
             class, subclass, family, allele, distal);
   
   FindFields("IGKV1-16D*01",
              class, subclass, family, &allele, distal);
   testPrint("IGKV1-16D*01",
             class, subclass, family, allele, distal);

   FindFields("IGHV1-NL1*01",
              class, subclass, family, &allele, distal);
   testPrint("IGHV1-NL1*01",
             class, subclass, family, allele, distal);
   
   FindFields("IGHV1-38-4*01",
              class, subclass, family, &allele, distal);
   testPrint("IGHV1-38-4*01",
             class, subclass, family, allele, distal);

   FindFields("TRDV1S1*01",
              class, subclass, family, &allele, distal);
   testPrint("TRDV1S1*01",
             class, subclass, family, allele, distal);

   FindFields("IGHV3-48*02",
              class, subclass, family, &allele, distal);
   testPrint("IGHV3-48*02",
             class, subclass, family, allele, distal);

   FindFields("IGKV12-e*01",
              class, subclass, family, &allele, distal);
   testPrint("IGKV12-e*01",
             class, subclass, family, allele, distal);

   FindFields("IGHG2*01",
              class, subclass, family, &allele, distal);
   testPrint("IGHG2*01",
             class, subclass, family, allele, distal);

   FindFields("IGHA*01",
              class, subclass, family, &allele, distal);
   testPrint("IGHA*01",
             class, subclass, family, allele, distal);

   FindFields("IGHD1/OR15-1a*01",
              class, subclass, family, &allele, distal);
   testPrint("IGHD1/OR15-1a*01",
             class, subclass, family, allele, distal);

   FindFields("IGLJ-C/OR18*01",
              class, subclass, family, &allele, distal);
   testPrint("IGLJ-C/OR18*01",
             class, subclass, family, allele, distal);

   FindFields("TRAV14/DV4*01",
              class, subclass, family, &allele, distal);
   testPrint("TRAV14/DV4*01",
             class, subclass, family, allele, distal);

   FindFields("TRAV38-2/DV8*01",
              class, subclass, family, &allele, distal);
   testPrint("TRAV38-2/DV8*01",
             class, subclass, family, allele, distal);

   FindFields("TRAV15-1/DV6-1*01",
              class, subclass, family, &allele, distal);
   testPrint("TRAV15-1/DV6-1*01",
             class, subclass, family, allele, distal);
   
   FindFields("TRAV15D-1/DV6D-1*01",
              class, subclass, family, &allele, distal);
   testPrint("TRAV15D-1/DV6D-1*01",
             class, subclass, family, allele, distal);

   FindFields("TRAV14D-3/DV8*01",
              class, subclass, family, &allele, distal);
   testPrint("TRAV14D-3/DV8*01",
             class, subclass, family, allele, distal);
   
   FindFields("TRAV16D/DV11*01",
              class, subclass, family, &allele, distal);
   testPrint("TRAV16D/DV11*01",
             class, subclass, family, allele, distal);

   printf("\n\nExtras\n\n");
   
   FindFields("IGHV4-4*08",
              class, subclass, family, &allele, distal);
   testPrint("IGHV4-4*08",
             class, subclass, family, allele, distal);

   FindFields("IGHV4-59*01",
              class, subclass, family, &allele, distal);
   testPrint("IGHV4-59*01",
             class, subclass, family, allele, distal);

   return(0);
}
#endif


/************************************************************************/
/*>static void CopyMatch(char *out, int maxbuff, char *in,
                         regmatch_t *pMatches, int matchNum)
   -------------------------------------------------------
*//**
   \param[out]  out        Output string
   \param[in]   maxbuff    Max size of output
   \param[in]   in         Input string
   \param[in]   pMatches   POSIX regex structure pointer
   \param[in]   matchNum   Match number we wish to extract

   Extracts the specified match from the input string

   - 31.03.20 Original   By: ACRM
*/
static void CopyMatch(char *out, int maxbuff, char *in,
                      regmatch_t *pMatches, int matchNum)
{
   regmatch_t *match = pMatches + matchNum;
   int        matchLen,
              nCopy;

   out[0] = '\0';
   if(match->rm_so == (-1))
      return;

   matchLen = match->rm_eo - match->rm_so;
   nCopy    = MIN(matchLen, maxbuff-1);
   strncpy(out, in+match->rm_so, nCopy);
   out[nCopy] = '\0';
}


/************************************************************************/
/*>static int IntMatch(char *in, regmatch_t *pMatches, int matchNum)
   -----------------------------------------------------------------
*//**
   \param[in]   in         Input string
   \param[in]   pMatches   POSIX regex structure pointer
   \param[in]   matchNum   Match number we wish to extract
   \return                 Integer extracted from string

   Extracts the specified integer match from the input string

   - 31.03.20 Original   By: ACRM
*/
static int IntMatch(char *in, regmatch_t *pMatches, int matchNum)
{
   regmatch_t *match = pMatches + matchNum;
   int         matchLen,
               nCopy;
   char        buffer[SMALLBUFF];

   buffer[0] = '\0';
   if(match->rm_so == (-1))
      return(0);

   matchLen = match->rm_eo - match->rm_so;
   nCopy    = MIN(matchLen, SMALLBUFF-1);
   strncpy(buffer, in+match->rm_so, nCopy);
   buffer[nCopy] = '\0';

   return(atoi(buffer));
}


/************************************************************************/
/*>void FindFields(char *id, char *class, char *subclass, char *family,
                   int *pAllele, char *distal)
   --------------------------------------------------------------------
*//**
   \param[in]  id         The IMGT gene ID
   \param[out] class      Class
   \param[out] subclass   Subclass
   \param[out] family     Family
   \param[out] pAllele    Allele number
   \param[out] distal     Distal flag

   Parses an IMGT gene ID to extract class, etc.

   - 31.03.20 Original   By: ACRM
*/
void FindFields(char *id, char *class, char *subclass, char *family,
                int *pAllele, char *distal)
{
   regex_t    re;
   regmatch_t matches[LABELBUFF];
   int        err;

   class[0]    = '\0';
   subclass[0] = '\0';
   distal[0]   = '\0';
   family[0]   = '\0';
   *pAllele    = 0;
   
   /* CCCCSSDD-FF*AA (e.g. IGKV1D-16*01)                                */
   err = regcomp(&re, "^([A-Z]+)([0-9]+)([A-Z0-9]+)-([0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      CopyMatch(distal,   SMALLBUFF, id, matches, 3);
      CopyMatch(family,   SMALLBUFF, id, matches, 4);
      *pAllele = IntMatch(id, matches, 5);
      regfree(&re);
      return;
   }

   /* CCCCSS-FFDD*AA (e.g. IGKV1-16D*01)                                */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([0-9]+)([A-Z]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      CopyMatch(distal,   SMALLBUFF, id, matches, 4);
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele = IntMatch(id, matches, 5);
      regfree(&re);
      return;
   }

   /* CCCCSS-FF*AA (e.g. IGHV1-NL1*01)                                  */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([A-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSS-FF-FF*AA (e.g. IGHV1-38-4*01)                              */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([0-9]+-[0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSDD*AA (e.g. TRDV1S1*01)                                     */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)([A-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS-FF*AA (e.g. IGHV3-48*02)                                  */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS-FF*AA (e.g. IGKV12-e*01)                                  */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([a-zA-Z]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS*AA (e.g. IGHG2*01)                                        */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      family[0] = '\0';
      *pAllele  = IntMatch(id, matches, 3);
      regfree(&re);
      return;
   }

   /*  CCCC*AA (e.g. IGHA*01)                                           */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      subclass[0] = '\0';
      distal[0]   = '\0';
      family[0]   = '\0';
      *pAllele    = IntMatch(id, matches, 2);
      regfree(&re);
      return;
   }

   /* CCCCSS/DD-FF*AA (e.g. IGHD1/OR15-1a*01)                           */
   regfree(&re);
   err = regcomp(&re,
                 "^([A-Z]+)([0-9]*)/([A-Z0-9]+)-([a-zA-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      CopyMatch(distal,   SMALLBUFF, id, matches, 3);
      CopyMatch(family,   SMALLBUFF, id, matches, 4);
      *pAllele = IntMatch(id, matches, 5);
      regfree(&re);
      return;
   }

   /*  CCCC-SS/FFFF*AA (e.g. IGLJ-C/OR18*01)                            */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)-([A-Z]+)/([A-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS/FFF*AA (e.g. TRAV14/DV4*01)                               */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)/([A-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS-SS/FFF*AA (e.g. TRAV38-2/DV8*01)                          */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+-[0-9]+)/([A-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSS-SS/FFF-F*AA (e.g. TRAV15-1/DV6-1*01)                       */
   regfree(&re);
   err = regcomp(&re,
                 "^([A-Z]+)([0-9]+-[0-9]+)/([A-Z0-9]+-[0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele  = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSD-SS/FFF-F*AA (e.g. TRAV15D-1/DV6D-1*01)                    */
   regfree(&re);
   err = regcomp(&re,
                 "^([A-Z]+)([0-9]+D-[0-9]+)/([A-Z0-9]+-[0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      strcpy(distal, "D");
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSD-SS/FFF*AA (e.g. TRAV14D-3/DV8*01)                         */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+D-[0-9]+)/([A-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      strcpy(distal, "D");
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSD/FFF*AA (e.g. TRAV16D/DV11*01)                             */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)D/([A-Z0-9]+)\\*([0-9]+)",
                 REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]),
              (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class,    SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      strcpy(distal, "D");
      CopyMatch(family,   SMALLBUFF, id, matches, 3);
      *pAllele = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }
   regfree(&re);
}
