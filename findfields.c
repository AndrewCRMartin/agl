#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <assert.h>
#include "bioplib/macros.h"
#include "findfields.h"
#define SMALLBUFF 80


int test_findfields(void);
static void CopyMatch(char *out, int maxbuff, char *in, regmatch_t *pMatches, int matchNum);
static int IntMatch(char *in, regmatch_t *pMatches, int matchNum);

#ifdef TEST
int main(int argc, char **argv)
{
   return(test_findfields());
}
#endif

int test_findfields(void)
{
   char class[SMALLBUFF],
      distal[SMALLBUFF],
      family[SMALLBUFF],
      subclass[SMALLBUFF];
   int allele;
   
   FindFields("IGKV1D-16*01", class, subclass, family, &allele, distal);
   printf("IGKV1D-16*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);
   
   FindFields("IGKV1-16D*01", class, subclass, family, &allele, distal);
   printf("IGKV1-16D*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("IGHV1-NL1*01", class, subclass, family, &allele, distal);
   printf("IGHV1-NL1*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);
   
   FindFields("IGHV1-38-4*01", class, subclass, family, &allele, distal);
   printf("IGHV1-38-4*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("TRDV1S1*01", class, subclass, family, &allele, distal);
   printf("TRDV1S1*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("IGHV3-48*02", class, subclass, family, &allele, distal);
   printf("IGHV3-48*02 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("IGKV12-e*01", class, subclass, family, &allele, distal);
   printf("IGKV12-e*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("IGHG2*01", class, subclass, family, &allele, distal);
   printf("IGHG2*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("IGHA*01", class, subclass, family, &allele, distal);
   printf("IGHA*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("IGHD1/OR15-1a*01", class, subclass, family, &allele, distal);
   printf("IGHD1/OR15-1a*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("IGLJ-C/OR18*01", class, subclass, family, &allele, distal);
   printf("IGLJ-C/OR18*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("TRAV14/DV4*01", class, subclass, family, &allele, distal);
   printf("TRAV14/DV4*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("TRAV38-2/DV8*01", class, subclass, family, &allele, distal);
   printf("TRAV38-2/DV8*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("TRAV15-1/DV6-1*01", class, subclass, family, &allele, distal);
   printf("TRAV15-1/DV6-1*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);
   
   FindFields("TRAV15D-1/DV6D-1*01", class, subclass, family, &allele, distal);
   printf("TRAV15D-1/DV6D-1*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   FindFields("TRAV14D-3/DV8*01", class, subclass, family, &allele, distal);
   printf("TRAV14D-3/DV8*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);
   
   FindFields("TRAV16D/DV11*01", class, subclass, family, &allele, distal);
   printf("TRAV16D/DV11*01 Class: %s, Subclass: %s, Family: %s, Allele: %d, Distal: %s\n",
          class, subclass, family, allele, distal);

   return(0);
}

/************************************************************************/
static void CopyMatch(char *out, int maxbuff, char *in, regmatch_t *pMatches, int matchNum)
{
   regmatch_t *match = pMatches + matchNum;
   int matchLen, nCopy;
   out[0] = '\0';
   if(match->rm_so == (-1))
      return;
   matchLen = match->rm_eo - match->rm_so;
   nCopy = MIN(matchLen, maxbuff-1);
   strncpy(out, in+match->rm_so, nCopy);
   out[nCopy] = '\0';
}

/************************************************************************/
static int IntMatch(char *in, regmatch_t *pMatches, int matchNum)
{
   regmatch_t *match = pMatches + matchNum;
   int matchLen, nCopy;
   char buffer[SMALLBUFF];

   buffer[0] = '\0';
   if(match->rm_so == (-1))
      return(0);
   matchLen = match->rm_eo - match->rm_so;
   nCopy = MIN(matchLen, SMALLBUFF-1);
   strncpy(buffer, in+match->rm_so, nCopy);
   buffer[nCopy] = '\0';

   return(atoi(buffer));
}

/************************************************************************/
void FindFields(char *id, char *class, char *subclass, char *family, int *pAllele, char *distal)
{
   regex_t re;
   regmatch_t matches[16];
   int err;

   class[0]    = '\0';
   subclass[0] = '\0';
   distal[0]   = '\0';
   family[0]   = '\0';
   *pAllele    = 0;
   
   /* CCCCSSDD-FF*AA (e.g. IGKV1D-16*01) */   
   err = regcomp(&re, "^([A-Z]+)([0-9]+)([A-Z0-9]+)-([0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      CopyMatch(distal, SMALLBUFF, id, matches, 3);
      CopyMatch(family, SMALLBUFF, id, matches, 4);
      *pAllele   = IntMatch(id, matches, 5);
      regfree(&re);
      return;
   }

   /* CCCCSS-FFDD*AA (e.g. IGKV1-16D*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([0-9]+)([A-Z]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      CopyMatch(distal, SMALLBUFF, id, matches, 4);
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 5);
      regfree(&re);
      return;
   }

   /* CCCCSS-FF*AA (e.g. IGHV1-NL1*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([A-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0]  = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSS-FF-FF*AA (e.g. IGHV1-38-4*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([0-9]+-[0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0]  = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSDD*AA (e.g. TRDV1S1*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)([A-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0]  = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS-FF*AA (e.g. IGHV3-48*02) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0]  = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS-FF*AA (e.g. IGKV12-e*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)-([a-zA-Z]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0]  = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS*AA (e.g. IGHG2*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0]  = '\0';
      family[0]  = '\0';
      *pAllele   = IntMatch(id, matches, 3);
      regfree(&re);
      return;
   }

   /*  CCCC*AA (e.g. IGHA*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      subclass[0] = '\0';
      distal[0]  = '\0';
      family[0]  = '\0';
      *pAllele   = IntMatch(id, matches, 2);
      regfree(&re);
      return;
   }

   /* CCCCSS/DD-FF*AA (e.g. IGHD1/OR15-1a*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]*)/([A-Z0-9]+)-([a-zA-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      CopyMatch(distal, SMALLBUFF, id, matches, 3);
      CopyMatch(family, SMALLBUFF, id, matches, 4);
      *pAllele   = IntMatch(id, matches, 5);
      regfree(&re);
      return;
   }

   /*  CCCC-SS/FFFF*AA (e.g. IGLJ-C/OR18*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)-([A-Z]+)/([A-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS/FFF*AA (e.g. TRAV14/DV4*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)/([A-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*  CCCCSS-SS/FFF*AA (e.g. TRAV38-2/DV8*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+-[0-9]+)/([A-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /*   CCCCSS-SS/FFF-F*AA (e.g. TRAV15-1/DV6-1*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+-[0-9]+)/([A-Z0-9]+-[0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      distal[0] = '\0';
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSD-SS/FFF-F*AA (e.g. TRAV15D-1/DV6D-1*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+D-[0-9]+)/([A-Z0-9]+-[0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      strcpy(distal, "D");
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSD-SS/FFF*AA (e.g. TRAV14D-3/DV8*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+D-[0-9]+)/([A-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      strcpy(distal, "D");
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   /* CCCCSSD/FFF*AA (e.g. TRAV16D/DV11*01) */
   regfree(&re);
   err = regcomp(&re, "^([A-Z]+)([0-9]+)D/([A-Z0-9]+)\\*([0-9]+)", REG_EXTENDED);
   assert(err == 0);
   if(regexec(&re, id, sizeof(matches)/sizeof(matches[0]), (regmatch_t *)&matches, 0) == 0)
   {
      CopyMatch(class, SMALLBUFF, id, matches, 1);
      CopyMatch(subclass, SMALLBUFF, id, matches, 2);
      strcpy(distal, "D");
      CopyMatch(family, SMALLBUFF, id, matches, 3);
      *pAllele   = IntMatch(id, matches, 4);
      regfree(&re);
      return;
   }

   
   
   
   
}
