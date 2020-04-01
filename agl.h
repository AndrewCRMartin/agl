/************************************************************************/
/**

   Program:    agl (Assign Germ Line)
   \file       agl.h
   
   \version    V1.0    
   \date       31.03.20   
   \brief      Assigns IMGT germline
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2020
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

*************************************************************************/
/* Includes
*/
#define SMALLBUFF         80
#define MAXBUFF           256
#define HUGEBUFF          1024
#define CHAINTYPE_UNKNOWN 0
#define CHAINTYPE_LIGHT   1
#define CHAINTYPE_HEAVY   2
#define AGLDATADIR        "AGLDATADIR"
#define THRESHOLD_LV      0.5
#define THRESHOLD_LJ      0.5
#define THRESHOLD_LC      0.5
#define THRESHOLD_HV      0.5
#define THRESHOLD_HJ      0.5
#define THRESHOLD_HC      0.5
