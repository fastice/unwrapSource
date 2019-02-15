#include "stdio.h"
#include"string.h"
#include "unwrapInclude.h"
#include <math.h>
#include <stdlib.h>

    static void initLabel( unwrapImageStructure *unwrapImage );
    static void labelPass1( unwrapImageStructure *unwrapImage );
    static void labelPass2( unwrapImageStructure *unwrapImage );  

    static fpNeighbors firstPassNeighborsL( unwrapImageStructure *unwrapImage,
                                    int i, int j);
    static spNeighbors secondPassNeighborsL( unwrapImageStructure *unwrapImage,
                                     int i, int j);
    static void forwardPropagateL(unwrapImageStructure *unwrapImage,
                         fpNeighbors fpN,int row,int column,int *newLabel);

    static void addFPTableL(unwrapImageStructure *unwrapImage,fpNeighbors fpN, 
                           int row,int column,eqTableElement *eqTable, int *nEQ);
    static void addSPTableL(unwrapImageStructure *unwrapImage,spNeighbors spN, 
                           int row,int column,eqTableElement *eqTable, int *nEQ);

    static void resolveEQL(eqTableElement *eqTable, int nEQ,
                          eqTableElement *labelList,int *nLabels);
    static void dfsL(int node,eqTableElement *eqTable, int nEQ, 
                    eqTableElement *labelList,int nLabels,
                    nodeListElement *nodesFound,int *nFound);    
    static void closeNodeL(int node, eqTableElement *labelList, int nLabels );
    static int nodeClosedL(int node, eqTableElement *labelList, int nLabels );
    static int equivalentLabelL(int label, eqTableElement *labelList,int nLabels);

 
    void labelRegions( unwrapImageStructure *unwrapImage )
{
    int maxlabel,i,j;
    initLabel( unwrapImage );         /* Init variables for labelint */
    labelPass1( unwrapImage );      /* First pass, top-left to bottom right*/
    labelPass2( unwrapImage );      /* Second pass, bottom-right to top-left */
    return;
}


    static void initLabel( unwrapImageStructure *unwrapImage )
{
    int i, j;
    for( i=0; i < unwrapImage->azimuthSize; i++) {
        for( j=0; j < unwrapImage->rangeSize; j++) { /* Compute differences */
            unwrapImage->labels[i][j] = LARGEINT;
        }
    }
    return;
}




    static void labelPass1( unwrapImageStructure *unwrapImage )
{
    fpNeighbors fpN;
    int newLabel = 1;
    int i,j;
    eqTableElement *eqList, *eqTable;
    int nLabels, nEQ;
    int eqIndex;
    nLabels=0;
/*
   malloc mem for equiv table.
*/
    eqTable = (eqTableElement *)
                      malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
    eqList =  (eqTableElement *) 
                      malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
/*
   First pass
*/
    for( i=0; i < unwrapImage->azimuthSize; i++) {
        nEQ = 0;
        for( j=0; j < unwrapImage->rangeSize; j++) { /* Label row */ 
            if( (unwrapImage->bCuts[i][j] & REGION) &&
                validPixel(unwrapImage->bCuts[i][j]) ) { /* If not branch cut */
/*
    Get neighbors 
*/
                fpN = firstPassNeighborsL( unwrapImage,i,j ); 
/*
    Propagate label and phase
*/
                forwardPropagateL(unwrapImage,fpN,i,j,&newLabel);
/*
    Update equivalence table 
*/
                addFPTableL(unwrapImage,fpN,i,j,eqTable,&nEQ ); 
             }
        }
/*
   Pixels in row have now all been labeled. Now get equivalence classes for    
   labels and relabel to remove equivalent labels.
*/
       resolveEQL(eqTable,nEQ,eqList,&nLabels);
/*
   Relabel
*/
       for( j=0; j < unwrapImage->rangeSize-1; j++) {
           eqIndex = equivalentLabelL(unwrapImage->labels[i][j],eqList,nLabels);
           if(eqIndex >= 0) 
               unwrapImage->labels[i][j] = eqList[eqIndex].b;              
       }
    }
/*
    Free used memory and return
*/ 
    free(eqTable);
    free(eqList);
    return;
}



    static void labelPass2( unwrapImageStructure *unwrapImage )
{
    spNeighbors spN;
    int newLabel = 1;
    int i,j;
    eqTableElement *eqList;
    int nLabels;
    eqTableElement *eqTable;
    int nEQ;
    int eqIndex;
/*
   malloc mem for equiv table.
*/
    eqTable = (eqTableElement *)
                      malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
    eqList =  (eqTableElement *) 
                      malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
/*
   Do second pass
*/
    for( i=unwrapImage->azimuthSize-1; i >= 0; i--) {
        nEQ = 0;
        for( j=unwrapImage->rangeSize-1; j >= 0; j--) { /* Label row */ 
            if( (unwrapImage->bCuts[i][j] & REGION) &&
                validPixel(unwrapImage->bCuts[i][j]) ) { /* If not branch cut */
/*
    Get neighbors 
*/
                spN = secondPassNeighborsL( unwrapImage,i,j ); 
/*
    Update equivalence table 
*/
                addSPTableL(unwrapImage,spN,i,j,eqTable,&nEQ ); 
             }
        }
/*
   Pixels in row have now all been labeled. Now get equivalence classes for    
   labels and relabel to remove equivalent labels.
*/
       resolveEQL(eqTable,nEQ,eqList,&nLabels);
/*
   Relabel
*/
       for( j=0; j < unwrapImage->rangeSize; j++) {
           eqIndex = equivalentLabelL(unwrapImage->labels[i][j],eqList,nLabels);
           if(eqIndex >= 0) 
               unwrapImage->labels[i][j] = eqList[eqIndex].b;  
       }
    }
/*
    Free used memory and return
*/ 
    free(eqTable);
    free(eqList);
    return;
}

  
    static fpNeighbors firstPassNeighborsL( unwrapImageStructure *unwrapImage,
                                           int i, int j)
{
/*
    Return structure with upper and left labels.
*/
    fpNeighbors fpN;
    fpN.left = LARGEINT;
    fpN.upper = LARGEINT;

    if( i != 0 && (unwrapImage->bCuts[i-1][j] & REGION) 
               && validPixel(unwrapImage->bCuts[i-1][j] ) ) 
        fpN.upper = unwrapImage->labels[i-1][j];
    if( j != 0 && ( unwrapImage->bCuts[i][j-1] & REGION )
               && validPixel(unwrapImage->bCuts[i][j-1] ) )
        fpN.left = unwrapImage->labels[i][j-1];
    return fpN;
}


    static spNeighbors secondPassNeighborsL( unwrapImageStructure *unwrapImage,
                                            int i, int j)
{
/*
    Return structure with lower and right labels.
*/
    spNeighbors spN;
    spN.right = LARGEINT;
    spN.lower = LARGEINT;

    if(i==unwrapImage->azimuthSize-1) return spN; 
    if( (i < (unwrapImage->azimuthSize-1)) && (unwrapImage->bCuts[i+1][j] & REGION )
               && validPixel(unwrapImage->bCuts[i+1][j] ) ) 
        spN.lower = unwrapImage->labels[i+1][j];
    if(j==unwrapImage->rangeSize-1) return spN;
    if( (j < (unwrapImage->rangeSize-1)) && ( unwrapImage->bCuts[i][j+1] & REGION)
               && validPixel(unwrapImage->bCuts[i][j+1] ) )
        spN.right = unwrapImage->labels[i][j+1];

    return spN;
}


    static void forwardPropagateL(unwrapImageStructure *unwrapImage,
                                fpNeighbors fpN,int row,int column,int *newLabel)
{
    if( noFPNeighbors(fpN) == TRUE) { /* New label */
	unwrapImage->labels[row][column] =(*newLabel)++;
	return;
    }
    if(fpN.left <= fpN.upper)
        unwrapImage->labels[row][column] = fpN.left;
    else 
        unwrapImage->labels[row][column] = fpN.upper;
}


    static void addFPTableL(unwrapImageStructure *unwrapImage,fpNeighbors fpN,
                           int row,int column,eqTableElement *eqTable, int *nEQ )
/*
   If new pair, add to table.
*/
{
   int i, new;
   int a,b;
   int currentLabel;

   eqTable[*nEQ].a = -1;
   eqTable[*nEQ].b = -1;
   currentLabel = unwrapImage->labels[row][column];
   new = FALSE;
/*
   Determine if upper or left is possible equivalency.
*/
   if( (fpN.left < LARGEINT) && (currentLabel < LARGEINT) && 
       (currentLabel != fpN.left) ) {    
       a = currentLabel;  b = fpN.left;      
       new = TRUE; 
   } else if( (fpN.upper < LARGEINT) && (currentLabel < LARGEINT) && 
              (currentLabel != fpN.upper) ) {      
       a = currentLabel;  b = fpN.upper;      
       new = TRUE; 
   }
/*
   Determine if equivalency is already in the table.
*/
   i = 0;
   while( (new  == TRUE) && (i < *nEQ) ) { /* Check if already in table */
      if( (a == eqTable[i].a) && (b == eqTable[i].b) ) new = FALSE;
      i++;
   }
/*
   If it is new, add to the table.
*/
   if( new == TRUE) { 
       eqTable[*nEQ].a = a;  eqTable[*nEQ].b = b;
       (*nEQ)++; 
   }  
   return;
}


   static void addSPTableL(unwrapImageStructure *unwrapImage,spNeighbors spN,
                          int row,int column,eqTableElement *eqTable, int *nEQ )
{
   int i, new;
   int a,b;
   int currentLabel;

   eqTable[*nEQ].a = -1;
   eqTable[*nEQ].b = -1;
   currentLabel = unwrapImage->labels[row][column];
   new = FALSE;
/*
   Determine if upper or left is possible equivalency.
*/
   if( (spN.right < LARGEINT) && (currentLabel < LARGEINT) && 
       (currentLabel != spN.right) ) {  
       b = currentLabel;  a = spN.right;      
       new = TRUE; 
   } else if( (spN.lower < LARGEINT) && (currentLabel < LARGEINT) && 
              (currentLabel != spN.lower) ) {     
       b = currentLabel;  a = spN.lower;      
       new = TRUE; 
   }
/*
   Determine if equivalency is already in the table.
*/
   i = 0;
   while( (new  == TRUE) && (i < *nEQ) ) { /* Check if already in table */
      if( (a == eqTable[i].a) && (b == eqTable[i].b) ) new = FALSE;
      i++;
   }
/*
   If it is new, add to the table.
*/
   if( new == TRUE) { 
       eqTable[*nEQ].a = a;  eqTable[*nEQ].b = b;
       (*nEQ)++; 
   }  
   return;
}


    static void resolveEQL(eqTableElement *eqTable, int nEQ,
                          eqTableElement *labelList,int *nLabels)
{
    int i, j, k;
    int node, nFound;
    nodeListElement labelsFound[10000];
    int a,b, minLabel;
  
    for(i = 0; i < 2 * nEQ; i++) labelList[i].a = -1;
    *nLabels = 0;
    for(i = 0; i < nEQ; i++) {     /* Build list of labels */
       a = eqTable[i].a;
       b = eqTable[i].b;
       for(j = 0; j < *nLabels; j++) {
           if(labelList[j].a == a) a = -1;
           if(labelList[j].a == b) b = -1;
       }
       if(a > 0) { 
           labelList[*nLabels].a = a; 
           labelList[*nLabels].b = -1;
           (*nLabels)++;
       }
       if(b > 0) { 
           labelList[*nLabels].a = b; 
           labelList[*nLabels].b = -1;
           (*nLabels)++; 
       }
    }    
    for( i = 0; i < *nLabels; i++ ) {  /* Resolve labels */
        nFound = 0;
        minLabel = LARGEINT;
        for(k = 0; k < *nLabels; k++) /* Find min open node */
            if(labelList[k].b < 0) minLabel = min(minLabel,labelList[k].a);
        node = minLabel;
        if(node == LARGEINT) return; /* Return if all nodes closed */

	dfsL(node,eqTable,nEQ,labelList,*nLabels,
		      labelsFound,&nFound); /*Search node */

	for(j = 0; j < nFound; j++) minLabel = min(minLabel,labelList[j].a);
	for(j = 0; j < nFound; j++) /* Mark labels with equivalent labels */
	    for(k = 0; k < *nLabels; k++) {
		if(labelList[k].a == labelsFound[j].node) {
		    labelList[k].b = labelsFound[0].node;
		}
	}
    }
}



    static void dfsL(int node,eqTableElement *eqTable, int nEQ, 
                    eqTableElement *labelList,int nLabels,
                    nodeListElement *nodesFound,int *nFound)
{
    int i;
    int newNode;
    float newDistance;
    if( nodeClosedL(node,labelList,nLabels) == TRUE ) 
        return; /* Node closed return */

    closeNodeL( node,labelList,nLabels );   /* Close node */
    nodesFound[*nFound].node = node;       /* Add node to found list */
    (*nFound)++;
/*
   Search node
*/
    for(i = 0; i < nEQ; i++ ) {  /* Search nodes */
        newNode = -1;
        if( node == eqTable[i].a ) { 
            newNode = eqTable[i].b; 
        } else if( node == eqTable[i].b ) {
            newNode = eqTable[i].a;
        }
        if(newNode > -1)
            dfsL(newNode,eqTable,nEQ,labelList,nLabels,
               nodesFound,nFound);
    }
}



    static void closeNodeL(int node, eqTableElement *labelList, int nLabels )
{
    int i;
    for( i = 0; i < nLabels; i++ ) 
        if( node == labelList[i].a ) {
            labelList[i].b = 0; 
            return;
        }
    error("*** closeNode: Invalid node ***\n");
}


    static int nodeClosedL(int node, eqTableElement *labelList, int nLabels )
{
    int i;
    for( i = 0; i < nLabels; i++ ) 
        if( node == labelList[i].a ) { 
            if(labelList[i].b > -1) return TRUE;
            else return FALSE;        
        }
    error("*** nodeClosed: Invalid node ***\n");
}




    static int equivalentLabelL(int currentlabel,eqTableElement *labelList,
                               int nLabels)
{
/* 
   Called by unwrapPass1 & 2. Find label in table and return equivalent label. 
*/
    int i;   
    for( i = 0; i < nLabels; i++ ) /* Find label and return equivalent */
        if( currentlabel == labelList[i].a ) return i;
    return -1;  /* Label unique */
}

