#include "unwrapInclude.h"
#include <math.h>
#include <stdlib.h>
/*
   Draw branch cuts for phase unwrapping.
*/
   static void residueDensity(unwrapImageStructure *unwrapImage);
    static void lowDensity( unwrapImageStructure *unwrapImage );
    static void balanceResidues( unwrapImageStructure *unwrapImage );
    static void dilateResidues( unwrapImageStructure *unwrapImage );
    static regionList *regionData( unwrapImageStructure *unwrapImage );
/****/
    static void maskRegionPoints(regionList *regions,
                                 unwrapImageStructure *unwrapImage);
/****/
static void balanceRegionsPass1(regionList *regions, regionList *topRegion, 
                                unwrapImageStructure *unwrapImage);
    static void subsumeRegions(regionList *topRegion,
                               regionList *region1,regionList *region2);
    static void balanceRegionsPass2( regionList *regions,regionList *topRegion, 
                              unwrapImageStructure *unwrapImage,int netCharge );
    static void regionsDrawLines( regionList *regions, 
                                  unwrapImageStructure *unwrapImage );
    static void freeRegions( regionList *regions );

    static void cancelBox(unwrapImageStructure *unwrapImage,
                          int r,int a,int boxSize);
    static void cancelResidues(unwrapImageStructure *unwrapImage,int nMatch,
                              int nMax,int aResidues[][2],int bResidues[][2]);
    static regionList *buildRegionList( regionList *regions, 
                              unwrapImageStructure *unwrapImage,int r,int a);
    static void addRegionPoints( regionList *regions, 
                              unwrapImageStructure *unwrapImage,int r,int a);
    static regionList *nearestRegion(regionList *regions1, 
                                     regionList *regions2, 
                                     int *rMin1,int *aMin1,
                                     int *rMin2,int *aMin2, int *dMin,int type);
    static int distanceBetweenRegions(regionList *regions1, 
                                      regionList *regions2, 
                                      int *rMin1,int *aMin1,
                                      int *rMin2,int *aMin2 );
    static int distanceToBorder( regionList *regions,
                                 unwrapImageStructure *unwrapImage,
                                 int *rMin, int *aMin, int *rB,int *aB );
    static void connectRegionPts( regionList *regions, 
                                  unwrapImageStructure *unwrapImage,
                                  int *charged1,int pt,int dMax);
    static void sortRegions( regionList *regions);
    static void sortRegionPts(int n, rList *ra);
    static int comparePts(rList r1,rList r2);
    static void drawBranchCut( unwrapImageStructure *unwrapImage,
		              int startX, int startY, int endX, int endY );
    static void connectLargeRegionPts( regionList *regions, 
                                       unwrapImageStructure *unwrapImage);
static void printLabels(regionList *regions)
{
/*
   Debug routine
*/
    fprintf(stderr,"%i %i\n",regions->label,regions->nTotal);
    if(regions->next != NULL) printLabels(regions->next);
    return;
}



   void branchCuts( unwrapImageStructure *unwrapImage )
{
    int netCharge = 0;
    int i,j;
    regionList *regions;
    FILE *fp;
/*
    Compute residue density
*/
    residueDensity( unwrapImage );
/*
    Connect low density residues
*/
    lowDensity( unwrapImage );
    fprintf(stderr,"nLeft after lowDensity %i \n",unwrapImage->nLeft);
/*
    Balance number of residues
*/
    balanceResidues( unwrapImage ); 
    fprintf(stderr,"nLeft after balanceResidues %i \n",unwrapImage->nLeft);
/*
    Dilate residues.
*/
    dilateResidues( unwrapImage );
/*
    Region label
*/
    labelRegions( unwrapImage );
    fprintf(stderr,"After final pass %i\n",unwrapImage->nLeft);
/*
   Region data.
*/
   regions = regionData( unwrapImage ); 
   fprintf(stderr,"After regionData\n");
/*
   Connect residues within a region
*/
    regionsDrawLines( regions, unwrapImage );
    fprintf(stderr,"After regionsDrawLines\n");  
/*
   balance regions, pass 1
*/
   balanceRegionsPass1( regions, regions,unwrapImage );
   fprintf(stderr,"After balanceRegionsPass1\n");
/*
   Pass 2
*/
   balanceRegionsPass2( regions, regions,unwrapImage, netCharge );
   fprintf(stderr,"After balanceRegionsPass2\n");
/*
   Fill in mask region as branch cuts
*/
   if(unwrapImage->mask != NULL) 
       for( i=0; i < unwrapImage->azimuthSize; i++)
           for( j = 0; j < unwrapImage->rangeSize; j++ ) {
                if( unwrapImage->bCuts[i][j] & MASK) {
                    unwrapImage->bCuts[i][j] |= BRANCHCUT;
                }
           }
   if(unwrapImage->cutFile != NULL) {
      fp = fopen(unwrapImage->cutFile,"w");
      for( i=0; i < unwrapImage->azimuthSize; i++) 
          fwriteBS(unwrapImage->bCuts[i],1,unwrapImage->rangeSize, fp,BYTEFLAG);

   }
   
/*
   freeRegion memory
*/

   freeRegions( regions );
   return;
}



/*
    Divide residues into high density and low density residues.
*/
    static void residueDensity(unwrapImageStructure *unwrapImage)
{
     int i,j,a,r;     
     int count=0;
     int iMin,iMax,jMin,jMax;
     int boxSize = 5;
     for( a=0; a < unwrapImage->azimuthSize-1; a++)
         for( r = 0; r < unwrapImage->rangeSize-1; r++ ) { 
	     if(r - boxSize < 0) jMin = 0; else jMin = -boxSize;
	     if(r+boxSize > unwrapImage->rangeSize-2) 
		 jMax= unwrapImage->rangeSize-2 - r;
	     else jMax=boxSize;
	     if(a-boxSize < 0) iMin = 0; else iMin = -boxSize;
	     if(a+boxSize > unwrapImage->azimuthSize-2) 
	         iMax= unwrapImage->azimuthSize-2-a;
	     else iMax=boxSize;
             if(unwrapImage->bCuts[a][r] & RESIDUE) {
                 count = 0;
		 for(i=iMin; i <= iMax; i++)
		     for(j=jMin; j <= jMax; j++) 
			 if( (unwrapImage->bCuts[a+i][r+j] & RESIDUE) ) count++;
                 if(count <= unwrapImage->residueDensity)
                     unwrapImage->bCuts[a][r] |= LOWDENSITY;
                 else unwrapImage->bCuts[a][r] &= HIGHDENSITY;
             }
         }
 }



/*
    Connect low density residues where possible
*/
   static void lowDensity( unwrapImageStructure *unwrapImage )
{
    int i,j;
    int charge, boxSize;
    int boxInc = 2;
/* 
  Connect low density residues 
*/
    for(boxSize = 2; boxSize < 12; boxSize+=boxInc) {
        for( i=0; i < unwrapImage->azimuthSize-1; i++)
             for( j = 0; j < unwrapImage->rangeSize-1; j++ ) 
                if( (unwrapImage->bCuts[i][j] & RESIDUE) &&
                    (unwrapImage->bCuts[i][j] & LOWDENSITY) &&
                     charged(unwrapImage->bCuts[i][j]) ) 
                     cancelBox(unwrapImage,j, i, boxSize);
         fprintf(stderr,"nLeft = %i\n",unwrapImage->nLeft);
         if(unwrapImage->nLeft == 0) return; /* Done ? */
    }   
    return;
}


/*
    Balance surplus residues so nPositive = nMinus. Surplus residues are 
    connected to the border with a branch cut.
*/
   static void balanceResidues( unwrapImageStructure *unwrapImage )
{
     typedef int (*coordinates)[3];
     int netResidues, nResidues; 
     int residue;
     coordinates coord;
     int d, minD;
     int i, j, r;
     int i1,i2,j1,j2,iMin,jMin,cMin;
     netResidues = unwrapImage->nPlus - unwrapImage->nMinus;
     if(netResidues < 0) residue = NRESIDUE;
     else if(netResidues > 0) residue = PRESIDUE;    
     else return;     

     nResidues = (unwrapImage->nLeft - abs(netResidues))/2 + abs(netResidues);
     coord = (coordinates) malloc(3 * nResidues * sizeof(int));    
/*
   Get coordinates of all residues of the unbalanced type.
*/
     r = 0;
     for( i=0; i < unwrapImage->azimuthSize-1; i++)
         for( j = 0; j < unwrapImage->rangeSize-1; j++ ) { 
                if( (unwrapImage->bCuts[i][j] & residue) && 
/*                    (unwrapImage->bCuts[i][j] & LOWDENSITY) && */
                    charged(unwrapImage->bCuts[i][j]) ) {
                 d = min(min(i,j),
                     min(unwrapImage->azimuthSize-i,unwrapImage->rangeSize-j));
                 coord[r][0] = i;
                 coord[r][1] = j;
                 coord[r][2] = d;
                 r++;
             }
         }
     if( r != nResidues ) 
         error("balanceResidues 1 %i  %i %i",r,nResidues,residue);

/*
    Now search for residues with minimum distance and draw branch cut off 
    the image.
*/

     for( r=0; r < abs(netResidues); r++ ) {
          minD = 30000;
          for(i = 0; i <  nResidues; i++ ) { /* Find min */
              if( coord[i][2] < minD ) {
                  minD = coord[i][2];
                  iMin = coord[i][0];
                  jMin = coord[i][1];
                  cMin = i;
              }
          }
          coord[cMin][2] = 30000;

          if( iMin == minD ) { /* Compute coordinates of branch cut */
               i1 = 0; i2 = iMin; j1 = jMin; j2 = jMin;
          } else if( (unwrapImage->azimuthSize-iMin) == minD) {
               i1 = iMin; 
               i2 = unwrapImage->azimuthSize-1; 
               j1 = jMin; j2 = jMin;
          } else if( jMin == minD ) {
               i1 = iMin; i2 = iMin; j1 = 0; j2 = jMin;
          } else if( (unwrapImage->rangeSize-jMin) == minD) {
               j1 = jMin; 
               j2 = unwrapImage->rangeSize-1; 
               i1 = iMin; i2 = iMin;
          } else error("balanceResidues 2 %i %i %i\n",iMin,jMin,minD);
         drawBranchCut( unwrapImage,j1, i1, j2, i2 );
         unwrapImage->bCuts[iMin][jMin] |= UNCHARGED;     
     }

     unwrapImage->nLeft -= abs(netResidues);
fprintf(stderr,"nLeft, nBalanced %i  %i\n",unwrapImage->nLeft,netResidues);
     free(coord);
     return;  
} 


/*
    Morphological dilation to connect regions of high density residues
*/
  static void dilateResidues( unwrapImageStructure *unwrapImage )
{
     int i,j,a,r;     
     int count=0;
     int iMin,iMax,jMin,jMax;
     int boxLeft, boxRight;
     int blockSize;
     blockSize = unwrapImage->blockSize;
     boxRight = (int)(blockSize/2);
     boxLeft = boxRight;
     if( (blockSize % 2) == 0 ) boxLeft--;
/*
     for( a=0; a < unwrapImage->azimuthSize-1; a++)
         for( r = 0; r < unwrapImage->rangeSize-1; r++ ) { 
*/
     for( a=0; a < unwrapImage->azimuthSize; a++)
         for( r = 0; r < unwrapImage->rangeSize; r++ ) { 

/***/
             if( (unwrapImage->mask !=NULL) &&
                 (unwrapImage->bCuts[a][r] & MASK) )
                 unwrapImage->bCuts[a][r] |= REGION;      
/***/
             if( (unwrapImage->bCuts[a][r] & RESIDUE) &&
                 charged(unwrapImage->bCuts[a][r]) ) {

		 if(r - boxLeft < 0) jMin = 0; else jMin = -boxLeft;
		 if(r+boxRight > unwrapImage->rangeSize-2) 
		     jMax= unwrapImage->rangeSize-2 - r;
		 else jMax=boxRight;
		 if(a-boxLeft < 0) iMin = 0; else iMin = -boxLeft;
		 if(a+boxRight > unwrapImage->azimuthSize-2) 
		     iMax= unwrapImage->azimuthSize-2-a;
		 else iMax=boxRight;
		 for(i=iMin; i <= iMax; i++)
		     for(j=jMin; j <= jMax; j++) 
			 unwrapImage->bCuts[a+i][r+j] |= REGION;
             }
         }
}



/*
    Build list of regions. Each element of the list contains information  
    describing corresponding region.
*/
    static regionList *regionData( unwrapImageStructure *unwrapImage )
{
     int i,j, c;
     int done=FALSE;
     regionList *regions,*tmp,first,*tmp1;
     int r1,a1,r2,a2;
     int maxLabel;
     regions = NULL;
/*
   Build list for masked points. This avoids larged masked regions getting
   stuck at the end of the list
*/
     regions = &first; /* this step is to avoid large label at start */
     first.next=NULL;
     first.label = -1;
     if(unwrapImage->mask != NULL) {
         for( i=unwrapImage->azimuthSize-1; i >= 0; i--) {
            if(i % 1000 == 0)  fprintf(stderr,"%i\n",i);
            for( j = unwrapImage->rangeSize -1; j >= 0; j-- ) {
                 if(unwrapImage->bCuts[i][j] & MASK)
                     regions = buildRegionList(regions,unwrapImage,j,i);
            } 
         }
     }
/*
  Now for residue points making sure not to do mask points again.
*/
     for( i=unwrapImage->azimuthSize-1; i >= 0; i--) {
        if(i % 1000 == 0)  fprintf(stderr,"%i\n",i);
        for( j = unwrapImage->rangeSize -1; j >= 0; j-- ) {
             if( ( (unwrapImage->bCuts[i][j] & RESIDUE) &&
                    charged(unwrapImage->bCuts[i][j]) ) &&
                    ((unwrapImage->bCuts[i][j] & MASK) == 0) )
                      regions = buildRegionList(regions,unwrapImage,j,i);
        } 
     }
     regions = regions->next;/*this step is to avoid large label at start */
fprintf(stderr,"After pass 1 of regionData\n");
     for( i=0; i < unwrapImage->azimuthSize; i++)
        for( j = 0; j < unwrapImage->rangeSize; j++ ) 
             if( (unwrapImage->bCuts[i][j] & RESIDUE) &&
                  charged(unwrapImage->bCuts[i][j])  )                
                      addRegionPoints(regions,unwrapImage,j,i);
fprintf(stderr,"After pass 2 of regionData\n");
/****/
if(unwrapImage->mask != NULL) fprintf(stderr,"MASK\n");
    if( unwrapImage->mask != NULL) maskRegionPoints(regions,unwrapImage);
fprintf(stderr,"After maskRegionPoints of regionData\n");
/****/
    sortRegions(regions);     
fprintf(stderr,"After sortRegions of regionData\n");
/*
    precompute some distances
*/
    tmp = regions;
    while (tmp->next != NULL) tmp = tmp->next;
    maxLabel = tmp->label;
    fprintf(stderr,"Max label = %i \n",maxLabel);
 
    tmp = regions;
    while(tmp != NULL) {
        if(tmp->nTotal > 500) {
             tmp->distances = (dist *) malloc((maxLabel+1)*sizeof(dist));
             tmp1 = regions;
             fprintf(stderr,"nTotal,label %i %i\n",tmp->nTotal,tmp->label);
             for(tmp1=regions; tmp1 != NULL;  tmp1=tmp1->next) { 
                 (tmp->distances[tmp1->label]).d = -1;
                 (tmp->distances[tmp1->label]).d = 
                    distanceBetweenRegions(tmp,tmp1,&r1,&a1,&r2,&a2); 
                 (tmp->distances[tmp1->label]).r1 =r1;
                 (tmp->distances[tmp1->label]).r2 =r2;
                 (tmp->distances[tmp1->label]).a1 =a1;
                 (tmp->distances[tmp1->label]).a2 =a2;
             }
        }
        tmp = tmp->next;
    }
    return regions;           
}


/*****/

/*
   Search through region list and find mask regions. Let borderpixels of 
   mask region act as list of points.
*/

    static void maskRegionPoints(regionList *regions,
                                 unwrapImageStructure *unwrapImage)

{
     int i,j;  
     int label,nTotal,n;

     if( regions == NULL) return;
     else if(regions->mask == TRUE) {
/*
    Pass one count border points
*/      
         nTotal = 0;
         label =regions->label;
/* fprintf(stderr,"%i",label); */

         for( i=1; i < unwrapImage->azimuthSize-1; i++)
             for( j = 1; j < unwrapImage->rangeSize-1; j++ ) {
                 unwrapImage->bCuts[i][j] &= NOTBORDER;
                 if(unwrapImage->labels[i][j] == label) {
/*
   Add area of same region as mask to mask
*/
                    unwrapImage->bCuts[i][j] |= MASK;
                    if( unwrapImage->labels[i+1][j-1] !=label ||
                        unwrapImage->labels[i+1][j+1] !=label ||
                        unwrapImage->labels[i-1][j-1] !=label ||
                        unwrapImage->labels[i-1][j+1] !=label  ) {
                        nTotal++;
                        unwrapImage->bCuts[i][j] |=BORDER;
                     } /* End if unwrapimage->labels[i+1]... */
                 } /* End if unwrapimage->labels[i] */
         }
         regions->pts = (rList *)malloc((nTotal+1) * sizeof(rList));
         regions->pts[nTotal].r = nTotal;
         regions->nTotal = nTotal;

/*
    Pass two add points
*/
         n = 0;
         for( i=1; i < unwrapImage->azimuthSize-1; i++)
             for( j = 1; j < unwrapImage->rangeSize-1; j++ ) 
                 if( (unwrapImage->labels[i][j] == label) &&
                     (unwrapImage->bCuts[i][j] & BORDER ) ) {
                     regions->pts[n].r = j;
                     regions->pts[n].a = i;
                     n++;
                 } /* End if label and Border */
    fprintf(stderr,"Mask region %i, nTotal = %i %i, charge=%i\n",
          label,nTotal,n,regions->charge);

     } /* End if regionList.... */
 
     maskRegionPoints(regions->next,unwrapImage); 
     return;
}





/*****/
/*
   First pass to balance regions. Connects pairs that cancel and are mutually
   nearest neighbors
*/
    static void balanceRegionsPass1(regionList *regions, regionList *topRegion,unwrapImageStructure *unwrapImage )
{
    int dBorder1,dBorder2, dNearest;
    int rBorder, aBorder, rBorder1, aBorder1;
    int rTmp, aTmp, rTmp1, aTmp1;
    int r1,a1,r2,a2;
    regionList *tmp, *tmp1;
/*
    This pass through list balance regions that are each the closest region
    to each other.
*/
/* 
IRJ converted from recursive to loop 7/14/2011 
yields no major differences, though a couple small ones that
could be due to ordering of lables. This appears to work,
but if problems crop up, this would be somewhere to check.
*/  
    for(regions=regions; regions != NULL; regions=regions->next) {
        if(regions->charge != 0) {
            dBorder1 = distanceToBorder(regions,unwrapImage,&rTmp,  &aTmp,&rTmp1,&aTmp1);
            tmp = nearestRegion(regions,topRegion,&r1,&a1,&r2,&a2,&dNearest, IGNORECHARGE);
            tmp1 = nearestRegion(tmp,topRegion,&r1,&a1,&r2,&a2,&dNearest, IGNORECHARGE);
/*
   If mutually nearest, continue
*/
            if(tmp1 == regions) { 
                dBorder2 = distanceToBorder(tmp,unwrapImage,&rTmp,  &aTmp,&rTmp1,&aTmp1);
                fprintf(stderr,"l1,l2,n,c,c2,d,db %7i %7i %7i %7i %7i %7i %7i\n",
                       regions->label,tmp->label,regions->nTotal,regions->charge,tmp->charge,dNearest,max(dBorder1,dBorder2));
                if(dNearest < (max(dBorder1,dBorder2)+25)) {
                    if( abs(tmp->charge) > abs(regions->charge) ) {
                    subsumeRegions(topRegion,tmp,regions);
                    } else {
                        subsumeRegions(topRegion,regions,tmp);
                    }
                    drawBranchCut( unwrapImage, r1, a1, r2, a2 );
                }
	    }
	}
    }
}


static void subsumeRegions(regionList *topRegion,regionList *region1,
                           regionList *region2)
{
   int newTotal;
   rList *newList;
   regionList *tmp;
   int r1,r2,a1,a2;
   int i,d1,d2;
/*
   Merge precomputed distances
*/
   if(region1->distances != NULL) {
       for(tmp = topRegion; tmp != NULL; tmp=tmp->next) {
           d1 = (region1->distances[tmp->label]).d;
           if(tmp->nTotal >= 0) 
               d2 = distanceBetweenRegions(region2,tmp,&r1,&a1,&r2,&a2); 
           else d2 = LARGEINT;
           if(d2 < d1) { /* the new region is closer so update list */
               (region1->distances[tmp->label]).d=d2;
               (region1->distances[tmp->label]).r1=r1;
               (region1->distances[tmp->label]).r2=r2;
               (region1->distances[tmp->label]).a1=a1;
               (region1->distances[tmp->label]).a2=a2;
           }
       }              
   } else if(region1->distances == NULL && region2->distances !=NULL) {
       for(tmp = topRegion; tmp != NULL; tmp=tmp->next) {
           d2 = (region2->distances[tmp->label]).d;
           if(tmp->nTotal >= 0) 
               d1 = distanceBetweenRegions(region1,tmp,&r1,&a1,&r2,&a2); 
           else d1 = LARGEINT;
           if(d1 < d2) { /* the new region is closer so update list */
               (region2->distances[tmp->label]).d=d1;
               (region2->distances[tmp->label]).r1=r1;
               (region2->distances[tmp->label]).r2=r2;
               (region2->distances[tmp->label]).a1=a1;
               (region2->distances[tmp->label]).a2=a2;
           }
       }     
       region1->distances = region2->distances;
   } /* end else if*/
/*
    Merge list of points
*/
   newTotal = region1->nTotal + region2->nTotal;
   newList = (rList *)malloc((newTotal+1) * sizeof(rList));

   for(i=0; i < region1->nTotal; i++) {
       newList[i].r = region1->pts[i].r;
       newList[i].a = region1->pts[i].a;
   }
   for(i=0; i < region2->nTotal; i++) {
       newList[i+region1->nTotal].r = region2->pts[i].r;
       newList[i+region1->nTotal].a = region2->pts[i].a;
   }
   newList[newTotal].r = newTotal;
/*
   Free mem from orginal lists
*/
   free(region2->pts);
   free(region1->pts);
/*
   Update regions
*/
   region1->pts = newList;
   region1->nTotal = newTotal;
   region1->charge += region2->charge;
   region2->pts = NULL;
   region2->nTotal = -1;
   region2->charge = 0;
   return;
}


/*
   Second pass to balance regions. Connects pairs of regions or connects a    
   single region to the border 
*/
    static void balanceRegionsPass2( regionList *regions,regionList *topRegion, 
                                unwrapImageStructure *unwrapImage,
                                int netCharge )
{
    int dBorder = LARGEINT, d, dNearest;
    int rBorder, aBorder, rBorder1, aBorder1;
    int rTmp, aTmp, rTmp1, aTmp1;
    int r1,a1,r2,a2;
    int charge;
    int tmpCharge;
    regionList *tmp;
/* 
IRJ converted from recursive to loop 7/14/2011 
*/  
for(regions=regions; regions != NULL; regions=regions->next) {
/*
   Loop until region is balanced.
*/
    dBorder = LARGEINT;
    while( regions->charge != 0) {
/*
   Distance to border
*/
        d = distanceToBorder(regions,unwrapImage,&rTmp, &aTmp,&rTmp1,&aTmp1);
        dBorder = min(d,dBorder);
        if(d == dBorder) { 
            rBorder = rTmp; aBorder = aTmp; 
            rBorder1 = rTmp1; aBorder1 = aTmp1; 
        }
/*
    Distance to nearest region
*/
        tmp = nearestRegion(regions,topRegion,&r1,&a1,&r2,&a2,&dNearest,IGNORECHARGE);
        if(tmp == NULL) fprintf(stderr,"Null tmp in balanceRegions pass 2\n");
        fprintf(stderr,"2: l1,l2,n,c,c2,d,db %7i %7i %7i %7i %7i %7i %7i\n",
                   regions->label,tmp->label,regions->nTotal,
                   regions->charge,tmp->charge,dNearest,dBorder);
        if(dNearest < (dBorder + 10)) { 
/*
            if( abs(tmp->charge) > abs(regions->charge) ) {
                subsumeRegions(topRegion,tmp,regions);
            } else {
                subsumeRegions(topRegion,regions,tmp);
            }
*/
            subsumeRegions(topRegion,regions,tmp);
            drawBranchCut( unwrapImage, r1, a1, r2, a2 );
        } else {                 /* Connect to border */
            drawBranchCut(unwrapImage, rBorder,aBorder,rBorder1,aBorder1);
            netCharge -= regions->charge;
            regions->charge = 0;
        }
    }
}
    /*
    if(regions->next == NULL) return;
    balanceRegionsPass2( regions->next, topRegion, unwrapImage, netCharge );
    */
}


/*
    Connect points for each element in region list.
*/
    static void regionsDrawLines ( regionList *regions, 
                                  unwrapImageStructure *unwrapImage )
{   
    int i;
    int *charged1;
    int dMax;
/*fprintf(stderr,"regions->nLeft  %i %i\n",regions->nLeft,regions->label);*/
    if((regions->nTotal > 1) && (regions->mask==FALSE) ) {
        dMax = unwrapImage->blockSize -1;
        dMax = dMax * dMax + unwrapImage->blockSize * unwrapImage->blockSize;
/*fprintf(stderr,"%i\n",regions->nTotal * sizeof(int));*/
        charged1 = (int *) malloc(regions->nTotal * sizeof(int));
        for( i = 0; i < regions->nTotal; i++) {
/*fprintf(stderr,"%i\n",i);*/
             charged1[i] = TRUE;
        }
        charged1[0] = FALSE;
        (regions->nLeft)--;
 
        if(regions->nLeft < 8000) {
/*fprintf(stderr,"Connectsmall\n");*/
            connectRegionPts(regions,unwrapImage, charged1,0,dMax);
/*fprintf(stderr,"FinshConnectsmall\n");*/
        } else {
/*
   Above recursive routine takes to long for very large regions, so dilate
   regions to connect.
*/

/*fprintf(stderr,"ConnectLarge\n");*/
            connectLargeRegionPts(regions,unwrapImage);
/*fprintf(stderr,"FinshConnectLarge\n");*/

        }
        free(charged1);
    }
    if(regions->next != NULL) regionsDrawLines(regions->next,unwrapImage);
}



/*
    Free regionsList memory, recursively.
*/
    static void freeRegions( regionList *regions )
{
/* 
    Free malloced memory 
*/
/*    fprintf(stderr,"%d %d %d %d\n",regions->pts,regions->distances,regions,NULL);*/
    if(regions->next != NULL) freeRegions( regions->next );
    if(regions->pts != NULL) free(regions->pts);  
/*    if(regions->distances != NULL) free(regions->distances); IRJ removed to avoid erroroneous free's 9/23/04*/
    free(regions);
    return;
}


/*
    Cancel as many residues possible in a box around a point.
*/
    static void cancelBox(unwrapImageStructure *unwrapImage,
                          int r,int a,int boxSize)
{
     int pResidue[1000][2], nResidue[1000][2];
     int i,j;     
     int nPlus=0,nMinus=0, charge, nMatch, nMax;
     int iMin,iMax,jMin,jMax;
/*
   Set border for box of boxSize
*/
     if(r - boxSize < 0) jMin = 0; else jMin = -boxSize;
     if(r+boxSize > unwrapImage->rangeSize-2) 
         jMax= unwrapImage->rangeSize-2 - r;
     else jMax=boxSize;
     if(a-boxSize < 0) iMin = 0; else iMin = -boxSize;
     if(a+boxSize > unwrapImage->azimuthSize-2) 
        iMax= unwrapImage->azimuthSize-2-a;
     else iMax=boxSize;
/*
   Loop over box finding residues in box. Make list of pos and neg residues.
*/
     for(i=iMin; i <= iMax; i++)
         for(j=jMin; j <= jMax; j++) {
             if( (unwrapImage->bCuts[a+i][r+j] & RESIDUE) &&
                 (unwrapImage->bCuts[a+i][r+j] & LOWDENSITY) &&
                 charged(unwrapImage->bCuts[a+i][r+j]) ) {    
                    if( unwrapImage->bCuts[a+i][r+j] & PRESIDUE ) {
                         pResidue[nPlus][0] = a+i;
                         pResidue[nPlus][1] = r+j;
                         nPlus++;
                         if(nPlus > 999) error("To many residues in box");
                    } else {
                         nResidue[nMinus][0] = a+i;
                         nResidue[nMinus][1] = r+j;
                         nMinus++;
                         if(nPlus > 999) error("To many residues in box");
                    }                     
                 }
         }
        charge = nPlus - nMinus;
/* 
   Only balance equal number of + and -.
*/
        if(charge > 0) {
            nMatch = nMinus;
            nMax = nPlus;
        } else {
             nMatch = nPlus;
             nMax = nMinus;
        }
        if(nMatch == 0) return;
/*
    Cancel residues.
*/
        if(charge > 0) /* Cancel residues */
            cancelResidues(unwrapImage,nMatch,nMax,nResidue,pResidue);
        else 
            cancelResidues(unwrapImage,nMatch,nMax,pResidue,nResidue);
        return;
}


/*
    Given list of + and - residues, cancel as many as possible.
*/
    static void cancelResidues(unwrapImageStructure *unwrapImage,int nMatch,
                              int nMax,
                              int aResidues[][2],int bResidues[][2])
{
    int i, j;
    int nLeft, dMin, jMin,d1;
/*
   Simple search. Start with element in  first list then find closest uncharged
   element in second list. Not perfect search. But shouldn't be to bad for a
   small number of residues as in low density cases.
*/
    nLeft = nMatch;
    for(i = 0; i < nMatch; i++) {
         dMin = 32000;
         for(j = 0; j < nMax; j++) { /* Find minimum distance */
             if(charged(unwrapImage->bCuts[bResidues[j][0]][bResidues[j][1]])) {
		 d1 = (aResidues[i][0]-bResidues[j][0]) *
		      (aResidues[i][0]-bResidues[j][0]) + 
		      (aResidues[i][1]-bResidues[j][1]) *
		      (aResidues[i][1]-bResidues[j][1]);
		 dMin = min(dMin,d1);
                 if(dMin == d1) jMin = j;
             }
        }    
        unwrapImage->bCuts[aResidues[i][0]][aResidues[i][1]] |= UNCHARGED;
        (unwrapImage->nLeft)--;

        unwrapImage->bCuts[bResidues[jMin][0]][bResidues[jMin][1]] |= UNCHARGED;
        (unwrapImage->nLeft)--;

        drawBranchCut( unwrapImage,aResidues[i][1],aResidues[i][0],
                       bResidues[jMin][1],bResidues[jMin][0] ); /* Draw line */

    }
    return;
}


/*
    Add new elements to region list as they are found.
*/
      static regionList *buildRegionList( regionList *regions, 
                              unwrapImageStructure *unwrapImage,int r,int a)
{
     regionList *tmp;
     int label;

     label = unwrapImage->labels[a][r];

     if(regions == NULL) {
         regions = (regionList *)malloc(sizeof(regionList)); /* First element */
         regions->next = NULL;
         regions->mask = FALSE;
         regions->charge = 0;
         regions->nTotal = 0;
         regions->nLeft = 0;
         regions->pts = NULL;
         regions->label = label;
         regions->distances = NULL;
     }
     if(label == regions->label) { /* If found, update */
         regions->nTotal++;
         if(unwrapImage->bCuts[a][r] & MASK) {
             regions->mask=TRUE;
         } else {
            if(unwrapImage->bCuts[a][r] & PRESIDUE) (regions->charge)++;
            else (regions->charge)--;
         }
         return regions;
     } 
     tmp = NULL;
     if(regions->next == NULL) {  /* End of list, add element */
         tmp = (regionList *)malloc(sizeof(regionList));
         tmp->pts=NULL; /* IRJ added 9/23/04 */
         tmp->distances=NULL; /* IRJ added 9/23/04 */
         tmp->next = NULL;
     }
     else if((regions->next)->label > label) { /* Next label > so add here */
         tmp = (regionList *)malloc(sizeof(regionList));
         tmp->pts=NULL; /* IRJ added 9/23/04 */
         tmp->distances=NULL; /* IRJ added 9/23/04 */
         tmp->next = regions->next;
     }
     if(tmp != NULL) {     /* Create new element */
         tmp->nTotal = 0;
         tmp->pts = NULL;
         tmp->charge = 0;
         tmp->nLeft = 0;
         tmp->mask = FALSE;
         tmp->label = unwrapImage->labels[a][r];
         tmp->distances= NULL;
         regions->next = tmp;
     }
     buildRegionList(regions->next,unwrapImage,r,a); /* Continue */
     return regions;
}


/*
    Create array for each region in list and fill with region residues.
*/
    static void addRegionPoints( regionList *regions, 
                              unwrapImageStructure *unwrapImage,int r,int a)
{
     int n;
     regionList *tmp;
     /* Run through list */
     for(tmp=regions; tmp != NULL; tmp=tmp->next ) {

       /* If found do whats needed and return */
         
         if(unwrapImage->labels[a][r] == tmp->label) {
/***HANDLE MASK REGIONS POINTS SEPARATELY*/
         if(tmp->mask == TRUE)  return;
/****/
            if(tmp->pts == NULL) {
                 tmp->pts = (rList *)malloc((tmp->nTotal+1) * 
                                             sizeof(rList));
                 tmp->pts[tmp->nTotal].r = 0;
                 tmp->nLeft = tmp->nTotal;
	    }
             n = tmp->pts[tmp->nTotal].r;
             tmp->pts[n].r = r;
             tmp->pts[n].a = a;
             tmp->pts[tmp->nTotal].r += 1;
             return;
         } /* End if */
     }  /* endfor */
    
/*
    Error eol before label found 
*/
     error("addRegionPoints: Invalid label\n");
}

/* MODIFIED TO CREATE NONRECURSIVE VERSION 1/05/07 */
static regionList *nearestRegion(regionList *region1, 
                                      regionList *region2,
                                      int *rMin1,int *aMin1,
                                     int *rMin2,int *aMin2, int *dMin,int type )
{
    regionList *regionTmp, *closestRegion;
    int r1,r2,a1,a2, d;
    int lookUp=FALSE;
    *dMin=LARGEINT;
    /*
      Loop through to find nearest region.
    */
    /* IRJ added 7/14 to fix segfault */
    closestRegion=NULL;
    /* end fix */
    for(regionTmp=region2; regionTmp != NULL; regionTmp=regionTmp->next) {
        if( ( (((region1->charge * regionTmp->charge) < 0) && (type == OPPOSITE)) ||
             (((region1->charge * regionTmp->charge) == 0) && (type == NEUTRAL)) ||
             (((region1->charge * regionTmp->charge) > 0) && (type == SAME)) ||
             (type==IGNORECHARGE) ) && (region1 != regionTmp) && 
             (region1->nTotal >= 0) && (regionTmp->nTotal >=0) ) {
             d = distanceBetweenRegions(region1,regionTmp,&r1,&a1,&r2,&a2);    
        } else d = LARGEINT;
        /* found closer value */
        if(d < *dMin) {
	  closestRegion=regionTmp;   
          *rMin1 = r1; *aMin1 = a1; *rMin2 = r2; *aMin2 = a2;
          *dMin=d;
        }
    }
    if(*dMin >= LARGEINT) return(region2);
    return(closestRegion);
}

/*
   Find nearest region to region 1. Either same sign, opposite sign, or neutral.
*/
    static regionList *nearestRegionOld(regionList *region1,       regionList *region2,
                                      int *rMin1,int *aMin1,   int *rMin2,int *aMin2, int *dMin,int type )
{   
    regionList *tmp;
    int r1,r2,a1,a2, d;
    int lookUp=FALSE;
    /*    fprintf(stderr,"a %i\n",region2);*/
    if( ( (((region1->charge * region2->charge) < 0) && (type == OPPOSITE)) ||
          (((region1->charge * region2->charge) == 0) && (type == NEUTRAL)) ||
          (((region1->charge * region2->charge) > 0) && (type == SAME)) ||
          (type==IGNORECHARGE) ) && (region1 != region2) && 
          (region1->nTotal >= 0) && (region2->nTotal >=0) ) {
          d = distanceBetweenRegions(region1,region2,&r1,&a1,&r2,&a2);   
 
    } else d = LARGEINT;

    if(region2->next == NULL) { /* End of list return */
        *rMin1 = r1; *aMin1 = a1; *rMin2 = r2; *aMin2 = a2;
        *dMin = d; 
        return region2;
    }
        fprintf(stderr,"c %i %i\n",(int)region1,(int)region2->next);
    tmp = 
        nearestRegion(region1,region2->next,rMin1,aMin1,rMin2,aMin2,dMin,type);
    /*    fprintf(stderr,"d\n");*/
    if( d < *dMin ) { /* Is this element closer than those further down list */
/*
    Yes, so return this 1 
*/
        *rMin1 = r1; *aMin1 = a1; *rMin2 = r2; *aMin2 = a2;        
        *dMin = d; 
	    fprintf(stderr,"e\n");
        return region2;
    }
        fprintf(stderr,"f\n");
    return tmp; /* No, return previous mininum */
}


/*
    Compute shortest distance between two regions.
    Nonrecursive.
*/
    static int distanceBetweenRegions(regionList *regions1, 
                                      regionList *regions2, 
                                      int *rMin1,int *aMin1,
                                      int *rMin2,int *aMin2 )
{
    int dMin = LARGEINT, d;
    int i,j;
    int r1,r2,a1,a2;
    int np1,np2;
    int iInc,jInc;
    rList *tmp1, *tmp2;

    np1 = regions1->nTotal;
    np2 = regions2->nTotal;
/*
   Use distance look up table if it exists.
*/
    if(regions1->distances !=NULL) {
         if((regions1->distances[regions2->label]).d >= 0) {
             dMin = (regions1->distances[regions2->label]).d;
             *rMin1 = (regions1->distances[regions2->label]).r1;
             *rMin2 = (regions1->distances[regions2->label]).r2;
             *aMin1 = (regions1->distances[regions2->label]).a1;
             *aMin2 = (regions1->distances[regions2->label]).a2;
             return dMin;
          }
    } else if(regions2->distances != NULL) {
        if((regions2->distances[regions1->label]).d >= 0) {
            dMin = (regions2->distances[regions1->label]).d;
            *rMin1 = (regions2->distances[regions1->label]).r1;
            *rMin2 = (regions2->distances[regions1->label]).r2;
            *aMin1 = (regions2->distances[regions1->label]).a1;
            *aMin2 = (regions2->distances[regions1->label]).a2;
            return dMin;
        }
    }
    if(np1 < 0 || np2 < 0) error("distanceBetweenRegions: invalid region");

    if(np1 < 20) iInc = 1;
    else if(np1 < 200) iInc=2;
    else if(np1 < 1000) iInc=8;
    else if(np1 < 5000) iInc=20;
    else if(np1 < 10000) iInc=40; 
    else if(regions1->mask == TRUE) iInc=25;
    else iInc = 80;
    
    if(np2 < 20) jInc = 1; 
    else if(np2 < 200) jInc=2;
    else if(np2 < 1000) jInc=8;
    else if(np2 < 5000) jInc = 20;
    else if(np2 < 10000) jInc=40; 
    else if(regions2->mask == TRUE) jInc=25;
    else jInc=80;
    
    tmp1 = regions1->pts;
    tmp2 = regions2->pts;

    for(i = 0; i < np1; i+=iInc) {  /* Distance of 1 to 2's */
        r1 = tmp1[i].r;   a1 = tmp1[i].a;
        for(j = 0; j < np2; j+=jInc) {
            r2= tmp2[j].r; a2 = tmp2[j].a;
            d = (r1-r2)*(r1-r2) + (a1-a2)*(a1-a2);
            if( (dMin = min(d,dMin)) == d) {
               *rMin2 = r2; *aMin2 = a2; *rMin1 = r1; *aMin1 = a1;
            }
        }
    }
    dMin = (int)sqrt((double)dMin);
    return dMin;
}



/*
    Compute shortest distance of a region to the border and points to connect
    it.
*/
    static int distanceToBorder(regionList *regions,
                                unwrapImageStructure *unwrapImage,
                                int *rMin, int *aMin,int *rB,int *aB )
{
    int i;
    int dMin = LARGEINT, dLeft, dRight, dUp, dDown, d;
    int left,right,up,down;
    int r,a;
    int np;
    rList *tmp;

    np = regions->nTotal;
    if(np < 0) error("distanceToBorder: Invalid region %i %i %i",
                     regions->label,regions->charge,regions->nTotal);
    left = unwrapImage->left; right = unwrapImage->right;
    up = unwrapImage->up; down = unwrapImage->down;
    tmp = regions->pts;
    for(i = 0; i < np; i++ ) { /* Distance from + to border */
        r = tmp[i].r; a = tmp[i].a;
        dLeft = abs(r - left); dRight = abs(r - right);
        dUp = abs(a - up);     dDown = abs(a - down);
        d = min( min(dLeft,dRight),min(dUp,dDown) );
        if( (dMin = min(d,dMin)) == d ) { 
            *rMin = r; *aMin = a; 
             if(dLeft == dMin) { *rB = left; *aB = a; }
             else if(dRight == dMin) { *rB = right; *aB = a; }
             else if(dUp == dMin) { *rB = r; *aB = up; }
             else { *rB = r; *aB = down; }
        }  
    }
    return dMin;
}


/*
   Recursively, connect all points in a regions with branch cuts.
*/
    static void connectRegionPts( regionList *regions, 
                                  unwrapImageStructure *unwrapImage,
                                  int *charged1,int pt,int dMax)
{
    int i;
    int r,a;
    int r1,a1,r2,a2;
    int d;
    if(regions->nTotal == 0) return;
    if(regions->mask == TRUE) return;
    r = regions->pts[pt].r;
    a = regions->pts[pt].a;
    for(i = 0; i < regions->nTotal; i++) {
        d = (r-regions->pts[i].r) * (r-regions->pts[i].r) +
            (a-regions->pts[i].a) * (a-regions->pts[i].a);
        if(d <= dMax && charged1[i] == TRUE) {
             r1 = r; a1 = a; r2 = regions->pts[i].r; a2 = regions->pts[i].a;
             if(r1 > r2) r2++; else if(r2 < r1) r1++;
             if(a1 > a2) a2++; else if(a2 < a1) a1++; 
             drawBranchCut(unwrapImage,r1,a1,r2,a2);
             (regions->nLeft)--;
             charged1[i] = FALSE;
             if(regions->nLeft <= 0) return;
             connectRegionPts(regions,unwrapImage,charged1,i,dMax);
        }
    }    
}
/*
  Nonrecursively connect points in large region via dilation
*/
    static void connectLargeRegionPts( regionList *regions, 
                                       unwrapImageStructure *unwrapImage)
{

     int i,j,k,a,r;     
     int count=0;
     int iMin,iMax,jMin,jMax;
     int boxLeft, boxRight;
     int blockSize;
     
     if(regions->mask == TRUE) return;
     blockSize = unwrapImage->blockSize;
     boxRight = (int)(blockSize/2);
     boxLeft = boxRight;
     if( (blockSize % 2) == 0 ) boxLeft--;

     for(k = 0; k < regions->nTotal; k++) {
         r = regions->pts[k].r;
         a = regions->pts[k].a;
         if(r - boxLeft < 0) jMin = 0; else jMin = -boxLeft;
         if(r+boxRight > unwrapImage->rangeSize-2) 
	     jMax= unwrapImage->rangeSize-2 - r;
	 else jMax=boxRight;
         if(a-boxLeft < 0) iMin = 0; else iMin = -boxLeft;
         if(a+boxRight > unwrapImage->azimuthSize-2) 
             iMax= unwrapImage->azimuthSize-2-a;
         else iMax=boxRight;
         for(i=iMin; i <= iMax; i++)
            for(j=jMin;j <= jMax;j++) unwrapImage->bCuts[a+i][r+j] |= BRANCHCUT;
    }
    return;
}


/*
    Recursively work through list sortintg points in each region in order of 
    ascending azimuth and when azimuth same ascending range.
*/
    static void sortRegions( regionList *regions )
{
     if(regions->nTotal > 1) sortRegionPts(regions->nTotal,regions->pts);
     if(regions->next == NULL) return;
     sortRegions(regions->next);
}


/*
   Sort points within a region.
*/
    static void sortRegionPts(int n, rList *ra)
{
    int l,j,ir,i;
    rList rra;

    ra--; /* This takes of numerical recipes indexing */
   
     l=(n >> 1)+1;
     ir=n;
     for (;;) {
         if (l > 1)
            rra=ra[--l];
	 else {
             rra=ra[ir];
             ra[ir]=ra[1];
             if (--ir == 1) {
                ra[1]=rra;
		return;
             }
         }
         i=l;
         j=l << 1;
         while (j <= ir) {
            if (j < ir && (comparePts(ra[j],ra[j+1])== TRUE) ) ++j; /* *** */
                if(comparePts(rra ,ra[j]) == TRUE) {                /* *** */
                    ra[i]=ra[j];
                    j += (i=j);
		} else j=ir+1;
            }
            ra[i]=rra;
         }
}



/*
   Compare to points for sortRegionPts.
*/
    static int comparePts(rList r1,rList r2)
{
    if( r1.a < r2.a) return TRUE;
    else if(r1.a == r2.a) if(r1.r < r2.r) return TRUE;
    return FALSE;
}



/*
    Draw a branch cut between specified points.
*/
     static void drawBranchCut( unwrapImageStructure *unwrapImage,
                               int startX, int startY, int endX, int endY )
                               
{
    int t, distance;
    int xerr=0, yerr=0, deltaX,deltaY;
    int incX, incY;
/*
   Compute the distances in both directions 
*/
     deltaX = endX - startX;
     deltaY = endY - startY;
/*
     Compute the direction  of the increment, an increment of 0 means either
     a vertical or horizontal line.
*/
    if( deltaX > 0) incX = 1;
    else if(deltaX == 0) incX = 0;
    else incX = -1;
    
    if( deltaY > 0) incY = 1;
    else if(deltaY == 0) incY = 0;
    else incY = -1;
/*
    Determine which distance is greater
*/    
     deltaX = abs(deltaX);
     deltaY = abs(deltaY);
     if(deltaX > deltaY) distance = deltaX;
     else distance = deltaY;
/*
    Draw the line
*/
     for(t = 0; t <= distance+1; t++) {
         unwrapImage->bCuts[startY][startX] |= BRANCHCUT;
/*
unwrapImage->labels[startY][startX] = 10;
*/
         xerr += deltaX;
         yerr += deltaY;
         if( xerr > distance ) {
             xerr -= distance;
             startX += incX;
         }
         if( yerr > distance ) {
             yerr -= distance;
             startY += incY;
         }
     }
     return;
}
