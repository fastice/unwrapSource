#include "stdio.h"
#include "string.h"
#include "unwrapInclude.h"
#include <math.h>
#include <stdlib.h>

static void initUnwrap(unwrapImageStructure *unwrapImage);
static void phaseDifference(unwrapImageStructure *unwrapImage);
static void unwrapPass1(unwrapImageStructure *unwrapImage, int32_t *maxLabel);
static void unwrapPass2(unwrapImageStructure *unwrapImage, int32_t *maxLabel, int32_t **labelCount);
static fpNeighbors firstPassNeighbors(unwrapImageStructure *unwrapImage,
                                      int32_t i, int32_t j);
static void secondPassNeighbors(unwrapImageStructure *unwrapImage,
                                int32_t i, int32_t j, spNeighbors *spN);
static void forwardPropagate(unwrapImageStructure *unwrapImage,
                             fpNeighbors fpN, int32_t row, int32_t column, int32_t *newLabel);
static void addFPTable(unwrapImageStructure *unwrapImage, fpNeighbors fpN,
                       int32_t row, int32_t column, eqTableElement *eqTable, int32_t *nEQ);
static void addSPTable(unwrapImageStructure *unwrapImage, spNeighbors spN,
                       int32_t row, int32_t column, eqTableElement *eqTable, int32_t *nEQ);
static void resolveEQ(eqTableElement *eqTable, int32_t nEQ,
                      eqTableElement *labelList, int32_t *nLabels);
static void dfs(int32_t node, eqTableElement *eqTable, int32_t nEQ, float distance,
                eqTableElement *labelList, int32_t nLabels,
                nodeListElement *nodesFound, int32_t *nFound);
static void closeNode(int32_t node, eqTableElement *labelList, int32_t nLabels);
static int32_t nodeClosed(int32_t node, eqTableElement *labelList, int32_t nLabels);
static int32_t equivalentLabel(int32_t label, eqTableElement *labelList, int32_t nLabels);
static void finishUnwrap(unwrapImageStructure *unwrapImage,
                         int32_t maxLabel, int32_t *labelCount);

void unwrapPhase(unwrapImageStructure *unwrapImage)
{
    int32_t maxLabel, i;
    int32_t *labelCount;

    phaseDifference(unwrapImage);                     /* Compute phasedifferences */
    initUnwrap(unwrapImage);                          /* Init variables for unwrapping */
    unwrapPass1(unwrapImage, &maxLabel);              /* pass 1, top-left to bottom right*/
    unwrapPass2(unwrapImage, &maxLabel, &labelCount); /*pass 2,bottom-right to
                                                       top-left */
    finishUnwrap(unwrapImage, maxLabel, labelCount);
    /* BEGIN IRJ added 10/31/00 for resolving ambiguities*/
    unwrapImage->labelCount = labelCount;
    /* END IRJ added 10/31/00 */
    for (i = 0; i < unwrapImage->azimuthSize; i++)
    {
        free(unwrapImage->dpdr[i]);
        free(unwrapImage->dpda[i]);
    }
    return;
}

static void finishUnwrap(unwrapImageStructure *unwrapImage, int32_t maxLabel,
                         int32_t *labelCount)
{
    int32_t i, j;

    if (unwrapImage->thresh < 0)
        unwrapImage->thresh = labelCount[maxLabel];
    fprintf(stderr, "Thresh %i\n", unwrapImage->thresh);

    for (i = 0; i < unwrapImage->azimuthSize; i++)
    {
        for (j = 0; j < unwrapImage->rangeSize; j++)
        {
            if (unwrapImage->labels[i][j] >= 0 &&
                unwrapImage->labels[i][j] < LARGEINT)
            {
                if (labelCount[unwrapImage->labels[i][j]] < unwrapImage->thresh)
                    unwrapImage->phase[i][j] = (float)-LARGEINT;
            }
            else
                unwrapImage->phase[i][j] = (float)-LARGEINT;
        }
    }
}

static void initUnwrap(unwrapImageStructure *unwrapImage)
{
    int32_t i, j;
    for (i = 0; i < unwrapImage->azimuthSize; i++)
    {
        for (j = 0; j < unwrapImage->rangeSize; j++)
        { /* Compute differences */
            /*            unwrapImage->phase[i][j] = 0.0;  IRJ 11/1/00 */
            unwrapImage->labels[i][j] = LARGEINT;
        }
    }
    return;
}

static void phaseDifference(unwrapImageStructure *unwrapImage)
{
    int32_t i, j;
    float dpdr, dpda;
    float twoPi;
    twoPi = 2 * PI;
    unwrapImage->dpdr = (float **)
        malloc(unwrapImage->azimuthSize * sizeof(float *));
    unwrapImage->dpda = (float **)
        malloc(unwrapImage->azimuthSize * sizeof(float *));

    for (i = 0; i < unwrapImage->azimuthSize; i++)
    {
        unwrapImage->dpdr[i] = (float *)
            malloc(unwrapImage->rangeSize * sizeof(float));
        unwrapImage->dpda[i] = (float *)
            malloc(unwrapImage->rangeSize * sizeof(float));
    }

    for (i = 0; i < unwrapImage->azimuthSize - 1; i++)
    {
        for (j = 0; j < unwrapImage->rangeSize - 1; j++)
        { /* Compute differences */
            dpdr = unwrapImage->phase[i][j + 1] - unwrapImage->phase[i][j];
            if (dpdr > PI)
                dpdr -= twoPi;
            else if (dpdr <= -PI)
                dpdr += twoPi;
            unwrapImage->dpdr[i][j + 1] = dpdr;
            dpda = unwrapImage->phase[i + 1][j] - unwrapImage->phase[i][j];
            if (dpda > PI)
                dpda -= twoPi;
            else if (dpda <= -PI)
                dpda += twoPi;
            unwrapImage->dpda[i + 1][j] = dpda;
        }
    }
    /*
    for( i=0; i < unwrapImage->azimuthSize-1; i++) unwrapImage->dpdr[i][0];
    for( j=0; j < unwrapImage->rangeSize-1; j++) unwrapImage->dpda[0][j];
    */
    return;
}

static void unwrapPass1(unwrapImageStructure *unwrapImage, int32_t *maxLabel)
{
    fpNeighbors fpN;
    int32_t newLabel = 1;
    int32_t i, j;
    eqTableElement *eqList, *eqTable;
    int32_t nLabels, nEQ;
    int32_t eqIndex;
    /*
       malloc mem for equiv table.
    */
    eqTable = (eqTableElement *)
        malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
    eqList = (eqTableElement *)
        malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
    /*
       First pass
    */
    for (i = 0; i < unwrapImage->azimuthSize; i++)
    {
        nEQ = 0;
        for (j = 0; j < unwrapImage->rangeSize; j++)
        { /* Label row */
            if (notBranch(unwrapImage->bCuts[i][j]) &&
                validPixel(unwrapImage->bCuts[i][j]))
            { /* If not branch cut */
                /*
                    Get neighbors
                */
                fpN = firstPassNeighbors(unwrapImage, i, j);
                /*
                    Propagate label and phase
                */
                forwardPropagate(unwrapImage, fpN, i, j, &newLabel);
                /*
                    Update equivalence table
                */
                addFPTable(unwrapImage, fpN, i, j, eqTable, &nEQ);
            }
        }
        /*
           Pixels in row have now all been labeled. Now get equivalence classes for
           labels and relabel to remove equivalent labels.
        */
        resolveEQ(eqTable, nEQ, eqList, &nLabels);
        /*
           Relabel
        */
        for (j = 0; j < unwrapImage->rangeSize - 1; j++)
        {
            eqIndex = equivalentLabel(unwrapImage->labels[i][j], eqList, nLabels);
            if (eqIndex >= 0)
            {
                unwrapImage->labels[i][j] = eqList[eqIndex].b;
                unwrapImage->phase[i][j] += eqList[eqIndex].dp;
            }
        }
    }
    *maxLabel = newLabel;
    /*
        Free used memory and return
    */
    free(eqTable);
    free(eqList);
    return;
}

static void unwrapPass2(unwrapImageStructure *unwrapImage, int32_t *maxLabel,
                        int32_t **labelCount1)
{
    spNeighbors spN;
    int32_t i, j;
    eqTableElement *eqList;
    int32_t nLabels;
    int32_t *labelCount, nCount, maxCount;
    eqTableElement *eqTable;
    int32_t nEQ;
    int32_t eqIndex;
    /*
       malloc mem for equiv table.
    */
    eqTable = (eqTableElement *)
        malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
    eqList = (eqTableElement *)
        malloc((unwrapImage->rangeSize) * sizeof(eqTableElement));
    nCount = *maxLabel + 1;
    labelCount = (int32_t *)malloc(nCount * sizeof(int));
    unwrapImage->nLabels = nCount; /* IRJ added 11/1/00 for resolveAmb. */
    *labelCount1 = labelCount;
    for (i = 0; i < nCount; i++)
        labelCount[i] = 0;
    /*
       Do second pass
    */
    for (i = unwrapImage->azimuthSize - 1; i >= 0; i--)
    {
        nEQ = 0;
        for (j = unwrapImage->rangeSize - 1; j >= 0; j--)
        { /* Label row */
            if (notBranch(unwrapImage->bCuts[i][j]) &&
                validPixel(unwrapImage->bCuts[i][j]))
            { /* If not branch cut */
                /*
                    Get neighbors
                */
                secondPassNeighbors(unwrapImage, i, j, &spN);
                /*
                    Update equivalence table
                */
                addSPTable(unwrapImage, spN, i, j, eqTable, &nEQ);
            }
        }
        /*
           Pixels in row have now all been labeled. Now get equivalence classes for
           labels and relabel to remove equivalent labels.
        */

        resolveEQ(eqTable, nEQ, eqList, &nLabels);
        /*
           Relabel
        */
        for (j = 0; j < unwrapImage->rangeSize; j++)
        {
            eqIndex = equivalentLabel(unwrapImage->labels[i][j], eqList, nLabels);
            if (eqIndex >= 0)
            {
                unwrapImage->labels[i][j] = eqList[eqIndex].b;
                unwrapImage->phase[i][j] += eqList[eqIndex].dp;
            }
            if (unwrapImage->labels[i][j] < LARGEINT)
                (labelCount[unwrapImage->labels[i][j]])++;
        }
    }
    /*
        Free used memory and return
    */
    maxCount = 0;
    *maxLabel = -1;
    if (unwrapImage->fpReport != NULL)
    {
        fprintf(unwrapImage->fpReport, "; Components > 40,000 pixels\n");
        fprintf(unwrapImage->fpReport, "; Component     size    relative size\n");
    }
    for (i = 0; i < nCount; i++)
    {
        if (labelCount[i] > 300)
            fprintf(stderr, "%i %i\n", i, labelCount[i]);
        if (unwrapImage->fpReport != NULL)
        {
            if (labelCount[i] > 40000)
            {
                fprintf(unwrapImage->fpReport, "%10i %10i %6.1f %% \n",
                        i, labelCount[i],
                        (float)labelCount[i] / (float)(unwrapImage->rangeSize * unwrapImage->azimuthSize) * 100);
            }
        }
        if (labelCount[i] > maxCount)
        {
            maxCount = labelCount[i];
            *maxLabel = i;
        }
    }
    fprintf(stderr, "%i %i\n", *maxLabel, maxCount);
    free(eqTable);
    free(eqList);
    return;
}

static fpNeighbors firstPassNeighbors(unwrapImageStructure *unwrapImage,
                                      int32_t i, int32_t j)
{
    /*
        Return structure with upper and left labels.
    */
    fpNeighbors fpN;
    fpN.left = LARGEINT;
    fpN.upper = LARGEINT;

    if (i != 0 && notBranch(unwrapImage->bCuts[i - 1][j]) && validPixel(unwrapImage->bCuts[i - 1][j]))
        fpN.upper = unwrapImage->labels[i - 1][j];
    if (j != 0 && notBranch(unwrapImage->bCuts[i][j - 1]) && validPixel(unwrapImage->bCuts[i][j - 1]))
        fpN.left = unwrapImage->labels[i][j - 1];

    return fpN;
}

static void secondPassNeighbors(unwrapImageStructure *unwrapImage,
                                int32_t i, int32_t j, spNeighbors *spN)
{
    /*
        Return structure with lower and right labels.
    */
    spN->right = LARGEINT;
    spN->lower = LARGEINT;

    if (i < unwrapImage->azimuthSize - 1)
    {
        if (notBranch(unwrapImage->bCuts[i + 1][j]) && validPixel(unwrapImage->bCuts[i + 1][j]))
            spN->lower = unwrapImage->labels[i + 1][j];
    }

    if (j < unwrapImage->rangeSize - 1)
    {
        if (notBranch(unwrapImage->bCuts[i][j + 1]) && validPixel(unwrapImage->bCuts[i][j + 1]))
            spN->right = unwrapImage->labels[i][j + 1];
    }
}

static void forwardPropagate(unwrapImageStructure *unwrapImage,
                             fpNeighbors fpN, int32_t row, int32_t column, int32_t *newLabel)
{
    if (noFPNeighbors(fpN) == TRUE)
    { /* New label */
        unwrapImage->labels[row][column] = (*newLabel)++;
        /*	unwrapImage->phase[row][column] = 0.0; IRJ 11/1/00 */
        return;
    }
    if (fpN.left <= fpN.upper)
    {
        unwrapImage->labels[row][column] = fpN.left;
        unwrapImage->phase[row][column] =
            unwrapImage->phase[row][column - 1] +
            unwrapImage->dpdr[row][column];
    }
    else
    {
        unwrapImage->labels[row][column] = fpN.upper;
        unwrapImage->phase[row][column] =
            unwrapImage->phase[row - 1][column] +
            unwrapImage->dpda[row][column];
    }
}

static void addFPTable(unwrapImageStructure *unwrapImage, fpNeighbors fpN,
                       int32_t row, int32_t column, eqTableElement *eqTable, int32_t *nEQ)
/*
   If new pair, add to table.
*/
{
    int32_t i, new;
    int32_t a, b;
    int32_t currentLabel;
    float phase, phaseOffset;

    eqTable[*nEQ].a = -1;
    eqTable[*nEQ].b = -1;
    currentLabel = unwrapImage->labels[row][column];
    new = FALSE;
    /*
       Determine if upper or left is possible equivalency.
    */
    if ((fpN.left < LARGEINT) && (currentLabel < LARGEINT) &&
        (currentLabel != fpN.left))
    {
        phase = unwrapImage->phase[row][column] -
                unwrapImage->dpdr[row][column];
        phaseOffset = phase - unwrapImage->phase[row][column - 1];
        a = currentLabel;
        b = fpN.left;
        new = TRUE;
    }
    else if ((fpN.upper < LARGEINT) && (currentLabel < LARGEINT) &&
             (currentLabel != fpN.upper))
    {
        phase = unwrapImage->phase[row][column] -
                unwrapImage->dpda[row][column];
        phaseOffset = phase - unwrapImage->phase[row - 1][column];
        a = currentLabel;
        b = fpN.upper;
        new = TRUE;
    }
    /*
       Determine if equivalency is already in the table.
    */
    i = 0;
    while ((new == TRUE) && (i < *nEQ))
    { /* Check if already in table */
        if ((a == eqTable[i].a) && (b == eqTable[i].b))
            new = FALSE;
        i++;
    }
    /*
       If it is new, add to the table.
    */
    if (new == TRUE)
    {
        eqTable[*nEQ].a = a;
        eqTable[*nEQ].b = b;
        eqTable[*nEQ].dp = phaseOffset;
        (*nEQ)++;
    }
    return;
}

static void addSPTable(unwrapImageStructure *unwrapImage, spNeighbors spN,
                       int32_t row, int32_t column, eqTableElement *eqTable, int32_t *nEQ)
{
    int32_t i, new;
    int32_t a, b;
    int32_t currentLabel;
    float phase, phaseOffset;

    eqTable[*nEQ].a = -1;
    eqTable[*nEQ].b = -1;
    currentLabel = unwrapImage->labels[row][column];
    new = FALSE;
    /*
       Determine if upper or left is possible equivalency.
    */
    if ((spN.right < LARGEINT) && (currentLabel < LARGEINT) &&
        (currentLabel != spN.right))
    {
        phase = unwrapImage->phase[row][column + 1] -
                unwrapImage->dpdr[row][column + 1];
        phaseOffset = phase - unwrapImage->phase[row][column];
        b = currentLabel;
        a = spN.right;
        new = TRUE;
    }
    else if ((spN.lower < LARGEINT) && (currentLabel < LARGEINT) &&
             (currentLabel != spN.lower))
    {
        phase = unwrapImage->phase[row + 1][column] -
                unwrapImage->dpda[row + 1][column];
        phaseOffset = phase - unwrapImage->phase[row][column];
        b = currentLabel;
        a = spN.lower;
        new = TRUE;
    }
    /*
       Determine if equivalency is already in the table.
    */
    i = 0;
    while ((new == TRUE) && (i < *nEQ))
    { /* Check if already in table */
        if ((a == eqTable[i].a) && (b == eqTable[i].b))
            new = FALSE;
        i++;
    }
    /*
       If it is new, add to the table.
    */
    if (new == TRUE)
    {
        eqTable[*nEQ].a = a;
        eqTable[*nEQ].b = b;
        eqTable[*nEQ].dp = phaseOffset;
        (*nEQ)++;
    }
    return;
}

static void resolveEQ(eqTableElement *eqTable, int32_t nEQ,
                      eqTableElement *labelList, int32_t *nLabels)
{
    int32_t i, j, k;
    int32_t node, nFound;
    nodeListElement labelsFound[10000];
    int32_t a, b, minLabel;

    for (i = 0; i < 2 * nEQ; i++)
        labelList[i].a = -1;
    *nLabels = 0;
    for (i = 0; i < nEQ; i++)
    { /* Build list of labels */
        a = eqTable[i].a;
        b = eqTable[i].b;
        for (j = 0; j < *nLabels; j++)
        {
            if (labelList[j].a == a)
                a = -1;
            if (labelList[j].a == b)
                b = -1;
        }
        if (a > 0)
        {
            labelList[*nLabels].a = a;
            labelList[*nLabels].b = -1;
            (*nLabels)++;
        }
        if (b > 0)
        {
            labelList[*nLabels].a = b;
            labelList[*nLabels].b = -1;
            (*nLabels)++;
        }
    }
    for (i = 0; i < *nLabels; i++)
    { /* Resolve labels */
        nFound = 0;
        minLabel = LARGEINT;
        for (k = 0; k < *nLabels; k++) /* Find min open node */
            if (labelList[k].b < 0)
                minLabel = min(minLabel, labelList[k].a);
        node = minLabel;
        if (node == LARGEINT)
            return; /* Return if all nodes closed */

        dfs(node, eqTable, nEQ, 0.0, labelList, *nLabels,
            labelsFound, &nFound); /*Search node */

        for (j = 0; j < nFound; j++)
            minLabel = min(minLabel, labelList[j].a);
        for (j = 0; j < nFound; j++) /* Mark labels with equivalent labels */
            for (k = 0; k < *nLabels; k++)
            {
                if (labelList[k].a == labelsFound[j].node)
                {
                    labelList[k].b = labelsFound[0].node;
                    labelList[k].dp = labelsFound[j].distance;
                }
            }
    }
}

static void dfs(int32_t node, eqTableElement *eqTable, int32_t nEQ, float distance,
                eqTableElement *labelList, int32_t nLabels,
                nodeListElement *nodesFound, int32_t *nFound)
{
    int32_t i;
    int32_t newNode;
    float newDistance;
    if (nodeClosed(node, labelList, nLabels) == TRUE)
        return; /* Node closed return */

    closeNode(node, labelList, nLabels); /* Close node */
    nodesFound[*nFound].node = node;     /* Add node to found list */
    nodesFound[*nFound].distance = distance;
    (*nFound)++;
    /*
       Search node
    */
    for (i = 0; i < nEQ; i++)
    { /* Search nodes */
        newNode = -1;
        if (node == eqTable[i].a)
        {
            newNode = eqTable[i].b;
            newDistance = distance + eqTable[i].dp;
        }
        else if (node == eqTable[i].b)
        {
            newNode = eqTable[i].a;
            newDistance = distance - eqTable[i].dp;
        }
        if (newNode > -1)
            dfs(newNode, eqTable, nEQ, newDistance, labelList, nLabels,
                nodesFound, nFound);
    }
}

static void closeNode(int32_t node, eqTableElement *labelList, int32_t nLabels)
{
    int32_t i;
    for (i = 0; i < nLabels; i++)
        if (node == labelList[i].a)
        {
            labelList[i].b = 0;
            return;
        }
    error("*** closeNode: Invalid node ***\n");
}

static int32_t nodeClosed(int32_t node, eqTableElement *labelList, int32_t nLabels)
{
    int32_t i;
    for (i = 0; i < nLabels; i++)
        if (node == labelList[i].a)
        {
            if (labelList[i].b > -1)
                return TRUE;
            else
                return FALSE;
        }
    error("*** nodeClosed: Invalid node ***\n");
    return 0;
}

static int32_t equivalentLabel(int32_t currentlabel, eqTableElement *labelList,
                               int32_t nLabels)
{
    /*
       Called by unwrapPass1 & 2. Find label in table and return equivalent label.
    */
    int32_t i;
    for (i = 0; i < nLabels; i++) /* Find label and return equivalent */
        if (currentlabel == labelList[i].a)
            return i;
    return -1; /* Label unique */
}
