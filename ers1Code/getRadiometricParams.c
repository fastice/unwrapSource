#include "ers1.h"
#include "string.h"

/*
    Get radiometric params
*/

    void getRadiometricParams(char *radiometricFile,float *a1,
                               float *a2,float *a3) 
{
    FILE *fp;
    char tmp[128],*string,*line;
    int notDone,i,lineLength;  
    fp = openInputFile(radiometricFile);  

    line = tmp; 
    i = 0;
    notDone = TRUE;
    while(notDone == TRUE) {
        lineLength = fgetline(fp,line,LINEMAX); /* Read line */
        if( (string = strstr(line,"NOISE_SCALE_FACTOR")) != NULL) {
             string = strstr(string,"=");
             string++; 
             sscanf(string,"%f",a1);                 
       } else if( (string = strstr(line,"LINEAR_CONVERSION_FACTOR")) != NULL ){
             string = strstr(string,"=");  
             string++;
             sscanf(string,"%f",a2);    
       } else if( (string = strstr(line,"OFFSET_CONVERSION_FACTOR")) != NULL){
             string = strstr(string,"=");
             string++;
             sscanf(string,"%f",a3);
             notDone = FALSE;
       } else if(i > 500) 
            error("*** getRadiometricParams : Data not found ***\n");
       i++;
    }
    
}
