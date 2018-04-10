#include <stdio.h>
#include <stdlib.h>
#include "X11/xpm.h"

int main(int argc, char *argv[]){
    XpmImage XPM;
    XpmInfo  XPM_info;
    FILE *fp;
    char* fin, fout;
    int i, j, o;
    int *acorr;
    
    fin = argv[1];
    fout = argv[2];
    
    printf("Start Reading file : %s\n", fin);
    
    XpmReadFileToXpmImage(fin, &XPM, &XPM_info);
    
    printf("Done Reading!\n");
    printf("****************\nXpmImage Contents\n");
    printf("width  : %10d\n", XPM.width);
    printf("height : %10d\n", XPM.height);
    printf("cpp    : %10d\n", XPM.cpp);
    printf("ncolors: %10d\n", XPM.ncolors);
    
    for (i=0;i<XPM.ncolors;i++){
        printf("%7d)%7s:%7s\n",i, XPM.colorTable[i].string, XPM.colorTable[i].c_color); 
    }
    
    printf("****************\nXpmInfo Contents\n");
    printf("valuemask      : %20d\n", XPM_info.valuemask);
    printf("hints  comment : %20s\n", XPM_info.hints_cmt);
    printf("colors comment : %20s\n", XPM_info.colors_cmt);
    printf("pixels comment : %20s\n", XPM_info.pixels_cmt);
    printf("x_hotspot      : %20d\n", XPM_info.x_hotspot);
    printf("y_hotspot      : %20d\n", XPM_info.y_hotspot);
    //printf("nextensions    : %20d\n", XPM_info.nextensions);
    
    //for (i=0;i<XPM_info.nextensions;i++){
        //printf("%7d)%7s:%7d\n",i, XPM_info.extensions[i].name, XPM_info.extensions[i].nlines);
    //}
    
    printf("Start writing file : %s\n", fout);
    acorr = (int *) calloc(XPM.height*XPM.width/2, sizeof(int));
    
    //fp = fopen(fout,'w');
    for (i=0;i<XPM.height;i++){
        for (o=0;o<(XPM.width/2);o++){
            acorr[o+i*(XPM.width/2)] = 0;
            printf ("Calculating delay %3d for set %3d\n", o,i);
            for (j=0;j<(XPM.width-o);j++){
                acorr[o+i*XPM.width/2] += XPM.data[j+i*XPM.width] * XPM.data[(j+o)+i*XPM.width];
            }
            printf("%6d ",acorr[o+i*XPM.width/2]);
        }
        printf("\n");
    }
    //fclose(fp);
    
    return 0;
}
