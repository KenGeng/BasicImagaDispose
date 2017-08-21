//
//  main.cpp
//  图像信息处理——双边滤波
//
//  Created by apple on 2017/5/10.
//  Copyright © 2017年 apple. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "DisposeImage.hpp"

int main(){
    
    
    FILE *fp;
    char BMPfilename[50];//origin
    char Greyfilename[50];//grey
    char Resultfilename1[50];//bilateral filter
    
    int picwidth,picheight;
    long int biSizeImage;// amount of bytes of the data
    int *mingrey=(int*)malloc(sizeof(int));
    int *maxgrey=(int*)malloc(sizeof(int));
    int *mingrey1=(int*)malloc(sizeof(int));
    int *maxgrey1=(int*)malloc(sizeof(int));
    printf("Please input the BMP file path (the image must be a 24-bit BMP image):\n");
    scanf("%s",BMPfilename);
    
    printf("Please input the save path of the generated grayscale BMP file:\n");
    scanf("%s",Greyfilename);
    printf("Please input the path to the BMP file generated after the mean filter operation:\n");
    scanf("%s",Resultfilename1);
    printf("Pleasr input the save path of the BMP file generated after histogram equalization:\n");
    scanf("%s",Resultfilename2);
    
    
    if ((fp = fopen(BMPfilename, "rb")) == NULL){
        printf("\nthe file does not exist！\n");
        return 1;
    }
    else{
        
        FILE *outfp0=fopen(Greyfilename,"w+");//grey
        FILE *outfp1=fopen(Resultfilename1,"w+");//bilateral_filter
        
        
        short ind;
        char bm[2];
        char buf1[54];
        
        //get the infomation of the file header and the image header to judge whether it is a 24 bmp image
        fseek(fp,0L,0);
        fread(bm,2,1,fp);
        
        //get the infomation of the file header and the image header and write them into the result image
        fseek(fp,0L,0);
        fread(buf1, 54, 1, fp);
        fwrite(buf1, 1, 54, outfp0);
        fwrite(buf1, 1, 54, outfp1);
        fwrite(buf1, 1, 54, outfp2);
        
        //if the bmp is not a 24-bit,print the error information
        fseek(fp,28L,0);
        fread(&ind,2,1,fp);
        if(bm[0]!='B'||bm[1]!='M'||ind!=24){
            printf("\nnot a 24-bit BMP picture!\n");
            return 1;
        }
        
        fseek(fp, 18L, 0);
        fread(&picwidth, 4, 1, fp);
        fread(&picheight, 4, 1, fp);
        fseek(fp,34L,0);
        
        //biSizeImage = picheight * picwidth * 3 ;
        fread(&biSizeImage,4,1,fp);//size of the image data
        
        fseek(fp, 54L, 0);
        printf("size of the image:%d*%d %ld\n",picwidth,picheight,biSizeImage);
        
        *mingrey = 255;
        *maxgrey = 0;
        *mingrey1 = 255;
        *maxgrey1 = 0;
        grey(fp, picheight, picwidth, mingrey, maxgrey, outfp0);
        
        int window_r = 10;
        double space_singa=2;
        double range_singa=10;
        fseek(fp,54L,0);
        //bilateral_filter
        bilateral_filter(fp, picheight, picwidth, window_r, space_singa, range_singa, outfp1);
        
        
        fclose(fp);
        fclose(outfp0);
        fclose(outfp1);
        
    }
    
    
    printf("Conversion is complete!\n");
    return 0;
}

