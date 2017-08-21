//
//  DisposeImage.hpp
//  图像信息处理_a2_二值化
//
//  Created by apple on 2017/4/5.
//  Copyright © 2017年 apple. All rights reserved.
//

#ifndef DisposeImage_hpp
#define DisposeImage_hpp

#include <stdio.h>

int FindThreshold(int mingrey,int maxgrey,int picheight,int picwidth,FILE* infp);

void readRGB(unsigned char *Ptr_b,unsigned char *Ptr_g,unsigned char *Ptr_r,FILE*infp);

void writeRGB(void *Ptr_b,void *Ptr_g,void *Ptr_r,FILE*outfp);
void rgb2yuv(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp);
void bmp_binarization(int picheight, int picwidth,int threshold,FILE* infp,FILE* outfp);

void erosion_delation(unsigned char* data, int width, int height,unsigned char* buffer_re,int model);

void delation(unsigned char* data, int width, int height,unsigned char* data2);

void opening(unsigned char* data, int width, int height,unsigned char* data2);

void closing(unsigned char* data, int width, int height,unsigned char* data2);

void grey(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp);

void logarithmic_operation(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp);

void Histogram_equalization(FILE *infp,int picheight,int picwidth,int*mingrey,int*maxgrey,FILE*outfp);

void Histogram_equalization_color(FILE *infp,int picheight,int picwidth,int*mingrey,int*maxgrey,FILE*outfp);

void logarithmic_operation_color(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp);

void translation(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,int dx,int dy);
void mirror(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp);
void rotation(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,double angle);
void scale(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,double ratio_x,double ration_y);

void shear(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,double dx,double dy);
unsigned char bilinear_interpolation(unsigned char x0, unsigned char y0,unsigned char x1,unsigned char y1,double dx, double dy);
void yuv2rgb(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp);
void yuv2rgb_value(unsigned char* y,unsigned char*u,unsigned char*v,unsigned char*r,unsigned char*g,unsigned char*b);
void mean_filter(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,int* mask);

void laplace_filter(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,int* mask);

void bilateral_filter(FILE* infp,int picheight,int picwidth,int window_r,double space_singa,double range_singa,FILE* outfp);

#endif /* DisposeImage_hpp */
