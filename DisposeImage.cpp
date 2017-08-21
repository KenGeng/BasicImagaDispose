//
//  DisposeImage.cpp
//  图像信息处理_a2_二值化
//
//  Created by apple on 2017/4/5.
//  Copyright © 2017年 apple. All rights reserved.
//

#include "DisposeImage.hpp"
#include <iostream>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#define PI 3.1415926
//to simplify the read code
void readRGB(unsigned char *Ptr_b,unsigned char *Ptr_g,unsigned char *Ptr_r,FILE*infp){
    fread(Ptr_b, 1, 1, infp);
    fread(Ptr_g, 1, 1, infp);
    fread(Ptr_r, 1, 1, infp);
    
}

//to simplify the write code
void writeRGB(void *Ptr_b,void *Ptr_g,void *Ptr_r,FILE*outfp){
    fwrite(Ptr_b, 1, 1, outfp);
    fwrite(Ptr_g, 1, 1, outfp);
    fwrite(Ptr_r, 1, 1, outfp);
}


void bilateral_filter(FILE* infp,int picheight,int picwidth,int window_r,double space_singa,double range_singa,FILE* outfp){
    
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);
    //acquire the gery image to handle
    int *mingrey=(int*)malloc(sizeof(int));
    int *maxgrey=(int*)malloc(sizeof(int));
    grey(infp, picheight, picwidth, mingrey, maxgrey, outfp);
    
    fseek(outfp, 54L, 0);
    //read the grey value into the buffer
    fread(buffer, 3*picwidth*picheight, 1, outfp);
    
    
    if (window_r<0) {
        printf("wrong input of window radius");
    }
    
    double space_weight[2*window_r+1][2*window_r+1];

    //bulid the globl space kernal for use to save time
    for (int i = -window_r; i<=window_r; i++) {
        for (int j = -window_r; j<=window_r; j++) {
            space_weight[i+window_r][j+window_r]=exp(-0.5*(i*i+j*j)/space_singa/space_singa);
        }
    }
    double ty=0;
    fseek(outfp, 54L, 0);
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            int y = buffer[3*i*picwidth+3*j];
            if (1/*i>=window_r&&j>=window_r&&i+window_r<picheight&&j+window_r<picwidth*/) {
                double w_sum=0;
                double w_cof_sum=0;
               
                for (int w_col = -window_r; w_col<=window_r; w_col++){
                    for (int w_row = -window_r; w_row<=window_r; w_row++){
                        int row,col;
                        row = w_row+j;
                        col = w_col+i;
//                        if (row<0) {
//                            row=-row;
//                        }else if(row >= picwidth)	row = 2 * picwidth - row;
//                        
//                        if (col<0) {
//                            col=-col;
//                        }else if(col >= picheight)	col = 2 * picheight - col;
                        int ty;
                        if(row<0||col<0||row >= picwidth||col >= picheight)
                         ty = 0;
                        else ty=buffer[3*(col)*picwidth+3*(row)];
                        double range_cof = exp(-0.5*(ty-y)*(ty-y)/range_singa/range_singa);
                        double total_cof = range_cof*space_weight[w_row+window_r][w_col+window_r];
                        w_cof_sum += total_cof;
                        w_sum += ty*total_cof;
                        
                    }
                }
                
                 ty = w_sum/w_cof_sum;
                 y=ty;
                writeRGB(&y,&y,&y,outfp);
            }else{
                y = buffer[3*i*picwidth+3*j];

                writeRGB(&y,&y,&y,outfp);
            }
            
            
        }
}

void laplace_filter(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,int* mask){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);
    //acquire the gery image to handle
    rgb2yuv(infp, picheight, picwidth, mingrey, maxgrey, outfp);
    
    fseek(outfp, 54L, 0);
    //read the grey value into the buffer
    fread(buffer, 3*picwidth*picheight, 1, outfp);
    unsigned char y = 0;
    unsigned char u,v;
    unsigned char y3 = 0;
    fseek(outfp, 54L, 0);
    
    //int mask[9]={1,1,1,1,-8,1,1,1,1};
    *mingrey =255;
    *maxgrey = 0;
    int ty;
    //find the largest and smallest grey value after laplace enhancement
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            if (i>=1&&j>=1&&i<=picheight-2&&j<=picwidth-2) {
                ty=0;
                for (int m =i-1; m<=i+1; m++) {
                    for (int n=j-1; n<=j+1; n++) {
                        ty+=buffer[3*m*picwidth+3*n]*mask[(m-i+1)*3+n-j+1];
                    }
                }
                if (*mingrey>ty) {
                    *mingrey = ty;
                }
                if (*maxgrey<ty) {
                    *maxgrey =ty;
                }
            }
            
        }
    
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            ty=0;
            if (i>=1&&j>=1&&i<=picheight-2&&j<=picwidth-2) {
                for (int m =i-1; m<=i+1; m++) {
                    for (int n=j-1; n<=j+1; n++) {
                        ty+=buffer[3*m*picwidth+3*n]*mask[(m-i+1)*3+n-j+1];
                    }
                }
                //ty = (buffer[3*(i-1)*picwidth+3*(j-1)]*mask[0]+buffer[3*(i-1)*picwidth+3*j]*mask[1]+buffer[3*(i-1)*picwidth+3*(j+1)]*mask[2]+buffer[3*i*picwidth+3*(j-1)]*mask[3]+buffer[3*i*picwidth+3*j]*mask[4]+buffer[3*i*picwidth+3*(j+1)]*mask[5]+buffer[3*(i+1)*picwidth+3*(j-1)]*mask[6]+buffer[3*(i+1)*picwidth+3*j]*mask[7]+buffer[3*(i+1)*picwidth+3*(j+1)]*mask[8]);
                y3 =buffer[3*i*picwidth+3*j];

                int y2=ty;
                y2 = 1.0*(ty-*mingrey)/(*maxgrey-*mingrey)*255;//rearrange
                //write
                int dy = 1.0*(y2*mask[1]*-1)/9;
                y = y3+dy;
                if ((y3+dy)>255) {
                    y= 255;
                }
                if ((y3+dy)<=0) {
                    y=0;
                }
            }
            
            else  y = buffer[3*i*picwidth+3*j];//border pixels don't change
            u = buffer[3*i*picwidth+3*j+1];
            v = buffer[3*i*picwidth+3*j+2];

            unsigned char r,g,b;
            yuv2rgb_value(&y, &u, &v, &r, &g, &b);
            writeRGB(&b,&g,&r,outfp);
            
        }
    printf("%d %d\n",*mingrey,*maxgrey);
    free(buffer);

    

}
void mean_filter(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,int* mask){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);
    //acquire the gery image to handle
    rgb2yuv(infp, picheight, picwidth, mingrey, maxgrey, outfp);
    
    fseek(outfp, 54L, 0);
    //read the grey value into the buffer
    fread(buffer, 3*picwidth*picheight, 1, outfp);
    unsigned char y = 0;
    unsigned char u,v;
    fseek(outfp, 54L, 0);
    
    //int mask[9]={1,1,1,1,1,1,1,1,1,};
    int ty=0;
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            if (i>=1&&j>=1&&i<=picheight-2&&j<=picwidth-2) {
                for (int m =i-1; m<=i+1; m++) {
                    for (int n=j-1; n<=j+1; n++) {
                        ty+=buffer[3*m*picwidth+3*n]*mask[(m-i+1)*3+n-j+1];
                    }
                }
                ty/=9;
                y=ty;
                if (ty>255) {
                    y=255;
                }
                if (ty<0) {
                    y=0;
                }
               //get the new y
            }
           
           else  y = buffer[3*i*picwidth+3*j];
            
            //write
            u = buffer[3*i*picwidth+3*j+1];
            v = buffer[3*i*picwidth+3*j+2];
            
            unsigned char r,g,b;
            yuv2rgb_value(&y, &u, &v, &r, &g, &b);
            writeRGB(&b,&g,&r,outfp);
        }
    free(buffer);

}

unsigned char bilinear_interpolation(unsigned char d0, unsigned char d1,unsigned char d2,unsigned char d3,double dx, double dy){
    int x1 = dx*d0+(1-dx)*d1;
    int x2 = dx*d2+(1-dx)*d3;
    unsigned char re = dy*x1+(1-dy)*x2;
//    unsigned char result =d0-dx*d0+dx*d1+dy*(d2+dx*-d2+dx*d3-d0+dx*d0-dx*d1);
    return re;
}

void shear(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,double dx,double dy){
    unsigned char buffer[picheight][picwidth][3];
    fseek(infp, 54L, 0);
    fread(buffer, 3*picwidth*picheight, 1, infp);
    fseek(outfp, 18L, 0);
    int new_picwidth = picwidth*1+picheight*fabs(dx)+3;
    new_picwidth = (new_picwidth)/4*4;//4字节对齐
    int new_picheight = picheight*1+picwidth*(dy)+3;
    new_picheight = (new_picheight)/4*4;
    fwrite(&new_picwidth, 4, 1, outfp);
    fwrite(&new_picheight, 4, 1, outfp);
    fseek(outfp,34L,0);//图像数据区的大小
    int biSizeImage = 3*new_picheight*new_picwidth;
    fwrite(&biSizeImage, 4, 1, outfp);
    unsigned int bfSize = biSizeImage+54;
    fseek(outfp, 2L, 0);
    fwrite(&bfSize, 4, 1, outfp);
    
    unsigned char re_buffer[new_picheight][new_picwidth][3];
    memset(re_buffer, 255, new_picheight*new_picwidth*3);
    int xoffset = picwidth/2;
    int yoffset = picheight/2;
    int new_xoffset = new_picwidth/2;
    int new_yoffset = new_picheight/2;
    int new_x,new_y;
    int x,y;
    
    for (int i =-new_yoffset; i<new_yoffset; i++) {
        for (int j = -new_xoffset; j<new_xoffset; j++) {
            
            new_x = j;
            new_y = i;
//            if (dx==0&&dy!=0) {
//                
//            }
            double tx,ty;
            tx = 1.0*(new_x-dx*new_y);
            ty = 1.0*(new_y-dy*new_x);
            x = (new_x-dx*new_y);
            y = (new_y-dy*new_x);
            //bilinear interpolation method
            double dx = fabs(tx-x);
            double dy = fabs(ty-y);
            if (x==xoffset-1||y==yoffset-1) {
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
            }
            if (x<xoffset-1&&x>=-xoffset&&y<yoffset-1&&y>=-yoffset) {
                //judeg whether the source x,y is an integer
                if ( fabs(tx - (int)tx) <= FLT_EPSILON&&fabs(ty - (int)ty) <= FLT_EPSILON){
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
                }
                else if(fabs(tx - (int)tx) <= FLT_EPSILON&&fabs(ty - (int)ty) > FLT_EPSILON){
                    
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+xoffset][0], buffer[y+1+yoffset][x+xoffset][0], buffer[y+1+yoffset][x+xoffset][0], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] =bilinear_interpolation(buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+xoffset][1], buffer[y+1+yoffset][x+xoffset][1], buffer[y+1+yoffset][x+xoffset][1], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+xoffset][2], buffer[y+1+yoffset][x+xoffset][2], buffer[y+1+yoffset][x+xoffset][2], dx, dy );
                }
                else if(fabs(tx - (int)tx) > FLT_EPSILON&&fabs(ty - (int)ty) <= FLT_EPSILON){
                    
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+xoffset+1][0], buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+1+xoffset][0], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] =bilinear_interpolation(buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+xoffset+1][1], buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+1+xoffset][1], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+xoffset+1][2], buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+1+xoffset][2], dx, dy );
                }
                else {
                    
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+xoffset+1][0], buffer[y+1+yoffset][x+xoffset][0], buffer[y+1+yoffset][x+1+xoffset][0], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] =bilinear_interpolation(buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+xoffset+1][1], buffer[y+1+yoffset][x+xoffset][1], buffer[y+1+yoffset][x+1+xoffset][1], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+xoffset+1][2], buffer[y+1+yoffset][x+xoffset][2], buffer[y+1+yoffset][x+1+xoffset][2], dx, dy );
                }
                
            }
//            x = (new_x-dx*new_y)/(1-dx*dy);
//            y = (new_y-dy*new_x)/(1-dx*dy);
//            if (x<xoffset&&x>=-xoffset&&y<yoffset&&y>=-yoffset) {
//                
//                re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
//                re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
//                re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
//            }
            
        }
    }
    
    fseek(outfp, 54L, 0);
    fwrite(re_buffer, 1, new_picheight*new_picwidth*3, outfp);
    

}
void scale(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,double ratio_x,double ration_y){
    
    unsigned char buffer[picheight][picwidth][3];
    fseek(infp, 54L, 0);
    fread(buffer, 3*picwidth*picheight, 1, infp);
    fseek(outfp, 18L, 0);
    int new_picwidth = picwidth*ratio_x+3;
    new_picwidth = (int)((int)new_picwidth)/4*4;//4字节对齐
    int new_picheight = picheight*ration_y+3;
    new_picheight = (int)((int)new_picheight)/4*4;
    fwrite(&new_picwidth, 4, 1, outfp);
    fwrite(&new_picheight, 4, 1, outfp);
    fseek(outfp,34L,0);//图像数据区的大小
    int biSizeImage = 3*new_picheight*new_picwidth;
    fwrite(&biSizeImage, 4, 1, outfp);
    unsigned int bfSize = biSizeImage+54;
    fseek(outfp, 2L, 0);
    fwrite(&bfSize, 4, 1, outfp);
    //    fread(&new_picwidth, 4, 1, outfp);
    //    fread(&new_picheight, 4, 1, outfp);
    //    std::cout<<std::endl<<"$"<<new_picwidth<<" "<<new_picheight;
    
    unsigned char re_buffer[(int)new_picheight][(int)new_picwidth][3];
    memset(re_buffer, 255, new_picheight*new_picwidth*3);
    int xoffset = picwidth/2;
    int yoffset = picheight/2;
    int new_xoffset = new_picwidth/2;
    int new_yoffset = new_picheight/2;
    int new_x,new_y;
    int x,y;

    for (int i =-new_yoffset; i<new_yoffset; i++) {
        for (int j = -new_xoffset; j<new_xoffset; j++) {
//            x2 = x*ratio_w
//            y2 = y*ratio_h
//            逆变换：
//            x = x2/ratio_w
//            y = y2/ratio_h
          
            double tx,ty;
            
            new_x = j;
            new_y = i;
            tx = 1.0*new_x/ratio_x;
            ty = 1.0*new_y/ration_y;
            x = new_x/ratio_x;
            y = new_y/ration_y;
            //bilinear interpolation method
            double dx = fabs(tx-x);
            double dy = fabs(ty-y);
            if (x==xoffset-1||y==yoffset-1) {
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
            }
            if (x<xoffset-1&&x>=-xoffset&&y<yoffset-1&&y>=-yoffset) {
                //judeg whether the source x,y is an integer
                if ( fabs(tx - (int)tx) <= FLT_EPSILON&&fabs(ty - (int)ty) <= FLT_EPSILON){
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
                }
                //if x is an integer,then just handle y
                else if(fabs(tx - (int)tx) <= FLT_EPSILON&&fabs(ty - (int)ty) > FLT_EPSILON){
                    
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+xoffset][0], buffer[y+1+yoffset][x+xoffset][0], buffer[y+1+yoffset][x+xoffset][0], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] =bilinear_interpolation(buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+xoffset][1], buffer[y+1+yoffset][x+xoffset][1], buffer[y+1+yoffset][x+xoffset][1], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+xoffset][2], buffer[y+1+yoffset][x+xoffset][2], buffer[y+1+yoffset][x+xoffset][2], dx, dy );
                }
                //if y is an integer,then just handle x
                else if(fabs(tx - (int)tx) > FLT_EPSILON&&fabs(ty - (int)ty) <= FLT_EPSILON){
                    
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+xoffset+1][0], buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+1+xoffset][0], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] =bilinear_interpolation(buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+xoffset+1][1], buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+1+xoffset][1], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+xoffset+1][2], buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+1+xoffset][2], dx, dy );
                }
                else {
                    
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+xoffset+1][0], buffer[y+1+yoffset][x+xoffset][0], buffer[y+1+yoffset][x+1+xoffset][0], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] =bilinear_interpolation(buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+xoffset+1][1], buffer[y+1+yoffset][x+xoffset][1], buffer[y+1+yoffset][x+1+xoffset][1], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+xoffset+1][2], buffer[y+1+yoffset][x+xoffset][2], buffer[y+1+yoffset][x+1+xoffset][2], dx, dy );
                }
                
            }
//Nearest neighbor based interpolation method:
//            x = new_x/ratio_x;
//            y = new_y/ration_y;
//            if (x<xoffset&&x>=-xoffset&&y<yoffset&&y>=-yoffset) {
//
//                re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
//                re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
//                re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
//            }
//            
        }
    }
    
    fseek(outfp, 54L, 0);
    fwrite(re_buffer, 1, new_picheight*new_picwidth*3, outfp);

}

void rotation(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,double angle){
    angle = angle/180*PI;//转成弧度制
    
    //unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);//dynamic array version
    unsigned char buffer[picheight][picwidth][3];
    fseek(infp, 54L, 0);
    fread(buffer, 3*picwidth*picheight, 1, infp);
    fseek(outfp, 18L, 0);
    int new_picwidth = picwidth*fabs(cos(angle)) + picheight*fabs(sin(angle))+3;//+3是为了下面的对齐
    new_picwidth = new_picwidth/4*4;//4字节对齐
    int new_picheight = picwidth*fabs(sin(angle)) + picheight*fabs(cos(angle))+3;
    new_picheight = new_picheight/4*4;
    fwrite(&new_picwidth, 4, 1, outfp);
    fwrite(&new_picheight, 4, 1, outfp);
    fseek(outfp,34L,0);//图像数据区的大小
    int biSizeImage = 3*new_picheight*new_picwidth;
    fwrite(&biSizeImage, 4, 1, outfp);
    unsigned int bfSize = biSizeImage+54;
    fseek(outfp, 2L, 0);
    fwrite(&bfSize, 4, 1, outfp);
    //    fread(&new_picwidth, 4, 1, outfp);
    //    fread(&new_picheight, 4, 1, outfp);
    //    std::cout<<std::endl<<"$"<<new_picwidth<<" "<<new_picheight;

    unsigned char re_buffer[new_picheight][new_picwidth][3];
    memset(re_buffer, 255, new_picheight*new_picwidth*3);
    int xoffset = picwidth/2;
    int yoffset = picheight/2;
    int new_xoffset = new_picwidth/2;
    int new_yoffset = new_picheight/2;
    int new_x,new_y;
    int x,y;
//    for (int i =-yoffset; i<yoffset; i++) {
//        for (int j = -xoffset; j<xoffset; j++) {
//            x2 = x cosa - ysina
//            y2 = xsina + ycosa
//            逆变换：
//            x = x2cosa+y2sina
//            y = y2cosa-x2sina
//            new_x = j*cos(angle)-i*sin(angle);
//            new_y = j*sin(angle)+i*cos(angle);
//            re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[i+yoffset][j+xoffset][0];
//            re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[i+yoffset][j+xoffset][1];
//            re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[i+yoffset][j+xoffset][2];
//        }
//    }
    for (int i =-new_yoffset; i<new_yoffset; i++) {
        for (int j = -new_xoffset; j<new_xoffset; j++) {
            //x2 = x cosa - ysina
            //y2 = xsina + ycosa
            //逆变换：
            //x = x2cosa+y2sina
            //y = y2cosa-x2sina
           
            double tx,ty;
            //bilinear interpolation
            new_x = j;
            new_y = i;
            tx = new_x*cos(angle)+new_y*sin(angle);
            ty = -new_x*sin(angle)+new_y*cos(angle);
            new_x = j;
            new_y = i;
            x = new_x*cos(angle)+new_y*sin(angle);
            y = -new_x*sin(angle)+new_y*cos(angle);
            double dx = fabs(tx-x);
            double dy = fabs(ty-y);
            if (x==xoffset-1||y==yoffset-1) {
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
                re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
            }

            if (x<xoffset-1&&x>=-xoffset&&y<yoffset-1&&y>=-yoffset) {
                //judeg whether the source x,y is an integer
                if ( fabs(tx - (int)tx)<= FLT_EPSILON&&fabs(ty - (int)ty) <= FLT_EPSILON){
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[y+yoffset][x+xoffset][0];
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[y+yoffset][x+xoffset][1];
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[y+yoffset][x+xoffset][2];
                }
                else{
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][0], buffer[y+yoffset][x+xoffset+1][0], buffer[y+1+yoffset][x+xoffset][0], buffer[y+1+yoffset][x+1+xoffset][0], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] =bilinear_interpolation(buffer[y+yoffset][x+xoffset][1], buffer[y+yoffset][x+xoffset+1][1], buffer[y+1+yoffset][x+xoffset][1], buffer[y+1+yoffset][x+1+xoffset][1], dx, dy );
                    re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = bilinear_interpolation(buffer[y+yoffset][x+xoffset][2], buffer[y+yoffset][x+xoffset+1][2], buffer[y+1+yoffset][x+xoffset][2], buffer[y+1+yoffset][x+1+xoffset][2], dx, dy );
                    
                }
                
            }
//Nearest neighbor based interpolation:
//this version is for dynamic array
//if (x<xoffset&&x>=-xoffset&&y<yoffset&&y>=-yoffset) {
//   re_buffer[new_y+new_yoffset][new_x+new_xoffset][0] = buffer[3*picwidth*(y+yoffset)+3*(x+xoffset)];
//// re_buffer[new_y+new_yoffset][new_x+new_xoffset][1] = buffer[3*picwidth* (y+yoffset)+3*(x+xoffset)+1];
//// re_buffer[new_y+new_yoffset][new_x+new_xoffset][2] = buffer[3*picwidth* (y+yoffset)+3*(x+xoffset)+2];
//            }
//            
        }
    }

    fseek(outfp, 54L, 0);
    fwrite(re_buffer, 1, new_picheight*new_picwidth*3, outfp);
    
//   free(buffer);
//   free(re_buffer);
}
void mirror(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);
    
    fseek(infp, 54L, 0);
    fread(buffer, 3*picwidth*picheight, 1, infp);
    //    fread(&new_picwidth, 4, 1, outfp);
    //    fread(&new_picheight, 4, 1, outfp);
    fseek(outfp, 54L, 0);
    int r1,g1,b1;
    r1 = b1 = g1 = 255;
    for (int i=0;i< picheight;i++)
        for (int j = 0; j < picwidth; j++){
            //r1 = b1 = g1 = 105;
            
            b1 = buffer[3*(i)*picwidth+3*(picwidth - j)];
            g1 = buffer[3*(i)*picwidth+3*(picwidth - j)+1];
            r1 = buffer[3*(i)*picwidth+3*(picwidth - j)+2];
            
                        //write the new color bmp
            writeRGB(&b1, &g1, &r1, outfp);
            
        }
    free(buffer);

}
void translation(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp,int dx,int dy){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);
    
    fseek(infp, 54L, 0);
    fread(buffer, 3*picwidth*picheight, 1, infp);
    fseek(outfp, 18L, 0);
    int new_picwidth = picwidth + abs(dx)+3;//+3是为了下面的对齐
    new_picwidth = new_picwidth/4*4;//4字节对齐
    int new_picheight = picheight + abs(dy)+3;//+3是为了下面的对齐
    new_picheight = new_picheight/4*4;//4字节对齐
    fwrite(&new_picwidth, 4, 1, outfp);
    fwrite(&new_picheight, 4, 1, outfp);
    fseek(outfp,34L,0);//图像数据区的大小
    int biSizeImage = 3*new_picheight*new_picwidth;
    fwrite(&biSizeImage, 4, 1, outfp);
    unsigned int bfSize = biSizeImage+54;
    fseek(outfp, 2L, 0);
    fwrite(&bfSize, 4, 1, outfp);
//    fread(&new_picwidth, 4, 1, outfp);
//    fread(&new_picheight, 4, 1, outfp);
//    std::cout<<std::endl<<"$"<<new_picwidth<<" "<<new_picheight;
    fseek(outfp, 54L, 0);
    int r1,g1,b1;
    r1 = b1 = g1 = 255;
    for (int i=0;i<new_picheight;i++)
        for (int j = 0; j < new_picwidth; j++){
            
            int op = 0;
            if (dx>=0&&dy>=0) {
                op=0;
            }else if(dx<0&&dy>0){
                op = 1;
            }else if(dx>0&&dy<0){
                op = 2;
            }else if(dx<0&&dy<0){
                op = 3;
            }
            
            switch (op) {
                case 0:
                    if (i<dy||j<dx) {
                        b1 = r1 = g1 = 0;
                    }
                    else {
                        b1 = buffer[3*(i-dy)*picwidth+3*(j-dx)];
                        g1 = buffer[3*(i-dy)*picwidth+3*(j-dx)+1];
                        r1 = buffer[3*(i-dy)*picwidth+3*(j-dx)+2];
                    }

                    break;
                case 1:
                    if (i<dy||j>=picwidth) {
                        b1 = r1 = g1 = 0;
                    }
                    else {
                        b1 = buffer[3*(i-dy)*picwidth+3*(j)];
                        g1 = buffer[3*(i-dy)*picwidth+3*(j)+1];
                        r1 = buffer[3*(i-dy)*picwidth+3*(j)+2];
                    }
                    break;
                case 2:
                    if (i>=picheight||j<dx) {
                        b1 = r1 = g1 = 0;
                    }
                    else {
                        b1 = buffer[3*(i)*picwidth+3*(j-dx)];
                        g1 = buffer[3*(i)*picwidth+3*(j-dx)+1];
                        r1 = buffer[3*(i)*picwidth+3*(j-dx)+2];
                    }
                    break;
                case 3:
                    if (i>picheight||j>=picwidth) {
                        b1 = r1 = g1 = 0;
                    }
                    else {
                        b1 = buffer[3*(i)*picwidth+3*(j)];
                        g1 = buffer[3*(i)*picwidth+3*(j)+1];
                        r1 = buffer[3*(i)*picwidth+3*(j)+2];
                    }

                    break;
                default:
                    break;
            }

            //write the new color bmp
            writeRGB(&b1, &g1, &r1, outfp);
            
        }
    free(buffer);
    
}
//Image logarithmic operation
//this version is for the grey bmp image;
void logarithmic_operation(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);
    //acquire the gery image to handle
    grey(infp, picheight, picwidth, mingrey, maxgrey, outfp);
    
    fseek(outfp, 54L, 0);
    fseek(infp, 54L, 0);
    //read the grey value into the buffer
    fread(buffer, 3*picwidth*picheight, 1, outfp);
    unsigned char y = 0;
    fseek(outfp, 54L, 0);
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){

            y = buffer[3*i*picwidth+3*j];
             y = 255 * (log2(y+1)-log2(*mingrey+1))/(log2(*maxgrey+1)-log2(*mingrey+1));//get the new y
            //write
            if (y>255) {
                y=255;
            }
            if (y<0) {
                y=0;
            }
            
            fwrite(&y, 1, 1, outfp);
            fwrite(&y, 1, 1, outfp);
            fwrite(&y, 1, 1, outfp);
        }
    free(buffer);
    
}

//Image logarithmic operation
//this version is for the color bmp image;
void logarithmic_operation_color(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picwidth*picheight);
    //get the yuv image to handle
    rgb2yuv(infp, picheight, picwidth, mingrey, maxgrey, outfp);
    
    fseek(outfp, 54L, 0);
    
    //read the grey value into the buffer
    fread(buffer, 3*picwidth*picheight, 1, outfp);
    int y, u,v;
    y=u=v=0;
    fseek(outfp, 54L, 0);
    int r1,g1,b1;
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            
            y = buffer[3*i*picwidth+3*j];
            u = buffer[3*i*picwidth+3*j+1];
            v = buffer[3*i*picwidth+3*j+2];
            
            y = 255 * (log2(y+1)-log2(*mingrey+1))/(log2(*maxgrey+1)-log2(*mingrey+1));//get the new y
            //change back to rgb
            r1= y+1.4075 *(v-128);
            g1= y-0.3455*(u-128)-0.7169*(v-128);
            b1= y+1.779 * (u-128);

            //handle overflow
            if (r1>255) r1=255;
            if (g1>255) g1=255;
            if (b1>255) b1=255;
            if (r1<0)   r1=0;
            if (g1<0)   g1=0;
            if (b1<0)   b1=0;
                
            //write the new color bmp
            writeRGB(&b1, &g1, &r1, outfp);

        }
    free(buffer);
}

//Histogram equalization
//this version is for the grey bmp image;
void Histogram_equalization(FILE *infp,int picheight,int picwidth,int* mingrey,int*maxgrey,FILE*outfp){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picheight*picwidth);
    grey(infp, picheight, picwidth, mingrey, maxgrey, outfp);
    
    fseek(outfp, 54L, 0);
    fread(buffer, 1, 3*picheight*picwidth, outfp);
    double OriginGrey[256]={0};
    double TransformGray[256]={0};
    double temp[256]={0};
    for (int j = 0; j<picheight; j++) {
        for (int k=0 ; k<picwidth; k++) {
            OriginGrey[ buffer[3*j*picwidth+3*k] ]++;//0~255 each level's amount
        }
    }
    for (int i = 0; i<256; i++) {
        OriginGrey[ i ] /= (picwidth*picheight);//0~255 each level's probability
    }
    TransformGray[0] = OriginGrey[0];
    for (int i =1; i<256; i++) {
        TransformGray[i]=TransformGray[i-1]+OriginGrey[i];//0~255 each level's 累加 probability
    }
    for (int i = 0; i<256; i++) {
        temp[i] = round(255*TransformGray[i]);//get the result grey value
    }
    for (int j = 0; j<picheight; j++) {
        for (int k=0 ; k<picwidth; k++) {
            buffer[3*j*picwidth+3*k] = temp[buffer[3*j*picwidth+3*k]];//change the image's grey to the result grey value
        }
    }
    int y = 0;
    fseek(outfp, 54L, 0);
    for (int j=0;j<picheight;j++)
        for (int k = 0; k < picwidth; k++){
            
            y = buffer[3*j*picwidth+3*k];
            writeRGB(&y, &y, &y, outfp);
        }
    free(buffer);
    
}
//Histogram equalization
//this version is for the color bmp image;
void Histogram_equalization_color(FILE *infp,int picheight,int picwidth,int*mingrey,int*maxgrey,FILE*outfp){
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*picheight*picwidth);

    int r1,g1,b1;
    int y,u,v;
    
    rgb2yuv(infp, picheight, picwidth, mingrey, maxgrey, outfp);

    fseek(outfp, 54L, 0);
    fread(buffer, 1, 3*picheight*picwidth, outfp);
    double OringeGray[256]={0};
    double TransformGray[256]={0};
    double temp[256]={0};
    for (int j = 0; j<picheight; j++) {
        for (int k=0 ; k<picwidth; k++) {
            OringeGray[ buffer[3*j*picwidth+3*k] ]++;//0~255 each level's amount
        }
    }
    for (int i = 0; i<256; i++) {
        OringeGray[ i ] /= (picwidth*picheight);//0~255 each level's probability
    }
    TransformGray[0] = OringeGray[0];
    for (int i =1; i<256; i++) {
        TransformGray[i]=TransformGray[i-1]+OringeGray[i];//0~255 each level's 累加 probability
    }
    for (int i = 0; i<256; i++) {
        temp[i] = round(255*TransformGray[i]);//get the result grey value
    }
    for (int j = 0; j<picheight; j++) {
        for (int k=0 ; k<picwidth; k++) {
            buffer[3*j*picwidth+3*k] = temp[buffer[3*j*picwidth+3*k]];//change the image's grey
        }
    }
    fseek(outfp, 54L, 0);
    for (int j=0;j<picheight;j++)
        for (int k = 0; k < picwidth; k++){
            
            y = buffer[3*j*picwidth+3*k];
            u = buffer[3*j*picwidth+3*k+1];
            v = buffer[3*j*picwidth+3*k+2];
            //change back to rgb
            r1= y+1.4075 *(v-128);
            g1= y-0.3455*(u-128)-0.7169*(v-128);
            b1= y+1.779 * (u-128);
            //handle overflow
            if (r1>255) r1=255;
            if (g1>255) g1=255;
            if (b1>255) b1=255;
            if (r1<0)   r1=0;
            if (g1<0)   g1=0;
            if (b1<0)   b1=0;
            writeRGB(&b1, &g1, &r1, outfp);
        }
    free(buffer);

}


void bmp_binarization(int picheight, int picwidth,int threshold,FILE* infp,FILE* outfp)
{
    int i,j;
    unsigned char y_b;
    fseek(outfp, 54L, 0);
    fseek(infp, 54L, 0);

    for (i=0;i<picheight;i++)
        for (j = 0; j < picwidth; j++){
//            fread(&y_b, 1, 1, infp);
//            fread(&y_b, 1, 1, infp);
//            fread(&y_b, 1, 1, infp);
            readRGB(&y_b, &y_b, &y_b, infp);
            if (y_b>threshold) {
                y_b = 255;
            }else y_b = 0;
            writeRGB(&y_b, &y_b, &y_b, outfp);
        }
}

//the strture unit is a black cross
//the model = 1 :erosion model=0:delation [in blcak]
void erosion_delation(unsigned char* data, int height, int width,unsigned char* data2,int model)
{
    int i, j,sum, flag;
    for (int i =0; i<3*width; i++) {
        data2[i]=data[i];
        data2[3*(height-1)*width+i]=data[3*(height-1)*width+i];
    }
    
    memset(data2, 0, height*width*3);
    
    for(i = 0;i < height ;i++)
    {
        for(j = 0;j < width;j++)
        {
            if (i==0||i==height-1||j==0||j==width-1) {
                //handle the pixels around the boundary
                data2[3*i * width + 3*j] = 255;
                data2[3*i * width + 3*j+1] = 255;
                data2[3*i * width + 3*j+2] =255;
            }
            else{
                flag = 1;
                //if all pixels in the cross centered with the target pixel are black ,the target pixel will be black; else be white
                //use black to erode <=> use white to delation
                //use white to erode <=> use black to delation
                //the model = 1 :erosion model=0:delation both use black
                if (model==1) {
                    sum =!!(data[3*i*width+3*(j-1)]==0 && data[3*i*width+3*(j+1)]==0 &&data[3*i*width+3*j]==0&&data[3*(i-1)*width+3*j]==0 && data[3*(i+1)*width+3*j]==0);
                }
                else  {sum =!!(data[3*i*width+3*(j-1)]==0||data[3*i*width+3*(j+1)]==0||data[3*(i-1)*width+3*(j)]==0||data[3*(i+1)*width+3*(j)]==0||data[3*i*width+3*j]==0);
                }
                
                
                if (sum) {
                    flag = 0;
                }
                if(flag == 0 )
                {   //be black
                    data2[3*i * width + 3*j] = 0;
                    data2[3*i * width + 3*j+1] = 0;
                    data2[3*i * width + 3*j+2] =0;
                    
                }
                else
                {   //be white
                    data2[3*i * width + 3*j] = 255;
                    data2[3*i * width + 3*j+1] = 255;
                    data2[3*i * width + 3*j+2] = 255;
                }
            }
        }
    }
}

void opening(unsigned char* data, int height, int width,unsigned char* data2){
    //erosion first then delation
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*width*height);
    erosion_delation(data, height, width,buffer ,1);
    erosion_delation(buffer, height, width,data2,0);
    free(buffer);
}

void closing(unsigned char* data, int height, int width,unsigned char* data2){
    //delation first then erosion
    unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned char)*3*width*height);
    erosion_delation(data, height, width,buffer,0);
    erosion_delation(buffer, height, width,data2 ,1);
    
    free(buffer);
}

//the way to find the suitable threshold:
//1.find minimum grey value t1;
//2.find the average grey m1 value of the pixels whose grey value is smaller than t1;the average grey m2 value of the pixels whose grey value is greater than t1
//3.t1' = (m1+m2) /2
//4.t1'?t1 equal then end; not equal ,then go back to 2

int FindThreshold(int mingrey,int maxgrey,int picheight,int picwidth,FILE* infp){
    //fseek(infp, 54L, 0);
    unsigned char ty;
    int a1,a2;//a1 is the result
    double t1,t2;// use double to get more accurate result and avoid oveeflow
    t1 =t2 =0;
    t1 = (mingrey+maxgrey)/2;//set the initial value to be the middle number in order to iterate more conveniently
    t2 = maxgrey;
    a1 = t1;
    a2 = t2;
    while (a1!=a2) {
        a2 = a1;
        t2 = t1;
        // use double to get more accurate result and avoid overflow
        double totalmin,totalmax;
        double countmin,countmax;
        double m1,m2;
        countmax = countmin = 0;
        totalmax = totalmin = 0;
        fseek(infp, 54L, 0);
        for (int i=0;i<picheight;i++)
            for (int j = 0; j < picwidth; j++){
                readRGB(&ty, &ty, &ty, infp);
                
                if (ty>a1) {
                    totalmax+=ty;
                    countmax++;
                }else {
                    totalmin+=ty;
                    countmin++;
                }
            }
        m1 = totalmin/countmin;
        m2 = totalmax/countmax;
        t1 = (m1+m2)/2;//iteration
        a1 = (int)t1;//get the int value
    }
    return a1;
}

void grey(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp){
    
    int high = 0;
    int low = 255;
    int i, j;
    unsigned char r,g,b;
    unsigned char y;
    
    fseek(infp, 54L, 0);
    //find the max and min grey value to rearrange the grey value
    for (i=0;i<picheight;i++)
        for (j = 0; j < picwidth; j++){
            readRGB(&b, &g, &r, infp);
            
            //rgb->y
            y=1.0* 0.2990*r+0.5869*g+0.1140*b;
//            if (r==0&&g==0&&b==0) {
//                printf(">\n");
//            }
  
            if (y>high) {
                high = y;
            }
            if (y<low) {
                low = y;

            }
            
        }
    *mingrey = low;
    *maxgrey = high;
    fseek(infp, 54L, 0);
    for (i=0;i<picheight;i++)
        for (j = 0; j < picwidth; j++){
            readRGB(&b, &g, &r, infp);
            //rgb->yuv
            y= 0.2990*r+0.5869*g+0.1140*b;
            
            //rearrange y to project into [0,255]
            y=1.0*(y-low)/(high-low)*255;
            //write the grey image
            if (y>=255) {
                y=255;
            }
            writeRGB(&y, &y, &y, outfp);
        }
}


void rgb2yuv(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp){
    
 
    unsigned char r,g,b;
    int r1,g1,b1;
    int y,u,v;
    int high = 0;
    int low = 255;
    // WriteBmp(outfp, fp , r1, g1, b1, picheight, picwidth);
    fseek(infp, 54L, 0);
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            
            readRGB(&b, &g, &r, infp);
            //rgb->yuv
            y= 0.2990*r+0.5869*g+0.1140*b;
            if (y>high) high = y;
            if (y<low)  low = y;
            
        }
    fseek(infp, 54L, 0);
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            readRGB(&b, &g, &r, infp);
            //to avoid rgb overflows，use  r1,b1,g1 of int type to store their value  for the calculation
            r1 = r;
            b1 = b;
            g1 = g;
            
            //rgb->yuv
            y= 0.2990*r1+0.5869*g1+0.1140*b1;
            u= -0.1687*r1-0.3313*g1+0.5000*b1 +128;
            v= 0.5000*r1-0.4187*g1-0.0813*b1 + 128;
            
            if (y<*mingrey) {
                *mingrey = y;
            }
            if (y>*maxgrey) {
                *maxgrey = y;
            }
            //rearrange y to project into [0,255]
            y=1.0*(y-low)/(high-low)*255;
            if (y>=255) {
                y=255;
            }
            //write the new color bmp in yuv
            writeRGB(&y, &u, &v, outfp);
            
        }
    

}

void yuv2rgb(FILE* infp,int picheight,int picwidth,int* mingrey,int* maxgrey,FILE* outfp){
    
    
    unsigned char r,g,b;
    int r1,g1,b1;
    unsigned char y,u,v;

    // WriteBmp(outfp, fp , r1, g1, b1, picheight, picwidth);
    fseek(infp, 54L, 0);
    for (int i=0;i<picheight;i++)
        for (int j = 0; j < picwidth; j++){
            readRGB(&y, &u, &v, infp);
            //to avoid rgb overflows，use  r1,b1,g1 of int type to store their value  for the calculation
            //yuv->rgb
            r1= y+1.4075 *(v-128);
            g1= y-0.3455*(u-128)-0.7169*(v-128);
            b1= y+1.779 * (u-128);
            if (r1>255) r1=255;
            if (g1>255) g1=255;
            if (b1>255) b1=255;
            if (r1<0)   r1=0;
            if (g1<0)   g1=0;
            if (b1<0)   b1=0;
            r=r1;
            b=b1;
            g=g1;
            
            //write the new color bmp in yuv
            writeRGB(&b, &g, &r, outfp);
            
        }
    
}

void yuv2rgb_value(unsigned char* y,unsigned char*u,unsigned char*v,unsigned char*r,unsigned char*g,unsigned char*b){
    int r1,b1,g1;
    r1= *y+1.4075 *(*v-128);
    g1= *y-0.3455*(*u-128)-0.7169*(*v-128);
    b1= *y+1.779 * (*u-128);
    //            y= 0.2990*r1+0.5869*g1+0.1140*b1;
    //            u= -0.1687*r1-0.3313*g1+0.5000*b1 +128;
    //            v= 0.5000*r1-0.4187*g1-0.0813*b1 + 128;

    if (r1>255) r1=255;
    if (g1>255) g1=255;
    if (b1>255) b1=255;
    if (r1<0)   r1=0;
    if (g1<0)   g1=0;
    if (b1<0)   b1=0;
    *r=r1;
    *b=b1;
    *g=g1;
}
