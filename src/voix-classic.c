/*********************************************************************
  
  VoiX -- An Open-Source Vocal Eliminator
  (http://vocaleliminator.sourceforge.net)
  Version 1.0.0 Classic
  Written by Zhang Zhichao(zzcworld@hotmail.com)
  License: GNU General Public License (GPL) Ver 2.0 & Above
  
  VoiX is a voal elimination software and can eliminate vocals from
  standard PCM Wave files. VoiX uses partial code from "dsp_centercut"
  plugin (http://www.moitah.net), which is distributed under the GNU
  General Public License(GPL). You can redistribute and/or modify this
  software under the GPL license.
  
  * Due to my lack of time, I will not continue developing VoiX. I
    hope my software is useful for you, but I don't garantee that it
    works well. For any questions please send me an e-mail.
    
  Acknowledgement
  * dsp_centercut (http://www.moitah.net)
  * FFT code from DSP Dimension (http://www.dspdimension.com)
  * LP and HP Filter code from Robert Bristow-Johnson
    (rbj@audioimagination.com)
  * The "Center Cut" Algorithm used in Virtual Dub
  
*********************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define pi 3.141592653589

#define VERSION "1.0.0 Classic"
#define kWindowSize 32768
#define kOverlapCount 4

typedef struct{
        unsigned short FormatTag;         
		unsigned short Channels;         
		unsigned long SampleRate;   
		unsigned long AvgBytesPerSec;  
		unsigned short BlockAlign;
        unsigned short Bits;
        } fmtChunk;

typedef struct{
        char ID[4];
        unsigned long Length;
        } wHeader;

/*
  Set two biquad filters: a low-pass filter and a high-pass filter
  See "Cookbook formulae for audio EQ biquad filter coefficients"
  (audioeq.txt) for details.
*/
void SetLP(float samplerate, float cutoff, double *a,double *b){
     double w0, Q, alpha;
     w0 = 2*pi*cutoff/samplerate;
     Q = 1.0;
     //Q = 1.0/(2.0*sinh(ln2/2.0*BW*w0/sin(w0)));
     alpha = sin(w0)/(2.00*Q);
     *b = (1.00-cos(w0))/2.00;
     *(b+1) = 1.00-cos(w0);
     *(b+2) = (1.00-cos(w0))/2.00;
     *a = 1.00+alpha;
     *(a+1) = -2.00*cos(w0);
     *(a+2) = 1.00-alpha;
}

void SetHP(float samplerate,float cutoff, double *a,double *b){
     double w0, Q, alpha;
     w0 = 2*pi*cutoff/samplerate;
     Q = 1.0;
     //Q = 1.0/(2.0*sinh(ln2/2.0*BW*w0/sin(w0)));
     alpha = sin(w0)/(2.00*Q);
     *b = (1.00+cos(w0))/2.00;
     *(b+1) = -(1.00+cos(w0));
     *(b+2) = (1.00+cos(w0))/2.00;
     *a = 1.00+alpha;
     *(a+1) = -2.00*cos(w0);
     *(a+2) = 1.00-alpha;
}

inline int BlockRead(double *blockdata, unsigned short bits, unsigned short FormatTag, FILE *fp){
    unsigned char block8;
    signed short block16;
    unsigned char block24[3];
    float block32f;
    signed long block32d;
    int stat;
    
    if ((bits==8)&&(FormatTag==1)) {
                 stat = fread(&block8,sizeof(block8),1,fp);
                 *blockdata = (double)block8;
                 *blockdata = *blockdata - 128.00;
                 return(stat);
                 }                
    if ((bits==16)&&(FormatTag==1)) {
                  stat = fread(&block16,sizeof(block16),1,fp);
                  *blockdata = (double)block16;
                  return(stat);
                  }
    if ((bits==24)&&(FormatTag==1)) {
                  stat = fread(&block24,sizeof(block24),1,fp);
                  *blockdata = block24[2]*65536.00 + block24[1]*256.00 + block24[0];
                  if (*blockdata >= 8388608.00) *blockdata -= 16777216.00;
                  return(stat);
                  }
    if ((bits==32)&&(FormatTag==1)) {
                  stat = fread(&block32d,sizeof(block32d),1,fp);
                  *blockdata = (double)block32d;
                  return(stat);
                  }                 
    if ((bits==32)&&(FormatTag==3)) {
                  stat = fread(&block32f,sizeof(block32f),1,fp);
                  *blockdata = (double)block32f;
                  return(stat);
                  }
    return(0);
}

inline int BlockWrite(double blockdata, short bits, unsigned short FormatTag, FILE *fp){
    unsigned char block8;
    signed short block16;
    unsigned char block24[3];
    float block32f;
    signed long block32d;
    int stat;

    if ((bits==8)&&(FormatTag==1)) {
                 if (blockdata >= 127.00) blockdata = 127.00;
                 if (blockdata <= -128.00) blockdata = -128.00;
                 blockdata = blockdata + 128.00;
                 block8 = (char)blockdata;
                 stat = fwrite(&block8,sizeof(block8),1,fp);
                 return(stat);
                 }
    if ((bits==16)&&(FormatTag==1)) {
                  if (blockdata >= 32767.00) blockdata = 32767.00;
                  if (blockdata <= -32768.00) blockdata = -32768.00;
                  block16 = (short)blockdata;
                  stat = fwrite(&block16,sizeof(block16),1,fp);
                  return(stat);
                  }
    if ((bits==24)&&(FormatTag==1)) {
                  if (blockdata >= 8388607.00) blockdata = 8388607.00;
                  if (blockdata <= -8388608.00) blockdata = -8388608.00;
                  if (blockdata <= 0.00) blockdata += 16777216.00;
                  block24[0] = ((long)blockdata)%256;
                  block24[1] = (((long)blockdata)/256)%256;
                  block24[2] = ((long)blockdata)/65536;
                  stat = fwrite(&block24,sizeof(block24),1,fp);
                  return(stat);
                  }
    if ((bits==32)&&(FormatTag==1)) {
                  if (blockdata >= 2147483647.00) blockdata = 2147483647.00;
                  if (blockdata <= -2147483648.00) blockdata = -2147483648.00;
                  block32d = (long)blockdata;
                  stat = fwrite(&block32d,sizeof(block32d),1,fp);
                  return(stat);
                  }
    if ((bits==32)&&(FormatTag==3)) {
                  block32f = (float)blockdata;
                  stat = fwrite(&block32f,sizeof(block32f),1,fp);
                  return(stat);
                  }
    return(0);
}
int main(int argc, char *argv[]){

    FILE *infile, *outfile;
    char infilename[260], outfilename[260];
    char WaveId[4];
    wHeader infileHeader;
    fmtChunk infileChunk;
    unsigned long Count, DataSize, tmpDataSize, tmpRIFFSize, FilePos;

    float LPCutoff, HPCutoff, Invert, VocalPan, Gain;
    float ChMix[4];
    double a1[3],b1[3],a2[3],b2[3];
    double LP_L_In[3], LP_L_Out[3], HP_L_In[3], HP_L_Out[3];
    double LP_R_In[3], LP_R_Out[3], HP_R_In[3], HP_R_Out[3];    
    double Left, Right, NewLeft, NewRight;
    unsigned long i, j, k;
    short percentage = 0;
    
    /* Display version and usage */
    printf("VoiX Version %s\n(http://vocaleliminator.sourceforge.net)\n", VERSION);

    if (argc==1) {printf("\nUsage: %s <infile> [outfile] [options]\n", argv[0]);
                  printf("\nOptions: <lowpass> <highpass> <invert> <pan> <gain>\n");
                  printf("\n<lowpass>: Cutoff frequency (in Hz) for low-pass filter\n");
                  printf("<highpass>: Cutoff frequency (in Hz) for high-pass filter\n");
                  printf("<invert>: Invert right channel to make it mono, 0=off 1=on\n");                  
                  printf("<pan>: Pan for Vocal Track: 0=center -100=left 100=right\n");
                  printf("<gain>: Gain (in dB) for output file\n");                  
                  printf("\nExample: %s input.wav output.wav 200 8000\n", argv[0]);
                  return(0);}
    strcpy(infilename,argv[1]); /* Copy input filename from first parameter */
    
    /* For Windows Drag & Drop*/
    char *FileExt = ".voix.wav";
    if (argc==2) {
                 strcpy(outfilename, infilename); strcat(outfilename,FileExt); 
                 } else strcpy(outfilename,argv[2]);
    /* Windows Drag & Drop End*/
    
    printf("Input: %s ", infilename);
    if (argc>=4) LPCutoff = atof(argv[3]); else LPCutoff = 200.0;
    if (argc>=5) HPCutoff = atof(argv[4]); else HPCutoff = 8000.0;
    if (argc>=6) Invert = atof(argv[5]); else Invert = 0.0;
    if (argc>=7) VocalPan = atof(argv[6]); else VocalPan = 0.0;
    if ((VocalPan<=-100.0)||(VocalPan>=100.0)) VocalPan = 0.0;
    if (argc>=8) Gain = atof(argv[7]); else Gain = 0.0;
    if ((Gain<-20.0)||(Gain>20.0)) Gain = 0.0;

    /* Read Header from Input Wave File*/
    if ((infile = fopen(infilename,"rb+")) == NULL) {
        printf("\n\nError(1): Could not open input file \"%s\".\n",infilename); return(1);}
    
    fread(&infileHeader, sizeof(infileHeader), 1, infile);
    if ((infileHeader.ID[0]!='R')||(infileHeader.ID[1]!='I')
       ||(infileHeader.ID[2]!='F')||(infileHeader.ID[3]!='F')) {
         printf("\n\nError(3): \"%s\" is not a standrad RIFF file.\n",infilename);
           fclose(infile); return(3);}

    fread(&WaveId, sizeof(WaveId), 1, infile);
    if ((WaveId[0]!='W')||(WaveId[1]!='A')||(WaveId[2]!='V')||(WaveId[3]!='E')) {
        printf("\n\nError(4): Could not find correct WAVE header form \"%s\".\n", infilename);
           fclose(infile); return(4);}

    fread(&infileHeader, sizeof(infileHeader), 1, infile);
    if ((infileHeader.ID[0]!='f')||(infileHeader.ID[1]!='m')
      ||(infileHeader.ID[2]!='t')||(infileHeader.ID[3]!=' ')) {
        printf("\n\nError(5): Could not find correct chunk header form \"%s\".\n", infilename);
           fclose(infile); return(5);}
    
    fread(&infileChunk, sizeof(infileChunk), 1, infile);
    fseek(infile, ftell(infile) + infileHeader.Length - 16, 0);
    /* Display the input file format*/
    if (infileChunk.FormatTag == 1) printf("(Windows PCM");
    if (infileChunk.FormatTag == 3) printf("(IEEE float");
    printf(", %d Hz, %d bit, ", infileChunk.SampleRate, infileChunk.Bits);
    if (infileChunk.Channels == 1) printf("mono)\n"); else
    if (infileChunk.Channels == 2) printf("stereo)\n"); else
    printf("%d channels)\n", infileChunk.Channels);
    /* end */
    if ((infileChunk.FormatTag!=0x0001)&&(infileChunk.FormatTag!=0x0003)) {
        printf("\nError(6): \"%s\" is not a PCM or IEEE float wave file.\n",infilename);
          fclose(infile); return(6);}
    if (infileChunk.Channels!=2) {
        printf("\nError(7): \"%s\" is not stereo.\n",infilename);
          fclose(infile); return(7);}
    if ((infileChunk.Bits!=8)&&(infileChunk.Bits!=16)&&(infileChunk.Bits!=24)&&(infileChunk.Bits!=32))
        { printf("\nError(8): VoiX could not process %d bit wave files.\n", infileChunk.Bits);
          fclose(infile); return(8);}
    if (infileChunk.Bits/8*infileChunk.Channels!=infileChunk.BlockAlign)
        { printf("\nError(9): VoiX could not process non-standrad wave files.\n");
          fclose(infile); return(9);}
    fread(&infileHeader, sizeof(infileHeader), 1, infile);
    while (((infileHeader.ID[0]!='d')||(infileHeader.ID[1]!='a')||
          (infileHeader.ID[2]!='t')||(infileHeader.ID[3]!='a'))
    &&(!feof(infile))) {
    fseek(infile, ftell(infile) + infileHeader.Length, 0); 
    fread(&infileHeader, sizeof(infileHeader), 1, infile);
    }
    if (feof(infile)) {
       printf("\nError(10): Could not find correct data header form \"%s\".\n", infilename);
         fclose(infile); return(10);}
    printf("Output: %s\n", outfilename);
    if ((outfile = fopen(outfilename,"wb+")) == NULL) {
        printf("\nError(2): Could not open output file \"%s\".\n",outfilename);
          fclose(infile); return(2);}
    
    /* Write Header to Output Wave File*/
    DataSize = infileHeader.Length;
    tmpDataSize = 0;
    strcpy(infileHeader.ID,"RIFF"); infileHeader.Length = tmpDataSize+36;
    fwrite(&infileHeader,sizeof(infileHeader),1,outfile);
    strcpy(WaveId,"WAVE"); fwrite(WaveId,sizeof(WaveId),1,outfile);
    strcpy(infileHeader.ID,"fmt "); infileHeader.Length = 16;
    fwrite(&infileHeader,sizeof(infileHeader),1,outfile);
    fwrite(&infileChunk,sizeof(infileChunk),1,outfile);
    strcpy(infileHeader.ID,"data"); infileHeader.Length = tmpDataSize;
    fwrite(&infileHeader,sizeof(infileHeader),1,outfile);    
    
    if ((LPCutoff < 1.00) || (LPCutoff > infileChunk.SampleRate/2 - 1)) LPCutoff = 1.00;
    if ((HPCutoff < 1.00) || (HPCutoff > infileChunk.SampleRate/2 - 1)) HPCutoff = infileChunk.SampleRate/2 - 1;
    if (HPCutoff < LPCutoff) HPCutoff = LPCutoff;
    printf("(i)Lowpass: %.1f Hz, Highpass: %.1f Hz\n", LPCutoff, HPCutoff);
    SetLP(infileChunk.SampleRate,LPCutoff,a1,b1);
    SetHP(infileChunk.SampleRate,HPCutoff,a2,b2);
    for (i=0;i<=2;i++) {
      LP_L_In[i]=0; LP_L_Out[i]=0;
      HP_L_In[i]=0; HP_L_Out[i]=0;
      LP_R_In[i]=0; LP_R_Out[i]=0;
      HP_R_In[i]=0; HP_R_Out[i]=0;
    }
    printf("(i)Invert: ");
    if (Invert == 1.00) printf("On"); else printf("Off");
    printf(" Vocal Pan: "); printf("%.1f", VocalPan);
    VocalPan = (VocalPan+100.00)/2.00;
    printf(" Gain: %.1f dB\n", Gain);
    Gain = pow(10.0,Gain/20.0);
    
    if (VocalPan <= 50.0) {
    ChMix[0] = 1.0;
    ChMix[1] = -(VocalPan/(100.0-VocalPan));
    ChMix[2] = 1.0;
    ChMix[3] = -(VocalPan/(100.0-VocalPan));
    }
    if (VocalPan > 50.0) {
    ChMix[0] = -((100.0-VocalPan)/VocalPan);
    ChMix[1] = 1.0;
    ChMix[2] = -((100.0-VocalPan)/VocalPan);
    ChMix[3] = 1.0;
    }
    if (Invert == 0.00) {
      if (VocalPan<=50.00) {
         ChMix[2] = -ChMix[2]; ChMix[3] = -ChMix[3];
      }else {
         ChMix[0] = -ChMix[0]; ChMix[1] = -ChMix[1];         
      }         
    }

    printf("(!)Press [Ctrl+C] to exit\n(i)Process: ");
    Count = 0;
    for (i = 1; i <= DataSize/infileChunk.BlockAlign; i++) {
        
        if (feof(infile)) {
                        printf("\nError(11): EOF(End of file) of input file detected.\n");
                        fclose(infile); fclose(outfile); return(11);
                        }
        BlockRead(&Left, infileChunk.Bits, infileChunk.FormatTag, infile);
        BlockRead(&Right, infileChunk.Bits, infileChunk.FormatTag, infile);      
      
        
        if (percentage<i*20/(DataSize/infileChunk.BlockAlign)) {
          tmpDataSize = Count*infileChunk.BlockAlign;
          FilePos = ftell(outfile);
          tmpRIFFSize = tmpDataSize + 36;
          fseek(outfile,4,0); fwrite(&tmpRIFFSize,sizeof(tmpRIFFSize),1,outfile);
          fseek(outfile,40,0); fwrite(&tmpDataSize,sizeof(tmpDataSize),1,outfile);        
          fseek(outfile,FilePos,0);
          percentage= i*20/(DataSize/infileChunk.BlockAlign);
          printf("."); if (percentage%5==0) printf("%d%c",percentage*5,'%');

        }
        
        Count++;
        LP_L_In[2]=Left; LP_R_In[2]=Right;
        LP_L_Out[2]=(b1[0]/a1[0])*LP_L_In[2]+(b1[1]/a1[0])*LP_L_In[1]+(b1[2]/a1[0])*LP_L_In[0]-(a1[1]/a1[0])*LP_L_Out[1]-(a1[2]/a1[0])*LP_L_Out[0];
        LP_R_Out[2]=(b1[0]/a1[0])*LP_R_In[2]+(b1[1]/a1[0])*LP_R_In[1]+(b1[2]/a1[0])*LP_R_In[0]-(a1[1]/a1[0])*LP_R_Out[1]-(a1[2]/a1[0])*LP_R_Out[0];

        HP_L_In[2]=Left; HP_R_In[2]=Right;
        HP_L_Out[2]=(b2[0]/a2[0])*HP_L_In[2]+(b2[1]/a2[0])*HP_L_In[1]+(b2[2]/a2[0])*HP_L_In[0]-(a2[1]/a2[0])*HP_L_Out[1]-(a2[2]/a2[0])*HP_L_Out[0];
        HP_R_Out[2]=(b2[0]/a2[0])*HP_R_In[2]+(b2[1]/a2[0])*HP_R_In[1]+(b2[2]/a2[0])*HP_R_In[0]-(a2[1]/a2[0])*HP_R_Out[1]-(a2[2]/a2[0])*HP_R_Out[0];
           
        NewLeft=(ChMix[0]*(Left-LP_L_Out[2]-HP_L_Out[2])+ChMix[1]*(Right-LP_R_Out[2]-HP_R_Out[2]))*Gain+LP_L_Out[2]+HP_L_Out[2];
        NewRight=(ChMix[2]*(Left-LP_L_Out[2]-HP_L_Out[2])+ChMix[3]*(Right-LP_R_Out[2]-HP_R_Out[2]))*Gain+LP_R_Out[2]+HP_R_Out[2];
      
        LP_L_In[0]=LP_L_In[1]; LP_L_In[1]=LP_L_In[2];
        LP_R_In[0]=LP_R_In[1]; LP_R_In[1]=LP_R_In[2];
        LP_L_Out[0]=LP_L_Out[1]; LP_L_Out[1]=LP_L_Out[2];
        LP_R_Out[0]=LP_R_Out[1]; LP_R_Out[1]=LP_R_Out[2]; 
      
        HP_L_In[0]=HP_L_In[1]; HP_L_In[1]=HP_L_In[2];
        HP_R_In[0]=HP_R_In[1]; HP_R_In[1]=HP_R_In[2];
        HP_L_Out[0]=HP_L_Out[1]; HP_L_Out[1]=HP_L_Out[2];
        HP_R_Out[0]=HP_R_Out[1]; HP_R_Out[1]=HP_R_Out[2];  
      
        BlockWrite(NewLeft, infileChunk.Bits, infileChunk.FormatTag, outfile);
        BlockWrite(NewRight,infileChunk.Bits, infileChunk.FormatTag, outfile);       
    }
    tmpDataSize = DataSize;
    FilePos = ftell(outfile);
    tmpRIFFSize = tmpDataSize + 36;
    fseek(outfile,4,0); fwrite(&tmpRIFFSize,sizeof(tmpRIFFSize),1,outfile);
    fseek(outfile,40,0); fwrite(&tmpDataSize,sizeof(tmpDataSize),1,outfile);        
    fseek(outfile,FilePos,0);    
    printf("\n(i)Done processing file!\n");
    fclose(infile); fclose(outfile); return(0);
}
