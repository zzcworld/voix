/*********************************************************************
  
  VoiX -- An Open-Source Vocal Eliminator
  (http://vocaleliminator.sourceforge.net)
  Version 1.0.0 Beta 5
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

#define VERSION "1.0.0 Beta 5"
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

void DoFFT(float *fftBuffer, long fftFrameSize, long sign)
/* 
	FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
	Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
	time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
	and returns the cosine and sine parts in an interleaved manner, ie.
	fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
	must be a power of 2. It expects a complex input signal (see footnote 2),
	ie. when working with 'common' audio signals our input signal has to be
	passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
	of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
{
	float wr, wi, arg, *p1, *p2, temp;
	float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
	long i, bitm, j, le, le2, k;

	for (i = 2; i < 2*fftFrameSize-2; i += 2) {
		for (bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1) {
			if (i & bitm) j++;
			j <<= 1;
		}
		if (i < j) {
			p1 = fftBuffer+i; p2 = fftBuffer+j;
			temp = *p1; *(p1++) = *p2;
			*(p2++) = temp; temp = *p1;
			*p1 = *p2; *p2 = temp;
		}
	}
	for (k = 0, le = 2; k < (long)(log(fftFrameSize)/log(2.)+.5); k++) {
		le <<= 1;
		le2 = le>>1;
		ur = 1.0;
		ui = 0.0;
		arg = pi / (le2>>1);
		wr = cos(arg);
		wi = sign*sin(arg);
		for (j = 0; j < le2; j += 2) {
			p1r = fftBuffer+j; p1i = p1r+1;
			p2r = p1r+le2; p2i = p2r+1;
			for (i = j; i < 2*fftFrameSize; i += le) {
				tr = *p2r * ur - *p2i * ui;
				ti = *p2r * ui + *p2i * ur;
				*p2r = *p1r - tr; *p2i = *p1i - ti;
				*p1r += tr; *p1i += ti;
				p1r += le; p1i += le;
				p2r += le; p2i += le;
			}
			tr = ur*wr - ui*wi;
			ui = ur*wi + ui*wr;
			ur = tr;
		}
	}
}

inline void CenterCut(float *lBuff, float *rBuff, float *cBuff, long windowsize){
   	 double lR, lI, rR, rI, cR, cI;
	 double A, B, C, D, alpha;
     long k;
     /*
         The "Center Cut" Algorithm used in Virtual Dub
         See "The "center cut" algorithm"(center.htm) for details.
     */
     *(cBuff) = 0.0; *(cBuff+1) = 0.0;
     for (k = 1; k < windowsize; k++){
         lR = *(lBuff+k*2); lI = *(lBuff+k*2+1);
         rR = *(rBuff+k*2); rI = *(rBuff+k*2+1);

         cR = lR + rR;
         cI = lI + rI;
         A = cR*cR + cI*cI;
         B = -cR*(lR+rR)-cI*(lI+rI);
         C = lR*rR+lI*rI;
         D = B*B-4*A*C;
         if ((D>=0.0) && (A!=0)) {
                      alpha = (-B-sqrt(D))/(2*A);
                      cR *= alpha;
		              cI *= alpha;
                      } else cR = cI = 0.0;
         lR = lR - cR; lI = lI - cI;
         rR = rR - cR; rI = rI - cI;
         
         *(lBuff+k*2) = lR;
         *(lBuff+k*2+1) = lI;
         *(rBuff+k*2) = rR;
         *(rBuff+k*2+1) = rI;
         *(cBuff+k*2) = cR;
         *(cBuff+k*2+1) = cI;
     }
}

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

    float LPCutoff, HPCutoff;
    double a1[3],b1[3],a2[3],b2[3];
    double LPInput[3], LPOutput[3], HPInput[3], HPOutput[3];
    float lBuf[kWindowSize], rBuf[kWindowSize], cBuf[kWindowSize];
    float CosWindow[kWindowSize], OverlapBuf[kOverlapCount-1][kWindowSize/kOverlapCount];
    float lFFTBuf[kWindowSize], rFFTBuf[kWindowSize], cFFTBuf[kWindowSize];
    double Left, Right, NewLeft, NewRight;
    unsigned long i, j, k;
    short percentage = 0;
    
    /* Display version and usage */
    printf("VoiX Version %s\n(http://vocaleliminator.sourceforge.net)\n", VERSION);

    if (argc==1) {printf("\nUsage: %s <infile> [outfile] [options]\n", argv[0]);
                  printf("\nOptions: <lowpass> <highpass>\n");
                  printf("\n<lowpass>: Cutoff frequency (in Hz) for low-pass filter\n");
                  printf("<highpass>: Cutoff frequency (in Hz) for high-pass filter\n");
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
      LPInput[i]=0; LPOutput[i]=0; HPInput[i]=0; HPOutput[i]=0;
    }
    
    memset(OverlapBuf, 0, sizeof(OverlapBuf));
    memset(cBuf, 0, sizeof(cBuf));
    
    for (i = 0; i < kWindowSize; i++) {
        CosWindow[i] = (1.0-cos(2*pi/(float)kWindowSize*(i+0.5)))/(kWindowSize/2.0)/kOverlapCount;    
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
      
        /* Display the percentage */
        if (percentage<i*20/(DataSize/infileChunk.BlockAlign)) {
          percentage= i*20/(DataSize/infileChunk.BlockAlign);
          printf("."); if (percentage%5==0) printf("%d%c",percentage*5,'%');
        }
        
        Count++;
        lBuf[Count-1] = Left;
        rBuf[Count-1] = Right;
       
        if (Count == kWindowSize) {
                             Count = kWindowSize*(kOverlapCount-1)/kOverlapCount;
                             for (k = 0; k < kWindowSize/2; k++) {
                                 lFFTBuf[2*k] = lBuf[k]*CosWindow[k];
                                 lFFTBuf[2*k+1] = lBuf[kWindowSize-k-1]*CosWindow[kWindowSize-k-1];
                                 rFFTBuf[2*k] = rBuf[k]*CosWindow[k];
                                 rFFTBuf[2*k+1] = rBuf[kWindowSize-k-1]*CosWindow[kWindowSize-k-1];
                                 }              
                             DoFFT(lFFTBuf, kWindowSize/2, -1);
                             DoFFT(rFFTBuf, kWindowSize/2, -1);
                             CenterCut(lFFTBuf, rFFTBuf, cFFTBuf, kWindowSize/2);                             
                             DoFFT(cFFTBuf, kWindowSize/2, 1);

			                 for (k=0; k < kWindowSize/2; k++) {
				                      cBuf[k] = cFFTBuf[2*k]; cBuf[kWindowSize-k-1] = cFFTBuf[2*k+1]; }                          
		                     for (k=0; k < kWindowSize/kOverlapCount; k++) {                                      

                                      LPInput[2] = cBuf[k] + OverlapBuf[0][k];
                                      HPInput[2] = cBuf[k] + OverlapBuf[0][k];
                                      LPOutput[2] = (b1[0]/a1[0])*LPInput[2]+(b1[1]/a1[0])*LPInput[1]+(b1[2]/a1[0])*LPInput[0]-(a1[1]/a1[0])*LPOutput[1]-(a1[2]/a1[0])*LPOutput[0];
                                      HPOutput[2] = (b2[0]/a2[0])*HPInput[2]+(b2[1]/a2[0])*HPInput[1]+(b2[2]/a2[0])*HPInput[0]-(a2[1]/a2[0])*HPOutput[1]-(a2[2]/a2[0])*HPOutput[0];                                      
                                      NewLeft = lBuf[k] - cBuf[k] - OverlapBuf[0][k] + LPOutput[2] + HPOutput[2];
                                      NewRight = rBuf[k] - cBuf[k] - OverlapBuf[0][k] + LPOutput[2] + HPOutput[2];
                                      
                                      BlockWrite(NewLeft, infileChunk.Bits, infileChunk.FormatTag, outfile);
                                      BlockWrite(NewRight,infileChunk.Bits, infileChunk.FormatTag, outfile);    

                                      LPInput[0] = LPInput[1]; LPInput[1]=LPInput[2];
                                      HPInput[0] = HPInput[1]; HPInput[1]=HPInput[2];
                                      LPOutput[0] = LPOutput[1]; LPOutput[1]=LPOutput[2];
                                      HPOutput[0] = HPOutput[1]; HPOutput[1]=HPOutput[2]; 
                                      
                                      for (j=0; j<kOverlapCount-2; j++)
                                      OverlapBuf[j][k] = OverlapBuf[j+1][k] + cBuf[k+kWindowSize*(j+1)/kOverlapCount];
			                          OverlapBuf[kOverlapCount-2][k] = cBuf[k+kWindowSize*(kOverlapCount-1)/kOverlapCount];
		                     }

                            for (k = 0; k < kWindowSize*(kOverlapCount-1)/kOverlapCount; k++) {
                                lBuf[k] = lBuf[k+kWindowSize/kOverlapCount];
                                rBuf[k] = rBuf[k+kWindowSize/kOverlapCount];                                
                                }                                
        tmpDataSize += kWindowSize/kOverlapCount*infileChunk.BlockAlign;
        FilePos = ftell(outfile);
        tmpRIFFSize = tmpDataSize + 36;
        fseek(outfile,4,0); fwrite(&tmpRIFFSize,sizeof(tmpRIFFSize),1,outfile);
        fseek(outfile,40,0); fwrite(&tmpDataSize,sizeof(tmpDataSize),1,outfile);        
        fseek(outfile,FilePos,0);
        }
    }
    for (i = Count; i<= kWindowSize; i++) lBuf[i-1] = rBuf[i-1] = 0.0;
    for (i = 1; i < kOverlapCount; i++) {
                             for (k = 0; k < kWindowSize/2; k++) {
                                 lFFTBuf[2*k] = lBuf[k]*CosWindow[k];
                                 lFFTBuf[2*k+1] = lBuf[kWindowSize-k-1]*CosWindow[kWindowSize-k-1];
                                 rFFTBuf[2*k] = rBuf[k]*CosWindow[k];
                                 rFFTBuf[2*k+1] = rBuf[kWindowSize-k-1]*CosWindow[kWindowSize-k-1];
                                 }              
                             DoFFT(lFFTBuf, kWindowSize/2, -1);
                             DoFFT(rFFTBuf, kWindowSize/2, -1);
                             CenterCut(lFFTBuf, rFFTBuf, cFFTBuf, kWindowSize/2);                             
                             DoFFT(cFFTBuf, kWindowSize/2, 1);

			                 for (k=0; k < kWindowSize/2; k++) {
				                      cBuf[k] = cFFTBuf[2*k]; cBuf[kWindowSize-k-1] = cFFTBuf[2*k+1]; }                          
		                     for (k=0; k < kWindowSize/kOverlapCount; k++) {                                      

                                      LPInput[2] = cBuf[k] + OverlapBuf[0][k];
                                      HPInput[2] = cBuf[k] + OverlapBuf[0][k];
                                      LPOutput[2] = (b1[0]/a1[0])*LPInput[2]+(b1[1]/a1[0])*LPInput[1]+(b1[2]/a1[0])*LPInput[0]-(a1[1]/a1[0])*LPOutput[1]-(a1[2]/a1[0])*LPOutput[0];
                                      HPOutput[2] = (b2[0]/a2[0])*HPInput[2]+(b2[1]/a2[0])*HPInput[1]+(b2[2]/a2[0])*HPInput[0]-(a2[1]/a2[0])*HPOutput[1]-(a2[2]/a2[0])*HPOutput[0];                                      
                                      NewLeft = lBuf[k] - cBuf[k] - OverlapBuf[0][k] + LPOutput[2] + HPOutput[2];
                                      NewRight = rBuf[k] - cBuf[k] - OverlapBuf[0][k] + LPOutput[2] + HPOutput[2];
                                      
                                      BlockWrite(NewLeft, infileChunk.Bits, infileChunk.FormatTag, outfile);
                                      BlockWrite(NewRight,infileChunk.Bits, infileChunk.FormatTag, outfile);    

                                      LPInput[0] = LPInput[1]; LPInput[1]=LPInput[2];
                                      HPInput[0] = HPInput[1]; HPInput[1]=HPInput[2];
                                      LPOutput[0] = LPOutput[1]; LPOutput[1]=LPOutput[2];
                                      HPOutput[0] = HPOutput[1]; HPOutput[1]=HPOutput[2]; 
                                      
                                      for (j=0; j<kOverlapCount-2; j++)
                                      OverlapBuf[j][k] = OverlapBuf[j+1][k] + cBuf[k+kWindowSize*(j+1)/kOverlapCount];
			                          OverlapBuf[kOverlapCount-2][k] = cBuf[k+kWindowSize*(kOverlapCount-1)/kOverlapCount];
		                     }

                            for (k = 0; k < kWindowSize*(kOverlapCount-1)/kOverlapCount; k++) {
                                lBuf[k] = lBuf[k+kWindowSize/kOverlapCount];
                                rBuf[k] = rBuf[k+kWindowSize/kOverlapCount];                                
                                }                                
                            for (k = kWindowSize*(kOverlapCount-1)/kOverlapCount; k < kWindowSize; k++)
                                lBuf[k] = rBuf[k] = 0.0;
                            tmpDataSize += kWindowSize/kOverlapCount*infileChunk.BlockAlign;
                            FilePos = ftell(outfile);
                            tmpRIFFSize = tmpDataSize + 36;
                            fseek(outfile,4,0); fwrite(&tmpRIFFSize,sizeof(tmpRIFFSize),1,outfile);
                            fseek(outfile,40,0); fwrite(&tmpDataSize,sizeof(tmpDataSize),1,outfile);        
                            fseek(outfile,FilePos,0);
    }
    if (Count != kWindowSize*(kOverlapCount-1)/kOverlapCount) {
                             for (k = 0; k < kWindowSize/2; k++) {
                                 lFFTBuf[2*k] = lBuf[k]*CosWindow[k];
                                 lFFTBuf[2*k+1] = lBuf[kWindowSize-k-1]*CosWindow[kWindowSize-k-1];
                                 rFFTBuf[2*k] = rBuf[k]*CosWindow[k];
                                 rFFTBuf[2*k+1] = rBuf[kWindowSize-k-1]*CosWindow[kWindowSize-k-1];
                                 }              
                             DoFFT(lFFTBuf, kWindowSize/2, -1);
                             DoFFT(rFFTBuf, kWindowSize/2, -1);
                             CenterCut(lFFTBuf, rFFTBuf, cFFTBuf, kWindowSize/2);                             
                             DoFFT(cFFTBuf, kWindowSize/2, 1);

			                 for (k=0; k < kWindowSize/2; k++) {
				                      cBuf[k] = cFFTBuf[2*k]; cBuf[kWindowSize-k-1] = cFFTBuf[2*k+1]; }                          
		                     for (k=0; k < Count - kWindowSize*(kOverlapCount-1)/kOverlapCount; k++) {                                      

                                      LPInput[2] = cBuf[k] + OverlapBuf[0][k];
                                      HPInput[2] = cBuf[k] + OverlapBuf[0][k];
                                      LPOutput[2] = (b1[0]/a1[0])*LPInput[2]+(b1[1]/a1[0])*LPInput[1]+(b1[2]/a1[0])*LPInput[0]-(a1[1]/a1[0])*LPOutput[1]-(a1[2]/a1[0])*LPOutput[0];
                                      HPOutput[2] = (b2[0]/a2[0])*HPInput[2]+(b2[1]/a2[0])*HPInput[1]+(b2[2]/a2[0])*HPInput[0]-(a2[1]/a2[0])*HPOutput[1]-(a2[2]/a2[0])*HPOutput[0];                                      
                                      NewLeft = lBuf[k] - cBuf[k] - OverlapBuf[0][k] + LPOutput[2] + HPOutput[2];
                                      NewRight = rBuf[k] - cBuf[k] - OverlapBuf[0][k] + LPOutput[2] + HPOutput[2];
                                      
                                      BlockWrite(NewLeft, infileChunk.Bits, infileChunk.FormatTag, outfile);
                                      BlockWrite(NewRight,infileChunk.Bits, infileChunk.FormatTag, outfile);    

                                      LPInput[0] = LPInput[1]; LPInput[1]=LPInput[2];
                                      HPInput[0] = HPInput[1]; HPInput[1]=HPInput[2];
                                      LPOutput[0] = LPOutput[1]; LPOutput[1]=LPOutput[2];
                                      HPOutput[0] = HPOutput[1]; HPOutput[1]=HPOutput[2]; 
                                                             			      
                                      for (j=0; j<kOverlapCount-2; j++)
                                      OverlapBuf[j][k] = OverlapBuf[j+1][k] + cBuf[k+kWindowSize*(j+1)/kOverlapCount];
			                          OverlapBuf[kOverlapCount-2][k] = cBuf[k+kWindowSize*(kOverlapCount-1)/kOverlapCount];
		                     }
                              
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
