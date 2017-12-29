#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include<time.h>
#include<sys/stat.h>

__constant__ int d_NumL[2];
__constant__ char d_Pchar[1084];
__constant__ double parameter[5916];
__constant__ float d_stab[4096];
__constant__ float d_deltah[16];
__constant__ float d_deltas[16];

char str2int_CPU(char c)
{
        switch (c)
        {
                case 'A':
                        return 0;
                case 'C':
                        return 1;
                case 'G':
                        return 2;              
                case 'T':  
                        return 3;       
        }
        return 4;
}

__device__ void str2int(char c,char *d_numSeq,int id)
{
        switch (c)
        {
                case 'A':
                        d_numSeq[id]=0;
                        break;
                case 'C':
                        d_numSeq[id]=1;
                        break;
                case 'G':
                        d_numSeq[id]=2;
                        break;
                case 'T':
                        d_numSeq[id]=3;
                        break;
                default:
                        d_numSeq[id]=4;
                        break;
        }
}

__device__ void str2int_rev(char c,char *d_numSeq,int id)
{
        switch (c)
        {
                case 'T':
                        d_numSeq[id]=0;
                        break;
                case 'G':
                        d_numSeq[id]=1;
                        break;
                case 'C':
                        d_numSeq[id]=2;
                        break;                 
                case 'A':               
                        d_numSeq[id]=3;
                        break;
                default:
                        d_numSeq[id]=4;
                        break;
        }
}

void readLoop(FILE *file,double *v1,double *v2,double *v3)
{
        char *line,*p,*q;
        
        line=(char *)malloc(200);
        memset(line,'\0',200);
        fgets(line,200,file);

        p = line;
        while (*p==' '||*p=='\t')
                p++;
        while (*p=='0'||*p=='1'||*p=='2'||*p=='3'||*p=='4'||*p=='5'||*p=='6'||*p=='7'||*p=='8'||*p=='9') 
                p++;
        while (*p==' '||*p=='\t') 
                p++;

        q = p;
        while (!(*q==' '||*q=='\t')) 
                q++;
        *q = '\0';
        q++;
        if (!strcmp(p, "inf"))
                *v1 =1.0*INFINITY;
        else 
                sscanf(p, "%lf", v1);
        while (*q==' '||*q=='\t')
                q++;

        p = q;
        while (!(*p==' '||*p=='\t'))
                p++;
        *p = '\0';
        p++;
        if (!strcmp(q, "inf"))
                *v2 =1.0*INFINITY;
        else 
                sscanf(q, "%lf", v2);
        while (*p==' '||*p=='\t')
                p++;

        q = p;
        while (!(*q==' '||*q=='\t') && (*q != '\0'))
                q++;
        *q = '\0';
        if (!strcmp(p, "inf"))
                *v3 =1.0*INFINITY;
        else 
                sscanf(p, "%lf", v3);
}

void getStack(char *path,double *parameter)
{
        int i, j, ii, jj;
        FILE *sFile, *hFile;
        char *line;

        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stack.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stack.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i = 0; i < 5; ++i)
        {
                for (ii = 0; ii < 5; ++ii)
                {
                        for (j = 0; j < 5; ++j)
                        {
                                for (jj = 0; jj < 5; ++jj)
                                {
                                        if (i == 4 || j == 4 || ii == 4 || jj == 4) //N 
                                        {
                                                parameter[i*125+ii*25+j*5+jj] = -1.0;
                                                parameter[625+i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                        }
                                        else 
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        parameter[i*125+ii*25+j*5+jj] = atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[625+i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        parameter[625+i*125+ii*25+j*5+jj] = atof(line);

                                                if (fabs(parameter[i*125+ii*25+j*5+jj])>999999999 ||fabs(parameter[625+i*125+ii*25+j*5+jj])>999999999) 
                                                {
                                                        parameter[i*125+ii*25+j*5+jj] = -1.0;
                                                        parameter[625+i*125+ii*25+j*5+jj] =1.0*INFINITY;
                                                }
                                        }
                                }
                        }
                }
        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getStackint2(char *path,double *parameter)
{
        int i, j, ii, jj;
        FILE *sFile, *hFile;
        char *line;

        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stackmm.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stackmm.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i = 0; i < 5; ++i)
        {
                for (ii = 0; ii < 5; ++ii)
                {
                        for (j = 0; j < 5; ++j)
                        {
                                for (jj = 0; jj < 5; ++jj)
                                {
                                        if (i == 4 || j == 4 || ii == 4 || jj == 4)
                                        {
                                                parameter[1250+i*125+ii*25+j*5+jj] = -1.0;
                                                parameter[1875+i*125+ii*25+j*5+jj] =1.0*INFINITY;
                                        } 
                                        else 
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStackint2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[1250+i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        parameter[1250+i*125+ii*25+j*5+jj] = atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStackint2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[1875+i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        parameter[1875+i*125+ii*25+j*5+jj] = atof(line);

                                                if(fabs(parameter[1250+i*125+ii*25+j*5+jj])>999999999||fabs(parameter[1875+i*125+ii*25+j*5+jj])>999999999)
                                                {
                                                        parameter[1250+i*125+ii*25+j*5+jj] = -1.0;
                                                        parameter[1875+i*125+ii*25+j*5+jj] =1.0*INFINITY;
                                                }
                                        }
                                }
                        }
                }
        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getDangle(char *path,double *parameter)
{
        int i, j, k;
        FILE *sFile, *hFile;
        char *line;
        
        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"dangle.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"dangle.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        for (k = 0; k < 5; ++k) 
                        {
                                if (i == 4 || j == 4) 
                                {
                                        parameter[2500+i*25+k*5+j] = -1.0;
                                        parameter[2625+i*25+k*5+j] =1.0*INFINITY;
                                }
                                else if (k == 4)
                                {
                                        parameter[2500+i*25+k*5+j] = -1.0;
                                        parameter[2625+i*25+k*5+j] =1.0*INFINITY;
                                } 
                                else
                                {
                                        if(fgets(line,20,sFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");
                                                exit(1);
                                        }
                                        if(strncmp(line, "inf", 3)==0)
                                                parameter[2500+i*25+k*5+j]=1.0*INFINITY;
                                        else
                                                parameter[2500+i*25+k*5+j]=atof(line);

                                        if(fgets(line,20,hFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");        
                                                exit(1);        
                                        }
                                        if(strncmp(line, "inf", 3)==0)        
                                                parameter[2625+i*25+k*5+j]=1.0*INFINITY;           
                                        else        
                                                parameter[2625+i*25+k*5+j]=atof(line);

                                        if(fabs(parameter[2500+i*25+k*5+j])>999999999||fabs(parameter[2625+i*25+k*5+j])>999999999) 
                                        {
                                                parameter[2500+i*25+k*5+j] = -1.0;
                                                parameter[2625+i*25+k*5+j] =1.0*INFINITY;
                                        }
                                }
                        }

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        for (k = 0; k < 5; ++k) 
                        {
                                if (i == 4 || j == 4)
                                {
                                        parameter[2750+i*25+j*5+k] = -1.0;
                                        parameter[2875+i*25+j*5+k] =1.0*INFINITY;
                                } 
                                else if (k == 4) 
                                {
                                        parameter[2750+i*25+j*5+k] = -1.0;
                                        parameter[2875+i*25+j*5+k] =1.0*INFINITY;
                                }
                                else
                                {
                                        if(fgets(line,20,sFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");
                                                exit(1);
                                        }
                                        if(strncmp(line, "inf", 3)==0)
                                                parameter[2750+i*25+j*5+k]=1.0*INFINITY;
                                        else
                                                parameter[2750+i*25+j*5+k]=atof(line);

                                        if(fgets(line,20,hFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");        
                                                exit(1);        
                                        }
                                        if(strncmp(line, "inf", 3)==0)        
                                                parameter[2875+i*25+j*5+k]=1.0*INFINITY;           
                                        else        
                                                parameter[2875+i*25+j*5+k]=atof(line);

                                        if(fabs(parameter[2750+i*25+j*5+k])>999999999||fabs(parameter[2875+i*25+j*5+k])>999999999)
                                        {
                                                parameter[2750+i*25+j*5+k] = -1.0;
                                                parameter[2875+i*25+j*5+k] =1.0*INFINITY;
                                        }
                                }
                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getLoop(char *path,double *parameter)
{
        int k;
        FILE *sFile, *hFile;
        char *line;

        k=strlen(path)+20;
        line=(char *)malloc(k);
        memset(line,'\0',k);
        strcpy(line,path);
        strcat(line,"loops.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',k);
        strcpy(line,path);
        strcat(line,"loops.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        for (k = 0; k < 30; ++k)
        {
                readLoop(sFile, &parameter[3030+k], &parameter[3060+k], &parameter[3000+k]);
                readLoop(hFile, &parameter[3120+k], &parameter[3150+k], &parameter[3090+k]);
        }
        fclose(sFile);
        fclose(hFile);
}

void getTstack(char *path,double *parameter)
{
        int i1, j1, i2, j2;
        FILE *sFile, *hFile;
        char *line;

        i1=strlen(path)+20;
        line=(char *)malloc(i1);
        memset(line,'\0',i1);
        strcpy(line,path);
        strcat(line,"tstack_tm_inf.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i1);
        strcpy(line,path);      
        strcat(line,"tstack.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }             
        hFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);   
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i1 = 0; i1 < 5; ++i1)
                for (i2 = 0; i2 < 5; ++i2)
                        for (j1 = 0; j1 < 5; ++j1)
                                for (j2 = 0; j2 < 5; ++j2)
                                        if (i1 == 4 || j1 == 4)
                                        {
                                                parameter[3805+i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                parameter[3180+i1*125+i2*25+j1*5+j2] = -1.0;
                                        }
                                        else if (i2 == 4 || j2 == 4)
                                        {
                                                parameter[3180+i1*125+i2*25+j1*5+j2] = 0.00000000001;
                                                parameter[3805+i1*125+i2*25+j1*5+j2] = 0.0;
                                        }
                                        else
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[3180+i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        parameter[3180+i1*125+i2*25+j1*5+j2]=atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[3805+i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        parameter[3805+i1*125+i2*25+j1*5+j2]=atof(line);

                                                if (fabs(parameter[3180+i1*125+i2*25+j1*5+j2])>999999999||fabs(parameter[3805+i1*125+i2*25+j1*5+j2])>999999999)
                                                {
                                                        parameter[3180+i1*125+i2*25+j1*5+j2] = -1.0;
                                                        parameter[3805+i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                }
                                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getTstack2(char *path,double *parameter)
{
        int i1, j1, i2, j2;
        FILE *sFile, *hFile;
        char *line;

        i1=strlen(path)+20;
        line=(char *)malloc(i1);
        memset(line,'\0',i1);
        strcpy(line,path);
        strcat(line,"tstack2.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i1);
        strcpy(line,path);      
        strcat(line,"tstack2.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }             
        hFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);   
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i1 = 0; i1 < 5; ++i1)
                for (i2 = 0; i2 < 5; ++i2)
                        for (j1 = 0; j1 < 5; ++j1)
                                for (j2 = 0; j2 < 5; ++j2)
                                        if (i1 == 4 || j1 == 4)
                                        {
                                                parameter[5055+i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                parameter[4430+i1*125+i2*25+j1*5+j2] = -1.0;
                                        }
                                        else if (i2 == 4 || j2 == 4)
                                        {
                                                parameter[4430+i1*125+i2*25+j1*5+j2] = 0.00000000001;
                                                parameter[5055+i1*125+i2*25+j1*5+j2] = 0.0;
                                        }
                                        else
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[4430+i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        parameter[4430+i1*125+i2*25+j1*5+j2]=atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        parameter[5055+i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        parameter[5055+i1*125+i2*25+j1*5+j2]=atof(line);


                                                if (fabs(parameter[4430+i1*125+i2*25+j1*5+j2])>999999999||fabs(parameter[5055+i1*125+i2*25+j1*5+j2])>999999999)
                                                {
                                                        parameter[4430+i1*125+i2*25+j1*5+j2] = -1.0;
                                                        parameter[5055+i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                }
                                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

int get_num_line(char *path,int flag)
{
	FILE *fp;
	int i,size;
	char *line;

	i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
	if(flag==0)
	        strcat(line,"triloop.ds");
	else
		strcat(line,"tetraloop.ds");

        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        fp=fopen(line,"r");
        if(fp==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

	size=0;
	while(fgets(line,i,fp)!=NULL)
		size++;
	return size;
}

void getTriloop(char *path,double *parameter,char *Pchar,int NumL[])
{
        FILE *sFile, *hFile;
        int i,turn;
        char *line,seq[10],value[10];
        
        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"triloop.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
	
	turn=0;
        while(fscanf(sFile,"%s\t%s\n",seq,value)!=EOF)
        {
		for (i=0;i<5;i++)
			Pchar[5*turn+i]=str2int_CPU(seq[i]);
		if(value[0]=='i')
			parameter[5730+turn]=1.0*INFINITY;
		else
			parameter[5730+turn]=atof(value);
		turn++;
        }
        fclose(sFile);

	i=strlen(path)+20;
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"triloop.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

	turn=0;
        while(fscanf(hFile,"%s\t%s\n",seq,value)!=EOF)
        {
		for(i=0;i<5;i++)
			Pchar[5*NumL[0]+turn*5+i]=str2int_CPU(seq[i]);
		if(value[0]=='i')
			parameter[5730+NumL[0]+turn]=1.0*INFINITY;
		else
			parameter[5730+NumL[0]+turn]=atof(value);
		turn++;
        }
        fclose(hFile);
}

void getTetraloop(char *path,double *parameter,char *Pchar,int NumL[])
{
        FILE *sFile, *hFile;
        int i, turn;
        char *line,seq[10],value[10];

        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"tetraloop.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

	turn=0;
        while(fscanf(sFile,"%s\t%s\n",seq,value)!=EOF)
        {
		for(i=0;i<6;i++)
			Pchar[10*NumL[0]+turn*6+i]=str2int_CPU(seq[i]);
		if(value[0]=='i')
			parameter[5730+2*NumL[0]+turn]=1.0*INFINITY;
		else
			parameter[5730+2*NumL[0]+turn]=atof(value);
		turn++;
        }
        fclose(sFile);

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"tetraloop.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);
        
	turn=0;
        while(fscanf(hFile,"%s\t%s\n",seq,value)!=EOF)
        {
		for(i=0;i<6;i++)
			Pchar[10*NumL[0]+6*NumL[1]+6*turn+i]=str2int_CPU(seq[i]);
		if(value[0]=='i')
			parameter[5730+2*NumL[0]+NumL[1]+turn]=1.0*INFINITY;
		else
			parameter[5730+2*NumL[0]+NumL[1]+turn]=atof(value);
		turn++;
        }
        fclose(hFile);
}

void tableStartATS(double atp_value,double parameter[] )
{
        int i, j;

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        parameter[5680+i*5+j] = 0.00000000001;
        parameter[5680+3] = parameter[5680+15] = atp_value;
}

void tableStartATH(double atp_value,double parameter[])
{
        int i, j;

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        parameter[5705+i*5+j] = 0.0;
        parameter[5705+3] = parameter[5705+15] = atp_value;
}

//end read parameter
__device__ void initMatrix2(int length,double *d_DPT,char *d_numSeq,int id)
{
	int i,j;
	for(i=1;i<=length;++i)
		for(j=i;j<=length;++j)
			if(j-i<4 || (d_numSeq[id*54+i]+d_numSeq[id*54+j]!=3))
			{
				d_DPT[id*1340+(i-1)*(length-1)+j-1]=1.0*INFINITY;
				d_DPT[id*1340+625+(i-1)*(length-1)+j-1]=-1.0;
			}
			else
			{
				d_DPT[id*1340+(i-1)*(length-1)+j-1]=0.0;
				d_DPT[id*1340+625+(i-1)*(length-1)+j-1]=-3224.0;
			}
}

__device__ double Ss(int i,int j,int k,int length,char *d_numSeq,int id)
{
	if(k==2)
	{
		if(i>=j)
			return -1.0;
		if(i==length||j==length+1)
			return -1.0;

		if(i>length)
			i-=length;
		if(j>length)
			j-=length;
		return parameter[d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j-1]];
	}
	else
		return parameter[d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
}

__device__ double Hs(int i,int j,int k,int length,char *d_numSeq,int id)
{
	if(k==2)
	{
		if(i>= j)
			return 1.0*INFINITY;
		if(i==length||j==length+1)
			return 1.0*INFINITY;

		if(i>length)
			i-=length;
		if(j>length)
			j-=length;
		if(fabs(parameter[625+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j-1]])<999999999)
			return parameter[625+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j-1]];
		else
			return 1.0*INFINITY;
	}
	else
		return parameter[625+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
}

__device__ void maxTM2(int i,int j,int length,double *d_DPT,char *d_numSeq,int id)
{
	double T0,T1,S0,S1,H0,H1;

	S0=d_DPT[id*1340+625+(i-1)*(length-1)+j-1];
	H0=d_DPT[id*1340+(i-1)*(length-1)+j-1];
	T0=(H0+d_DPT[id*1340+1302])/(S0+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
	if(fabs(d_DPT[id*1340+(i-1)*(length-1)+j-1])<999999999)
	{
		S1=(d_DPT[id*1340+625+i*(length-1)+j-2]+Ss(i,j,2,length,d_numSeq,id));
		H1=(d_DPT[id*1340+i*(length-1)+j-2]+Hs(i,j,2,length,d_numSeq,id));
	}
	else
	{
		S1=-1.0;
		H1=1.0*INFINITY;
	}
	T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
	if(S1<-2500.0)
	{
		S1=-3224.0;
		H1=0.0;
	}
	if(S0<-2500.0)
	{
		S0=-3224.0;
		H0=0.0;
 	}

	if(T1>T0)
	{
		d_DPT[id*1340+625+(i-1)*(length-1)+j-1]=S1;
		d_DPT[id*1340+(i-1)*(length-1)+j-1]= H1;
	}
	else
	{
		d_DPT[id*1340+625+(i-1)*(length-1)+j-1]=S0;
		d_DPT[id*1340+(i-1)*(length-1)+j-1]=H0;
	}
}

__device__ void calc_bulge_internal2(int i,int j,int ii,int jj,int pos,int traceback,int length,double *d_DPT,char *d_numSeq,int id)
{
	int loopSize1,loopSize2,loopSize;
	double T1,T2,S,H;

	S=-3224.0;
	H=0.0;
	loopSize1=ii-i-1;
	loopSize2=j-jj-1;
	if(loopSize1+loopSize2>30)
	{
		d_DPT[id*1340+pos]=-1.0;
		d_DPT[id*1340+pos+1]=1.0*INFINITY;
		return;
	}

	loopSize=loopSize1+loopSize2-1;
	if((loopSize1==0&&loopSize2>0)||(loopSize2==0&&loopSize1>0))
	{
		if(loopSize2==1||loopSize1==1)
		{ 
			if((loopSize2==1&&loopSize1==0)||(loopSize2==0&&loopSize1==1))
			{
				H=parameter[3150+loopSize]+parameter[625+d_numSeq[id*54+i]*125+d_numSeq[id*54+ii]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+jj]];
				S=parameter[3060+loopSize]+parameter[d_numSeq[id*54+i]*125+d_numSeq[id*54+ii]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+jj]];
 			}
			if(traceback!=1)
			{
				H+=d_DPT[id*1340+(ii-1)*(length-1)+jj-1];
				S+=d_DPT[id*1340+625+(ii-1)*(length-1)+jj-1];
			}

			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
			T2=(d_DPT[id*1340+(i-1)*(length-1)+j-1]+d_DPT[id*1340+1302])/((d_DPT[id*1340+625+(i-1)*(length-1)+j-1])+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if((T1>T2)||((traceback&&T1>=T2)||traceback==1))
			{
				d_DPT[id*1340+pos]=S;
				d_DPT[id*1340+pos+1]=H;
			}
		}
		else
		{
			H=parameter[3150+loopSize]+parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[5705+d_numSeq[id*54+ii]*5+d_numSeq[id*54+27+jj]];
			if(traceback!=1)
				H+=d_DPT[id*1340+(ii-1)*(length-1)+jj-1];

			S=parameter[3060+loopSize]+parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[5680+d_numSeq[id*54+ii]*5+d_numSeq[id*54+27+jj]];
			if(traceback!=1)
				S+=d_DPT[id*1340+625+(ii-1)*(length-1)+jj-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
			T2=(d_DPT[id*1340+(i-1)*(length-1)+j-1]+d_DPT[id*1340+1302])/(d_DPT[id*1340+625+(i-1)*(length-1)+j-1]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
			{
				d_DPT[id*1340+pos]=S;
				d_DPT[id*1340+pos+1]=H;
			}
		}
	}
	else if(loopSize1==1&&loopSize2==1)
	{
		S=parameter[1250+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j-1]]+parameter[1250+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj+1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]];
		if(traceback!=1)
			S+=d_DPT[id*1340+625+(ii-1)*(length-1)+jj-1];

		H=parameter[1875+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j-1]]+parameter[1875+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj+1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]];
		if(traceback!=1)
			H+=d_DPT[id*1340+(ii-1)*(length-1)+jj-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
		T2=(d_DPT[id*1340+(i-1)*(length-1)+j-1]+d_DPT[id*1340+1302])/(d_DPT[id*1340+625+(i-1)*(length-1)+j-1]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if((T1-T2>=0.000001)||traceback)
		{
			if((T1>T2)||((traceback&&T1>= T2)||traceback==1))
			{
				d_DPT[id*1340+pos]=S;
				d_DPT[id*1340+pos+1]=H;
			}
		}
		return;
	}
	else
	{
		H=parameter[3120+loopSize]+parameter[3805+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j-1]]+parameter[3805+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj+1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]];
		if(traceback!=1)
			H+=d_DPT[id*1340+(ii-1)*(length-1)+jj-1];

		S=parameter[3030+loopSize]+parameter[3180+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j-1]]+parameter[3180+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj+1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
		if(traceback!=1)
			S+=d_DPT[id*1340+625+(ii-1)*(length-1)+jj-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}

		T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
		T2=(d_DPT[id*1340+(i-1)*(length-1)+j-1]+d_DPT[id*1340+1302])/((d_DPT[id*1340+625+(i-1)*(length-1)+j-1])+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
		{
			d_DPT[id*1340+pos]=S;
			d_DPT[id*1340+pos+1]=H;
		}
	}
	return;
}

__device__ void CBI(int i,int j,int pos,int traceback,int length,double *d_DPT,char *d_numSeq,int id)
{
	int d,ii,jj;

	for(d=j-i-3;d>=4&&d>=j-i-32;--d)
		for(ii=i+1;ii<j-d&&ii<=length;++ii)
		{
			jj=d+ii;
			if(traceback==0)
			{
				d_DPT[id*1340+pos]=-1.0;
				d_DPT[id*1340+pos+1]=1.0*INFINITY;
			}
			if(fabs(d_DPT[id*1340+(ii-1)*(length-1)+jj-1])<999999999)
			{
				calc_bulge_internal2(i,j,ii,jj,pos,traceback,length,d_DPT,d_numSeq,id);
				if(fabs(d_DPT[id*1340+pos+1])<999999999)
				{
					if(d_DPT[id*1340+pos] <-2500.0)
					{
						d_DPT[id*1340+pos+1]=-3224.0;
						d_DPT[id*1340+pos+1]=0.0;
					}
					if(traceback==0)
					{
						d_DPT[id*1340+(i-1)*(length-1)+j-1]=d_DPT[id*1340+pos+1];
						d_DPT[id*1340+625+(i-1)*(length-1)+j-1]=d_DPT[id*1340+pos];
					}
				}
			}
		}
	return;
}

__device__ int find_pos(char *ref,int ref_start,int start,int length,int num)
{
	int flag,i,j;

	for(i=0;i<num;i++)
	{
		flag=0;
		for(j=0;j<length;j++)
		{
			if(ref[ref_start+j]!=d_Pchar[start+i*length+j])
			{
				flag++;
				break;
			}
		}
		if(flag==0)
			return i;
	}
	return -1;
}

__device__ void calc_hairpin(int i,int j,int pos_start,int traceback,int length,double *d_DPT,char *d_numSeq,int id)
{
	int pos,loopSize=j-i-1;
	double T1,T2;
	
	if(loopSize < 3)
	{
		d_DPT[id*1340+pos_start]=-1.0;
		d_DPT[id*1340+pos_start+1]=1.0*INFINITY;
		return;
	}
	if(i<=length&&length<j)
	{
		d_DPT[id*1340+pos_start]=-1.0;
		d_DPT[id*1340+pos_start+1]=1.0*INFINITY;
		return;
	}
	else if(i>length)
	{
		i-= length;
		j-= length;
	}
	if(loopSize<=30)
	{
		d_DPT[id*1340+pos_start+1]=parameter[3090+loopSize-1];
		d_DPT[id*1340+pos_start]=parameter[3000+loopSize-1];
	}
	else
	{
		d_DPT[id*1340+pos_start+1]=parameter[3090+29];
		d_DPT[id*1340+pos_start]=parameter[3000+29];
	}

	if(loopSize>3) // for loops 4 bp and more in length, terminal mm are accounted
	{
		d_DPT[id*1340+pos_start+1]+=parameter[5055+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+j]*5+d_numSeq[id*54+j-1]];
		d_DPT[id*1340+pos_start]+=parameter[4430+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+j]*5+d_numSeq[id*54+j-1]];
	}
	else if(loopSize == 3) // for loops 3 bp in length at-penalty is considered
	{
		d_DPT[id*1340+pos_start+1]+=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+j]];
		d_DPT[id*1340+pos_start]+=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+j]];
	}

	if(loopSize==3) // closing AT-penalty (+), triloop bonus, hairpin of 3 (+) 
	{
		pos=find_pos(d_numSeq,(id*54+i),5*d_NumL[0],5,d_NumL[0]);
		if(pos!=-1)
			d_DPT[id*1340+pos_start+1]+=parameter[5730+d_NumL[0]+pos];

		pos=find_pos(d_numSeq,(id*54+i),0,5,d_NumL[0]);
		if(pos!=-1)
			d_DPT[id*1340+pos_start]+=parameter[5730+pos];
	}
	else if (loopSize == 4) // terminal mismatch, tetraloop bonus, hairpin of 4
	{
		pos=find_pos(d_numSeq,(id*54+i),10*d_NumL[0]+6*d_NumL[1],6,d_NumL[1]);
		if(pos!=-1)
			d_DPT[id*1340+pos_start+1]+=parameter[5730+2*d_NumL[0]+d_NumL[1]+pos];

		pos=find_pos(d_numSeq,(id*54+i),10*d_NumL[0],6,d_NumL[1]);
		if(pos!=-1)
			d_DPT[id*1340+pos_start]+=parameter[5730+2*d_NumL[0]+pos];
	}
	if(fabs(d_DPT[id*1340+pos_start+1])>999999999)
	{
		d_DPT[id*1340+pos_start+1] =1.0*INFINITY;
		d_DPT[id*1340+pos_start] = -1.0;
	}
	T1 = (d_DPT[id*1340+pos_start+1] +d_DPT[id*1340+1302]) / ((d_DPT[id*1340+pos_start] +d_DPT[id*1340+1303]+ d_DPT[id*1340+1304]));
	T2 = (d_DPT[id*1340+(i-1)*(length-1)+j-1] +d_DPT[id*1340+1302]) / ((d_DPT[id*1340+625+(i-1)*(length-1)+j-1]) +d_DPT[id*1340+1303]+ d_DPT[id*1340+1304]);
	if(T1 < T2 && traceback == 0)
	{
		d_DPT[id*1340+pos_start] =d_DPT[id*1340+625+(i-1)*(length-1)+j-1];
		d_DPT[id*1340+pos_start+1] =d_DPT[id*1340+(i-1)*(length-1)+j-1];
	}
	return;
}

__device__ void fillMatrix2(int length,double *d_DPT,char *d_numSeq,int id)
{
	int i, j;

	for (j = 2; j <= length; ++j)
		for (i = j - 3 - 1; i >= 1; --i)
		{
			if (fabs(d_DPT[id*1340+(i-1)*(length-1)+j-1])<999999999)
			{
				d_DPT[id*1340+1310] = -1.0;
				d_DPT[id*1340+1311] =1.0*INFINITY;
				maxTM2(i,j,length,d_DPT,d_numSeq,id);
				CBI(i,j,1310,0,length,d_DPT,d_numSeq,id);

				d_DPT[id*1340+1310] = -1.0;
				d_DPT[id*1340+1311]=1.0*INFINITY;
				calc_hairpin(i,j,1310,0,length,d_DPT,d_numSeq,id);
				if(fabs(d_DPT[id*1340+1311])<999999999)
				{
					if(d_DPT[id*1340+1310]<-2500.0) /* to not give dH any value if dS is unreasonable */
					{
						d_DPT[id*1340+1310]=-3224.0;
						d_DPT[id*1340+1311]= 0.0;
					}
					d_DPT[id*1340+625+(i-1)*(length-1)+j-1]=d_DPT[id*1340+1310];
					d_DPT[id*1340+(i-1)*(length-1)+j-1]=d_DPT[id*1340+1311];
				}
			}
		}
}

__device__ int max5(double a,double b,double c,double d,double e)
{
	if(a>b&&a>c&&a>d&&a>e)
		return 1;
	else if(b>c&&b>d&&b>e)
		return 2;
	else if(c>d&&c>e)
		return 3;
	else if(d>e)
		return 4;
	else
		return 5;
}

__device__ double Sd5(int i,int j,char *d_numSeq,int id)
{
	return parameter[2750+d_numSeq[id*54+i]*25+d_numSeq[id*54+j]*5+d_numSeq[id*54+j-1]];
}

__device__ double Hd5(int i,int j,char *d_numSeq,int id)
{
	return parameter[2875+d_numSeq[id*54+i]*25+d_numSeq[id*54+j]*5+d_numSeq[id*54+j-1]];
}

__device__ double Sd3(int i,int j,char *d_numSeq,int id)
{
	return parameter[2500+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+j]];
}

__device__ double Hd3(int i,int j,char *d_numSeq,int id)
{
	return parameter[2625+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+j]];
}

__device__ double Ststack(int i,int j,char *d_numSeq,int id)
{
	return parameter[4430+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+j]*5+d_numSeq[id*54+j-1]];
}

__device__ double Htstack(int i,int j,char *d_numSeq,int id)
{
	return parameter[5055+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+j]*5+d_numSeq[id*54+j-1]];
}

__device__ double END5_1(int i,int hs,int length,double *d_DPT,char *d_numSeq,int id)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	max_tm=-1.0*INFINITY;
	H_max=1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-5;++k)
	{
		T1=(d_DPT[id*1340+1276+k]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1250+k]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		T2=d_DPT[id*1340+1302]/(d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(T1>=T2)
		{
			H=d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+k*(length-1)+i-1];
			S=d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+625+k*(length-1)+i-1];
			if(fabs(H)>999999999||H>0||S>0)  // H and S must be greater than 0 to avoid BS
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		}
		else
		{
			H=parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+k*(length-1)+i-1];
			S=parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+625+k*(length-1)+i-1];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		}

		if(max_tm<T1)
		{
			if(S>-2500.0)
			{
				H_max=H;
				S_max=S;
				max_tm=T1;
			}
		}
	}
	if(hs==1)
		return H_max;
	return S_max;
}

__device__ double END5_2(int i,int hs,int length,double *d_DPT,char *d_numSeq,int id)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=1.0*INFINITY;
	max_tm=-1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-6;++k)
	{
		T1=(d_DPT[id*1340+1276+k]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1250+k]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		T2=d_DPT[id*1340+1302]/(d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(T1>=T2)
		{
			H=d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Hd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-1];
			S=d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Sd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-1];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		}
		else
		{
			H=parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Hd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-1];
			S=parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Sd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-1];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		}

		if(max_tm<T1)
		{
			if(S>-2500.0)
			{
				H_max=H;
				S_max=S;
				max_tm=T1;
			}
		}
	}
	if(hs==1)
		return H_max;
	return S_max;
}

__device__ double END5_3(int i,int hs,int length,double *d_DPT,char *d_numSeq,int id)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=1.0*INFINITY;
	max_tm=-1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-6;++k)
	{
		T1=(d_DPT[id*1340+1276+k]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1250+k]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		T2=d_DPT[id*1340+1302]/(d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(T1>=T2)
		{
			H=d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Hd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+k*(length-1)+i-2];
			S=d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Sd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+625+k*(length-1)+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		}
		else
		{
			H=parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Hd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+k*(length-1)+i-2];
			S=parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Sd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+625+k*(length-1)+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		}

		if(max_tm<T1)
		{
			if(S>-2500.0)
			{
				H_max=H;
				S_max=S;
				max_tm=T1;
			}
		}
	}
	if(hs==1)
		return H_max;
	return S_max;
}

__device__ double END5_4(int i,int hs,int length,double *d_DPT,char *d_numSeq,int id)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=1.0*INFINITY;
	max_tm=-1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-7;++k)
	{
		T1=(d_DPT[id*1340+1276+k]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1250+k]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		T2=d_DPT[id*1340+1302]/(d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(T1>=T2)
		{
			H=d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Htstack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-2];
			S=d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Ststack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		}
		else
		{
			H=parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Htstack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-2];
			S=parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Ststack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/(S+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
 		}

		if(max_tm<T1)
		{
			if(S>-2500.0)
			{
				H_max=H;
				S_max=S;
				max_tm=T1;
			}
		}
	}
	if(hs==1)
		return H_max;
	return S_max;
}

__device__ void calc_terminal_bp(double temp,int length,double *d_DPT,char *d_numSeq,int id)
{
	int i,max;
//back here
	
	d_DPT[id*1340+1250+0]=d_DPT[id*1340+1250+1]= -1.0;
	d_DPT[id*1340+1276+0]=d_DPT[id*1340+1276+1]=1.0*INFINITY;

	for(i=2;i<=length;i++)
	{
		d_DPT[id*1340+1250+i]=-3224.0;
		d_DPT[id*1340+1276+i]=0;
	}

// adding terminal penalties to 3' end and to 5' end 
	for(i=2;i<=length;++i)
	{
		max=0;
		d_DPT[id*1340+1310]=(d_DPT[id*1340+1276+i-1]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1250+i-1]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		d_DPT[id*1340+1315]=END5_1(i,1,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1316]=END5_1(i,2,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1311]=(d_DPT[id*1340+1315]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1316]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		d_DPT[id*1340+1317]=END5_2(i,1,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1318]=END5_2(i,2,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1312]=(d_DPT[id*1340+1317]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1318]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		d_DPT[id*1340+1319]=END5_3(i,1,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1320]=END5_3(i,2,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1313]=(d_DPT[id*1340+1319]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1320]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		d_DPT[id*1340+1321]=END5_4(i,1,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1322]=END5_4(i,2,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1314]=(d_DPT[id*1340+1321]+d_DPT[id*1340+1302])/(d_DPT[id*1340+1322]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);

		max=max5(d_DPT[id*1340+1310],d_DPT[id*1340+1311],d_DPT[id*1340+1312],d_DPT[id*1340+1313],d_DPT[id*1340+1314]);
		switch(max)
		{
			case 1:
				d_DPT[id*1340+1250+i]=d_DPT[id*1340+1250+i-1];
				d_DPT[id*1340+1276+i]=d_DPT[id*1340+1276+i-1];
				break;
			case 2:
				if(d_DPT[id*1340+1315]<temp*d_DPT[id*1340+1316])
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1316];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1315];
				}
				else
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1250+i-1];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1276+i-1];
				}
				break;
			case 3:
				if(d_DPT[id*1340+1317]<temp*d_DPT[id*1340+1318])
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1318];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1317];
				}
				else
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1250+i-1];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1276+i-1];
				}
				break;
			case 4:
				if(d_DPT[id*1340+1319]<temp*d_DPT[id*1340+1320])
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1320];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1319];
				}
				else
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1250+i-1];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1276+i-1];
				}
				break;
			case 5:
				if(d_DPT[id*1340+1321]<temp*d_DPT[id*1340+1322])
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1322];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1321];
				}
				else
				{
					d_DPT[id*1340+1250+i]=d_DPT[id*1340+1250+i-1];
					d_DPT[id*1340+1276+i]=d_DPT[id*1340+1276+i-1];
				}
				break;
			default:
				break;
		}
	}
}

__device__ int newpush(int store[],int i,int j,int mtrx,int total,int next)
{
        int k;
        for(k=total-1;k>=next;k--)
        {
                store[(k+1)*3]=store[k*3];
                store[(k+1)*3+1]=store[k*3+1];
                store[(k+1)*3+2]=store[k*3+2];
        }
        store[next*3]=i;                  
        store[next*3+1]=j;
        store[next*3+2]=mtrx;

        return total+1;           
}

__device__ int equal(double a,double b)
{
	if(fabs(a)>999999999||fabs(b)>999999999)
		return 0;
	return fabs(a-b)<1e-5;
}

__device__ void tracebacku(int *d_ps,int length,double *d_DPT,char *d_numSeq,int id)
{
	int i,j,store[50],total,now,ii,jj,k,d,done;

        total=newpush(store,length,0,1,0,0);
        now=0;
        while(now<total)
        {
                i=store[3*now]; // top->i;
                j=store[3*now+1]; // top->j;
                if(store[now*3+2]==1)
                {
                        while(equal(d_DPT[id*1340+1250+i],d_DPT[id*1340+1250+i-1])&&equal(d_DPT[id*1340+1276+i],d_DPT[id*1340+1276+i-1])) // if previous structure is the same as this one
                                --i;
                        if(i==0)
                                continue;
                        if(equal(d_DPT[id*1340+1250+i],END5_1(i,2,length,d_DPT,d_numSeq,id))&&equal(d_DPT[id*1340+1276+i],END5_1(i,1,length,d_DPT,d_numSeq,id)))
                        {
                                for(k=0;k<=i-5;++k)
                                        if(equal(d_DPT[id*1340+1250+i],parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+625+k*(length-1)+i-1])&&equal(d_DPT[id*1340+1276+i],parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+k*(length-1)+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);                    
                                                break;
                                        }
                                        else if(equal(d_DPT[id*1340+1250+i],d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+625+k*(length-1)+i-1])&&equal(d_DPT[id*1340+1276+i],d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i]]+d_DPT[id*1340+k*(length-1)+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(d_DPT[id*1340+1250+i],END5_2(i,2,length,d_DPT,d_numSeq,id))&&equal(d_DPT[id*1340+1276+i],END5_2(i,1,length,d_DPT,d_numSeq,id)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(d_DPT[id*1340+1250+i],parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Sd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-1])&&equal(d_DPT[id*1340+1276+i],parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Hd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(d_DPT[id*1340+1250+i],d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Sd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-1])&&equal(d_DPT[id*1340+1276+i],d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i]]+Hd5(i,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(d_DPT[id*1340+1250+i],END5_3(i,2,length,d_DPT,d_numSeq,id))&&equal(d_DPT[id*1340+1276+i],END5_3(i,1,length,d_DPT,d_numSeq,id)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(d_DPT[id*1340+1250+i],parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Sd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+625+k*(length-1)+i-2])&&equal(d_DPT[id*1340+1276+i],parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Hd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+k*(length-1)+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(d_DPT[id*1340+1250+i],d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Sd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+625+k*(length-1)+i-2])&&equal(d_DPT[id*1340+1276+i],d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+1]*5+d_numSeq[id*54+i-1]]+Hd3(i-1,k+1,d_numSeq,id)+d_DPT[id*1340+k*(length-1)+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(d_DPT[id*1340+1250+i],END5_4(i,2,length,d_DPT,d_numSeq,id))&&equal(d_DPT[id*1340+1276+i],END5_4(i,1,length,d_DPT,d_numSeq,id)))
                        {
                                for (k=0;k<=i-7;++k)
                                        if(equal(d_DPT[id*1340+1250+i],parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Ststack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-2])&&equal(d_DPT[id*1340+1276+i],parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Htstack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-2]))
                                        {
                                                total=newpush(store,k+2,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(d_DPT[id*1340+1250+i],d_DPT[id*1340+1250+k]+parameter[5680+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Ststack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+625+(k+1)*(length-1)+i-2])&&equal(d_DPT[id*1340+1276+i],d_DPT[id*1340+1276+k]+parameter[5705+d_numSeq[id*54+k+2]*5+d_numSeq[id*54+i-1]]+Htstack(i-1,k+2,d_numSeq,id)+d_DPT[id*1340+(k+1)*(length-1)+i-2]))
                                        {
                                                total=newpush(store,k+2,i-1,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                }
                else if(store[3*now+2]==0)
                {
                        d_ps[id*50+i-1]=j;
                        d_ps[id*50+j-1]=i;
                        d_DPT[id*1340+1310]=-1.0;
                        d_DPT[id*1340+1311]=1.0*INFINITY;
                        calc_hairpin(i,j,1310,1,length,d_DPT,d_numSeq,id);

                        d_DPT[id*1340+1312]=-1.0;
                        d_DPT[id*1340+1313]=1.0*INFINITY;
                        CBI(i,j,1312,2,length,d_DPT,d_numSeq,id);

                        if (equal(d_DPT[id*1340+625+(i-1)*(length-1)+j-1],Ss(i,j,2,length,d_numSeq,id)+d_DPT[id*1340+625+i*(length-1)+j-2])&&equal(d_DPT[id*1340+(i-1)*(length-1)+j-1],Hs(i,j,2,length,d_numSeq,id)+d_DPT[id*1340+i*(length-1)+j-2]))
                                total=newpush(store,i+1,j-1,0,total,now+1);
                        else if(equal(d_DPT[id*1340+625+(i-1)*(length-1)+j-1],d_DPT[id*1340+1312])&&equal(d_DPT[id*1340+(i-1)*(length-1)+j-1],d_DPT[id*1340+1313]))
                        {
                                for (done=0,d=j-i-3;d>=4&&d>=j-i-32&&!done;--d)
                                        for (ii=i+1;ii<j-d;++ii)
                                        {
                                                jj=d+ii;
                                                d_DPT[id*1340+1314]=-1.0;
                                                d_DPT[id*1340+1315]=1.0*INFINITY;
                                                calc_bulge_internal2(i,j,ii,jj,1314,1,length,d_DPT,d_numSeq,id);

                                                if (equal(d_DPT[id*1340+625+(i-1)*(length-1)+j-1],d_DPT[id*1340+1314]+d_DPT[id*1340+625+(ii-1)*(length-1)+jj-1])&&equal(d_DPT[id*1340+(i-1)*(length-1)+j-1],d_DPT[id*1340+1315]+d_DPT[id*1340+(ii-1)*(length-1)+jj-1]))
                                                {
                                                        total=newpush(store,ii,jj,0,total,now+1);
                                                        ++done;
                                                        break;
                                                }
                                        }
                        }
                }
                now++;
        }
}

__device__ double drawHairpin(int *d_ps,int id,double mh,double ms,int length)
{
        int i,N;

        N=0;
        if(fabs(ms)>999999999||fabs(mh)>999999999)
        {
		return 0.0;
        }
        else
        {
		for(i=1;i<length;++i)
		{
			if(d_ps[id*50+i-1]>0)
				N++;
                }
                return mh/(ms+(((N/2)-1)*-0.51986))-273.15;
        }
}

__device__ void initMatrix(int length,double *d_DPT,char *d_numSeq,int id)
{
	int i,j;

	for(i=1;i<=length;++i)
	{
		for(j=1;j<=length;++j)
		{
			if(d_numSeq[id*54+i]+d_numSeq[id*54+27+j]!=3)
			{
				d_DPT[id*1340+(i-1)*length+j-1]=1.0*INFINITY;
				d_DPT[id*1340+625+(i-1)*length+j-1]=-1.0;
			}
			else
			{
				d_DPT[id*1340+(i-1)*length+j-1]=0.0;
				d_DPT[id*1340+625+(i-1)*length+j-1]=-3224.0;
			}
		}
	}
}

__device__ void LSH(int i,int j,double *EntropyEnthalpy,int length,double *d_DPT,char *d_numSeq,int id)
{
	double S1,H1,T1,S2,H2,T2;

	if(d_numSeq[id*54+i]+d_numSeq[id*54+27+j]!=3)
	{
		d_DPT[id*1340+625+(i-1)*length+j-1]=-1.0;
		d_DPT[id*1340+(i-1)*length+j-1]=1.0*INFINITY;
		return;
	}

	S1=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[4430+d_numSeq[id*54+27+j]*125+d_numSeq[id*54+27+j-1]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]];
	H1=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[5055+d_numSeq[id*54+27+j]*125+d_numSeq[id*54+27+j-1]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]];
	if(fabs(H1)>999999999)
	{
		H1=1.0*INFINITY;
		S1=-1.0;
	}
// If there is two dangling ends at the same end of duplex
	if(fabs(parameter[2625+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+27+j-1]*5+d_numSeq[id*54+i]])<999999999&&fabs(parameter[2875+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]])<999999999)
	{
		S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2500+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+27+j-1]*5+d_numSeq[id*54+i]]+parameter[2750+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]];
		H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2625+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+27+j-1]*5+d_numSeq[id*54+i]]+parameter[2875+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if(T1<T2)
			{
				S1=S2;
				H1=H2;
				T1=T2;
			}
		}
		else
		{
			S1=S2;
			H1=H2;
			T1=T2;
		}
	}
	else if(fabs(parameter[2625+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+27+j-1]*5+d_numSeq[id*54+i]])<999999999)
	{
		S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2500+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+27+j-1]*5+d_numSeq[id*54+i]];
		H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2625+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+27+j-1]*5+d_numSeq[id*54+i]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if(T1<T2)
			{
				S1=S2;
				H1=H2;
				T1=T2;
			}
		}
		else
		{
			S1=S2;
			H1=H2;
			T1=T2;
		}
	}
	else if(fabs(parameter[2875+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]])<999999999)
	{
		S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2750+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]];
		H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2875+d_numSeq[id*54+27+j]*25+d_numSeq[id*54+i]*5+d_numSeq[id*54+i-1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if(T1<T2)
			{
				S1=S2;
				H1=H2;
				T1=T2;
			}
		}
		else
		{
			S1=S2;
			H1=H2;
			T1=T2;
		}
	}

	S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]];
	H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]];
	T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
	if(fabs(H1)<999999999)
	{
		if(T1<T2)
		{
			EntropyEnthalpy[0]=S2;
			EntropyEnthalpy[1]=H2;
		}
		else
		{
			EntropyEnthalpy[0]=S1;
			EntropyEnthalpy[1]=H1;
		}
	}
	else
	{
		EntropyEnthalpy[0]=S2;
		EntropyEnthalpy[1]=H2;
	}
	return;
}

__device__ void maxTM(int i,int j,int length,double *d_DPT,char *d_numSeq,int id)
{
	double T0,T1,S0,S1,H0,H1;

	S0=d_DPT[id*1340+625+(i-1)*length+j-1];
	H0=d_DPT[id*1340+(i-1)*length+j-1];
	T0=(H0+d_DPT[id*1340+1302])/(S0+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]); // at current position 
	if(fabs(d_DPT[id*1340+(i-2)*length+j-2])<999999999&&fabs(Hs(i-1,j-1,1,length,d_numSeq,id))<999999999)
	{
		S1=(d_DPT[id*1340+625+(i-2)*length+j-2]+Ss(i-1,j-1,1,length,d_numSeq,id));
		H1=(d_DPT[id*1340+(i-2)*length+j-2]+Hs(i-1,j-1,1,length,d_numSeq,id));
	}
	else
	{
		S1=-1.0;
		H1=1.0*INFINITY;
	}
	T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);

	if(S1<-2500.0)
	{
// to not give dH any value if dS is unreasonable
		S1=-3224.0;
		H1=0.0;
	}
	if(S0<-2500.0)
	{
// to not give dH any value if dS is unreasonable
		S0=-3224.0;
		H0=0.0;
	}
	if((T1>T0)||(S0>0&&H0>0)) // T1 on suurem 
	{
		d_DPT[id*1340+625+(i-1)*length+j-1]=S1;
		d_DPT[id*1340+(i-1)*length+j-1]=H1;
	}
	else if(T0>=T1)
	{
		d_DPT[id*1340+625+(i-1)*length+j-1]=S0;
		d_DPT[id*1340+(i-1)*length+j-1]=H0;
	}
}

__device__ void calc_bulge_internal(int i,int j,int ii,int jj,double* EntropyEnthalpy,int traceback,int length,double *d_DPT,char *d_numSeq,int id)
{
	int loopSize1,loopSize2,loopSize,N,N_loop;
	double T1,T2,S,H;

	S=-3224.0;
	H=0;
	loopSize1=ii-i-1;
	loopSize2=jj-j-1;
	if(ii<jj)
	{
		N=i;
		N_loop=N;
		if(loopSize1>2)
			N_loop-=(loopSize1-2);
		if(loopSize2>2)
			N_loop-=(loopSize2-2);
	}
	else
	{
		N=j;
		N_loop=2*jj;
		if(loopSize1>2)
			N_loop-=(loopSize1-2);
		if(loopSize2>2)
			N_loop-=(loopSize2-2);
		N_loop=(N_loop/2)-1;
	}

	loopSize=loopSize1+loopSize2-1;
	if((loopSize1==0&&loopSize2>0)||(loopSize2==0&&loopSize1>0))// only bulges have to be considered
	{
		if(loopSize2==1||loopSize1==1) // bulge loop of size one is treated differently the intervening nn-pair must be added
		{
			if((loopSize2==1&&loopSize1==0)||(loopSize2==0&&loopSize1==1))
			{
				H=parameter[3150+loopSize]+parameter[625+d_numSeq[id*54+i]*125+d_numSeq[id*54+ii]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+jj]];
				S=parameter[3060+loopSize]+parameter[d_numSeq[id*54+i]*125+d_numSeq[id*54+ii]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+jj]];
			}
			H+=d_DPT[id*1340+(i-1)*length+j-1];
			S+=d_DPT[id*1340+625+(i-1)*length+j-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}

			T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
			T2=(d_DPT[id*1340+(ii-1)*length+jj-1]+d_DPT[id*1340+1302])/((d_DPT[id*1340+625+(ii-1)*length+jj-1])+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
		else // we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30
		{
			H=parameter[3150+loopSize]+parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[5705+d_numSeq[id*54+ii]*5+d_numSeq[id*54+27+jj]];
			H+=d_DPT[id*1340+(i-1)*length+j-1];

			S=parameter[3060+loopSize]+parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[5680+d_numSeq[id*54+ii]*5+d_numSeq[id*54+27+jj]];
			S+=d_DPT[id*1340+625+(i-1)*length+j-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
			T2=(d_DPT[id*1340+(ii-1)*length+jj-1]+d_DPT[id*1340+1302])/(d_DPT[id*1340+625+(ii-1)*length+jj-1]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
	}
	else if(loopSize1==1&&loopSize2==1)
	{
		S=parameter[1250+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]]+parameter[1250+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj-1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]];
		S+=d_DPT[id*1340+625+(i-1)*length+j-1];

		H=parameter[1875+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]]+parameter[1875+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj-1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]];
		H+=d_DPT[id*1340+(i-1)*length+j-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
		T2=(d_DPT[id*1340+(ii-1)*length+jj-1]+d_DPT[id*1340+1302])/(d_DPT[id*1340+625+(ii-1)*length+jj-1]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if((T1-T2>=0.000001)||traceback==1)
		{
			if((T1>T2)||(traceback&&T1>=T2))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
		return;
	}
	else // only internal loops
	{
		H=parameter[3120+loopSize]+parameter[3805+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]]+parameter[3805+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj-1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]];
		H+=d_DPT[id*1340+(i-1)*length+j-1];

		S=parameter[3030+loopSize]+parameter[3180+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]]+parameter[3180+d_numSeq[id*54+27+jj]*125+d_numSeq[id*54+27+jj-1]*25+d_numSeq[id*54+ii]*5+d_numSeq[id*54+ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
		S+=d_DPT[id*1340+625+(i-1)*length+j-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+d_DPT[id*1340+1302])/((S+d_DPT[id*1340+1303])+d_DPT[id*1340+1304]);
		T2=(d_DPT[id*1340+(ii-1)*length+jj-1]+d_DPT[id*1340+1302])/((d_DPT[id*1340+625+(ii-1)*length+jj-1])+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
		{
			EntropyEnthalpy[0]=S;
			EntropyEnthalpy[1]=H;
		}
	}
	return;
}

__device__ void fillMatrix(int length,double *d_DPT,char *d_numSeq,int id)
{
	int d,i,j,ii,jj;
	double SH[2];

	for(i=1;i<=length;++i)
	{
		for(j=1;j<=length;++j)
		{
			if(fabs(d_DPT[id*1340+(i-1)*length+j-1])<999999999)
			{
				SH[0]=-1.0;
				SH[1]=1.0*INFINITY;
				LSH(i,j,SH,length,d_DPT,d_numSeq,id);

				if(fabs(SH[1])<999999999)
				{
					d_DPT[id*1340+625+(i-1)*length+j-1]=SH[0];
					d_DPT[id*1340+(i-1)*length+j-1]=SH[1];
				}
				if(i>1&&j>1)
				{
					maxTM(i,j,length,d_DPT,d_numSeq,id);
					for(d=3;d<=32;d++)
					{
						ii=i-1;
						jj=-ii-d+(j+i);
						if(jj<1)
						{
							ii-=abs(jj-1);
							jj=1;
						}
						for(;ii>0&&jj<j;--ii,++jj)
						{
							if(fabs(d_DPT[id*1340+(ii-1)*length+jj-1])<999999999)
							{
								SH[0]=-1.0;
								SH[1]=1.0*INFINITY;
								calc_bulge_internal(ii,jj,i,j,SH,0,length,d_DPT,d_numSeq,id);

								if(SH[0]<-2500.0)
								{
									SH[0] =-3224.0;
									SH[1] = 0.0;
								}
								if(fabs(SH[1])<999999999)
								{
									d_DPT[id*1340+(i-1)*length+j-1]=SH[1];
									d_DPT[id*1340+625+(i-1)*length+j-1]=SH[0];
								}
							}
						}
					}
				} // if 
			}
		} // for 
	} //for
}

__device__ void RSH(int i,int j,double *d_DPT,char *d_numSeq,int id)
{
	double S1,S2,H1,H2,T1,T2;

	if(d_numSeq[id*54+i]+d_numSeq[id*54+27+j]!=3)
	{
		d_DPT[id*1340+1306]=-1.0;
		d_DPT[id*1340+1307]=1.0*INFINITY;
		return;
	}
	S1=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[4430+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
	H1=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[5055+d_numSeq[id*54+i]*125+d_numSeq[id*54+i+1]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
	if(fabs(H1)>999999999)
	{
		H1=1.0*INFINITY;
		S1=-1.0;
	}
	if(fabs(parameter[2625+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+27+j]])<999999999&&fabs(parameter[2875+d_numSeq[id*54+i]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]])<999999999)
	{
		S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2500+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+27+j]]+parameter[2750+d_numSeq[id*54+i]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
		H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2625+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+27+j]]+parameter[2875+d_numSeq[id*54+i]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if(T1<T2)
			{
				S1=S2;
				H1=H2;
				T1=T2;
			}
		}
		else
		{
			S1=S2;
			H1=H2;
			T1=T2;
		}
	}

	if(fabs(parameter[2625+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+27+j]])<999999999)
	{
		S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2500+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+27+j]];
		H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2625+d_numSeq[id*54+i]*25+d_numSeq[id*54+i+1]*5+d_numSeq[id*54+27+j]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if(T1<T2)
			{
				S1=S2;
				H1=H2;
				T1=T2;
			}
		}
		else
		{
			S1=S2;
			H1=H2;
			T1=T2;
		}
	}

	if(fabs(parameter[2875+d_numSeq[id*54+i]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]])<999999999)
	{
		S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2750+d_numSeq[id*54+i]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
		H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]]+parameter[2875+d_numSeq[id*54+i]*25+d_numSeq[id*54+27+j]*5+d_numSeq[id*54+27+j+1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+d_DPT[id*1340+1302])/(S1+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
			if(T1<T2)
			{
				S1=S2;
				H1=H2;
				T1=T2;
			}
		}
		else
		{
			S1=S2;
			H1=H2;
			T1=T2;
		}
	}
	S2=parameter[5680+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]];
	H2=parameter[5705+d_numSeq[id*54+i]*5+d_numSeq[id*54+27+j]];
	T2=(H2+d_DPT[id*1340+1302])/(S2+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]);
	if(fabs(H1)<999999999)
	{
		if(T1<T2)
		{
			d_DPT[id*1340+1306]=S2;
			d_DPT[id*1340+1307]=H2;
		}
		else
		{
			d_DPT[id*1340+1306]=S1;
			d_DPT[id*1340+1307]=H1;
		}
	}
	else
	{
		d_DPT[id*1340+1306]=S2;
		d_DPT[id*1340+1307]=H2;
	}
	return;
}

__device__ void traceback(int i,int j,int *d_ps,int length,double *d_DPT,char *d_numSeq,int id)
{
	int d,ii,jj,done;
	double SH[2];

	d_ps[id*50+i-1]=j;
	d_ps[id*50+25+j-1]=i;
	while(1)
	{
		SH[0]=-1.0;
		SH[1]=1.0*INFINITY;
		LSH(i,j,SH,length,d_DPT,d_numSeq,id);
		if(equal(d_DPT[id*1340+625+(i-1)*length+j-1],SH[0])&&equal(d_DPT[id*1340+(i-1)*length+j-1],SH[1]))
			break;

		done = 0;
		if(i>1&&j>1&&equal(d_DPT[id*1340+625+(i-1)*length+j-1],Ss(i-1,j-1,1,length,d_numSeq,id)+d_DPT[id*1340+625+(i-2)*length+j-2]))
		{
			i=i-1;
			j=j-1;
			d_ps[id*50+i-1]=j;
			d_ps[id*50+25+j-1]=i;
			done=1;
		}
		for(d=3;!done&&d<=32;++d)
		{
			ii=i-1;
			jj=-ii-d+(j+i);
			if(jj<1)
			{
				ii-=abs(jj-1);
				jj=1;
			}
			for(;!done&&ii>0&&jj<j;--ii,++jj)
			{
				SH[0]=-1.0;
				SH[1]=1.0*INFINITY;
				calc_bulge_internal(ii,jj,i,j,SH,1,length,d_DPT,d_numSeq,id);
				if(equal(d_DPT[id*1340+625+(i-1)*length+j-1],SH[0])&&equal(d_DPT[id*1340+(i-1)*length+j-1],SH[1]))
				{
					i=ii;
					j=jj;
					d_ps[id*50+i-1]=j;
					d_ps[id*50+25+j-1]=i;
					done=1;
					break;
				}
			}
		}
	}
}

__device__ double drawDimer(int *d_ps,int id,double H,double S,double *d_DPT,int length)
{
        int i,N;

        if(fabs(d_DPT[id*1340+1305])>999999999)
                return (double)0.0;
        else
        {
                N=0;
                for(i=0;i<length;i++)
                {
                        if(d_ps[id*50+i]>0)
                                ++N;
                }
                for(i=0;i<length;i++)
                {
                        if(d_ps[id*50+25+i]>0)
                                ++N;
                }
                N=(N/2)-1;
                return (double)(H/(S+(N*-0.51986)+d_DPT[id*1340+1304])-273.15);
        }
}

__device__ int symmetry_thermo(char *d_seq,int start,int length)
{
	int i = 0;
	if(length%2==1)
		return 0;

	while(i<length/2)
	{
		if((d_seq[start+i]=='A'&&d_seq[start+length-1-i]!='T')||(d_seq[start+i]=='T'&&d_seq[start+length-1-i]!='A')||(d_seq[start+length-1-i]=='A'&&d_seq[start+i]!='T')||(d_seq[start+length-1-i]=='T'&&d_seq[start+i]!='A'))
			return 0;
		if((d_seq[start+i]=='C'&&d_seq[start+length-1-i]!='G')||(d_seq[start+i]=='G'&&d_seq[start+length-1-i]!='C')||(d_seq[start+length-1-i]=='C'&&d_seq[start+i]!='G')||(d_seq[start+length-1-i]=='G'&&d_seq[start+i]!='C'))
			return 0;
		i++;
	}
	return 1;
}

__device__ void thal(char *d_seq,int start,int length,int strand_flag,int type,char *d_numSeq,int id,double *d_DPT,int *d_ps)
{
	int Initint[5]; //0 is len1, 1 is len2, 2 is len3, 3 is bestI, 4 is bestJ
	int i, j;

/*** INIT values for unimolecular and bimolecular structures ***/
	if (type==4) /* unimolecular folding */
	{
		d_DPT[id*1340+1302]= 0.0;
		d_DPT[id*1340+1303] = -0.00000000001;
		d_DPT[id*1340+1304]=0;
	}
	else /* hybridization of two oligos */
	{
		d_DPT[id*1340+1302]= 200;
		d_DPT[id*1340+1303]= -5.7;
		if(symmetry_thermo(d_seq,start,length))
			d_DPT[id*1340+1304]=1.9872* log(38/1000000000.0);
		else
			d_DPT[id*1340+1304]=1.9872* log(38/4000000000.0);
	}
/* convert nucleotides to numbers */
	if(type==1 || type==2)
	{
		if(strand_flag==0) //plus
		{
	 		for(i=1;i<=length;++i)
			{
				str2int(d_seq[start+i-1],d_numSeq,(id*54+i));
				str2int(d_seq[start+length-i],d_numSeq,(id*54+27+i));
			}
		}
		else
		{
			for(i=1;i<=length;++i)
			{
				str2int_rev(d_seq[start+length-i],d_numSeq,(id*54+i));
				str2int_rev(d_seq[start+i-1],d_numSeq,(id*54+27+i));
			}
		}
	}
	else
	{
		if(strand_flag==0)
		{
                	for(i=1;i<=length;++i)
			{
				str2int(d_seq[start+i-1],d_numSeq,(id*54+i));
				d_numSeq[id*54+27+i]=d_numSeq[id*54+i];
			}
		}
		else
		{
			for(i=1;i<=length;++i)
			{
				str2int_rev(d_seq[start+length-i],d_numSeq,(id*54+i));
				d_numSeq[id*54+27+i]=d_numSeq[id*54+i];
			}
		}
	}
	d_numSeq[id*54+0]=d_numSeq[id*54+length+1]=d_numSeq[id*54+27+0]=d_numSeq[id*54+27+length+1]=4; /* mark as N-s */

	d_DPT[id*1340+1309]=0;
	if (type==4) /* calculate structure of monomer */
	{
		initMatrix2(length,d_DPT,d_numSeq,id);
		fillMatrix2(length,d_DPT,d_numSeq,id);
		calc_terminal_bp(310.15,length,d_DPT,d_numSeq,id);
		d_DPT[id*1340+1306]=d_DPT[id*1340+1276+length];
		d_DPT[id*1340+1307]=d_DPT[id*1340+1250+length];
		for (i=0;i<length;i++)
			d_ps[id*50+i]=0;
		if(fabs(d_DPT[id*1340+1306])<999999999)
		{
			tracebacku(d_ps,length,d_DPT,d_numSeq,id);
			d_DPT[id*1340+1309]=drawHairpin(d_ps,id,d_DPT[id*1340+1306],d_DPT[id*1340+1307],length);
			d_DPT[id*1340+1309]=(int)(d_DPT[id*1340+1309]*100+0.5)/100.0;
		}
	}
	else  /* Hybridization of two moleculs */
	{
		initMatrix(length,d_DPT,d_numSeq,id);
		fillMatrix(length,d_DPT,d_numSeq,id);

		d_DPT[id*1340+1305]=-1.0*INFINITY;
	/* calculate terminal basepairs */
		Initint[3]=Initint[4]=0;
		if(type==1)
			for (i=1;i<=length;i++)
			{
				for (j=1;j<=length;j++)
				{
					RSH(i,j,d_DPT,d_numSeq,id);
					d_DPT[id*1340+1306]=d_DPT[id*1340+1306]+0.000001; /* this adding is done for compiler, optimization -O2 vs -O0 */
					d_DPT[id*1340+1307]=d_DPT[id*1340+1307]+0.000001;
					d_DPT[id*1340+1308]=((d_DPT[id*1340+(i-1)*length+j-1]+d_DPT[id*1340+1307]+d_DPT[id*1340+1302]) / ((d_DPT[id*1340+625+(i-1)*length+j-1])+d_DPT[id*1340+1306]+d_DPT[id*1340+1303] + d_DPT[id*1340+1304])) -273.15;
					if(d_DPT[id*1340+1308]>d_DPT[id*1340+1305]&&((d_DPT[id*1340+625+(i-1)*length+j-1]+d_DPT[id*1340+1306])<0&&(d_DPT[id*1340+1307]+d_DPT[id*1340+(i-1)*length+j-1])<0))
					{
						d_DPT[id*1340+1305]=d_DPT[id*1340+1308];
						Initint[3]=i;
						Initint[4]=j;
					}
				}
			}
		if(type==2)
		{
		 //THAL_END1
			Initint[4]=0;
			Initint[3]=length;
			i=length;
			d_DPT[id*1340+1305]=-1.0*INFINITY;
			for (j=1;j<=length;++j)
			{
				RSH(i,j,d_DPT,d_numSeq,id);
				d_DPT[id*1340+1306]=d_DPT[id*1340+1306]+0.000001; // this adding is done for compiler, optimization -O2 vs -O0,that compiler could understand that SH is changed in this cycle 
				d_DPT[id*1340+1307]=d_DPT[id*1340+1307]+0.000001;
				d_DPT[id*1340+1308]=((d_DPT[id*1340+(i-1)*length+j-1]+d_DPT[id*1340+1307]+d_DPT[id*1340+1302])/((d_DPT[id*1340+625+(i-1)*length+j-1])+d_DPT[id*1340+1306]+d_DPT[id*1340+1303]+d_DPT[id*1340+1304]))-273.15;
				if (d_DPT[id*1340+1308]>d_DPT[id*1340+1305]&&((d_DPT[id*1340+1306]+d_DPT[id*1340+625+(i-1)*length+j-1])<0&&(d_DPT[id*1340+1307]+d_DPT[id*1340+(i-1)*length+j-1])<0))
				{
					d_DPT[id*1340+1305]=d_DPT[id*1340+1308];
					Initint[4]=j;
				}
			}
		}
		if(fabs(d_DPT[id*1340+1305])>999999999)
			Initint[3]=Initint[4]=1;
		RSH(Initint[3],Initint[4],d_DPT,d_numSeq,id);
	 // tracebacking 
		for (i=0;i<length;++i)
			d_ps[id*50+i]=0;
		for (j=0;j<length;++j)
			d_ps[id*50+25+j] = 0;
		if(fabs(d_DPT[id*1340+(Initint[3]-1)*length+Initint[4]-1])<999999999)
		{
			traceback(Initint[3],Initint[4],d_ps,length,d_DPT,d_numSeq,id);
			d_DPT[id*1340+1309]=drawDimer(d_ps,id,(d_DPT[id*1340+(Initint[3]-1)*length+Initint[4]-1]+d_DPT[id*1340+1307]+d_DPT[id*1340+1302]),(d_DPT[id*1340+625+(Initint[3]-1)*length+Initint[4]-1]+d_DPT[id*1340+1306]+d_DPT[id*1340+1303]),d_DPT,length);
			d_DPT[id*1340+1309]=(int)(d_DPT[id*1340+1309]*100+0.5)/100.0;
		}
	}
}

///function in gpu, check the GC-content; int length: the length of read
__device__ int gc(char *d_seq,int start,int length)
{
	int i,number;
	float gc;

	number=0;
	for(i=0;i<length;i++)
	{
		if(d_seq[start+i]=='C')
		{
			number++;
			continue;
		}
	
		if(d_seq[start+i]=='G')
		{
			number++;
		}
	}

	gc=1.0*number/length*100;
	if((gc<40)||(gc>65))
	{
		return 0;
	}
	return 1;
}

///function in gpu, translate A...G to int
__device__ int translate(char a)
{
	if(a=='A')
		return 0;
	if(a=='T')
		return 1;
	if(a=='C')
		return 2;
	return 3;
}

__device__ int translate_rev(char a)
{
        if(a=='T')
                return 0;
        if(a=='A')
                return 1;
        if(a=='G')
                return 2;
        return 3;
}
//function in gpu, caculate tm
__device__ int tm(char *d_seq,int start,int length,float max_tm,float min_tm)
{
	int i,pos;
	float deltah,deltas,result;

	deltah=0;
	deltas=0;
	for(i=0;i<length-1;i++)
	{
		pos=translate(d_seq[start+i]);
		pos=pos*4+translate(d_seq[start+i+1]);
		deltah+=d_deltah[pos];
		deltas+=d_deltas[pos];
	}

	deltah=(-1.0)*deltah;
	deltas=(-1.0)*deltas;
	if((d_seq[start]=='A')||(d_seq[start]=='T'))
	{
		deltah+=2.3;
		deltas+=4.1;
	}
	else
	{
		deltah+=0.1;
		deltas-=2.8;
	}
        if((d_seq[start+length-1]=='A')||(d_seq[start+length-1]=='T'))
        {
                deltah+=2.3;
                deltas+=4.1;
        }
        else
        {
                deltah+=0.1;
                deltas-=2.8;
        }
	result=1000.0*deltah/(deltas-0.51986*(length-1)-36.70381)-273.15;
	if((result<min_tm)||(result>max_tm))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

///function in gpu, caculate stability, int strand: 0 is 5' and 1 is 3'
__device__ int stability(char *d_seq,int start,int flag,int length,int strand)//flag=0: plus
{
	int i,pos;
	
	pos=0;
	for(i=0;i<6;i++)
	{
		if(flag==0)
		{
			if(strand==0)
			{
				pos=pos*4+translate(d_seq[start+i]);
			}
			else
			{
				pos=pos*4+translate(d_seq[start+i+length-6]);
			}
		}
		else
		{
			if(strand==0)
                        {
                                pos=pos*4+translate_rev(d_seq[start+length-1-i]);
                        }
                        else
                        {
                                pos=pos*4+translate_rev(d_seq[start+5-i]);
                        }
		}
	}
	
	if(d_stab[pos]<4)
	{
		return 0;
	}
//the other part
        pos=0;
        for(i=0;i<6;i++)
        {
		if(flag==0)
		{
                	if(strand==1)
                	{
				pos=pos*4+translate(d_seq[start+i]);
                	}
                	else
                	{
				pos=pos*4+translate(d_seq[start+i+length-6]);
                	}
		}
		else //minus
		{
			if(strand==1)
                        {
                                pos=pos*4+translate_rev(d_seq[start+length-1-i]);      
                        }
                        else
                        {
                                pos=pos*4+translate_rev(d_seq[start+5-i]);
                        }
		}
        }

        if(d_stab[pos]<3)
        {
                return 0;
        }

	return 1;
}

//function in gpu: whether species chars in reads
__device__ int words(char *d_seq,int position,int length)
{
	int i;
	
	for(i=0;i<length;i++)
	{
		if(d_seq[position+i]=='N')
		{
			return 0;
		}
	}
	return 1;
}

__device__ int check_long_ploy(char *d_seq,int start,int length)
{
        int i,same;
        char ref;

        same=1;
        ref=d_seq[start];
        for(i=1;i<length;i++)
        {
                if(d_seq[start+i]==ref)
                        same++;
                else
                {
                        if(same>=6)
                                return 0;
                        same=1;
                        ref=d_seq[start+i];
                }
        }
        if(same>=6)
                return 0;
        return 1;
}

///function: int length: the length of genome
__global__ void candidate_primer(char *d_seq,int *d_pos,int *d_len,int *d_rev_len,int strand,float max_tm,float min_tm,int length,int check_flag,char *d_numSeq,double *d_DPT,int *d_ps)
{
	int id,i,circle,check,plus,minus;

	id=threadIdx.x+blockIdx.x*blockDim.x;
	for(circle=id;circle<length;circle=circle+blockDim.x*gridDim.x)
	{
		for(i=0;i<8;i++)   //primer length is from 18 to 25
		{
			d_len[8*circle+i]=0;
			d_rev_len[8*circle+i]=0;
		}
		d_pos[circle]=0;
	
		for(i=18;i<=25;i++)  //read length is from 18 to 25
		{
			if(circle+i>length)
				break;
			check=words(d_seq,circle,i);
			if(check==0)
                                break;

			check=gc(d_seq,circle,i);
			if(check==0)
				continue;

			check=check_long_ploy(d_seq,circle,i);
			if(check==0)
                                continue;

			check=tm(d_seq,circle,i,max_tm,min_tm);
			if(check==0)
				continue;

                        check=stability(d_seq,circle,0,i,strand);
                        if(check==1)     //+ strand
				plus=1;
			else
				plus=0;
			
		//secondary structure
			if(check_flag&&plus)
			{
				thal(d_seq,circle,i,0,1,d_numSeq,id,d_DPT,d_ps);
				if(d_DPT[id*1340+1309]>min_tm-10)
					plus=0;	
			}
			if(check_flag&&plus)
                        {
                                thal(d_seq,circle,i,0,2,d_numSeq,id,d_DPT,d_ps);
				if(d_DPT[id*1340+1309]>min_tm-10)  
                                        plus=0;
                        }
			if(check_flag&&plus)
                        {
                                thal(d_seq,circle,i,0,4,d_numSeq,id,d_DPT,d_ps);
				if(d_DPT[id*1340+1309]>min_tm-10)
                                        plus=0;         
                        }
			if(plus)
                                d_len[circle*8+i-18]=1;
	//reverse
			check=stability(d_seq,circle,1,i,strand);
			if(check==1)
				minus=1;
			else
				minus=0;
		//secondary structure      
                        if(check_flag&&minus)
                        {                
                                thal(d_seq,circle,i,1,1,d_numSeq,id,d_DPT,d_ps);
				if(d_DPT[id*1340+1309]>min_tm-10)
                                        minus=0;         
                        }           
                        if(check_flag&&minus)
                        {
                                thal(d_seq,circle,i,1,2,d_numSeq,id,d_DPT,d_ps);
				if(d_DPT[id*1340+1309]>min_tm-10)
                                        minus=0;
                        }                
                        if(check_flag&&minus)
                        {
                                thal(d_seq,circle,i,1,4,d_numSeq,id,d_DPT,d_ps);
				if(d_DPT[id*1340+1309]>min_tm-10)
                                        minus=0;
                        }
                        if(minus)
				d_rev_len[circle*8+i-18]=1;
		}
		
		for(i=0;i<8;i++)
		{
			d_pos[circle]+=(d_len[circle*8+i]+d_rev_len[8*circle+i]);
		}
	}
	__syncthreads();
}

void usage()
{
        printf("Usage:\n");
        printf("    single  -in <fasta_file>  -out <primers_file_name>  -high[-low] [options]*\n\n");
        printf("    -in   <string>:  the reference sequence file, fasta formate\n");
        printf("    -out  <string>:  the prefix of output files, those files store candidate single primers\n");
        printf("    -dir  <string>:  the directory to store candidate single primers. default is current directory\n");
        printf("    -stab <string>:  the parameter file used in calculating the primers' stability. default is stab_parameter.txt in Par/ directory\n");
        printf("    -tm   <string>:  the parameter file used in calcalating Tm and second structure. default is stab_parameter.txt in Par/ directory\n");
	printf("    -check   <int>:  0: don't check primers' secondary structure; !=0: check, default is 1\n");
        printf("    -par  <string>:  the directory of storing parameter files used to check primers' secondary structure, default is Par/\n");
        printf("    -high/-low:      design candidate single primers in high/low GC region. high: the GC content>=45%%; low: the GC content <=45%%.\n");
        printf("    -loop:           design candidate loop single primers\n");
        printf("    -h/-help:        print usage\n");
}

int create_file(char *prefix,char *dir,char *seq,int *pos,int *len,int *rev_len,int length)
{
	char *file;
	int total,i,j;
	FILE *OUT;

	total=0;
	i=strlen(dir)+strlen(prefix)+20;
	file=(char *)malloc(i);
        memset(file,'\0',i);
        strcpy(file,dir);
        strcat(file,prefix);
        OUT=fopen(file,"w");
        if(OUT==NULL)
        {
                printf("Error! Can't create the %s file!\n",file);
                exit(1);
        }
	
        for(i=0;i<length;i++)
        {
                if(pos[i]==0)
                        continue;
                for(j=0;j<8;j++)
                {
                        if((len[8*i+j]+rev_len[8*i+j])==0)
                                continue;
                       	fprintf(OUT,"pos:%d\tlength:%d\t+:%d\t-:%d\n",i,(j+18),len[8*i+j],rev_len[8*i+j]);
			total++;
                }
        }
	fclose(OUT);
	free(file);
	return total;
}

main(int argc, char **argv)
{
	double *H_parameter,*d_DPT;
	int *pos,*d_pos,*len,*d_len,length,flag[10],i,*rev_len,*d_rev_len,num_outer,num_inner,num_loop,NumL[2],thread,block,*d_ps;
	float deltah[16],deltas[16],stab[4096],temp1,temp2;
	char *seq,*d_seq,*store_path,*prefix,*stab_path,*tm_path,*curren_path,*input,*outer,*inner,*loop,*par_path,*temp,*Pchar,*d_numSeq;
	FILE *fp;
	time_t start,end;
        struct stat statbuf;
//flag: 0:input; 1: out_prefix; 2: dir; 3: stab; 4: tm; 5: high; 6: low; 7: loop; 8: secondary structure; 9: path for secondary structure

	start=time(NULL);
	thread=200;
	block=200;
//get input
        for(i=0;i<10;i++)
        {
                flag[i]=0;
        }
	flag[8]=1;
        for(i=1;i<argc;)
        {
                if(strcmp(argv[i],"-in")==0)
                {
                        flag[0]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-in\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
			length=strlen(argv[i+1]);
			input=(char *)malloc(length+1);
			memset(input,'\0',length+1);
                        strcpy(input,argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-out")==0)
                {
                        flag[1]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-out\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
			length=strlen(argv[i+1]);
                        prefix=(char *)malloc(length+1);
                        memset(prefix,'\0',length+1);
                        strcpy(prefix,argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-dir")==0)
                {
                        flag[2]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-dir\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
			length=strlen(argv[i+1]);
			if(argv[i+1][length-1]=='/')
			{
                        	store_path=(char *)malloc(length+1);
                        	memset(store_path,'\0',length+1);
                        	strcpy(store_path,argv[i+1]);
			}
			else
			{
				store_path=(char *)malloc(length+2);
				memset(store_path,'\0',length+2);
				strcpy(store_path,argv[i+1]);
				store_path[length]='/';
			}
                        i=i+2;
                }
                else if(strcmp(argv[i],"-stab")==0)
                {
                        flag[3]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-stab\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
			length=strlen(argv[i+1]);
                        stab_path=(char *)malloc(length+1);
                        memset(stab_path,'\0',length+1);
                        strcpy(stab_path,argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-tm")==0)
                {
                        flag[4]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-tm\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
			length=strlen(argv[i+1]);
                        tm_path=(char *)malloc(length+1);
                        memset(tm_path,'\0',length+1);
                        strcpy(tm_path,argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-high")==0)
                {
                        flag[5]=1;
                        i++;
                }
                else if(strcmp(argv[i],"-low")==0)
                {
                        flag[6]=1;
                        i++;
                }
                else if(strcmp(argv[i],"-loop")==0) 
                {
                        flag[7]=1;
                        i++;
                }
                else if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-help")==0)
                {
                        usage();
                        exit(1);
                }
		else if(strcmp(argv[i],"-check")==0)
                {
                        if(i+1==argc)
                        {
                                printf("Error! The \"-check\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
                        flag[8]=atoi(argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-par")==0)
                {
                        flag[9]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-par\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
                        length=strlen(argv[i+1]);
                        if(argv[i+1][length-1]=='/')
                        {
                                par_path=(char *)malloc(length+1);
                                strcpy(par_path,argv[i+1]);
                                par_path[length]='\0';
                        }
                        else
                        {
                                par_path=(char *)malloc(length+2);
                                strcpy(par_path,argv[i+1]);
                                par_path[length]='/';
                                par_path[length+1]='\0';
                        }
                        i=i+2;
                }		
                else
                {
                        printf("Error: don't have the parameter: %s\n",argv[i]);
                        usage();
                        exit(1);
                }
        }
//check paramters
        if(flag[5]+flag[6]!=1)
        {
                printf("Error! The input parameter must contain one of -high and -low!\n");
                usage();
                exit(1);
        }
        if(flag[0]==0)
        {
                printf("Error! Users must input the reference sequence file with -in!\n");
                usage();
                exit(1);
        }
        if(flag[1]==0)
        {
                printf("Error! Users must supply the prefix name for output file with -out!\n");
                usage();
                exit(1);
        }
        for(i=0;i<strlen(prefix);i++)
        {
                if(prefix[i]=='/')
                {
                        printf("Error! the -out parameter couldn't contain any directory!\n");
                        usage();
                        exit(1);
                }
        }
//prepare
	inner=(char *)malloc(4096);
        memset(inner,'\0',4096);
        getcwd(inner,4096);
        length=strlen(inner);
        curren_path=(char *)malloc(length+1);
        memset(curren_path,'\0',length+1);
        strcpy(curren_path,inner);
        if(flag[2]==0)
        {
                store_path=(char *)malloc(length+2);
                memset(store_path,'\0',length+2);
                strcpy(store_path,curren_path);
                store_path[length]='/';
        }
        free(inner);

        length=strlen(store_path)+12;
        outer=(char *)malloc(length);
        memset(outer,'\0',length);
        strcpy(outer,store_path);

        inner=(char *)malloc(length);
        memset(inner,'\0',length);
        strcpy(inner,store_path);

        if(flag[7]==1)
        {
                loop=(char *)malloc(length);
                memset(loop,'\0',length);
                strcpy(loop,store_path);
        }
        if(flag[5]==1)
        {
                strcat(outer,"high-outer/");
                strcat(inner,"high-inner/");
                if(flag[7]==1)
                        strcat(loop,"high-loop/");
        }
        else          
        {                
                strcat(outer,"low-outer/");
                strcat(inner,"low-inner/");
                if(flag[7]==1)
                        strcat(loop,"low-loop/");
        }
        mkdir(outer,0755);
        mkdir(inner,0755);        
        if(flag[7]==1)
                mkdir(loop,0755);        

//stability parameter file
        if(flag[3]==0)
        {
		length=strlen(curren_path);
                stab_path=(char *)malloc(length+30);
                memset(stab_path,'\0',length+30);
                strcpy(stab_path,curren_path);
                i=length-1;
                while(stab_path[i]!='/'&&i>=0)
                {
                        stab_path[i]='\0';
                        i--;
                }
                strcat(stab_path,"Par/stab_parameter.txt");
        }
//tm parameter file
        if(flag[4]==0)
        {
		length=strlen(curren_path);
                tm_path=(char *)malloc(length+30);
                memset(tm_path,'\0',length+30);
                strcpy(tm_path,curren_path);
                i=length-1;
                while(tm_path[i]!='/'&&i>=0)
                {
                        tm_path[i]='\0';
                        i--;
                }
                strcat(tm_path,"Par/tm_nn_parameter.txt");
        }
//secondary structure
	if(flag[8]&&flag[9]==0)
        {
                length=strlen(curren_path);
                par_path=(char *)malloc(length+10);
                memset(par_path,'\0',length+10);
                strcpy(par_path,curren_path);
                i=length-1;
                while(par_path[i]!='/'&&i>=0)
                {
                        par_path[i]='\0';
                        i--;
                }
                strcat(par_path,"Par/");
        }
	if(flag[8])
	{
		NumL[0]=get_num_line(par_path,0);
	        NumL[1]=get_num_line(par_path,1);
	        H_parameter=(double *)malloc((5730+2*NumL[0]+2*NumL[1])*sizeof(double));
	        memset(H_parameter,'\0',(5730+2*NumL[0]+2*NumL[1])*sizeof(double));
	        Pchar=(char *)malloc(10*NumL[0]+12*NumL[1]);
	        memset(Pchar,'\0',10*NumL[0]+12*NumL[1]);

		getStack(par_path,H_parameter);
	        getStackint2(par_path,H_parameter);
	        getDangle(par_path,H_parameter);
	        getLoop(par_path,H_parameter);
	        getTstack(par_path,H_parameter);
	        getTstack2(par_path,H_parameter);
	        getTriloop(par_path,H_parameter,Pchar,NumL);
	        getTetraloop(par_path,H_parameter,Pchar,NumL);
	        tableStartATS(6.9,H_parameter);
	        tableStartATH(2200.0,H_parameter);

		cudaMemcpyToSymbol(d_NumL,NumL,2*sizeof(int));
		cudaMemcpyToSymbol(d_Pchar,Pchar,10*NumL[0]+12*NumL[1]);
		cudaMemcpyToSymbol(parameter,H_parameter,(5730+2*NumL[0]+2*NumL[1])*sizeof(double));

		cudaMalloc((void **)&d_numSeq,54*thread*block*sizeof(char));
		cudaMalloc((void **)&d_DPT,1340*thread*block*sizeof(double));
		cudaMalloc((void **)&d_ps,50*thread*block*sizeof(int));
	}

//input reference sequence
        if(access(input,0)==-1)
        {
                printf("Error! Don't have the %s file.\n",input);
                exit(1);
        }
        stat(input,&statbuf);
        length=statbuf.st_size;
        length=length+100;
        temp=(char *)malloc(length);
        memset(temp,'\0',length);
        seq=(char *)malloc(length*sizeof(char));
        memset(seq,'\0',length*sizeof(char));

        fp=fopen(input,"r");   //open the sequence file
        if(fp==NULL)
        {
                printf("Error! can't open the %s file!\n",input);
                exit(1);
        }
        fread(temp,length*sizeof(char),1,fp);
        fclose(fp); 

        length=0;
        i=0;
        while(temp[i]!='\n')
        {
                i++;
        }
        i++;
        while(temp[i]!='\0')
        {
                if(temp[i]=='\n')
                {
                        i++;
                        continue;
                }
		if(temp[i]=='a'||temp[i]=='A')
                        seq[length]='A';
                else if(temp[i]=='t'||temp[i]=='T')
                        seq[length]='T';
                else if(temp[i]=='c'||temp[i]=='C')
                        seq[length]='C';
                else if(temp[i]=='g'||temp[i]=='G')
                        seq[length]='G';
                else
                        seq[length]='N';
                i++;
                length++;
        }
        free(temp);
        length=strlen(seq);

//input Tm parameter
        fp=fopen(tm_path,"r");  //read the paramter of deltah and deltas
        if(fp==NULL)
        {
                printf("Error: can't open the %s file!\n",tm_path);
                exit(1);
        }
        while(fscanf(fp,"%d\t%f\t%f",&i,&temp1,&temp2)!=EOF)
        {
                deltah[i]=temp1;
                deltas[i]=temp2;
        }
        fclose(fp);

//input stability parameter
        fp=fopen(stab_path,"r");  //read the parameters of stability
        if(fp==NULL)
        {
                printf("Error: can't open the %s file!\n",stab_path);
                exit(1);
        }
        while(fscanf(fp,"%d\t%f",&i,&temp1)!=EOF)
        {
                stab[i]=temp1;
        }
        fclose(fp);

	cudaMalloc((void **)&d_seq,length*sizeof(char));
	cudaMemset(d_seq,'\0',length*sizeof(char));

	/////from cpu to gpu
	cudaMemcpy(d_seq,seq,length*sizeof(char),cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_deltah,deltah,16*sizeof(float));
	cudaMemcpyToSymbol(d_deltas,deltas,16*sizeof(float));
	cudaMemcpyToSymbol(d_stab,stab,4096*sizeof(float));

	cudaMalloc((void **)&d_pos,length*sizeof(int));
	cudaMemset(d_pos,'\0',length*sizeof(int));
	cudaMalloc((void **)&d_len,8*length*sizeof(int));
	cudaMemset(d_len,'\0',8*length*sizeof(int));
	cudaMalloc((void **)&d_rev_len,8*length*sizeof(int));
        cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
	pos=(int *)malloc(length*sizeof(int));
	memset(pos,'\0',length*sizeof(int));
	len=(int *)malloc(8*length*sizeof(int));
	memset(len,'\0',8*length*sizeof(int));
        rev_len=(int *)malloc(8*length*sizeof(int));
        memset(rev_len,'\0',8*length*sizeof(int));

	end=time(NULL);
	printf("It takes %d seconds to prepare.\n",(int)difftime(end,start));
	start=time(NULL);
	if(flag[5]==1)
        {
		cudaMemset(d_pos,'\0',length*sizeof(int));
		cudaMemset(d_len,'\0',8*length*sizeof(int));
		cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
		candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,1,61,59,length,flag[8],d_numSeq,d_DPT,d_ps);
		cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
        	cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
        	cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_outer=create_file(prefix,outer,seq,pos,len,rev_len,length);

		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,0,66,64,length,flag[8],d_numSeq,d_DPT,d_ps);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_inner=create_file(prefix,inner,seq,pos,len,rev_len,length);

                if(flag[7]==1)
		{
			cudaMemset(d_pos,'\0',length*sizeof(int));
                	cudaMemset(d_len,'\0',8*length*sizeof(int));
                	cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                	candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,1,66,64,length,flag[8],d_numSeq,d_DPT,d_ps);
                	cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	num_loop=create_file(prefix,loop,seq,pos,len,rev_len,length);
		}
        }
        else
        {
		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,1,56,54,length,flag[8],d_numSeq,d_DPT,d_ps);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_outer=create_file(prefix,outer,seq,pos,len,rev_len,length);

		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,0,61,59,length,flag[8],d_numSeq,d_DPT,d_ps);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_inner=create_file(prefix,inner,seq,pos,len,rev_len,length);
                if(flag[7]==1)
		{
			cudaMemset(d_pos,'\0',length*sizeof(int));
                	cudaMemset(d_len,'\0',8*length*sizeof(int));
                	cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                	candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,1,61,59,length,flag[8],d_numSeq,d_DPT,d_ps);
                	cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	num_loop=create_file(prefix,loop,seq,pos,len,rev_len,length);
		}
        }
	cudaFree(d_pos);
	cudaFree(d_len);
	cudaFree(d_rev_len);
	cudaFree(d_seq);
	free(pos);
        free(len);
        free(rev_len);
	free(seq);

	printf("There ara %d candidate primers used as F3/F2/B2/B3.\n",num_outer);
        printf("There are %d candidate primers used as F1c/B1c.\n",num_inner);
        if(flag[7]==1)
                printf("There are %d candidate primers used as LF/LB.\n",num_loop);
        //check
        if(num_outer<4)
                printf("Warning: there don't have enough primers(>=4) used as F3/F2/B2/B3.\n");
        if(num_inner<2)
                printf("Warning: there don't have enough primers(>=2) used as F1c/B1c.\n");
        if(flag[7]==1 && num_loop<1)
                printf("Warning: there don't have enough primers(>=1) used as LF/LB. But you can design LAMP primers without loop primer.\n");
	end=time(NULL);
        printf("It takes %d seconds to design candidate single primers.\n",(int)difftime(end,start));

	free(store_path);
	free(prefix);
	free(stab_path);
	free(tm_path);
	free(curren_path);
	free(input);
	free(outer);
	free(inner);
	if(flag[7])
		free(loop);

	if(flag[8])
	{
		free(Pchar);
		free(H_parameter);
		cudaFree(d_numSeq);
		cudaFree(d_DPT);
		cudaFree(d_ps);
	}
	if(flag[8]||flag[9])
		free(par_path);
}
