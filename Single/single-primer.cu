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

__device__ char str2int(char c)
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

__device__ int seq_length(char seq[])
{
	int i=0;
	while(seq[i]!='\0')
		i++;
	return i;
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
__device__ void initMatrix2(int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[])
{
	int i,j;
	for(i=1;i<=Initint[0];++i)
		for(j=i;j<=Initint[1];++j)
			if(j-i<4 || (numSeq1[i]+numSeq1[j]!=3))
			{
				enthalpyDPT[(i-1)*Initint[2]+j-1]=1.0*INFINITY;
				entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
			}
			else
			{
				enthalpyDPT[(i-1)*Initint[2]+j-1]=0.0;
				entropyDPT[(i-1)*Initint[2]+j-1]=-3224.0;
			}
}

__device__ double Ss(int i,int j,int k,int Initint[],char numSeq1[],char numSeq2[],double parameter[])
{
	if(k==2)
	{
		if(i>=j)
			return -1.0;
		if(i==Initint[0]||j==Initint[1]+1)
			return -1.0;

		if(i>Initint[0])
			i-=Initint[0];
		if(j>Initint[1])
			j-=Initint[1];
		return parameter[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]];
	}
	else
		return parameter[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
}

__device__ double Hs(int i,int j,int k,int Initint[],char numSeq1[],char numSeq2[],double parameter[])
{
	if(k==2)
	{
		if(i>= j)
			return 1.0*INFINITY;
		if(i==Initint[0]||j==Initint[1]+1)
			return 1.0*INFINITY;

		if(i>Initint[0])
			i-=Initint[0];
		if(j>Initint[1])
			j-=Initint[1];
		if(fabs(parameter[625+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]])<999999999)
			return parameter[625+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]];
		else
			return 1.0*INFINITY;
	}
	else
		return parameter[625+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
}

__device__ void maxTM2(int i,int j,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double parameter[])
{
	double T0,T1,S0,S1,H0,H1;

	S0=entropyDPT[(i-1)*Initint[2]+j-1];
	H0=enthalpyDPT[(i-1)*Initint[2]+j-1];
	T0=(H0+Initdouble[0])/(S0+Initdouble[1]+Initdouble[2]);
	if(fabs(enthalpyDPT[(i-1)*Initint[2]+j-1])<999999999)
	{
		S1=(entropyDPT[i*Initint[2]+j-2]+Ss(i,j,2,Initint,numSeq1,numSeq2,parameter));
		H1=(enthalpyDPT[i*Initint[2]+j-2]+Hs(i,j,2,Initint,numSeq1,numSeq2,parameter));
	}
	else
	{
		S1=-1.0;
		H1=1.0*INFINITY;
	}
	T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
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
		entropyDPT[(i-1)*Initint[2]+j-1]=S1;
		enthalpyDPT[(i-1)*Initint[2]+j-1]= H1;
	}
	else
	{
		entropyDPT[(i-1)*Initint[2]+j-1]=S0;
		enthalpyDPT[(i-1)*Initint[2]+j-1]=H0;
	}
}

__device__ void calc_bulge_internal2(int i,int j,int ii,int jj,double *EntropyEnthalpy,int traceback,double Initdouble[0],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double parameter[])
{
	int loopSize1,loopSize2,loopSize;
	double T1,T2,S,H;

	S=-3224.0;
	H=0.0;
	loopSize1=ii-i-1;
	loopSize2=j-jj-1;
	if(loopSize1+loopSize2>30)
	{
		EntropyEnthalpy[0]=-1.0;
		EntropyEnthalpy[1]=1.0*INFINITY;
		return;
	}

	loopSize=loopSize1+loopSize2-1;
	if((loopSize1==0&&loopSize2>0)||(loopSize2==0&&loopSize1>0))
	{
		if(loopSize2==1||loopSize1==1)
		{ 
			if((loopSize2==1&&loopSize1==0)||(loopSize2==0&&loopSize1==1))
			{
				H=parameter[3150+loopSize]+parameter[625+numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
				S=parameter[3060+loopSize]+parameter[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
 			}
			if(traceback!=1)
			{
				H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];
				S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
			}

			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
			T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/((entropyDPT[(i-1)*Initint[2]+j-1])+Initdouble[1]+Initdouble[2]);
			if((T1>T2)||((traceback&&T1>=T2)||traceback==1))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
		else
		{
			H=parameter[3150+loopSize]+parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[5705+numSeq1[ii]*5+numSeq2[jj]];
			if(traceback!=1)
				H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];

			S=parameter[3060+loopSize]+parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[5680+numSeq1[ii]*5+numSeq2[jj]];
			if(traceback!=1)
				S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
			T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/(entropyDPT[(i-1)*Initint[2]+j-1]+Initdouble[1]+Initdouble[2]);
			if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
	}
	else if(loopSize1==1&&loopSize2==1)
	{
		S=parameter[1250+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+parameter[1250+numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		if(traceback!=1)
			S+=entropyDPT[(ii-1)*Initint[2]+jj-1];

		H=parameter[1875+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+parameter[1875+numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		if(traceback!=1)
			H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
		T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/(entropyDPT[(i-1)*Initint[2]+j-1]+Initdouble[1]+Initdouble[2]);
		if((T1-T2>=0.000001)||traceback)
		{
			if((T1>T2)||((traceback&&T1>= T2)||traceback==1))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
		return;
	}
	else
	{
		H=parameter[3120+loopSize]+parameter[3805+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+parameter[3805+numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		if(traceback!=1)
			H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];

		S=parameter[3030+loopSize]+parameter[3180+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+parameter[3180+numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
		if(traceback!=1)
			S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}

		T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
		T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/((entropyDPT[(i-1)*Initint[2]+j-1])+Initdouble[1]+Initdouble[2]);
		if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
		{
			EntropyEnthalpy[0]=S;
			EntropyEnthalpy[1]=H;
		}
	}
	return;
}

__device__ void CBI(int i,int j,double* EntropyEnthalpy,int traceback,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double parameter[])
{
	int d,ii,jj;

	for(d=j-i-3;d>=4&&d>=j-i-32;--d)
		for(ii=i+1;ii<j-d&&ii<=Initint[0];++ii)
		{
			jj=d+ii;
			if(traceback==0)
			{
				EntropyEnthalpy[0]=-1.0;
				EntropyEnthalpy[1]=1.0*INFINITY;
			}
			if(fabs(enthalpyDPT[(ii-1)*Initint[2]+jj-1])<999999999)
			{
				calc_bulge_internal2(i,j,ii,jj,EntropyEnthalpy,traceback,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);
				if(fabs(EntropyEnthalpy[1])<999999999)
				{
					if(EntropyEnthalpy[0] <-2500.0)
					{
						EntropyEnthalpy[0]=-3224.0;
						EntropyEnthalpy[1]=0.0;
					}
					if(traceback==0)
					{
						enthalpyDPT[(i-1)*Initint[2]+j-1]=EntropyEnthalpy[1];
						entropyDPT[(i-1)*Initint[2]+j-1]=EntropyEnthalpy[0];
					}
				}
			}
		}
	return;
}

__device__ int find_pos(char *ref,int ref_start,char *source,int start,int length,int num)
{
	int flag,i,j;

	for(i=0;i<num;i++)
	{
		flag=0;
		for(j=0;j<length;j++)
		{
			if(ref[ref_start+j]!=source[start+i*length+j])
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

__device__ void calc_hairpin(int i,int j,double *EntropyEnthalpy,int traceback,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],double parameter[],char *d_Pchar,int *d_NumL)
{
	int pos,loopSize=j-i-1;
	double T1,T2;
	
	if(loopSize < 3)
	{
		EntropyEnthalpy[0]=-1.0;
		EntropyEnthalpy[1]=1.0*INFINITY;
		return;
	}
	if(i<=Initint[0]&&Initint[1]<j)
	{
		EntropyEnthalpy[0]=-1.0;
		EntropyEnthalpy[1]=1.0*INFINITY;
		return;
	}
	else if(i>Initint[1])
	{
		i-= Initint[0];
		j-= Initint[1];
	}
	if(loopSize<=30)
	{
		EntropyEnthalpy[1]=parameter[3090+loopSize-1];
		EntropyEnthalpy[0]=parameter[3000+loopSize-1];
	}
	else
	{
		EntropyEnthalpy[1]=parameter[3090+29];
		EntropyEnthalpy[0]=parameter[3000+29];
	}

	if(loopSize>3) // for loops 4 bp and more in length, terminal mm are accounted
	{
		EntropyEnthalpy[1]+=parameter[5055+numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
		EntropyEnthalpy[0]+=parameter[4430+numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
	}
	else if(loopSize == 3) // for loops 3 bp in length at-penalty is considered
	{
		EntropyEnthalpy[1]+=parameter[5705+numSeq1[i]*5+numSeq1[j]];
		EntropyEnthalpy[0]+=parameter[5680+numSeq1[i]*5+numSeq1[j]];
	}

	if(loopSize==3) // closing AT-penalty (+), triloop bonus, hairpin of 3 (+) 
	{
		pos=find_pos(numSeq1,i,d_Pchar,5*d_NumL[0],5,d_NumL[0]);
		if(pos!=-1)
			EntropyEnthalpy[1]+=parameter[5730+d_NumL[0]+pos];

		pos=find_pos(numSeq1,i,d_Pchar,0,5,d_NumL[0]);
		if(pos!=-1)
			EntropyEnthalpy[0]+=parameter[5730+pos];
	}
	else if (loopSize == 4) // terminal mismatch, tetraloop bonus, hairpin of 4
	{
		pos=find_pos(numSeq1,i,d_Pchar,10*d_NumL[0]+6*d_NumL[1],6,d_NumL[1]);
		if(pos!=-1)
			EntropyEnthalpy[1]+=parameter[5730+2*d_NumL[0]+d_NumL[1]+pos];

		pos=find_pos(numSeq1,i,d_Pchar,10*d_NumL[0],6,d_NumL[1]);
		if(pos!=-1)
			EntropyEnthalpy[0]+=parameter[5730+2*d_NumL[0]+pos];
	}
	if(fabs(EntropyEnthalpy[1])>999999999)
	{
		EntropyEnthalpy[1] =1.0*INFINITY;
		EntropyEnthalpy[0] = -1.0;
	}
	T1 = (EntropyEnthalpy[1] +Initdouble[0]) / ((EntropyEnthalpy[0] +Initdouble[1]+ Initdouble[2]));
	T2 = (enthalpyDPT[(i-1)*Initint[2]+j-1] +Initdouble[0]) / ((entropyDPT[(i-1)*Initint[2]+j-1]) +Initdouble[1]+ Initdouble[2]);
	if(T1 < T2 && traceback == 0)
	{
		EntropyEnthalpy[0] =entropyDPT[(i-1)*Initint[2]+j-1];
		EntropyEnthalpy[1] =enthalpyDPT[(i-1)*Initint[2]+j-1];
	}
	return;
}

__device__ void fillMatrix2(double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double *parameter,char *d_Pchar,int *d_NumL)
{
	int i, j;
	double SH[2];

	for (j = 2; j <= Initint[1]; ++j)
		for (i = j - 3 - 1; i >= 1; --i)
		{
			if (fabs(enthalpyDPT[(i-1)*Initint[2]+j-1])<999999999)
			{
				SH[0] = -1.0;
				SH[1] =1.0*INFINITY;
				maxTM2(i,j,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);
				CBI(i,j,SH,0,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);

				SH[0] = -1.0;
				SH[1] =1.0*INFINITY;
				calc_hairpin(i, j, SH, 0,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,parameter,d_Pchar,d_NumL);
				if(fabs(SH[1])<999999999)
				{
					if(SH[0] <-2500.0) /* to not give dH any value if dS is unreasonable */
					{
						SH[0] =-3224.0;
						SH[1] = 0.0;
					}
					entropyDPT[(i-1)*Initint[2]+j-1]= SH[0];
					enthalpyDPT[(i-1)*Initint[2]+j-1]= SH[1];
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

__device__ double Sd5(int i,int j,char numSeq1[],double parameter[])
{
	return parameter[2750+numSeq1[i]*25+numSeq1[j]*5+numSeq1[j-1]];
}

__device__ double Hd5(int i,int j,char numSeq1[],double parameter[])
{
	return parameter[2875+numSeq1[i]*25+numSeq1[j]*5+numSeq1[j-1]];
}

__device__ double Sd3(int i,int j,char numSeq1[],double parameter[])
{
	return parameter[2500+numSeq1[i]*25+numSeq1[i+1]*5+numSeq1[j]];
}

__device__ double Hd3(int i,int j,char numSeq1[],double parameter[])
{
	return parameter[2625+numSeq1[i]*25+numSeq1[i+1]*5+numSeq1[j]];
}

__device__ double Ststack(int i,int j,char numSeq1[],double parameter[])
{
	return parameter[4430+numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
}

__device__ double Htstack(int i,int j,char numSeq1[],double parameter[])
{
	return parameter[5055+numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
}

__device__ double END5_1(int i,int hs,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[],double parameter[])
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	max_tm=-1.0*INFINITY;
	H_max=1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-5;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+parameter[5705+numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1];
			S=send5[k]+parameter[5680+numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1];
			if(fabs(H)>999999999||H>0||S>0)  // H and S must be greater than 0 to avoid BS
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=parameter[5705+numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1];
			S=parameter[5680+numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
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

__device__ double END5_2(int i,int hs,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[],double parameter[])
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=1.0*INFINITY;
	max_tm=-1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-6;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+parameter[5705+numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-1];
			S=send5[k]+parameter[5680+numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-1];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=parameter[5705+numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-1];
			S=parameter[5680+numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-1];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
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

__device__ double END5_3(int i,int hs,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[],double parameter[])
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=1.0*INFINITY;
	max_tm=-1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-6;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+parameter[5705+numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,numSeq1,parameter)+enthalpyDPT[k*Initint[2]+i-2];
			S=send5[k]+parameter[5680+numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,numSeq1,parameter)+entropyDPT[k*Initint[2]+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=parameter[5705+numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,numSeq1,parameter)+enthalpyDPT[k*Initint[2]+i-2];
			S=parameter[5680+numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,numSeq1,parameter)+entropyDPT[k*Initint[2]+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
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

__device__ double END5_4(int i,int hs,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[],double parameter[])
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=1.0*INFINITY;
	max_tm=-1.0*INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-7;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+parameter[5705+numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-2];
			S=send5[k]+parameter[5680+numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=parameter[5705+numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-2];
			S=parameter[5680+numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-2];
			if(fabs(H)>999999999||H>0||S>0)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
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

__device__ void calc_terminal_bp(double temp,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[],double parameter[])
{
	int i,max;
	double T1,T2,T3,T4,T5,G,end5_11,end5_12,end5_21,end5_22,end5_31,end5_32,end5_41,end5_42;
	
	send5[0]=send5[1]= -1.0;
	hend5[0]=hend5[1]=1.0*INFINITY;

	for(i=2;i<=Initint[0];i++)
	{
		send5[i]=-3224.0;
		hend5[i]=0;
	}

// adding terminal penalties to 3' end and to 5' end 
	for(i=2;i<=Initint[0];++i)
	{
		max=0;
		T1=(hend5[i-1]+Initdouble[0])/(send5[i-1]+Initdouble[1]+Initdouble[2]);
		end5_11=END5_1(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		end5_12=END5_1(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		T2=(end5_11+Initdouble[0])/(end5_12+Initdouble[1]+Initdouble[2]);
		end5_21=END5_2(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		end5_22=END5_2(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		T3=(end5_21+Initdouble[0])/(end5_22+Initdouble[1]+Initdouble[2]);
		end5_31=END5_3(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		end5_32=END5_3(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		T4=(end5_31+Initdouble[0])/(end5_32+Initdouble[1]+Initdouble[2]);
		end5_41=END5_4(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		end5_42=END5_4(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		T5=(end5_41+Initdouble[0])/(end5_42+Initdouble[1]+Initdouble[2]);

		max=max5(T1,T2,T3,T4,T5);
		switch(max)
		{
			case 1:
				send5[i]=send5[i-1];
				hend5[i]=hend5[i-1];
				break;
			case 2:
				G=end5_11-temp*end5_12;
				if(G<0.0)
				{
					send5[i]=end5_12;
					hend5[i]=end5_11;
				}
				else
				{
					send5[i]=send5[i-1];
					hend5[i]=hend5[i-1];
				}
				break;
			case 3:
				G=end5_21-temp*end5_22;
				if(G<0.0)
				{
					send5[i]=end5_22;
					hend5[i]=end5_21;
				}
				else
				{
					send5[i]=send5[i-1];
					hend5[i]=hend5[i-1];
				}
				break;
			case 4:
				G=end5_31-temp*end5_32;
				if(G<0.0)
				{
					send5[i]=end5_32;
					hend5[i]=end5_31;
				}
				else
				{
					send5[i]=send5[i-1];
					hend5[i]=hend5[i-1];
				}
				break;
			case 5:
				G=end5_41-temp*end5_42;
				if(G<0.0)
				{
					send5[i]=end5_42;
					hend5[i]=end5_41;
				}
				else
				{
					send5[i]=send5[i-1];
					hend5[i]=hend5[i-1];
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

__device__ void tracebacku(int bp[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[],char numSeq2[],double *parameter,char *d_Pchar,int *d_NumL)
{
	int i,j,store[50],total,now,ii,jj,k,d,done;
	double SH1[2],SH2[2],EntropyEnthalpy[2];

        total=newpush(store,Initint[0],0,1,0,0);
        now=0;
        while(now<total)
        {
                i=store[3*now]; // top->i;
                j=store[3*now+1]; // top->j;
                if(store[now*3+2]==1)
                {
                        while(equal(send5[i],send5[i-1])&&equal(hend5[i],hend5[i-1])) // if previous structure is the same as this one
                                --i;
                        if(i==0)
                                continue;
                        if(equal(send5[i],END5_1(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter))&&equal(hend5[i],END5_1(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter)))
                        {
                                for(k=0;k<=i-5;++k)
                                        if(equal(send5[i],parameter[5680+numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1])&&equal(hend5[i],parameter[5705+numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);                    
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+parameter[5680+numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1])&&equal(hend5[i],hend5[k]+parameter[5705+numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_2(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter))&&equal(hend5[i],END5_2(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(send5[i],parameter[5680+numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-1])&&equal(hend5[i],parameter[5705+numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+parameter[5680+numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-1])&&equal(hend5[i],hend5[k]+parameter[5705+numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_3(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter))&&equal(hend5[i],END5_3(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(send5[i],parameter[5680+numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,numSeq1,parameter)+entropyDPT[k*Initint[2]+i-2])&&equal(hend5[i],parameter[5705+numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,numSeq1,parameter)+enthalpyDPT[k*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+parameter[5680+numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,numSeq1,parameter)+entropyDPT[k*Initint[2]+i-2])&&equal(hend5[i],hend5[k]+parameter[5705+numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,numSeq1,parameter)+enthalpyDPT[k*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_4(i,2,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter))&&equal(hend5[i],END5_4(i,1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter)))
                        {
                                for (k=0;k<=i-7;++k)
                                        if(equal(send5[i],parameter[5680+numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-2])&&equal(hend5[i],parameter[5705+numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+2,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+parameter[5680+numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,numSeq1,parameter)+entropyDPT[(k+1)*Initint[2]+i-2])&&equal(hend5[i],hend5[k]+parameter[5705+numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,numSeq1,parameter)+enthalpyDPT[(k+1)*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+2,i-1,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                }
                else if(store[3*now+2]==0)
                {
                        bp[i-1]=j;
                        bp[j-1]=i;
                        SH1[0]=-1.0;
                        SH1[1]=1.0*INFINITY;
                        calc_hairpin(i,j,SH1,1,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,parameter,d_Pchar,d_NumL);

                        SH2[0]=-1.0;
                        SH2[1]=1.0*INFINITY;
                        CBI(i,j,SH2,2,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);

                        if (equal(entropyDPT[(i-1)*Initint[2]+j-1],Ss(i,j,2,Initint,numSeq1,numSeq2,parameter)+entropyDPT[i*Initint[2]+j-2])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],Hs(i,j,2,Initint,numSeq1,numSeq2,parameter)+enthalpyDPT[i*Initint[2]+j-2]))
                                total=newpush(store,i+1,j-1,0,total,now+1);
                        else if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH1[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH1[1]));
                        else if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH2[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH2[1]))
                        {
                                for (done=0,d=j-i-3;d>=4&&d>=j-i-32&&!done;--d)
                                        for (ii=i+1;ii<j-d;++ii)
                                        {
                                                jj=d+ii;
                                                EntropyEnthalpy[0]=-1.0;
                                                EntropyEnthalpy[1]=1.0*INFINITY;
                                                calc_bulge_internal2(i,j,ii,jj,EntropyEnthalpy,1,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);

                                                if (equal(entropyDPT[(i-1)*Initint[2]+j-1],EntropyEnthalpy[0]+entropyDPT[(ii-1)*Initint[2]+jj-1])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],EntropyEnthalpy[1]+enthalpyDPT[(ii-1)*Initint[2]+jj-1]))
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

__device__ double drawHairpin(int bp[],double mh,double ms,int Initint[])
{
        int i,N;

        N=0;
        if(fabs(ms)>999999999||fabs(mh)>999999999)
        {
		return 0.0;
        }
        else
        {
		for(i=1;i<Initint[0];++i)
		{
			if(bp[i-1]>0)
				N++;
                }
                return mh/(ms+(((N/2)-1)*-0.51986))-273.15;
        }
}

__device__ void initMatrix(int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
	int i,j;

	for(i=1;i<=Initint[0];++i)
	{
		for(j=1;j<=Initint[1];++j)
		{
			if(numSeq1[i]+numSeq2[j]!=3)
			{
				enthalpyDPT[(i-1)*Initint[2]+j-1]=1.0*INFINITY;
				entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
			}
			else
			{
				enthalpyDPT[(i-1)*Initint[2]+j-1]=0.0;
				entropyDPT[(i-1)*Initint[2]+j-1]=-3224.0;
			}
		}
	}
}

__device__ void LSH(int i,int j,double *EntropyEnthalpy,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double parameter[])
{
	double S1,H1,T1,S2,H2,T2;

	if(numSeq1[i]+numSeq2[j]!=3)
	{
		entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
		enthalpyDPT[(i-1)*Initint[2]+j-1]=1.0*INFINITY;
		return;
	}

	S1=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[4430+numSeq2[j]*125+numSeq2[j-1]*25+numSeq1[i]*5+numSeq1[i-1]];
	H1=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[5055+numSeq2[j]*125+numSeq2[j-1]*25+numSeq1[i]*5+numSeq1[i-1]];
	if(fabs(H1)>999999999)
	{
		H1=1.0*INFINITY;
		S1=-1.0;
	}
// If there is two dangling ends at the same end of duplex
	if(fabs(parameter[2625+numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]])<999999999&&fabs(parameter[2875+numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]])<999999999)
	{
		S2=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[2500+numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]]+parameter[2750+numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		H2=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[2625+numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]]+parameter[2875+numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
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
	else if(fabs(parameter[2625+numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]])<999999999)
	{
		S2=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[2500+numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]];
		H2=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[2625+numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
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
	else if(fabs(parameter[2875+numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]])<999999999)
	{
		S2=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[2750+numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		H2=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[2875+numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
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

	S2=parameter[5680+numSeq1[i]*5+numSeq2[j]];
	H2=parameter[5705+numSeq1[i]*5+numSeq2[j]];
	T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
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

__device__ void maxTM(int i,int j,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double parameter[])
{
	double T0,T1,S0,S1,H0,H1;

	S0=entropyDPT[(i-1)*Initint[2]+j-1];
	H0=enthalpyDPT[(i-1)*Initint[2]+j-1];
	T0=(H0+Initdouble[0])/(S0+Initdouble[1]+Initdouble[2]); // at current position 
	if(fabs(enthalpyDPT[(i-2)*Initint[2]+j-2])<999999999&&fabs(Hs(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter))<999999999)
	{
		S1=(entropyDPT[(i-2)*Initint[2]+j-2]+Ss(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter));
		H1=(enthalpyDPT[(i-2)*Initint[2]+j-2]+Hs(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter));
	}
	else
	{
		S1=-1.0;
		H1=1.0*INFINITY;
	}
	T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);

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
		entropyDPT[(i-1)*Initint[2]+j-1]=S1;
		enthalpyDPT[(i-1)*Initint[2]+j-1]=H1;
	}
	else if(T0>=T1)
	{
		entropyDPT[(i-1)*Initint[2]+j-1]=S0;
		enthalpyDPT[(i-1)*Initint[2]+j-1]=H0;
	}
}

__device__ void calc_bulge_internal(int i,int j,int ii,int jj,double* EntropyEnthalpy,int traceback,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double parameter[])
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
				H=parameter[3150+loopSize]+parameter[625+numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
				S=parameter[3060+loopSize]+parameter[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
			}
			H+=enthalpyDPT[(i-1)*Initint[2]+j-1];
			S+=entropyDPT[(i-1)*Initint[2]+j-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}

			T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
			T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/((entropyDPT[(ii-1)*Initint[2]+jj-1])+Initdouble[1]+Initdouble[2]);
			if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
		else // we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30
		{
			H=parameter[3150+loopSize]+parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[5705+numSeq1[ii]*5+numSeq2[jj]];
			H+=enthalpyDPT[(i-1)*Initint[2]+j-1];

			S=parameter[3060+loopSize]+parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[5680+numSeq1[ii]*5+numSeq2[jj]];
			S+=entropyDPT[(i-1)*Initint[2]+j-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
			T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/(entropyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[1]+Initdouble[2]);
			if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
	}
	else if(loopSize1==1&&loopSize2==1)
	{
		S=parameter[1250+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+parameter[1250+numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		S+=entropyDPT[(i-1)*Initint[2]+j-1];

		H=parameter[1875+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+parameter[1875+numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		H+=enthalpyDPT[(i-1)*Initint[2]+j-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
		T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/(entropyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[1]+Initdouble[2]);
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
		H=parameter[3120+loopSize]+parameter[3805+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+parameter[3805+numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		H+=enthalpyDPT[(i-1)*Initint[2]+j-1];

		S=parameter[3030+loopSize]+parameter[3180+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+parameter[3180+numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
		S+=entropyDPT[(i-1)*Initint[2]+j-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
		T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/((entropyDPT[(ii-1)*Initint[2]+jj-1])+Initdouble[1]+Initdouble[2]);
		if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
		{
			EntropyEnthalpy[0]=S;
			EntropyEnthalpy[1]=H;
		}
	}
	return;
}

__device__ void fillMatrix(double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double *parameter)
{
	int d,i,j,ii,jj;
	double SH[2];

	for(i=1;i<=Initint[0];++i)
	{
		for(j=1;j<=Initint[1];++j)
		{
			if(fabs(enthalpyDPT[(i-1)*Initint[2]+j-1])<999999999)
			{
				SH[0]=-1.0;
				SH[1]=1.0*INFINITY;
				LSH(i,j,SH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);

				if(fabs(SH[1])<999999999)
				{
					entropyDPT[(i-1)*Initint[2]+j-1]=SH[0];
					enthalpyDPT[(i-1)*Initint[2]+j-1]=SH[1];
				}
				if(i>1&&j>1)
				{
					maxTM(i,j,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);
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
							if(fabs(enthalpyDPT[(ii-1)*Initint[2]+jj-1])<999999999)
							{
								SH[0]=-1.0;
								SH[1]=1.0*INFINITY;
								calc_bulge_internal(ii,jj,i,j,SH,0,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);

								if(SH[0]<-2500.0)
								{
									SH[0] =-3224.0;
									SH[1] = 0.0;
								}
								if(fabs(SH[1])<999999999)
								{
									enthalpyDPT[(i-1)*Initint[2]+j-1]=SH[1];
									entropyDPT[(i-1)*Initint[2]+j-1]=SH[0];
								}
							}
						}
					}
				} // if 
			}
		} // for 
	} //for
}

__device__ void RSH(int i,int j,double EntropyEnthalpy[],double Initdouble[],char numSeq1[],char numSeq2[],double *parameter)
{
	double S1,S2,H1,H2,T1,T2;

	if(numSeq1[i]+numSeq2[j]!=3)
	{
		EntropyEnthalpy[0]=-1.0;
		EntropyEnthalpy[1]=1.0*INFINITY;
		return;
	}
	S1=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[4430+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
	H1=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[5055+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
	if(fabs(H1)>999999999)
	{
		H1=1.0*INFINITY;
		S1=-1.0;
	}
	if(fabs(parameter[2625+numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]])<999999999&&fabs(parameter[2875+numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]])<999999999)
	{
		S2=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[2500+numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]]+parameter[2750+numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
		H2=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[2625+numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]]+parameter[2875+numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
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

	if(fabs(parameter[2625+numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]])<999999999)
	{
		S2=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[2500+numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]];
		H2=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[2625+numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
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

	if(fabs(parameter[2875+numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]])<999999999)
	{
		S2=parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[2750+numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
		H2=parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[2875+numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
		if(fabs(H2)>999999999)
		{
			H2=1.0*INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(fabs(H1)<999999999)
		{
			T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
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
	S2=parameter[5680+numSeq1[i]*5+numSeq2[j]];
	H2=parameter[5705+numSeq1[i]*5+numSeq2[j]];
	T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
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

__device__ void traceback(int i,int j,int* ps1,int* ps2,double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[],double *parameter)
{
	int d,ii,jj,done;
	double SH[2];

	ps1[i-1]=j;
	ps2[j-1]=i;
	while(1)
	{
		SH[0]=-1.0;
		SH[1]=1.0*INFINITY;
		LSH(i,j,SH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);
		if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH[1]))
			break;

		done = 0;
		if(i>1&&j>1&&equal(entropyDPT[(i-1)*Initint[2]+j-1],Ss(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter)+entropyDPT[(i-2)*Initint[2]+j-2]))
		{
			i=i-1;
			j=j-1;
			ps1[i-1]=j;
			ps2[j-1]=i;
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
				calc_bulge_internal(ii,jj,i,j,SH,1,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);
				if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH[1]))
				{
					i=ii;
					j=jj;
					ps1[i-1]=j;
					ps2[j-1]=i;
					done=1;
					break;
				}
			}
		}
	}
}

__device__ double drawDimer(int *ps1,int *ps2,double H,double S,double Initdouble[],int Initint[])
{
        int i,N;

        if(fabs(Initdouble[3])>999999999)
                return (double)0.0;
        else
        {
                N=0;
                for(i=0;i<Initint[0];i++)
                {
                        if(ps1[i]>0)
                                ++N;
                }
                for(i=0;i<Initint[1];i++)
                {
                        if(ps2[i]>0)
                                ++N;
                }
                N=(N/2)-1;
                return (double)(H/(S+(N*-0.51986)+Initdouble[2])-273.15);
        }
}

__device__ int symmetry_thermo(char seq[])
{
	int i = 0;
	int seq_len=seq_length(seq);
	if(seq_len%2==1)
		return 0;

	while(i<seq_len/2)
	{
		if((seq[i]=='A'&&seq[seq_len-1-i]!='T')||(seq[i]=='T'&&seq[seq_len-1-i]!='A')||(seq[seq_len-1-i]=='A'&&seq[i]!='T')||(seq[seq_len-1-i]=='T'&&seq[i]!='A'))
			return 0;
		if((seq[i]=='C'&&seq[seq_len-1-i]!='G')||(seq[i]=='G'&&seq[seq_len-1-i]!='C')||(seq[seq_len-1-i]=='C'&&seq[i]!='G')||(seq[seq_len-1-i]=='G'&&seq[i]!='C'))
			return 0;
		i++;
	}
	return 1;
}

__device__ double thal(char oligo_f[],char oligo_r[],int type,double *parameter,char *d_Pchar,int *d_NumL)
{
	double SH[2],Initdouble[4];//0 is dplx_init_H, 1 is dplx_init_S, 2 is RC, 3 is SHleft
	int Initint[5]; //0 is len1, 1 is len2, 2 is len3, 3 is bestI, 4 is bestJ
	int i, j;
	double T1,enthalpyDPT[625],entropyDPT[625],send5[26],hend5[26],result_TH;
	int ps1[25],ps2[25];
	char numSeq1[27],numSeq2[27];
	double mh, ms;

/*** INIT values for unimolecular and bimolecular structures ***/
	if (type==4) /* unimolecular folding */
	{
		Initdouble[0]= 0.0;
		Initdouble[1] = -0.00000000001;
		Initdouble[2]=0;
	}
	else /* hybridization of two oligos */
	{
		Initdouble[0]= 200;
		Initdouble[1]= -5.7;
		if(symmetry_thermo(oligo_f) && symmetry_thermo(oligo_r))
			Initdouble[2]=1.9872* log(38/1000000000.0);
		else
			Initdouble[2]=1.9872* log(38/4000000000.0);
	}
/* convert nucleotides to numbers */
	if(type==1 || type==2)
	{
		Initint[0]=seq_length(oligo_f);
		Initint[1]=seq_length(oligo_r);
	 	for(i=1;i<=Initint[0];++i)
			numSeq1[i]=str2int(oligo_f[i-1]);
		for(i=1;i<=Initint[1];++i)
			numSeq2[i]=str2int(oligo_r[Initint[1]-i]);
	}
	else if(type==3)
	{
		Initint[0]=seq_length(oligo_r);
		Initint[1]=seq_length(oligo_f);
		for(i=1;i<=Initint[0];++i)
			numSeq1[i]=str2int(oligo_r[i-1]);
		for(i=1;i<=Initint[1];++i)
			numSeq2[i]=str2int(oligo_f[Initint[1]-i]);
	}
	else
	{
		Initint[0]=seq_length(oligo_f);
                Initint[1]=seq_length(oligo_r);
		Initint[2]=Initint[1]-1;
                for(i=1;i<=Initint[0];++i)      
                        numSeq1[i]=str2int(oligo_f[i-1]);   
                for(i=1;i<=Initint[1];++i)      
                        numSeq2[i]=str2int(oligo_r[i-1]);
	}
	numSeq1[0]=numSeq1[Initint[0]+1]=numSeq2[0]=numSeq2[Initint[1]+1]=4; /* mark as N-s */

	result_TH=0;
	if (type==4) /* calculate structure of monomer */
	{
		initMatrix2(Initint,enthalpyDPT,entropyDPT,numSeq1);
		fillMatrix2(Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter,d_Pchar,d_NumL);
		calc_terminal_bp(310.15,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,parameter);
		mh=hend5[Initint[0]];
		ms=send5[Initint[0]];
		for (i=0;i<Initint[0];i++)
			ps1[i]=0;
		if(fabs(mh)<999999999)
		{
			tracebacku(ps1,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,numSeq2,parameter,d_Pchar,d_NumL);
			result_TH=drawHairpin(ps1,mh,ms,Initint);
			result_TH=(int)(result_TH*100+0.5)/100.0;
		}
	}
	else  /* Hybridization of two moleculs */
	{
		Initint[2]=Initint[1];
		initMatrix(Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
		fillMatrix(Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);

		Initdouble[3]=-1.0*INFINITY;
	/* calculate terminal basepairs */
		Initint[3]=Initint[4]=0;
		if(type==1)
			for (i=1;i<=Initint[0];i++)
			{
				for (j=1;j<=Initint[1];j++)
				{
					RSH(i,j,SH,Initdouble,numSeq1,numSeq2,parameter);
					SH[0]=SH[0]+0.000001; /* this adding is done for compiler, optimization -O2 vs -O0 */
					SH[1]=SH[1]+0.000001;
					T1=((enthalpyDPT[(i-1)*Initint[2]+j-1]+ SH[1] +Initdouble[0]) / ((entropyDPT[(i-1)*Initint[2]+j-1]) + SH[0] +Initdouble[1] + Initdouble[2])) -273.15;
					if(T1>Initdouble[3]&&((entropyDPT[(i-1)*Initint[2]+j-1]+SH[0])<0&&(SH[1]+enthalpyDPT[(i-1)*Initint[2]+j-1])<0))
					{
						Initdouble[3]=T1;
						Initint[3]=i;
						Initint[4]=j;
					}
				}
			}
		if(type==2||type==3)
		{
		 //THAL_END1
			Initint[4]=0;
			Initint[3]=Initint[0];
			i=Initint[0];
			Initdouble[3]=-1.0*INFINITY;
			for (j=1;j<=Initint[1];++j)
			{
				RSH(i,j,SH,Initdouble,numSeq1,numSeq2,parameter);
				SH[0]=SH[0]+0.000001; // this adding is done for compiler, optimization -O2 vs -O0,that compiler could understand that SH is changed in this cycle 
				SH[1]=SH[1]+0.000001;
				T1=((enthalpyDPT[(i-1)*Initint[2]+j-1]+SH[1]+Initdouble[0])/((entropyDPT[(i-1)*Initint[2]+j-1])+SH[0]+Initdouble[1]+Initdouble[2]))-273.15;
				if (T1>Initdouble[3]&&((SH[0]+entropyDPT[(i-1)*Initint[2]+j-1])<0&&(SH[1]+enthalpyDPT[(i-1)*Initint[2]+j-1])<0))
				{
					Initdouble[3]=T1;
					Initint[4]=j;
				}
			}
		}
		if(fabs(Initdouble[3])>999999999)
			Initint[3]=Initint[4]=1;
		RSH(Initint[3],Initint[4],SH,Initdouble,numSeq1,numSeq2,parameter);
	 // tracebacking 
		for (i=0;i<Initint[0];++i)
			ps1[i]=0;
		for (j=0;j<Initint[1];++j)
			ps2[j] = 0;
		if(fabs(enthalpyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1])<999999999)
		{
			traceback(Initint[3],Initint[4],ps1,ps2,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2,parameter);
			result_TH=drawDimer(ps1,ps2,(enthalpyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1]+SH[1]+Initdouble[0]),(entropyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1]+SH[0]+Initdouble[1]),Initdouble,Initint);
			result_TH=(int)(result_TH*100+0.5)/100.0;
		}
	}
        return result_TH;
}

///function in gpu, generate a read; int length: the length of reads
__device__ void generate(char *d_seq,char seq[],int pos,int length)
{
	int i;
	for(i=0;i<length;i++)
	{
		seq[i]=d_seq[pos+i];
	}
	seq[i]='\0';
}

///function in gpu, check the GC-content; int length: the length of read
__device__ int gc(char seq[],int length)
{
	int i,number;
	float gc;

	number=0;
	for(i=0;i<length;i++)
	{
		if(seq[i]=='C')
		{
			number++;
			continue;
		}
	
		if(seq[i]=='G')
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

//function in gpu, caculate tm
__device__ int tm(char seq[],float *d_deltah,float *d_deltas,int length,float max_tm,float min_tm)
{
	int i,pos;
	float deltah,deltas,result;

	deltah=0;
	deltas=0;
	for(i=0;i<length-1;i++)
	{
		pos=translate(seq[i]);
		pos=pos*4+translate(seq[i+1]);
		deltah+=d_deltah[pos];
		deltas+=d_deltas[pos];
	}

	deltah=(-1.0)*deltah;
	deltas=(-1.0)*deltas;
	if((seq[0]=='A')||(seq[0]=='T'))
	{
		deltah+=2.3;
		deltas+=4.1;
	}
	else
	{
		deltah+=0.1;
		deltas-=2.8;
	}
        if((seq[length-1]=='A')||(seq[length-1]=='T'))
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
__device__ int stability(char seq[],float *d_stab,int length,int strand)
{
	int i,pos;
	
	pos=0;
	for(i=0;i<6;i++)
	{
		if(strand==0)
		{
			pos=pos*4+translate(seq[i]);
		}
		else
		{
			pos=pos*4+translate(seq[i+length-6]);
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
                if(strand==1)
                {
                        pos=pos*4+translate(seq[i]);
                }
                else
                {
                        pos=pos*4+translate(seq[i+length-6]);
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

//function in gpu, reverse the strand,+ to - strand
__device__ void reverse(char seq[],char rev[],int length)
{
	int i;
	
	for(i=0;i<length;i++)
	{
		if(seq[length-1-i]=='A')
		{
			rev[i]='T';
			continue;
		}
                if(seq[length-1-i]=='T')
                {
                        rev[i]='A';
                        continue;
                }
                if(seq[length-1-i]=='C')
                {
                        rev[i]='G';
                        continue;
                }
		rev[i]='C';
	}
	rev[i]='\0';
}

__device__ int check_long_ploy(char primer[],int length)
{
        int i,same;
        char ref;

        same=1;
        ref=primer[0];
        for(i=1;i<length;i++)
        {
                if(primer[i]==ref)
                        same++;
                else
                {
                        if(same>=6)
                                return 0;
                        same=1;
                        ref=primer[i];
                }
        }
        if(same>=6)
                return 0;
        return 1;
}

///function: int length: the length of genome
__global__ void candidate_primer(char *d_seq,int *d_pos,int *d_len,int *d_rev_len,float *d_stab,float *d_deltah,float *d_deltas,int strand,float max_tm,float min_tm,int length,int check_flag,double *parameter,char *d_Pchar,int *d_NumL)
{
	int id,i,circle,check,plus,minus;
	char primer[30],rev[30];

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

			generate(d_seq,primer,circle,i);
			check=gc(primer,i);
			if(check==0)
				continue;

			check=check_long_ploy(primer,i);
			if(check==0)
                                continue;

			check=tm(primer,d_deltah,d_deltas,i,max_tm,min_tm);
			if(check==0)
				continue;

                        check=stability(primer,d_stab,i,strand);
                        if(check==1)     //+ strand
				plus=1;
			else
				plus=0;
			
		//secondary structure
			if(check_flag&&plus)
			{
				if(thal(primer,primer,1,parameter,d_Pchar,d_NumL)>min_tm-10)
					plus=0;	
			}
			if(check_flag&&plus)
                        {
                                if(thal(primer,primer,2,parameter,d_Pchar,d_NumL)>min_tm-10)  
                                        plus=0;
                        }
			if(check_flag&&plus)
                        {                
                                if(thal(primer,primer,4,parameter,d_Pchar,d_NumL)>min_tm-10)
                                        plus=0;         
                        }
			if(plus)
                                d_len[circle*8+i-18]=1;
	
			reverse(primer,rev,i);  //generate - strand
			check=stability(rev,d_stab,i,strand);
			if(check==1)
				minus=1;
			else
				minus=0;
		//secondary structure      
                        if(check_flag&&minus)
                        {                
                                if(thal(rev,rev,1,parameter,d_Pchar,d_NumL)>min_tm-10)
                                        minus=0;         
                        }           
                        if(check_flag&&minus)
                        {
                                if(thal(rev,rev,2,parameter,d_Pchar,d_NumL)>min_tm-10)
                                        minus=0;
                        }                
                        if(check_flag&&minus)
                        {
                                if(thal(rev,rev,4,parameter,d_Pchar,d_NumL)>min_tm-10)
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
	double *H_parameter,*parameter;
	int *pos,*d_pos,*len,*d_len,length,flag[10],i,*rev_len,*d_rev_len,num_outer,num_inner,num_loop,NumL[2],*d_NumL;
	float deltah[16],deltas[16],stab[4096],*d_deltah,*d_deltas,*d_stab,temp1,temp2;
	char *seq,*d_seq,*store_path,*prefix,*stab_path,*tm_path,*curren_path,*input,*outer,*inner,*loop,*par_path,*temp,*Pchar,*d_Pchar;
	FILE *fp;
	time_t start,end;
        struct stat statbuf;
//flag: 0:input; 1: out_prefix; 2: dir; 3: stab; 4: tm; 5: high; 6: low; 7: loop; 8: secondary structure; 9: path for secondary structure

	start=time(NULL);
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
	        cudaMalloc((void **)&d_Pchar,10*NumL[0]+12*NumL[1]);
	        cudaMemset(d_Pchar,'\0',10*NumL[0]+12*NumL[1]);
		cudaMalloc((void **)&parameter,(5730+2*NumL[0]+2*NumL[1])*sizeof(double));
	        cudaMemset(parameter,'\0',(5730+2*NumL[0]+2*NumL[1])*sizeof(double));
		cudaMalloc((void **)&d_NumL,2*sizeof(int));
	        cudaMemset(d_NumL,'\0',2*sizeof(int));

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

	cudaMalloc((void **)&d_deltah,16*sizeof(float));
	cudaMemset(d_deltah,'\0',16*sizeof(float));
	cudaMalloc((void **)&d_deltas,16*sizeof(float));
	cudaMemset(d_deltas,'\0',16*sizeof(float));
	cudaMalloc((void **)&d_stab,4096*sizeof(float));
	cudaMemset(d_stab,'\0',4096*sizeof(float));

	/////from cpu to gpu
	cudaMemcpy(d_seq,seq,length*sizeof(char),cudaMemcpyHostToDevice);
	cudaMemcpy(d_deltah,deltah,16*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_deltas,deltas,16*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_stab,stab,4096*sizeof(float),cudaMemcpyHostToDevice);

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

//secondary structure
	if(flag[8])
	{
		cudaMemcpy(parameter,H_parameter,(5730+2*NumL[0]+2*NumL[1])*sizeof(double),cudaMemcpyHostToDevice);
        	cudaMemcpy(d_Pchar,Pchar,10*NumL[0]+12*NumL[1],cudaMemcpyHostToDevice);
		cudaMemcpy(d_NumL,NumL,2*sizeof(int),cudaMemcpyHostToDevice);
	}
	end=time(NULL);
	printf("It takes %d seconds to prepare.\n",(int)difftime(end,start));
	start=time(NULL);

	if(flag[5]==1)
        {
		cudaMemset(d_pos,'\0',length*sizeof(int));
		cudaMemset(d_len,'\0',8*length*sizeof(int));
		cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
		candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,61,59,length,flag[8],parameter,d_Pchar,d_NumL);
		cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
        	cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
        	cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_outer=create_file(prefix,outer,seq,pos,len,rev_len,length);

		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,0,66,64,length,flag[8],parameter,d_Pchar,d_NumL);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_inner=create_file(prefix,inner,seq,pos,len,rev_len,length);

                if(flag[7]==1)
		{
			cudaMemset(d_pos,'\0',length*sizeof(int));
                	cudaMemset(d_len,'\0',8*length*sizeof(int));
                	cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                	candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,66,64,length,flag[8],parameter,d_Pchar,d_NumL);
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
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,56,54,length,flag[8],parameter,d_Pchar,d_NumL);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_outer=create_file(prefix,outer,seq,pos,len,rev_len,length);

		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,0,61,59,length,flag[8],parameter,d_Pchar,d_NumL);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_inner=create_file(prefix,inner,seq,pos,len,rev_len,length);
                if(flag[7]==1)
		{
			cudaMemset(d_pos,'\0',length*sizeof(int));
                	cudaMemset(d_len,'\0',8*length*sizeof(int));
                	cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                	candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,61,59,length,flag[8],parameter,d_Pchar,d_NumL);
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
	cudaFree(d_stab);
	cudaFree(d_deltah);
	cudaFree(d_deltas);
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
		cudaFree(parameter);
        	cudaFree(d_Pchar);
        	cudaFree(d_NumL);
	}
	if(flag[8]||flag[9])
		free(par_path);
}
