#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include<sys/stat.h>
#include<cuda.h>
#include<cuda_runtime.h>

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
		i++;
		if((seq[i]=='A'&&seq[seq_len-1-i]!='T')||(seq[i]=='T'&&seq[seq_len-1-i]!='A')||(seq[seq_len-1-i]=='A'&&seq[i]!='T')||(seq[seq_len-1-i]=='T'&&seq[i]!='A'))
			return 0;
		if((seq[i]=='C'&&seq[seq_len-1-i]!='G')||(seq[i]=='G'&&seq[seq_len-1-i]!='C')||(seq[seq_len-1-i]=='C'&&seq[i]!='G')||(seq[seq_len-1-i]=='G'&&seq[i]!='C'))
			return 0;
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
			result_TH=(int)(100*result_TH+0.5)/100.0;
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
			result_TH=(int)(100*result_TH+0.5)/100.0;
		}
	}
        return result_TH;
}

struct Node
{
	int pos;
	int gi;
	int plus;  //as a flag, 1 is OK, 0 is no
	int minus; //as a flag
	struct Node *next;
};

struct Primer
{
	int pos;
	int len;
	int plus;
	int minus;
	int total_common;
	int total_special;
	int total; //common number
	struct Primer *next;
	struct Node *common;
	struct Node *special;
};

struct INFO
{
        char name[301];
        int turn;
        struct INFO *next;
};

int check_add(int F3_pos,int *par,int have)
{
        int i,dis;

        for(i=0;i<have;i++)
        {
                dis=par[i]-F3_pos;
                if(abs(dis)<300)              
                        return 0;
        }
        return 1;        
}

void generate_primer(char *seq,char primer[],int start,int length,int flag)
{
        int i;
        if(flag==0)
        {
                for(i=0;i<length;i++)
                	primer[i]=seq[start+i];
        }
        else
        {
                for(i=0;i<length;i++)
                {
                        if(seq[start+length-1-i]=='A')
                                primer[i]='T';
                        else if(seq[start+length-1-i]=='T')
                                primer[i]='A';
                        else if(seq[start+length-1-i]=='C')
                                primer[i]='G';
                        else
                                primer[i]='C';
                }
        }
	primer[length]='\0';
}

__device__ void get_primer(char *d_seq,char primer[],int start,int length,int flag) //flag=0,plus
{
	int i;

	for(i=0;i<length;i++)
	{
		if(flag==0)
			primer[i]=d_seq[i+start];
		else
		{
			if(d_seq[start+length-1-i]=='A')
				primer[i]='T';
			else if(d_seq[start+length-1-i]=='T')
				primer[i]='A';
			else if(d_seq[start+length-1-i]=='C')
				primer[i]='G';
			else
				primer[i]='C';
		}

	}
	primer[i]='\0';
}

__device__ int check_structure(char *d_seq,int *d_S,int *d_L,int turn[],int high_flag,double *parameter,char *d_Pchar,int *d_NumL)
{
	double TH;
	char P_F3[26],P_F2[26],P_F1c[26],P_B1c[26],P_B2[26],P_B3[26],r_F3[26],r_F2[26],r_F1c[26],r_B1c[26],r_B2[26],r_B3[26],*list[6],*rev[6];
	int i,j;

	get_primer(d_seq,P_F3,d_S[4*turn[0]],d_S[4*turn[0]+1],0);
	get_primer(d_seq,P_F2,d_S[4*turn[1]],d_S[4*turn[1]+1],0);
	get_primer(d_seq,P_F1c,d_L[4*turn[2]],d_L[4*turn[2]+1],1);
	get_primer(d_seq,P_B1c,d_L[4*turn[3]],d_L[4*turn[3]+1],0);
	get_primer(d_seq,P_B2,d_S[4*turn[4]],d_S[4*turn[4]+1],1);
	get_primer(d_seq,P_B3,d_S[4*turn[5]],d_S[4*turn[5]+1],1);

	get_primer(d_seq,r_F3,d_S[4*turn[0]],d_S[4*turn[0]+1],1);
        get_primer(d_seq,r_F2,d_S[4*turn[1]],d_S[4*turn[1]+1],1);
        get_primer(d_seq,r_F1c,d_L[4*turn[2]],d_L[4*turn[2]+1],0);
        get_primer(d_seq,r_B1c,d_L[4*turn[3]],d_L[4*turn[3]+1],1);
        get_primer(d_seq,r_B2,d_S[4*turn[4]],d_S[4*turn[4]+1],0);                               
        get_primer(d_seq,r_B3,d_S[4*turn[5]],d_S[4*turn[5]+1],0);

	list[0]=P_F3;
	list[1]=P_F2;
	list[2]=P_F1c;
	list[3]=P_B1c;
	list[4]=P_B2;
	list[5]=P_B3;
	rev[0]=r_F3;
        rev[1]=r_F2;
        rev[2]=r_F1c;
        rev[3]=r_B1c;
        rev[4]=r_B2;
        rev[5]=r_B3;
	for(i=0;i<5;i++)
	{
		for(j=i+1;j<6;j++)
		{
			TH=thal(list[i],list[j],1,parameter,d_Pchar,d_NumL);
			if(TH>44+5*high_flag)
                                return 0;
			TH=thal(list[i],list[j],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*high_flag)
                                return 0;
			TH=thal(list[i],list[j],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*high_flag)
                                return 0;
			TH=thal(rev[j],rev[i],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*high_flag)
                                return 0;
                        TH=thal(rev[j],rev[i],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*high_flag)
                                return 0;
		}
	}
	return 1;
}

void how_many(struct Primer *head,int common)
{
        struct Primer *p_primer;
        struct Node *p_node;
        int i,num,*list;

	list=(int *)malloc(common*sizeof(int));
        p_primer=head;
        while(p_primer)
        {
                p_node=p_primer->common;
		for(i=0;i<common;i++)
                {
                        list[i]=0;
                }
                i=0;
                while(p_node)
                {
                        i++;
			list[p_node->gi]=1;
                        p_node=p_node->next;
                }
		p_primer->total_common=i;

        //special
                p_node=p_primer->special;
                i=0;
                while(p_node)       
                {
                        i++;
                        p_node=p_node->next;
                }
		p_primer->total_special=i;

		num=0;
                for(i=0;i<common;i++)
                {
                        num=num+list[i];
                }
		p_primer->total=num;
                p_primer=p_primer->next;
        }
	free(list);
}

//get the file size
int file_size2(char* filename)
{
        struct stat statbuf;
        stat(filename,&statbuf);
        int size=statbuf.st_size;
        return size;
}

////function read primer informatin and align information 
struct Primer *read_par(char *path,int common_flag,int special_flag)
{
        char *in;
        int pos,len,gi,position,plus,minus,size,i,flag;
        struct Primer *new_primer,*p_primer,*head;
        struct Node *new_node,*p_node;
        FILE *fp;

///read the  primer file
        if(access(path,0)==-1)
        {
                printf("Error! Don't have the %s file!\n",path);
                exit(1);
        }
        fp=fopen(path,"r");
        if(fp==NULL)
        {
                printf("Error: can't open the %s file!\n",path);
                exit(1);
        }
        
        size=sizeof(struct Primer);
        i=0;
        while(fscanf(fp,"pos:%d\tlength:%d\t+:%d\t-:%d\n",&pos,&len,&plus,&minus)!=EOF)
        {
                new_primer=(struct Primer *)malloc(size);
                new_primer->pos=pos;
                new_primer->len=len;
                new_primer->total=1;
		new_primer->total_common=0;
		new_primer->total_special=0;
                new_primer->plus=plus;
                new_primer->minus=minus;
                new_primer->next=NULL;
                new_primer->common=NULL;
                new_primer->special=NULL;

                if(i==0)
                {
                        head=new_primer;
                        p_primer=new_primer;
                        i++;
                }
                else
                {
                        p_primer->next=new_primer;
                        p_primer=new_primer;
                }
        }
        fclose(fp);
        if(i==0)
        {
                printf("Sorry! Don't have any candidate single primers in %s!\n",path);
                exit(1);
        }

//parameter of common
        if(common_flag==1)
        {
                i=strlen(path);
                in=(char *)malloc(i+20);
                memset(in,'\0',i+20);
                strcpy(in,path);
                strcat(in,"-common.txt"); //suffix of parameter
                if(access(in,0)==-1)
                {
                        printf("Error! Don't have the %s file!\n",in);
                        exit(1);
                }

                fp=fopen(in,"r");
                if(fp==NULL)
                {
                        printf("Error: can't open the %s file!\n",in);
                        exit(1);
                }

                p_primer=head;
                size=sizeof(struct Node);
                while(fscanf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",&pos,&len,&gi,&position,&plus,&minus)!=EOF)
                {
                        new_node=(struct Node *)malloc(size);
                        new_node->pos=position;
                        new_node->gi=gi;
                        new_node->plus=plus;
                        new_node->minus=minus;

        //find the primer
                        flag=0;
                        while((p_primer->pos!=pos||p_primer->len!=len)&&flag<2)
                        {
                                if((p_primer->next==NULL)||(p_primer->pos>pos))
                                {
                                        flag++;
                                        p_primer=head;
                                }
                                else
                                {
                                        p_primer=p_primer->next;
                                }
                        }
                        if(flag==2)
                        {
                                printf("Warning: the single primer(pos is %d, length is %d) is not in %s!\n",pos,len,path);
                                free(new_node);
                                continue;
                        } 
                        p_node=p_primer->common;
                        p_primer->common=new_node;
			p_primer->total_common++;
                        new_node->next=p_node;
                }
                fclose(fp);
                free(in);
        }
//paramter for special
        if(special_flag==1)
        {
                i=strlen(path);
                in=(char *)malloc(i+20);
                memset(in,'\0',i+20);
                strcpy(in,path);
                strcat(in,"-special.txt"); //suffix of parameter
                if(access(in,0)==-1)
                {
                        printf("Error! Don't have the %s file!\n",in);
                        exit(1);
                }

                fp=fopen(in,"r");
                if(fp==NULL)
                {
                        printf("Error: can't open the %s file!\n",in);
                        exit(1);
                }
        
                p_primer=head;
                while(fscanf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",&pos,&len,&gi,&position,&plus,&minus)!=EOF)
                {
                        new_node=(struct Node *)malloc(size);
                        new_node->pos=position;
                        new_node->gi=gi;
                        new_node->plus=plus;
                        new_node->minus=minus;
        
                        //find the primer
                        flag=0;
                        while((p_primer->pos!=pos||p_primer->len!=len)&&flag<2)
                        {
                                if((p_primer->next==NULL)||(p_primer->pos>pos))
                                {
                                        flag++;
                                        p_primer=head;
                                }
                                else
                                        p_primer=p_primer->next;
                        }
                        if(flag==2)
                        {
                                printf("Warning: the single primer(pos is %d, length is %d) is not in %s!\n",pos,len,path);
                                free(new_node);
                                continue;
                        }
                        p_node=p_primer->special;
                        p_primer->special=new_node;
			p_primer->total_special++;
                        new_node->next=p_node;
                }
                fclose(fp);
                free(in);
        }
        return head;
}

//function: check how many GIs this primer can be used for
__device__ int check_common(int *d_L,int *d_S,int *d_cL,int *d_cS,int *d_scL,int *d_scS,int *d_ecL,int *d_ecS,int turn[],int common,int *d_apply)
{
        int pos[6],i,dis;

        for(i=0;i<common;i++)
        {
                d_apply[common*turn[0]+i]=0;
        }
//plus
        for(pos[0]=d_scS[turn[0]];pos[0]<d_ecS[turn[0]];pos[0]++)
        {
                if(d_cS[4*pos[0]+2]!=1)
                        continue;
		i=d_cS[4*pos[0]];
                if(d_apply[common*turn[0]+i]==1)
                        continue;
                for(pos[1]=d_scS[turn[1]];pos[1]<d_ecS[turn[1]];pos[1]++)
                {
                        if(d_cS[4*pos[1]]!=i)
                                continue;
                        if(d_cS[4*pos[1]+2]!=1)
                                continue;
                        for(pos[2]=d_scL[turn[2]];pos[2]<d_ecL[turn[2]];pos[2]++)
                        {
                                if(d_cL[4*pos[2]]!=i)
                                        continue;
                                if(d_cL[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=d_scL[turn[3]];pos[3]<d_ecL[turn[3]];pos[3]++)
                                {
                                        if(d_cL[4*pos[3]]!=i)
                                                continue;
                                        if(d_cL[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=d_scS[turn[4]];pos[4]<d_ecS[turn[4]];pos[4]++)
                                        {
                                                if(d_cS[4*pos[4]]!=i)
                                                        continue;
                                                if(d_cS[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=d_scS[turn[5]];pos[5]<d_ecS[turn[5]];pos[5]++)
                                                {
                                                        if(d_cS[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_cS[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_cS[4*pos[1]+1]-(d_cS[4*pos[0]+1]+d_S[4*turn[0]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=d_cL[4*pos[2]+1]-d_cS[4*pos[1]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_cL[4*pos[3]+1]-(d_cL[4*pos[2]+1]+d_L[4*turn[2]+1]);
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=(d_cS[4*pos[4]+1]+d_S[4*turn[4]+1]-1)-(d_cL[4*pos[3]+1]+d_L[4*turn[3]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_cS[4*pos[4]+1]+d_S[4*turn[4]+1]-1-d_cS[4*pos[1]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_cS[4*pos[5]+1]-(d_cS[4*pos[4]+1]+d_S[4*turn[4]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                        d_apply[common*turn[0]+i]=1;
                                                }
                                        }
                                }
                        }
                }
        }
//minus
        for(pos[0]=d_scS[turn[0]];pos[0]<d_ecS[turn[0]];pos[0]++)
        {
                if(d_cS[4*pos[0]+3]!=1)
                        continue;
		i=d_cS[4*pos[0]];
                if(d_apply[turn[0]*common+i]==1)
                        continue;  //this GI can common

                for(pos[1]=d_scS[turn[1]];pos[1]<d_ecS[turn[1]];pos[1]++)
                {
                        if(d_cS[4*pos[1]]!=i)
                                continue;
                        if(d_cS[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=d_scL[turn[2]];pos[2]<d_ecL[turn[2]];pos[2]++)
                        {
                                if(d_cL[4*pos[2]]!=i)
                                        continue;
                                if(d_cL[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=d_scL[turn[3]];pos[3]<d_ecL[turn[3]];pos[3]++)
                                {
                                        if(d_cL[4*pos[3]]!=i)
                                                continue;
                                        if(d_cL[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=d_scS[turn[4]];pos[4]<d_ecS[turn[4]];pos[4]++)
                                        {
                                                if(d_cS[4*pos[4]]!=i)
                                                        continue;
                                                if(d_cS[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=d_scS[turn[5]];pos[5]<d_ecS[turn[5]];pos[5]++)
                                                {
                                                        if(d_cS[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_cS[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_cS[4*pos[0]+1]-(d_cS[4*pos[1]+1]+d_S[4*turn[1]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=(d_cS[4*pos[1]+1]+d_S[4*turn[1]+1]-1)-(d_cL[4*pos[2]+1]+d_L[4*turn[2]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_cL[4*pos[2]+1]-(d_cL[4*pos[3]+1]+d_L[4*turn[3]+1]);
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=d_cL[4*pos[3]+1]-d_cS[4*pos[4]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_cS[4*pos[1]+1]+d_S[4*turn[1]+1]-1-d_cS[4*pos[4]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_cS[4*pos[4]+1]-(d_cS[4*pos[5]+1]+d_S[4*turn[5]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                        d_apply[common*turn[0]+i]=1;
                                                }
                                        }
                                }
                        }
                }
        }
        dis=0;
        for(i=0;i<common;i++)
        {
                dis=dis+d_apply[common*turn[0]+i];
        }
        return dis;
}

//check this LAMP primers are uniq or not
//return=0: stop and return=1: go on
__device__ int check_uniq(int *d_L,int *d_S,int *d_sL,int *d_sS,int *d_ssL,int *d_ssS,int *d_esL,int *d_esS,int turn[])
{
        int pos[6],gi;

//plus
        for(pos[0]=d_ssS[turn[0]];pos[0]<d_esS[turn[0]];pos[0]++)
        {
                if(d_sS[4*pos[0]+2]!=1)
                        continue;
		gi=d_sS[4*pos[0]];
                for(pos[1]=d_ssS[turn[1]];pos[1]<d_esS[turn[1]];pos[1]++)
                {
			if(d_sS[4*pos[1]]!=gi)
                                continue;
                        if(d_sS[4*pos[1]+2]!=1)
				continue;
                        for(pos[2]=d_ssL[turn[2]];pos[2]<d_esL[turn[2]];pos[2]++) //F1c
                        {
                                if(d_sL[4*pos[2]]!=gi)
                                        continue;
                                if(d_sL[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=d_ssL[turn[3]];pos[3]<d_esL[turn[3]];pos[3]++) //B1c
                                {
                                        if(d_sL[pos[3]*4]!=gi)
                                                continue;
                                        if(d_sL[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=d_ssS[turn[4]];pos[4]<d_esS[turn[4]];pos[4]++) //B2
                                        {
                                                if(d_sS[4*pos[4]]!=gi)
                                                        continue;
                                                if(d_sS[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=d_ssS[turn[5]];pos[5]<d_esS[turn[5]];pos[5]++)
                                                {
                                                        if(d_sS[4*pos[5]]!=gi)
                                                                continue;
                                                        if(d_sS[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(d_sS[4*pos[1]+1]<d_sS[4*pos[0]+1])
                                                                continue;
                                                //F2-F1c
                                                        if(d_sL[4*pos[2]+1]<d_sS[4*pos[1]+1]+d_S[4*turn[1]+1])
                                                                continue;
                                                //F1c-B1c
                                                        if(d_sL[4*pos[3]+1]<d_sL[4*pos[2]+1]+d_L[4*turn[2]+1])
                                                                continue;
                                                //B1c-B2
                                                        if(d_sS[4*pos[4]+1]<d_sL[4*pos[3]+1]+d_L[4*turn[3]+1])
                                                                continue;
                                                //B2-B3
                                                        if(d_sS[4*pos[5]+1]<d_sS[4*pos[4]+1])
                                                                continue;
                                                //whole
                                                        if(d_sS[4*pos[5]+1]-d_sS[4*pos[0]+1]>1000)
                                                                continue;
                                                        return 0;
                                                }//B3
                                        }
                                }//B1c
                        }
                }//F2
        }

//minus
        for(pos[0]=d_ssS[turn[0]];pos[0]<d_esS[turn[0]];pos[0]++)
        {
                if(d_sS[4*pos[0]+3]!=1)
                        continue;
		gi=d_sS[4*pos[0]];
                for(pos[1]=d_ssS[turn[1]];pos[1]<d_esS[turn[1]];pos[1]++)
                {
                        if(d_sS[4*pos[1]]!=gi)
                                continue;
                        if(d_sS[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=d_ssL[turn[2]];pos[2]<d_esL[turn[2]];pos[2]++)
                        {
                                if(d_sL[4*pos[2]]!=gi)
                                        continue;
                                if(d_sL[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=d_ssL[turn[3]];pos[3]<d_esL[turn[3]];pos[3]++)
                                {
                                        if(d_sL[4*pos[3]]!=gi)
                                                continue;
                                        if(d_sL[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=d_ssS[turn[4]];pos[4]<d_esS[turn[4]];pos[4]++)
                                        {
                                                if(d_sS[4*pos[4]]!=gi)
                                                        continue;
                                                if(d_sS[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=d_ssS[turn[5]];pos[5]<d_esS[turn[5]];pos[5]++)
                                                {
                                                        if(d_sS[4*pos[5]]!=gi)
                                                                continue;
                                                        if(d_sS[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(d_sS[4*pos[0]+1]<d_sS[4*pos[1]+1])
                                                                continue;
                                                //F2-F1c
                                                        if(d_sS[4*pos[1]+1]<d_sL[4*pos[2]+1]+d_L[4*turn[2]+1])
                                                                continue;
                                                //F1c-B1c
                                                        if(d_sL[4*pos[2]+1]<d_sL[4*pos[3]+1]+d_L[4*turn[3]+1])
                                                                continue;
                                                //B1c-B2
                                                        if(d_sL[4*pos[3]+1]<d_sS[4*pos[4]+1]+d_S[4*turn[4]+1])
                                                                continue;
                                                //B2-B3
                                                        if(d_sS[4*pos[4]+1]<d_sS[4*pos[5]+1])
                                                                continue;
                                                //whole
                                                        if(d_sS[4*pos[0]+1]-d_sS[4*pos[5]+1]>1000)
                                                                continue;
                                                        return 0;
                                                }
                                        }
                                }
                        }
                }
        }
        return 1;
}

//from first to second
__global__ void next_one(int *first,int *second,int *next,int num_first,int num_second)
{
        int id=blockDim.x*blockIdx.x+threadIdx.x;
	int i;

	while(id<num_first)
	{
		i=id;
		if(i>=num_second)
		{
			i=num_second-1;
		}
		if(second[4*i]>=first[id*4]+first[id*4+1])
		{
			while((i>=0)&&(second[4*i]>=first[id*4]+first[id*4+1]))
			{
				next[id]=i;
				i--;
			}
		}
		else
		{
			while((i<num_second)&&(second[i*4]<first[id*4]+first[id*4+1]))
				i++;
			next[id]=i;
			if(i==num_second)
				next[id]=-1;
		}
		id=id+blockDim.x*gridDim.x;
	}
	__syncthreads();
}

__device__ int check_gc(char *d_seq,int start,int end,int flag)
{
        int i,total=0;
        float gc;

        for(i=start;i<end;i++)
        {
                if(d_seq[i]=='C'||d_seq[i]=='G')
                        total++;
        }
        gc=total*100.0/(end-start);
        if(flag==1&&gc>=45)
                return 1;
        if(flag==0&&gc<=45)
                return 1;
        return 0;
}

__device__ int check_common_loop(int *d_S,int *d_L,int *d_Lp,int *d_cS,int *d_cL,int *d_cLp,int *d_scS,int *d_ecS,int *d_scL,int *d_ecL,int *d_scLp,int *d_ecLp,int turn[],int LF,int LB,int *d_int,int *d_apply)
{
        int dis,i,pos[7];
//plus
        for(pos[0]=d_scS[turn[0]];pos[0]<d_ecS[turn[0]];pos[0]++)
        {
                if(d_cS[4*pos[0]+2]!=1)
                        continue;
		i=d_cS[4*pos[0]];
                if(d_apply[turn[0]*d_int[7]+i]==0||d_apply[turn[0]*d_int[7]+i]==2)
                        continue;
                for(pos[1]=d_scS[turn[1]];pos[1]<d_ecS[turn[1]];pos[1]++)
                {
                        if(d_cS[4*pos[1]]!=i)
                                continue;
                        if(d_cS[4*pos[1]+2]!=1)
                                continue;
                        for(pos[2]=d_scL[turn[2]];pos[2]<d_ecL[turn[2]];pos[2]++)
                        {
                                if(d_cL[4*pos[2]]!=i)
                                        continue;
                                if(d_cL[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=d_scL[turn[3]];pos[3]<d_ecL[turn[3]];pos[3]++)
                                {
                                        if(d_cL[4*pos[3]]!=i)
                                                continue;
                                        if(d_cL[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=d_scS[turn[4]];pos[4]<d_ecS[turn[4]];pos[4]++)
                                        {
                                                if(d_cS[4*pos[4]]!=i)
                                                        continue;
                                                if(d_cS[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=d_scS[turn[5]];pos[5]<d_ecS[turn[5]];pos[5]++)
                                                {
                                                        if(d_cS[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_cS[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_cS[4*pos[1]+1]-(d_cS[4*pos[0]+1]+d_S[4*turn[0]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=d_cL[4*pos[2]+1]-d_cS[4*pos[1]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_cL[4*pos[3]+1]-(d_cL[4*pos[2]+1]+d_L[4*turn[2]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=(d_cS[4*pos[4]+1]+d_S[4*turn[4]+1]-1)-(d_cL[4*pos[3]+1]+d_L[4*turn[3]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_cS[4*pos[4]+1]+d_S[4*turn[4]+1]-1-d_cS[4*pos[1]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_cS[4*pos[5]+1]-(d_cS[4*pos[4]+1]+d_S[4*turn[4]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //LF
                                                        if(LF!=-1)
                                                        {
                                                                dis=0;
                                                                for(pos[6]=d_scLp[LF];pos[6]<d_ecLp[LF];pos[6]++)
                                                                {
                                                                        if(d_cLp[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+3]!=1)
                                                                                continue;
                                                                        if(d_cS[4*pos[1]+1]+d_S[4*turn[1]+1]>d_cLp[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+1]+d_Lp[4*LF+1]>d_cL[4*pos[2]+1])
                                                                                continue;
                                                                        dis=1;
                                                                        break;
                                                                }
                                                                if(dis==0)
                                                                        continue;
                                                        }
                                                //LB
                                                        if(LB!=-1)
                                                        {
                                                                dis=0;
                                                                for(pos[6]=d_scLp[LB];pos[6]=d_ecLp[LB];pos[6]++)
                                                                {
                                                                        if(d_cLp[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+2]!=1)
                                                                                continue;
                                                                        if(d_cL[4*pos[3]+1]+d_L[4*turn[3]+1]>d_cLp[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+1]+d_Lp[4*LB+1]>d_cS[4*pos[4]+1])
                                                                                continue;
                                                                        dis=1;
                                                                        break;
                                                                }
                                                                if(dis==0)
                                                                        continue;
                                                        }
                                                        d_apply[turn[0]*d_int[7]+i]=2;
                                                }
                                        }
                                }
                        }
                }
        }
//minus
	for(pos[0]=d_scS[turn[0]];pos[0]<d_ecS[turn[0]];pos[0]++)
        {
                if(d_cS[4*pos[0]+3]!=1)
                        continue;
                i=d_cS[4*pos[0]];
                if(d_apply[turn[0]*d_int[7]+i]==0||d_apply[turn[0]*d_int[7]+i]==2)
                        continue;
                for(pos[1]=d_scS[turn[1]];pos[1]<d_ecS[turn[1]];pos[1]++)
                {
                        if(d_cS[4*pos[1]]!=i)
                                continue;
                        if(d_cS[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=d_scL[turn[2]];pos[2]<d_ecL[turn[2]];pos[2]++)
                        {
                                if(d_cL[4*pos[2]]!=i)
                                        continue;
                                if(d_cL[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=d_scL[turn[3]];pos[3]<d_ecL[turn[3]];pos[3]++)
                                {
                                        if(d_cL[4*pos[3]]!=i)
                                                continue;
                                        if(d_cL[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=d_scS[turn[4]];pos[4]<d_ecS[turn[4]];pos[4]++)
                                        {
                                                if(d_cS[4*pos[4]]!=i)
                                                        continue;
                                                if(d_cS[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=d_scS[turn[5]];pos[5]<d_ecS[turn[5]];pos[5]++)
                                                {
                                                        if(d_cS[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_cS[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_cS[4*pos[0]+1]-(d_cS[4*pos[1]+1]+d_S[4*turn[1]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=(d_cS[4*pos[1]+1]+d_S[4*turn[1]+1]-1)-(d_cL[4*pos[2]+1]+d_L[4*turn[2]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_cL[4*pos[2]+1]-(d_cL[4*pos[3]+1]+d_L[4*turn[3]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=d_cL[4*pos[3]+1]-d_cS[4*pos[4]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_cS[4*pos[1]+1]+d_S[4*turn[1]+1]-1-d_cS[4*pos[4]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_cS[4*pos[4]+1]-(d_cS[4*pos[5]+1]+d_S[4*turn[5]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //LF
                                                        if(LF!=-1)
                                                        {
                                                                dis=0;
                                                                for(pos[6]=d_scLp[LF];pos[6]<d_ecLp[LF];pos[6]++)
                                                                {
                                                                        if(d_cLp[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+2]!=1)
                                                                                continue;
                                                                        if(d_cL[4*pos[2]+1]+d_L[4*turn[2]+1]>d_cLp[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+1]+d_Lp[4*LF+1]>d_cS[4*pos[1]+1])
                                                                                continue;
                                                                        dis=1;
                                                                        break;
                                                                }
                                                                if(dis==0)
                                                                        continue;
                                                        }
                                                //LB
                                                        if(LB!=-1)
                                                        {
                                                                dis=0;
                                                                for(pos[6]=d_scLp[LB];pos[6]=d_ecLp[LB];pos[6]++)
                                                                {
                                                                        if(d_cLp[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+3]!=1)
                                                                                continue;
                                                                        if(d_cS[4*pos[4]+1]+d_S[4*turn[4]+1]>d_cLp[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_cLp[4*pos[6]+1]+d_Lp[4*LB+1]>d_cL[4*pos[3]+1])
                                                                                continue;
                                                                        dis=1;
                                                                        break;
                                                                }
                                                                if(dis==0)
                                                                        continue;
                                                        }
							d_apply[turn[0]*d_int[7]+i]=2;
                                                }
                                        }
                                }
                        }
                }
        }
        for(i=0;i<d_int[7];i++)
        {
                if(d_apply[turn[0]*d_int[7]+i]==1)
                        return 0;
        }
        return 1;
}

__device__ int check_structure_loop(char *list[],char *rev[],int flag,double *parameter,char *d_Pchar,int *d_NumL)
{
        int i;
        double TH;

        if(list[6]!=NULL)
        {
                for(i=0;i<=1;i++)
                {
                        TH=thal(list[i],list[6],1,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[i],list[6],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[i],list[6],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;

                        TH=thal(rev[6],rev[i],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(rev[6],rev[i],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                }

                for(i=2;i<6;i++)
                {
                        TH=thal(list[6],list[i],1,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[6],list[i],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[6],list[i],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;

                        TH=thal(rev[i],rev[6],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag) 
                                return 0;
                        TH=thal(rev[i],rev[6],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                }
        }
        if(list[7]!=NULL)
        {
                for(i=0;i<4;i++)
                {
                        TH=thal(list[i],list[7],1,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[i],list[7],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[i],list[7],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;

                        TH=thal(rev[7],rev[i],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(rev[7],rev[i],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                }
                for(i=4;i<6;i++)
                {
                        TH=thal(list[7],list[i],1,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[7],list[i],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(list[7],list[i],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;

                        TH=thal(rev[i],rev[7],2,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                        TH=thal(rev[i],rev[7],3,parameter,d_Pchar,d_NumL);
                        if(TH>44+5*flag)
                                return 0;
                }
        }
        if(list[6]!=NULL&&list[7]!=NULL)
        {
                TH=thal(list[6],list[7],1,parameter,d_Pchar,d_NumL);
                if(TH>44+5*flag)
                        return 0;
                TH=thal(list[6],list[7],2,parameter,d_Pchar,d_NumL);
                if(TH>44+5*flag)
                        return 0;
                TH=thal(list[6],list[7],3,parameter,d_Pchar,d_NumL);
                if(TH>44+5*flag)
                        return 0;

                TH=thal(rev[7],rev[6],2,parameter,d_Pchar,d_NumL);
                if(TH>44+5*flag)
                        return 0;
                TH=thal(rev[7],rev[6],3,parameter,d_Pchar,d_NumL);
                if(TH>44+5*flag)
                        return 0;
        }
        return 1;
}

__device__ int design_loop(int *d_S,int *d_L,int *d_Lp,int *d_SLp,int *d_LLp,int *d_LpLp,char *d_seq,int *d_cS,int *d_cL,int *d_cLp,int *d_scS,int *d_ecS,int *d_scL,int *d_ecL,int *d_scLp,int *d_ecLp,int turn[],int *d_int,int *d_apply,int *d_par,double *parameter,char *d_Pchar,int *d_NumL)
{
        int success,LF,LB;
        char P_F3[26],P_F2[26],P_F1c[26],P_B1c[26],P_B2[26],P_B3[26],P_LF[26],P_LB[26],*list[8];
	char r_F3[26],r_F2[26],r_F1c[26],r_B1c[26],r_B2[26],r_B3[26],r_LF[26],r_LB[26],*rev[8];

        if(d_int[6])
        {
		get_primer(d_seq,P_F3,d_S[4*turn[0]],d_S[4*turn[0]+1],0);
	        get_primer(d_seq,P_F2,d_S[4*turn[1]],d_S[4*turn[1]+1],0);
	        get_primer(d_seq,P_F1c,d_L[4*turn[2]],d_L[4*turn[2]+1],1);
	        get_primer(d_seq,P_B1c,d_L[4*turn[3]],d_L[4*turn[3]+1],0);
	        get_primer(d_seq,P_B2,d_S[4*turn[4]],d_S[4*turn[4]+1],1);
	        get_primer(d_seq,P_B3,d_S[4*turn[5]],d_S[4*turn[5]+1],1);

		get_primer(d_seq,r_F3,d_S[4*turn[0]],d_S[4*turn[0]+1],1);
                get_primer(d_seq,r_F2,d_S[4*turn[1]],d_S[4*turn[1]+1],1);
                get_primer(d_seq,r_F1c,d_L[4*turn[2]],d_L[4*turn[2]+1],0);
                get_primer(d_seq,r_B1c,d_L[4*turn[3]],d_L[4*turn[3]+1],1);
                get_primer(d_seq,r_B2,d_S[4*turn[4]],d_S[4*turn[4]+1],0);
                get_primer(d_seq,r_B3,d_S[4*turn[5]],d_S[4*turn[5]+1],0);
	
		list[0]=P_F3;
		list[1]=P_F2;
		list[2]=P_F1c;
		list[3]=P_B1c;
		list[4]=P_B2;
		list[5]=P_B3;
		rev[0]=r_F3;
		rev[1]=r_F2;
		rev[2]=r_F1c;
		rev[3]=r_B1c;
		rev[4]=r_B2;
		rev[5]=r_B3;
        }
//LF and LB 
        success=0;
	LF=d_SLp[turn[1]];
        while(LF<d_int[2])
        {
		if(LF==-1)
			break;
		if(d_Lp[4*LF+3]!=1)
		{
			LF++;
			continue;
		}
                if(d_Lp[4*LF]+18>d_L[4*turn[2]])
                        break;
                LB=d_LLp[turn[3]];
		if(LB==-1||d_Lp[4*LB]+18>d_S[4*turn[4]])
			break;
                if(d_int[6])
		{
			get_primer(d_seq,P_LF,d_Lp[4*LF],d_Lp[4*LF+1],1);
			get_primer(d_seq,r_LF,d_Lp[4*LF],d_Lp[4*LF+1],0);
			list[6]=P_LF;
			rev[6]=r_LF;
		}
                while(LB<d_int[2])
                {
			if(d_Lp[4*LB+2]!=1)
			{
				LB++;
				continue;
			}
                        if(d_Lp[4*LB]+18>d_S[4*turn[4]])
                                break;
                //check_common
                        if(d_int[3])
                        {
                                success=check_common_loop(d_S,d_L,d_Lp,d_cS,d_cL,d_cLp,d_scS,d_ecS,d_scL,d_ecL,d_scLp,d_ecLp,turn,LF,LB,d_int,d_apply);
                                if(success==0)
                                {
                                        LB++;
                                        continue;
                                }
                        }
                //check_structure
                        if(d_int[6])
                        {
				get_primer(d_seq,P_LB,d_Lp[4*LB],d_Lp[4*LB+1],0);
				get_primer(d_seq,r_LB,d_Lp[4*LB],d_Lp[4*LB+1],1);
				list[7]=P_LB;
				rev[7]=r_LB;
                                success=check_structure_loop(list,rev,d_int[9],parameter,d_Pchar,d_NumL);
                                if(success==0)
                                {
                                        LB++;
                                        continue;
                                }
                        }
                        d_par[16*turn[0]+12]=d_Lp[4*LF];
			d_par[16*turn[0]+13]=d_Lp[4*LF+1];
			d_par[16*turn[0]+14]=d_Lp[4*LB];
			d_par[16*turn[0]+15]=d_Lp[4*LB+1];
                        success=1;
                        break;
                }
                if(success==1)
                        break;
                else
                        LF++;
        }
        if(success==1)
                return success;
//only LF
        LF=d_SLp[turn[1]];
        while(LF<d_int[2])
        {
		if(LF==-1)
			break;
                if(d_Lp[4*LF]+18>d_L[4*turn[2]])
                        break;
		if(d_Lp[4*LF+3]!=1)
		{
			LF++;
			continue;
		}
        //check_common
                if(d_int[3])
                {
                        success=check_common_loop(d_S,d_L,d_Lp,d_cS,d_cL,d_cLp,d_scS,d_ecS,d_scL,d_ecL,d_scLp,d_ecLp,turn,LF,-1,d_int,d_apply);
                        if(success==0)
                        {
                                LF++;
                                continue;
                        }
                }
        //check_structure
                if(d_int[6])
                {
                        get_primer(d_seq,P_LF,d_Lp[4*LF],d_Lp[4*LF+1],1);
			get_primer(d_seq,r_LF,d_Lp[4*LF],d_Lp[4*LF+1],0);
			list[6]=P_LF;
			rev[6]=r_LF;
			list[7]=NULL;
			rev[7]=NULL;
                        success=check_structure_loop(list,rev,d_int[9],parameter,d_Pchar,d_NumL);
                        if(success==0)
                        {
                                LF++;
                                continue;
                        }
                }
		d_par[16*turn[0]+12]=d_Lp[4*LF];
		d_par[16*turn[0]+13]=d_Lp[4*LF+1];
		d_par[16*turn[0]+14]=0;
		d_par[16*turn[0]+15]=0;
                success=1;
                break;
        }
        if(success==1)
                return success;
//only LB
        LB=d_LLp[turn[3]];
        while(LB<d_int[2])
        {
		if(LB==-1)
			break;
                if(d_Lp[LB*4]+18>d_S[4*turn[4]])
                        break;
		if(d_Lp[LB*4+2]!=1)
		{
			LB++;
			continue;
		}
        //check_common
                if(d_int[3])
                {
                        success=check_common_loop(d_S,d_L,d_Lp,d_cS,d_cL,d_cLp,d_scS,d_ecS,d_scL,d_ecL,d_scLp,d_ecLp,turn,-1,LB,d_int,d_apply);
                        if(success==0)
                        {
                                LB++;
                                continue;
                        }
                }
        //check_structure
                if(d_int[6])
                {
                        get_primer(d_seq,P_LB,d_Lp[4*LB],d_Lp[4*LB+1],0);
			get_primer(d_seq,r_LB,d_Lp[4*LB],d_Lp[4*LB+1],1);
			list[6]=NULL;
			rev[6]=NULL;
			list[7]=P_LB;
			rev[7]=r_LB;
                        success=check_structure_loop(list,rev,d_int[9],parameter,d_Pchar,d_NumL);
                        if(success==0)
                        {
                                LB++;
                                continue;
                        }
                }
		d_par[16*turn[0]+12]=0;
		d_par[16*turn[0]+13]=0;
		d_par[16*turn[0]+14]=d_Lp[4*LB];
		d_par[16*turn[0]+15]=d_Lp[4*LB+1];
                success=1;
                break;
        }
        return success;
}

//caculate
__global__ void LAMP(char *d_seq,int *d_S,int *d_L,int *d_Lp,int *d_cS,int *d_cL,int *d_cLp,int *d_sS,int *d_sL,int *d_scS,int *d_scL,int *d_scLp,int *d_ecS,int *d_ecL,int *d_ecLp,int *d_ssS,int *d_ssL,int *d_esS,int *d_esL,int *d_SS,int *d_SL,int *d_SLp,int *d_LL,int *d_LS,int *d_LLp,int *d_LpLp,int *d_int,int *d_par,int *d_apply,double *parameter,char *d_Pchar,int *d_NumL,int *d_pos)
//d_int: 0:numS,1:numL,2:numLp,3:common_flag,4:special_flag,5:loop_flag,6:secondary_flag,7:common_num,8:this turn common_num,9:high_GC_flag; 10:expect
{
	int id=blockDim.x*blockIdx.x+threadIdx.x;
	int turn[6],flag;

	while(id<d_int[0])
	{
		d_par[16*id+1]=0;//not LAMP, as a flag
//check add by F3'pos
		turn[0]=0;
		for(turn[1]=0;turn[1]<d_int[10];turn[1]++)
		{
			if(d_pos[turn[1]]==-1)
				break;
			if(d_S[4*id]-d_pos[turn[1]]<300&&d_S[4*id]-d_pos[turn[1]]>-300)
			{
				turn[0]++;
				break;
			}
		}
		if(turn[0])
		{
			id=id+blockDim.x*gridDim.x;
			continue;
		}
		if(d_S[4*id+2]!=1)
		{
			id=id+blockDim.x*gridDim.x;	
			continue;
		}
	//combine
		turn[0]=id; //one thread, one F3
		flag=0;
		for(turn[1]=d_SS[turn[0]];turn[1]<d_int[0];turn[1]++) //F2
		{
			if(turn[1]==-1)
				break;
			if(flag!=0)
				break; //have find one LAMP primer
			if(d_S[4*turn[1]+2]!=1)
				continue;
			if(d_S[4*turn[1]]-(d_S[4*turn[0]]+d_S[4*turn[0]+1])>20)
				break;
			for(turn[2]=d_SL[turn[1]];turn[2]<d_int[1];turn[2]++) //F1c
			{
				if(turn[2]==-1)
					break;
				if(flag!=0)
					break;
				if(d_L[4*turn[2]+3]!=1)
					continue;
				if(d_L[4*turn[2]]-d_S[4*turn[1]]-1<40)
					continue;
                                if(d_L[4*turn[2]]-d_S[4*turn[1]]-1>60)
                                	break;
                                for(turn[3]=d_LL[turn[2]];turn[3]<d_int[1];turn[3]++)   //B1c
                                {
                                        if(turn[3]==-1)
                                        	break;
					if(flag!=0)
						break;
					if(d_L[4*turn[3]+2]!=1)
						continue;
                                        if(d_L[4*turn[3]]-d_L[4*turn[2]]>85)
                                        	break;
                                        for(turn[4]=d_LS[turn[3]];turn[4]<d_int[0];turn[4]++)   //B2
                                        {
                                                if(turn[4]==-1)
                                                	break;
						if(flag!=0)
							break;
						if(d_S[4*turn[4]+3]!=1)
							continue;
                                                if((d_S[4*turn[4]]+d_S[4*turn[4]+1]-1)-(d_L[4*turn[3]]+d_L[4*turn[3]+1])<40)
                                                	continue;
                                                if((d_S[4*turn[4]]+d_S[4*turn[4]+1]-1)-(d_L[4*turn[3]]+d_L[4*turn[3]+1])>60)
                                                	break;
                                                if(d_S[4*turn[4]]+d_S[4*turn[4]+1]-1-d_S[turn[1]*4]-1<120)
                                                	continue;
                                                if(d_S[4*turn[4]]+d_S[4*turn[4]+1]-1-d_S[turn[1]*4]-1>180)
                                                	break;
						if(d_int[5]&&(d_SLp[turn[1]]==-1||(d_Lp[4*d_SLp[turn[1]]]+18>d_L[4*turn[2]]))&&(d_LLp[turn[3]]==-1||(d_Lp[4*d_LLp[turn[3]]]+18>d_S[4*turn[4]])))
							continue;
                                                for(turn[5]=d_SS[turn[4]];turn[5]<d_int[0];turn[5]++)  //B3
                                                {
                                                        if(turn[5]==-1)
                                                        	break;
							if(d_S[4*turn[5]+3]!=1)
								continue;
                                                        if(d_S[turn[5]*4]-(d_S[4*turn[4]]+d_S[4*turn[4]+1])>20)
                                                        	break;
							flag=check_gc(d_seq,d_S[4*turn[0]],(d_S[4*turn[5]]+d_S[4*turn[5]+1]),d_int[9]);
							if(flag==0)
								continue;
							if(d_int[4])
							{
								flag=check_uniq(d_L,d_S,d_sL,d_sS,d_ssL,d_ssS,d_esL,d_esS,turn);
								if(flag==0)
									continue;
							}
							if(d_int[3])
							{
								flag=check_common(d_L,d_S,d_cL,d_cS,d_scL,d_scS,d_ecL,d_ecS,turn,d_int[7],d_apply);
								if(flag<d_int[8])
									continue;
							}

							if(d_int[6])
							{
								flag=check_structure(d_seq,d_S,d_L,turn,d_int[9],parameter,d_Pchar,d_NumL);
								if(flag==0)
									continue;
							}
	
							if(d_int[5])
							{
								flag=design_loop(d_S,d_L,d_Lp,d_SLp,d_LLp,d_LpLp,d_seq,d_cS,d_cL,d_cLp,d_scS,d_ecS,d_scL,d_ecL,d_scLp,d_ecLp,turn,d_int,d_apply,d_par,parameter,d_Pchar,d_NumL);
								if(flag==0)
									continue;
							}
							d_par[id*16]=d_S[4*turn[0]];
							d_par[id*16+1]=d_S[4*turn[0]+1];
							d_par[id*16+2]=d_S[4*turn[1]];
							d_par[id*16+3]=d_S[4*turn[1]+1];
							d_par[id*16+4]=d_L[4*turn[2]];
							d_par[id*16+5]=d_L[4*turn[2]+1];
							d_par[id*16+6]=d_L[4*turn[3]];
							d_par[id*16+7]=d_L[4*turn[3]+1];
							d_par[id*16+8]=d_S[4*turn[4]];
							d_par[id*16+9]=d_S[4*turn[4]+1];
							d_par[id*16+10]=d_S[4*turn[5]];
							d_par[id*16+11]=d_S[4*turn[5]+1];
							break;
						}
					}
				}
			}
		}
		id=id+blockDim.x*gridDim.x;
	}
	__syncthreads();
}

void usage()
{
        printf("Usage:\n");
        printf("    LAMP_GPU  -in <name>  -out <result>  -high[-low] [options]*\n\n");
        printf("    -in   <string>:  the name of candidate single primers file\n");
        printf("    -out  <string>:  the result file of LAMP primers\n");
        printf("    -dir  <string>:  the directory to store candidate single primers, default is current directory\n");
        printf("    -ref  <string>:  the reference sequence file used in single program, fasta format\n");
        printf("    -expect  <int>:  the number of LAMP primers needed to be design, default is 10\n"); 
        printf("    -common:         design common LAMP primers\n");
        printf("    -special:        design special LAMP primers\n");
        printf("    -check   <int>:  0: don't check tendency of the left primer to bind to the right primer; !=0: check, default is 1\n");
        printf("    -par  <string>:  the directory of storing parameter files used to check the tendency of two primers binding, default is Par/\n");
        printf("    -high/-low:      design candidate single primers in high/low GC region, high: the GC content>=45%%; low: the GC content <=45%%.\n");
        printf("    -loop:           design LAMP primer with loop primers\n");
        printf("    -h/-help:        usage\n");
}

struct INFO *read_list(char *path,int common_num[])  
{
        char *in,name[301];
        int turn,i,size;
        struct INFO *new_primer,*p_primer,*head;
        FILE *fp;

        i=strlen(path);
        in=(char *)malloc(i+20);              
        memset(in,'\0',i+20);
        strcpy(in,path);
        strcat(in,"-common_list.txt");
        if(access(in,0)==-1)  
        {
                printf("Error! Don't have the %s file!\n",in);
                exit(1);           
        }
        fp=fopen(in,"r");
        if(fp==NULL)
        {          
                printf("Error: can't open the %s file!\n",in);
                exit(1);
        }

        size=sizeof(struct INFO);
        i=0;
        memset(name,'\0',301);
        while(fscanf(fp,"%s\t%d\n",name,&turn)!=EOF)
        {
                new_primer=(struct INFO *)malloc(size);
                new_primer->turn=turn;
                strcpy(new_primer->name,name);
                new_primer->next=NULL;

                if(i==0)
                {
                        head=new_primer;
                        p_primer=new_primer;
                        i++;
                }
                else
                {
                        p_primer->next=new_primer;
                        p_primer=new_primer;
                }
                memset(name,'\0',301);
        }
        fclose(fp);
        common_num[0]=turn;
        free(in);
        return head;
}

main(int argc,char **argv)
{
	int i,j,flag[12],expect,circle,have,common_num[1],NumL[2],*d_NumL,num[11],max_loop,min_loop,count[3],block;
	char *output,*prefix,*store_path,*path_fa,*inner,*outer,*loop,*par_path,*Pchar,*d_Pchar,*temp,*seq,*d_seq,primer[26];
	FILE *fp;
	struct Primer *headL,*headS,*headLoop,*tempL,*tempS,*tempLoop,*storeL,*storeS,*storeLoop; 
	struct Node *p_node,*p_temp;
	struct INFO *headList,*p_list;
	time_t start,end;
	double *H_parameter,*parameter;	
	long int memory;
	cudaDeviceProp prop;
	int *d_S,*d_L,*d_Lp,*d_cS,*d_cL,*d_cLp,*d_scS,*d_ecS,*d_scL,*d_ecL,*d_scLp,*d_ecLp,*d_sS,*d_sL,*d_ssS,*d_esS,*d_ssL,*d_esL,*d_SS,*d_SL,*d_SLp,*d_LL,*d_LS,*d_LLp,*d_LpLp,*d_apply,*d_par,*d_int,*d_pos;
	int *h_S,*h_L,*h_Lp,*h_cS,*h_cL,*h_cLp,*h_scS,*h_ecS,*h_scL,*h_ecL,*h_scLp,*h_ecLp,*h_sS,*h_sL,*h_ssS,*h_esS,*h_ssL,*h_esL,*h_apply,*h_par,h_int[11],*h_pos;
	
	expect=10; //default output max 10 LAMP primers
	start=time(NULL);
/////read the parameters
        for(i=0;i<=11;i++)
                flag[i]=0;
	for(i=0;i<10;i++)
		h_int[i]=0;
        flag[7]=1;
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
                        j=strlen(argv[i+1]);
                        prefix=(char *)malloc(j+1);
                        memset(prefix,'\0',j+1);
                        strcpy(prefix,argv[i+1]);
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
                        j=strlen(argv[i+1]);
                        output=(char *)malloc(j+1);
                        memset(output,'\0',j+1);
                        strcpy(output,argv[i+1]);
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
                        j=strlen(argv[i+1]);
                        if(argv[i+1][j-1]=='/')
                        {
                                store_path=(char *)malloc(j+1);
                                memset(store_path,'\0',j+1);
                                strcpy(store_path,argv[i+1]);
                        }
                        else
                        {
                                store_path=(char *)malloc(j+2);
                                memset(store_path,'\0',j+2);
                                strcpy(store_path,argv[i+1]);
                                store_path[j]='/';
                        }
                        i=i+2;
                }
                else if(strcmp(argv[i],"-ref")==0)
                {
                        flag[3]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-ref\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
                        j=strlen(argv[i+1]);
                        path_fa=(char *)malloc(j+1);
                        memset(path_fa,'\0',j+1);
                        strcpy(path_fa,argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-expect")==0)
                {
                        flag[4]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-tm\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
                        expect=atoi(argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-high")==0)
                {
                        flag[8]=1;
                        i++;
                }
                else if(strcmp(argv[i],"-low")==0)
                {
                        flag[9]=1;
                        i++;
                }
                else if(strcmp(argv[i],"-loop")==0) 
                {
                        flag[10]=1;
			h_int[5]=1;
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
                        flag[7]=atoi(argv[i+1]);
                        i=i+2;
                }
                else if(strcmp(argv[i],"-par")==0)
                {
                        flag[11]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-par\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
                        j=strlen(argv[i+1]);
                        if(argv[i+1][j-1]=='/')
                        {
                                par_path=(char *)malloc(j+1);
                                strcpy(par_path,argv[i+1]);
                                par_path[j]='\0';
                        }
                        else
                        {
                                par_path=(char *)malloc(j+2);
                                strcpy(par_path,argv[i+1]);
                                par_path[j]='/';
                                par_path[j+1]='\0';
                        }
                        i=i+2;
                }
                else if(strcmp(argv[i],"-common")==0)
                {
                        flag[5]=1;
			h_int[3]=1;
                        i++;
                }
                else if(strcmp(argv[i],"-special")==0)
                {
                        flag[6]=1;
			h_int[4]=1;
                        i++;
                }
                else
                {
                        printf("Warning! The parameter of %s is invalid.\n\n",argv[i]);
                        i++;
                }
        }
//check parameters
        if(flag[0]==0)
        {
                printf("Error! Users must supply the name of candidate single primers file with -in!\n");
                usage();
                exit(1);
        }
        if(flag[1]==0)
        {
                printf("Error! Users must supply the name of output file with -out!\n");
                usage();
                exit(1);
        }
        if(flag[3]==0)
        {
                printf("Error! Users must supply the reference sequence file with -ref!\n");
                usage();
                exit(1);
        }
        if(flag[8]+flag[9]!=1)
        {
                printf("Error! The input parameter must contain one of -high and -low!\n");
                usage();
                exit(1);
        }
//prepare
	if(flag[7]!=0)
		h_int[6]=1;
	else
		h_int[6]=0;
	h_int[9]=flag[8];
	h_int[10]=expect;
        if(flag[2]==0)
        {
                temp=(char *)malloc(4096);
                memset(temp,'\0',4096);
                getcwd(temp,4096);
                j=strlen(temp);
                store_path=(char *)malloc(j+2);
                memset(store_path,'\0',j+2);
                strcpy(store_path,temp);
                free(temp);
                store_path[j]='/';
        }
//secondary
        if(flag[7]&&flag[11]==0)
        {
                temp=(char *)malloc(4096);
                memset(temp,'\0',4096);
                getcwd(temp,4096);
                j=strlen(temp);
                par_path=(char *)malloc(j+10);
                memset(par_path,'\0',j+10);
                strcpy(par_path,temp);
                free(temp);

                j--;
                while(par_path[j]!='/'&&j>=0)
                {
                        par_path[j]='\0';
                        j--;
                }
                strcat(par_path,"Par/");
        }
        if(flag[7])
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
//F3's pos 
	h_pos=(int *)malloc(expect*sizeof(int));
	for(i=0;i<expect;i++)
		h_pos[i]=-1;
	cudaMalloc((void **)&d_pos,expect*sizeof(int));
//d_int
	cudaMalloc((void **)&d_int,11*sizeof(int));
//directory for single primers
        j=strlen(store_path)+strlen(prefix)+12;
        outer=(char *)malloc(j);
        memset(outer,'\0',j);
        strcpy(outer,store_path);

        inner=(char *)malloc(j);
        memset(inner,'\0',j);
        strcpy(inner,store_path);

        if(flag[10]==1)
        {
                loop=(char *)malloc(j);
                memset(loop,'\0',j);
                strcpy(loop,store_path);
        }
        if(flag[8]==1)
        {                       
                strcat(outer,"high-outer/");
                strcat(outer,prefix);
                strcat(inner,"high-inner/");
                strcat(inner,prefix);
                if(flag[10]==1)
                {
                        strcat(loop,"high-loop/");
                        strcat(loop,prefix);
                }
        }
        else          
        {                
                strcat(outer,"low-outer/");
                strcat(outer,prefix);
                strcat(inner,"low-inner/");
                strcat(inner,prefix);
                if(flag[10]==1)
                {
                        strcat(loop,"low-loop/");
                        strcat(loop,prefix);
                }
        }

//reference sequence fa
        if(access(path_fa,0)==-1)
        {
                printf("Error! Don't have the %s file!\n",path_fa);
                exit(1);
        }
        i=file_size2(path_fa);
        i=i+100;
        temp=(char *)malloc(i*sizeof(char));
        memset(temp,'\0',i*sizeof(char));
        fp=fopen(path_fa,"r");
        if(fp==NULL)
        {
                printf("Error! Can't open the sequence file %s\n",path_fa);
                exit(1);
        }

        fread(temp,i*sizeof(char),1,fp);
        fclose(fp); 
        seq=(char *)malloc(i*sizeof(char));
        memset(seq,'\0',i*sizeof(char));
        
        j=0;
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
                        seq[j]='A';
                else if(temp[i]=='t'||temp[i]=='T')
                        seq[j]='T';
                else if(temp[i]=='c'||temp[i]=='C')
                        seq[j]='C';
                else if(temp[i]=='g'||temp[i]=='G')
                        seq[j]='G';
                else
                        seq[j]='N';
                i++;
                j++;
        }
        free(temp);
	num[0]=j; //the length of genome
	cudaMalloc((void **)&d_seq,num[0]);
	cudaMemset(d_seq,'\0',num[0]);
	cudaMemcpy(d_seq,seq,num[0],cudaMemcpyHostToDevice);

//common-list
        if(flag[5])
        {
                headList=read_list(inner,common_num);
                common_num[0]++;
        }
        else
                common_num[0]=1;
	h_int[7]=common_num[0];
//read parameters
        headS=read_par(outer,flag[5],flag[6]);
        headL=read_par(inner,flag[5],flag[6]);
        if(flag[10])
        {
                headLoop=read_par(loop,flag[5],0);
                tempLoop=headLoop;
                while(tempLoop->next!=NULL)
                        tempLoop=tempLoop->next;
                max_loop=tempLoop->pos;
		min_loop=headLoop->pos;
        }

//common statistics
        if(flag[5])
        {
                how_many(headL,common_num[0]);
		how_many(headS,common_num[0]);
                if(flag[10])
                        how_many(headLoop,common_num[0]);
        }

	cudaGetDeviceProperties(&prop,0); //read parameters
	fp=fopen(output,"w");
        if(fp==NULL)
        {
                printf("Error: can't create the %s file!\n",output);
                exit(1);
        }
	have=1;
        end=time(NULL);
        printf("The prepare time is %0.1f seconds!\n",difftime(end,start));
        start=time(NULL);

//LAMP-GPU
	for(circle=common_num[0];circle>=1;circle--)
	{
		if(have>expect)
			break;
		storeL=headL;
		storeS=headS;
		while((storeL->pos<=storeS->pos+18)&&storeL!=NULL)
			storeL=storeL->next;
		if(flag[10])
		{
                        storeLoop=headLoop;
			while((storeLoop->pos<=storeS->pos+18)&&storeLoop!=NULL)
				storeLoop=storeLoop->next;
		}
		num[10]=0;
		while(storeS)
		{
			if(have>expect)
				break;
			if(num[10]==1)  //don't have enough primers
				break;
			for(i=1;i<9;i++)
			{
				num[i]=0;
			}
			memory=num[0]/2;
		//statistics	
			tempL=storeL;
			tempS=storeS;
			if(flag[10])
				tempLoop=storeLoop;
			while(tempS&&(memory<prop.totalGlobalMem/6))
			{
				if(flag[10]&&(tempS->pos+200)<min_loop)
					continue;
				if(flag[10]&&(tempS->pos-200)>max_loop)
					break;
				if(tempS->total<circle)
				{
					tempS=tempS->next;
					continue;
				}
				while(tempL&&(tempL->pos<tempS->pos))
				{
					if(tempL->total<circle)
					{
						tempL=tempL->next;
						continue;
					}
					num[1]++;
					if(flag[5])
					{
						num[3]=num[3]+tempL->total_common;
						memory=memory+2+4*tempL->total_common; //2 is common_start and common_end
					}
					if(flag[6])
					{
						num[4]=num[4]+tempL->total_special;
						memory=memory+2+4*tempL->total_special;
					}
					memory=memory+6+flag[10]; //6=4(primers)+2(next_to),flag[10]: next_to_loop
					tempL=tempL->next;
				}
				
				while(flag[10]&&tempLoop&&(tempLoop->pos<tempS->pos))
                                {
                                        if(tempLoop->total<circle)
                                        {
                                                tempLoop=tempLoop->next;
                                                continue;
                                        }
                                        num[7]++;
                                        if(flag[5])
                                        {
                                                num[8]=num[8]+tempLoop->total_common;
                                                memory=memory+2+4*tempLoop->total_common; //2 is common_start and common_end
                                        }
                                        memory=memory+5; //5=4(primers)+1(next_to_self)
                                        tempLoop=tempLoop->next;
                                }

				num[2]++;
				if(flag[5])
				{
					num[5]=num[5]+tempS->total_common;
					memory=memory+2+4*tempS->total_common;
				}
				if(flag[6])
				{
					num[6]=num[6]+tempS->total_special;
					memory=memory+2+4*tempS->total_special;
				}
				memory=memory+22+flag[10]+common_num[0];//14=4(primers)+16(result_par,include loop)+2(next_to)
				tempS=tempS->next;
			}
			if(num[2]<4||num[1]<2||(flag[10]&&num[7]<1)) //don't have enough primers
			{
				num[10]=1;
				break;
			}
			if(tempS==NULL)  //check all primers
				num[10]=1;	

		//malloc small primer
			h_S=(int *)malloc(4*num[2]*sizeof(int));
			cudaMalloc((void **)&d_S,4*num[2]*sizeof(int));
			if(flag[5])
			{
				h_cS=(int *)malloc(4*num[5]*sizeof(int));
				cudaMalloc((void **)&d_cS,4*num[5]*sizeof(int));
				h_scS=(int *)malloc(num[2]*sizeof(int));
				cudaMalloc((void **)&d_scS,num[2]*sizeof(int));
	                        h_ecS=(int *)malloc(num[2]*sizeof(int));
	                        cudaMalloc((void **)&d_ecS,num[2]*sizeof(int));
			}
			if(flag[6])
			{
				h_sS=(int *)malloc(4*num[6]*sizeof(int));
				cudaMalloc((void **)&d_sS,4*num[6]*sizeof(int));
				h_ssS=(int *)malloc(num[2]*sizeof(int));
	                        cudaMalloc((void **)&d_ssS,num[2]*sizeof(int));
	                        h_esS=(int *)malloc(num[2]*sizeof(int));
	                        cudaMalloc((void **)&d_esS,num[2]*sizeof(int));
			}
			cudaMalloc((void **)&d_SS,num[2]*sizeof(int));
			cudaMalloc((void **)&d_SL,num[2]*sizeof(int));
			if(flag[10])
				cudaMalloc((void **)&d_SLp,num[2]*sizeof(int));
		
			tempS=storeS;
			for(i=0;i<3;i++)
				count[i]=0;
			while(count[0]<num[2])
			{
				if(tempS->total<circle)
				{
					tempS=tempS->next;
					continue;
				}
		//primer info
				h_S[4*count[0]]=tempS->pos;
				h_S[4*count[0]+1]=tempS->len;
				h_S[4*count[0]+2]=tempS->plus;
				h_S[4*count[0]+3]=tempS->minus;
			//common
				if(flag[5])
				{
					h_scS[count[0]]=count[1];
					if(tempS->total_common==0)
						h_ecS[count[0]]=-1;
					else
					{
						p_node=tempS->common;
						while(p_node)
						{
							h_cS[4*count[1]]=p_node->gi;
							h_cS[4*count[1]+1]=p_node->pos;
							h_cS[4*count[1]+2]=p_node->plus; 
                                		        h_cS[4*count[1]+3]=p_node->minus;
							count[1]++;
							p_node=p_node->next;
						}
						h_ecS[count[0]]=count[1];
					}
				}
			//special
				if(flag[6])
				{
					h_ssS[count[0]]=count[2];
                	        	if(tempS->total_special==0)
                	        	        h_esS[count[0]]=-1;
                	        	else
                	        	{
                	        	        p_node=tempS->special;
                	        	        while(p_node)
                	        	        {
                	        	                h_sS[4*count[2]]=p_node->gi;
                	        	                h_sS[4*count[2]+1]=p_node->pos;
                	        	                h_sS[4*count[2]+2]=p_node->plus;
                	        	                h_sS[4*count[2]+3]=p_node->minus;
                	        	                count[2]++;
                	        	                p_node=p_node->next;
                	        	        }
                	        	        h_esS[count[0]]=count[2];
                	        	}
				}
				count[0]++;
				tempS=tempS->next;
			}
		//copy
			h_int[0]=num[2];
			cudaMemcpy(d_S,h_S,num[2]*4*sizeof(int),cudaMemcpyHostToDevice);
			free(h_S);
			if(flag[5])
			{
				cudaMemcpy(d_cS,h_cS,num[5]*4*sizeof(int),cudaMemcpyHostToDevice);
				free(h_cS);
				cudaMemcpy(d_scS,h_scS,num[2]*sizeof(int),cudaMemcpyHostToDevice);
				free(h_scS);
				cudaMemcpy(d_ecS,h_ecS,num[2]*sizeof(int),cudaMemcpyHostToDevice);
        		        free(h_ecS);
			}
			if(flag[6])
			{
				cudaMemcpy(d_sS,h_sS,num[6]*4*sizeof(int),cudaMemcpyHostToDevice);
                        	free(h_sS);
				cudaMemcpy(d_ssS,h_ssS,num[2]*sizeof(int),cudaMemcpyHostToDevice);
        	        	free(h_ssS);
        	        	cudaMemcpy(d_esS,h_esS,num[2]*sizeof(int),cudaMemcpyHostToDevice);    
        	        	free(h_esS);
			}
		//large primer
			h_L=(int *)malloc(4*num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_L,4*num[1]*sizeof(int));
			if(flag[5])
			{
        	        	h_cL=(int *)malloc(4*num[3]*sizeof(int));
        	        	cudaMalloc((void **)&d_cL,4*num[3]*sizeof(int));
				h_scL=(int *)malloc(num[1]*sizeof(int));
                        	cudaMalloc((void **)&d_scL,num[1]*sizeof(int));
                        	h_ecL=(int *)malloc(num[1]*sizeof(int));
                        	cudaMalloc((void **)&d_ecL,num[1]*sizeof(int));
			}
			if(flag[6])
			{
        	        	h_sL=(int *)malloc(4*num[4]*sizeof(int));
 				cudaMalloc((void **)&d_sL,4*num[4]*sizeof(int));
				h_ssL=(int *)malloc(num[1]*sizeof(int));
                        	cudaMalloc((void **)&d_ssL,num[1]*sizeof(int));
                        	h_esL=(int *)malloc(num[1]*sizeof(int));
                        	cudaMalloc((void **)&d_esL,num[1]*sizeof(int));
			}
        	        cudaMalloc((void **)&d_LL,num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_LS,num[1]*sizeof(int));
			if(flag[10])
				cudaMalloc((void **)&d_LLp,num[1]*sizeof(int));

        	        tempL=storeL;
        	        for(i=0;i<3;i++)
        	                count[i]=0;
        	        while(count[0]<num[1])
        	        {
				if(tempL->total<circle)
				{
					tempL=tempL->next;
					continue;
				}
                	//primer info
                	        h_L[4*count[0]]=tempL->pos;
                	        h_L[4*count[0]+1]=tempL->len;
                	        h_L[4*count[0]+2]=tempL->plus;
                	        h_L[4*count[0]+3]=tempL->minus;
                	//common
				if(flag[5])
				{
                	        	h_scL[count[0]]=count[1];
                	        	if(tempL->total_common==0)
                	        	        h_ecL[count[0]]=-1;
                	        	else
                	        	{
                	        	        p_node=tempL->common;
                	        	        while(p_node)
                	        	        {
                	        	                h_cL[4*count[1]]=p_node->gi;
                	        	                h_cL[4*count[1]+1]=p_node->pos;
                	        	                h_cL[4*count[1]+2]=p_node->plus; 
                        		                h_cL[4*count[1]+3]=p_node->minus;
                        		                count[1]++;
                        		                p_node=p_node->next;
                        		        }
                        		        h_ecL[count[0]]=count[1];
					}
                        	}
                	//special
				if(flag[6])
				{
                        		h_ssL[count[0]]=count[2];
                        		if(tempL->total_special==0)
                        		        h_esL[count[0]]=-1;
                        		else
                        		{
                        		        p_node=tempL->special;
                        		        while(p_node)
                        		        {
                        		                h_sL[4*count[2]]=p_node->gi;
                        		                h_sL[4*count[2]+1]=p_node->pos;
                        		                h_sL[4*count[2]+2]=p_node->plus;
                        		                h_sL[4*count[2]+3]=p_node->minus;
                        		                count[2]++;
                        		                p_node=p_node->next;
                        		        }
                        		        h_esL[count[0]]=count[2];
                        		}
				}
                        	count[0]++;
                        	tempL=tempL->next;
                	}
        	//copy
			h_int[1]=num[1];
        	        cudaMemcpy(d_L,h_L,num[1]*4*sizeof(int),cudaMemcpyHostToDevice);
			free(h_L);
			if(flag[5])
			{
                		cudaMemcpy(d_cL,h_cL,num[3]*4*sizeof(int),cudaMemcpyHostToDevice);
                		free(h_cL);
	                	cudaMemcpy(d_scL,h_scL,num[1]*sizeof(int),cudaMemcpyHostToDevice);
	                	free(h_scL);
	                	cudaMemcpy(d_ecL,h_ecL,num[1]*sizeof(int),cudaMemcpyHostToDevice);
	                	free(h_ecL);
			}
			if(flag[6])
			{
				cudaMemcpy(d_sL,h_sL,num[4]*4*sizeof(int),cudaMemcpyHostToDevice);
                        	free(h_sL);
                		cudaMemcpy(d_ssL,h_ssL,num[1]*sizeof(int),cudaMemcpyHostToDevice);
                		free(h_ssL);
                		cudaMemcpy(d_esL,h_esL,num[1]*sizeof(int),cudaMemcpyHostToDevice);    
                		free(h_esL);
			}

                //loop primer
			if(flag[10])
			{
                        	h_Lp=(int *)malloc(4*num[7]*sizeof(int));
                        	cudaMalloc((void **)&d_Lp,4*num[7]*sizeof(int));
                        	if(flag[5])
                        	{
                        	        h_cLp=(int *)malloc(4*num[8]*sizeof(int));
                        	        cudaMalloc((void **)&d_cLp,4*num[8]*sizeof(int));
                        	        h_scLp=(int *)malloc(num[7]*sizeof(int));
                        	        cudaMalloc((void **)&d_scLp,num[7]*sizeof(int));
                        	        h_ecLp=(int *)malloc(num[7]*sizeof(int));
                        	        cudaMalloc((void **)&d_ecLp,num[7]*sizeof(int));
                        	}
	                        cudaMalloc((void **)&d_LpLp,num[7]*sizeof(int));

	                        tempLoop=storeLoop;
                        	for(i=0;i<2;i++)
                        	        count[i]=0;
	                        while(count[0]<num[7])
                        	{
                                	if(tempLoop->total<circle)
                                	{
                                        	tempLoop=tempLoop->next;
                                        	continue;
                                	}
                        	//primer info
                                	h_Lp[4*count[0]]=tempLoop->pos;
                                	h_Lp[4*count[0]+1]=tempLoop->len;
                                	h_Lp[4*count[0]+2]=tempLoop->plus;
                                	h_Lp[4*count[0]+3]=tempLoop->minus;
                        	//common
                                	if(flag[5])
                                	{
                                	        h_scLp[count[0]]=count[1];
                                	        if(tempLoop->total_common==0)
                                	                h_ecLp[count[0]]=-1;
                                	        else
                                	        {
                                	                p_node=tempLoop->common;
                                	                while(p_node)
                                	                {
                                	                        h_cLp[4*count[1]]=p_node->gi;
                                	                        h_cLp[4*count[1]+1]=p_node->pos;
                                	                        h_cLp[4*count[1]+2]=p_node->plus; 
                                	                        h_cLp[4*count[1]+3]=p_node->minus;
                                	                        count[1]++;
                                	                        p_node=p_node->next;
                                	                }
                                	                h_ecLp[count[0]]=count[1];
                                	        }
                                	}
	                                count[0]++;
	                                tempLoop=tempLoop->next;
	                        }
	                //copy
				h_int[2]=num[7];
                        	cudaMemcpy(d_Lp,h_Lp,num[7]*4*sizeof(int),cudaMemcpyHostToDevice);
                        	free(h_Lp);
                        	if(flag[5])
                        	{
                        	        cudaMemcpy(d_cLp,h_cLp,num[8]*4*sizeof(int),cudaMemcpyHostToDevice);
                        	        free(h_cLp);
                        	        cudaMemcpy(d_scLp,h_scLp,num[7]*sizeof(int),cudaMemcpyHostToDevice);
                        	        free(h_scLp);
                        	        cudaMemcpy(d_ecLp,h_ecLp,num[7]*sizeof(int),cudaMemcpyHostToDevice);
                        	        free(h_ecLp);
                        	}
                        }
	//run
			if(num[2]%prop.maxThreadsPerBlock==0)
                                block=num[2]/prop.maxThreadsPerBlock;
                        else
                                block=(num[2]-num[2]%prop.maxThreadsPerBlock)/prop.maxThreadsPerBlock+1;

			if(block>prop.maxGridSize[0])
				block=prop.maxGridSize[0];

		//next primer	
			next_one<<<block,prop.maxThreadsPerBlock>>>(d_S,d_S,d_SS,num[2],num[2]);
			next_one<<<block,prop.maxThreadsPerBlock>>>(d_L,d_L,d_LL,num[1],num[1]);
			next_one<<<block,prop.maxThreadsPerBlock>>>(d_S,d_L,d_SL,num[2],num[1]);
        	        next_one<<<block,prop.maxThreadsPerBlock>>>(d_L,d_S,d_LS,num[1],num[2]);
			if(flag[10])
			{
				next_one<<<block,prop.maxThreadsPerBlock>>>(d_Lp,d_Lp,d_LpLp,num[7],num[7]);
                        	next_one<<<block,prop.maxThreadsPerBlock>>>(d_S,d_Lp,d_SLp,num[2],num[7]);
                        	next_one<<<block,prop.maxThreadsPerBlock>>>(d_L,d_Lp,d_LLp,num[1],num[7]);
			}

		//calculate
			cudaMalloc((void **)&d_apply,common_num[0]*num[2]*sizeof(int));
			cudaMemcpy(d_pos,h_pos,expect*sizeof(int),cudaMemcpyHostToDevice);
			h_int[8]=circle;
			cudaMemcpy(d_int,h_int,11*sizeof(int),cudaMemcpyHostToDevice);
			cudaMemcpy(d_pos,h_pos,expect*sizeof(int),cudaMemcpyHostToDevice);
			cudaMalloc((void **)&d_par,num[2]*16*sizeof(int));
			LAMP<<<block,prop.maxThreadsPerBlock>>>(d_seq,d_S,d_L,d_Lp,d_cS,d_cL,d_cLp,d_sS,d_sL,d_scS,d_scL,d_scLp,d_ecS,d_ecL,d_ecLp,d_ssS,d_ssL,d_esS,d_esL,d_SS,d_SL,d_SLp,d_LL,d_LS,d_LLp,d_LpLp,d_int,d_par,d_apply,parameter,d_Pchar,d_NumL,d_pos);

			h_par=(int *)malloc(16*num[2]*sizeof(int));
			cudaMemcpy(h_par,d_par,16*num[2]*sizeof(int),cudaMemcpyDeviceToHost);
			cudaFree(d_par);
			h_apply=(int *)malloc(common_num[0]*num[2]*sizeof(int));
			cudaMemcpy(h_apply,d_apply,common_num[0]*num[2]*sizeof(int),cudaMemcpyDeviceToHost);
			cudaFree(d_apply);
	//free
			cudaFree(d_S);
			cudaFree(d_L);
			cudaFree(d_SS);
			cudaFree(d_LL);
			cudaFree(d_SL);
			cudaFree(d_LS);
			if(flag[10])
			{
				cudaFree(d_Lp);
				cudaFree(d_LpLp);
				cudaFree(d_SLp);
				cudaFree(d_LLp);
				if(flag[5])
				{
					cudaFree(d_cLp);
					cudaFree(d_scLp);
					cudaFree(d_ecLp);
				}
			}
			if(flag[5])
			{
				cudaFree(d_cS);
				cudaFree(d_scS);
				cudaFree(d_ecS);
				cudaFree(d_cL);
				cudaFree(d_scL);
				cudaFree(d_ecL);
			}
			if(flag[6])
			{
                                cudaFree(d_sS);
                                cudaFree(d_ssS);
                                cudaFree(d_esS);
                                cudaFree(d_sL);
                                cudaFree(d_ssL);
                                cudaFree(d_esL);
                        }
		//LAMP primers, output
			for(i=0;i<num[2];i++)
			{
				if(have>expect)
					break;
				if(h_par[i*16+1]==0)
					continue;
				if(check_add(h_par[i*16],h_pos,have)==0)
					continue;
				fprintf(fp,"The %d LAMP primers:\n",have);
		                generate_primer(seq,primer,h_par[i*16],h_par[i*16+1],0);
		                fprintf(fp,"  F3: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16],h_par[i*16+1],primer);
		                generate_primer(seq,primer,h_par[i*16+2],h_par[i*16+3],0);
		                fprintf(fp,"  F2: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16+2],h_par[i*16+3],primer);
		                generate_primer(seq,primer,h_par[i*16+4],h_par[i*16+5],1);
		                fprintf(fp,"  F1c: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16+4],h_par[i*16+5],primer);
		                generate_primer(seq,primer,h_par[i*16+6],h_par[i*16+7],0);
		                fprintf(fp,"  B1c: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16+6],h_par[i*16+7],primer);
		                generate_primer(seq,primer,h_par[i*16+8],h_par[i*16+9],1);
		                fprintf(fp,"  B2: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16+8],h_par[i*16+9],primer);
		                generate_primer(seq,primer,h_par[i*16+10],h_par[i*16+11],1);
		                fprintf(fp,"  B3: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16+10],h_par[i*16+11],primer);
                		if(flag[10])
                		{
                        		if(h_par[i*16+13]==0)
                        		        fprintf(fp,"  LF: NULL\n");
                        		else
                        		{
                        		        generate_primer(seq,primer,h_par[i*16+12],h_par[i*16+13],1);
                        		        fprintf(fp,"  LF: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16+12],h_par[i*16+13],primer);
                        		}
		                        if(h_par[i*16+15]==0)
		                                fprintf(fp,"  LB: NULL\n");
		                        else
		                        {
		                                generate_primer(seq,primer,h_par[i*16+14],h_par[i*16+15],0);
		                                fprintf(fp,"  LB: pos:%d,length:%d bp, primer(5'-3'):%s\n",h_par[i*16+14],h_par[i*16+15],primer);
		                        }
		                }
		                if(flag[5])
		                {
					j=0;
		                        fprintf(fp,"  This set of LAMP primers could be used in %d genomes, there are: ",circle);
		                        p_list=headList;
                        		for(j=0;j<common_num[0];j++)
                        		{
                                		if(h_apply[i*common_num[0]+j]==0)
                                		        continue;
                                		while(p_list)
                                		{
                                		        if(p_list->turn==j)
                                		                break;
                                		        else
                                		                p_list=p_list->next;
                                		}
                                		if(j==0)
                                		        fprintf(fp,"%s",p_list->name);
                                		else
                                		        fprintf(fp,", %s",p_list->name);
                                		j++;
	       	                 	}
		                        fprintf(fp,"\n");
				}
				h_pos[have-1]=h_par[i*16];
				have++;
			}
			free(h_apply);
			free(h_par);
			if(have>expect)
				break;
		//new primer start
			if(tempS==NULL)
				storeS=tempS;
			else
			{
				while(tempS->pos-storeS->pos>300)
				{
					storeS=storeS->next;
				}
			}
		}//one circle
	}
	cudaFree(d_int);
	cudaFree(d_pos);
	cudaFree(d_seq);
	free(h_pos);
	free(seq);
	free(output);
        free(prefix);
        free(store_path);
        free(inner);
        free(outer);
        free(path_fa);
	fclose(fp);
//free struct list
        while(headL)
        {
                p_node=headL->common;
                while(p_node)
                {
                        p_temp=p_node->next;
                        free(p_node);
                        p_node=p_temp;
                }
                p_node=headL->special;
                while(p_node)  
                {
                        p_temp=p_node->next;  
                        free(p_node);  
                        p_node=p_temp;  
                }
                
                storeL=headL->next;
                free(headL);
                headL=storeL;
        }
        while(headS)
        {
                p_node=headS->common;  
                while(p_node)  
                {
                        p_temp=p_node->next;  
                        free(p_node);  
                        p_node=p_temp;  
                }
                p_node=headS->special;
                while(p_node)
                {               
                        p_temp=p_node->next;
                        free(p_node);
                        p_node=p_temp;
                }

                storeS=headS->next;
                free(headS);
                headS=storeS;
        }

        if(flag[5])
        {
                while(headList)
                {
                        p_list=headList->next;
                        free(headList);
                        headList=p_list;
                }
        }

        if(flag[7])
        {
		free(Pchar);
                free(H_parameter);
                cudaFree(parameter);
                cudaFree(d_Pchar);
                cudaFree(d_NumL);
        }
        if(flag[7]||flag[11])
                free(par_path);

        if(flag[10])
        {
                free(loop);
                while(headLoop)
                {
                        p_node=headLoop->common;
                        while(p_node)
                        {
                                p_temp=p_node->next;
                                free(p_node);
                                p_node=p_temp;
                        }
                        p_node=headLoop->special;
                        while(p_node)
                        {
                                p_temp=p_node->next;
                                free(p_node);
                                p_node=p_temp;
                        }

                        storeLoop=headLoop->next;
                        free(headLoop);
                        headLoop=storeLoop;
                }
        }	
	end=time(NULL);
        printf("the time for design is %0.1f seconds!\n",difftime(end,start));
}