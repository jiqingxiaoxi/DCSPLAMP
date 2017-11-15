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

__device__ char str2int_rev(char c)
{
        switch (c)
        {
                case 'T':
                        return 0;
                case 'G':
                        return 1;
                case 'C':
                        return 2;                 
                case 'A':               
                        return 3;         
        }
        return 4;
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

__device__ int equal(double a,double b)
{
	if(fabs(a)>999999999||fabs(b)>999999999)
		return 0;
	return fabs(a-b)<1e-5;
}

__device__ void initMatrix(int Initint[],double *d_DPT,int id,char numSeq1[],char numSeq2[])
{
	int i,j;

	for(i=1;i<=Initint[0];++i)
	{
		for(j=1;j<=Initint[1];++j)
		{
			if(numSeq1[i]+numSeq2[j]!=3)
			{
				d_DPT[id*1250+(i-1)*Initint[2]+j-1]=1.0*INFINITY;
				d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]=-1.0;
			}
			else
			{
				d_DPT[id*1250+(i-1)*Initint[2]+j-1]=0.0;
				d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]=-3224.0;
			}
		}
	}
}

__device__ void LSH(int i,int j,double *EntropyEnthalpy,double Initdouble[],int Initint[],double *d_DPT,int id,char numSeq1[],char numSeq2[],double parameter[])
{
	double S1,H1,T1,S2,H2,T2;

	if(numSeq1[i]+numSeq2[j]!=3)
	{
		d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]=-1.0;
		d_DPT[id*1250+(i-1)*Initint[2]+j-1]=1.0*INFINITY;
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

__device__ void maxTM(int i,int j,double Initdouble[],int Initint[],double *d_DPT,int id,char numSeq1[],char numSeq2[],double parameter[])
{
	double T0,T1,S0,S1,H0,H1;

	S0=d_DPT[id*1250+625+(i-1)*Initint[2]+j-1];
	H0=d_DPT[id*1250+(i-1)*Initint[2]+j-1];
	T0=(H0+Initdouble[0])/(S0+Initdouble[1]+Initdouble[2]); // at current position 
	if(fabs(d_DPT[id*1250+(i-2)*Initint[2]+j-2])<999999999&&fabs(Hs(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter))<999999999)
	{
		S1=(d_DPT[id*1250+625+(i-2)*Initint[2]+j-2]+Ss(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter));
		H1=(d_DPT[id*1250+(i-2)*Initint[2]+j-2]+Hs(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter));
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
		d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]=S1;
		d_DPT[id*1250+(i-1)*Initint[2]+j-1]=H1;
	}
	else if(T0>=T1)
	{
		d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]=S0;
		d_DPT[id*1250+(i-1)*Initint[2]+j-1]=H0;
	}
}

__device__ void calc_bulge_internal(int i,int j,int ii,int jj,double* EntropyEnthalpy,int traceback,double Initdouble[],int Initint[],double *d_DPT,int id,char numSeq1[],char numSeq2[],double parameter[])
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
			H+=d_DPT[id*1250+(i-1)*Initint[2]+j-1];
			S+=d_DPT[id*1250+625+(i-1)*Initint[2]+j-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}

			T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
			T2=(d_DPT[id*1250+(ii-1)*Initint[2]+jj-1]+Initdouble[0])/((d_DPT[id*1250+625+(ii-1)*Initint[2]+jj-1])+Initdouble[1]+Initdouble[2]);
			if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
			{
				EntropyEnthalpy[0]=S;
				EntropyEnthalpy[1]=H;
			}
		}
		else // we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30
		{
			H=parameter[3150+loopSize]+parameter[5705+numSeq1[i]*5+numSeq2[j]]+parameter[5705+numSeq1[ii]*5+numSeq2[jj]];
			H+=d_DPT[id*1250+(i-1)*Initint[2]+j-1];

			S=parameter[3060+loopSize]+parameter[5680+numSeq1[i]*5+numSeq2[j]]+parameter[5680+numSeq1[ii]*5+numSeq2[jj]];
			S+=d_DPT[id*1250+625+(i-1)*Initint[2]+j-1];
			if(fabs(H)>999999999)
			{
				H=1.0*INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
			T2=(d_DPT[id*1250+(ii-1)*Initint[2]+jj-1]+Initdouble[0])/(d_DPT[id*1250+625+(ii-1)*Initint[2]+jj-1]+Initdouble[1]+Initdouble[2]);
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
		S+=d_DPT[id*1250+625+(i-1)*Initint[2]+j-1];

		H=parameter[1875+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+parameter[1875+numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		H+=d_DPT[id*1250+(i-1)*Initint[2]+j-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
		T2=(d_DPT[id*1250+(ii-1)*Initint[2]+jj-1]+Initdouble[0])/(d_DPT[id*1250+625+(ii-1)*Initint[2]+jj-1]+Initdouble[1]+Initdouble[2]);
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
		H+=d_DPT[id*1250+(i-1)*Initint[2]+j-1];

		S=parameter[3030+loopSize]+parameter[3180+numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+parameter[3180+numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
		S+=d_DPT[id*1250+625+(i-1)*Initint[2]+j-1];
		if(fabs(H)>999999999)
		{
			H=1.0*INFINITY;
			S=-1.0;
		}
		T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
		T2=(d_DPT[id*1250+(ii-1)*Initint[2]+jj-1]+Initdouble[0])/((d_DPT[id*1250+625+(ii-1)*Initint[2]+jj-1])+Initdouble[1]+Initdouble[2]);
		if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
		{
			EntropyEnthalpy[0]=S;
			EntropyEnthalpy[1]=H;
		}
	}
	return;
}

__device__ void fillMatrix(double Initdouble[],int Initint[],double *d_DPT,int id,char numSeq1[],char numSeq2[],double *parameter)
{
	int d,i,j,ii,jj;
	double SH[2];

	for(i=1;i<=Initint[0];++i)
	{
		for(j=1;j<=Initint[1];++j)
		{
			if(fabs(d_DPT[id*1250+(i-1)*Initint[2]+j-1])<999999999)
			{
				SH[0]=-1.0;
				SH[1]=1.0*INFINITY;
				LSH(i,j,SH,Initdouble,Initint,d_DPT,id,numSeq1,numSeq2,parameter);

				if(fabs(SH[1])<999999999)
				{
					d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]=SH[0];
					d_DPT[id*1250+(i-1)*Initint[2]+j-1]=SH[1];
				}
				if(i>1&&j>1)
				{
					maxTM(i,j,Initdouble,Initint,d_DPT,id,numSeq1,numSeq2,parameter);
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
							if(fabs(d_DPT[id*1250+(ii-1)*Initint[2]+jj-1])<999999999)
							{
								SH[0]=-1.0;
								SH[1]=1.0*INFINITY;
								calc_bulge_internal(ii,jj,i,j,SH,0,Initdouble,Initint,d_DPT,id,numSeq1,numSeq2,parameter);

								if(SH[0]<-2500.0)
								{
									SH[0] =-3224.0;
									SH[1] = 0.0;
								}
								if(fabs(SH[1])<999999999)
								{
									d_DPT[id*1250+(i-1)*Initint[2]+j-1]=SH[1];
									d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]=SH[0];
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

__device__ void traceback(int i,int j,int *d_ps,double Initdouble[],int Initint[],double *d_DPT,int id,char numSeq1[],char numSeq2[],double *parameter)
{
	int d,ii,jj,done;
	double SH[2];

	d_ps[id*50+i-1]=j;
	d_ps[id*50+25+j-1]=i;
	while(1)
	{
		SH[0]=-1.0;
		SH[1]=1.0*INFINITY;
		LSH(i,j,SH,Initdouble,Initint,d_DPT,id,numSeq1,numSeq2,parameter);
		if(equal(d_DPT[id*1250+625+(i-1)*Initint[2]+j-1],SH[0])&&equal(d_DPT[id*1250+(i-1)*Initint[2]+j-1],SH[1]))
			break;

		done = 0;
		if(i>1&&j>1&&equal(d_DPT[id*1250+625+(i-1)*Initint[2]+j-1],Ss(i-1,j-1,1,Initint,numSeq1,numSeq2,parameter)+d_DPT[id*1250+625+(i-2)*Initint[2]+j-2]))
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
				calc_bulge_internal(ii,jj,i,j,SH,1,Initdouble,Initint,d_DPT,id,numSeq1,numSeq2,parameter);
				if(equal(d_DPT[id*1250+625+(i-1)*Initint[2]+j-1],SH[0])&&equal(d_DPT[1250*id+(i-1)*Initint[2]+j-1],SH[1]))
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

__device__ double drawDimer(int *d_ps,int id,double H,double S,double Initdouble[],int Initint[])
{
        int i,N;

        if(fabs(Initdouble[3])>999999999)
                return (double)0.0;
        else
        {
                N=0;
                for(i=0;i<Initint[0];i++)
                {
                        if(d_ps[id*50+i]>0)
                                ++N;
                }
                for(i=0;i<Initint[1];i++)
                {
                        if(d_ps[id*50+25+i]>0)
                                ++N;
                }
                N=(N/2)-1;
                return (double)(H/(S+(N*-0.51986)+Initdouble[2])-273.15);
        }
}

__device__ int symmetry_thermo(char *d_seq,int start,int length )
{
	int i = 0;
	if(length%2==1)
		return 0;

	while(i<length/2)
	{
		if((d_seq[i+start]=='A'&&d_seq[start+length-1-i]!='T')||(d_seq[i+start]=='T'&&d_seq[start+length-1-i]!='A')||(d_seq[start+length-1-i]=='A'&&d_seq[i+start]!='T')||(d_seq[start+length-1-i]=='T'&&d_seq[i+start]!='A'))
			return 0;
		if((d_seq[i+start]=='C'&&d_seq[start+length-1-i]!='G')||(d_seq[i+start]=='G'&&d_seq[start+length-1-i]!='C')||(d_seq[start+length-1-i]=='C'&&d_seq[i+start]!='G')||(d_seq[start+length-1-i]=='G'&&d_seq[i+start]!='C'))
			return 0;
		i++;
	}
	return 1;
}

__device__ double thal(char *d_seq,int *d_primer,int one_turn,int two_turn,int one_flag,int two_flag,int type,double *parameter,double *d_DPT,int id,int *d_ps)
{
	double SH[2],Initdouble[4];//0 is dplx_init_H, 1 is dplx_init_S, 2 is RC, 3 is SHleft
	int Initint[5]; //0 is len1, 1 is len2, 2 is len3, 3 is bestI, 4 is bestJ
	int i, j;
	double T1,result_TH;
	char numSeq1[27],numSeq2[27];

/*** INIT values for unimolecular and bimolecular structures ***/
	Initdouble[0]= 200;
	Initdouble[1]= -5.7;
	if(symmetry_thermo(d_seq,d_primer[4*one_turn],d_primer[4*one_turn+1])&&symmetry_thermo(d_seq,d_primer[4*two_turn],d_primer[4*two_turn+1]))
		Initdouble[2]=1.9872* log(38/1000000000.0);
	else
		Initdouble[2]=1.9872* log(38/4000000000.0);
/* convert nucleotides to numbers */
	if(type==1 || type==2)
	{
		Initint[0]=d_primer[4*one_turn+1];
		Initint[1]=d_primer[4*two_turn+1];
		if(one_flag==0) //plus
		{
	 		for(i=1;i<=Initint[0];++i)
				numSeq1[i]=str2int(d_seq[d_primer[4*one_turn]+i-1]);
		}
		else
		{
			for(i=1;i<=Initint[0];++i)
				numSeq1[i]=str2int_rev(d_seq[d_primer[4*one_turn]+d_primer[4*one_turn+1]-i]);
		}

		if(two_flag==0)
		{
			for(i=1;i<=Initint[1];++i)
				numSeq2[i]=str2int(d_seq[d_primer[4*two_turn]+d_primer[4*two_turn+1]-i]);
		}
		else
		{
			for(i=1;i<=Initint[1];++i)
				numSeq2[i]=str2int_rev(d_seq[d_primer[4*two_turn]+i-1]);
		}
	}
	else if(type==3)
	{
		Initint[0]=d_primer[4*two_turn+1];
		Initint[1]=d_primer[4*one_turn+1];
		if(two_flag==0)
		{
			for(i=1;i<=Initint[0];++i)
				numSeq1[i]=str2int(d_seq[d_primer[4*two_turn]+i-1]);
		}
		else
		{
			for(i=1;i<=Initint[0];++i)
				numSeq1[i]=str2int_rev(d_seq[d_primer[4*two_turn]+d_primer[4*two_turn+1]-i]);
		}
		if(one_flag==0)
		{
			for(i=1;i<=Initint[1];++i)
				numSeq2[i]=str2int(d_seq[d_primer[4*one_turn]+d_primer[4*one_turn+1]-i]);
		}
		else
		{
			for(i=1;i<=Initint[1];++i)
				numSeq2[i]=str2int_rev(d_seq[d_primer[4*one_turn]+i-1]);
		}
	}
	numSeq1[0]=numSeq1[Initint[0]+1]=numSeq2[0]=numSeq2[Initint[1]+1]=4; /* mark as N-s */

	result_TH=0;
	Initint[2]=Initint[1];
	initMatrix(Initint,d_DPT,id,numSeq1,numSeq2);
	fillMatrix(Initdouble,Initint,d_DPT,id,numSeq1,numSeq2,parameter);

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
				T1=((d_DPT[id*1250+(i-1)*Initint[2]+j-1]+ SH[1] +Initdouble[0]) / ((d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]) + SH[0] +Initdouble[1] + Initdouble[2])) -273.15;
				if(T1>Initdouble[3]&&((d_DPT[id*1250+625+(i-1)*Initint[2]+j-1]+SH[0])<0&&(SH[1]+d_DPT[id*1250+(i-1)*Initint[2]+j-1])<0))
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
			T1=((d_DPT[id*1250+(i-1)*Initint[2]+j-1]+SH[1]+Initdouble[0])/((d_DPT[id*1250+625+(i-1)*Initint[2]+j-1])+SH[0]+Initdouble[1]+Initdouble[2]))-273.15;
			if (T1>Initdouble[3]&&((SH[0]+d_DPT[id*1250+625+(i-1)*Initint[2]+j-1])<0&&(SH[1]+d_DPT[id*1250+(i-1)*Initint[2]+j-1])<0))
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
		d_ps[id*50+i]=0;
	for (j=0;j<Initint[1];++j)
		d_ps[id*50+25+j] = 0;
	if(fabs(d_DPT[id*1250+(Initint[3]-1)*Initint[2]+Initint[4]-1])<999999999)
	{
		traceback(Initint[3],Initint[4],d_ps,Initdouble,Initint,d_DPT,id,numSeq1,numSeq2,parameter);
		result_TH=drawDimer(d_ps,id,(d_DPT[id*1250+(Initint[3]-1)*Initint[2]+Initint[4]-1]+SH[1]+Initdouble[0]),(d_DPT[id*1250+625+(Initint[3]-1)*Initint[2]+Initint[4]-1]+SH[0]+Initdouble[1]),Initdouble,Initint);
		result_TH=(int)(100*result_TH+0.5)/100.0;
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
		if(par[i]==-1)
			return 1;
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

__device__ int check_structure(char *d_seq,int *d_primer,int turn[],int *d_int,double *parameter,double *d_TH,int id,double *d_DPT,int *d_ps)
{
	double TH;
	int i,j;

	for(i=0;i<5;i++)
	{
		for(j=i+1;j<6;j++)
		{
		if(i!=2||j!=3)
			continue;
			TH=thal(d_seq,d_primer,turn[i],turn[j],d_int[11+i],d_int[11+j],1,parameter,d_DPT,id,d_ps);
			if(TH>44+5*d_int[9])
                                return 0;
		d_TH[id*2]=TH;
			TH=thal(d_seq,d_primer,turn[i],turn[j],d_int[11+i],d_int[11+j],2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
		d_TH[id*2+1]=TH;
			TH=thal(d_seq,d_primer,turn[i],turn[j],d_int[11+i],d_int[11+j],3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
		if(TH>d_TH[id*2+1])
			d_TH[id*2+1]=TH;
			TH=thal(d_seq,d_primer,turn[j],turn[i],(1-d_int[11+j]),(1-d_int[11+i]),2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
		if(TH>d_TH[id*2+1])
                        d_TH[id*2+1]=TH;
                        TH=thal(d_seq,d_primer,turn[j],turn[i],(1-d_int[11+j]),(1-d_int[11+i]),3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
		if(TH>d_TH[id*2+1])
                        d_TH[id*2+1]=TH;
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
__device__ int check_common(int *d_primer,int *d_common,int *d_sc,int *d_ec,int turn[],int common,int *d_apply)
{
        int pos[6],i,dis;

        for(i=0;i<common;i++)
        {
                d_apply[common*turn[0]+i]=0;
        }
//plus
        for(pos[0]=d_sc[turn[0]];pos[0]<d_ec[turn[0]];pos[0]++)
        {
                if(d_common[4*pos[0]+2]!=1)
                        continue;
		i=d_common[4*pos[0]];
                if(d_apply[common*turn[0]+i]==1)
                        continue;
                for(pos[1]=d_sc[turn[1]];pos[1]<d_ec[turn[1]];pos[1]++)
                {
                        if(d_common[4*pos[1]]!=i)
                                continue;
                        if(d_common[4*pos[1]+2]!=1)
                                continue;
                        for(pos[2]=d_sc[turn[2]];pos[2]<d_ec[turn[2]];pos[2]++)
                        {
                                if(d_common[4*pos[2]]!=i)
                                        continue;
                                if(d_common[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=d_sc[turn[3]];pos[3]<d_ec[turn[3]];pos[3]++)
                                {
                                        if(d_common[4*pos[3]]!=i)
                                                continue;
                                        if(d_common[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=d_sc[turn[4]];pos[4]<d_ec[turn[4]];pos[4]++)
                                        {
                                                if(d_common[4*pos[4]]!=i)
                                                        continue;
                                                if(d_common[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=d_sc[turn[5]];pos[5]<d_ec[turn[5]];pos[5]++)
                                                {
                                                        if(d_common[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_common[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_common[4*pos[1]+1]-(d_common[4*pos[0]+1]+d_primer[4*turn[0]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=d_common[4*pos[2]+1]-d_common[4*pos[1]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_common[4*pos[3]+1]-(d_common[4*pos[2]+1]+d_primer[4*turn[2]+1]);
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=(d_common[4*pos[4]+1]+d_primer[4*turn[4]+1]-1)-(d_common[4*pos[3]+1]+d_primer[4*turn[3]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_common[4*pos[4]+1]+d_primer[4*turn[4]+1]-1-d_common[4*pos[1]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_common[4*pos[5]+1]-(d_common[4*pos[4]+1]+d_primer[4*turn[4]+1]);
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
        for(pos[0]=d_sc[turn[0]];pos[0]<d_ec[turn[0]];pos[0]++)
        {
                if(d_common[4*pos[0]+3]!=1)
                        continue;
		i=d_common[4*pos[0]];
                if(d_apply[turn[0]*common+i]==1)
                        continue;  //this GI can common

                for(pos[1]=d_sc[turn[1]];pos[1]<d_ec[turn[1]];pos[1]++)
                {
                        if(d_common[4*pos[1]]!=i)
                                continue;
                        if(d_common[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=d_sc[turn[2]];pos[2]<d_ec[turn[2]];pos[2]++)
                        {
                                if(d_common[4*pos[2]]!=i)
                                        continue;
                                if(d_common[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=d_sc[turn[3]];pos[3]<d_ec[turn[3]];pos[3]++)
                                {
                                        if(d_common[4*pos[3]]!=i)
                                                continue;
                                        if(d_common[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=d_sc[turn[4]];pos[4]<d_ec[turn[4]];pos[4]++)
                                        {
                                                if(d_common[4*pos[4]]!=i)
                                                        continue;
                                                if(d_common[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=d_sc[turn[5]];pos[5]<d_ec[turn[5]];pos[5]++)
                                                {
                                                        if(d_common[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_common[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_common[4*pos[0]+1]-(d_common[4*pos[1]+1]+d_primer[4*turn[1]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=(d_common[4*pos[1]+1]+d_primer[4*turn[1]+1]-1)-(d_common[4*pos[2]+1]+d_primer[4*turn[2]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_common[4*pos[2]+1]-(d_common[4*pos[3]+1]+d_primer[4*turn[3]+1]);
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=d_common[4*pos[3]+1]-d_common[4*pos[4]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_common[4*pos[1]+1]+d_primer[4*turn[1]+1]-1-d_common[4*pos[4]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_common[4*pos[4]+1]-(d_common[4*pos[5]+1]+d_primer[4*turn[5]+1]);
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
__device__ int check_uniq(int *d_primer,int *d_special,int *d_ss,int *d_es,int turn[])
{
        int pos[6],gi;

//plus
        for(pos[0]=d_ss[turn[0]];pos[0]<d_es[turn[0]];pos[0]++)
        {
                if(d_special[4*pos[0]+2]!=1)
                        continue;
		gi=d_special[4*pos[0]];
                for(pos[1]=d_ss[turn[1]];pos[1]<d_es[turn[1]];pos[1]++)
                {
			if(d_special[4*pos[1]]!=gi)
                                continue;
                        if(d_special[4*pos[1]+2]!=1)
				continue;
                        for(pos[2]=d_ss[turn[2]];pos[2]<d_es[turn[2]];pos[2]++) //F1c
                        {
                                if(d_special[4*pos[2]]!=gi)
                                        continue;
                                if(d_special[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=d_ss[turn[3]];pos[3]<d_es[turn[3]];pos[3]++) //B1c
                                {
                                        if(d_special[pos[3]*4]!=gi)
                                                continue;
                                        if(d_special[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=d_ss[turn[4]];pos[4]<d_es[turn[4]];pos[4]++) //B2
                                        {
                                                if(d_special[4*pos[4]]!=gi)
                                                        continue;
                                                if(d_special[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=d_ss[turn[5]];pos[5]<d_es[turn[5]];pos[5]++)
                                                {
                                                        if(d_special[4*pos[5]]!=gi)
                                                                continue;
                                                        if(d_special[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(d_special[4*pos[1]+1]<d_special[4*pos[0]+1])
                                                                continue;
                                                //F2-F1c
                                                        if(d_special[4*pos[2]+1]<d_special[4*pos[1]+1]+d_primer[4*turn[1]+1])
                                                                continue;
                                                //F1c-B1c
                                                        if(d_special[4*pos[3]+1]<d_special[4*pos[2]+1]+d_primer[4*turn[2]+1])
                                                                continue;
                                                //B1c-B2
                                                        if(d_special[4*pos[4]+1]<d_special[4*pos[3]+1]+d_primer[4*turn[3]+1])
                                                                continue;
                                                //B2-B3
                                                        if(d_special[4*pos[5]+1]<d_special[4*pos[4]+1])
                                                                continue;
                                                //whole
                                                        if(d_special[4*pos[5]+1]-d_special[4*pos[0]+1]>1000)
                                                                continue;
                                                        return 0;
                                                }//B3
                                        }
                                }//B1c
                        }
                }//F2
        }

//minus
        for(pos[0]=d_ss[turn[0]];pos[0]<d_es[turn[0]];pos[0]++)
        {
                if(d_special[4*pos[0]+3]!=1)
                        continue;
		gi=d_special[4*pos[0]];
                for(pos[1]=d_ss[turn[1]];pos[1]<d_es[turn[1]];pos[1]++)
                {
                        if(d_special[4*pos[1]]!=gi)
                                continue;
                        if(d_special[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=d_ss[turn[2]];pos[2]<d_es[turn[2]];pos[2]++)
                        {
                                if(d_special[4*pos[2]]!=gi)
                                        continue;
                                if(d_special[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=d_ss[turn[3]];pos[3]<d_es[turn[3]];pos[3]++)
                                {
                                        if(d_special[4*pos[3]]!=gi)
                                                continue;
                                        if(d_special[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=d_ss[turn[4]];pos[4]<d_es[turn[4]];pos[4]++)
                                        {
                                                if(d_special[4*pos[4]]!=gi)
                                                        continue;
                                                if(d_special[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=d_ss[turn[5]];pos[5]<d_es[turn[5]];pos[5]++)
                                                {
                                                        if(d_special[4*pos[5]]!=gi)
                                                                continue;
                                                        if(d_special[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(d_special[4*pos[0]+1]<d_special[4*pos[1]+1])
                                                                continue;
                                                //F2-F1c
                                                        if(d_special[4*pos[1]+1]<d_special[4*pos[2]+1]+d_primer[4*turn[2]+1])
                                                                continue;
                                                //F1c-B1c
                                                        if(d_special[4*pos[2]+1]<d_special[4*pos[3]+1]+d_primer[4*turn[3]+1])
                                                                continue;
                                                //B1c-B2
                                                        if(d_special[4*pos[3]+1]<d_special[4*pos[4]+1]+d_primer[4*turn[4]+1])
                                                                continue;
                                                //B2-B3
                                                        if(d_special[4*pos[4]+1]<d_special[4*pos[5]+1])
                                                                continue;
                                                //whole
                                                        if(d_special[4*pos[0]+1]-d_special[4*pos[5]+1]>1000)
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
__global__ void next_one(int *d_primer,int *next,int one_start,int one_end,int two_start,int two_end)
{
        int id=blockDim.x*blockIdx.x+threadIdx.x;
	int i;

	while(one_start+id<one_end)
	{
		i=id+two_start;
		if(i>=two_end)
		{
			i=two_end-1;
		}
		if(d_primer[4*i]>=d_primer[(id+one_start)*4]+d_primer[(id+one_start)*4+1])
		{
			while((i>=two_start)&&(d_primer[4*i]>=d_primer[(id+one_start)*4]+d_primer[(id+one_start)*4+1]))
			{
				next[id]=i;
				i--;
			}
		}
		else
		{
			while((i<two_end)&&(d_primer[4*i]<d_primer[(id+one_start)*4]+d_primer[(id+one_start)*4+1]))
				i++;
			next[id]=i;
			if(i==two_end)
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

__device__ int check_common_loop(int *d_primer,int *d_common,int *d_sc,int *d_ec,int turn[],int LF,int LB,int *d_int,int *d_apply)
{
        int dis,i,pos[7];
//plus
        for(pos[0]=d_sc[turn[0]];pos[0]<d_ec[turn[0]];pos[0]++)
        {
                if(d_common[4*pos[0]+2]!=1)
                        continue;
		i=d_common[4*pos[0]];
                if(d_apply[turn[0]*d_int[7]+i]==0||d_apply[turn[0]*d_int[7]+i]==2)
                        continue;
                for(pos[1]=d_sc[turn[1]];pos[1]<d_ec[turn[1]];pos[1]++)
                {
                        if(d_common[4*pos[1]]!=i)
                                continue;
                        if(d_common[4*pos[1]+2]!=1)
                                continue;
                        for(pos[2]=d_sc[turn[2]];pos[2]<d_ec[turn[2]];pos[2]++)
                        {
                                if(d_common[4*pos[2]]!=i)
                                        continue;
                                if(d_common[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=d_sc[turn[3]];pos[3]<d_ec[turn[3]];pos[3]++)
                                {
                                        if(d_common[4*pos[3]]!=i)
                                                continue;
                                        if(d_common[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=d_sc[turn[4]];pos[4]<d_ec[turn[4]];pos[4]++)
                                        {
                                                if(d_common[4*pos[4]]!=i)
                                                        continue;
                                                if(d_common[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=d_sc[turn[5]];pos[5]<d_ec[turn[5]];pos[5]++)
                                                {
                                                        if(d_common[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_common[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_common[4*pos[1]+1]-(d_common[4*pos[0]+1]+d_primer[4*turn[0]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=d_common[4*pos[2]+1]-d_common[4*pos[1]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_common[4*pos[3]+1]-(d_common[4*pos[2]+1]+d_primer[4*turn[2]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=(d_common[4*pos[4]+1]+d_primer[4*turn[4]+1]-1)-(d_common[4*pos[3]+1]+d_primer[4*turn[3]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_common[4*pos[4]+1]+d_primer[4*turn[4]+1]-1-d_common[4*pos[1]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_common[4*pos[5]+1]-(d_common[4*pos[4]+1]+d_primer[4*turn[4]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //LF
                                                        if(LF!=-1)
                                                        {
                                                                dis=0;
                                                                for(pos[6]=d_sc[LF];pos[6]<d_ec[LF];pos[6]++)
                                                                {
                                                                        if(d_common[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_common[4*pos[6]+3]!=1)
                                                                                continue;
                                                                        if(d_common[4*pos[1]+1]+d_primer[4*turn[1]+1]>d_common[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_common[4*pos[6]+1]+d_primer[4*LF+1]>d_common[4*pos[2]+1])
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
                                                                for(pos[6]=d_sc[LB];pos[6]=d_ec[LB];pos[6]++)
                                                                {
                                                                        if(d_common[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_common[4*pos[6]+2]!=1)
                                                                                continue;
                                                                        if(d_common[4*pos[3]+1]+d_primer[4*turn[3]+1]>d_common[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_common[4*pos[6]+1]+d_primer[4*LB+1]>d_common[4*pos[4]+1])
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
	for(pos[0]=d_sc[turn[0]];pos[0]<d_ec[turn[0]];pos[0]++)
        {
                if(d_common[4*pos[0]+3]!=1)
                        continue;
                i=d_common[4*pos[0]];
                if(d_apply[turn[0]*d_int[7]+i]==0||d_apply[turn[0]*d_int[7]+i]==2)
                        continue;
                for(pos[1]=d_sc[turn[1]];pos[1]<d_ec[turn[1]];pos[1]++)
                {
                        if(d_common[4*pos[1]]!=i)
                                continue;
                        if(d_common[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=d_sc[turn[2]];pos[2]<d_ec[turn[2]];pos[2]++)
                        {
                                if(d_common[4*pos[2]]!=i)
                                        continue;
                                if(d_common[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=d_sc[turn[3]];pos[3]<d_ec[turn[3]];pos[3]++)
                                {
                                        if(d_common[4*pos[3]]!=i)
                                                continue;
                                        if(d_common[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=d_sc[turn[4]];pos[4]<d_ec[turn[4]];pos[4]++)
                                        {
                                                if(d_common[4*pos[4]]!=i)
                                                        continue;
                                                if(d_common[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=d_sc[turn[5]];pos[5]<d_ec[turn[5]];pos[5]++)
                                                {
                                                        if(d_common[4*pos[5]]!=i)
                                                                continue;
                                                        if(d_common[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=d_common[4*pos[0]+1]-(d_common[4*pos[1]+1]+d_primer[4*turn[1]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=(d_common[4*pos[1]+1]+d_primer[4*turn[1]+1]-1)-(d_common[4*pos[2]+1]+d_primer[4*turn[2]+1]-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=d_common[4*pos[2]+1]-(d_common[4*pos[3]+1]+d_primer[4*turn[3]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=d_common[4*pos[3]+1]-d_common[4*pos[4]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=d_common[4*pos[1]+1]+d_primer[4*turn[1]+1]-1-d_common[4*pos[4]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=d_common[4*pos[4]+1]-(d_common[4*pos[5]+1]+d_primer[4*turn[5]+1]-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //LF
                                                        if(LF!=-1)
                                                        {
                                                                dis=0;
                                                                for(pos[6]=d_sc[LF];pos[6]<d_ec[LF];pos[6]++)
                                                                {
                                                                        if(d_common[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_common[4*pos[6]+2]!=1)
                                                                                continue;
                                                                        if(d_common[4*pos[2]+1]+d_primer[4*turn[2]+1]>d_common[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_common[4*pos[6]+1]+d_primer[4*LF+1]>d_common[4*pos[1]+1])
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
                                                                for(pos[6]=d_sc[LB];pos[6]=d_ec[LB];pos[6]++)
                                                                {
                                                                        if(d_common[4*pos[6]]!=i)
                                                                                continue;
                                                                        if(d_common[4*pos[6]+3]!=1)
                                                                                continue;
                                                                        if(d_common[4*pos[4]+1]+d_primer[4*turn[4]+1]>d_common[4*pos[6]+1])
                                                                                continue;
                                                                        if(d_common[4*pos[6]+1]+d_primer[4*LB+1]>d_common[4*pos[3]+1])
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

__device__ int check_structure_loop(char *d_seq,int *d_primer,int turn[],int *d_int,int LF,int LB,double *parameter,double *d_DPT,int id,int *d_ps)
{
        int i;
        double TH;

        if(LF!=-1)
        {
                for(i=0;i<=1;i++)
                {
                        TH=thal(d_seq,d_primer,turn[i],LF,d_int[11+i],1,1,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,turn[i],LF,d_int[11+i],1,2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,turn[i],LF,d_int[11+i],1,3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;

                        TH=thal(d_seq,d_primer,LF,turn[i],0,(1-d_int[11+i]),2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,LF,turn[i],0,(1-d_int[11+i]),3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                }

                for(i=2;i<6;i++)
                {
                        TH=thal(d_seq,d_primer,LF,turn[i],1,d_int[11+i],1,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,LF,turn[i],1,d_int[11+i],2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,LF,turn[i],1,d_int[11+i],3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;

                        TH=thal(d_seq,d_primer,turn[i],LF,(1-d_int[11+i]),0,2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9]) 
                                return 0;
                        TH=thal(d_seq,d_primer,turn[i],LF,(1-d_int[11+i]),0,3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                }
        }
        if(LB!=-1)
        {
                for(i=0;i<4;i++)
                {
                        TH=thal(d_seq,d_primer,turn[i],LB,d_int[11+i],0,1,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,turn[i],LB,d_int[11+i],0,2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,turn[i],LB,d_int[11+i],0,3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;

                        TH=thal(d_seq,d_primer,LB,turn[i],1,(1-d_int[11+i]),2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,LB,turn[i],1,(1-d_int[11+i]),3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                }
                for(i=4;i<6;i++)
                {
                        TH=thal(d_seq,d_primer,LB,turn[i],0,d_int[11+i],1,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,LB,turn[i],0,d_int[11+i],2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,LB,turn[i],0,d_int[11+i],3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;

                        TH=thal(d_seq,d_primer,turn[i],LB,(1-d_int[11+i]),1,2,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                        TH=thal(d_seq,d_primer,turn[i],LB,(1-d_int[11+i]),1,3,parameter,d_DPT,id,d_ps);
                        if(TH>44+5*d_int[9])
                                return 0;
                }
        }
        if(LF!=-1&&LB!=-1)
        {
                TH=thal(d_seq,d_primer,LF,LB,1,0,1,parameter,d_DPT,id,d_ps);
                if(TH>44+5*d_int[9])
                        return 0;
                TH=thal(d_seq,d_primer,LF,LB,1,0,2,parameter,d_DPT,id,d_ps);
                if(TH>44+5*d_int[9])
                        return 0;
                TH=thal(d_seq,d_primer,LF,LB,1,0,3,parameter,d_DPT,id,d_ps);
                if(TH>44+5*d_int[9])
                        return 0;

                TH=thal(d_seq,d_primer,LB,LF,1,0,2,parameter,d_DPT,id,d_ps);
                if(TH>44+5*d_int[9])
                        return 0;
                TH=thal(d_seq,d_primer,LB,LF,1,0,3,parameter,d_DPT,id,d_ps);
                if(TH>44+5*d_int[9])
                        return 0;
        }
        return 1;
}

__device__ int design_loop(int *d_primer,int *d_SLp,int *d_LLp,int *d_LpLp,char *d_seq,int *d_common,int *d_sc,int *d_ec,int turn[],int *d_int,int *d_apply,int *d_par,double *parameter,double *d_DPT,int id,int *d_ps)
{
        int success,LF,LB;

//LF and LB 
        success=0;
	LF=d_SLp[turn[1]];
        while(LF<d_int[2]+d_int[0]+d_int[1])
        {
		if(LF==-1)
			break;
		if(d_primer[4*LF+3]!=1)
		{
			LF++;
			continue;
		}
                if(d_primer[4*LF]+18>d_primer[4*turn[2]])
                        break;
                LB=d_LLp[turn[3]-d_int[0]];
		if(LB==-1||d_primer[4*LB]+18>d_primer[4*turn[4]])
			break;
                while(LB<d_int[2]+d_int[0]+d_int[1])
                {
			if(d_primer[4*LB+2]!=1)
			{
				LB++;
				continue;
			}
                        if(d_primer[4*LB]+18>d_primer[4*turn[4]])
                                break;
                //check_common
                        if(d_int[3])
                        {
                                success=check_common_loop(d_primer,d_common,d_sc,d_ec,turn,LF,LB,d_int,d_apply);
                                if(success==0)
                                {
                                        LB++;
                                        continue;
                                }
                        }
                //check_structure
                        if(d_int[6])
                        {
                                success=check_structure_loop(d_seq,d_primer,turn,d_int,LF,LB,parameter,d_DPT,id,d_ps);
                                if(success==0)
                                {
                                        LB++;
                                        continue;
                                }
                        }
                        d_par[16*turn[0]+12]=d_primer[4*LF];
			d_par[16*turn[0]+13]=d_primer[4*LF+1];
			d_par[16*turn[0]+14]=d_primer[4*LB];
			d_par[16*turn[0]+15]=d_primer[4*LB+1];
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
        while(LF<d_int[2]+d_int[1]+d_int[0])
        {
		if(LF==-1)
			break;
                if(d_primer[4*LF]+18>d_primer[4*turn[2]])
                        break;
		if(d_primer[4*LF+3]!=1)
		{
			LF++;
			continue;
		}
        //check_common
                if(d_int[3])
                {
                        success=check_common_loop(d_primer,d_common,d_sc,d_ec,turn,LF,-1,d_int,d_apply);
                        if(success==0)
                        {
                                LF++;
                                continue;
                        }
                }
        //check_structure
                if(d_int[6])
                {
                        success=check_structure_loop(d_seq,d_primer,turn,d_int,LF,-1,parameter,d_DPT,id,d_ps);
                        if(success==0)
                        {
                                LF++;
                                continue;
                        }
                }
		d_par[16*turn[0]+12]=d_primer[4*LF];
		d_par[16*turn[0]+13]=d_primer[4*LF+1];
		d_par[16*turn[0]+14]=0;
		d_par[16*turn[0]+15]=0;
                success=1;
                break;
        }
        if(success==1)
                return success;
//only LB
        LB=d_LLp[turn[3]-d_int[0]];
        while(LB<d_int[2]+d_int[0]+d_int[1])
        {
		if(LB==-1)
			break;
                if(d_primer[LB*4]+18>d_primer[4*turn[4]])
                        break;
		if(d_primer[LB*4+2]!=1)
		{
			LB++;
			continue;
		}
        //check_common
                if(d_int[3])
                {
                        success=check_common_loop(d_primer,d_common,d_sc,d_ec,turn,-1,LB,d_int,d_apply);
                        if(success==0)
                        {
                                LB++;
                                continue;
                        }
                }
        //check_structure
                if(d_int[6])
                {
                        success=check_structure_loop(d_seq,d_primer,turn,d_int,-1,LB,parameter,d_DPT,id,d_ps);
                        if(success==0)
                        {
                                LB++;
                                continue;
                        }
                }
		d_par[16*turn[0]+12]=0;
		d_par[16*turn[0]+13]=0;
		d_par[16*turn[0]+14]=d_primer[4*LB];
		d_par[16*turn[0]+15]=d_primer[4*LB+1];
                success=1;
                break;
        }
        return success;
}

//caculate
__global__ void LAMP(char *d_seq,int *d_primer,int *d_common,int *d_special,int *d_sc,int *d_ec,int *d_ss,int *d_es,int *d_SS,int *d_SL,int *d_SLp,int *d_LL,int *d_LS,int *d_LLp,int *d_LpLp,int *d_int,int *d_par,int *d_apply,double *parameter,int *d_pos,double *d_TH,double *d_DPT,int *d_ps)
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
			if(d_primer[4*id]-d_pos[turn[1]]<300&&d_primer[4*id]-d_pos[turn[1]]>-300)
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
		if(d_primer[4*id+2]!=1)
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
			if(d_primer[4*turn[1]+2]!=1)
				continue;
			if(d_primer[4*turn[1]]-(d_primer[4*turn[0]]+d_primer[4*turn[0]+1])>20)
				break;
		//	d_par[16*id]=1;
		//	d_par[16*id+1]=1;
		//	flag=1;
		//	break;
			for(turn[2]=d_SL[turn[1]];turn[2]<d_int[1]+d_int[0];turn[2]++) //F1c
			{
				if(turn[2]==-1)
					break;
				if(flag!=0)
					break;
				if(d_primer[4*turn[2]+3]!=1)
					continue;
				if(d_primer[4*turn[2]]-d_primer[4*turn[1]]-1<40)
					continue;
                                if(d_primer[4*turn[2]]-d_primer[4*turn[1]]-1>60)
                                	break;
		//	d_par[16*id]=2;
                  //      d_par[16*id+1]=1;
                    //    flag=1;
                      //  break;
                                for(turn[3]=d_LL[turn[2]-d_int[0]];turn[3]<d_int[1]+d_int[0];turn[3]++)   //B1c
                                {
                                        if(turn[3]==-1)
                                        	break;
					if(flag!=0)
						break;
					if(d_primer[4*turn[3]+2]!=1)
						continue;
                                        if(d_primer[4*turn[3]]-d_primer[4*turn[2]]>85)
                                        	break;
	//		d_par[16*id]=3;
          //              d_par[16*id+1]=1;
            //            flag=1;
              //          break;
                                        for(turn[4]=d_LS[turn[3]-d_int[0]];turn[4]<d_int[0];turn[4]++)   //B2
                                        {
                                                if(turn[4]==-1)
                                                	break;
						if(flag!=0)
							break;
						if(d_primer[4*turn[4]+3]!=1)
							continue;
                                                if((d_primer[4*turn[4]]+d_primer[4*turn[4]+1]-1)-(d_primer[4*turn[3]]+d_primer[4*turn[3]+1])<40)
                                                	continue;
                                                if((d_primer[4*turn[4]]+d_primer[4*turn[4]+1]-1)-(d_primer[4*turn[3]]+d_primer[4*turn[3]+1])>60)
                                                	break;
                                                if(d_primer[4*turn[4]]+d_primer[4*turn[4]+1]-1-d_primer[turn[1]*4]-1<120)
                                                	continue;
                                                if(d_primer[4*turn[4]]+d_primer[4*turn[4]+1]-1-d_primer[turn[1]*4]-1>180)
                                                	break;
						if(d_int[5]&&(d_SLp[turn[1]]==-1||(d_primer[4*d_SLp[turn[1]]]+18>d_primer[4*turn[2]]))&&(d_LLp[turn[3]-d_int[0]]==-1||(d_primer[4*d_LLp[turn[3]-d_int[0]]]+18>d_primer[4*turn[4]])))
							continue;
			//		d_par[16*id]=4;
			//		d_par[16*id+1]=id;
			//		flag=1;
			//		break;
                                                for(turn[5]=d_SS[turn[4]];turn[5]<d_int[0];turn[5]++)  //B3
                                                {
                                                        if(turn[5]==-1)
                                                        	break;
							if(d_primer[4*turn[5]+3]!=1)
								continue;
                                                        if(d_primer[turn[5]*4]-(d_primer[4*turn[4]]+d_primer[4*turn[4]+1])>20)
                                                        	break;
							flag=check_gc(d_seq,d_primer[4*turn[0]],(d_primer[4*turn[5]]+d_primer[4*turn[5]+1]),d_int[9]);
							if(flag==0)
								continue;

							if(d_int[4]!=0)
							{
								flag=check_uniq(d_primer,d_special,d_ss,d_es,turn);
								if(flag==0)
									continue;
							}

				//	d_par[16*id]=5;
                                  //      d_par[16*id+1]=id;
                                    //    flag=1;
                                      //  break;
							if(d_int[3])
							{
								flag=check_common(d_primer,d_common,d_sc,d_ec,turn,d_int[7],d_apply);
								if(flag<d_int[8])
								{
									flag=0;
									continue;
								}
							}

							if(d_int[6])
							{
								flag=check_structure(d_seq,d_primer,turn,d_int,parameter,d_TH,id,d_DPT,d_ps);
								if(flag==0)
									continue;
							}
	
							if(d_int[5])
							{
								flag=design_loop(d_primer,d_SLp,d_LLp,d_LpLp,d_seq,d_common,d_sc,d_ec,turn,d_int,d_apply,d_par,parameter,d_DPT,id,d_ps);
								if(flag==0)
									continue;
							}
							d_par[id*16]=d_primer[4*turn[0]];
							d_par[id*16+1]=d_primer[4*turn[0]+1];
							d_par[id*16+2]=d_primer[4*turn[1]];
							d_par[id*16+3]=d_primer[4*turn[1]+1];
							d_par[id*16+4]=d_primer[4*turn[2]];
							d_par[id*16+5]=d_primer[4*turn[2]+1];
							d_par[id*16+6]=d_primer[4*turn[3]];
							d_par[id*16+7]=d_primer[4*turn[3]+1];
							d_par[id*16+8]=d_primer[4*turn[4]];
							d_par[id*16+9]=d_primer[4*turn[4]+1];
							d_par[id*16+10]=d_primer[4*turn[5]];
							d_par[id*16+11]=d_primer[4*turn[5]+1];
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
	int i,j,flag[12],expect,circle,have,common_num[1],num[11],max_loop,min_loop,count[3],block,thread;
	char *output,*prefix,*store_path,*path_fa,*inner,*outer,*loop,*par_path,*temp,*seq,*d_seq,primer[26];
	FILE *fp;
	struct Primer *headL,*headS,*headLoop,*tempL,*tempS,*tempLoop,*storeL,*storeS,*storeLoop; 
	struct Node *p_node,*p_temp;
	struct INFO *headList,*p_list;
	time_t start,end;
	double *H_parameter,*parameter;	
	long int memory;
	cudaDeviceProp prop;
	int *d_primer,*d_common,*d_special,*d_sc,*d_ec,*d_ss,*d_es,*d_SS,*d_SL,*d_SLp,*d_LL,*d_LS,*d_LLp,*d_LpLp,*d_apply,*d_par,*d_int,*d_pos;
	int *h_primer,*h_common,*h_special,*h_sc,*h_ec,*h_ss,*h_es,*h_apply,*h_par,h_int[17],*h_pos,*d_ps;
	double *h_TH,*d_TH,*d_DPT;
	
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
                H_parameter=(double *)malloc(5730*sizeof(double));
                memset(H_parameter,'\0',5730*sizeof(double));
                cudaMalloc((void **)&parameter,5730*sizeof(double));
                cudaMemset(parameter,'\0',5730*sizeof(double));

                getStack(par_path,H_parameter);
                getStackint2(par_path,H_parameter);
                getDangle(par_path,H_parameter);
                getLoop(par_path,H_parameter);
                getTstack(par_path,H_parameter);
                getTstack2(par_path,H_parameter);
                tableStartATS(6.9,H_parameter);
                tableStartATH(2200.0,H_parameter);
		cudaMemcpy(parameter,H_parameter,5730*sizeof(double),cudaMemcpyHostToDevice);

		h_int[11]=0; //F3,plus
		h_int[12]=0;
		h_int[13]=1; //F1c,minus
		h_int[14]=0;
		h_int[15]=1;
		h_int[16]=1;
        }
//F3's pos 
	h_pos=(int *)malloc(expect*sizeof(int));
	for(i=0;i<expect;i++)
		h_pos[i]=-1;
	cudaMalloc((void **)&d_pos,expect*sizeof(int));
//d_int
	cudaMalloc((void **)&d_int,17*sizeof(int));
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
	thread=1;
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
			while(tempS&&(memory<prop.totalGlobalMem/6)&&num[2]<20000)
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
				if(flag[7])
					memory=memory+5000+50; //one double=4 int, DPT; 50: ps1+ps2
				tempS=tempS->next;
			}
			if(num[2]<4||num[1]<2||(flag[10]&&num[7]<1)) //don't have enough primers
			{
				num[10]=1;
				break;
			}
			if(tempS==NULL)  //check all primers
				num[10]=1;	

			printf("memory is %ld\n",2*memory);
		//malloc
			h_primer=(int *)malloc(4*(num[2]+num[1]+num[7])*sizeof(int));
			cudaMalloc((void **)&d_primer,4*(num[2]+num[1]+num[7])*sizeof(int));
			if(flag[5])
			{
				h_common=(int *)malloc(4*(num[5]+num[3]+num[8])*sizeof(int));
				cudaMalloc((void **)&d_common,4*(num[5]+num[3]+num[8])*sizeof(int));
				h_sc=(int *)malloc((num[2]+num[1]+num[7])*sizeof(int));
				cudaMalloc((void **)&d_sc,(num[2]+num[1]+num[7])*sizeof(int));
	                        h_ec=(int *)malloc((num[2]+num[1]+num[7])*sizeof(int));
	                        cudaMalloc((void **)&d_ec,(num[2]+num[1]+num[7])*sizeof(int));
			}
			if(flag[6])
			{
				h_special=(int *)malloc(4*(num[6]+num[4])*sizeof(int));
				cudaMalloc((void **)&d_special,4*(num[6]+num[4])*sizeof(int));
				h_ss=(int *)malloc((num[2]+num[1])*sizeof(int));
	                        cudaMalloc((void **)&d_ss,(num[2]+num[1])*sizeof(int));
	                        h_es=(int *)malloc((num[2]+num[1])*sizeof(int));
	                        cudaMalloc((void **)&d_es,(num[2]+num[1])*sizeof(int));
			}
		//small
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
				h_primer[4*count[0]]=tempS->pos;
				h_primer[4*count[0]+1]=tempS->len;
				h_primer[4*count[0]+2]=tempS->plus;
				h_primer[4*count[0]+3]=tempS->minus;
			//common
				if(flag[5])
				{
					h_sc[count[0]]=count[1];
					if(tempS->total_common==0)
						h_ec[count[0]]=-1;
					else
					{
						p_node=tempS->common;
						while(p_node)
						{
							h_common[4*count[1]]=p_node->gi;
							h_common[4*count[1]+1]=p_node->pos;
							h_common[4*count[1]+2]=p_node->plus; 
                                		        h_common[4*count[1]+3]=p_node->minus;
							count[1]++;
							p_node=p_node->next;
						}
						h_ec[count[0]]=count[1];
					}
				}
			//special
				if(flag[6])
				{
					h_ss[count[0]]=count[2];
                	        	if(tempS->total_special==0)
                	        	        h_es[count[0]]=-1;
                	        	else
                	        	{
                	        	        p_node=tempS->special;
                	        	        while(p_node)
                	        	        {
                	        	                h_special[4*count[2]]=p_node->gi;
                	        	                h_special[4*count[2]+1]=p_node->pos;
                	        	                h_special[4*count[2]+2]=p_node->plus;
                	        	                h_special[4*count[2]+3]=p_node->minus;
                	        	                count[2]++;
                	        	                p_node=p_node->next;
                	        	        }
                	        	        h_es[count[0]]=count[2];
                	        	}
				}
				count[0]++;
				tempS=tempS->next;
			}
			h_int[0]=num[2];

		//large primer
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
                	        h_primer[4*num[2]+4*count[0]]=tempL->pos;
                	        h_primer[4*num[2]+4*count[0]+1]=tempL->len;
                	        h_primer[4*num[2]+4*count[0]+2]=tempL->plus;
                	        h_primer[4*num[2]+4*count[0]+3]=tempL->minus;
                	//common
				if(flag[5])
				{
                	        	h_sc[num[2]+count[0]]=count[1];
                	        	if(tempL->total_common==0)
                	        	        h_ec[num[2]+count[0]]=-1;
                	        	else
                	        	{
                	        	        p_node=tempL->common;
                	        	        while(p_node)
                	        	        {
                	        	                h_common[4*num[5]+4*count[1]]=p_node->gi;
                	        	                h_common[4*num[5]+4*count[1]+1]=p_node->pos;
                	        	                h_common[4*num[5]+4*count[1]+2]=p_node->plus; 
                        		                h_common[4*num[5]+4*count[1]+3]=p_node->minus;
                        		                count[1]++;
                        		                p_node=p_node->next;
                        		        }
                        		        h_ec[num[2]+count[0]]=count[1];
					}
                        	}
                	//special
				if(flag[6])
				{
                        		h_ss[num[2]+count[0]]=count[2];
                        		if(tempL->total_special==0)
                        		        h_es[num[2]+count[0]]=-1;
                        		else
                        		{
                        		        p_node=tempL->special;
                        		        while(p_node)
                        		        {
                        		                h_special[4*num[6]+4*count[2]]=p_node->gi;
                        		                h_special[4*num[6]+4*count[2]+1]=p_node->pos;
                        		                h_special[4*num[6]+4*count[2]+2]=p_node->plus;
                        		                h_special[4*num[6]+4*count[2]+3]=p_node->minus;
                        		                count[2]++;
                        		                p_node=p_node->next;
                        		        }
                        		        h_es[num[2]+count[0]]=count[2];
                        		}
				}
                        	count[0]++;
                        	tempL=tempL->next;
                	}
			h_int[1]=num[1];

                //loop primer
			if(flag[10])
			{
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
                                	h_primer[4*(num[1]+num[2])+4*count[0]]=tempLoop->pos;
                                	h_primer[4*(num[1]+num[2])+4*count[0]+1]=tempLoop->len;
                                	h_primer[4*(num[1]+num[2])+4*count[0]+2]=tempLoop->plus;
                                	h_primer[4*(num[1]+num[2])+4*count[0]+3]=tempLoop->minus;
                        	//common
                                	if(flag[5])
                                	{
                                	        h_sc[num[1]+num[2]+count[0]]=count[1];
                                	        if(tempLoop->total_common==0)
                                	                h_ec[num[1]+num[2]+count[0]]=-1;
                                	        else
                                	        {
                                	                p_node=tempLoop->common;
                                	                while(p_node)
                                	                {
                                	                        h_common[4*(num[5]+num[3])+4*count[1]]=p_node->gi;
                                	                        h_common[4*(num[5]+num[3])+4*count[1]+1]=p_node->pos;
                                	                        h_common[4*(num[5]+num[3])+4*count[1]+2]=p_node->plus; 
                                	                        h_common[4*(num[5]+num[3])+4*count[1]+3]=p_node->minus;
                                	                        count[1]++;
                                	                        p_node=p_node->next;
                                	                }
                                	                h_ec[num[1]+num[2]+count[0]]=count[1];
                                	        }
                                	}
	                                count[0]++;
	                                tempLoop=tempLoop->next;
	                        }
				h_int[2]=num[7];
                        }
	//run
			if(num[2]%thread==0)
                                block=num[2]/thread;
                        else
                                block=(num[2]-num[2]%thread)/thread+1;

			if(block>prop.maxGridSize[0]/2)
				block=prop.maxGridSize[0]/2;

			cudaMemcpy(d_primer,h_primer,4*(num[1]+num[2]+num[7])*sizeof(int),cudaMemcpyHostToDevice);
			free(h_primer);
			if(flag[5])
			{
				cudaMemcpy(d_common,h_common,4*(num[3]+num[5]+num[8])*sizeof(int),cudaMemcpyHostToDevice);
				free(h_common);
				cudaMemcpy(d_sc,h_sc,(num[1]+num[2]+num[7])*sizeof(int),cudaMemcpyHostToDevice);
				free(h_sc);
				cudaMemcpy(d_ec,h_ec,(num[1]+num[2]+num[7])*sizeof(int),cudaMemcpyHostToDevice);
				free(h_ec);
			}
			if(flag[6])
			{
				cudaMemcpy(d_special,h_special,4*(num[4]+num[6])*sizeof(int),cudaMemcpyHostToDevice);
				free(h_special);
				cudaMemcpy(d_ss,h_ss,(num[1]+num[2])*sizeof(int),cudaMemcpyHostToDevice);
				free(h_ss);
				cudaMemcpy(d_es,h_es,(num[1]+num[2])*sizeof(int),cudaMemcpyHostToDevice);
				free(h_es);
			}
		//next primer
		printf("block is %d,thread is %d\n",block,thread);
			next_one<<<block,thread>>>(d_primer,d_SS,0,num[2],0,num[2]);
			next_one<<<block,thread>>>(d_primer,d_LL,num[2],(num[1]+num[2]),num[2],(num[1]+num[2]));
			next_one<<<block,thread>>>(d_primer,d_SL,0,num[2],num[2],(num[1]+num[2]));
        	        next_one<<<block,thread>>>(d_primer,d_LS,num[2],(num[1]+num[2]),0,num[2]);
			if(flag[10])
			{
				next_one<<<block,thread>>>(d_primer,d_LpLp,(num[1]+num[2]),(num[1]+num[2]+num[7]),(num[1]+num[2]),(num[1]+num[2]+num[7]));
                        	next_one<<<block,thread>>>(d_primer,d_SLp,0,num[2],(num[1]+num[2]),(num[1]+num[2]+num[7]));
                        	next_one<<<block,thread>>>(d_primer,d_LLp,num[2],(num[1]+num[2]),(num[1]+num[2]),(num[1]+num[2]+num[7]));
			}
		//calculate
			if(flag[5])
				cudaMalloc((void **)&d_apply,common_num[0]*num[2]*sizeof(int));
			cudaMemcpy(d_pos,h_pos,expect*sizeof(int),cudaMemcpyHostToDevice);
			h_int[8]=circle;
			cudaMemcpy(d_int,h_int,17*sizeof(int),cudaMemcpyHostToDevice);
			cudaMemcpy(d_pos,h_pos,expect*sizeof(int),cudaMemcpyHostToDevice);
			cudaMalloc((void **)&d_par,num[2]*16*sizeof(int));
		cudaMalloc((void **)&d_TH,2*num[2]*sizeof(double));
		h_TH=(double *)malloc(2*num[2]*sizeof(double));
			cudaMalloc((void **)&d_DPT,num[2]*1250*sizeof(double));
			cudaMalloc((void **)&d_ps,num[2]*50*sizeof(int));
			LAMP<<<block,thread>>>(d_seq,d_primer,d_common,d_special,d_sc,d_ec,d_ss,d_es,d_SS,d_SL,d_SLp,d_LL,d_LS,d_LLp,d_LpLp,d_int,d_par,d_apply,parameter,d_pos,d_TH,d_DPT,d_ps);
			cudaFree(d_DPT);
			cudaFree(d_ps);
		cudaMemcpy(h_TH,d_TH,2*num[2]*sizeof(double),cudaMemcpyDeviceToHost);
		cudaFree(d_TH);
		printf("%lf\t%lf\n",h_TH[0],h_TH[1]);
		free(h_TH);
			h_par=(int *)malloc(16*num[2]*sizeof(int));
			memset(h_par,'\0',16*num[2]*sizeof(int));
			cudaMemcpy(h_par,d_par,16*num[2]*sizeof(int),cudaMemcpyDeviceToHost);
			cudaFree(d_par);
			if(flag[5])
			{
				h_apply=(int *)malloc(common_num[0]*num[2]*sizeof(int));
				cudaMemcpy(h_apply,d_apply,common_num[0]*num[2]*sizeof(int),cudaMemcpyDeviceToHost);
				cudaFree(d_apply);
			}
	//free
			cudaFree(d_primer);
			cudaFree(d_SS);
			cudaFree(d_LL);
			cudaFree(d_SL);
			cudaFree(d_LS);
			if(flag[10])
			{
				cudaFree(d_LpLp);
				cudaFree(d_SLp);
				cudaFree(d_LLp);
			}
			if(flag[5])
			{
				cudaFree(d_common);
				cudaFree(d_sc);
				cudaFree(d_ec);
			}
			if(flag[6])
			{
                                cudaFree(d_special);
                                cudaFree(d_ss);
                                cudaFree(d_es);
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
			if(flag[5])
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
                free(H_parameter);
                cudaFree(parameter);
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
