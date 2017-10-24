#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

const double  _INFINITY=INFINITY;
#define isFinite(x) finite(x)

char str2int(char c)
{
        switch (c)
        {
                case 'A': case '0':
                        return 0;
                case 'C': case '1':
                        return 1;
                case 'G': case '2':
                        return 2;
                case 'T': case '3':
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
        while (isspace(*p)) 
                p++;
        while (isdigit(*p)) 
                p++;
        while (isspace(*p)) 
                p++;

        q = p;
        while (!isspace(*q)) 
                q++;
        *q = '\0';
        q++;
        if (!strcmp(p, "inf"))
                *v1 =1.0*INFINITY;
        else 
                sscanf(p, "%lf", v1);
        while (isspace(*q))
                q++;

        p = q;
        while (!isspace(*p))
                p++;
        *p = '\0';
        p++;
        if (!strcmp(q, "inf"))
                *v2 =1.0*INFINITY;
        else 
                sscanf(q, "%lf", v2);
        while (isspace(*p))
                p++;

        q = p;
        while (!isspace(*q) && (*q != '\0'))
                q++;
        *q = '\0';
        if (!strcmp(p, "inf"))
                *v3 =1.0*INFINITY;
        else 
                sscanf(p, "%lf", v3);
}

int readTLoop(FILE *file, char *s,double *v, int triloop)
{
        char *line,*p,*q;

        line=(char *)malloc(200);
        memset(line,'\0',200);
        fgets(line,200,file);
        if (!isalpha(line[0]))
                return -1;

        p = line;
/* skip first spaces */
        while (isspace(*p))
                p++;
/* read the string */
        q = p;
        while (isalpha(*q))
                q++;
        *q = '\0';
        q++;
        if (triloop)
        {
                strncpy(s, p, 5);   /*triloop string has 5 characters*/
                s[5] = '\0';
        }
        else
        {
                strncpy(s, p, 6);   /*tetraloop string has 6 characters*/
                s[6] = '\0';
        }

/* skip all spaces */
        while (isspace(*q))
                q++;
        p = q;
        while (!isspace(*p) && (*p != '\0'))
                p++;
        *p = '\0';
        if (!strcmp(q, "inf"))
                *v =1.0*INFINITY;
        else 
                sscanf(q, "%lf", v);
        return 0;
}

void getStack(double stackEntropies[],double stackEnthalpies[], char *path)
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
                                                stackEntropies[i*125+ii*25+j*5+jj] = -1.0;
                                                stackEnthalpies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                        }
                                        else 
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackEntropies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackEntropies[i*125+ii*25+j*5+jj] = atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackEnthalpies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackEnthalpies[i*125+ii*25+j*5+jj] = atof(line);

                                                if (!isFinite(stackEntropies[i*125+ii*25+j*5+jj]) || !isFinite(stackEnthalpies[i*125+ii*25+j*5+jj])) 
                                                {
                                                        stackEntropies[i*125+ii*25+j*5+jj] = -1.0;
                                                        stackEnthalpies[i*125+ii*25+j*5+jj] =1.0*INFINITY;
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

void getStackint2(double stackint2Entropies[],double stackint2Enthalpies[], char *path)
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
                                                stackint2Entropies[i*125+ii*25+j*5+jj] = -1.0;
                                                stackint2Enthalpies[i*125+ii*25+j*5+jj] =1.0*INFINITY;
                                        } 
                                        else 
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStackint2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackint2Entropies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackint2Entropies[i*125+ii*25+j*5+jj] = atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStackint2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackint2Enthalpies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackint2Enthalpies[i*125+ii*25+j*5+jj] = atof(line);

                                                if(!isFinite(stackint2Entropies[i*125+ii*25+j*5+jj]) || !isFinite(stackint2Enthalpies[i*125+ii*25+j*5+jj]))
                                                {
                                                        stackint2Entropies[i*125+ii*25+j*5+jj] = -1.0;
                                                        stackint2Enthalpies[i*125+ii*25+j*5+jj] =1.0*INFINITY;
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

void getDangle(double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],char *path)
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
                                        dangleEntropies3[i*25+k*5+j] = -1.0;
                                        dangleEnthalpies3[i*25+k*5+j] =1.0*INFINITY;
                                }
                                else if (k == 4)
                                {
                                        dangleEntropies3[i*25+k*5+j] = -1.0;
                                        dangleEnthalpies3[i*25+k*5+j] =1.0*INFINITY;
                                } 
                                else
                                {
                                        if(fgets(line,20,sFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");
                                                exit(1);
                                        }
                                        if(strncmp(line, "inf", 3)==0)
                                                dangleEntropies3[i*25+k*5+j]=1.0*INFINITY;
                                        else
                                                dangleEntropies3[i*25+k*5+j]=atof(line);

                                        if(fgets(line,20,hFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");        
                                                exit(1);        
                                        }
                                        if(strncmp(line, "inf", 3)==0)        
                                                dangleEnthalpies3[i*25+k*5+j]=1.0*INFINITY;           
                                        else        
                                                dangleEnthalpies3[i*25+k*5+j]=atof(line);

                                        if(!isFinite(dangleEntropies3[i*25+k*5+j]) || !isFinite(dangleEnthalpies3[i*25+k*5+j])) 
                                        {
                                                dangleEntropies3[i*25+k*5+j] = -1.0;
                                                dangleEnthalpies3[i*25+k*5+j] =1.0*INFINITY;
                                        }
                                }
                        }

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        for (k = 0; k < 5; ++k) 
                        {
                                if (i == 4 || j == 4)
                                {
                                        dangleEntropies5[i*25+j*5+k] = -1.0;
                                        dangleEnthalpies5[i*25+j*5+k] =1.0*INFINITY;
                                } 
                                else if (k == 4) 
                                {
                                        dangleEntropies5[i*25+j*5+k] = -1.0;
                                        dangleEnthalpies5[i*25+j*5+k] =1.0*INFINITY;
                                }
                                else
                                {
                                        if(fgets(line,20,sFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");
                                                exit(1);
                                        }
                                        if(strncmp(line, "inf", 3)==0)
                                                dangleEntropies5[i*25+j*5+k]=1.0*INFINITY;
                                        else
                                                dangleEntropies5[i*25+j*5+k]=atof(line);

                                        if(fgets(line,20,hFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");        
                                                exit(1);        
                                        }
                                        if(strncmp(line, "inf", 3)==0)        
                                                dangleEnthalpies5[i*25+j*5+k]=1.0*INFINITY;           
                                        else        
                                                dangleEnthalpies5[i*25+j*5+k]=atof(line);

                                        if(!isFinite(dangleEntropies5[i*25+j*5+k]) || !isFinite(dangleEnthalpies5[i*25+j*5+k]))
                                        {
                                                dangleEntropies5[i*25+j*5+k] = -1.0;
                                                dangleEnthalpies5[i*25+j*5+k] =1.0*INFINITY;
                                        }
                                }
                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getLoop(double hairpinLoopEntropies[30],double interiorLoopEntropies[30],double bulgeLoopEntropies[30],double hairpinLoopEnthalpies[30],double interiorLoopEnthalpies[30],double bulgeLoopEnthalpies[30],char *path)
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
                readLoop(sFile, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k]);
                readLoop(hFile, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k]);
        }
        fclose(sFile);
        fclose(hFile);
}

void getTstack(double tstackEntropies[],double tstackEnthalpies[],char *path)
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
                                                tstackEnthalpies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                tstackEntropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                        }
                                        else if (i2 == 4 || j2 == 4)
                                        {
                                                tstackEntropies[i1*125+i2*25+j1*5+j2] = 0.00000000001;
                                                tstackEnthalpies[i1*125+i2*25+j1*5+j2] = 0.0;
                                        }
                                        else
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstackEntropies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstackEntropies[i1*125+i2*25+j1*5+j2]=atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstackEnthalpies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstackEnthalpies[i1*125+i2*25+j1*5+j2]=atof(line);

                                                if (!isFinite(tstackEntropies[i1*125+i2*25+j1*5+j2]) || !isFinite(tstackEnthalpies[i1*125+i2*25+j1*5+j2]))
                                                {
                                                        tstackEntropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                                        tstackEnthalpies[i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                }
                                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getTstack2(double tstack2Entropies[],double tstack2Enthalpies[],char *path)
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
                                                tstack2Enthalpies[i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                tstack2Entropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                        }
                                        else if (i2 == 4 || j2 == 4)
                                        {
                                                tstack2Entropies[i1*125+i2*25+j1*5+j2] = 0.00000000001;
                                                tstack2Enthalpies[i1*125+i2*25+j1*5+j2] = 0.0;
                                        }
                                        else
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstack2Entropies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstack2Entropies[i1*125+i2*25+j1*5+j2]=atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstack2Enthalpies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstack2Enthalpies[i1*125+i2*25+j1*5+j2]=atof(line);


                                                if (!isFinite(tstack2Entropies[i1*125+i2*25+j1*5+j2]) || !isFinite(tstack2Enthalpies[i1*125+i2*25+j1*5+j2]))
                                                {
                                                        tstack2Entropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                                        tstack2Enthalpies[i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
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

void getTriloop(char *triloopEntropies1,char *triloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,char *path)
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
			triloopEntropies1[5*turn+i]=str2int(seq[i]);
		if(value[0]=='i')
			triloopEntropies2[turn]=1.0*INFINITY;
		else
			triloopEntropies2[turn]=atof(value);
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
			triloopEnthalpies1[turn*5+i]=str2int(seq[i]);
		if(value[0]=='i')
			triloopEnthalpies2[turn]=1.0*INFINITY;
		else
			triloopEnthalpies2[turn]=atof(value);
		turn++;
        }
        fclose(hFile);
}

void getTetraloop(char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *tetraloopEntropies2,double *tetraloopEnthalpies2,char *path)
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
			tetraloopEntropies1[turn*6+i]=str2int(seq[i]);
		if(value[0]=='i')
			tetraloopEntropies2[turn]=1.0*INFINITY;
		else
			tetraloopEntropies2[turn]=atof(value);
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
			tetraloopEnthalpies1[6*turn+i]=str2int(seq[i]);
		if(value[0]=='i')
			tetraloopEnthalpies2[turn]=1.0*INFINITY;
		else
			tetraloopEnthalpies2[turn]=atof(value);
		turn++;
        }
        fclose(hFile);
}

void tableStartATS(double atp_value, double atpS[])
{
        int i, j;

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        atpS[i*5+j] = 0.00000000001;
        atpS[3] = atpS[15] = atp_value;
}

void tableStartATH(double atp_value,double atpH[])
{
        int i, j;

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        atpH[i*5+j] = 0.0;
        atpH[3] = atpH[15] = atp_value;
}

void initMatrix2(int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1)
{
	int i,j;
	for(i=1;i<=Initint[0];++i)
		for(j=i;j<=Initint[1];++j)
			if(j-i<4 || (numSeq1[i]+numSeq1[j]!=3))
			{
				enthalpyDPT[(i-1)*Initint[2]+j-1]=_INFINITY;
				entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
			}
			else
			{
				enthalpyDPT[(i-1)*Initint[2]+j-1]=0.0;
				entropyDPT[(i-1)*Initint[2]+j-1]=-3224.0;
			}
}

double Ss(int i,int j,int k,double stackEntropies[],int Initint[],char *numSeq1,char *numSeq2)
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
		return stackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]];
	}
	else
		return stackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
}

double Hs(int i,int j,int k,double stackEnthalpies[],int Initint[],char *numSeq1,char *numSeq2)
{
	if(k==2)
	{
		if(i>= j)
			return _INFINITY;
		if(i==Initint[0]||j==Initint[1]+1)
			return _INFINITY;

		if(i>Initint[0])
			i-=Initint[0];
		if(j>Initint[1])
			j-=Initint[1];
		if(isFinite(stackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]))
			return stackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]];
		else
			return _INFINITY;
	}
	else
		return stackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
}

void maxTM2(int i,int j,double stackEntropies[],double stackEnthalpies[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	double T0,T1,S0,S1,H0,H1;

	S0=entropyDPT[(i-1)*Initint[2]+j-1];
	H0=enthalpyDPT[(i-1)*Initint[2]+j-1];
	T0=(H0+Initdouble[0])/(S0+Initdouble[1]+Initdouble[2]);
	if(isFinite(enthalpyDPT[(i-1)*Initint[2]+j-1]))
	{
		S1=(entropyDPT[i*Initint[2]+j-2]+Ss(i,j,2,stackEntropies,Initint,numSeq1,numSeq2));
		H1=(enthalpyDPT[i*Initint[2]+j-2]+Hs(i,j,2,stackEnthalpies,Initint,numSeq1,numSeq2));
	}
	else
	{
		S1=-1.0;
		H1=_INFINITY;
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

void calc_bulge_internal2(int i,int j,int ii,int jj,double *EntropyEnthalpy,int traceback,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double atpS[],double atpH[],double Initdouble[0],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	int loopSize1,loopSize2,loopSize;
	double T1,T2,S,H;

	T1 = T2 = -_INFINITY;
	S=-3224.0;
	H=0.0;
	loopSize1=ii-i-1;
	loopSize2=j-jj-1;
	if(loopSize1+loopSize2>30)
	{
		EntropyEnthalpy[0]=-1.0;
		EntropyEnthalpy[1]=_INFINITY;
		return;
	}

	loopSize=loopSize1+loopSize2-1;
	if((loopSize1==0&&loopSize2>0)||(loopSize2==0&&loopSize1>0))
	{
		if(loopSize2==1||loopSize1==1)
		{ 
			if((loopSize2==1&&loopSize1==0)||(loopSize2==0&&loopSize1==1))
			{
				H=bulgeLoopEnthalpies[loopSize]+stackEnthalpies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
				S=bulgeLoopEntropies[loopSize]+stackEntropies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
 			}
			if(traceback!=1)
			{
				H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];
				S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
			}

			if(!isFinite(H))
			{
				H=_INFINITY;
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
			H=bulgeLoopEnthalpies[loopSize]+atpH[numSeq1[i]*5+numSeq2[j]]+atpH[numSeq1[ii]*5+numSeq2[jj]];
			if(traceback!=1)
				H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];

			S=bulgeLoopEntropies[loopSize]+atpS[numSeq1[i]*5+numSeq2[j]]+atpS[numSeq1[ii]*5+numSeq2[jj]];
			if(traceback!=1)
				S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
			if(!isFinite(H))
			{
				H=_INFINITY;
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
		S=stackint2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+stackint2Entropies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		if(traceback!=1)
			S+=entropyDPT[(ii-1)*Initint[2]+jj-1];

		H=stackint2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+stackint2Enthalpies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		if(traceback!=1)
			H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];
		if(!isFinite(H))
		{
			H=_INFINITY;
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
		H=interiorLoopEnthalpies[loopSize]+tstackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+tstackEnthalpies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		if(traceback!=1)
			H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];

		S=interiorLoopEntropies[loopSize]+tstackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+tstackEntropies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
		if(traceback!=1)
			S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
		if(!isFinite(H))
		{
			H=_INFINITY;
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

void CBI(int i,int j,double* EntropyEnthalpy,int traceback,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	int d,ii,jj;

	for(d=j-i-3;d>=4&&d>=j-i-32;--d)
		for(ii=i+1;ii<j-d&&ii<=Initint[0];++ii)
		{
			jj=d+ii;
			if(traceback==0)
			{
				EntropyEnthalpy[0]=-1.0;
				EntropyEnthalpy[1]=_INFINITY;
			}
			if(isFinite(enthalpyDPT[(ii-1)*Initint[2]+jj-1]))
			{
				calc_bulge_internal2(i,j,ii,jj,EntropyEnthalpy,traceback,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
				if(isFinite(EntropyEnthalpy[1]))
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

int find_pos(char *ref,int ref_start,char *source,int length,int num)
{
	int flag,i,j;

	for(i=0;i<num;i++)
	{
		flag=0;
		for(j=0;j<length;j++)
		{
			if(ref[ref_start+j]!=source[i*length+j])
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

void calc_hairpin(int i,int j,double *EntropyEnthalpy,int traceback,double hairpinLoopEntropies[],double hairpinLoopEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1)
{
	int pos,loopSize=j-i-1;
	double T1,T2;
	
	T1=T2=-_INFINITY;
	if(loopSize < 3)
	{
		EntropyEnthalpy[0]=-1.0;
		EntropyEnthalpy[1]=_INFINITY;
		return;
	}
	if(i<=Initint[0]&&Initint[1]<j)
	{
		EntropyEnthalpy[0]=-1.0;
		EntropyEnthalpy[1]=_INFINITY;
		return;
	}
	else if(i>Initint[1])
	{
		i-= Initint[0];
		j-= Initint[1];
	}
	if(loopSize<=30)
	{
		EntropyEnthalpy[1]=hairpinLoopEnthalpies[loopSize-1];
		EntropyEnthalpy[0]=hairpinLoopEntropies[loopSize-1];
	}
	else
	{
		EntropyEnthalpy[1]=hairpinLoopEnthalpies[29];
		EntropyEnthalpy[0]=hairpinLoopEntropies[29];
	}

	if(loopSize>3) // for loops 4 bp and more in length, terminal mm are accounted
	{
		EntropyEnthalpy[1]+=tstack2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
		EntropyEnthalpy[0]+=tstack2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
	}
	else if(loopSize == 3) // for loops 3 bp in length at-penalty is considered
	{
		EntropyEnthalpy[1]+=atpH[numSeq1[i]*5+numSeq1[j]];
		EntropyEnthalpy[0]+=atpS[numSeq1[i]*5+numSeq1[j]];
	}

	if(loopSize==3) // closing AT-penalty (+), triloop bonus, hairpin of 3 (+) 
	{
		pos=find_pos(numSeq1,i,triloopEnthalpies1,5,numTriloops);
		if(pos!=-1)
			EntropyEnthalpy[1]+=triloopEnthalpies2[pos];

		pos=find_pos(numSeq1,i,triloopEntropies1,5,numTriloops);
		if(pos!=-1)
			EntropyEnthalpy[0]+=triloopEntropies2[pos];
	}
	else if (loopSize == 4) // terminal mismatch, tetraloop bonus, hairpin of 4
	{
		pos=find_pos(numSeq1,i,tetraloopEnthalpies1,6,numTetraloops);
		if(pos!=-1)
			EntropyEnthalpy[1]+=tetraloopEnthalpies2[pos];

		pos=find_pos(numSeq1,i,tetraloopEntropies1,6,numTetraloops);
		if(pos!=-1)
			EntropyEnthalpy[0]+=tetraloopEntropies2[pos];
	}
	if(!isFinite(EntropyEnthalpy[1]))
	{
		EntropyEnthalpy[1] = _INFINITY;
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

void fillMatrix2(double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	int i, j;
	double SH[2];

	for (j = 2; j <= Initint[1]; ++j)
		for (i = j - 3 - 1; i >= 1; --i)
		{
			if (isFinite(enthalpyDPT[(i-1)*Initint[2]+j-1]))
			{
				SH[0] = -1.0;
				SH[1] = _INFINITY;
				maxTM2(i,j,stackEntropies,stackEnthalpies,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
				CBI(i,j,SH,0,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

				SH[0] = -1.0;
				SH[1] = _INFINITY;
				calc_hairpin(i, j, SH, 0,hairpinLoopEntropies,hairpinLoopEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1);
				if(isFinite(SH[1]))
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

int max5(double a,double b,double c,double d,double e)
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

double Sd5(int i,int j,double dangleEntropies5[],char *numSeq1)
{
	return dangleEntropies5[numSeq1[i]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double Hd5(int i,int j,double dangleEnthalpies5[],char *numSeq1)
{
	return dangleEnthalpies5[numSeq1[i]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double Sd3(int i,int j,double dangleEntropies3[],char *numSeq1)
{
	return dangleEntropies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq1[j]];
}

double Hd3(int i,int j,double dangleEnthalpies3[],char *numSeq1)
{
	return dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq1[j]];
}

double Ststack(int i,int j,double tstack2Entropies[],char *numSeq1)
{
	return tstack2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double Htstack(int i,int j,double tstack2Enthalpies[],char *numSeq1)
{
	return tstack2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double END5_1(int i,int hs,double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,double *send5,double *hend5,char *numSeq1)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	max_tm=-_INFINITY;
	H_max=_INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-5;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1];
			S=send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1];
			if(!isFinite(H)||H>0||S>0)  // H and S must be greater than 0 to avoid BS
			{
				H=_INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1];
			S=atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1];
			if(!isFinite(H)||H>0||S>0)
			{
				H=_INFINITY;
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

double END5_2(int i,int hs,double dangleEntropies5[],double dangleEnthalpies5[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,double *send5,double *hend5,char *numSeq1)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=_INFINITY;
	max_tm=-_INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-6;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1];
			S=send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1];
			if(!isFinite(H)||H>0||S>0)
			{
				H=_INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1];
			S=atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1];
			if(!isFinite(H)||H>0||S>0)
			{
				H=_INFINITY;
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

double END5_3(int i,int hs,double dangleEntropies3[],double dangleEnthalpies3[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,double *send5,double *hend5,char *numSeq1)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=_INFINITY;
	max_tm=-_INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-6;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2];
			S=send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2];
			if(!isFinite(H)||H>0||S>0)
			{
				H=_INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2];
			S=atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2];
			if(!isFinite(H)||H>0||S>0)
			{
				H=_INFINITY;
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

double END5_4(int i,int hs,double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,double *send5,double *hend5,char *numSeq1)
{
	int k;
	double max_tm,T1,T2,H,S,H_max,S_max;

	H_max=_INFINITY;
	max_tm=-_INFINITY;
	S_max=-1.0;
	for(k=0;k<=i-7;++k)
	{
		T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
		T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
		if(T1>=T2)
		{
			H=hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2];
			S=send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2];
			if(!isFinite(H)||H>0||S>0)
			{
				H=_INFINITY;
				S=-1.0;
			}
			T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
		}
		else
		{
			H=atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2];
			S=atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2];
			if(!isFinite(H)||H>0||S>0)
			{
				H=_INFINITY;
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

void calc_terminal_bp(double temp,double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,double *send5,double *hend5,char *numSeq1)
{
	int i,max;
	double T1,T2,T3,T4,T5,G,end5_11,end5_12,end5_21,end5_22,end5_31,end5_32,end5_41,end5_42;
	
	send5[0]=send5[1]= -1.0;
	hend5[0]=hend5[1]=_INFINITY;

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
		end5_11=END5_1(i,1,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		end5_12=END5_1(i,2,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		T2=(end5_11+Initdouble[0])/(end5_12+Initdouble[1]+Initdouble[2]);
		end5_21=END5_2(i,1,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		end5_22=END5_2(i,2,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		T3=(end5_21+Initdouble[0])/(end5_22+Initdouble[1]+Initdouble[2]);
		end5_31=END5_3(i,1,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		end5_32=END5_3(i,2,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		T4=(end5_31+Initdouble[0])/(end5_32+Initdouble[1]+Initdouble[2]);
		end5_41=END5_4(i,1,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		end5_42=END5_4(i,2,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
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

int newpush(int store[],int i,int j,int mtrx,int total,int next)
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

int equal(double a,double b)
{
	if(!finite(a)||!finite(b))
		return 0;
	return fabs(a-b)<1e-5;
}

void tracebacku(int* bp,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,double *send5,double *hend5,char *numSeq1,char *numSeq2)
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
                        if(equal(send5[i],END5_1(i,2,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_1(i,1,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for(k=0;k<=i-5;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1])&&equal(hend5[i],atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);                    
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_2(i,2,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_2(i,1,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1])&&equal(hend5[i],atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_3(i,2,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_3(i,1,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2])&&equal(hend5[i],atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_4(i,2,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_4(i,1,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for (k=0;k<=i-7;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2])&&equal(hend5[i],atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+2,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2]))
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
                        SH1[1]=_INFINITY;
                        calc_hairpin(i,j,SH1,1,hairpinLoopEntropies,hairpinLoopEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1);

                        SH2[0]=-1.0;
                        SH2[1]=_INFINITY;
                        CBI(i,j,SH2,2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

                        if (equal(entropyDPT[(i-1)*Initint[2]+j-1],Ss(i,j,2,stackEntropies,Initint,numSeq1,numSeq2)+entropyDPT[i*Initint[2]+j-2])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],Hs(i,j,2,stackEnthalpies,Initint,numSeq1,numSeq2)+enthalpyDPT[i*Initint[2]+j-2]))
                                total=newpush(store,i+1,j-1,0,total,now+1);
                        else if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH1[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH1[1]));
                        else if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH2[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH2[1]))
                        {
                                for (done=0,d=j-i-3;d>=4&&d>=j-i-32&&!done;--d)
                                        for (ii=i+1;ii<j-d;++ii)
                                        {
                                                jj=d+ii;
                                                EntropyEnthalpy[0]=-1.0;
                                                EntropyEnthalpy[1]=_INFINITY;
                                                calc_bulge_internal2(i,j,ii,jj,EntropyEnthalpy,1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

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

double drawHairpin(int *bp,double mh,double ms,int Initint[])
{
        int i,N;
	double mg,t;

        N=0;
        if(!isFinite(ms)||!isFinite(mh))
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

void initMatrix(int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	int i,j;

	for(i=1;i<=Initint[0];++i)
	{
		for(j=1;j<=Initint[1];++j)
		{
			if(numSeq1[i]+numSeq2[j]!=3)
			{
				enthalpyDPT[(i-1)*Initint[2]+j-1]=_INFINITY;
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

void LSH(int i,int j,double *EntropyEnthalpy,double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	double S1,H1,T1,S2,H2,T2;

	if(numSeq1[i]+numSeq2[j]!=3)
	{
		entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
		enthalpyDPT[(i-1)*Initint[2]+j-1]=_INFINITY;
		return;
	}

	S1=atpS[numSeq1[i]*5+numSeq2[j]]+tstack2Entropies[numSeq2[j]*125+numSeq2[j-1]*25+numSeq1[i]*5+numSeq1[i-1]];
	H1=atpH[numSeq1[i]*5+numSeq2[j]]+tstack2Enthalpies[numSeq2[j]*125+numSeq2[j-1]*25+numSeq1[i]*5+numSeq1[i-1]];
	if(!isFinite(H1))
	{
		H1=_INFINITY;
		S1=-1.0;
	}
// If there is two dangling ends at the same end of duplex
	if(isFinite(dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]])&&isFinite(dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]]))
	{
		S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]]+dangleEntropies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]]+dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		if(!isFinite(H2))
		{
			H2=_INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(isFinite(H1))
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
	else if(isFinite(dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]]))
	{
		S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]];
		H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]];
		if(!isFinite(H2))
		{
			H2=_INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(isFinite(H1))
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
	else if(isFinite(dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]]))
	{
		S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
		if(!isFinite(H2))
		{
			H2=_INFINITY;
			S2=-1.0;
		}
		T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
		if(isFinite(H1))
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

	S2=atpS[numSeq1[i]*5+numSeq2[j]];
	H2=atpH[numSeq1[i]*5+numSeq2[j]];
	T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
	if(isFinite(H1))
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

void maxTM(int i,int j,double stackEntropies[],double stackEnthalpies[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	double T0,T1,S0,S1,H0,H1;

	S0=entropyDPT[(i-1)*Initint[2]+j-1];
	H0=enthalpyDPT[(i-1)*Initint[2]+j-1];
	T0=(H0+Initdouble[0])/(S0+Initdouble[1]+Initdouble[2]); // at current position 
	if(isFinite(enthalpyDPT[(i-2)*Initint[2]+j-2])&&isFinite(Hs(i-1,j-1,1,stackEnthalpies,Initint,numSeq1,numSeq2)))
	{
		S1=(entropyDPT[(i-2)*Initint[2]+j-2]+Ss(i-1,j-1,1,stackEntropies,Initint,numSeq1,numSeq2));
		H1=(enthalpyDPT[(i-2)*Initint[2]+j-2]+Hs(i-1,j-1,1,stackEnthalpies,Initint,numSeq1,numSeq2));
	}
	else
	{
		S1=-1.0;
		H1=_INFINITY;
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

void calc_bulge_internal(int i,int j,int ii,int jj,double* EntropyEnthalpy,int traceback,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
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
				H=bulgeLoopEnthalpies[loopSize]+stackEnthalpies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
				S=bulgeLoopEntropies[loopSize]+stackEntropies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
			}
			H+=enthalpyDPT[(i-1)*Initint[2]+j-1];
			S+=entropyDPT[(i-1)*Initint[2]+j-1];
			if(!isFinite(H))
			{
				H=_INFINITY;
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
			H=bulgeLoopEnthalpies[loopSize]+atpH[numSeq1[i]*5+numSeq2[j]]+atpH[numSeq1[ii]*5+numSeq2[jj]];
			H+=enthalpyDPT[(i-1)*Initint[2]+j-1];

			S=bulgeLoopEntropies[loopSize]+atpS[numSeq1[i]*5+numSeq2[j]]+atpS[numSeq1[ii]*5+numSeq2[jj]];
			S+=entropyDPT[(i-1)*Initint[2]+j-1];
			if(!isFinite(H))
			{
				H=_INFINITY;
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
		S=stackint2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+stackint2Entropies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		S+=entropyDPT[(i-1)*Initint[2]+j-1];

		H=stackint2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+stackint2Enthalpies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		H+=enthalpyDPT[(i-1)*Initint[2]+j-1];
		if(!isFinite(H))
		{
			H=_INFINITY;
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
		H=interiorLoopEnthalpies[loopSize]+tstackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+tstackEnthalpies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
		H+=enthalpyDPT[(i-1)*Initint[2]+j-1];

		S=interiorLoopEntropies[loopSize]+tstackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+tstackEntropies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
		S+=entropyDPT[(i-1)*Initint[2]+j-1];
		if(!isFinite(H))
		{
			H=_INFINITY;
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

void fillMatrix(double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
	int d,i,j,ii,jj;
	double SH[2];

	for(i=1;i<=Initint[0];++i)
	{
		for(j=1;j<=Initint[1];++j)
		{
			if(isFinite(enthalpyDPT[(i-1)*Initint[2]+j-1]))
			{
				SH[0]=-1.0;
				SH[1]=_INFINITY;
				LSH(i,j,SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

				if(isFinite(SH[1]))
				{
					entropyDPT[(i-1)*Initint[2]+j-1]=SH[0];
					enthalpyDPT[(i-1)*Initint[2]+j-1]=SH[1];
				}
				if(i>1&&j>1)
				{
					maxTM(i,j,stackEntropies,stackEnthalpies,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
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
							if(isFinite(enthalpyDPT[(ii-1)*Initint[2]+jj-1]))
							{
								SH[0]=-1.0;
								SH[1]=_INFINITY;
								calc_bulge_internal(ii,jj,i,j,SH,0,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

								if(SH[0]<-2500.0)
								{
									SH[0] =-3224.0;
									SH[1] = 0.0;
								}
								if(isFinite(SH[1]))
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

#include "thal.h"

#ifndef THAL_EXIT_ON_ERROR
#define THAL_EXIT_ON_ERROR 0
#endif

#define CHECK_ERROR(COND,MSG) if (COND) { strcpy(o->msg, MSG); errno = 0; longjmp(_jmp_buf, 1); }
#define THAL_OOM_ERROR { strcpy(o->msg, "Out of memory"); errno = ENOMEM; longjmp(_jmp_buf, 1); }
#define THAL_IO_ERROR(f) { sprintf(o->msg, "Unable to open file %s", f); longjmp(_jmp_buf, 1); }

#define STR(X) #X
#define LONG_SEQ_ERR_STR(MAX_LEN) "Target sequence length > maximum allowed (" STR(MAX_LEN) ") in thermodynamic alignment"
#define XSTR(X) STR(X)

/*
#ifdef INTEGER
# define isFinite(x) (x < _INFINITY / 2)
#else
# define isFinite(x) finite(x)
#endif
*/
/*** BEGIN CONSTANTS ***/
/*
# ifdef INTEGER
const double _INFINITY = 999999.0;
# else
# ifdef INFINITY
const double _INFINITY = INFINITY;
# else
const double _INFINITY = 1.0 / 0.0;
# endif
# endif
*/

/* get thermodynamic tables */

/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end (Left hand/3' end - Right end) */
static void RSH(int i,int j, double* EntropyEnthalpy,double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],char *numSeq1,char *numSeq2);

static void reverse(unsigned char *s);

/* Is sequence symmetrical */
static int symmetry_thermo(const unsigned char* seq);

/* traceback for dimers */
static void traceback(int i, int j, int* ps1, int* ps2, int maxLoop, thal_results* o,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2);

/* prints ascii output of dimer structure */
static void drawDimer(int*, int*,double, double, int, double, thal_results *,double Initdouble[],int Initint[],char *oligo1,char *oligo2);

static void strcatc(char*, char);


static jmp_buf _jmp_buf;

/* central method: execute all sub-methods for calculating secondary
   structure for dimer or for monomer */
void thal(const unsigned char *oligo_f,const unsigned char *oligo_r,const thal_args *a,thal_results *o,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[])
{
	double *SH,Initdouble[4];//0 is dplx_init_H, 1 is dplx_init_S, 2 is RC, 3 is SHleft
	int Initint[5]; //0 is len1, 1 is len2, 2 is len3, 3 is bestI, 4 is bestJ
	int i, j;
	int len_f, len_r;
	double T1,*enthalpyDPT,*entropyDPT,*send5,*hend5,result_TH;
	int k;
	int *bp;
	char *oligo1,*oligo2,*numSeq1,*numSeq2;
	unsigned char *oligo2_rev = NULL;
	double mh, ms;

	send5 = hend5 = NULL;
	enthalpyDPT = entropyDPT = NULL;
	numSeq1 = numSeq2 = NULL;
	oligo1 = oligo2 = NULL;
	strcpy(o->msg, "");
	o->temp = THAL_ERROR_SCORE;
	errno = 0; 

	len_f =strlen(oligo_f);
	len_r =strlen(oligo_r);

	o->align_end_1 = -1;
	o->align_end_2 = -1;

	if(a->type!=3)
	{
		oligo1 = (unsigned char*)malloc((len_f + 1) * sizeof(unsigned char));
		oligo2 = (unsigned char*)malloc((len_r + 1) * sizeof(unsigned char));
		strcpy((char*)oligo1,(const char*)oligo_f);
		strcpy((char*)oligo2,(const char*)oligo_r);
	}
	else
	{
		oligo1 = (unsigned char*)malloc((len_r + 1) * sizeof(unsigned char));
		oligo2 = (unsigned char*)malloc((len_f + 1) * sizeof(unsigned char));
		strcpy((char*)oligo1,(const char*)oligo_r);
		strcpy((char*)oligo2,(const char*)oligo_f);
	}
/*** INIT values for unimolecular and bimolecular structures ***/
	if (a->type==4) /* unimolecular folding */
	{
		Initint[1] =strlen(oligo2);
		Initint[2] = Initint[1] -1;
		Initdouble[0]= 0.0;
		Initdouble[1] = -0.00000000001;
		Initdouble[2]=0;
	}
	else /* hybridization of two oligos */
	{
		Initdouble[0]= 200;
		Initdouble[1]= -5.7;
		if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2))
			Initdouble[2]=1.9872* log(38/1000000000.0);
		else
			Initdouble[2]=1.9872* log(38/4000000000.0);
		
		if(a->type!=3)
		{
			oligo2_rev = (unsigned char*)malloc((strlen(oligo_r) + 1) * sizeof(unsigned char));
			strcpy((char*)oligo2_rev,(const char*)oligo_r);
		}
		else
		{
			oligo2_rev = (unsigned char*)malloc((strlen(oligo_f) + 1) * sizeof(unsigned char));
			strcpy((char*)oligo2_rev,(const char*)oligo_f);
		}
		reverse(oligo2_rev); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
		free(oligo2);
		oligo2=NULL;
		oligo2=&oligo2_rev[0];
	}
	Initint[0] =strlen(oligo1);
	Initint[1] =strlen(oligo2);
/* convert nucleotides to numbers */
	numSeq1 = (unsigned char*)realloc(numSeq1,Initint[0]+2);
	numSeq2 = (unsigned char*)realloc(numSeq2,Initint[1]+2);

	if(a->type == 4) /* monomer */
	{
	/* terminal basepairs */
		send5 = (double*)realloc(send5,(Initint[0]+1)*sizeof(double));
		hend5 = (double*)realloc(hend5,(Initint[0]+1)*sizeof(double));
	}
	for(i = 0; i < Initint[0]; i++)
		oligo1[i] = toupper(oligo1[i]);
	for(i = 0; i < Initint[1]; i++)
		oligo2[i] = toupper(oligo2[i]);
 	for(i = 1; i <= Initint[0]; ++i)
		numSeq1[i] = str2int(oligo1[i - 1]);
	for(i = 1; i <= Initint[1]; ++i)
		numSeq2[i] = str2int(oligo2[i - 1]);
	numSeq1[0] = numSeq1[Initint[0] + 1] = numSeq2[0] = numSeq2[Initint[1] + 1] = 4; /* mark as N-s */

	if (a->type==4) /* calculate structure of monomer */
	{
		enthalpyDPT =(double *)realloc(enthalpyDPT,Initint[0]*Initint[1]*sizeof(double));
		entropyDPT =(double *)realloc(entropyDPT,Initint[0]*Initint[1]*sizeof(double));
		initMatrix2(Initint,enthalpyDPT,entropyDPT,numSeq1);
		fillMatrix2(stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
		calc_terminal_bp(a->temp,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
		mh=hend5[Initint[0]];
		ms=send5[Initint[0]];
		o->align_end_1 = (int) mh;
		o->align_end_2 = (int) ms;
		bp=(int *)malloc(Initint[0]*sizeof(int));
		memset(bp,'\0',Initint[0]*sizeof(int));
		for (k = 0; k < Initint[0]; ++k)
			bp[k] = 0;
		if(isFinite(mh))
		{
			tracebacku(bp,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,numSeq2);
			o->temp=drawHairpin(bp,mh,ms,Initint);
		}

		if(o->temp==-_INFINITY && (!strcmp(o->msg, "")))
			o->temp=0.0;
		free(bp);
		free(enthalpyDPT);
		free(entropyDPT);
		free(numSeq1);
		free(numSeq2);
		free(send5);
		free(hend5);
		free(oligo1);
		free(oligo2);
		return;
	}
	else if(a->type!=4) /* Hybridization of two moleculs */
	{
		Initint[2]=Initint[1];
		enthalpyDPT =(double *)realloc(enthalpyDPT,Initint[0]*Initint[1]*sizeof(double)); /* dyn. programming table for dS and dH */
		entropyDPT =(double *)realloc(entropyDPT,Initint[0]*Initint[1]*sizeof(double)); /* enthalpyDPT is 3D array represented as 1D array */
		initMatrix(Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
		fillMatrix(stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

		Initdouble[3] = -_INFINITY;
		SH = (double*)malloc(2 * sizeof(double));
	/* calculate terminal basepairs */
		Initint[3] = Initint[4] = 0;
		if(a->type==1)
			for (i = 1; i <= Initint[0]; i++)
			{
				for (j = 1; j <= Initint[1]; j++)
				{
					RSH(i, j, SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,numSeq1,numSeq2);
					SH[0] = SH[0]+0.000001; /* this adding is done for compiler, optimization -O2 vs -O0 */
					SH[1] = SH[1]+0.000001;
					T1 = ((enthalpyDPT[(i-1)*Initint[2]+j-1]+ SH[1] +Initdouble[0]) / ((entropyDPT[(i-1)*Initint[2]+j-1]) + SH[0] +Initdouble[1] + Initdouble[2])) -273.15;
					if (T1 > Initdouble[3]  && ((entropyDPT[(i-1)*Initint[2]+j-1]+ SH[0])<0 && (SH[1] + enthalpyDPT[(i-1)*Initint[2]+j-1])<0))
					{
						Initdouble[3] = T1;
						Initint[3] = i;
						Initint[4] = j;
					}
				}
			}
		int *ps1, *ps2;
		ps1=(int *)malloc(Initint[0]*sizeof(int));
		ps2=(int *)malloc(Initint[1]*sizeof(int));
		for (i = 0; i < Initint[0]; ++i)
			ps1[i] = 0;
		for (j = 0; j < Initint[1]; ++j)
			ps2[j] = 0;
		if(a->type == 2 || a->type == 3)
		{
		 /* THAL_END1 */
			Initint[3] = Initint[4] = 0;
			Initint[3] = Initint[0];
			i = Initint[0];
			Initdouble[3] = -_INFINITY;
			for (j = 1; j <= Initint[1]; ++j)
			{
				RSH(i, j, SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,numSeq1,numSeq2);
				SH[0] = SH[0]+0.000001; /* this adding is done for compiler, optimization -O2 vs -O0,
					   that compiler could understand that SH is changed in this cycle */
				SH[1] = SH[1]+0.000001;
				T1 = ((enthalpyDPT[(i-1)*Initint[2]+j-1]+ SH[1] +Initdouble[0]) / ((entropyDPT[(i-1)*Initint[2]+j-1]) + SH[0] +Initdouble[1]+ Initdouble[2])) -273.15;
				if (T1 > Initdouble[3] && ((SH[0] +entropyDPT[(i-1)*Initint[2]+j-1])<0 && (SH[1] + enthalpyDPT[(i-1)*Initint[2]+j-1])<0))
				{
					Initdouble[3] = T1;
					Initint[4] = j;
				}
			}
		}
		if (!isFinite(Initdouble[3]))
			Initint[3] = Initint[4] = 1;
		double dH, dS;
		RSH(Initint[3], Initint[4], SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,numSeq1,numSeq2);

		dH = enthalpyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1]+ SH[1] +Initdouble[0];
		dS = (entropyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1] + SH[0] +Initdouble[1]);
	 /* tracebacking */
		for (i = 0; i < Initint[0]; ++i)
			ps1[i] = 0;
		for (j = 0; j < Initint[1]; ++j)
			ps2[j] = 0;
		if(isFinite(enthalpyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1]))
		{
			traceback(Initint[3], Initint[4],ps1, ps2, a->maxLoop, o,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
			drawDimer(ps1, ps2, dH, dS, a->temponly,a->temp, o,Initdouble,Initint,oligo1,oligo2);
			o->align_end_1=Initint[3];
			o->align_end_2=Initint[4];
		}
		else
		{
			o->temp = 0.0;
		}
		free(ps1);
		free(ps2);
		free(SH);
		free(oligo2_rev);
		free(enthalpyDPT);
		free(entropyDPT);
		free(numSeq1);
		free(numSeq2);
		free(oligo1);
		return;
	}
	return;
}


/* Set default args */
void 
set_thal_default_args(thal_args *a)
{
   memset(a, 0, sizeof(*a));
   a->debug = 0;
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop =30;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.8; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = 310.15; /* Kelvin */
   a->temponly = 1; /* return only melting temperature of predicted structure */
   a->dimer = 1; /* by default dimer structure is calculated */
}

/* Set default args for oligo */
void
set_thal_oligo_default_args(thal_args *a)    
{
   memset(a, 0, sizeof(*a));
   a->debug = 0;
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop =30;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.0; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = 310.15; /* Kelvin */
   a->temponly = 1; /* return only melting temperature of predicted structure */
   a->dimer = 1; /* by default dimer structure is calculated */
}


static void 
reverse(unsigned char *s)
{
   int i,j;
   char c;
   for (i = 0, j =strlen(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

#define INIT_BUF_SIZE 1024

static char*
p3_read_line(FILE *file, thal_results* o)
{
  static size_t ssz;
  static char *s = NULL;

  size_t remaining_size;
  char *p, *n;

  if (NULL == s) {
    ssz = INIT_BUF_SIZE;
    s = (char *) malloc(ssz);
  }
  p = s;
  remaining_size = ssz;
  while (1) {
    if (fgets(p, remaining_size, file) == NULL) /* End of file. */
      return p == s ? NULL : s;

    if ((n = strchr(p, '\n')) != NULL) {
      *n = '\0';
      return s;
    }

    /* We did not get the whole line. */

    if (ssz >= INT_MAX / 2)
      ssz = INT_MAX;
    else {
      ssz *= 2;
    }
    s = (char *)realloc(s,ssz);
    p = strchr(s, '\0');
    remaining_size = ssz - (p - s);
  }
}

static void 
RSH(int i, int j, double* EntropyEnthalpy,double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],char *numSeq1,char *numSeq2)
{
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   H1 = H2 = _INFINITY;
   T1 = T2 = -_INFINITY;
	if (numSeq1[i]+numSeq2[j]!=3) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 =atpS[numSeq1[i]*5+numSeq2[j]] + tstack2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
   H1 =atpH[numSeq1[i]*5+numSeq2[j]]+ tstack2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
   if(!isFinite(H1)) {
      H1 = _INFINITY;
      S1 = -1.0;
   }
   if(isFinite(dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]])) {
      S2 =atpS[numSeq1[i]*5+numSeq2[j]] + dangleEntropies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]] +
	dangleEntropies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
      H2 =atpH[numSeq1[i]*5+numSeq2[j]]+ dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]] +
	dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
      if(!isFinite(H2)) {
	 H2 = _INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 +Initdouble[0]) / (S2 +Initdouble[1]+ Initdouble[2]);
      if(isFinite(H1)) {
	 T1 = (H1 +Initdouble[0]) / (S1 +Initdouble[1]+ Initdouble[2]);
	 if(T1 < T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }

   if(isFinite(dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]])) {
      S2 =atpS[numSeq1[i]*5+numSeq2[j]] +dangleEntropies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]];
      H2 =atpH[numSeq1[i]*5+numSeq2[j]]+ dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]];
      if(!isFinite(H2)) {
	 H2 = _INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 +Initdouble[0]) / (S2 +Initdouble[1]+ Initdouble[2]);
      if(isFinite(H1)) {
	 T1 = (H1 +Initdouble[0]) / (S1 +Initdouble[1]+ Initdouble[2]);
	 if(T1 < T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }

   if(isFinite(dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]])) {
      S2 =atpS[numSeq1[i]*5+numSeq2[j]]+ dangleEntropies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
      H2 =atpH[numSeq1[i]*5+numSeq2[j]]+ dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
      if(!isFinite(H2)) {
	 H2 = _INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 +Initdouble[0]) / (S2 +Initdouble[1]+ Initdouble[2]);
      if(isFinite(H1)) {
	 T1 = (H1 +Initdouble[0]) / (S1 +Initdouble[1]+ Initdouble[2]);
	 if(T1 < T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }
   S2 =atpS[numSeq1[i]*5+numSeq2[j]];
   H2 =atpH[numSeq1[i]*5+numSeq2[j]];
   T2 = (H2 +Initdouble[0]) / (S2 +Initdouble[1]+ Initdouble[2]);
   if(isFinite(H1)) {
      if(T1 < T2) {
	 EntropyEnthalpy[0] = S2;
	 EntropyEnthalpy[1] = H2;
      } else {
	 EntropyEnthalpy[0] = S1;
	 EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}

/* Return 1 if string is symmetrical, 0 otherwise. */
static int 
symmetry_thermo(const unsigned char* seq)
{
   register char s;
   register char e;
   const unsigned char *seq_end=seq;
   int i = 0;
   int seq_len=strlen(seq);
   int mp = seq_len/2;
   if(seq_len%2==1) {
      return 0;
   }
   seq_end+=seq_len;
   seq_end--;
   while(i<mp) {
      i++;
      s=toupper(*seq);
      e=toupper(*seq_end);
      if ((s=='A' && e!='T')
	  || (s=='T' && e!='A')
	  || (e=='A' && s!='T')
	  || (e=='T' && s!='A')) {
	 return 0;
      }
      if ((s=='C' && e!='G')
	  || (s=='G' && e!='C')
	  || (e=='C' && s!='G')
	  || (e=='G' && s!='C')) {
	 return 0;
      }
      seq++;
      seq_end--;
   }
   return 1;
}


static void 
traceback(int i, int j, int* ps1, int* ps2, int maxLoop, thal_results* o,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double *enthalpyDPT,double *entropyDPT,char *numSeq1,char *numSeq2)
{
   int d, ii, jj, done;
   double* SH;
   SH = (double*)malloc(2 * sizeof(double));
   ps1[i - 1] = j;
   ps2[j - 1] = i;
   while(1) {
      SH[0] = -1.0;
      SH[1] = _INFINITY;
      LSH(i,j,SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

      if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH[0]) && equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH[1])) {
	 break;
      }
      done = 0;
      if (i > 1 && j > 1 && equal(entropyDPT[(i-1)*Initint[2]+j-1], Ss(i - 1, j - 1, 1,stackEntropies,Initint,numSeq1,numSeq2) +entropyDPT[(i-2)*Initint[2]+j-2])) {
	 i = i - 1;
	 j = j - 1;
	 ps1[i - 1] = j;
	 ps2[j - 1] = i;
	 done = 1;
      }
      for (d = 3; !done && d <= maxLoop + 2; ++d) {
	 ii = i - 1;
	 jj = -ii - d + (j + i);
	 if (jj < 1) {
	    ii -= abs(jj-1);
	    jj = 1;
	 }
	 for (; !done && ii > 0 && jj < j; --ii, ++jj) {
	    SH[0] = -1.0;
	    SH[1] = _INFINITY;
	    calc_bulge_internal(ii, jj, i, j, SH,1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

	    if (equal(entropyDPT[(i-1)*Initint[2]+j-1], SH[0]) && equal(enthalpyDPT[(i-1)*Initint[2]+j-1], SH[1])) {
	       i = ii;
	       j = jj;
	       ps1[i - 1] = j;
	       ps2[j - 1] = i;
	       done = 1;
	       break;
	    }
	 }
      }
   }
   free(SH);
}

static void 
drawDimer(int* ps1, int* ps2, double H, double S, int temponly, double t37, thal_results *o,double Initdouble[],int Initint[],char *oligo1,char *oligo2)
{
   int i, j, k, numSS1, numSS2, N;
   char* duplex[4];
   double G, t;
   t = G = 0;
   if (!isFinite(Initdouble[3])){
      if(temponly==0) {
	 printf("No predicted secondary structures for given sequences\n");
      }
      o->temp = 0.0; /* lets use generalization here; this should rather be very negative value */
      strcpy(o->msg, "No predicted sec struc for given seq");
      return;
   } else {
      N=0;
      for(i=0;i<Initint[0];i++){
	 if(ps1[i]>0) ++N;
      }
      for(i=0;i<Initint[1];i++) {
	 if(ps2[i]>0) ++N;
      }
      N = (N/2) -1;
	t = ((H) / (S + (N *-0.51986) + Initdouble[2])) -273.15;
      if(temponly==0) {
	G = (H) - (t37 * (S + (N *-0.51986)));
	S = S + (N *-0.51986);
	 o->temp = (double) t;
	 /* maybe user does not need as precise as that */
	 /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n",
		Initint[0], (double) S, (double) H, (double) G, (double) t, (int) N, saltCorrection, RC); */
	 printf("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
		(double) S, (double) H, (double) G, (double) t);
      } else {
	 o->temp = (double) t;
	 return;
      }
   }

   duplex[0] = (char*)malloc(Initint[0] + Initint[1] + 1);
   duplex[1] = (char*)malloc(Initint[0] + Initint[1] + 1);
   duplex[2] = (char*)malloc(Initint[0] + Initint[1] + 1);
   duplex[3] = (char*)malloc(Initint[0] + Initint[1] + 1);
   duplex[0][0] = duplex[1][0] = duplex[2][0] = duplex[3][0] = 0;

   i = 0;
   numSS1 = 0;
   while (ps1[i++] == 0) ++numSS1;
   j = 0;
   numSS2 = 0;
   while (ps2[j++] == 0) ++numSS2;

   if (numSS1 >= numSS2){
      for (i = 0; i < numSS1; ++i) {
	 strcatc(duplex[0], oligo1[i]);
	 strcatc(duplex[1], ' ');
	 strcatc(duplex[2], ' ');
      }
      for (j = 0; j < numSS1 - numSS2; ++j) strcatc(duplex[3], ' ');
      for (j = 0; j < numSS2; ++j) strcatc(duplex[3], oligo2[j]);
   } else {
      for (j = 0; j < numSS2; ++j) {
	 strcatc(duplex[3], oligo2[j]);
	 strcatc(duplex[1], ' ');
	 strcatc(duplex[2], ' ');
      }
      for (i = 0; i < numSS2 - numSS1; ++i)
	strcatc(duplex[0], ' ');
      for (i = 0; i < numSS1; ++i)
	strcatc(duplex[0], oligo1[i]);
   }
   i = numSS1 + 1;
   j = numSS2 + 1;

   while (i <= Initint[0]) {
      while (i <= Initint[0] && ps1[i - 1] != 0 && j <= Initint[1] && ps2[j - 1] != 0) {
	 strcatc(duplex[0], ' ');
	 strcatc(duplex[1], oligo1[i - 1]);
	 strcatc(duplex[2], oligo2[j - 1]);
	 strcatc(duplex[3], ' ');
	 ++i;
	 ++j;
      }
      numSS1 = 0;
      while (i <= Initint[0] && ps1[i - 1] == 0) {
	 strcatc(duplex[0], oligo1[i - 1]);
	 strcatc(duplex[1], ' ');
	 ++numSS1;
	 ++i;
      }
      numSS2 = 0;
      while (j <= Initint[1] && ps2[j - 1] == 0) {
	 strcatc(duplex[2], ' ');
	 strcatc(duplex[3], oligo2[j - 1]);
	 ++numSS2;
	 ++j;
      }
      if (numSS1 < numSS2)
	for (k = 0; k < numSS2 - numSS1; ++k) {
	   strcatc(duplex[0], '-');
	   strcatc(duplex[1], ' ');
	}
      else if (numSS1 > numSS2)
	for (k = 0; k < numSS1 - numSS2; ++k) {
	   strcatc(duplex[2], ' ');
	   strcatc(duplex[3], '-');
	}
   }
   printf("SEQ\t");
   printf("%s\n", duplex[0]);
   printf("SEQ\t");
   printf("%s\n", duplex[1]);
   printf("STR\t");
   printf("%s\n", duplex[2]);
   printf("STR\t");
   printf("%s\n", duplex[3]);

   free(duplex[0]);
   free(duplex[1]);
   free(duplex[2]);
   free(duplex[3]);

   return;
}

static void 
strcatc(char* str, char c)
{
   str[strlen(str) + 1] = 0;
   str[strlen(str)] = c;
}

main()
{
	double stackEntropies[625],stackEnthalpies[625],stackint2Entropies[625],stackint2Enthalpies[625];
	double dangleEntropies3[125],dangleEnthalpies3[125],dangleEntropies5[125],dangleEnthalpies5[125];
	double hairpinLoopEntropies[30],interiorLoopEntropies[30],bulgeLoopEntropies[30],hairpinLoopEnthalpies[30],interiorLoopEnthalpies[30],bulgeLoopEnthalpies[30];
	double tstackEntropies[625],tstackEnthalpies[625],tstack2Entropies[625],tstack2Enthalpies[625];
	char *triloopEntropies1,*triloopEnthalpies1,*tetraloopEntropies1,*tetraloopEnthalpies1;
	double *triloopEntropies2,*triloopEnthalpies2,*tetraloopEntropies2,*tetraloopEnthalpies2;
	int numTriloops,numTetraloops;
	double atpS[25],atpH[25];

        char path[100]="/home/bjia/GPU-LAMP/Create_Primer3/primer3-2.3.5/src/primer3_config/";
        char one[30]="GAGCTAGAGTCGTTAGCTAAACC";
        char two[30]="GAGCTAGAGTCGTTAGCTAAACC";
        thal_results *o;
        thal_args *a;

        o=(thal_results *)malloc(sizeof(thal_results));
        memset(o,'\0',sizeof(thal_results));
        a=(thal_args *)malloc(sizeof(thal_args));              
        memset(a,'\0',sizeof(thal_args));

        a->debug=0;
        a->maxLoop=30;
        a->mv=50;
        a->dv=4;
        a->dntp=1.4;
        a->dna_conc=38;
        a->temp=310.15; 
        a->temponly=1;

//read_parameter
	getStack(stackEntropies,stackEnthalpies,path);
	getStackint2(stackint2Entropies,stackint2Enthalpies,path);
	getDangle(dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,path);
	getLoop(hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,path);
	getTstack(tstackEntropies,tstackEnthalpies,path);
	getTstack2(tstack2Entropies,tstack2Enthalpies,path);

	numTriloops=get_num_line(path,0);
	triloopEntropies1=(char *)malloc(numTriloops*5);
	triloopEnthalpies1=(char *)malloc(numTriloops*5);
	triloopEntropies2=(double *)malloc(numTriloops*sizeof(double));
        triloopEnthalpies2=(double *)malloc(numTriloops*sizeof(double));
	getTriloop(triloopEntropies1,triloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,path);
	
	numTetraloops=get_num_line(path,1);
	tetraloopEntropies1=(char *)malloc(numTetraloops*6);
	tetraloopEnthalpies1=(char *)malloc(numTetraloops*6);
	tetraloopEntropies2=(double *)malloc(numTetraloops*sizeof(double));
	tetraloopEnthalpies2=(double *)malloc(numTetraloops*sizeof(double));
	getTetraloop(tetraloopEntropies1,tetraloopEnthalpies1,tetraloopEntropies2,tetraloopEnthalpies2,path);
	tableStartATS(6.9,atpS);
	tableStartATH(2200.0,atpH);

printf("single-any: real is 11.37, ours is ");
	a->type=1;
	a->dimer=1;	
        thal(one,two,a,o,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
	printf("%lf\n",o->temp);
printf("single-end: real is 0.65, ours is ");
        a->type=2;
        a->dimer=1;
        thal(one,two,a,o,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
        printf("%lf\n",o->temp);
printf("single-hairpin: real is 35.66, ours is ");
        a->type=4;
        a->dimer=0;
        thal(one,two,a,o,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
        printf("%lf\n",o->temp);

	memset(one,'\0',30);
	strcpy(one,"TTAGACTTCTAAGCGCTGTGAAC");
	memset(two,'\0',30);
	strcpy(two,"TGGTTCACGTATGCCTGC");
printf("two-any: real is 1.78, ours is ");
	a->type=1;
	a->dimer=1;
	thal(one,two,a,o,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
        printf("%lf\n",o->temp);

printf("two-end: real is 6.76, ours is ");
	a->type=3;
	a->dimer=1;
	thal(one,two,a,o,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
        printf("%lf\n",o->temp);
//when calculate the compl_end, a->type=2 and 3, inputs are (one, two) and (rev_one,rev_two :5'-3'), total four times, select the biggest
        free(o);
        free(a);
	free(triloopEntropies1);
	free(triloopEnthalpies1);
	free(tetraloopEntropies1);
	free(tetraloopEnthalpies1);
	free(triloopEntropies2);
	free(triloopEnthalpies2);
	free(tetraloopEntropies2);
	free(tetraloopEnthalpies2);
}
