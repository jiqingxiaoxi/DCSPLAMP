#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/stat.h>
#include<regex.h>

struct Node
{
	int pos;
	int gi;
	int plus;  //as a flag, 1 is OK, 0 is no
	int minus; //as a flag
	struct Node *next;
};

struct INFO
{
	char name[301];
	int turn;
	struct INFO *next;
};

struct Primer
{
	int pos;
	int len;
	int plus;
	int minus;
	int total;  //how many GIs can be used
	struct Primer *self;
	struct Primer *loop;
	struct Primer *notloop;
	struct Primer *next;
	struct Node *common;
	struct Node *special;
};

//get the file size
int file_size2(char *filename)
{
	int size;
        struct stat statbuf;

        stat(filename,&statbuf);
        size=statbuf.st_size;
        return size;
}

void usage()
{
        printf("Usage:\n");
        printf("    LAMP_CPU  -in <name>  -out <result>  -high[-low] [options]*\n\n");
        printf("    -in   <string>:  the name of candidate single primers file\n");
        printf("    -out  <string>:  the result file of LAMP primers\n");
        printf("    -dir  <string>:  the directory to store candidate single primers, default is current directory\n");
        printf("    -ref  <string>:  the reference sequence file used in single program, fasta format\n");
	printf("    -expect  <int>:  the number of LAMP primers needed to be design, default is 10\n"); 
        printf("    -common:         design common LAMP primers\n");
	printf("    -special:        design special LAMP primers\n");
        printf("    -P    <string>:  the primer3_core program in Primer3 software, use primer3 to check the complementarity among primers.\n");
        printf("    -high/-low:      design candidate single primers in high/low GC region, high: the GC content>=45%%; low: the GC content <=45%%.\n");
        printf("    -loop:           design LAMP primer with loop primers\n");
        printf("    -h/-help:        usage\n");
}

void generate_primer(char *seq,char primer[],int start,int length,int flag) //flag=0:plus
{
        int i;
        
        if(flag==0)
        {
                for(i=0;i<length;i++)
                        primer[i]=seq[start+i];
                primer[i]='\0';
        }
        else
        {
                for(i=0;i<length;i++)
                {
                        if((seq[start+length-1-i]=='A')||(seq[start+length-1-i]=='a'))
                                primer[i]='T';
                        else if((seq[start+length-1-i]=='T')||(seq[start+length-1-i]=='t'))
                                primer[i]='A';
                        else if((seq[start+length-1-i]=='C')||(seq[start+length-1-i]=='c'))
                                primer[i]='G';
                        else 
                                primer[i]='C';
                }
                primer[length]='\0';
        }
}

void take_regulate(regmatch_t pmatch[],int which,char *out,char *input)
{
        int i,j=0;
        for(i=pmatch[which].rm_so;i<pmatch[which].rm_eo;i++)
        {
                out[j]=input[i];
                j++;
        }
        out[j]='\0';
}

main(int argc,char **argv)
{
	int i,j,*result,new,*best_par,*record,expect,flag[11],common_num[1],turn,success;
	char *output,*prefix,*store_path,*path_fa,*inner,*outer,*loop,*primer3,*primer3_config;
	char *temp,*seq,F3[26],F2[26],F1c[26],B1c[26],B2[26],B3[26],LF[26],LB[26];
	FILE *fp;
	struct Primer *headL,*headS,*headLoop,*p_F3,*p_F2,*p_F1c,*p_B1c,*p_B2,*p_B3,*result_Loop[2];
	struct Node *p_node,*p_temp;
	struct INFO *headList,*p_list;
	time_t start,end;
	
	struct Primer *read_par();
	struct INFO *read_list();
	void next_one();
	int check_common();
	void how_many();
	void add();
	int check_add();
	int check_gc();
	int check_uniq();
	int check_structure();
	int design_loop();

	start=time(NULL);
        for(i=0;i<=10;i++)
        {
                flag[i]=0;
        }
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
                        i++;
                }
                else if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-help")==0)
                {
                        usage();
                        exit(1);
                }
                else if(strcmp(argv[i],"-P")==0)
                {
                        flag[7]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-P\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
			j=strlen(argv[i+1]);
                        primer3=(char *)malloc(j+1);
			memset(primer3,'\0',j+1);
                        strcpy(primer3,argv[i+1]);
                        i=i+2;

			if(access(primer3,0)==-1)
			{
				printf("Error! Don't have the %s program!\n",primer3);
				exit(1);
			}
                }
		else if(strcmp(argv[i],"-common")==0)
		{
			flag[5]=1;
			i++;
		}
		else if(strcmp(argv[i],"-special")==0)
		{
			flag[6]=1;
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
	if(flag[4]==0)
		expect=10;

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

        //primer3 program
        if(flag[7]==1)
        {
		j=strlen(primer3);
		j=j+20;
		primer3_config=(char *)malloc(j);
                memset(primer3_config,'\0',j);
                strcpy(primer3_config,primer3);
                i=j-20;
                i--;
                while(primer3_config[i]!='/'&&i>=0)
                {
                        primer3_config[i]='\0';
                        i--;
                }
                strcat(primer3_config,"primer3_config/");
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
		seq[j]=temp[i];
		i++;
		j++;
	}
       	free(temp);

//common-list
	if(flag[5])
	{
		headList=read_list(inner,common_num);
		common_num[0]++;
	}
	else
		common_num[0]=1;

//read parameters
	headS=read_par(outer,flag[5],flag[6]);
	headL=read_par(inner,flag[5],flag[6]);
	if(flag[10])
		headLoop=read_par(loop,flag[5],flag[6]);

//common statistics
	if(flag[5])
	{
		how_many(headL,common_num[0]);
		how_many(headS,common_num[0]);
		if(flag[10])
			how_many(headLoop,common_num[0]);
	}

//the next one 
	next_one(headL,headL,0);
	next_one(headS,headS,0);
	next_one(headS,headL,1);
	next_one(headL,headS,1);
	if(flag[10]==1)
	{
		next_one(headS,headLoop,2);
		next_one(headL,headLoop,2);
	}

	if(flag[5])
	{
		result=(int *)malloc(common_num[0]*sizeof(int));
		record=(int *)malloc(expect*(common_num[0]+1)*sizeof(int));
		for(i=0;i<expect*(1+common_num[0]);i++)
			record[i]=0;
	}
	best_par=(int *)malloc(expect*16*sizeof(int)); //include LF/LB
	for(i=0;i<expect*16;i++)    
        {
                best_par[i]=0;
        }
	end=time(NULL);
	printf("The prepare time is %0.1f seconds!\n",difftime(end,start));

//design LAMP primers
	start=time(NULL);
	turn=0;
	for(j=common_num[0];j>=1;j--)
	{
		if(turn>=expect)
			break;
        	for(p_F3=headS;p_F3;p_F3=p_F3->next)   //F3
        	{
			if(turn>=expect)
				break;
			if(p_F3->total<j)
				continue;
			if(p_F3->plus==0)
				continue;
			success=check_add(p_F3->pos,best_par,expect); //in best_par are sorted primers by common
			if(success==0)
				continue;
			new=0;  //whether find a new primer, if find, adjust F3 position
	                for(p_F2=p_F3->self;p_F2;p_F2=p_F2->next)  //F2
                	{
				if(new==1)
                                        break;
				if(p_F2->total<j)
                                        continue;
				if(p_F2->plus==0)
					continue;
	                        if(p_F2->pos-(p_F3->pos+p_F3->len)>20)
	                                break;
	                        for(p_F1c=p_F2->notloop;p_F1c;p_F1c=p_F1c->next)   //F1c
                        	{
					if(new==1)
						break;
					if(p_F1c->total<j)
						continue;
					if(p_F1c->minus==0)
						continue;
					if(p_F1c->pos-p_F2->pos-1<40)
						continue;
	                                if(p_F1c->pos-p_F2->pos-1>60)
                                        	break;
                                	for(p_B1c=p_F1c->self;p_B1c;p_B1c=p_B1c->next)   //B1c
                                	{
						if(new==1)
                                                        break;
                                                if(p_B1c->total<j)
                                                        continue;
						if(p_B1c->plus==0)
							continue;
	                                        if(p_B1c->pos-p_F1c->pos>85)
                                                	break;
                                        	for(p_B2=p_B1c->notloop;p_B2;p_B2=p_B2->next)   //B2
                                        	{
							if(new==1)
								break;
							if(p_B2->total<j)
								continue;
							if(p_B2->minus==0)
								continue;
							if((p_B2->pos+p_B2->len-1)-(p_B1c->pos+p_B1c->len-1)-1<40)
								continue;
	                                                if((p_B2->pos+p_B2->len-1)-(p_B1c->pos+p_B1c->len-1)-1>60)
                                                        	break;
                                                	if(p_B2->pos+p_B2->len-1-p_F2->pos-1<120)
                                                        	continue;
                                                	if(p_B2->pos+p_B2->len-1-p_F2->pos-1>180)
                                                        	break;
						//check whether has enough positions for loop
							if(flag[10]&&((p_F1c->pos-p_F2->pos-p_F2->len)<18)&&(p_B2->pos-p_B1c->pos-p_B1c->len<18))
								continue;
                                                	for(p_B3=p_B2->self;p_B3;p_B3=p_B3->next)  //B3
                                                	{
								if(p_B3->total<j)
                                                                        continue;
								if(p_B3->minus==0)
									continue;
	                                                        if(p_B3->pos-(p_B2->pos+p_B2->len)>20)
                                                                	break;
							//primer GC and target GC
								success=check_gc(seq,p_F3->pos,(p_B3->pos+p_B3->len),flag[8]);
								if(success==0)
									continue;
                                                        //check uniq and common
								if(flag[6])
								{
									success=check_uniq(p_F3,p_F2,p_F1c,p_B1c,p_B2,p_B3);
									if(success==0)
										continue;
								}
								if(flag[5])
								{
	                                                        	success=check_common(p_F3,p_F2,p_F1c,p_B1c,p_B2,p_B3,result,common_num[0]);
									if(success!=j)
										continue;
								}
							//check second structure
								if(flag[7])
								{
									generate_primer(seq,F3,p_F3->pos,p_F3->len,0);
									generate_primer(seq,F2,p_F2->pos,p_F2->len,0);
									generate_primer(seq,F1c,p_F1c->pos,p_F1c->len,1);
									generate_primer(seq,B1c,p_B1c->pos,p_B1c->len,0);
									generate_primer(seq,B2,p_B2->pos,p_B2->len,1);
									generate_primer(seq,B3,p_B3->pos,p_B3->len,1);
									success=check_structure(F3,F2,F1c,B1c,B2,B3,primer3,primer3_config,flag[8]);
									if(success==0)
										continue;
								}
							//design loop
								if(flag[10])
								{
									success=design_loop(p_F3,p_F2,p_F2->loop,p_F1c,p_B1c,p_B1c->loop,p_B2,p_B3,result_Loop,result,primer3,primer3_config,flag,common_num[0],seq);
									if(success==0)
										continue;
									add(p_F3,p_F2,p_F1c,p_B1c,p_B2,p_B3,result_Loop[0],result_Loop[1],best_par,turn);
								}
								else
									add(p_F3,p_F2,p_F1c,p_B1c,p_B2,p_B3,NULL,NULL,best_par,turn);

								if(flag[5])
								{
									for(i=0;i<common_num[0];i++)
										record[turn*(common_num[0]+1)+i]=result[i];
									record[turn*(common_num[0]+1)+i]=j;
								}
								new=1;
								turn++;
								break;
							} //B3
                                                }  //B2
                                        } //B1c
                                } //F1c
                        }  //F2
                }  //F3
        }
	end=time(NULL);
	printf("the time for design is %0.1f seconds!\n",difftime(end,start));

	fp=fopen(output,"w");
        if(fp==NULL)
        {
                printf("Error: can't create the %s file!\n",output);
                exit(1);
        }
        
        for(i=0;i<expect;i++)
        {
		if(best_par[i*16+1]==0)
			break;
		turn=i+1;
		fprintf(fp,"The %d LAMP primers:\n",turn);
		generate_primer(seq,F3,best_par[i*16],best_par[i*16+1],0);
		fprintf(fp,"  F3: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16],best_par[i*16+1],F3);
		generate_primer(seq,F2,best_par[i*16+2],best_par[i*16+3],0);
		fprintf(fp,"  F2: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16+2],best_par[i*16+3],F2);
		generate_primer(seq,F1c,best_par[i*16+4],best_par[i*16+5],1);
		fprintf(fp,"  F1c: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16+4],best_par[i*16+5],F1c);
		generate_primer(seq,B1c,best_par[i*16+6],best_par[i*16+7],0);
		fprintf(fp,"  B1c: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16+6],best_par[i*16+7],B1c);
		generate_primer(seq,B2,best_par[i*16+8],best_par[i*16+9],1);
		fprintf(fp,"  B2: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16+8],best_par[i*16+9],B2);
		generate_primer(seq,B3,best_par[i*16+10],best_par[i*16+11],1);
		fprintf(fp,"  B3: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16+10],best_par[i*16+11],B3);
		if(flag[10])
		{
			if(best_par[i*16+13]==0)
				fprintf(fp,"  LF: NULL\n");
			else
			{
				generate_primer(seq,LF,best_par[i*16+12],best_par[i*16+13],1);
				fprintf(fp,"  LF: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16+12],best_par[i*16+13],LF);
			}
			if(best_par[i*16+15]==0)
                                fprintf(fp,"  LB: NULL\n");
                        else
			{
				generate_primer(seq,LB,best_par[i*16+14],best_par[i*16+15],0);
                                fprintf(fp,"  LB: pos:%d,length:%d bp, primer(5'-3'):%s\n",best_par[i*16+14],best_par[i*16+15],LB);
			}
		}
		if(flag[5])
		{
			turn=0;
			fprintf(fp,"  This set of LAMP primers could be used in %d genomes, there are: ",record[i*(common_num[0]+1)+common_num[0]]);
			p_list=headList;
			for(j=0;j<common_num[0];j++)
			{
				if(record[i*(common_num[0]+1)+j]==0)
					continue;
				while(p_list)
				{
					if(p_list->turn==j)
						break;
					else
						p_list=p_list->next;
				}
				if(turn==0)
					fprintf(fp,"%s",p_list->name);
				else
					fprintf(fp,", %s",p_list->name);
				turn++;
			}
			fprintf(fp,"\n");
		}
        }
	fclose(fp);
	free(best_par);
	free(seq);
	free(output);
        free(prefix);
        free(store_path);
        free(inner);
        free(outer);
        free(path_fa);

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
		
		p_F3=headL->next;
		free(headL);
		headL=p_F3;
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

                p_F3=headS->next;
                free(headS);
                headS=p_F3;
        }

	if(flag[5])
	{
		while(headList)
		{
			p_list=headList->next;
			free(headList);
			headList=p_list;
		}
		free(record);
		free(result);
	}

	if(flag[7])
	{
		free(primer3);
		free(primer3_config);
	}

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

                	p_F3=headLoop->next;
                	free(headLoop);
                	headLoop=p_F3;
        	}
	}
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
		new_primer->plus=plus;
		new_primer->minus=minus;
		new_primer->next=NULL;
		new_primer->self=NULL;
		new_primer->loop=NULL;
		new_primer->notloop=NULL;
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
                	new_node->next=p_node;
		}
		fclose(fp);
		free(in);
	}
	return head;
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

//function: the next one
void next_one(struct Primer *first, struct Primer *second,int flag) //0:self,1:notloop;2:loop
{
	struct Primer *one,*two,*start;
	int pos=-1;

	one=first;
	start=second;
	two=start;

	while(one)
	{
		if(pos!=one->pos)
		{
			while(start)
			{
				if(start->pos+18<one->pos)
					start=start->next;
				else
					break;
			}
			pos=one->pos;
		}
		//move second
		two=start;
		while(two)
		{
			if(two->pos<one->pos+one->len)
				two=two->next;
			else
			{
				if(flag==0)
					one->self=two;
				else if(flag==1)
					one->notloop=two;
				else
					one->loop=two;
				break;
			}
		}
		one=one->next;
	}			
}

//function: check how many GIs this primer can be used for
int check_common(struct Primer *F3,struct Primer *F2,struct Primer *F1c,struct Primer *B1c,struct Primer *B2,struct Primer *B3,int *result,int common)
{
        int dis,num,i;
	struct Node *c_F3,*c_F2,*c_F1c,*c_B1c,*c_B2,*c_B3;

	for(i=0;i<common;i++)
		result[i]=0;
//plus
	for(c_F3=F3->common;c_F3;c_F3=c_F3->next)
	{
		if(c_F3->plus!=1)
			continue;
		if(result[c_F3->gi])
			continue;
		for(c_F2=F2->common;c_F2;c_F2=c_F2->next)
		{
			if(c_F2->gi!=c_F3->gi)
				continue;
			if(c_F2->plus!=1)
				continue;
			for(c_F1c=F1c->common;c_F1c;c_F1c=c_F1c->next)
			{
				if(c_F1c->gi!=c_F3->gi)
					continue;
				if(c_F1c->minus!=1)
					continue;
				for(c_B1c=B1c->common;c_B1c;c_B1c=c_B1c->next)
				{
					if(c_B1c->gi!=c_F3->gi)
						continue;
					if(c_B1c->plus!=1)
						continue;
					for(c_B2=B2->common;c_B2;c_B2=c_B2->next)
					{
						if(c_B2->gi!=c_F3->gi)
							continue;
						if(c_B2->minus!=1)
							continue;
						for(c_B3=B3->common;c_B3;c_B3=c_B3->next)
						{
							if(c_B3->gi!=c_F3->gi)
								continue;
							if(c_B3->minus!=1)
								continue;
						//F3-F2 
        						dis=c_F2->pos-(c_F3->pos+F3->len-1)-1;
						        if(dis<0)
                						continue;
        						if(dis>20)
                						continue;
						//F2-F1c
						        dis=c_F1c->pos-c_F2->pos-1;
        						if(dis<40)
                						continue;
        						if(dis>60)
                						continue;
						//F1c-B1c
        						dis=c_B1c->pos-(c_F1c->pos+F1c->len-1)-1;
        						if(dis<0)
                						continue;
						//B1c-B2
        						dis=(c_B2->pos+B2->len-1)-(c_B1c->pos+B1c->len-1)-1;
        						if(dis<40)
								continue;
						        if(dis>60)
                						continue;
						//F2-B2
        						dis=c_B2->pos+B2->len-1-c_F2->pos-1;
        						if(dis<120)
                						continue;
						        if(dis>180)
                						continue;
						//B2-B3
        						dis=c_B3->pos-(c_B2->pos+B2->len-1)-1;
        						if(dis<0)
                						continue;
						        if(dis>20)
                						continue;
							result[c_F3->gi]=1;
						}
					}
				}
			}
		}
	}
//minus
        for(c_F3=F3->common;c_F3;c_F3=c_F3->next)
        {
                if(c_F3->minus!=1)
                        continue;
                if(result[c_F3->gi])
                        continue;  //this GI can common

                for(c_F2=F2->common;c_F2;c_F2=c_F2->next)
                {
                        if(c_F2->gi!=c_F3->gi)
                                continue;
                        if(c_F2->minus!=1)
                                continue;
                        for(c_F1c=F1c->common;c_F1c;c_F1c=c_F1c->next)
                        {
                                if(c_F1c->gi!=c_F3->gi)
                                        continue;
                                if(c_F1c->plus!=1)
                                        continue;
                                for(c_B1c=B1c->common;c_B1c;c_B1c=c_B1c->next)
                                {
                                        if(c_B1c->gi!=c_F3->gi)
                                                continue;
                                        if(c_B1c->minus!=1)
                                                continue;
                                        for(c_B2=B2->common;c_B2;c_B2=c_B2->next)
                                        {
                                                if(c_B2->gi!=c_F3->gi)
                                                        continue;
                                                if(c_B2->plus!=1)
                                                        continue;
                                                for(c_B3=B3->common;c_B3;c_B3=c_B3->next)
                                                {
                                                        if(c_B3->gi!=c_F3->gi)
                                                                continue;
                                                        if(c_B3->plus!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=c_F3->pos-(c_F2->pos+F2->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=(c_F2->pos+F2->len-1)-(c_F1c->pos+F1c->len-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=c_F1c->pos-(c_B1c->pos+B1c->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=c_B1c->pos-c_B2->pos-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=c_F2->pos+F2->len-1-c_B2->pos-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=c_B2->pos-(c_B3->pos+B3->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                        result[c_F3->gi]=1;
                                                }
                                        }
                                }
                        }
                }
        }
	num=0;
	for(i=0;i<common;i++)
	{
		num=num+result[i];
	}
	return num;
}

int check_structure(char F3[],char F2[],char F1c[],char B1c[],char B2[],char B3[],char *primer3,char *primer3_config,int flag)
{
	int i,j,random,cflags,success;
	float TH;
	char *in,*out,*dir,*list[6],pattern[50],value[20],FIP[35],BIP[35],*line,*script;
	FILE *fp;
	regex_t reg;
        regmatch_t pmatch[2];

//prepare
	list[0]=F3;
	list[1]=F2;
	list[2]=F1c;
	list[3]=B1c;
	list[4]=B2;
	list[5]=B3;
	j=strlen(F1c);
	for(i=0;i<17;i++)
	{
		FIP[i]=F1c[j-17+i];
		FIP[i+17]=F2[i];
	}
	j=strlen(B1c);
	for(i=0;i<17;i++)
	{
		BIP[i]=B1c[j-17+i];
		BIP[i+17]=B2[i];
	}
	FIP[34]='\0';
	BIP[34]='\0';

	in=(char *)malloc(4096);
	memset(in,'\0',4096);
	getcwd(in,4096);
	i=strlen(in);
	dir=(char *)malloc(i+22);
	memset(dir,'\0',i+22);
	strcpy(dir,in);
	free(in);
        dir[i]='/';

	random=rand();
	memset(value,'\0',20);
	sprintf(value,"%d",random);
	strcat(dir,value);

	strcpy(pattern,"PRIMER_PAIR_0_COMPL.+TH=(.+)");
        cflags=REG_EXTENDED;
        regcomp(&reg,pattern,cflags);

	in=(char *)malloc(i+20);
	memset(in,'\0',i+20);
	strcpy(in,dir);
	strcat(in,"-par.txt");
	fp=fopen(in,"w");
        if(fp==NULL)
        {
                printf("Error: can't create the %s file!\n",in);
                exit(1);
        }
	
//generate
	fprintf(fp,"PRIMER_TASK=check_primers\nPRIMER_PICK_ANYWAY=1\nPRIMER_SALT_DIVALENT=4\nPRIMER_DNTP_CONC=1.4\nPRIMER_DNA_CONC=38\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n",primer3_config);
	for(i=0;i<5;i++)
	{
		for(j=i+1;j<6;j++)
			fprintf(fp,"SEQUENCE_ID=%d-%d\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",i,j,list[i],list[j]);
        }
	fprintf(fp,"SEQUENCE_ID=F3-FIP\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",list[0],FIP);
	fprintf(fp,"SEQUENCE_ID=F3-BIP\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",list[0],BIP);
	fprintf(fp,"SEQUENCE_ID=FIP-BIP\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",FIP,BIP);
	fprintf(fp,"SEQUENCE_ID=FIP-B3\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",FIP,list[5]);
        fprintf(fp,"SEQUENCE_ID=BIP-B3\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",BIP,list[5]);
	fclose(fp);

//run generage
	i=strlen(dir);
	out=(char *)malloc(i+20);
	memset(out,'\0',i+20);
        strcpy(out,dir);
        strcat(out,"-result.txt");

	i=strlen(primer3)+strlen(in)+strlen(out)+50;
	script=(char *)malloc(i);
	memset(script,'\0',i);
        sprintf(script,"%s -strict_tags %s > %s",primer3,in,out);
	system(script);
        remove(in);

	success=1;
	fp=fopen(out,"r");
	if(fp==NULL)
	{
		printf("Can't open %s file!\n",out);
		exit(1);
	}
	memset(value,'\0',20);
	j=strlen(primer3)+200;
	line=(char *)malloc(j);
	memset(line,'\0',j);
	while(fgets(line,j,fp)!=NULL)
        {
                if(regexec(&reg,line,2,pmatch,0)==0)  //begin
                {
                        take_regulate(pmatch,1,value,line);
                        TH=atof(value);
			if(flag==1&&TH>49)
			{
				success=0;
				break;
			}
			if(flag==0&&TH>44)
			{
				success=0;
				break;
			}
                }
	}
	fclose(fp);
	remove(out);
	free(in);
	free(out);
	free(script);
	regfree(&reg);
	free(line);
	return success;
}

int check_structure_loop(char F3[],char F2[],char F1c[],char B1c[],char B2[],char B3[],char LF[],char LB[],char *primer3,char *primer3_config,int flag)
{
        int i,j,random,cflags,success;
        float TH;
        char *in,*out,*dir,*list[10],pattern[50],value[20],FIP[35],BIP[35],*line,*script;
        FILE *fp;
        regex_t reg;
        regmatch_t pmatch[2];

//prepare
        list[0]=F3;
        list[1]=F2;
        list[3]=LF;
        list[4]=F1c;
        list[5]=B1c;
	list[6]=LB;
	list[8]=B2;
	list[9]=B3;
        j=strlen(F1c);
        for(i=0;i<17;i++)
        {
                FIP[i]=F1c[j-17+i];
                FIP[i+17]=F2[i];
        }
        j=strlen(B1c);
        for(i=0;i<17;i++)
        {
                BIP[i]=B1c[j-17+i];
                BIP[i+17]=B2[i];
        }
        FIP[34]='\0';
        BIP[34]='\0';
	list[2]=FIP;
	list[7]=BIP;

	in=(char *)malloc(4096);
        memset(in,'\0',4096);
        getcwd(in,4096);
	i=strlen(in);
	dir=(char *)malloc(i+22);
	memset(dir,'\0',i+22);
	strcpy(dir,in);
        dir[i]='/';
	free(in);

        random=rand();
        memset(value,'\0',20);
        sprintf(value,"%d",random);
        strcat(dir,value);

        strcpy(pattern,"PRIMER_PAIR_0_COMPL.+TH=(.+)");
        cflags=REG_EXTENDED;
        regcomp(&reg,pattern,cflags);

	i=strlen(dir);
	in=(char *)malloc(i+10);
        memset(in,'\0',i+10);
        strcpy(in,dir);
        strcat(in,"-par.txt");
        fp=fopen(in,"w");
        if(fp==NULL)
        {
                printf("Error: can't create the %s file!\n",in);
                exit(1);
        }
        
//generate
        fprintf(fp,"PRIMER_TASK=check_primers\nPRIMER_PICK_ANYWAY=1\nPRIMER_SALT_DIVALENT=4\nPRIMER_DNTP_CONC=1.4\nPRIMER_DNA_CONC=38\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n",primer3_config);
	if(list[3]!=NULL)
	{
        	for(i=0;i<=2;i++)
                        fprintf(fp,"SEQUENCE_ID=%d-LF\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",i,list[i],list[3]);
		for(i=4;i<10;i++)
		{
			if(i==6)
				continue;
			fprintf(fp,"SEQUENCE_ID=%d-LF\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",i,list[3],list[i]);
		}
        }
	if(list[6]!=NULL)
	{
		for(i=0;i<=5;i++)
		{
			if(i==3)
				continue;
			fprintf(fp,"SEQUENCE_ID=%d-LB\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",i,list[i],list[6]);
		}
		for(i=7;i<10;i++)
			fprintf(fp,"SEQUENCE_ID=%d-LB\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",i,list[6],list[i]);
	}
	if(list[3]!=NULL&&list[6]!=NULL)
	        fprintf(fp,"SEQUENCE_ID=LF-LB\nSEQUENCE_PRIMER=%s\nSEQUENCE_PRIMER_REVCOMP=%s\n=\n",list[3],list[6]);
        fclose(fp);

//run generage
	i=strlen(dir);
	out=(char *)malloc(i+20);
        memset(out,'\0',i+20);
        strcpy(out,dir);
        strcat(out,"-result.txt");

	i=strlen(primer3)+strlen(in)+strlen(out)+50;
	script=(char *)malloc(i);
        memset(script,'\0',i);
        sprintf(script,"%s -strict_tags %s > %s",primer3,in,out);
        system(script);
        remove(in);

        success=1;
        fp=fopen(out,"r");
        if(fp==NULL)
        {
                printf("Can't open %s file!\n",out);
                exit(1);
        }
	memset(value,'\0',20);
	j=strlen(primer3)+200;
	line=(char *)malloc(j);
	memset(line,'\0',j);
        while(fgets(line,j,fp)!=NULL)
        {
                if(regexec(&reg,line,2,pmatch,0)==0)  //begin
                {
                        take_regulate(pmatch,1,value,line);
                        TH=atof(value);
                        if(flag==1&&TH>49)
                        {
                                success=0;
                                break;
                        }
                        if(flag==0&&TH>44)
                        {
                                success=0;
                                break;
                        }
                }
        }
        fclose(fp);
        remove(out);
	free(in);
	free(out);
	free(script);
	regfree(&reg);
	free(line);
        return success;
}

int check_common_loop(struct Primer *F3,struct Primer *F2,struct Primer *LF,struct Primer *F1c,struct Primer *B1c,struct Primer *LB,struct Primer *B2,struct Primer *B3,int *result,int common_num)
{
        int dis,i,*temp,success;
        struct Node *c_F3,*c_F2,*c_LF,*c_F1c,*c_B1c,*c_LB,*c_B2,*c_B3;

	temp=(int *)malloc(common_num*sizeof(int));
	for(i=0;i<common_num;i++)
		temp[i]=0;
//plus
        for(c_F3=F3->common;c_F3;c_F3=c_F3->next)
        {
                if(c_F3->plus!=1)
                        continue;
                if(result[c_F3->gi]==0)
                        continue;
		if(temp[c_F3->gi])
			continue;
                for(c_F2=F2->common;c_F2;c_F2=c_F2->next)
                {
                        if(c_F2->gi!=c_F3->gi)
                                continue;
                        if(c_F2->plus!=1)
                                continue;
                        for(c_F1c=F1c->common;c_F1c;c_F1c=c_F1c->next)
                        {
                                if(c_F1c->gi!=c_F3->gi)
                                        continue;
                                if(c_F1c->minus!=1)
                                        continue;
                                for(c_B1c=B1c->common;c_B1c;c_B1c=c_B1c->next)
                                {
                                        if(c_B1c->gi!=c_F3->gi)
                                                continue;
                                        if(c_B1c->plus!=1)
                                                continue;
                                        for(c_B2=B2->common;c_B2;c_B2=c_B2->next)
                                        {
                                                if(c_B2->gi!=c_F3->gi)
                                                        continue;
                                                if(c_B2->minus!=1)
                                                        continue;
                                                for(c_B3=B3->common;c_B3;c_B3=c_B3->next)
                                                {
                                                        if(c_B3->gi!=c_F3->gi)
                                                                continue;
                                                        if(c_B3->minus!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=c_F2->pos-(c_F3->pos+F3->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=c_F1c->pos-c_F2->pos-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=c_B1c->pos-(c_F1c->pos+F1c->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=(c_B2->pos+B2->len-1)-(c_B1c->pos+B1c->len-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=c_B2->pos+B2->len-1-c_F2->pos-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=c_B3->pos-(c_B2->pos+B2->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
						//LF
							if(LF)
							{
								success=0;
								for(c_LF=LF->common;c_LF;c_LF=c_LF->next)
								{
									if(c_LF->gi!=c_F3->gi)
										continue;
									if(c_LF->minus!=1)
										continue;
									if(c_F2->pos+F2->len>c_LF->pos)
										continue;
									if(c_LF->pos+LF->len>c_F1c->pos)
										continue;
									success=1;
									break;
								}
								if(success==0)
									continue;
							}
						//LB
							if(LB)
                                                        {
                                                                success=0;
                                                                for(c_LB=LB->common;c_LB;c_LB=c_LB->next)
                                                                {
                                                                        if(c_LB->gi!=c_F3->gi)
                                                                                continue;
                                                                        if(c_LB->plus!=1)
                                                                                continue;
                                                                        if(c_B1c->pos+B1c->len>c_LB->pos)
                                                                                continue;
                                                                        if(c_LB->pos+LB->len>c_B2->pos)
                                                                                continue;
                                                                        success=1;
                                                                        break;
                                                                }
                                                                if(success==0)
                                                                        continue;
                                                        }
							temp[c_F3->gi]=1;
                                                }
                                        }
                                }
                        }
                }
        }
//minus
        for(c_F3=F3->common;c_F3;c_F3=c_F3->next)
        {
                if(c_F3->minus!=1)
                        continue;
                if(result[c_F3->gi]==0)
                        continue;  
		if(temp[c_F3->gi])
			continue;
                for(c_F2=F2->common;c_F2;c_F2=c_F2->next)
                {
                        if(c_F2->gi!=c_F3->gi)
                                continue;
                        if(c_F2->minus!=1)
                                continue;
                        for(c_F1c=F1c->common;c_F1c;c_F1c=c_F1c->next)
                        {
                                if(c_F1c->gi!=c_F3->gi)
                                        continue;
                                if(c_F1c->plus!=1)
                                        continue;
                                for(c_B1c=B1c->common;c_B1c;c_B1c=c_B1c->next)
                                {
                                        if(c_B1c->gi!=c_F3->gi)
                                                continue;
                                        if(c_B1c->minus!=1)
                                                continue;
                                        for(c_B2=B2->common;c_B2;c_B2=c_B2->next)
                                        {
                                                if(c_B2->gi!=c_F3->gi)
                                                        continue;
                                                if(c_B2->plus!=1)
                                                        continue;
                                                for(c_B3=B3->common;c_B3;c_B3=c_B3->next)
                                                {
                                                        if(c_B3->gi!=c_F3->gi)
                                                                continue;
                                                        if(c_B3->plus!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=c_F3->pos-(c_F2->pos+F2->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=(c_F2->pos+F2->len-1)-(c_F1c->pos+F1c->len-1)-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=c_F1c->pos-(c_B1c->pos+B1c->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=c_B1c->pos-c_B2->pos-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=c_F2->pos+F2->len-1-c_B2->pos-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=c_B2->pos-(c_B3->pos+B3->len-1)-1;
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
						//LF
                                                        if(LF)
                                                        {
                                                                success=0;
                                                                for(c_LF=LF->common;c_LF;c_LF=c_LF->next)
                                                                {
                                                                        if(c_LF->gi!=c_F3->gi)
                                                                                continue;
                                                                        if(c_LF->plus!=1)
                                                                                continue;
                                                                        if(c_F1c->pos+F1c->len>c_LF->pos)
                                                                                continue;
                                                                        if(c_LF->pos+LF->len>c_F2->pos)
                                                                                continue;
                                                                        success=1;
                                                                        break;
                                                                }
                                                                if(success==0)
                                                                        continue;
                                                        }
                                                //LB
                                                        if(LB)
                                                        {
                                                                success=0;
                                                                for(c_LB=LB->common;c_LB;c_LB=c_LB->next)
                                                                {
                                                                        if(c_LB->gi!=c_F3->gi)
                                                                                continue;
                                                                        if(c_LB->minus!=1)
                                                                                continue;
                                                                        if(c_B2->pos+B2->len>c_LB->pos)
                                                                                continue;
                                                                        if(c_LB->pos+LB->len>c_B1c->pos)
                                                                                continue;
                                                                        success=1;
                                                                        break;
                                                                }
                                                                if(success==0)
                                                                        continue;
                                                        }
                                                        temp[c_F3->gi]=1;
                                                }
                                        }
                                }
                        }
                }
        }
        for(i=0;i<common_num;i++)
        {
                if(result[i]&&temp[i]==0)
			return 0;
        }
        return 1;
}

int design_loop(struct Primer *p_F3,struct Primer *p_F2,struct Primer *p_LF,struct Primer *p_F1c,struct Primer *p_B1c,struct Primer *p_LB,struct Primer *p_B2,struct Primer *p_B3,struct Primer *result_Loop[], int *result,char *primer3,char *primer3_config,int flag[],int common_num,char *seq)
{
	int success;
	struct Primer *LF,*LB;
	struct Node *c_LF,*c_LB,*s_LF,*s_LB;
	char primer_F3[26],primer_F2[26],primer_F1c[26],primer_B1c[26],primer_B2[26],primer_B3[26],primer_LF[26],primer_LB[26];

	if(flag[7])
	{
		generate_primer(seq,primer_F3,p_F3->pos,p_F3->len,0);
		generate_primer(seq,primer_F2,p_F2->pos,p_F2->len,0);
		generate_primer(seq,primer_F1c,p_F1c->pos,p_F1c->len,1);
		generate_primer(seq,primer_B1c,p_B1c->pos,p_B1c->len,0);
		generate_primer(seq,primer_B2,p_B2->pos,p_B2->len,1);
		generate_primer(seq,primer_B3,p_B3->pos,p_B3->len,1);
	}
//LF and LB 
	success=0;
	LF=p_LF;
	while(LF)
	{
		if(LF->pos+LF->len>p_F1c->pos)
			break;
		LB=p_LB;
		if(flag[7])
			generate_primer(seq,primer_LF,LF->pos,LF->len,1);
		while(LB)
		{
			if(LB->pos+LB->len>p_B2->pos)
				break;
		//check_common
			if(flag[5])
			{
				success=check_common_loop(p_F3,p_F2,LF,p_F1c,p_B1c,LB,p_B2,p_B3,result,common_num);
				if(success==0)
				{
					LB=LB->next;
					continue;
				}
			}
		//check_structure
			if(flag[7])
			{
				generate_primer(seq,primer_LB,LB->pos,LB->len,0);
				success=check_structure_loop(primer_F3,primer_F2,primer_F1c,primer_B1c,primer_B2,primer_B3,primer_LF,primer_LB,primer3,primer3_config,flag[8]);
				if(success==0)
                        	{
                        	        LB=LB->next;
                        	        continue;
                        	}
			}
			result_Loop[0]=LF;
			result_Loop[1]=LB;
			success=1;
			break;
		}
		if(success==1)
			break;
		else
			LF=LF->next;
	}
	if(success==1)
		return success;
//only LF
	LF=p_LF;
	result_Loop[1]=NULL;
        while(LF)
        {
                if(LF->pos+LF->len>p_F1c->pos)
                        break;
	//check_common
		if(flag[5])
		{
			success=check_common_loop(p_F3,p_F2,LF,p_F1c,p_B1c,NULL,p_B2,p_B3,result,common_num);
			if(success==0)
			{
				LF=LF->next;
				continue;
			}
		}
	//check_structure
		if(flag[7])
		{
			generate_primer(seq,primer_LF,LF->pos,LF->len,1);
			success=check_structure_loop(primer_F3,primer_F2,primer_F1c,primer_B1c,primer_B2,primer_B3,primer_LF,NULL,primer3,primer3_config,flag[8]);
			if(success==0)
			{
				LF=LF->next;
				continue;
			}
		}
		result_Loop[0]=LF;
		success=1;
		break;
        }
        if(success==1)
                return success;
//only LB
	LB=p_LB;
	result_Loop[0]=NULL;
        while(LB)
        {
                if(LB->pos+LB->len>p_B2->pos)
                        break;
	//check_common
		if(flag[5])
		{
			success=check_common_loop(p_F3,p_F2,NULL,p_F1c,p_B1c,LB,p_B2,p_B3,result,common_num);
			if(success==0)
			{
				LB=LB->next;
				continue;
			}
		}
	//check_structure
		if(flag[7])
		{
			generate_primer(seq,primer_LB,LB->pos,LB->len,0);
			success=check_structure_loop(primer_F3,primer_F2,primer_F1c,primer_B1c,primer_B2,primer_B3,NULL,primer_LB,primer3,primer3_config,flag[8]);
			if(success==0)
			{
				LB=LB->next;
				continue;
			}
		}
		result_Loop[1]=LB;
		success=1;
		break;
        }
	return success;
}			
		
int check_add(int F3_pos,int *best_par,int expect)
{
        int i,dis;

	for(i=0;i<expect;i++)
	{
		if(best_par[i*16+1]==0) //the empty record
			return 1;
		dis=best_par[i*16]-F3_pos;
		if(abs(dis)<300)
			return 0;
	}
        return 1;
}

int check_gc(char *seq,int start,int end,int flag)//p_F3->pos,(p_B3->pos+p_B3->len),flag[8])
{
	int i,total=0;
	float gc;

	for(i=start;i<end;i++)
	{
		if(seq[i]=='C'||seq[i]=='G'||seq[i]=='c'||seq[i]=='g')
			total++;
	}
	gc=total*100.0/(end-start);
	if(flag==1&&gc>=45)
		return 1;
	if(flag==0&&gc<=45)
		return 1;
	return 0;
}
	
//function: add new one and return the pos
void add(struct Primer *F3,struct Primer *F2,struct Primer *F1c,struct Primer *B1c,struct Primer *B2,struct Primer *B3,struct Primer *LF,struct Primer *LB,int *best_par,int turn)
{

        best_par[turn*16]=F3->pos;
        best_par[turn*16+1]=F3->len;
        best_par[turn*16+2]=F2->pos;
        best_par[turn*16+3]=F2->len;
        best_par[turn*16+4]=F1c->pos;
        best_par[turn*16+5]=F1c->len;
        best_par[turn*16+6]=B1c->pos;
        best_par[turn*16+7]=B1c->len;
        best_par[turn*16+8]=B2->pos;
        best_par[turn*16+9]=B2->len;
        best_par[turn*16+10]=B3->pos;
        best_par[turn*16+11]=B3->len;
	if(LF!=NULL)
	{
		best_par[turn*16+12]=LF->pos;
		best_par[turn*16+13]=LF->len;
	}
	if(LB!=NULL)
	{
		best_par[turn*16+14]=LB->pos;
		best_par[turn*16+15]=LB->len;
	}
}

///how many GIs this primer can be used in, return the biggest common number
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
			list[i]=0;
		while(p_node)
		{
			list[p_node->gi]=1;
			p_node=p_node->next;
		}
		num=0;
                for(i=0;i<common;i++)
                        num=num+list[i];
		p_primer->total=num;
		p_primer=p_primer->next;
	}
	free(list);
}

//check this LAMP primers are uniq or not
//return=0: stop and return=1: go on
int check_uniq(struct Primer *F3,struct Primer *F2,struct Primer *F1c,struct Primer *B1c,struct Primer *B2,struct Primer *B3)
{
        struct Node *s_F3,*s_F2,*s_F1c,*s_B1c,*s_B2,*s_B3;

//plus
        for(s_F3=F3->special;s_F3;s_F3=s_F3->next)
        {
                if(s_F3->plus!=1)
                        continue;
                for(s_F2=F2->special;s_F2;s_F2=s_F2->next)
                {
                        if(s_F2->gi!=s_F3->gi)
                                continue;
                        if(s_F2->plus!=1)
                                continue;
                        for(s_F1c=F1c->special;s_F1c;s_F1c=s_F1c->next)
                        {
                                if(s_F1c->gi!=s_F3->gi)
                                        continue;
                                if(s_F1c->minus!=1)
                                        continue;
                                for(s_B1c=B1c->special;s_B1c;s_B1c=s_B1c->next)
                                {
                                        if(s_B1c->gi!=s_F3->gi)
                                                continue;
                                        if(s_B1c->plus!=1)
                                                continue;
                                        for(s_B2=B2->special;s_B2;s_B2=s_B2->next)
                                        {
                                                if(s_B2->gi!=s_F3->gi)
                                                        continue;
                                                if(s_B2->minus!=1)
                                                        continue;
                                                for(s_B3=B3->special;s_B3;s_B3=s_B3->next)
                                                {
                                                        if(s_B3->gi!=s_F3->gi)
                                                                continue;
                                                        if(s_B3->minus!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(s_F2->pos<s_F3->pos)
                                                                continue;
                                                //F2-F1c
                                                        if(s_F1c->pos<s_F2->pos+F2->len)
                                                                continue;
                                                //F1c-B1c
                                                        if(s_B1c->pos<s_F1c->pos+F1c->len)
                                                                continue;
                                                //B1c-B2
                                                        if(s_B2->pos<s_B1c->pos+B1c->len)
                                                                continue;
                                                //B2-B3
                                                        if(s_B3->pos<s_B2->pos)
                                                                continue;
						//whole
							if(s_B3->pos-s_F3->pos>1000)
								continue;
							return 0;
                                                }
                                        }
                                }
                        }
                }
        }
//minus
        for(s_F3=F3->special;s_F3;s_F3=s_F3->next)
        {
                if(s_F3->minus!=1)
                        continue;
                for(s_F2=F2->special;s_F2;s_F2=s_F2->next)
                {
                        if(s_F2->gi!=s_F3->gi)
                                continue;
                        if(s_F2->minus!=1)
                                continue;
                        for(s_F1c=F1c->special;s_F1c;s_F1c=s_F1c->next)
                        {
                                if(s_F1c->gi!=s_F3->gi)
                                        continue;
                                if(s_F1c->plus!=1)
                                        continue;
                                for(s_B1c=B1c->special;s_B1c;s_B1c=s_B1c->next)
                                {
                                        if(s_B1c->gi!=s_F3->gi)
                                                continue;
                                        if(s_B1c->minus!=1)
                                                continue;
                                        for(s_B2=B2->special;s_B2;s_B2=s_B2->next)
                                        {
                                                if(s_B2->gi!=s_F3->gi)
                                                        continue;
                                                if(s_B2->plus!=1)
                                                        continue;
                                                for(s_B3=B3->special;s_B3;s_B3=s_B3->next)
                                                {
                                                        if(s_B3->gi!=s_F3->gi)
                                                                continue;
                                                        if(s_B3->plus!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(s_F3->pos<s_F2->pos)
                                                                continue;
                                                //F2-F1c
                                                        if(s_F2->pos<s_F1c->pos+F1c->len)
                                                                continue;
                                                //F1c-B1c
                                                        if(s_F1c->pos<s_B1c->pos+B1c->len)
                                                                continue;
                                                //B1c-B2
                                                        if(s_B1c->pos<s_B2->pos+B2->len)
                                                                continue;
                                                //B2-B3
                                                        if(s_B2->pos<s_B3->pos)
                                                                continue;
						//whole
							if(s_F3->pos-s_B3->pos>1000)
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
