//gpu version
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<regex.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<cuda_runtime.h>
#include<cuda.h>

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

void primer(char *seq,char primer[],int start,int length,int flag)
{
        char read[270];
        int i;
        
        for(i=start;i<start+length;i++)
                read[i-start]=seq[i];


        if(flag==0)
        {
                for(i=0;i<length;i++)
                {
                        if(read[i]=='A'||read[i]=='a')
                                primer[i]='A';
                        else if(read[i]=='T'||read[i]=='t')
                                primer[i]='T';
                        else if(read[i]=='C'||read[i]=='c')
                                primer[i]='C';    
                        else if(read[i]=='G'||read[i]=='g')
                                primer[i]='G';
                        else
                                primer[i]='N';
                }
                primer[i]='\0';
        }
        else
        {
                for(i=0;i<length;i++)
                {
                        if(read[length-1-i]=='A'||read[length-1-i]=='a')
                                primer[i]='T';
                        else if(read[length-1-i]=='T'||read[length-1-i]=='t')
                                primer[i]='A';
                        else if(read[length-1-i]=='C'||read[length-1-i]=='c')
                                primer[i]='G';
                        else if(read[length-1-i]=='G'||read[length-1-i]=='g')
                                primer[i]='C';
                        else 
                                primer[i]='N';
                }
                primer[length]='\0';
        }
}

int secondary_check(int expect,int have,int circle,int *par,char *file,FILE *LAMP,int *h_result,int *apply,char *seq,int common,int *list)
{
	FILE *fp;       
        char line[1000],result[20],pattern1[100],pattern2[100],pattern3[50],F3[30],F2[30],F1c[30],B1c[30],B2[30],B3[30];
        int flag,cflags,value[3],count,i,j;
        regex_t reg1,reg2,reg3;
        regmatch_t pmatch[4];
        float TH,max_TH;
        
        strcpy(pattern1,"SEQUENCE.+\\=(\\w+)\\-(\\w+)\\-(\\w+)");
        cflags=REG_EXTENDED;
        regcomp(&reg1,pattern1,cflags);
        strcpy(pattern2,"PRIMER\\_LEFT\\_0\\_TM\\=(.+)$");
        regcomp(&reg2,pattern2,cflags);
        strcpy(pattern3,"PRIMER.PAIR.+COMPL.+\\=(.+)$");
        regcomp(&reg3,pattern3,cflags);

        fp=fopen(file,"r");
        if(fp==NULL)
        {
                printf("can't open temp-result.txt file\n");
                exit(1);
        }
	count=have;
        while(fgets(line,1000,fp)!=NULL)
        {
                if(regexec(&reg1,line,4,pmatch,0)==0)  //begin
                {
                        take_regulate(pmatch,3,result,line);
                        value[2]=atoi(result); //Primer_turn in LAMP
                        if(value[2]==1) //a new LAMP primers
                        {
                                flag=1;
                        //check F3 pos
                                take_regulate(pmatch,2,result,line);
                                value[1]=atoi(result); //F3_pos
                                flag=check_add(value[1],par,count);
				if(flag==0)
					continue;
                        }

                        if(value[2]==45&&flag)
                        {
                                take_regulate(pmatch,1,result,line);
                                value[0]=atoi(result);
                        }
                        continue;
                }

                if(regexec(&reg2,line,2,pmatch,0)==0&&value[2]==1&&flag) //the max TH
                {
                        take_regulate(pmatch,1,result,line);
                        max_TH=atof(result);
                        continue;
                }
                if(regexec(&reg3,line,2,pmatch,0)==0)
                {
                        take_regulate(pmatch,1,result,line);
                        TH=atof(result);
                        if(TH>max_TH-10)
                                flag=0;
                        continue;                       
                }
                if(line[0]=='='&&value[2]==45&&flag==1)
                {
		//add
			primer(seq,F3,h_result[12*value[0]],h_result[12*value[0]+1],0);
			primer(seq,F2,h_result[12*value[0]+2],h_result[12*value[0]+3],0);
			primer(seq,F1c,h_result[12*value[0]+4],h_result[12*value[0]+5],1);
                        primer(seq,B1c,h_result[12*value[0]+6],h_result[12*value[0]+7],0);
                        primer(seq,B2,h_result[12*value[0]+8],h_result[12*value[0]+9],1);
                        primer(seq,B3,h_result[12*value[0]+10],h_result[12*value[0]+11],1);

			fprintf(LAMP,"The %d LAMP primer can be used in %d genomes:\n",(count+1),circle);
			fprintf(LAMP,"this LAMP primer can be used in:");
			i=0;
			for(j=0;j<common;j++)
			{
				if(apply[common*value[0]+j]==0)
					continue;
				i++;
				if(i==circle)
					fprintf(LAMP,"%d\n",list[j]);
				else
					fprintf(LAMP,"%d,",list[j]);
			}
			if(circle==common)
				fprintf(LAMP,"this LAMP primer can't be use in:None\n");
			else
			{
				i=0;
				for(j=0;j<common;j++)
				{
					if(apply[common*value[0]+j]==1)
						continue;
					i++;
					if(i==(common-circle))
						fprintf(LAMP,"%d\n",list[j]);
					else
						fprintf(LAMP,"%d,",list[j]);
				}
			}
			fprintf(LAMP,"  F3 start pos:%d, length:%d bp, sequence:%s\n",h_result[12*value[0]],h_result[12*value[0]+1],F3);
			fprintf(LAMP,"  F2 start pos:%d, length:%d bp, sequence:%s\n",h_result[12*value[0]+2],h_result[12*value[0]+3],F2);
			fprintf(LAMP,"  F1c start pos:%d, length:%d bp, sequence:%s\n",h_result[12*value[0]+4],h_result[12*value[0]+5],F1c);
			fprintf(LAMP," FIP sequence:%s-%s\n",F1c,F2);
			fprintf(LAMP," BIP sequence:%s-%s\n",B1c,B2);
			fprintf(LAMP,"  B1c start pos:%d, length:%d bp, sequence:%s\n",h_result[12*value[0]+6],h_result[12*value[0]+7],B1c);
			fprintf(LAMP,"  B2 start pos:%d, length:%d bp, sequence:%s\n",h_result[12*value[0]+8],h_result[12*value[0]+9],B2);
			fprintf(LAMP,"  B3 start pos:%d, length:%d bp, sequence:%s\n\n",h_result[12*value[0]+10],h_result[12*value[0]+11],B3);

			par[count]=h_result[12*value[0]];
			count++;
			if(count==expect)
				return count;
                }
        }
        fclose(fp);
        regfree(&reg1);
        regfree(&reg2);
        regfree(&reg3);
	return count;
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
struct Primer *read_par(char path[],int common)
{
        char in[2000];
        int pos,len,gi,position,plus,minus,size,i;
        struct Primer *new_primer,*p_primer,*head;
        struct Node *new_node,*p_node;
        FILE *fp;

///read the  primer file
        memset(in,'\0',2000*sizeof(char));       
        strcpy(in,path);
        strcat(in,".txt");  //suffix of primer candidate file
        fp=fopen(in,"r");
        if(fp==NULL)
        {
                printf("Error: can't open the %s file!\n",in);
                exit(1);
        }
        
        size=sizeof(struct Primer);
        i=0;
        while(fscanf(fp,"pos:%d\tlength:%d\t+:%d\t-:%d\n",&pos,&len,&plus,&minus)!=EOF)
        {
                new_primer=(struct Primer *)malloc(size);
                new_primer->pos=pos;
                new_primer->len=len;
                new_primer->total_common=0;
		new_primer->total_special=0;
		new_primer->total=0;
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
                printf("Sorry! Don't have candidate primer!\n");
                exit(1);
        }

//parameter of common
        memset(in,'\0',2000*sizeof(char));
        strcpy(in,path);
        strcat(in,"-common.txt"); //suffix of parameter
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
        //check the common
                if(gi>=common)
                {
                        gi++;
                        printf("ERROR!\n The primer(pos is %d and length is %d) from %s can be used in more than %d GIs, but the common parameter is %d! You should check the common parameter!\n",pos,len,in,gi,common);
                        exit(1);
                } 
                new_node=(struct Node *)malloc(size);
                new_node->pos=position;
                new_node->gi=gi;
                new_node->plus=plus;
                new_node->minus=minus;

        //find the primer
                while(p_primer->pos!=pos||p_primer->len!=len)
                {
                        if(p_primer->next==NULL)
                        {
                                p_primer=head;
                        }
                        else
                        {
                                p_primer=p_primer->next;
                        }
                }
                p_node=p_primer->common;
                p_primer->common=new_node;
                new_node->next=p_node;
        }
        fclose(fp);

//paramter for special
        memset(in,'\0',2000*sizeof(char));
        strcpy(in,path);
        strcat(in,"-special.txt"); //suffix of parameter
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
                while(p_primer->pos!=pos||p_primer->len!=len)
                {
                        if(p_primer->next==NULL)
                        {
                                p_primer=head;
                        }
                        else
                        {
                                p_primer=p_primer->next;
                        }
                }
                
                p_node=p_primer->special;
                p_primer->special=new_node;
                new_node->next=p_node;
        }
        fclose(fp);
        return head;
}

//function: check how many GIs this primer can be used for
__device__ int check_common(int *primerL,int *primerS,int *commonL,int *commonS,int *start_commonL,int *start_commonS,int *end_commonL,int *end_commonS,int turn[],int common,int *apply)
{
        int pos[6],i,dis,num,start[6],end[6],gi;

        for(i=0;i<common;i++)
        {
                apply[common*turn[0]+i]=0;
        }
//plus
	start[0]=start_commonS[turn[0]];
	end[0]=end_commonS[turn[0]];
	start[1]=start_commonS[turn[1]];
        end[1]=end_commonS[turn[1]];
	start[2]=start_commonL[turn[2]];
        end[2]=end_commonL[turn[2]];
        start[3]=start_commonL[turn[3]];
        end[3]=end_commonL[turn[3]];
	start[4]=start_commonS[turn[4]];
        end[4]=end_commonS[turn[4]];
        start[5]=start_commonS[turn[5]];
        end[5]=end_commonS[turn[5]];
        for(pos[0]=start[0];pos[0]<end[0];pos[0]++)
        {
                if(commonS[4*pos[0]+2]!=1)
                        continue;
		gi=commonS[4*pos[0]];
                if(apply[common*turn[0]+gi]==1)
                        continue;
                for(pos[1]=start[1];pos[1]<end[1];pos[1]++)
                {
                        if(commonS[4*pos[1]]!=gi)
                                continue;
                        if(commonS[4*pos[1]+2]!=1)
                                continue;
                        for(pos[2]=start[2];pos[2]<end[2];pos[2]++)
                        {
                                if(commonL[4*pos[2]]!=gi)
                                        continue;
                                if(commonL[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=start[3];pos[3]<end[3];pos[3]++)
                                {
                                        if(commonL[4*pos[3]]!=gi)
                                                continue;
                                        if(commonL[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=start[4];pos[4]<end[4];pos[4]++)
                                        {
                                                if(commonS[4*pos[4]]!=gi)
                                                        continue;
                                                if(commonS[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=start[5];pos[5]<end[5];pos[5]++)
                                                {
                                                        if(commonS[4*pos[5]]!=gi)
                                                                continue;
                                                        if(commonS[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=commonS[4*pos[1]+1]-(commonS[4*pos[0]+1]+primerS[4*turn[0]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=commonL[4*pos[2]+1]-commonS[4*pos[1]+1];
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=commonL[4*pos[3]+1]-(commonL[4*pos[2]+1]+primerL[4*turn[2]+1]);
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=(commonS[4*pos[4]+1]+primerS[4*turn[4]+1])-(commonL[4*pos[3]+1]+primerL[4*turn[3]+1]);
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=commonS[4*pos[4]+1]+primerS[4*turn[4]+1]-1-commonS[4*pos[1]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=commonS[4*pos[5]+1]-(commonS[4*pos[4]+1]+primerS[4*turn[4]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                        apply[common*turn[0]+gi]=1;
                                                }
                                        }
                                }
                        }
                }
        }
//minus
        for(pos[0]=start[0];pos[0]<end[0];pos[0]++)
        {
                if(commonS[4*pos[0]+3]!=1)
                        continue;
		gi=commonS[4*pos[0]];
                if(apply[turn[0]*common+gi]==1)
                        continue;  //this GI can common

                for(pos[1]=start[1];pos[1]<end[1];pos[1]++)
                {
                        if(commonS[4*pos[1]]!=gi)
                                continue;
                        if(commonS[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=start[2];pos[2]<end[2];pos[2]++)
                        {
                                if(commonL[4*pos[2]]!=gi)
                                        continue;
                                if(commonL[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=start[3];pos[3]<end[3];pos[3]++)
                                {
                                        if(commonL[4*pos[3]]!=gi)
                                                continue;
                                        if(commonL[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=start[4];pos[4]<end[4];pos[4]++)
                                        {
                                                if(commonS[4*pos[4]]!=gi)
                                                        continue;
                                                if(commonS[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=start[5];pos[5]<end[5];pos[5]++)
                                                {
                                                        if(commonS[4*pos[5]]!=gi)
                                                                continue;
                                                        if(commonS[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        dis=commonS[4*pos[0]+1]-(commonS[4*pos[1]+1]+primerS[4*turn[1]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                //F2-F1c
                                                        dis=(commonS[4*pos[1]+1]+primerS[4*turn[1]+1])-(commonL[4*pos[2]+1]+primerL[4*turn[2]+1]);
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F1c-B1c
                                                        dis=commonL[4*pos[2]+1]-(commonL[4*pos[3]+1]+primerL[4*turn[3]+1]);
                                                        if(dis<0)
                                                                continue;
                                                //B1c-B2
                                                        dis=commonL[4*pos[3]+1]-commonS[4*pos[4]+1]-1;
                                                        if(dis<40)
                                                                continue;
                                                        if(dis>60)
                                                                continue;
                                                //F2-B2
                                                        dis=commonS[4*pos[1]+1]+primerS[4*turn[1]+1]-1-commonS[4*pos[4]+1]-1;
                                                        if(dis<120)
                                                                continue;
                                                        if(dis>180)
                                                                continue;
                                                //B2-B3
                                                        dis=commonS[4*pos[4]+1]-(commonS[4*pos[5]+1]+primerS[4*turn[5]+1]);
                                                        if(dis<0)
                                                                continue;
                                                        if(dis>20)
                                                                continue;
                                                        apply[common*turn[0]+gi]=1;
                                                }
                                        }
                                }
                        }
                }
        }
        num=0;
        for(i=0;i<common;i++)
        {
                num=num+apply[common*turn[0]+i];
        }
        return num;
}

//check this LAMP primers are uniq or not
//return=0: stop and return=1: go on
__device__ int check_uniq(int *primerL,int *primerS,int *specialL,int *specialS,int *start_specialL,int *start_specialS,int *end_specialL,int *end_specialS,int turn[])
{
        int pos[6],start[6],end[6],gi;

	start[0]=start_specialS[turn[0]];
	end[0]=end_specialS[turn[0]];
	start[1]=start_specialS[turn[1]];
	end[1]=end_specialS[turn[1]];
	start[2]=start_specialL[turn[2]];
        end[2]=end_specialL[turn[2]];
        start[3]=start_specialL[turn[3]];
        end[3]=end_specialL[turn[3]];
	start[4]=start_specialS[turn[4]];
        end[4]=end_specialS[turn[4]];
        start[5]=start_specialS[turn[5]];
        end[5]=end_specialS[turn[5]];
//plus
        for(pos[0]=start[0];pos[0]<end[0];pos[0]++)
        {
                if(specialS[4*pos[0]+2]!=1)
                        continue;
		gi=specialS[4*pos[0]];
                for(pos[1]=start[1];pos[1]<end[1];pos[1]++)
                {
			if(specialS[4*pos[1]]!=gi)
                                continue;
                        if(specialS[4*pos[1]+2]!=1)
				continue;
                        for(pos[2]=start[2];pos[2]<end[2];pos[2]++) //F1c
                        {
                                if(specialL[4*pos[2]]!=gi)
                                        continue;
                                if(specialL[4*pos[2]+3]!=1)
                                        continue;
                                for(pos[3]=start[3];pos[3]<end[3];pos[3]++) //B1c
                                {
                                        if(specialL[pos[3]*4]!=gi)
                                                continue;
                                        if(specialL[4*pos[3]+2]!=1)
                                                continue;
                                        for(pos[4]=start[4];pos[4]<end[4];pos[4]++) //B2
                                        {
                                                if(specialS[4*pos[4]]!=gi)
                                                        continue;
                                                if(specialS[4*pos[4]+3]!=1)
                                                        continue;
                                                for(pos[5]=start[5];pos[5]<end[5];pos[5]++)
                                                {
                                                        if(specialS[4*pos[5]]!=gi)
                                                                continue;
                                                        if(specialS[4*pos[5]+3]!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(specialS[4*pos[1]+1]<specialS[4*pos[0]+1])
                                                                continue;
                                                //F2-F1c
                                                        if(specialL[4*pos[2]+1]<specialS[4*pos[1]+1]+primerS[4*turn[1]+1])
                                                                continue;
                                                //F1c-B1c
                                                        if(specialL[4*pos[3]+1]<specialL[4*pos[2]+1])
                                                                continue;
                                                //B1c-B2
                                                        if(specialS[4*pos[4]+1]<specialL[4*pos[3]+1]+primerL[4*turn[3]+1])
                                                                continue;
                                                //B2-B3
                                                        if(specialS[4*pos[5]+1]<specialS[4*pos[4]+1])
                                                                continue;
                                                //whole
                                                        if(specialS[4*pos[5]+1]-specialS[4*pos[0]+1]>1000)
                                                                continue;
                                                        return 0;
                                                }//B3
                                        }
                                }//B1c
                        }
                }//F2
        }

//minus
        for(pos[0]=start[0];pos[0]<end[0];pos[0]++)
        {
                if(specialS[4*pos[0]+3]!=1)
                        continue;
		gi=specialS[4*pos[0]];
                for(pos[1]=start[1];pos[1]<end[1];pos[1]++)
                {
                        if(specialS[4*pos[1]]!=gi)
                                continue;
                        if(specialS[4*pos[1]+3]!=1)
                                continue;
                        for(pos[2]=start[2];pos[2]<end[2];pos[2]++)
                        {
                                if(specialL[4*pos[2]]!=gi)
                                        continue;
                                if(specialL[4*pos[2]+2]!=1)
                                        continue;
                                for(pos[3]=start[3];pos[3]<end[3];pos[3]++)
                                {
                                        if(specialL[4*pos[3]]!=gi)
                                                continue;
                                        if(specialL[4*pos[3]+3]!=1)
                                                continue;
                                        for(pos[4]=start[4];pos[4]<end[4];pos[4]++)
                                        {
                                                if(specialS[4*pos[4]]!=gi)
                                                        continue;
                                                if(specialS[4*pos[4]+2]!=1)
                                                        continue;
                                                for(pos[5]=start[5];pos[5]<end[5];pos[5]++)
                                                {
                                                        if(specialS[4*pos[5]]!=gi)
                                                                continue;
                                                        if(specialS[4*pos[5]+2]!=1)
                                                                continue;
                                                //F3-F2 
                                                        if(specialS[4*pos[0]+1]<specialS[4*pos[1]+1])
                                                                continue;
                                                //F2-F1c
                                                        if(specialS[4*pos[1]+1]<specialL[4*pos[2]+1]+primerL[4*turn[2]+1])
                                                                continue;
                                                //F1c-B1c
                                                        if(specialL[4*pos[2]+1]<specialL[4*pos[3]+1])
                                                                continue;
                                                //B1c-B2
                                                        if(specialL[4*pos[3]+1]<specialS[4*pos[4]+1]+primerS[4*turn[4]+1])
                                                                continue;
                                                //B2-B3
                                                        if(specialS[4*pos[4]+1]<specialS[4*pos[5]+1])
                                                                continue;
                                                //whole
                                                        if(specialS[4*pos[0]+1]-specialS[4*pos[5]+1]>1000)
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

	if(id<num_first)
	{
		next[id]=-1;
		i=id;
		if(i>=num_second)
		{
			i=num_second-1;
		}
		if(second[4*i]>first[id*4]+first[id*4+1])
		{
			while((i>=0)&&(second[4*i]>first[id*4]+first[id*4+1]))
			{
				next[id]=i;
				i--;
			}
		}
		else
		{
			while((i<num_second)&&(second[i*4]<=first[id*4]+first[id*4+1]))
			{
				next[id]=i;
				i++;
			}
		}
	}
	__syncthreads();
}

__device__ float cal_GC(char *d_seq,int begin,int length)
{
	int i,total=0;
	float result;

	for(i=0;i<length;i++)
	{
		if(d_seq[begin+i]=='C'||d_seq[begin+i]=='c'||d_seq[begin+i]=='G'||d_seq[begin+i]=='g')
			total++;
	}
	result=(float)total/length;
	return result;
}
//caculate
__global__ void find_primer(int *primerL,int *primerS,int *specialL,int *specialS,int *commonL,int *commonS,int *start_specialL,int *start_specialS,int *start_commonL,int *start_commonS,int *end_specialL,int *end_specialS,int *end_commonL,int *end_commonS,int *nextS,int *nextL,int *other_toL,int *other_toS,int *result,int numL,int numS,int common,int *apply,int *best_par,int circle,int expect,char *d_seq,int begin) //apply is temp
{
	int id=blockDim.x*blockIdx.x+threadIdx.x;
	int turn[6],flag,i;
	float gc[2];

	if(id<numS)
	{
		result[12*id+1]=0;//not LAMP, as a flag
	//check overlap
		for(i=0;i<expect;i++)
		{
			if(best_par[i]==-1)
				break;
			flag=primerS[4*id]-best_par[i];
			if(flag<300&&flag>-300)
				return;
		}
	//combine
		turn[0]=id; //one thread, one F3
		for(turn[1]=nextS[turn[0]];turn[1]<numS;turn[1]++) //F2
		{
			if(turn[1]==-1)
				break;
			if(primerS[turn[1]*4]<primerS[turn[0]*4]+primerS[turn[0]*4+1])
				continue;
			if(primerS[4*turn[1]]-(primerS[4*turn[0]]+primerS[4*turn[0]+1])>20)
				break;
			for(turn[2]=other_toL[turn[1]];turn[2]<numL;turn[2]++) //F1c
			{
				if(turn[2]==-1)
					break;
				if(primerL[4*turn[2]]-primerS[4*turn[1]]-1<40)
					continue;
                                if(primerL[4*turn[2]]-primerS[4*turn[1]]-1>60)
                                	break;
                                for(turn[3]=nextL[turn[2]];turn[3]<numL;turn[3]++)   //B1c
                                {
                                        if(turn[3]==-1)
                                        	break;
                                        if(primerL[4*turn[3]]<primerL[4*turn[2]]+primerL[4*turn[2]+1])
                                        	continue;
                                        if(primerL[4*turn[3]]-primerL[4*turn[2]]>85)
                                        	break;
                                        for(turn[4]=other_toS[turn[3]];turn[4]<numS;turn[4]++)   //B2
                                        {
                                                if(turn[4]==-1)
                                                	break;
                                                if((primerS[4*turn[4]]+primerS[4*turn[4]+1]-1)-(primerL[4*turn[3]]+primerL[4*turn[3]+1])<40)
                                                	continue;
                                                if((primerS[4*turn[4]]+primerS[4*turn[4]+1]-1)-(primerL[4*turn[3]]+primerL[4*turn[3]+1])>60)
                                                	break;
                                                if(primerS[4*turn[4]]+primerS[4*turn[4]+1]-1-primerS[turn[1]*4]-1<120)
                                                	continue;
                                                if(primerS[4*turn[4]]+primerS[4*turn[4]+1]-1-primerS[turn[1]*4]-1>180)
                                                	break;
                                                for(turn[5]=nextS[turn[4]];turn[5]<numS;turn[5]++)  //B3
                                                {
                                                        if(turn[5]==-1)
                                                        	break;
                                                        if(primerS[turn[5]*4]<primerS[4*turn[4]]+primerS[4*turn[4]+1])
                                                        	continue;
                                                        if(primerS[turn[5]*4]-(primerS[4*turn[4]]+primerS[4*turn[4]+1])>20)
                                                        	break;
						//plus-minus
							flag=primerS[4*turn[0]+2]&&primerS[4*turn[1]+2]&&primerL[4*turn[2]+3]&&primerL[4*turn[3]+2]&&primerS[4*turn[4]+3]&&primerS[4*turn[5]+3];
                                                        if(flag==0)
                                                        	continue;
		
							gc[0]=cal_GC(d_seq,(primerS[4*turn[0]]-begin),primerS[4*turn[0]+1]);
							gc[1]=cal_GC(d_seq,(primerS[4*turn[0]]-begin),(primerS[4*turn[5]]+primerS[4*turn[5]+1]-primerS[4*turn[0]]));
							if(gc[1]<0.45&&gc[0]>0.5)
								continue;
							if(gc[1]>0.45&&gc[0]<0.5)
								continue;
							flag=check_uniq(primerL,primerS,specialL,specialS,start_specialL,start_specialS,end_specialL,end_specialS,turn);
							if(flag==0)
								continue;
							flag=check_common(primerL,primerS,commonL,commonS,start_commonL,start_commonS,end_commonL,end_commonS,turn,common,apply);

							if(flag<circle)
								continue;
							result[id*12]=primerS[4*turn[0]];
							result[id*12+1]=primerS[4*turn[0]+1];
							result[id*12+2]=primerS[4*turn[1]];
							result[id*12+3]=primerS[4*turn[1]+1];
							result[id*12+4]=primerL[4*turn[2]];
							result[id*12+5]=primerL[4*turn[2]+1];
							result[id*12+6]=primerL[4*turn[3]];
							result[id*12+7]=primerL[4*turn[3]+1];
							result[id*12+8]=primerS[4*turn[4]];
							result[id*12+9]=primerS[4*turn[4]+1];
							result[id*12+10]=primerS[4*turn[5]];
							result[id*12+11]=primerS[4*turn[5]+1];
							return;
						}
					}
				}
			}
		}
	}
	__syncthreads();
}

main(int argc,char **argv)
{
	int i,j,k,common,flag,expect,circle,have,*list,begin,stop;
	char out[2000],path_small[2000],path_large[2000],path_fa[2000],script[6000],directory[2000],file1[2000],file2[2000],path_common[2000],path_primer3[2000],check[2000];
	char *temp,*seq,F3[26],F2[26],F1c[26],B1c[26],B2[26],B3[26],*h_seq,*d_seq;
	FILE *fp,*LAMP;
	struct Primer *headL,*headS,*tempL,*tempS,*storeL,*storeS; 
	struct Node *p_node,*p_temp;
	char usage[200]="./a.out -small small_primers_path -large large_primer3_path -out out_file -fa fna_path -common common_GI_file -primer3 primer3_path -expect(optional,default=10) expect_output_number\n";
	time_t start,end;
	
	long int memory;
	cudaDeviceProp prop;
	int *h_primerL,*d_primerL,*h_primerS,*d_primerS; //primer info
	int *h_commonL,*d_commonL,*h_commonS,*d_commonS,*h_specialL,*d_specialL,*h_specialS,*d_specialS; //common and special info
	int *h_result,*d_result,*d_apply,*h_apply,*h_par,*d_par;
	int *d_other_to_L,*d_other_to_S,*d_nextL,*d_nextS;
	int *h_start_commonL,*d_start_commonL,*h_end_commonL,*d_end_commonL,*h_start_specialL,*d_start_specialL,*h_end_specialL,*d_end_specialL,*h_start_commonS,*d_start_commonS,*h_end_commonS,*d_end_commonS,*h_start_specialS,*d_start_specialS,*h_end_specialS,*d_end_specialS;
	int num[7],count[3],block,m,n;
	char *p_primer[6];
	
	expect=10; //default output max 10 LAMP primers
	start=time(NULL);
/////read the parameters
	if(argc<11)
	{
		printf("Error!\n%s\n",usage);
		exit(1);
	}

	j=0;
	for(i=1;i<argc-1;i=i+2)
	{
		if(strcmp(argv[i],"-small")==0)
		{
			strcpy(path_small,argv[i+1]);
			j++;
		}
                else if(strcmp(argv[i],"-large")==0)
                {
                        strcpy(path_large,argv[i+1]);
                        j++;
                }
                else if(strcmp(argv[i],"-out")==0)
                {
                        strcpy(out,argv[i+1]);
                        j++;
                }
		else if(strcmp(argv[i],"-fa")==0)
		{
			strcpy(path_fa,argv[i+1]);
			j++;
		}
		else if(strcmp(argv[i],"-common")==0)
		{
			strcpy(path_common,argv[i+1]);
			j++;
		}
		else if(strcmp(argv[i],"-primer3")==0)
		{
			strcpy(path_primer3,argv[i+1]);
			j++;
		}
		else if(strcmp(argv[i],"-expect")==0)
		{
			expect=atoi(argv[i+1]);
		}
		else
		{
			printf("Error!\n%s\n",usage);
			exit(1);
		}
	}

	if(j!=6)
	{
		printf("Error!\n%s\n",usage);
		exit(1);
	}
//the directory of program
	strcpy(directory,out);
	j=strlen(directory);
	j--;
	while(directory[j]!='/'&&j>=0)
	{
		directory[j]='\0';
		j--;
	}
	strcpy(file1,directory);
	strcat(file1,"temp-par.txt");
	strcpy(file2,directory);
	strcat(file2,"temp-result.txt");
//common_GI-list
	fp=fopen(path_common,"r");                 
        if(fp==NULL)
        {          
                printf("Error: can't open the %s file!\n",path_common);
                exit(1);                  
        }

	i=0;
        while(fscanf(fp,"%d\n",&j)!=EOF)
        {       
                i++;
        }
        list=(int *)malloc(i*sizeof(int));    
        rewind(fp);
        common=0;
        while(fscanf(fp,"%d\n",&j)!=EOF)
        {
                list[common]=j;
                common++;
        }
        fclose(fp);	

//primer3 path
	i=strlen(path_primer3);
	i--;
	if(path_primer3[i]=='/')
	{
		strcpy(check,path_primer3);
		strcat(check,"primer3_core");
		if(access(check,0)!=0)
		{
			printf("Can't fine primer3_core program in %s directory!\n",path_primer3);
			exit(1);
		}
	}
	else
	{
		for(j=0;j<12;j++)
			check[j]=path_primer3[i-j];
		check[j]='\0';
		if(strcmp(check,"eroc_3remirp")==0)
		{
			if(access(path_primer3,0)==0)
			{
				j=i;
				while(j>=0&&path_primer3[j]!='/')
				{
					path_primer3[j]='\0';
					j--;
				}
			}
			else
			{
				printf("Can't fine primer3_core program in %s directory!\n",path_primer3);
				exit(1);
			}
		}
		else
		{
			strcpy(check,path_primer3);
			strcat(check,"/primer3_core");
			if(access(check,0)==0)
			{
				path_primer3[i+1]='/';
				path_primer3[i+2]='\0';
			}
			else
			{
				printf("Can't fine primer3_core program in %s directory!\n",path_primer3);
				exit(1);
			}
		}
	}

	cudaGetDeviceProperties(&prop,0); //read parameters
	headS=read_par(path_small,common);
	headL=read_par(path_large,common);

//common statistics
	how_many(headL,common);
	how_many(headS,common);

//the genome sequence
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

        flag=0;
        j=0;
        for(i=0;temp[i]!='\0';i++)
        {
                if(temp[i]=='\n')
                {
                        flag+=1;
                        continue;
                }
                else
                {
                        if(flag!=0)
                        {
                                seq[j]=temp[i];
                                j++;
                        }
                }
        }
        free(temp);

	end=time(NULL);
	printf("The prepare time is %0.1f seconds!\n",difftime(end,start));
	start=time(NULL);

//LAMP-GPU
	LAMP=fopen(out,"w");
        if(LAMP==NULL)
        {
                printf("Error: can't create the %s file!\n",out);
                exit(1);
        }
	have=0;
	prop.maxThreadsPerBlock=prop.maxThreadsPerBlock/4;
	k=prop.maxThreadsPerBlock*60000; //max primers, max threads
	h_par=(int *)malloc(expect*sizeof(int));
	for(i=0;i<expect;i++)
	{
		h_par[i]=-1;
	}
	cudaMalloc((void **)&d_par,expect*sizeof(int));
	for(circle=common;circle>=1;circle--)
	{
		if(have>=expect)
			break;
		storeL=headL;
		storeS=headS;
		while(storeL->pos<=storeS->pos+storeS->len)
		{
			storeL=storeL->next;
		}
		flag=0;
		while(storeS)
		{
			if(have>=expect)
				break;
			if(flag==1)  //don't have enough primers
				break;
			for(i=0;i<7;i++)
			{
				num[i]=0;
			}
			memory=0;
		//statistics	
			tempL=storeL;
			tempS=storeS;
			while(tempS&&(memory<prop.totalGlobalMem/4)&&(num[0]<k))
			{
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
					num[0]++;
					num[1]++;
					num[3]=num[3]+tempL->total_common;
					num[4]=num[4]+tempL->total_special;
					memory=memory+10+4*(tempL->total_common+tempL->total_special);
					tempL=tempL->next;
				}
				if(num[0]==0)
				{
					begin=tempS->pos;
					stop=tempS->pos+tempS->len-1;
					memory=memory+tempS->len;
				}
				else
				{
					if(tempS->pos+tempS->len-1>stop)
					{
						memory=memory+(tempS->pos+tempS->len-1-stop);
						stop=tempS->pos+tempS->len-1;
					}
				}
				num[0]++;
				num[2]++;
				num[5]=num[5]+tempS->total_common;
				num[6]=num[6]+tempS->total_special;
				memory=memory+22+4*(tempS->total_common+tempS->total_special)+common;
				tempS=tempS->next;
			}
			if(num[2]<4||num[1]<2) //don't have enough primers
			{
				flag=1;
				break;
			}
			if(tempS==NULL)  //check all primers
				flag=1;	
			h_seq=(char *)malloc((stop-begin+1)*sizeof(char));
			cudaMalloc((void **)&d_seq,(stop-begin+1)*sizeof(char));
			for(i=begin;i<=stop;i++)
				h_seq[i-begin]=seq[i];

		//malloc small primer
			h_primerS=(int *)malloc(4*num[2]*sizeof(int));
			cudaMalloc((void **)&d_primerS,4*num[2]*sizeof(int));
			h_commonS=(int *)malloc(4*num[5]*sizeof(int));
			cudaMalloc((void **)&d_commonS,4*num[5]*sizeof(int));
			h_specialS=(int *)malloc(4*num[6]*sizeof(int));
			cudaMalloc((void **)&d_specialS,4*num[6]*sizeof(int));
			cudaMalloc((void **)&d_nextS,num[2]*sizeof(int));
			cudaMalloc((void **)&d_other_to_L,num[2]*sizeof(int));

			h_start_commonS=(int *)malloc(num[2]*sizeof(int));
			cudaMalloc((void **)&d_start_commonS,num[2]*sizeof(int));
			h_end_commonS=(int *)malloc(num[2]*sizeof(int));
        	        cudaMalloc((void **)&d_end_commonS,num[2]*sizeof(int));
			h_start_specialS=(int *)malloc(num[2]*sizeof(int));
        	        cudaMalloc((void **)&d_start_specialS,num[2]*sizeof(int));
        	        h_end_specialS=(int *)malloc(num[2]*sizeof(int));
        	        cudaMalloc((void **)&d_end_specialS,num[2]*sizeof(int));
		
			tempS=storeS;
			for(i=0;i<3;i++)
			{
				count[i]=0;
			}
			while(count[0]<num[2])
			{
				if(tempS->total<circle)
				{
					tempS=tempS->next;
					continue;
				}
		//primer info
				h_primerS[4*count[0]]=tempS->pos;
				h_primerS[4*count[0]+1]=tempS->len;
				h_primerS[4*count[0]+2]=tempS->plus;
				h_primerS[4*count[0]+3]=tempS->minus;
			//common
				h_start_commonS[count[0]]=count[1];
				if(tempS->total_common==0)
					h_end_commonS[count[0]]=-1;
				else
				{
					p_node=tempS->common;
					while(p_node)
					{
						h_commonS[4*count[1]]=p_node->gi;
						h_commonS[4*count[1]+1]=p_node->pos;
						h_commonS[4*count[1]+2]=p_node->plus; 
                                	        h_commonS[4*count[1]+3]=p_node->minus;
						count[1]++;
						p_node=p_node->next;
					}
					h_end_commonS[count[0]]=count[1];
				}
			//special
				h_start_specialS[count[0]]=count[2];
                	        if(tempS->total_special==0)
                	                h_end_specialS[count[0]]=-1;
                	        else
                	        {
                	                p_node=tempS->special;
                	                while(p_node)
                	                {
                	                        h_specialS[4*count[2]]=p_node->gi;
                	                        h_specialS[4*count[2]+1]=p_node->pos;
                	                        h_specialS[4*count[2]+2]=p_node->plus;
                	                        h_specialS[4*count[2]+3]=p_node->minus;
                	                        count[2]++;
                	                        p_node=p_node->next;
                	                }
                	                h_end_specialS[count[0]]=count[2];
                	        }
				count[0]++;
				tempS=tempS->next;
			}
		//copy
			cudaMemcpy(d_primerS,h_primerS,num[2]*4*sizeof(int),cudaMemcpyHostToDevice);
			free(h_primerS);
			cudaMemcpy(d_commonS,h_commonS,num[5]*4*sizeof(int),cudaMemcpyHostToDevice);
			free(h_commonS);
			cudaMemcpy(d_specialS,h_specialS,num[6]*4*sizeof(int),cudaMemcpyHostToDevice);
        	        free(h_specialS);
			cudaMemcpy(d_start_commonS,h_start_commonS,num[2]*sizeof(int),cudaMemcpyHostToDevice);
			free(h_start_commonS);
			cudaMemcpy(d_end_commonS,h_end_commonS,num[2]*sizeof(int),cudaMemcpyHostToDevice);
        	        free(h_end_commonS);
			cudaMemcpy(d_start_specialS,h_start_specialS,num[2]*sizeof(int),cudaMemcpyHostToDevice);
        	        free(h_start_specialS);
        	        cudaMemcpy(d_end_specialS,h_end_specialS,num[2]*sizeof(int),cudaMemcpyHostToDevice);    
        	        free(h_end_specialS);
		//large primer
			h_primerL=(int *)malloc(4*num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_primerL,4*num[1]*sizeof(int));
        	        h_commonL=(int *)malloc(4*num[3]*sizeof(int));
        	        cudaMalloc((void **)&d_commonL,4*num[3]*sizeof(int));
        	        h_specialL=(int *)malloc(4*num[4]*sizeof(int));
        	        cudaMalloc((void **)&d_specialL,4*num[4]*sizeof(int));
        	        cudaMalloc((void **)&d_nextL,num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_other_to_S,num[1]*sizeof(int));

        	        h_start_commonL=(int *)malloc(num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_start_commonL,num[1]*sizeof(int));
        	        h_end_commonL=(int *)malloc(num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_end_commonL,num[1]*sizeof(int));
        	        h_start_specialL=(int *)malloc(num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_start_specialL,num[1]*sizeof(int));
        	        h_end_specialL=(int *)malloc(num[1]*sizeof(int));
        	        cudaMalloc((void **)&d_end_specialL,num[1]*sizeof(int));
                
        	        tempL=storeL;
        	        for(i=0;i<3;i++)
        	        {
        	                count[i]=0;
        	        }
        	        while(count[0]<num[1])
        	        {
				if(tempL->total<circle)
				{
					tempL=tempL->next;
					continue;
				}
                	//primer info
                	        h_primerL[4*count[0]]=tempL->pos;
                	        h_primerL[4*count[0]+1]=tempL->len;
                	        h_primerL[4*count[0]+2]=tempL->plus;
                	        h_primerL[4*count[0]+3]=tempL->minus;
                	//common
                	        h_start_commonL[count[0]]=count[1];
                	        if(tempL->total_common==0)
                	                h_end_commonL[count[0]]=-1;
                	        else
                	        {
                	                p_node=tempL->common;
                	                while(p_node)
                	                {
                	                        h_commonL[4*count[1]]=p_node->gi;
                	                        h_commonL[4*count[1]+1]=p_node->pos;
                	                        h_commonL[4*count[1]+2]=p_node->plus; 
                        	                h_commonL[4*count[1]+3]=p_node->minus;
                        	                count[1]++;
                        	                p_node=p_node->next;
                        	        }
                        	        h_end_commonL[count[0]]=count[1];
                        	}
                	//special
                        	h_start_specialL[count[0]]=count[2];
                        	if(tempL->total_special==0)
                        	        h_end_specialL[count[0]]=-1;
                        	else
                        	{
                        	        p_node=tempL->special;
                        	        while(p_node)
                        	        {
                        	                h_specialL[4*count[2]]=p_node->gi;
                        	                h_specialL[4*count[2]+1]=p_node->pos;
                        	                h_specialL[4*count[2]+2]=p_node->plus;
                        	                h_specialL[4*count[2]+3]=p_node->minus;
                        	                count[2]++;
                        	                p_node=p_node->next;
                        	        }
                        	        h_end_specialL[count[0]]=count[2];
                        	}
                        	count[0]++;
                        	tempL=tempL->next;
                	}
        	//copy
        	        cudaMemcpy(d_primerL,h_primerL,num[1]*4*sizeof(int),cudaMemcpyHostToDevice);
			free(h_primerL);
                	cudaMemcpy(d_commonL,h_commonL,num[3]*4*sizeof(int),cudaMemcpyHostToDevice);
                	free(h_commonL);
                	cudaMemcpy(d_specialL,h_specialL,num[4]*4*sizeof(int),cudaMemcpyHostToDevice);
                	free(h_specialL);
                	cudaMemcpy(d_start_commonL,h_start_commonL,num[1]*sizeof(int),cudaMemcpyHostToDevice);
                	free(h_start_commonL);
                	cudaMemcpy(d_end_commonL,h_end_commonL,num[1]*sizeof(int),cudaMemcpyHostToDevice);
                	free(h_end_commonL);
                	cudaMemcpy(d_start_specialL,h_start_specialL,num[1]*sizeof(int),cudaMemcpyHostToDevice);
                	free(h_start_specialL);
                	cudaMemcpy(d_end_specialL,h_end_specialL,num[1]*sizeof(int),cudaMemcpyHostToDevice);    
                	free(h_end_specialL);
	//run
			if(num[2]%prop.maxThreadsPerBlock==0)
				block=num[2]/prop.maxThreadsPerBlock;
			else
				block=(num[2]-num[2]%prop.maxThreadsPerBlock)/prop.maxThreadsPerBlock+1;

		//next primer	
			next_one<<<block,prop.maxThreadsPerBlock>>>(d_primerS,d_primerS,d_nextS,num[2],num[2]);
			next_one<<<block,prop.maxThreadsPerBlock>>>(d_primerL,d_primerL,d_nextL,num[1],num[1]);
			next_one<<<block,prop.maxThreadsPerBlock>>>(d_primerS,d_primerL,d_other_to_L,num[2],num[1]);
        	        next_one<<<block,prop.maxThreadsPerBlock>>>(d_primerL,d_primerS,d_other_to_S,num[1],num[2]);

		//calculate
			cudaMalloc((void **)&d_apply,common*num[2]*sizeof(int));
			cudaMalloc((void **)&d_result,12*num[2]*sizeof(int));
			cudaMemset(d_result,'\0',12*num[2]*sizeof(int));
			cudaMemcpy(d_par,h_par,expect*sizeof(int),cudaMemcpyHostToDevice);
			cudaMemcpy(d_seq,h_seq,(stop-begin+1)*sizeof(char),cudaMemcpyHostToDevice);
			free(h_seq);
			find_primer<<<block,prop.maxThreadsPerBlock>>>(d_primerL,d_primerS,d_specialL,d_specialS,d_commonL,d_commonS,d_start_specialL,d_start_specialS,d_start_commonL,d_start_commonS,d_end_specialL,d_end_specialS,d_end_commonL,d_end_commonS,d_nextS,d_nextL,d_other_to_L,d_other_to_S,d_result,num[1],num[2],common,d_apply,d_par,circle,expect,d_seq,begin);
//	exit(0);

	//	cudaError_t result0 = cudaGetLastError();
	//	printf("%s\n",result0);
			h_result=(int *)malloc(12*num[2]*sizeof(int));
			cudaMemcpy(h_result,d_result,12*num[2]*sizeof(int),cudaMemcpyDeviceToHost);
			h_apply=(int *)malloc(common*num[2]*sizeof(int));

		//check secondary structure
			fp=fopen(file1,"w");
			if(fp==NULL)
			{
				printf("Can't create the temp-par.txt file!\n");
				exit(1);
			}
			j=0; //as a flag, how many primers
			for(i=0;i<num[2];i++)
			{
				if(h_result[12*i+1]==0)
					continue;
				primer(seq,F3,h_result[12*i],h_result[12*i+1],0);
				primer(seq,F2,h_result[12*i+2],h_result[12*i+3],0);
				primer(seq,F1c,h_result[12*i+4],h_result[12*i+5],1);
				primer(seq,B1c,h_result[12*i+6],h_result[12*i+7],0);
                                primer(seq,B2,h_result[12*i+8],h_result[12*i+9],1);
                                primer(seq,B3,h_result[12*i+10],h_result[12*i+11],1);
				p_primer[0]=F3;
				p_primer[1]=F2;
				p_primer[2]=F1c;
				p_primer[3]=B1c;
				p_primer[4]=B2;
				p_primer[5]=B3;
			//output
				for(m=0;m<5;m++)
				{
					for(n=m+1;n<6;n++)
					{
						fprintf(fp,"SEQUENCE_ID=%d-%d-%d\n",i,h_result[12*i],(m*10+n));
						fprintf(fp,"PRIMER_TASK=check_primers\nPRIMER_PICK_ANYWAY=1\nPRIMER_SALT_DIVALENT=4\nPRIMER_DNTP_CONC=1.4\nPRIMER_DNA_CONC=38\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=%sprimer3_config/\n",path_primer3);
						fprintf(fp,"SEQUENCE_PRIMER=%s\n",p_primer[m]);
						fprintf(fp,"SEQUENCE_PRIMER_REVCOMP=%s\n=\n",p_primer[n]);
					}
				}
				j++;
			}
			fclose(fp);
			if(j==0) //don't have any candidate LAMP primers
			{
				remove(file1);
				continue;
			}
			memset(script,'\0',6000*sizeof(char));
			sprintf(script,"%sprimer3_core -strict_tags %s > %s",path_primer3,file1,file2);
			system(script);
			remove(file1);

			cudaMemcpy(h_apply,d_apply,common*num[2]*sizeof(int),cudaMemcpyDeviceToHost);
			have=secondary_check(expect,have,circle,h_par,file2,LAMP,h_result,h_apply,seq,common,list);
			remove(file2);
	//free
			cudaFree(d_seq);
			cudaFree(d_primerS);
			cudaFree(d_primerL);
			cudaFree(d_commonS);
			cudaFree(d_commonL);
			cudaFree(d_specialS);
			cudaFree(d_specialL);
			cudaFree(d_nextS);
			cudaFree(d_nextL);
			cudaFree(d_other_to_S);
			cudaFree(d_other_to_L);
			cudaFree(d_start_commonS);
			cudaFree(d_start_commonL);
			cudaFree(d_start_specialS);
			cudaFree(d_start_specialL);
			cudaFree(d_end_commonS);
			cudaFree(d_end_specialS);
			cudaFree(d_end_commonL);
			cudaFree(d_end_specialL);
			cudaFree(d_apply);
			cudaFree(d_result);
			free(h_result);
			free(h_apply);
		//new primer start
			if(tempS==NULL)
				storeS=tempS;
			else
			{
				while(tempS->pos-storeS->pos>300)
				{
					storeS=storeS->next;
				}
				while(storeL&&(storeL->pos<=storeS->pos+storeS->len))
				{
					storeL=storeL->next;
				}
			}
		}//one circle
	}
	cudaFree(d_par);
	free(h_par);
	end=time(NULL);
	printf("the time for design is %0.1f seconds!\n",difftime(end,start));
	fclose(LAMP);
	free(seq);
	free(list);
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
		
		tempL=headL->next;
		free(headL);
		headL=tempL;
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

                tempS=headS->next;
                free(headS);
                headS=tempS;
        }
}
