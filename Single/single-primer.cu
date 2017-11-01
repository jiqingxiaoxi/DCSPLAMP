#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<cuda.h>
#include<time.h>
#include<sys/stat.h>
#include<regex.h>
#include<unistd.h>

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
		if(seq[i]=='C'||seq[i]=='c')
		{
			number++;
			continue;
		}
	
		if(seq[i]=='G'||seq[i]=='g')
		{
			number++;
		}
	}

	gc=1.0*number/length*100;
	if((gc<40)||(gc>60))
	{
		return 0;
	}
	return 1;
}

///function in gpu, translate A...G to int
__device__ int translate(char a)
{
	if(a=='A'||a=='a')
		return 0;
	if(a=='T'||a=='t')
		return 1;
	if(a=='C'||a=='c')
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
	if((seq[0]=='A')||(seq[0]=='T')||seq[0]=='a'||seq[0]=='t')
	{
		deltah+=2.3;
		deltas+=4.1;
	}
	else
	{
		deltah+=0.1;
		deltas-=2.8;
	}
        if((seq[length-1]=='A')||(seq[length-1]=='T')||seq[length-1]=='a'||seq[length-1]=='t')
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
	int i,flag;
	
	flag=1;
	for(i=0;i<length;i++)
	{
		if(d_seq[position+i]=='A'||d_seq[position+i]=='a')
		{
			continue;
		}
		if(d_seq[position+i]=='T'||d_seq[position+i]=='t')
		{
			continue;
		}
		if(d_seq[position+i]=='C'||d_seq[position+i]=='c')
                {
                        continue;
                }
                if(d_seq[position+i]=='G'||d_seq[position+i]=='g')
                {
                        continue;
                }
		flag--;
		break;
	}
	return flag;
}

//function in gpu, reverse the strand,+ to - strand
__device__ void reverse(char seq[],char rev[],int length)
{
	int i;
	
	for(i=0;i<length;i++)
	{
		if(seq[length-1-i]=='A'||seq[length-1-i]=='a')
		{
			rev[i]='T';
			continue;
		}
                if(seq[length-1-i]=='T'||seq[length-1-i]=='t')
                {
                        rev[i]='A';
                        continue;
                }
                if(seq[length-1-i]=='C'||seq[length-1-i]=='c')
                {
                        rev[i]='G';
                        continue;
                }
		rev[i]='C';
	}
}

///function: int length: the length of genome
__global__ void candidate_primer(char *d_seq,int *d_pos,int *d_len,int *d_rev_len,float *d_stab,float *d_deltah,float *d_deltas,int strand,float max_tm,float min_tm,int length)
{
	int id,i,circle,check;
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
                                continue;

			generate(d_seq,primer,circle,i);
			check=gc(primer,i);
			if(check==0)
				continue;

			check=tm(primer,d_deltah,d_deltas,i,max_tm,min_tm);
			if(check==0)
				continue;

                        check=stability(primer,d_stab,i,strand);
                        if(check==1)     //+ strand
                                d_len[circle*8+i-18]=1;
	
			reverse(primer,rev,i);  //generate - strand
			check=stability(rev,d_stab,i,strand);
			if(check==1)
				d_rev_len[circle*8+i-18]=1;
		}
		
		for(i=0;i<8;i++)
		{
			d_pos[circle]+=(d_len[circle*8+i]+d_rev_len[8*circle+i]);
		}
	}
	__syncthreads();
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

void cpu_reverse(char seq[],char rev[],int length)
{
        int i;
        
        for(i=0;i<length;i++)
        {
                if(seq[length-1-i]=='A'||seq[length-1-i]=='a')
                {
                        rev[i]='T';
                        continue;
                }
                if(seq[length-1-i]=='T'||seq[length-1-i]=='t')
                {
                        rev[i]='A';
                        continue;
                }
                if(seq[length-1-i]=='C'||seq[length-1-i]=='c')
                {
                        rev[i]='G';
                        continue;
                }
                rev[i]='C';
        }
        rev[i]='\0';
}

int second_structure(char *prefix,char *dir,float tm_threshold,char *primer3)
{
	char *in,*out,*script,*line,result[20],pattern1[50],pattern2[50],pattern3[50];
	FILE *fp1,*fp2;
	int pos,len,plus,minus,total,flag,len_now,pos_now,status,cflags,length;
	regex_t reg1,reg2,reg3;
	regmatch_t pmatch[3];
	float TH;

	len=strlen(prefix)+strlen(dir)+20;
	in=(char *)malloc(len*sizeof(char));
	memset(in,'\0',len*sizeof(char));
	strcpy(in,dir);
	strcat(in,prefix);
	strcat(in,"-primer3.txt");
//run primer3
	out=(char *)malloc(len*sizeof(char));
	memset(out,'\0',len*sizeof(char));
	strcpy(out,dir);
	strcat(out,prefix);
	strcat(out,"-result.txt");

	len=strlen(primer3)+strlen(in)+strlen(out)+50;
	script=(char *)malloc(len*sizeof(char));
	memset(script,'\0',len*sizeof(char));
	sprintf(script,"%s -strict_tags %s > %s",primer3,in,out);
	system(script);
	free(script);
	remove(in);

//check structure
	len=strlen(prefix)+strlen(dir)+20;
	memset(in,'\0',len*sizeof(char));
	strcpy(in,out);
	memset(out,'\0',len*sizeof(char));
	strcpy(out,dir);
	strcat(out,prefix);
	fp1=fopen(in,"r");
	if(fp1==NULL)
	{
		printf("Can't open %s file!\n",in);
		exit(1);
	}
	fp2=fopen(out,"w");
	if(fp2==NULL)
	{
		printf("Can't create %s file!\n",out);
		exit(1);
	}

//prepare regular
	strcpy(pattern1,"ID=(\\w+)-(\\w+)-(.)");
	cflags=REG_EXTENDED;
	regcomp(&reg1,pattern1,cflags);
	strcpy(pattern2,"_TH=(.+)");
	regcomp(&reg2,pattern2,cflags);
	strcpy(pattern3,"Contains too-long poly nucleotide tract");
	regcomp(&reg3,pattern3,cflags);

//read file
	pos=0;
	len=0;
	plus=0;
	minus=0;
	total=0;
	length=200+strlen(primer3);
	line=(char *)malloc(length*sizeof(char));
	memset(line,'\0',length*sizeof(char));
	while(fgets(line,length*sizeof(char),fp1)!=NULL)
	{
		if(regexec(&reg1,line,4,pmatch,0)==0)  //begin
		{
			take_regulate(pmatch,1,result,line);
			pos_now=atoi(result); //pos
			take_regulate(pmatch,2,result,line);
			len_now=atoi(result);

			if(pos_now!=pos||len_now!=len) //new one
			{
				if(plus+minus!=0)
				{
					fprintf(fp2,"pos:%d\tlength:%d\t+:%d\t-:%d\n",pos,len,plus,minus);
					total++;
				}
				pos=pos_now;
				len=len_now;
				plus=0;
				minus=0;
			}
			take_regulate(pmatch,3,result,line);
			if(result[0]=='+')
                        {
                                flag=1;
                        }
                        else
                        {
                                flag=0;
                        }
                        status=1;
			continue;
                }

                if(regexec(&reg2,line,2,pmatch,0)==0&&status==1) //the max TH
                {
                        take_regulate(pmatch,1,result,line);
                        TH=atof(result);
                        if(TH>tm_threshold)
                                status=0;
                        continue;
                }
                if(regexec(&reg3,line,1,pmatch,0)==0)
                {
                        status=0;
                        continue;                       
                }
                if(line[0]=='='&&status==1)
                {
                        if(flag==1)
                                plus=1;
                        else
                                minus=1;
                }
        }
        if(plus+minus!=0)
        {
                fprintf(fp2,"pos:%d\tlength:%d\t+:%d\t-:%d\n",pos,len,plus,minus);
                total++;
        }
        fclose(fp1);
        fclose(fp2);
	remove(in);
	free(out);
	free(in);
	free(line);
	regfree(&reg1);
	regfree(&reg2);
	regfree(&reg3);
        return total;
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
        printf("    -P    <string>:  the primer3_core program in Primer3 software. If add this parameter, the second structure would be checked using primer3\n");
        printf("    -high/-low:      design candidate single primers in high/low GC region. high: the GC content>=45%%; low: the GC content <=45%%.\n");
        printf("    -loop:           design candidate loop single primers\n");
        printf("    -h/-help:        print usage\n");
}

void cpu_generate(char *seq,char out[],int pos,int length)
{
        int i;
        for(i=0;i<length;i++)
        {
                out[i]=seq[pos+i];
        }
        out[i]='\0';
}

int create_file(char *prefix,char *dir,char *seq,int *pos,int *len,int *rev_len,int length,int P_flag,char *primer3_config)
{
	char primer[26],rev[26],*file;
	int total,i,j;
	FILE *OUT;

	total=0;
	i=strlen(dir)+strlen(prefix)+20;
	file=(char *)malloc(i);
        memset(file,'\0',i);
        strcpy(file,dir);
        strcat(file,prefix);
	if(P_flag==1)
		strcat(file,"-primer3.txt");
        OUT=fopen(file,"w");
        if(OUT==NULL)
        {
                printf("Error! Can't create the %s file!\n",file);
                exit(1);
        }
	
	if(P_flag)
	{
		fprintf(OUT,"PRIMER_TASK=check_primers\nPRIMER_PICK_ANYWAY=1\nPRIMER_SALT_DIVALENT=4\nPRIMER_DNTP_CONC=1.4\nPRIMER_DNA_CONC=38\n");
		fprintf(OUT,"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n",primer3_config);
	}
        for(i=0;i<length;i++)
        {
                if(pos[i]==0)
                        continue;
                for(j=0;j<8;j++)
                {
                        if((len[8*i+j]+rev_len[8*i+j])==0)
                                continue;
			if(P_flag==0)
                        	fprintf(OUT,"pos:%d\tlength:%d\t+:%d\t-:%d\n",i,(j+18),len[8*i+j],rev_len[8*i+j]);
			else
			{
				cpu_generate(seq,primer,i,(j+18));
				if(len[8*i+j]==1)
                		{
                        		fprintf(OUT,"SEQUENCE_ID=%d-%d-+\n",i,(j+18));
                		        fprintf(OUT,"SEQUENCE_PRIMER=%s\n=\n",primer);
                		}
                		if(rev_len[8*i+j]==1)
                		{
                        		cpu_reverse(primer,rev,(j+18));
                        		fprintf(OUT,"SEQUENCE_ID=%d-%d--\n",i,(j+18));
		                        fprintf(OUT,"SEQUENCE_PRIMER=%s\n=\n",rev);
        		        }
			}
			total++;
                }
        }
	fclose(OUT);
	free(file);
	return total;
}

main(int argc, char **argv)
{
	int *pos,*d_pos,*len,*d_len,length,flag[9],i,*rev_len,*d_rev_len,num_outer,num_inner,num_loop;
	float deltah[16],deltas[16],stab[4096],*d_deltah,*d_deltas,*d_stab,temp1,temp2;
	char *seq,*d_seq,*store_path,*prefix,*stab_path,*tm_path,*curren_path,*input,*outer,*inner,*loop,*primer3,*primer3_config,*temp;
	FILE *fp;
	time_t start,end;
        struct stat statbuf;
//flag: 0:input; 1: out_prefix; 2: dir; 3: stab; 4: tm; 5: high; 6: low; 7: loop; 8: primer3

	start=time(NULL);
//get input
        for(i=0;i<9;i++)
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
                else if(strcmp(argv[i],"-P")==0)
                {
                        flag[8]=1;
                        if(i+1==argc)
                        {
                                printf("Error! The \"-P\" parameter is not completed.\n");
                                usage();
                                exit(1);
                        }
			length=strlen(argv[i+1]);
                        primer3=(char *)malloc(length+1);
                        memset(primer3,'\0',length+1);
                        strcpy(primer3,argv[i+1]);
			if(access(primer3,0)==-1)
                        {
                                printf("Error! Don't have the %s program!\n",primer3);
                                exit(1);
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
//primer3 program
	if(flag[8]==1)
	{
		length=strlen(primer3);
                primer3_config=(char *)malloc(length+20);
                memset(primer3_config,'\0',length+20);
                strcpy(primer3_config,primer3);
                i=length-1;
        	while(primer3_config[i]!='/'&&i>=0)
        	{
        	        primer3_config[i]='\0';
        	        i--;
        	}
        	strcat(primer3_config,"primer3_config/");
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
                seq[length]=temp[i];
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

	end=time(NULL);
	printf("It takes %d seconds to prepare.\n",(int)difftime(end,start));
	start=time(NULL);

	if(flag[5]==1)
        {
		cudaMemset(d_pos,'\0',length*sizeof(int));
		cudaMemset(d_len,'\0',8*length*sizeof(int));
		cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
		candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,61,59,length);
		cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
        	cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
        	cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_outer=create_file(prefix,outer,seq,pos,len,rev_len,length,flag[8],primer3_config);

		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,0,66,64,length);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_inner=create_file(prefix,inner,seq,pos,len,rev_len,length,flag[8],primer3_config);

                if(flag[7]==1)
		{
			cudaMemset(d_pos,'\0',length*sizeof(int));
                	cudaMemset(d_len,'\0',8*length*sizeof(int));
                	cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                	candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,66,64,length);
                	cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	num_loop=create_file(prefix,loop,seq,pos,len,rev_len,length,flag[8],primer3_config);
		}
        }
        else
        {
		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,56,54,length);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_outer=create_file(prefix,outer,seq,pos,len,rev_len,length,flag[8],primer3_config);

		cudaMemset(d_pos,'\0',length*sizeof(int));
                cudaMemset(d_len,'\0',8*length*sizeof(int));
                cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,0,61,59,length);
                cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                num_inner=create_file(prefix,inner,seq,pos,len,rev_len,length,flag[8],primer3_config);
                if(flag[7]==1)
		{
			cudaMemset(d_pos,'\0',length*sizeof(int));
                	cudaMemset(d_len,'\0',8*length*sizeof(int));
                	cudaMemset(d_rev_len,'\0',8*length*sizeof(int));
                	candidate_primer<<<200,200>>>(d_seq,d_pos,d_len,d_rev_len,d_stab,d_deltah,d_deltas,1,61,59,length);
                	cudaMemcpy(pos,d_pos,length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(len,d_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	cudaMemcpy(rev_len,d_rev_len,8*length*sizeof(int),cudaMemcpyDeviceToHost);
                	num_loop=create_file(prefix,loop,seq,pos,len,rev_len,length,flag[8],primer3_config);
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
        printf("It takes %d seconds to design candidate single primers(without second structure check).\n",(int)difftime(end,start));
        start=time(NULL);

        if(flag[8]==0)
        {
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
                exit(1);
        }
//check the second structure
        if(flag[5]==1)
        {
		num_outer=second_structure(prefix,outer,49.0,primer3);
                num_inner=second_structure(prefix,inner,54.0,primer3);
                if(flag[7]==1)
                        num_loop=second_structure(prefix,loop,54.0,primer3);
        }
        else
        {
                num_outer=second_structure(prefix,outer,44.0,primer3);
                num_inner=second_structure(prefix,inner,49.0,primer3);
                if(flag[7]==1)
                        num_loop=second_structure(prefix,loop,49.0,primer3);
        }

        //check
        printf("After second structure check, there ara %d candidate primers used as F3/F2/B2/B3.\n",num_outer);
        printf("After second structure check, there are %d candidate primers used as F1c/B1c.\n",num_inner);
        if(flag[7]==1)
                printf("After second structure check, there are %d candidate primers used as LF/LB.\n",num_loop);
        //check
        if(num_outer<4)
                printf("Warning: there don't have enough primers(>=4) used as F3/F2/B2/B3.\n");
        if(num_inner<2)
                printf("Warning: there don't have enough primers(>=2) used as F1c/B1c.\n");
        if(flag[7]==1 && num_loop<1)
                printf("Warning: there don't have enough primers(>=1) used as LF/LB. But you can design LAMP primers without loop primer.\n");
	end=time(NULL);
	printf("It takes %d seconds to check single primers' second structure.\n",(int)difftime(end,start));

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
	free(primer3);
	free(primer3_config);
}
