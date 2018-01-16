use strict; use warnings;
use Getopt::Long;

##used in read parameters
my $prefix;
my $common_file;
my $special_file;
my $par_common;
my $par_special;
my $flag;
my $help;
my $path_bowtie;
my $index;
my $high;
my $low;
my $loop;
my $dir;
my $left;
my $check;
my $thread;
my $ref;
my $seq="";

my $line;
my @array;
my $pos;
my $len;
my $i;
my $j;
my $primer;
my $name;
my $mismatch;
my @list;
my @list_index;
my %common;
my %special;
my @list_common;
my @file;
my @status;
my $num;
my $start;
my $end;
my $take;
my $num_common=0;
my $num_special=0;

$help="usage is: perl $0 parameters\n
	--in <string>: The name of candidate primer files, generate from the single program.
	--ref <string>: The reference genome fasta file, used in designing single primers.
	--dir <string>: The directory to store candidate single primers. default is current directory.
	--common <string>: The file contains the name of genomes, the LAMP primers should be used in those genomes.
	--special <string>: The file contains the name of genomes, the LAMP primers should not be used in those genomes.
	--left: The LAMP primers should not be used in all the genomes beside in common file.
	--bowtie <string>: The bowtie program (include the directory).
	--index <string>: Comma-separated index files used in bowtie.
	--par_special <int>: The max mismatches in primer, 0-3, default is 2. The bigger, the more strict in specifity.
	--par_common <int>: The max mismatches in primer, 0-3, default is 0 and <= par_special. The smaller, the more strict in common.
	--threads <int>: Number of alignment threads to launch used in bowtie, default is 1.
	--help|--h Print help information.
	--high/--low: Design LAMP primers in high/low GC region. high: the GC content>=45\%; low: the GC content <=45\%.
	--loop: Design LAMP primers with loop single primers.\n";

$check=GetOptions("in=s"=>\$prefix,"ref=s"=>\$ref,"dir=s"=>\$dir,"common=s"=>\$common_file,"special=s"=>\$special_file,"left"=>\$left,"bowtie=s"=>\$path_bowtie,"index=s"=>\$index,"par_common=i"=>\$par_common,"par_special=i"=>\$par_special,"threads=i"=>\$thread,"help|h"=>\$flag,"high"=>\$high,"low"=>\$low,"loop"=>\$loop);

if($check==0)
{
	print $help;
	exit;
}
##check input paramters
if($flag)
{
        print $help;
        exit;
}
if(!$prefix)
{
	print "Error!\nDon't have the --in! The name of candidate primer files should be supplied.\n";
	print $help;
	exit;
}

if(!$ref)
{
        print "Error!\nDon't have the --ref! The reference genome file should be supplied.\n";
        print $help;
        exit;
}

if(!$dir)
{
	$dir=$ENV{'PWD'};
	$dir=$dir."/";
}
else
{
	if($dir!~/\/$/)
	{
		$dir=$dir."/";
	}
}

if((!$common_file) && (!$special_file) && (!$left))
{
	print "Error!\nOne of --common, --special and --left should be existed.\n";
	exit;
}

if(!$path_bowtie)
{
	print "Error!\nDon't have the --bowtie! The bowtie program should be supplied.\n";
	print $help;
	exit;
}
else
{
	if(!(-e $path_bowtie))
	{
		print "Don't have the bowtie program in $path_bowtie!\n";
		exit;
	}
}

if(!$index)
{
	print "Error!\nDon't have --index! Comma-separated index files used in bowtie.\n";
	print $help;
	exit;
}

if(!$par_common)
{
	$par_common=0;
}
else
{
	if($par_common<0 ||$par_common>3)
	{
		print "The argument of --par_common must be int between 0 and 3!\n";
		exit;
	}
}

if(!$par_special)
{
        $par_special=2;
}
else
{
        if($par_special<0 || $par_special>3) 
        {
                print "The argument of --par_special must be int between 0 and 3!\n";
                exit;
        }
}

if($par_common>$par_special)
{
	print "Error!\nThe argument of --par_common (now is $par_common) must <= the argument of --par_special (now is $par_special)!\n";
	exit;
}

if((!$special_file) && (!$left))
{
	$par_special=$par_common;
}##don't take care of special
if(!$thread)
{
	$thread=1;
}

if(!$high && !$low)
{
	print "Error!\nOne of --high or --low must be existed!\n";
	exit;
}
if($high && $low)
{
	print "Error!\nOnly one of --high or --low must be existed!\n";
	exit;
}

###################################prepare
if($common_file)
{
	open(IN,"<$common_file") or die "Can't open the $common_file inputed from --common!\n";
	while(<IN>)
	{
		chomp;
		$line=$_;
		if($line=~/^\>/)
		{
			($line)=$line=~/^\>(.+)$/;
		}
		@array=split(" ",$line);
		if(length($array[0])<=300)
		{
			if(exists $common{$array[0]})
			{
				print "Warnings: in $common_file, the name of \"$array[0]\" exists in different genomes.\n";
				next;
			}
			$common{$array[0]}=$num_common;
			$list_common[$num_common]=$array[0];
			$num_common++;
		}
		else
		{
			$flag=substr($array[0],0,300);
			if(exists $common{$flag})
			{
				print "Warnings: in $common_file, the name of \"$flag\" exists in different genomes.\n";
				next;
			}
			$common{$flag}=$num_common;
			$list_common[$num_common]=$flag;
			$num_common++;
		}
	}
	close IN;
}

if($special_file)
{
	open(IN,"<$special_file") or die "Can't open the $special_file file inputed from --special!\n";
	while(<IN>)
	{
		chomp;
		$line=$_;
		if($line=~/^\>/)
                {
                        ($line)=$line=~/^\>(.+)$/;
                }
                @array=split(" ",$line);
		if(length($array[0])<=300)
		{
			if(exists $special{$array[0]})
			{
				print "Warnings: in $special_file, the name of \"$array[0]\" exists in different genomes.\n";
				next;
			}
			$special{$array[0]}=$num_special;
			$num_special++
		}
		else
		{
			$flag=substr($array[0],0,300);
			if(exists $special{$flag})
			{
				print "Warnings: in $special_file, the name of \"$flag\" exists in different genomes.\n";
				next;
			}
			$special{$flag}=$num_special;
			$num_special++;
		}
	}
	close IN;
}

$num=2;
if($high)
{
	$list[0]=$dir."high-inner/".$prefix;
	$list[1]=$dir."high-outer/".$prefix;
	if($loop)
	{
		$num=3;
		$list[2]=$dir."high-loop/".$prefix;
	}
}
else
{
	$list[0]=$dir."low-inner/".$prefix;
        $list[1]=$dir."low-outer/".$prefix;
        if($loop)
        {
		$num=3;
                $list[2]=$dir."low-loop/".$prefix;
        }
}
@list_index=split(",",$index);

open(IN,"<$ref") or die "Can't open $ref file!\n";
while(<IN>)
{
        chomp;
        $line=$_;
        if($line=~/^\>/)
        {
                next;
        }
        $seq=$seq.$line;
}
close IN;

#####################################run bowtie and analysis
for($i=0;$i<@list;$i++)
{
	$flag=$i+1;
	print "Now the program is handling the $flag\-th file, total files is $num...\n";
	$start=time();
	open(IN,"<$list[$i]") or die "Can't open $list[$i] file!\n";
	$file[0]=$list[$i].".fa";
	open(OUT,">$file[0]") or die "Can't create $file[0] file!\n";
	while(<IN>)
	{
		chomp;
		$line=$_;
		($pos,$len,$status[0],$status[1])=$line=~/pos\:(\d+)\tlength\:(\d+)\t\+\:(\d)\t\-\:(\d)/;
        	$name=$pos."-".$len."-".$status[0]."-".$status[1];
		$primer=substr($seq,$pos,$len);
	        print OUT ">$name\n$primer\n";
	}
	close IN;
	close OUT;

	if($common_file)
	{
		if($i==0)
		{
			$file[1]=$list[$i]."-common_list.txt";
			open(LIST,">$file[1]") or die "Can't create $file[1] file!\n";
			for($j=0;$j<$num_common;$j++)
			{
				print LIST "$list_common[$j]\t$j\n";
			}
			close LIST;
		}
		$file[1]=$list[$i]."-common.txt";
		open(COMMON,">$file[1]") or die "Can't create $file[1] file!\n";
	}
	if($special_file || $left)
	{
		$file[2]=$list[$i]."-special.txt";
		open(SPECIAL,">$file[2]") or die "Can't create $file[2] file!\n";
	}

	##run bowtie
	for($j=0;$j<@list_index;$j++)
	{
		$file[3]=$list[$i]."-".$j.".bowtie";
		system("$path_bowtie -f --suppress 5,6,7 -v $par_special -p $thread -a $list_index[$j] $file[0] $file[3]");

		open(BOWTIE,"<$file[3]") or die "Can't open the $file[3] file!\n";
		while(<BOWTIE>)
		{
       			$line=$_;
			@array=split("\t",$line);

			if(length($array[2])<=300)
			{
				$flag=$array[2];
			}
			else
			{
				$flag=substr($array[2],0,300);
			}
			$mismatch=$array[4]=~tr/:/:/;
			($pos,$len,$status[0],$status[1])=$array[0]=~/^(\d+)\-(\d+)\-(\d)\-(\d)/;
			$status[2]=0;
			$status[3]=0;
			if($array[1] eq "+")
               	        {
               	                if($status[0]==1)
               	                {
               	                        $status[2]=1;
               	                }
               	                if($status[1]==1)
               	                {
               	                        $status[3]=1;
               	                }
               	        }
               	        else
               	        {
               	                if($status[0]==1)
               	                {
               	                        $status[3]=1;
               	                }
               	                if($status[1]==1)
               	                {
               	                        $status[2]=1;
               	                }
               	        }
	
			if($common_file)
			{
				if(exists $common{$flag}&&$mismatch<=$par_common)
				{
					print COMMON "$pos\t$len\t$common{$flag}\t$array[3]\t$status[2]\t$status[3]\n";
					next;
				}
			}
			if($special_file)
			{
				if(exists $special{$flag})
				{
					print SPECIAL "$pos\t$len\t$special{$flag}\t$array[3]\t$status[2]\t$status[3]\n";
				}
				next;
			}
	
			if($left)
			{
				if(exists $special{$flag})
				{
					print SPECIAL "$pos\t$len\t$special{$flag}\t$array[3]\t$status[2]\t$status[3]\n";
				}
				else
				{
					$special{$flag}=$num_special;
					print SPECIAL "$pos\t$len\t$num_special\t$array[3]\t$status[2]\t$status[3]\n";
					$num_special++;
				}
			}
		}
		close BOWTIE;	
		system("rm $file[3]");
	}
	system("rm $file[0]");
	if($common_file)
	{
		close COMMON;
	}
	if($special_file || $left)
	{
		close SPECIAL;
	}
	$end=time();
	$take=$end-$start;
	print "    In this step, it takes $take seconds.\n";
}
