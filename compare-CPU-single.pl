use strict; use warnings;

my $line;
my @array;
my @store;
my $id;
my $value;
my $temp;

if(@ARGV!=2)
{
	print "perl $0 ours  primer3\n";
	exit;
}

open(IN,"<$ARGV[0]") or die "$ARGV[0]\n";
while(<IN>)
{
	chomp;
	$line=$_;
	@array=split("\t",$line);
	$store[$array[0]][0]=$array[1]+0.005;
	$store[$array[0]][1]=$array[2]+0.005;
	$store[$array[0]][2]=$array[3]+0.005;
}
close IN;

open(IN,"<$ARGV[1]") or die "$ARGV[1]\n";
while(<IN>)
{
	chomp;
	$line=$_;
	if($line=~/SEQUENCE_ID/)
	{
		($id)=$line=~/\=(\d+)$/;
		$id=$id/20-1;
		next;
	}

	if($line=~/SELF_ANY_TH=/)
	{
		($value)=$line=~/\=(.+)$/;
		$temp=$store[$id][0]-$value;
		if($temp>0.01||$temp<-0.01)
		{
			$id=($id+1)*20;
			print "id is $id, self_any\n";
			$id=$id/20-1;
		}
		next;
	}
	if($line=~/SELF_END_TH\=/)
	{
		($value)=$line=~/\=(.+)$/;
		$temp=$store[$id][1]-$value;           
                if($temp>0.01||$temp<-0.01)
                {
			$id=($id+1)*20;
                        print "id is $id, self_end\n";
			$id=$id/20-1;
                }
                next;
        }
	if($line=~/HAIRPIN_TH\=/) 
        {
                ($value)=$line=~/\=(.+)$/;
                $temp=$store[$id][2]-$value;
                if($temp>0.01||$temp<-0.01)
                {
                        $id=($id+1)*20;
                        print "id is $id, hairpin\n";
                        $id=$id/20-1;
                }
        }
}
close IN;
