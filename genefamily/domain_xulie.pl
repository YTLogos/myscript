die "perl $0 <hmmoutfile> <fa> <OUT> <E-value>" unless(@ARGV==4 );
use Math::BigFloat;
use Bio::SeqIO;
use Bio::Seq;
$in = Bio::SeqIO -> new(-file => "$ARGV[1]",
                                  -format => 'Fasta');
$out = Bio::SeqIO -> new(-file => ">$ARGV[2]",
                                  -format => 'Fasta');
my %keep=() ;
open IN,"$ARGV[0]" or die "$!";
while (<IN>) {
		chomp;
		next if /^#/;
		
		my @a= split /\s+/;
		next if $a[6] > $ARGV[3] ;
		my @b=($a[17],$a[18]);
		my $keys = $a[0];
		if (!exists $keep{$keys} ) {
			$keep{$keys} = \@b ;
	
		}
		
}
close (IN);
while ( my $seq = $in->next_seq() ) {
     my($id,$sequence,$desc)=($seq->id,$seq->seq,$seq->desc);

     if(exists $keep{$id}){
               my$newseqobj=$seq->trunc($keep{$id}[0],$keep{$id}[1]);
     	
     	$out->write_seq($newseqobj);      	
     }
}
$in->close();
$out->close();
