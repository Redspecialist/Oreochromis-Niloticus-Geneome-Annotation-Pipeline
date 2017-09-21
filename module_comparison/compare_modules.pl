#Compares the contents of modules from two WGCNA runs and creates a table of the overlaps between those file outputs
#Vertical axis represents 
use strict;
use Fcntl;
use warnings;
use List::MoreUtils qw(uniq);
use List::Util qw(first);

my $full_genome = $ARGV[0];
my $partial_genome = $ARGV[1];
my $out = "$ARGV[2].module_comparison_matrix.out";
my $second_output = "$ARGV[2].comparison_sub_modules.out";

sysopen(my $FULL, $full_genome, O_RDONLY);

my %accession_to_fullmod;
my %all_mods;
<$FULL>;
#loads a mapping of each module to loc IDs based on the first inputs reference file
while(my $line = <$FULL>){
   
   my @terms = split(/\t+/, $line);
   if(!(scalar @terms < 5)){
       chomp $terms[4];
       $all_mods{$terms[4]} = "1"; 
       $accession_to_fullmod{$terms[0]} = $terms[4];
   }
}
close($FULL);

sysopen(my $PART, $partial_genome, O_RDONLY);

my %accession_to_partmod;
my %mod_to_access;
my %loc_to_info;
<$PART>;

#loads a mapping of each module to loc IDs based on the second inputs reference file  
while(my $line = <$PART>){
    chomp $line;
    my @terms = split(/\t+/, $line);
   
    if(!((scalar @terms) < 5)){
	$accession_to_partmod{$terms[0]} = $terms[4];
	push (@{$mod_to_access{$terms[4]}}, $terms[0]); 
	$loc_to_info{$terms[0]} = $terms[1]."\t".$terms[2]."\t".$terms[3];
    }
}

close $PART;

my @header = keys %all_mods;
sysopen(my $OUTPUT, $out, O_WRONLY | O_CREAT);
sysopen(my $SECOND_OUT, $second_output, O_WRONLY | O_CREAT);

print $OUTPUT "MODULE\tTOTAL\tUNKNOWN";

# outputs theremaining portions of the header consisting of modules in the first input 
foreach my $elem (@header){
    print $OUTPUT "\t$elem";
}
print $OUTPUT "\n";

my @locs = ("LOC100694761", "LOC100695018", "LOC100695287", "LOC100710676", "LOC100710942","LOC100711209", "LOC100710249","LOC100711443","LOC100702555","LOC100705124");
my %names = ("LWS" => "LOC100694761", "SWS2B" => "LOC100695018", "SWS2A" => "LOC100695287", "RH2A-beta" => "LOC100710676","RH2A-alpha" => "LOC100710942","RH2B" => "LOC100711209","SWS1" => "LOC100710249","Rhodopsin" => "LOC100711443","Tbx2a" => "LOC100702555","Rx1" => "LOC100705124");
my %intersections = map {$_ => "not_found"} @locs;

#parses each module from the first input to comparing them to each module in the second input. (Works in quadratic time)
foreach my $module (keys %mod_to_access){
    
    my $num_genes = @{$mod_to_access{$module}};
    my %mod_intersect;
    my $num_match = 0;

    
    foreach my $accession (@{$mod_to_access{$module}}) {
	if(exists $accession_to_fullmod{$accession}){    
	    my $match = $accession_to_fullmod{$accession};
	    push(@{$mod_intersect{$match}}, $accession);
	    $num_match += 1;
	} 
    }
    
    my $unk = $num_genes - $num_match;
    $unk = "-" if($unk == 0);
    chomp $module;
    print $OUTPUT "$module";
    print $OUTPUT "\t$num_genes\t$unk";
    
    foreach my $mod (@header){
	my $val = "-";
	
	if(exists $mod_intersect{$mod}){
	    
	    $val = scalar(@{$mod_intersect{$mod}});
	    print $SECOND_OUT ">$mod-X-$module\n";

	    #prints all of the contenta of a given module intersction
	    foreach my $accession (@{$mod_intersect{$mod}}){
		
		$intersections{$accession} = "$mod-X-$module" if(exists($intersections{$accession})); 		
		print $SECOND_OUT "$accession\t$loc_to_info{$accession}\t$mod-X-$module\n";
	    
	    }
	    print $SECOND_OUT "\n";
	
	}
	
	print $OUTPUT "\t$val";
	
    }
    print $OUTPUT "\n";
}
#prints out where each opsin and other important genes are located in the module intersections
foreach my $name (keys %names){
    print $OUTPUT "$name is found in the module intersection of $intersections{$names{$name}}\n";
}
    
close($OUTPUT);
close($SECOND_OUT);
