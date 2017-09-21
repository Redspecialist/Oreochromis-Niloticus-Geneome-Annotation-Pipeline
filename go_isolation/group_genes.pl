use strict;
use Fcntl;
use warnings;
use List::MoreUtils qw(uniq);

my $module_data = $ARGV[0];
my $GO_Data = "./go_isolation/library.go_table.txt";
my $out = $ARGV[1];

sysopen(my $MODREADER,$module_data,O_RDONLY);

readline $MODREADER;

my %loc_to_accession;
my %module_to_loc;
my %loc_to_info;

sub pointless{
    my $access = $_[0];

    if(not ((index $access, "ref\|") == -1)){
        $access = substr $access,(index $access, "ref\|") +4;
    }
    if(not ((index $access, ".") == -1)){
        $access = substr $access, 0, (index $access, ".");
    }else{
	$access = substr $access, 0, (index $access, "|");
    }
    return $access;
}

#ALL CLEAR
while(my $line = <$MODREADER>){
    
    chomp $line;
    my @terms = split(/\t+/,$line);
    if(!(scalar @terms < 5)){
	my @names = split(/,/,$terms[2]); 
    
	my @refined_names;
	foreach my $name (@names){
	    push @refined_names, (pointless $name);
	}
	
	push(@{$module_to_loc{$terms[4]}}, $terms[0]);
	push(@{$loc_to_accession{$terms[0]}}, @refined_names);
	push(@{$loc_to_info{$terms[0]}}, $terms[1]);
	push(@{$loc_to_info{$terms[0]}}, $terms[2]);
	push(@{$loc_to_info{$terms[0]}}, $terms[3]);
    }
}

close $MODREADER;
my %accession_to_go;

sysopen(my $GOREADER,$GO_Data,O_RDONLY);
#All Clear
while(my $line = <$GOREADER>){    
    chomp $line;
    my @terms = split(/\s+/,$line);
    my $gene;
    if ($terms[0] =~ /ref\|([^\.]+).{0,3}\|/){
	$gene = $1;
    }
    else{
	next;
    }
    my $go = $terms[1];
   
    push @{$accession_to_go{$gene}}, $go;
}
close $GOREADER;

#All clear
sysopen (my $GOFREQ, "$out.module_goterms.txt", O_WRONLY|O_CREAT);
foreach my $key (keys %module_to_loc){
    my %go_to_count;
    print $GOFREQ "$key\n";
    foreach my $gene (@{$module_to_loc{$key}}) {
	foreach my $accession (@{$loc_to_accession{$gene}}){
	    foreach my $goterm (@{$accession_to_go{$accession}}){
	    
		if(not (exists $go_to_count{$goterm})){
		    $go_to_count{$goterm} = 0; 
		}
		$go_to_count{$goterm} += 1;
	    }
	}
    }
    my @keys = sort{$go_to_count{$b} <=> $go_to_count{$a} } (keys (%go_to_count));
    for (my $i = 0; $i < 10 && $i < $#keys; $i++){
        print $GOFREQ "$keys[$i] - $go_to_count{$keys[$i]}\n";
    }   
}
close $GOFREQ;

#UNDER SCRUTINY
sysopen(my $OUTPUT, "$out.accession_to_go.txt", O_WRONLY|O_CREAT);
#parses each module
foreach my $module (keys %module_to_loc){
    #parses each loc in each module
    foreach my $loc (@{$module_to_loc{$module}}){
	foreach my $accession (@{$loc_to_accession{$loc}}){
	    my $go = join(';',uniq @{$accession_to_go{$accession}});
	    print $OUTPUT "$loc\t$accession\t$module\t$go\n";
	}
    }
}
close $OUTPUT;

#ALL CLEAR
sysopen($OUTPUT, "$out.loc_to_go.txt", O_WRONLY|O_CREAT);
foreach my $module (keys %module_to_loc){
    foreach my $loc (@{$module_to_loc{$module}}){
	my @go_terms;
	foreach my $accession (@{$loc_to_accession{$loc}}){
	    push @go_terms, (@{$accession_to_go{$accession}});
	}
	
	my $go = join(';', (uniq @go_terms));
	if(not $go){
	    $go = "-";
	}

	my $loc_info = join("\t",@{$loc_to_info{$loc}});
	print $OUTPUT "$loc\t$module\t$go\t$loc_info\n";
    }
}
close $OUTPUT;


sysopen(my $RESTRICTIONS,"./go_isolation/library.chosen_go_terms.txt",O_RDONLY);
sysopen(my $ORTHOLOGS,"./go_isolation/library.identified_suspect_TFs.txt", O_RDONLY);
sysopen($OUTPUT, "$ARGV[1].filtered.ACCEPTANCE_INFO.txt", O_WRONLY| O_CREAT);
sysopen(my $OUTPUT2, "$ARGV[1].filtered.reference.txt", O_WRONLY| O_CREAT);

my %weighted_go_terms;

my %exceptions = (
    
    "LOC100694761" => "LWS",
    "LOC100695018" => "SWS2B",
    "LOC100695287" => "SWS2A",
    "LOC100710676" => "RH2A-beta",
    "LOC100710942" => "RH2A-alpha",
    "LOC100711209" => "RH2B",
    "LOC100710249" => "SWS1",
    "LOC100711443" => "rho",
    "LOC100702555" => "Tbx2a",
    "LOC100705124" => "Rx1",
);

while(my $line = <$ORTHOLOGS>){

    chomp($line);
    $exceptions{$line} = $line;
    
}
while(my $line = <$RESTRICTIONS>){

    my @elements = split(/\t/,$line);
    $weighted_go_terms{$elements[0]} = $elements[1];
    
}



foreach my $module (keys %module_to_loc){
    foreach my $loc (@{$module_to_loc{$module}}){

	my $success = 0;
	my @vals;
	
	foreach my $go (keys %weighted_go_terms){

	    foreach my $accession (@{$loc_to_accession{$loc}}){
		if(grep { $_ eq "$go"} @{$accession_to_go{$accession}}){
		    $success += $weighted_go_terms{$go};
		    push @vals, $go;
		}
	    }
	    
	}

	if($success >= 5 || $exceptions{$loc}){

	    my $accept_info = ($exceptions{$loc} ? $exceptions{$loc} : join(",",@vals));
            my $gene_info = join("\t",@{$loc_to_info{$loc}});
	    my $xms = join(",",@{$loc_to_accession{$loc}});
	    
	    print $OUTPUT "$loc\t$accept_info\t$module\n";
	    
	    print $OUTPUT2 "$loc\t$gene_info\t$module\n";
		
	}
	
    }
}

close $OUTPUT;
close $OUTPUT2;
close $RESTRICTIONS;

#my %count;
#alternative means of filtering on GOTerm rules rather than GOTerm significance
#while(my $line = <$RESTRICTIONS>){
#    chomp $line;
#    my @requirements = split(/;/,$line);
#    print $OUTPUT ">$line\n";
#    
#    foreach my $module (keys %module_to_loc){
#	foreach my $loc (@{$module_to_loc{$module}}){
#	    my $success = 1;
#	    foreach my $term (@requirements){
#		my $has_term = 0;
#		foreach my $accession (@{$loc_to_accession{$loc}}){
#		    if(grep { $_ eq "$term"} @{$accession_to_go{$accession}}){
#			$has_term = 1;
#		    }
#		}
#		if(not $has_term){
#		    $success = 0;
#		    last;	    
#		}
#	    }
#	    if($success){ 
#		push @{$count{$loc}}, $line;
#		print $OUTPUT "$loc\t$module\n";
#	    }
#	}
#   }
#}

#close $OUTPUT;
#sysopen($OUTPUT, "$ARGV[1].filtered.threshold_limited.txt", O_WRONLY| O_CREAT);

#foreach my $module (keys %module_to_loc){
#    foreach my $loc (@{$module_to_loc{$module}}){
#	if(exists $count{$loc}){
#	    
#	    my @vals = @{$count{$loc}};
#	    my $size = @vals;
#	    if($size >= 3){
#		my $info = join(",",@vals);
#		my $xms = join(",",@{$loc_to_accession{$loc}});
#		print $OUTPUT "$loc\t@{$loc_to_info{$loc}}[0]\t$xms\t$info\t$module\n"
#	    }
#	}
#   }
#}

#close $RESTRICTIONS;
#close $OUTPUT;

