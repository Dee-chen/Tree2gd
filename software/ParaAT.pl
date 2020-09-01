#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

##################################################
my $Version     = "2.0";
my $ReleaseDate = "Oct. 4, 2014";
my $Author      = "Zhang Zhang";
my $Email       = "zhangzhang\@big.ac.cn";
my $Function    = "Parallel Alignment & Translation";
my $Reference   = "Zhang Z., et al. (2012) ParaAT: A parallel tool for constructing multiple protein-coding DNA alignments, Biochem Biophys Res Commun, 419, 779-781.";
my $WebSite     = "http://cbb.big.ac.cn";
##################################################

#***************************************************************
#                      The Genetic Codes
#   http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#***************************************************************
my %Genetic_Code = ();
$Genetic_Code{"1"}  = "The Standard Code";
$Genetic_Code{"2"}  = "The Vertebrate Mitochondrial Code";
$Genetic_Code{"3"}  = "The Yeast Mitochondrial Code";
$Genetic_Code{"4"}  = "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code";
$Genetic_Code{"5"}  = "The Invertebrate Mitochondrial Code";
$Genetic_Code{"6"}  = "The Ciliate, Dasycladacean and Hexamita Nuclear Code";
$Genetic_Code{"7"}  = "";
$Genetic_Code{"8"}  = "";
$Genetic_Code{"9"}  = "The Echinoderm and Flatworm Mitochondrial Code";
$Genetic_Code{"10"} = "The Euplotid Nuclear Code";
$Genetic_Code{"11"} = "The Bacterial, Archaeal and Plant Plastid Code";
$Genetic_Code{"12"} = "The Alternative Yeast Nuclear Code";
$Genetic_Code{"13"} = "The Ascidian Mitochondrial Code";
$Genetic_Code{"14"} = "The Alternative Flatworm Mitochondrial Code";
$Genetic_Code{"15"} = "Blepharisma Nuclear Code";
$Genetic_Code{"16"} = "Chlorophycean Mitochondrial Code";
$Genetic_Code{"17"} = "";
$Genetic_Code{"18"} = "";
$Genetic_Code{"19"} = "";
$Genetic_Code{"20"} = "";
$Genetic_Code{"21"} = "Trematode Mitochondrial Code";
$Genetic_Code{"22"} = "Scenedesmus Obliquus Mitochondrial Code";
$Genetic_Code{"23"} = "Thraustochytrium Mitochondrial Code";

#64 Codons and 20 Amino Acids in All Genetic Codes
my %Codon64 = ();
$Codon64{"AAA"} = "KKKKKK!!NKKKKNKK!!!!NKK";
$Codon64{"AAG"} = "KKKKKK!!KKKKKKKK!!!!KKK";
$Codon64{"AAT"} = "NNNNNN!!NNNNNNNN!!!!NNN";
$Codon64{"AAC"} = "NNNNNN!!NNNNNNNN!!!!NNN";
$Codon64{"TAA"} = "*****Q!!*****Y**!!!!***";
$Codon64{"TAG"} = "*****Q!!******QL!!!!*L*";
$Codon64{"TAT"} = "YYYYYY!!YYYYYYYY!!!!YYY";
$Codon64{"TAC"} = "YYYYYY!!YYYYYYYY!!!!YYY";
$Codon64{"ATA"} = "IMMIMI!!IIIIMIII!!!!MII";
$Codon64{"ATG"} = "MMMMMM!!MMMMMMMM!!!!MMM";
$Codon64{"ATT"} = "IIIIII!!IIIIIIII!!!!III";
$Codon64{"ATC"} = "IIIIII!!IIIIIIII!!!!III";
$Codon64{"TTA"} = "LLLLLL!!LLLLLLLL!!!!LL*";
$Codon64{"TTG"} = "LLLLLL!!LLLLLLLL!!!!LLL";
$Codon64{"TTT"} = "FFFFFF!!FFFFFFFF!!!!FFF";
$Codon64{"TTC"} = "FFFFFF!!FFFFFFFF!!!!FFF";
$Codon64{"GAA"} = "EEEEEE!!EEEEEEEE!!!!EEE";
$Codon64{"GAG"} = "EEEEEE!!EEEEEEEE!!!!EEE";
$Codon64{"GAT"} = "DDDDDD!!DDDDDDDD!!!!DDD";
$Codon64{"GAC"} = "DDDDDD!!DDDDDDDD!!!!DDD";
$Codon64{"CAA"} = "QQQQQQ!!QQQQQQQQ!!!!QQQ";
$Codon64{"CAG"} = "QQQQQQ!!QQQQQQQQ!!!!QQQ";
$Codon64{"CAT"} = "HHHHHH!!HHHHHHHH!!!!HHH";
$Codon64{"CAC"} = "HHHHHH!!HHHHHHHH!!!!HHH";
$Codon64{"GTA"} = "VVVVVV!!VVVVVVVV!!!!VVV";
$Codon64{"GTG"} = "VVVVVV!!VVVVVVVV!!!!VVV";
$Codon64{"GTT"} = "VVVVVV!!VVVVVVVV!!!!VVV";
$Codon64{"GTC"} = "VVVVVV!!VVVVVVVV!!!!VVV";
$Codon64{"CTA"} = "LLTLLL!!LLLLLLLL!!!!LLL";
$Codon64{"CTG"} = "LLTLLL!!LLLSLLLL!!!!LLL";
$Codon64{"CTT"} = "LLTLLL!!LLLLLLLL!!!!LLL";
$Codon64{"CTC"} = "LLTLLL!!LLLLLLLL!!!!LLL";
$Codon64{"AGA"} = "R*RRSR!!SRRRGSRR!!!!SRR";
$Codon64{"AGG"} = "R*RRSR!!SRRRGSRR!!!!SRR";
$Codon64{"AGT"} = "SSSSSS!!SSSSSSSS!!!!SSS";
$Codon64{"AGC"} = "SSSSSS!!SSSSSSSS!!!!SSS";
$Codon64{"TGA"} = "*WWWW*!!WC**WW**!!!!W**";
$Codon64{"TGG"} = "WWWWWW!!WWWWWWWW!!!!WWW";
$Codon64{"TGT"} = "CCCCCC!!CCCCCCCC!!!!CCC";
$Codon64{"TGC"} = "CCCCCC!!CCCCCCCC!!!!CCC";
$Codon64{"ACA"} = "TTTTTT!!TTTTTTTT!!!!TTT";
$Codon64{"ACG"} = "TTTTTT!!TTTTTTTT!!!!TTT";
$Codon64{"ACT"} = "TTTTTT!!TTTTTTTT!!!!TTT";
$Codon64{"ACC"} = "TTTTTT!!TTTTTTTT!!!!TTT";
$Codon64{"TCA"} = "SSSSSS!!SSSSSSSS!!!!S*S";
$Codon64{"TCG"} = "SSSSSS!!SSSSSSSS!!!!SSS";
$Codon64{"TCT"} = "SSSSSS!!SSSSSSSS!!!!SSS";
$Codon64{"TCC"} = "SSSSSS!!SSSSSSSS!!!!SSS";
$Codon64{"GGA"} = "GGGGGG!!GGGGGGGG!!!!GGG";
$Codon64{"GGG"} = "GGGGGG!!GGGGGGGG!!!!GGG";
$Codon64{"GGT"} = "GGGGGG!!GGGGGGGG!!!!GGG";
$Codon64{"GGC"} = "GGGGGG!!GGGGGGGG!!!!GGG";
$Codon64{"CGA"} = "RRRRRR!!RRRRRRRR!!!!RRR";
$Codon64{"CGG"} = "RRRRRR!!RRRRRRRR!!!!RRR";
$Codon64{"CGT"} = "RRRRRR!!RRRRRRRR!!!!RRR";
$Codon64{"CGC"} = "RRRRRR!!RRRRRRRR!!!!RRR";
$Codon64{"GCA"} = "AAAAAA!!AAAAAAAA!!!!AAA";
$Codon64{"GCG"} = "AAAAAA!!AAAAAAAA!!!!AAA";
$Codon64{"GCT"} = "AAAAAA!!AAAAAAAA!!!!AAA";
$Codon64{"GCC"} = "AAAAAA!!AAAAAAAA!!!!AAA";
$Codon64{"CCA"} = "PPPPPP!!PPPPPPPP!!!!PPP";
$Codon64{"CCG"} = "PPPPPP!!PPPPPPPP!!!!PPP";
$Codon64{"CCT"} = "PPPPPP!!PPPPPPPP!!!!PPP";
$Codon64{"CCC"} = "PPPPPP!!PPPPPPPP!!!!PPP";

#*********** Begin of Sub-functions Declaration **********#
sub AlignProteinSeqs;
sub AppendStr2Array;
sub BriefDesc;
sub CheckCMD;
sub ConvertFasta2AXT;
sub DelFile;
sub DelFiles;
sub ExecProgram;
sub GenerateSequenceFiles;
sub KaKsCommands;
sub PA2CACommands;
sub ProteinAlignmentCommands;
sub ReadFile;
sub ReadHomologList;
sub ReadSeqs;
sub ReadThreadNumber;
sub RunMultiProcesses;
sub Usage;
#************ End of Sub-functions Declaration ***********#

#Multiple sequence aligners for amino acid sequences
my @MSA = ("clustalw2", "t_coffee", "mafft", "muscle");
#KaKs_Calculator for estimating Ka and Ks
my $KaKsCMD = "KaKs_Calculator";
my $m="YN";
my $alignCMD="mafft";
#Epal2nal
my $Pal2NalCMD = "Epal2nal.pl";

#input parameters
my $HomologFile     =  "";
my $AAFile          =  "";
my $NUCFile         =  "";
my $MSACMD          =  "";
my $Process_num     =  "1";
my $OutputFolder    =  "";
my $OutputFormat    =  "fasta";
my $Code            =  "1";

my $KaKs_Calculator;
my $RemoveMismatch ;
my $RemoveGap      ;
my $Verbose        ;
my $Help           ;

#Parse input parameters
GetOptions(
    'h|homolog=s'    => \$HomologFile,
    'a|aminoacid=s'  => \$AAFile,
    'n|nuc=s'        => \$NUCFile,
    'p|processor=s'  => \$Process_num,
    'o|output=s'     => \$OutputFolder,
    'm|msa:s'        => \$MSACMD,
    'f|format:s'     => \$OutputFormat,
    'c|code:s'       => \$Code,
    'g|nogap!'       => \$RemoveGap,
    't|nomismatch!'  => \$RemoveMismatch,
    'k!'        => \$KaKs_Calculator,
    'v|verbose!'     => \$Verbose,
		'agcmd:s'        => \$alignCMD,
    'kaks|KaKs_Calculator:s'  => \$KaKsCMD,
    'p2n:s'  =>\$Pal2NalCMD,
    'mod|KaKs_Calculator_mod:s'=> \$m,
    'help|?!'        => \$Help,
) or die "Incorrect usage!\n\n" && Usage();

#Help or no parameter
if ($Help) {
    Usage();
}

BriefDesc();

print "Checking input parameters: ";
if ( $HomologFile eq "" || !(-e $HomologFile) ) {
    print "Required parameter for homologous group file or the file does not exist.\n";
    print "For help information: $0 -help.\n";
    exit(0);
}

if ( $AAFile eq "" || !(-e $AAFile) ) {
    print "Required parameter for amino acid sequence file or the file does not exist.\n";
    print "For help information: $0 -help.\n";
    exit(0);
}

if ( $NUCFile eq "" || !(-e $NUCFile) ) {
    print "Required parameter for nucleotide sequence file or the file does not exist.\n";
    print "For help information: $0 -help.\n";
    exit(0);
}

if ( $Process_num eq "0" ) {
    print "Required parameter for processor file or the file does not exist.\n";
    print "For help information: $0 -help.\n";
    exit(0);
}

if ($OutputFolder eq "") {
    print "Required parameter for output folder.\n";
    print "For help information: $0 -help.\n";
    exit(0);
}

if ($Genetic_Code{"$Code"} eq "") {
    print "Invalid parameter for the selected genetic code.\n";
    print "For help information: $0 -help.\n";
    exit(0);
}
print "[OK]\n";

print "Checking Epal2nal.pl: ";
if (CheckCMD($Pal2NalCMD)==1) {
    print "[OK]\n";
} else {
    print "[NA]\n";
    print "Please make sure Epal2nal.pl is executable (+X) and put into the global directoryy.\n";
    exit(0);
}

print "Checking multiple sequence aligner: ";
if ($MSACMD eq "") { #not specify by user
    for (my $i=0; $i<@MSA; $i++) {
        print "$MSA[$i]";
        if (CheckCMD($MSA[$i])==1) {
            print "[OK] ";
            if ($MSACMD eq "") {
                $MSACMD = $MSA[$i];
            }
        }
        else {
            print "[NA] ";
        }
    }
    print "\n";
    if ($MSACMD eq "") {
        print "No multiple sequence aligner is available. Please install any supported aligner first (".join("|",@MSA).").\n";
        print "For help information: $0 -help.\n";
        exit(0);
    }
}
else {
    print "$MSACMD";
    if ( CheckCMD($alignCMD) == 0 ) {
        print "[NA]\n$alignCMD is not accessible. Please make sure it is put into the global directory or specify its full path.\n";
        print "For help information: $0 -help.\n";
        exit(0);
    }
    else {
        print "[OK]\n";
    }
}

print "Checking output format: ";
if ( !($OutputFormat eq "axt" ||
        $OutputFormat eq "clustal" ||
        $OutputFormat eq "paml" ||
        $OutputFormat eq "fasta" ||
        $OutputFormat eq "codon") ) {
	print "The output format should be ‘axt‘, ‘fasta’, ‘paml’, ‘codon’, or ‘clustal’\n";
	print "For help information: $0 -help.\n";
	exit(0);
}
print "[OK]\n";

if ( $KaKs_Calculator) {
    print "Checking KaKs_Calculator: ";
    if ($OutputFormat ne "axt") {
	    print "Error in parameters’ setting. To calculate Ka and Ks by KaKs_Calculator, the output format should be axt.\n";
	    print "For help information: $0 -help.\n";
	    exit(0);
    }
    if (CheckCMD($KaKsCMD)==0) {
	    print "$KaKsCMD is not accessible. Please make sure it is put into the global directory or specify its full path.\n";
	    print "For help information: $0 -help.\n";
	    exit(0);
    }
    print "[OK]\n";
}

print "Checking output folder: ";
if (-e $OutputFolder) {
	print "Error in creating the folder $OutputFolder, since it already exists.\n";
	print "For help information: $0 -help.\n";
	exit(0);
}
else {
	mkdir($OutputFolder);
}
print "[OK]\n";

print "\nOutputing input parameters:\n";
print "\t Homologous groups = $HomologFile\n";
print "\t Amino acids sequences = $AAFile\n";
print "\t Nucleotide sequences = $NUCFile\n";
print "\t Process num = $Process_num\n";
print "\t Output folder = $OutputFolder\n";
print "\t Multiple sequence aligner = $MSACMD\n";
print "\t Output format = $OutputFormat\n";
print "\t Genetic code = $Code-$Genetic_Code{$Code} \n";
print "\t KaKs_Calculator used = ".($KaKs_Calculator?"TRUE":"FALSE")."\n";
print "\t Remove gap = ".($RemoveGap?"TRUE":"FALSE")."\n";
print "\t Remove mismatched codons = ".($RemoveMismatch?"TRUE":"FALSE")."\n";
print "\t Verbose output = ".($Verbose?"TRUE":"FALSE")."\n";
print "\n";

#Read Homolog List
print "Reading homologous groups from $HomologFile: ";
my @Homologs = ReadHomologList($HomologFile);
my $NumberOfHomologs = @Homologs;
print "$NumberOfHomologs groups\n";

#Read nucleotide sequences
print "Reading nucleotide sequences from $NUCFile: ";
my %hash_nuc = ReadSeqs($NUCFile);
my $NumberOfNUCSeqs = scalar(keys %hash_nuc);
print "$NumberOfNUCSeqs nucleotide sequences\n";

#Read amino acid sequences
print "Reading amino acid sequences from $AAFile: ";
my %hash_aa  = ReadSeqs($AAFile);
my $NumberOfAASeqs = scalar(keys %hash_aa);
print "$NumberOfAASeqs amino acid sequences\n";

#Check consistency of sequences


my $start_time = time();

#Change current working directory
chdir($OutputFolder);

#Generate sequence files
print "Generating homologous group files for alignments: ";
my @SeqFiles = GenerateSequenceFiles(\@Homologs, \%hash_aa, \%hash_nuc);
my $NumberOfFiles = @SeqFiles;
print "$NumberOfFiles from $NumberOfHomologs groups\n";

#Align Multiple Protein Sequences
my @Commands = ProteinAlignmentCommands(\@SeqFiles);
my $NumberOfCommands = @Commands;
print "Aligning amino acid sequences for $NumberOfCommands homologous groups:\n";
RunMultiProcesses(\@Commands);

#Translate Protein Alignments into Codon Alignments
@Commands = PA2CACommands(\@SeqFiles);
$NumberOfCommands = @Commands;
print "\nTranslating amino acid alignments into codon alignments:\n";
RunMultiProcesses(\@Commands);

my $end_time = time();

#Calculate Ka and Ks, if specifying -k or -kaks
if ($KaKs_Calculator) {
    @Commands = KaKsCommands(\@SeqFiles);
    $NumberOfCommands = @Commands;
    print "\nCalculating Ka and Ks by KaKs_Calculator:\n";
    RunMultiProcesses(\@Commands);
}

if ($Verbose) {
}
else {
    print "\nCleaning temporary files:\n";
    my @tmpFiles;
    @tmpFiles = AppendStr2Array(\@SeqFiles, ".pep");
    DelFiles(\@tmpFiles);
    @tmpFiles = AppendStr2Array(\@SeqFiles, ".cds");
    DelFiles(\@tmpFiles);
    @tmpFiles = AppendStr2Array(\@SeqFiles, ".dnd");
    DelFiles(\@tmpFiles);
    @tmpFiles = AppendStr2Array(\@SeqFiles, ".pep_aln");
    DelFiles(\@tmpFiles);

    DelFile("msg.msa");
    DelFile("msg.kaks");
}

print "\nMission Accomplished (Time used: ", $end_time - $start_time, "s).\n\n";

######## End of Main #########


#################### Begin of Sub functions #################

#Append a string to all elements of a array
sub AppendStr2Array {
	my ($array, $str) = @_;

	my @files = @$array;
	my $i = 0;
	for($i=0; $i<@files; $i++) {
	    $files[$i] .= $str;
	}

	return @files;
}

#Convert fasta file into axt file
sub ConvertFasta2AXT {
	my ($file) = @_;

	my @array = ReadFile($file);
	my $count = 0;
	my @seq = ();
	my @name = ();
	my $tmp_seq = "";
	while ($count < @array) {
		if(index($array[$count], ">")==0) {
			push(@name, substr($array[$count], 1, length($array[$count])-1));
	   	if ($tmp_seq ne "") {
				push(@seq, $tmp_seq);
				$tmp_seq = "";
			}
		}
		else {
			$tmp_seq .= $array[$count];
		}
		$count++;

		if ($count==@array) {
	  	push(@seq, $tmp_seq);
	  }
	}

	my $output  = join("-", @name)."\n";
	$output .= join("\n", @seq)."\n";
	$output .= "\n";

	my $outfile = $file.".axt";
	open(OUTFILE, ">$outfile");
	print OUTFILE ($output);

	return;
}

#*** Check the availability of command ***#
sub CheckCMD {
	my ($cmd) = @_;

	my $Exist = 0;
	if (-x $cmd) {
        $Exist = 1;
    }
    else {
        my $i = 0;
        #PATH variable
        my $string_path = $ENV{"PATH"};
        #Replace ~ as HOME variable
        $string_path =~ s/~/$ENV{"HOME"}/;
        #Obtain multiple global directories
        my @PATH = split(/:/, $string_path);
        chomp @PATH;
        #delete the last character ending with "/"
        for ($i=0; $i<@PATH; $i++) {
            $PATH[$i] =~ s/\/$//;
        }
        #Push current directory
        push(@PATH, $ENV{'PWD'});

        #Judge which directory the cmd can be accessed
        for ($i=0; $i<@PATH && $Exist==0; ) {
            #Found executive file
            if (-x $PATH[$i]."/".$cmd) {
                $Exist = 1;
            }
            else {
                $i++;
            }
        } #End of for
	} #End of else

	return $Exist;
}

#*** Run multiple processes ***#
sub RunMultiProcesses {
	my ($array) = @_;

	my @Commands = @$array;
	my $totalCmds = @Commands;
	#Display the process at every loop
	my $loop = (int($totalCmds / 100))==0?1:int($totalCmds / 100);

	my $delimit = "   ";
	#Init process count
	my $maxnum = 1;
	#num is the number of process defined in the file, initialized as maxnum
	my $num = $maxnum;
	my $CmdCount = 0;

	while (@Commands) {
        my $cmd = shift(@Commands);
        #init a new sub process,and CmdCount++
        #my $pid = fork();
        #if($pid == 0){
        unless (fork()){
        	ExecProgram($cmd);
        }
        $CmdCount++;

        #Read the current number of processors, which can be changeable at running time
        $num =$Process_num;

        #Compare the current number of processors with that obtained from the configration file
        if ($num > $maxnum) {
                my $count = 0;
                for(my $i = 0; $i < ($num-$maxnum); $i++) {
                    if(my $newcmd = shift(@Commands)) {
                  	    unless (fork()) {
                  		    ExecProgram($newcmd);
                  	    }
                  	    $count++;
                  	    $CmdCount++;
                    }
                    else {
                        last;
                    }
                }
                $maxnum += $count;
                print $delimit."Number of processors = $maxnum\n";
        }
        elsif ($num<$maxnum) {
            if($num == 0) { #suspend mode
                print $delimit."Enter suspend mode! Please wait the last $maxnum processors\n";
            }
            for (my $i=0; $i<($maxnum-$num);$i++) {
                wait();
            }

            $maxnum = $num;

            if($num == 0) { #suspend mode
                if(@Commands){
                    my $nextcmd = (shift @Commands);
                    print $delimit."Process Suspend!\nNext command will be: $nextcmd \n";
                    print $delimit."Total $CmdCount commands done\n";
                }
                else {
                    print $delimit."Process Suspend!\nBut all commands are done! :)\n";
                }
                exit(0);
            }
            print $delimit."\nCurrent processor number: $maxnum\n";
        }

        wait();

        my $percent = int(10000*$CmdCount/$totalCmds)/100;
        my $now = localtime(time());
        if ($CmdCount % $loop == 0 || $CmdCount == $totalCmds) {
            print $delimit."$CmdCount ($percent%) commands done at $now\n";
        }
	}

	for (my $i=0; $i<$maxnum-1; $i++) {
		wait();
	}
}


#*** Display the program's brief descritption ***#
sub BriefDesc {
	print "****************************************************************\n";
	print "Program: $0 [Version: $Version; $ReleaseDate] \n";
	print "Function: $Function \n";
	print "Reference: $Reference \n";
	print "****************************************************************\n";
}

#*** Show help information ***#
sub Usage{

    BriefDesc();

    print "USAGE: $0 <Options>\n";

    print "\t-h, -homolog\tHomolog group file [string, required]\n";
    print "\t-a, -aminoacid\tFile containing multiple amino acid sequences [string, required]\n";
    print "\t-n, -nuc\tFile containing multiple nucleotide sequences [string, required]\n";
    print "\t-p, -processor\tFile containing the number of processors (< number of physical cores), changeable during running [string, required]\n";
    print "\t-o, -output\tOutput folder [string, required]\n";

    print "\t-m, -msa\tMultiple sequence aligner or specify with its full path [" . join(" | ",@MSA) . ", optional]\n";
    print "\t-f, -format\tOutput file format [axt | fasta | paml | codon | clustal, optional], default = fasta\n";
    print "\t-c, -code\tGenetic Code used [integer, optional], default = 1-The Standard Code\n";
    my $code;
    foreach $code (sort{$a<=>$b} keys %Genetic_Code) {
        my $value = $Genetic_Code{$code};
        if ($value ne "") {
            print "\t\t\t$code-$value\n";
        }
    }
    print "\t\t\tSee all documented codes at <http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>\n";

    print "\t-g, -nogap\tRemove aligned codons with gaps\n";
    print "\t-t, -nomismatch\tRemove mismatched codons\n";
    print "\t-k, -kaks\tEnable using KaKs_Calculator for Ka and Ks estimation (requiring axt format)\n";
    print "\t-v, -verbose\tVerbose output\n";

    print "\t-?, -help\tDisplay help information\n";
    print "\n";

    print "NOTE:"."\n";
	print "\tSequence IDs in homologous group file, amino acid file and nucleotide file should be the same.\n";
	print "\n";

	print "DOWNLOAD:"."\n";
	print "\tPlease visit <$WebSite> for downloads and updates."."\n";
	print "\n";

	print "CONTACT:"."\n";
	print "\tPlease send suggestions or bugs to $Author at $Email."."\n";
	print "\n";

    exit(1);
}

#*** Read homologous relationship list ***#
sub ReadHomologList(){
    my ($file) = @_;
    return ReadFile($file);
}

#*** Read fasta sequences ***#
sub ReadSeqs(){
    my ($file) = @_;

	my @array = ReadFile($file);
	my %hash = ();

	while (@array) {
		my $name = shift(@array);
		if(index($name, ">") > -1) {
			#remove >
			$name =~ s/>//g;
			#remove blank character
			$name =~ s/\s//g;
			#remove version, e.g., NP_049009.1 to NP_049009
			#$name =~ s/\.[0-9]+$//g;
			#$name = substr($name, 0, index($name, "."));

			my $seq = "";
			my $tmp = shift(@array);
			while (@array && index($tmp, ">")<0) {
				$seq .= $tmp;
				$tmp = shift(@array);
			}
			if (@array) {
				unshift(@array, $tmp);
			}
			else {
				$seq .= $tmp;
			}

			$seq =~ s/\s//g;
			$hash{$name} = $seq;
		}
	}
	@array = ();

	return %hash;

}

#*** Read file text into array ***#
sub ReadFile() {
	my ($file) = @_;

	unless (open (MYFILE, $file)) {
        die ("Error in opening file: $file\n");
	}
	my @array = <MYFILE>;
	chomp @array;
	close(MYFILE);

	return @array;
}

#*** Read the number of thread from file ***#
sub ReadThreadNumber(){
	my ($current, $con_file) = @_;

	my $threads = $current;
	$con_file = "../".$con_file;

	my @config = &ReadFile($con_file);
	if (@config) {
		$threads = $config[0];
	  $threads =~ s/\s+//g;
	  if($threads < 0 || $threads!~/^\d+$/){
	  	warn("The nubmer of threads given in $con_file is invalid: $threads!\n");
	  }
	}

	return $threads;
}

#*** Delete single file ***#
sub DelFile(){
	my ($file) = @_;
	unlink($file);
}

#*** Delete multiple files stored in a array ***#
sub DelFiles(){
	my ($files) = @_;
	for (my $i=0; $i<@$files; $i++) {
		unlink(@$files[$i]);
	}
}

#*** Generate multiple sequence files by splitting ***#
sub GenerateSequenceFiles(){
	my ($array, $aa_hash, $nuc_hash) = @_;

	my @homologs = @$array;
	my %aa_seqs  = %$aa_hash;
	my %nuc_seqs = %$nuc_hash;

	my @files = ();
	for (my $i=0; $i<@homologs; $i++) {
		my @genes = split(/\s+/, $homologs[$i]);
		chomp @genes;

		my $file = join("-", @genes);
		my $aa_seq = "";
		my $nuc_seq = "";
		my $flag = 1;
		for (my $j=0; $j<@genes && $flag; $j++) {
		    #remove version, e.g., NP_049009.1 to NP_049009
			#$genes[$j] =~ s/\.[0-9]+$//g;
			if (defined $aa_seqs{$genes[$j]} && defined $nuc_seqs{$genes[$j]}) {
				$aa_seq .= ">".$genes[$j]."\n";
				$aa_seq .= $aa_seqs{$genes[$j]}."\n";
				#delete ( $aa_seqs{$genes[$j]} );

				$nuc_seq .= ">".$genes[$j]."\n";
				$nuc_seq .= $nuc_seqs{$genes[$j]}."\n";
				#delete ( $nuc_seqs{$genes[$j]} );
			}
			else {
				$flag = 0;
			}
		}
		if ($flag) {
			open(OUTFILE, ">$file.pep");
	  		print OUTFILE ("$aa_seq");

	  		open(OUTFILE, ">$file.cds");
	  		print OUTFILE ("$nuc_seq");

	  		push(@files, $file);
	  	}
	  	else {
	  	    #warn "\tThe sequence(s) is not available for homologous group ‘".join(" - ",@genes)."’.\n";
	  	}
	}
	chomp @files;

	%aa_seqs  = ();
	%nuc_seqs = ();

	return @files;
}

#*** Generate commands for protein alignments ***#
sub ProteinAlignmentCommands(){
	my ($array) = @_;
	my @files = @$array;

	my @cmds = ();
	for (my $i=0; $i<@files; $i++) {
	  push(@cmds, AlignProteinSeqs($files[$i]));
	}
	chomp @cmds;
	return @cmds;
}

#*** Generate a command for protein alignment ***#
sub AlignProteinSeqs(){
	my ($file) = @_;
	my $cmd = "";

	if ($MSACMD =~ /clustalw2/) {
	    $cmd = $alignCMD." -INFILE=".$file.".pep  -PWMATRIX=BLOSUM -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE=".$file.".pep_aln >> msg.msa";
	}
	elsif ($MSACMD =~ /mafft/) {
	    $cmd = 	$alignCMD."  --auto --quiet  ".$file.".pep  >  ".$file.".pep_aln";
	}
	elsif ($MSACMD =~ /muscle/) {
	    $cmd = $alignCMD." -quiet  -in ".$file.".pep  -out ".$file.".pep_aln >> msg.msa";
	}
	elsif ($MSACMD =~ /t_coffee/) { # not stable and greedy for memory, tested on MAC OSX
        $cmd = 	$alignCMD."  ".$file.".pep  -output=fasta  -outfile=".$file.".pep_aln -quiet >> msg.msa";
	}
	return $cmd;
}

#*** Generate commands for translating protein alignments into codon alignments ***#
sub PA2CACommands {
	my ($array) = @_;

	my @files = @$array;

	my $para = "";
	if ($RemoveGap) {
	    $para .=  " -nogap ";
	}
	if ($RemoveMismatch) {
	    $para .=  " -nomismatch ";
	}

	my @cmds = ();
	for (my $i=0; $i<@files; $i++) {
		my $file = $files[$i];
		my $cmd = $Pal2NalCMD."  ".$file.".pep_aln  ".$file.".cds  -output $OutputFormat  $para > ".$file.".cds_aln.".$OutputFormat;
		push(@cmds, $cmd);
	}
	chomp @cmds;
	return @cmds;
}

#*** Generate commands for calculating Ka and Ks ***#
sub KaKsCommands {
	my ($array) = @_;

	my @files = @$array;
	my @cmds = ();
	for (my $i=0; $i<@files; $i++) {
		my $file = $files[$i];
		my $cmd = $KaKsCMD." -m ".$m." -i ".$file.".cds_aln.axt -o  ".$file.".cds_aln.axt.kaks >> msg.kaks";
		push(@cmds, $cmd);
	}
	chomp @cmds;
	return @cmds;
}

#*** Execute a command ***#
sub ExecProgram(){
    my($cmd) = @_;
	exec($cmd);
}

#################### End of Sub functions #################
