#.....................License.................................
#	 MysiRNA-Designer a Software for Rational siRNA Design
#    Copyright (C) 2011  <M.Mysara et al>

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
	
#......................Packages Used..........................
use Getopt::Std;
use Cwd;
use strict;
use Bio::DB::GenBank;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::RemoteBlast;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Spreadsheet::ParseExcel;
use Win32::OLE qw(in with);
use Win32::OLE::Const 'Microsoft Excel';
use Modern::Perl;
use Bio::DB::EUtilities;
use Bio::SearchIO;
require LWP::UserAgent;
#.................adjusting the PATH to the current PATH...........................
my $dir = cwd;	use Cwd 'chdir';    chdir "/tmp";
#......................Creating directories ......................
mkdir "results";	
mkdir "RNAxsresult";	
mkdir "MysiRNA_Designer_Final_results";	system "del /Q /s MysiRNA_Designer_Final_results > temp.txt";
mkdir "blastresult";	
mkdir "Multi_score_Accepted";	system "del /Q /s Multi_score_Accepted > temp.txt";	
#...................Varriables...................
my $num=0;#wholemain
my  $trackernum=0;# fasta1
my %opts;
my $header;
my $database = new Bio::DB::GenBank;
my $trans_num=0;
#....................................
getopt('studowpei',\%opts);
if(defined($opts{s})&&$opts{o}&&$opts{w})
{
	my $off_target_filtration=1;# to perform off-target filtration or not as per used requirement
	if(!defined($opts{t})) {$opts{t}=0.01157; }     
	if(!defined($opts{u})) {$opts{u}=0.001002; }
	if(!defined($opts{p})) {$opts{p}="N"; } 	
	if(!defined($opts{e})) {$opts{e}="ENST_UTR.txt"; }  
	if(!defined($opts{i})) {$opts{i}="ENST.txt"; }
	if(!defined($opts{d})) {$off_target_filtration=0;$opts{d}="None"; }
	#............................................
	open Main, $opts{s};		
	while (my $acc=<Main>){
	$num++;
		print '
	||||||||||||||||||||||||||||||||||||||||||||||||
	||        Welcome To MysiRNA-Designer         ||
	||    A Software For Rational siRNA Design    ||
	||    Copyright (C) 2011  <M.Mysara et al>    ||
	||||||||||||||||||||||||||||||||||||||||||||||||
MysiRNA-Designer version 1, Copyright (C) 2011, M.Mysara et al
MysiRNA-Designer comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under certain conditions; type `notepad copying.txt\' for details.';
		#............Delete Buffer................
		system "del /Q /s results > temp.txt";	
		system "del /Q /s RNAxsresult > temp.txt";	
		system "del /Q /s blastresult > temp.txt";	
		system "del /Q /s Multi_score_Accepted > temp.txt";
		system "del /Q /s MysiRNA_Designer_Final_results > temp.txt";
		#..............................................................
		my $seq = $database->get_Seq_by_id($acc);
		print "\nSeq: ", $seq->accession_number(), " -- ", $seq->desc(), "\n\n";
		my $Desc;
		my $Desc_;
		my $habal=$seq->desc();
		if ($habal=~/(^[\w\W]*, transcript variant)/) {
			$Desc=$1; 
			
			}
			
		#make inputs for scoring (Iscore)
		
		my $out = Bio::SeqIO->newFh ( -file => ">".$ENV{'PWD'}."/results/fasta1result.fasta", -format => 'fasta');
		open FH,">",$ENV{'PWD'}.'\results\iscoreinput.fasta'; print FH "";close FH;
		open FH,">>",$ENV{'PWD'}.'\results\iscoreinput.fasta';
		print $out $seq;
		print FH $seq->seq;
		
		#Make SNPs Database search input
		
		my $out_ = Bio::SeqIO->newFh ( -file => ">".$ENV{'PWD'}."/results/fasta1result.gb", -format => 'GenBank');
		print $out_ $seq;
		close (FH);
		open FH,">",$ENV{'PWD'}.'\AlignIO_.fasta'; print FH "";close FH;
		open FH,">>",$ENV{'PWD'}.'\AlignIO_.fasta';
		print FH ">habal\n",$seq->seq;
		close (FH);
		print "fasta1 done \n";
		
		#Blast in case of Multitranscripts, and getting the sequence info for all transcripts
		open FILE, ">", $ENV{'PWD'}.'\results\blast1result.fasta';
		if ($seq->desc=~/transcript variant/){
		my $prog = "blastn";
		my $db = "nr";
		my $e_val = "1e-5";
		my $remoteBlast = Bio::Tools::Run::RemoteBlast->new(-prog => $prog,
									-data => $db,
									-expect => $e_val);
		# TO DEFINE SPECIFIC ORGANISM: $Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = 'Homo sapiens [ORGN]';
		my $r = $remoteBlast->submit_blast($ENV{'PWD'}.'\results\fasta1result.fasta');
		
		#Submit blast jobs to ncbi blast queue on sequence(s)
		
			while (my @reqIDs = $remoteBlast->each_rid ) {
				print STDERR join(" ", "\nINFO: RIDs: ", @reqIDs), "\n";
				foreach my $reqID (@reqIDs) {          
				
				# each search results
				
					my $rc = $remoteBlast->retrieve_blast($reqID);
					if (! ref ($rc)){
							if ($rc < 0) {                  # no match
								
								$remoteBlast->remove_rid($reqID);
							}
						
						# Search is not done yet, wait 10 sec, and try to retrieve again
						
						print STDERR ".";
						sleep (10);
					}
					
					# got some blast hit
				
					else {                            
						my $result = $rc->next_result;  # get the blast output
						while(my $hit = $result->next_hit) {
							$Desc_=0;
							my $haball=$hit->description();
							my $accc = $hit->accession();
							if (($haball=~/(^[\w\W]*, transcript variant)/)){
								$Desc_ =$1;}#&& $accc=~/[A-Z][A-Z]_/
								if(($Desc eq $Desc_)){
									# print out the accession etc of all hits
									print FILE $hit->accession,"\n"; #$hit->description, "\t", $hit->name, "\t ", $hit->description, "\t", $hit->significance, "\t", $hit->score, "\t", $hit->frac_identical, "\n";
								   $trackernum ++;
								}
							}
						print STDERR "\nINFO: removing $reqID\n";
						$remoteBlast->remove_rid($reqID);  # remove this RID since we 
														   # already  got the results  
						}
					}
				}
			}
		else{
			print FILE $seq->accession_number();
		}
		close FILE;
		#To get the number of transcripts
		my $file3=$ENV{'PWD'}.'\results\blast1result.fasta';
		open FH1,$file3;
		my @array=<FH1>;
		print @array,"\n";
		$trans_num=@array;
		print $trans_num,"\n";
		close FH1; 
		if ($trans_num>1){
			#Getting sequence fasta for all transcripts 
			open FH1, "<", $ENV{'PWD'}.'\results\blast1result.fasta';
			open FH2,">",$ENV{'PWD'}.'\results\fasta2result.fasta';
			print  FH2 "";
			my $kk;
			while (defined ($kk=<FH1>)){
			
				my $seq = $database->get_Seq_by_id($kk);
				print "Seq: ", $seq->accession_number(), " -- ", $seq->desc(), "\n\n";
				my $out = Bio::SeqIO->newFh ( -file => ">>".$ENV{'PWD'}.'\results\fasta2result.fasta', -format => 'fasta');
				print $out $seq;
				
			}	
			close (FH1);
			print "Multi varriant\n";
			my $file1=$ENV{'PWD'}.'\results\fasta2result.fasta';
			my$file2=$ENV{'PWD'}.'\results\clustalwresult';
			
			#Calling Local ClustalW (needs ClustalW to be installed in ur PC and varriables intialization in PATH and ClustPath)
			my $filename = $ENV{'PWD'}.'\results\clustalwresult.aln';
			
			#  Build a clustalw alignment factory
			my @params = ('output' => 'clustal','TYPE'=>'DNA','outfile'=>"$filename");
			my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
			
			#  Pass the factory a list of sequences to be aligned.	
			my $inputfilename = $ENV{'PWD'}.'\results\fasta2result.fasta';
			my $aln = $factory->align($inputfilename); # $aln is a SimpleAlign object.
			my $file4=$ENV{'PWD'}.'\results\clustalwresult.aln';
			my $file5=$ENV{'PWD'}.'\results\consresult.fasta';
			#use cons.exe tool
			system "cons.exe -sequence $file4 -outseq $file5 -identity $trans_num -setcase 0";
			#analysisng cons.exe results
			my $lines;
			open FH1,"<","$file5";
			while (my $line=<FH1>){
				if ($line=~/^>/){
				#Do Nothing "Header"
				}
				else{
					chomp $line;
					$lines=$lines.$line;
				}
			}
			close FH1;
			open FH1,">","$file5";
			print FH1 $lines;
			close FH1;
			
			open FH1,"<","$file5";
			open FH2,">",$ENV{'PWD'}.'\results\AlignIO.fasta';print FH2 ""; close FH2;
			open FH2,">>",$ENV{'PWD'}.'\results\AlignIO.fasta';
			while (my $seq=<FH1>){
				if ($seq=~/^>/){}
				else{
					for my $i(0..(length($seq)-19)){
						my $substring=substr($seq,$i,19);
						if ($substring=~/N/i){
						#if not 100% conserved Do notting
						}
						else {
							#if 100% conserved switch T => U
							$substring =~s/T/U/gi;print FH2 $substring,"\n" ;
						}
					}
				}
			}
			close FH1;
			close FH2;
			print "consresult.pl done!!!\n"; 			
		}
		else {
			print "one varriant\n";
			my $del=$ENV{'PWD'}.'\results\AlignIO.fasta'; 
			system "del $del > temp.txt";
		}
		my $file6=$ENV{'PWD'}.'\results\AlignIO_.fasta';
		system "move AlignIO_.fasta $file6 > temp.txt";
		system "copy RNAduplex.exe \.\\RNAxsresult > temp.txt";
		system "copy RNAfold.exe \.\\RNAxsresult > temp.txt";
		system "copy RNAplfold.exe \.\\RNAxsresult > temp.txt";
		system "copy cyggcc_s-1.dll \.\\RNAxsresult > temp.txt";
		system "copy cygwin1.dll \.\\RNAxsresult > temp.txt";
		system"perl RNAxs.pl -s $file6 -u $opts{u} -f 0 -g 0 -e 0 -t $opts{t} -a 0 -w 80 -l 40";
		my $RNAxsresult=$ENV{'PWD'}.'\\'.'RNAxsresult'.'\\output.csv';
		my $RNAxsresult_=$ENV{'PWD'}.'\\'.'results'.'\\Tar_acc_result.txt';
		#Analyzing RNAxs results
		open FH,">","$RNAxsresult_";
		open FH1,"<","$RNAxsresult";
		while (my $line=<FH1>){
		if ($line=~/^([A-Za-z]*),/){
			print FH ">habal\n",$1,"\n"; 
			}
		}
		close (FH);
		close (FH1);
		
		#Thermo-iscore evaluation
		my $Thermo_iscore = $ENV{'PWD'}.'\iscore_thermo.xls';
		my $Excel = Win32::OLE->GetActiveObject('Excel.Application')|| Win32::OLE->new('Excel.Application', 'Quit');  # use the Excel application if it's open, otherwise open new

		my $Book = $Excel->Workbooks->Open($Thermo_iscore); # open the file
		$Excel->Run("SeqInput");
		$Excel->Run("Secondary_Structure"); #macro_name
		$Book->Save; #optional - save any changes made by the macro
		$Book->Close;
		print "Iscore finished \n";
		#scores filtration
		scoring90 ($num,$trans_num);
		#snps filtration
		snps($num);
		#Off-target filtration
		if ($off_target_filtration==1){
			off_target($trans_num,$num,$opts{d},$opts{i},$opts{e});}
		if($off_target_filtration==0){
			No_off_target($trans_num,$num,$opts{d},$opts{i},$opts{e});}
		#MysiRNA-Model	
		MysiRNA_Model($num,$opts{w},$opts{s},$opts{d},$opts{o},$opts{t},$opts{u},$opts{p});
		my $file7=$ENV{'PWD'}.'\\'.'results'.'\\survey.log';
		chomp($acc);
		chomp($num);
		my $Outfile=$num."_MysiRNA_Designer_results_".$acc;
		system "del /Q /s $opts{o}\\$Outfile > temp.txt";	
		system "move survey.log $file7 > temp.txt";
		system "Xcopy /e/-Y/I/R .\\results $opts{o}\\$Outfile\\tools_results > temp.txt";
		system "Xcopy /e/-Y/I/R .\\blastresult $opts{o}\\$Outfile\\blastresult > temp.txt";
		system "Xcopy /e/-Y/I/R .\\Multi_score_Accepted $opts{o}\\$Outfile\\Multi_score_Accepted > temp.txt";
		system "Xcopy /e/-Y/I/R .\\RNAxsresult $opts{o}\\$Outfile\\RNAxsresult > temp.txt";
		system "Xcopy /e/-Y/I/R .\\MysiRNA_Designer_Final_results $opts{o}\\$Outfile\\MysiRNA_Designer_Final_results > temp.txt";
		
	}
close Main;
}
else{

	usage();
}
sub usage{
print
"	
	||||||||||||||||||||||||||||||||||||||||||||||||
	||        Welcome To MysiRNA-Designer         ||
	||    A Software For Rational siRNA Design    ||
	||    Copyright (C) 2011  <M.Mysara et al>    ||
	||||||||||||||||||||||||||||||||||||||||||||||||

 MysiRNA-Designer version 1, Copyright (C) 2011, M.Mysara et al
 MysiRNA-Designer comes with ABSOLUTELY NO WARRANTY.
 This is free software, and you are welcome to redistribute it under
 certain conditions; type `notepad copying.txt\' for details.;
   
	
 Use the Following Mandatory Options:\n
 -s File containing mRNA(s) Refseq Accession(s) one/line, 
  <Example, D:\\Acc.txt containing 'NM_001211' accession number >
 -o out put PATH.
  <Example, D:\\ >
 -w WEKA installation PATH.
  <Example, C:\\Weka-3-7 >
	
 Use the Following Non-Mandatory Options:\n
 -d To perform off-target filtration, insert mRNA refseq dataset name
 without its extension and complete PATH, <Example, D:\\refseq_rna >.
 -t RNAxs threshold on the 8 nts accessibility.
 -u RNAxs threshold on the 16 nts accessibility.
 -p Use MysiRNA-Model high specificity threshold (Y or N).
 -e Updated 3UTR Database, else use default
 -i Updated refseq2ensumble ID Database, else use default 
	
 For Queries about the installaion, type \'notepad README.txt\'
 For Queries about the Copy rights, type \'notepad COPYING.txt\'
	
";
}
sub scoring90{
	#Local varriable
		my $filename = $ENV{'PWD'}.'\Iscore_thermo.xls';
	
		my $trackernum=0;
		my $Tracker=2; # to include siRNA number at the header, START WITH TWO AS THE 1ST TWO siRNA are rejected
	#////////////////	
	#testing 
	open TEST, ">>",$ENV{'PWD'}.'\Multi_score_Accepted\\habal.csv';
	#///////////
		open FH1, ">",$ENV{'PWD'}.'\Multi_score_Accepted\\'.$_[0]."_Multi_score_Approved_siRNA.fasta";
		print FH1 "";
		close (FH1);
		open FH2, ">>",$ENV{'PWD'}.'\Multi_score_Accepted\\'.$_[0]."_Multi_score_Approved_siRNA.fasta";
		open FH6, ">>",$ENV{'PWD'}.'\blastresult\\'.$_[0]."_Whole_siRNA_off_targets.csv";
	#............... open the conserved siRNA and read it in array...........
	my @cons_sirna;
	my $sirna_arr_length;
	if ($_[1]>1){
	open FH3, "<",$ENV{'PWD'}.'\results\AlignIO.fasta';
	 @cons_sirna=<FH3>;
	$sirna_arr_length=scalar(@cons_sirna);
	print $sirna_arr_length,"\t";
	close (FH3);
	} 
	my $string;
	#............... open the target_Accessibility siRNA and read it in array...........
	open FH4, "<",$ENV{'PWD'}.'\results\Tar_acc_result.txt';
	my @tar_Acc=<FH4>;
	my $tar_acc_length=scalar(@tar_Acc);
	print $tar_acc_length,"\t";
	close (FH4);
	my $tar_acc_string;
	#..................open thermo-iscore results.......
		my $parser   = Spreadsheet::ParseExcel->new();
		my $workbook = $parser->parse($filename);
		my @Habal;
		my @Habal_;
		if ( !defined $workbook ) {
			die $parser->error(), ".\n";
		}

		for my $worksheet ($workbook->worksheet('Designer')) {

		   my ( $row_min, $row_max ) = $worksheet->row_range();
			my ( $col_min, $col_max ) = $worksheet->col_range();
			print "Row min: ",$row_min," ","Row max: ",$row_max," ","Col min: ",$col_min," ","Col max: ",$col_max,"\n";
			my $cell = $worksheet->get_cell( 0, 1 );
			 my $Maximum=  $cell->value();
			my $cell__ = $worksheet->get_cell( 15, 3 );
			print "Value       = ", $cell__->value(),       "\n";
			my $counter=1;
			
			my $i=0;
			for my $row ( 18 .. $Maximum-3 ) {
			$Tracker++; #add 1 to the current tracker
			print TEST $Tracker.",";
			my $track=0;
			my $track1=0;
			#.................................................................................................
			my $celll = $worksheet->get_cell( $row, 1); #to get the siRNA sense
			my $celll_AS=$worksheet->get_cell( $row, 2); #to get the anitsense
			my $sirna=$celll->value();
			my $sirna_AS=$celll_AS->value();
			my $sirna_=substr($sirna,0,19);# to remove the extra 2 nt
			print TEST "$sirna_AS,$sirna_,";
			#.......................................consider only siRNA conserved by comparing to the @cons_sirna
			if ($_[1]==1){
				$track=1;
				print TEST "OK_TRANSCRIPT,";
			}
			else{
			for my $numb(0..($sirna_arr_length-1)){
					$string=$cons_sirna[$numb];
					chomp ($string);
					$string=~s/T/U/gi;
					#print $string,"\t";
					$string=substr($string,0,19);#no need but to double check
					if ($sirna_ =~ /$string/i){
						$track=1;
					}
				}
			}
			# do the whole scoring else don't do anything
	#..........................Target access filteration.................................
			for my $numbb(0..($tar_acc_length-1)){
				$tar_acc_string=$tar_Acc[$numbb];
				chomp ($tar_acc_string);
				$tar_acc_string=~s/T/U/gi;
				#print $tar_acc_string,"\t";
				$tar_acc_string=substr($tar_acc_string,0,19);#no need but to double check
				if ($sirna_ =~ /$tar_acc_string/i){
					$track1=1;
					print TEST "ok_accessibility,";
				}
				else{
				print TEST "false_accessibility,";
				}
			}
			# do the whole scoring else dont do anything
			#...........ThermoComp filteration.....................................................................................
		if($track==1&&$track1==1){
		print TEST "enter,";
			#if (($Thermocomp[$row-20]>.6) && ($Thermocomp[$row-20]<1)){
			my $Thermocomp = $worksheet->get_cell( $row, 121); #Thermo21
			my $Thermocomp_ = $worksheet->get_cell( $row, 120); #Thermo19
			if ($Thermocomp->value() >=0.70){	 
				print TEST "ok_thermo,";
				 my $iscore=$worksheet->get_cell($row,116); #26*4+nn #iscore
					if ($iscore->value() >=70){
						print TEST "ok_iscore,";	
					  my $dG= $worksheet->get_cell( $row, 112);
						my $S_Bio= $worksheet->get_cell( $row, 115); #s-Biopredsi
						if ($S_Bio->value() >= 0.7){
								print TEST "ok_bio,";
							my $DISR = $worksheet->get_cell( $row, 119); #Disr
							if ($DISR->value() >= 70){
								print TEST "ok_Dsir,";
								my $rey = $worksheet->get_cell( $row, 117);#reynold
								if ($rey ->value() > 1.9 &&$rey ->value()< 9.15){
										print TEST "ok_rey,";
									my $Ui = $worksheet->get_cell( $row, 111);#Ui-Tie
									if ($Ui->value() =~/[Ia,Ib,II,III]/){	
											print TEST "ok_ui,";
										my $Ama = $worksheet->get_cell( $row, 112);#Amarzguioui
										if($Ama->value() > -1.21 && $Ama->value()<5.3){
											print TEST "ok_ama,";
											my $kat = $worksheet->get_cell( $row, 118);#katoh
											if ($kat->value() >42.03 &&$kat->value()  < 97.01){
													print TEST "ok_ka,";
												my $Hs = $worksheet->get_cell( $row, 113);#Hsieh	
												if ($Hs->value() > -1.11 && $Hs->value() < 3.11){
													print TEST "ok_hs,";
													my $Taka = $worksheet->get_cell( $row, 114);#Takasaki
													if ($Taka->value() >-10.22&&$Taka->value() <14.06){
												
													#print "SCORE!!! \n";
													#Getting Total G 
													my $G = $worksheet->get_cell( $row, 108);
													print FH2 ">XX_000000 siRNAnum:",$Tracker,"_iscore:",$iscore->value(),"_SBiopredsi:",$S_Bio->value(),"_DISR:",$DISR->value(),"_reynold:", $rey ->value(),"_UiTie:",$Ui->value(),"_Amarzguioui:",$Ama->value(),"_Katoh:",$kat->value(),"_Hiesh:",$Hs->value(),"_Takasaki:",$Taka->value(),"_ThermoComp21:",$Thermocomp->value(),"_Thermo19:",$Thermocomp_->value(),"_Gwhole:",$G->value(),"_siRNA_AntiSense:",$sirna_AS,"\n",$sirna_,"\n";
													print FH6 $Tracker,$iscore->value(),",",$S_Bio->value(),",",$DISR->value(),",",$rey ->value(),",",$Ui->value(),",",$Ama->value(),",",$kat->value(),",",$Taka->value(),",",$Thermocomp->value(),",",$Thermocomp_->value(),",",$G->value(),",",$sirna_,"\n";
														print TEST "accepted,";
													$trackernum++;	
														
													}
												}
											}
										}
									}
								}
							}
						}
				   }
				   
				   next unless $cell;
				   $i++;    
					$counter ++;
	#................ ...................          
						}
				#	}
	#.......................................
				}
			print TEST "\n";
			}
		}
		
		#.................................................................
	print FH6 '_,_,_,_,_,_,_,_,_,_,_,_'; print FH6 "\n";
	close FH6;
	close (FH2);
	print "\n SCORING DONE \n";
	close TEST;
}

sub snps{
	my $id;
	open  FH,$ENV{'PWD'}."/results/fasta1result.gb";
	while (my $line=<FH>){
		if ($line =~/^VERSION/){
			$line=~/GI:([0-9]+)/;#VERSION     NM_001211.5  GI:168229167
			$id = $1;
		}
	}

	print $id,"\n";
	my $eutil = Bio::DB::EUtilities->new(-eutil => 'elink',
										-id    => $id,
										-email  => 'setyourown@foo.bar',
										-verbose   => 1,
										-dbfrom => 'nuccore',
										-db  => 'snp',
										-cmd   => 'neighbor_history',
									);
	 
	my $hist = $eutil->next_History || die "No history data returned";
	 
	$eutil->set_parameters(-eutil => 'efetch',
						  -history   => $hist,
						  -retmode => 'text',
						  # 'chr', 'flt', 'brief', 'rsr', 'docset'
						  -rettype => 'brief'  
	);
	 
	$eutil->get_Response(-file =>$ENV{'PWD'}."/results/snps.txt",-format=>'fasta');
	 
	my $counter=0;
	open FH1, $ENV{'PWD'}."/results/snps.txt";
	open FH2,">",$ENV{'PWD'}."/results/snps_region.fasta"; print FH2 ""; close FH2;
	open FH2,">>",$ENV{'PWD'}."/results/snps_region.fasta";
	while (my $line_=<FH1>){
		$line_=~s/T/U/gi;
		#getting all possibilities of SNPs, and rejecting any siRNA with one or more snps
		if  ($line_=~/([A,G,C,U]{18})\[([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\]([A,G,C,U]{18})/i){print FH2 $1,$2,$8,"\n"; print FH2 $1,$3,$8,"\n";print FH2 $1,$4,$8,"\n";print FH2 $1,$6,$8,"\n";print FH2 $1,$5,$8,"\n";print FH2 $1,$7,$8,"\n";}#$line_=~/[A,G,C,U]{18}\[([A,U,C,G,-])*[\/[A,U,C,G,-]*]*\]([A,G,C,U]{18})/i
		if  ($line_=~/([A,G,C,U]{18})\[([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\]([A,G,C,U]{18})/i){print FH2 $1,$2,$7,"\n"; print FH2 $1,$3,$7,"\n";print FH2 $1,$4,$7,"\n";print FH2 $1,$5,$7,"\n";print$1,$6,$7,"\n";}#$line_=~/[A,G,C,U]{18}\[([A,U,C,G,-])*[\/[A,U,C,G,-]*]*\]([A,G,C,U]{18})/i
		if  ($line_=~/([A,G,C,U]{18})\[([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\]([A,G,C,U]{18})/i){print FH2 $1,$2,$6,"\n"; print FH2 $1,$3,$6,"\n";print FH2 $1,$4,$6,"\n";print FH2 $1,$5,$6,"\n";}#$line_=~/[A,G,C,U]{18}\[([A,U,C,G,-])*[\/[A,U,C,G,-]*]*\]([A,G,C,U]{18})/i
		if  ($line_=~/([A,G,C,U]{18})\[([A,U,C,G,-]*)\/([A,U,C,G,-]*)\/([A,U,C,G,-]*)\]([A,G,C,U]{18})/i){print FH2 $1,$2,$5,"\n"; print FH2 $1,$3,$5,"\n";print FH2 $1,$4,$5,"\n";}#$line_=~/[A,G,C,U]{18}\[([A,U,C,G,-])*[\/[A,U,C,G,-]*]*\]([A,G,C,U]{18})/i
		if  ($line_=~/([A,G,C,U]{18})\[([A,U,C,G,-]*)\/([A,U,C,G,-]*)\]([A,G,C,U]{18})/i){print FH2 $1,$2,$4,"\n"; print FH2 $1,$3,$4,"\n";}#$line_=~/[A,G,C,U]{18}\[([A,U,C,G,-])*[\/[A,U,C,G,-]*]*\]([A,G,C,U]{18})/i
	}
	close FH2; 
	open FH4,">",$ENV{'PWD'}."/results/SNPs_free_siRNAs.fasta"; print FH4 "";close FH4;
	my $out_ = Bio::SeqIO->newFh ( -file => ">".$ENV{'PWD'}."/results/SNPs_free_siRNAs.fasta", -format => 'fasta');
	my $quit=0;
	my $seq;
	my $str = Bio::SeqIO->new(-file=>$ENV{'PWD'}.'\Multi_score_Accepted\\'.$_[0].'_Multi_score_Approved_siRNA.fasta' , -format => 'fasta' );
	while (my $input = $str->next_seq()){
		$quit=0;
		$seq=$input->seq;
		open FH3,$ENV{'PWD'}."/results/snps_region.fasta";
		while (my $line1=<FH3>){
			last if $quit==1; 
			if ($line1=~/$seq/){
				$quit=1;
			}
		}
		if ($quit==0){
			print $out_ $input; 
			print $seq,"\n";
		}
	}
}
sub off_target{
	#....making the result folder
	my $num=1;
	my $quit=0;
	my $Tracer;
	my $flage;
	my $blast_db=$_[2];
	#////////////get the current path	
	my $filename=$ENV{'PWD'}.'\blastresult'.'\\';
	my $filename1=$ENV{'PWD'}.'\results\Blastquery.fasta';
	#..........geting the @ of all input acc
	open FH,"<",$ENV{'PWD'}.'\results\blast1result.fasta';
	my $t_acc;
	if ($_[0]>1){
		while (my$s_acc=<FH>){
			my $NN=chomp($s_acc);$t_acc=$t_acc.$s_acc;
		}
	}
	#..........to ensure that Acc of the mRNA and its transcript is rejected
	else{  
			my $str_ = Bio::SeqIO->new(-file=>$ENV{'PWD'}.'\results\fasta1result.fasta' , -format => 'fasta' );
			while (my $input_ = $str_->next_seq()){$t_acc=$input_->id; print "$t_acc\n";}
	}
	#...........to check if the snps worked and produced the snps_score_final.fasta or to shift to iscore_final.fasta
	 my $checkfile=$ENV{'PWD'}.'\results\SNPs_free_siRNAs.fasta';
	 my $checkfile_=$ENV{'PWD'}.'\Multi_score_Accepted\\'.$_[1]."_Multi_score_Approved_siRNA.fasta";
	 my $str;
	 if (-e "$checkfile"){ 
		$str = Bio::SeqIO->new(-file=>"$checkfile" , -format => 'fasta' );
	}
	else{
		$str = Bio::SeqIO->new(-file=>"$checkfile_" , -format => 'fasta' );
	}
	  
	while (my $input = $str->next_seq()){
		$Tracer=1;
		#.......... reseting 
		system"del Temp.fasta > temp.txt";
		#...............getting the sirna seq and 2-7 seq
		my $queryseq=$input->seq();
		my $seed=$input->subseq(2,7);
		#......................
		print $queryseq,"\n";
		open Blastquery, ">", $ENV{'PWD'}.'\results\Blastquery.fasta';
		print Blastquery ">$t_acc Homo sapiens, mRNA.\n",$queryseq;
		close Blastquery;
		system "blastall.exe -p blastn -d $blast_db -i $filename1 -o Temp.fasta -W 7 -e 1000 -q -1 -G 1 -E 2 -I T";
		print "blastall done\n";
		#...................reading the result file
		my $report_obj = new Bio::SearchIO(-format => 'blast',                                   
									  -file   => "Temp.fasta");   
		while( my $result = $report_obj->next_result ){     
			while( my $hit = $result->next_hit ){       
				#if the hit with one of the gene transcripts, Do nothing			
				my $hit_acc=$hit->accession; 		
				if ($t_acc=~/$hit_acc/){}		
				else{					
					if ($hit->description() =~/^Homo sapiens/ ){
						#Homosapiens only species 
						print $hit->accession()." 	".$result->query_accession(),"\n";
						while( my $hsp = $hit->next_hsp ){
							# searching the hists
							if ( ($hsp->length('query')==19 && $hsp->percent_identity > 93 )||($hsp->length('query')==18 && $hsp->percent_identity ==100)) { print "1st stage filteration\n"; $Tracer=0;}            
							else{
								my @filter=$hsp->seq_inds('query','identical');
								#to convert array to string			 
								my $String= join('',@filter);
								my @hit_range=$hsp->range('hit');
								my $hit_range_String=join('',@hit_range);
								if ( $String=~/([0-9]*)245678([0-9]*)/){
									print "Seed matching..check if it matches 3UTR \(true -ve\).. or not \(false -ve\)!!!! \n";
									my $hit_acc=$hit->accession();							
									#......................
									open FHH,$_[3];
									my $enst;
									my $hitid = $hit->accession();
									#to make sure no accession versions are included
									if($hitid=~/([\w]*)\./){
										print $1,"\n";
										$hitid=$1;
									}
									while (my $id=<FHH>){ 
										if ($id=~/^([\w]*),$hitid/){
											$enst=$1;#to take the ensumbl id corrosponding to refseq accession
											print"\.\.\.\. $enst \.\.\.\.\.";
											
											}
									}
									close FHH;
									if(!$enst){#ID still couldn't be found,  reject the siRNA
										$Tracer=0;
										print "siRNA's off-target id couldn't be found, rejecting the siRNA\n";
									}
									else{
										my $in  = Bio::SeqIO->new(-file => $_[4],
																   -format => 'fasta');
										$flage=0;#to know when to use ensembl id to get 3UTR in case we couldn't find in lacal database
										while ( my $seq = $in->next_seq() ){
										
											if ($seq->id=~/$enst/){
												$flage=1;
												print $seed,"\n";
												if ($seq->seq=~/$seed/){
													print "True -ve \n";
													$Tracer=0;
												}
												else{
													print "False -ve \n";
													
												}
											} 
										}
									}
									#..........................
									if ($flage==0){
												##Cant find ENSG ID in the ENSG_UTR database, reject the siRNA 
												$Tracer=0;
												print "True -ve \n";
												
									}			
						#............................			
								}
							}		
						last if $Tracer == 0;
						}
					last if $Tracer == 0;
					}
				last if $Tracer == 0;
				}	
				last if $Tracer == 0;
			}
		last if $Tracer == 0; 
		}	 
		if($Tracer==1){
			my $out_ = Bio::SeqIO->newFh ( -file => ">>".$ENV{'PWD'}."/MysiRNA_Designer_Final_results/".$_[1]."_Final_siRNA.fasta", -format => 'fasta');
			print $out_ $input;
			 #open FH6,">>", $ENV{'PWD'}.'\blastresult\Final.fasta';
			 #print FH6 $queryseq,"\n";
			 my $off_target=$filename."siRNA_num_".$num."_Off_Target\.fasta";
			 system "copy Temp.fasta $off_target > temp.txt"; 
		}	
		$num ++;
	}
	print "off_target done!!! \n !";
}
sub No_off_target{

	#....making the result folder
	my $num=1;
	my $quit=0;
	my $Tracer;
	my $flage;
	my $blast_db=$_[2];
	#...............get the current path	
	my $filename=$ENV{'PWD'}.'\blastresult'.'\\';
	my $filename1=$ENV{'PWD'}.'\results\Blastquery.fasta';
	#.............. geting the @ of all input acc
	open FH,"<",$ENV{'PWD'}.'\results\blast1result.fasta';
	my $t_acc;
	if ($_[0]>1){
		while (my$s_acc=<FH>){
			my $NN=chomp($s_acc);$t_acc=$t_acc.$s_acc;
		}
	}
	#................to ensure that Acc of the mRNA and its transcript is rejected
	else{  
			my $str_ = Bio::SeqIO->new(-file=>$ENV{'PWD'}.'\results\fasta1result.fasta' , -format => 'fasta' );
			while (my $input_ = $str_->next_seq()){$t_acc=$input_->id; print "$t_acc\n";}
	}
	#...................to check if the snps worked and produced the snps_score_final.fasta or to shift to iscore_final.fasta
	 my $checkfile=$ENV{'PWD'}.'\results\SNPs_free_siRNAs.fasta';
	 my $checkfile_=$ENV{'PWD'}.'\Multi_score_Accepted\\'.$_[1]."_Multi_score_Approved_siRNA.fasta";
	 my $str;
	 if (-e "$checkfile"){ 
		$str = Bio::SeqIO->new(-file=>"$checkfile" , -format => 'fasta' );
	}
	else{
		$str = Bio::SeqIO->new(-file=>"$checkfile_" , -format => 'fasta' );
	}
	  
	while (my $input = $str->next_seq()){ 
		{
			my $out_ = Bio::SeqIO->newFh ( -file => ">>".$ENV{'PWD'}."/MysiRNA_Designer_Final_results/".$_[1]."_Final_siRNA.fasta", -format => 'fasta');
			print $out_ $input;
			my $off_target=$filename."siRNA_num_".$num."_Off_Target\.fasta";
			system "copy Temp.fasta $off_target > temp.txt"; 
		}	
		$num ++;
	}
	print "off_target Skipped!!! \n !";
}
sub MysiRNA_Model{
	open FH1,">>",$ENV{'PWD'}.'\MysiRNA_Designer_Final_results\\'.$_[0].'_Final_MysiRNA.fasta';
	open FH2, $ENV{'PWD'}.'\results\blast1result.fasta';
	my $ID_total="";
	while (my $ID=<FH2>){
		chomp $ID;
		$ID_total=$ID_total.$ID."\n";
	}
	print FH1 "
	||||||||||||||||||||||||||||||||||||||||||||||||
	||        MysiRNA-Designer Results            ||
	||    A Software For Rational siRNA Design    ||
	||    Copyright (C) 2011  <M.Mysara et al>    ||
	||||||||||||||||||||||||||||||||||||||||||||||||

These siRNA designed to target mRNA(s) with the following Accession:
 $ID_total
Using the following runtime parameters:

mRNA refseq database PATH: $_[3] 
mRNA input accession number PATH: $_[2]
Output PATH: $_[4]
WEKA installation PATH: $_[1]
RNAxs threshold on the 8 nts accessibility: $_[5]
RNAxs threshold on the 16 nts accessibility: $_[6]
Use MysiRNA-Model High specificity Threshold: $_[7]


";
	my $inputfile=$ENV{'PWD'}.'\MysiRNA_Designer_Final_results\\'.$_[0].'_Final_siRNA.fasta';
	my $test=$ENV{'PWD'}.'\results\test.arff';
	my $str = Bio::SeqIO->new(-file=>$inputfile,-format => 'fasta' );
	my $Thermo=0;
	my $iscore=0;
	my $DSIR=0;
	my $DG=0;
	my $PATH_=$ENV{'PWD'};
	while (my $input = $str->next_seq()){
		open FH,">",$test;
		my $header= $input->desc;
		my $seq=$input->seq;
		if ($header=~/_iscore:([0-9]*[\.]?[0-9]*)_/){
			$iscore=$1;
		}
		if ($header=~/_Gwhole:(-[0-9]*[\.]?[0-9]*)_/){
			$DG=$1;
		}
		if ($header=~/_ThermoComp21:([0-9]*[\.]?[0-9]*)_/){
			$Thermo=$1;
		}
		print FH "\@relation sirnahabal\n\n\@attribute ThermoComposition21 numeric\n\@attribute i-Score numeric\n\@attribute Whole-DG numeric\n\@attribute ' Inhibition' numeric\n\n\@data\n$Thermo,$iscore,$DG,$DSIR\n";
		close FH;
		#.................
		my $model=$ENV{'PWD'}.'\MysiRNA_Model.model';
		my $result=$ENV{'PWD'}.'\results\Weka_result.txt';
		chdir "$_[1]";
		system "java weka.classifiers.functions.MultilayerPerceptron -l $model -T $test -p 0 > $result";
		chdir "$PATH_";
		#................
		open FH2,$result;
		my @mysirna=<FH2>;
		$mysirna[5]=~/ 1      0         ([0-9]*\.?[0-9]*)     [0-9]*\.?[0-9]* /;
	    my $MysiRNA_Score=$1;
	#	$header=~/(\>[\w\W]*)siRNA_AntiSense\:([A-Za-z]*)/;
		$header=~s/_siRNA_AntiSense:/_MysiRNA:$MysiRNA_Score\nAntiSense:	/gi;
		if ($_[7]=~/N/){
		print FH1 ">",$header,"\nSense:		$seq\n";
		}
		if ($_[7]=~/Y/){
			if ($MysiRNA_Score>=93){
			print FH1 ">",$header,"\nSense:		$seq\n";
			}
		}
		close FH2;
	}
	close FH1;
	print "MysiRNA_Model DONE!!!!\n";
}