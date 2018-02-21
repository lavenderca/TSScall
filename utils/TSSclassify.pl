#!/usr/bin/perl
#Written by Adam Burkholder, National Institute of Environmental Health Sciences, 2016-17

die "Error: detail and GTF files must be specified\nUsage: TSSclassify [detail file] [GTF]\n" if scalar(@ARGV)==0;
die "Error: GTF file must be specified\nUsage: TSSclassify [detail file] [GTF]\n" if scalar(@ARGV)==1;

open(IN,$ARGV[0]) || die "Error: could not open detail file \"$ARGV[0]\"\n";
while(1) {
	$char=getc(IN);
	if($char eq "\n") {
		seek(IN,0,0);
		last;
	}
	elsif($char eq "\r") {
		$nchar=getc(IN);
		if($nchar eq "\n") {
			seek(IN,0,0);
			last;
		}
		else {
			$/="\r";
			seek(IN,0,0);
			last;
		}
	}
}
$head=<IN>;
$head=~s/[\r\n]//g;
@fields=split(/\t/,$head);
$col_count=scalar(@fields);
$name_field=-1;
$type_field=-1;
for($i=0;$i<scalar(@fields);$i++) {
	if($fields[$i] eq "gene_name") {
		$name_field=$i;
	}
	if($fields[$i] eq "parent_type" || $fields[$i] eq "type") {
		$type_field=$i;
	}
}
die "Error: \"$ARGV[0]\" does not appear to be a properly formatted detail file\n" if $fields[0] ne "TSS ID" && scalar(@fields)!=15;
$head.="\tIntersects Gene Model?\tTranscript ID(s)\tGene ID(s)\tGene Name(s)\tType(s)\tSense\/Anti-sense\/Both?\tDistance from obsTSS\tIntersects Exon?\tNearest obsTSS\tTranscript ID(s)\tGene ID(s)\tGene Names(s)\tType(s)\tDistance from nearest\tNearest strand\tNearest obsTSS to [N] uTSS\'s\tSelected 5\' UTR\tObserved Exon 1 Length\tObserved 5\' UTR Length\tCDS Start Within Exon\n";
while(<IN>) {
	$_=~s/[\r\n]//g;
	@fields=split(/\t/,$_);
	push(@table,[@fields,0,"X","X","X","X","X","X",0,"X","X","X","X","X","X","X","X","X","X","X","X"]);	#store each line of input file in anonymous array, add default values for new columns, and push onto @table
	push(@prevs,-1);					#store dummy values for location of previous and next obsTSS
	push(@nexts,-1);
	@tx=split(/\;/,$fields[2]);
	foreach(@tx) {
		$obs{$_}=$fields[6];
	}
}
close(IN);
@table=sort { $a->[5] cmp $b->[5] || $a->[6] <=> $b->[6] || $a->[4] cmp $b->[4] } @table;	#sort by chromosome, position, and strand

open(IN,$ARGV[1]) || die "Error: could not open GTF file \"$ARGV[1]\"\n";
while(1) {
	$char=getc(IN);
	if($char eq "\n") {
		$/="\n";
		seek(IN,0,0);
		last;
	}
	elsif($char eq "\r") {
		$nchar=getc(IN);
		if($nchar eq "\n") {
			$/="\n";
			seek(IN,0,0);
			last;
		}
		else {
			$/="\r";
			seek(IN,0,0);
			last;
		}
	}
}
while(<IN>) {
	$_=~s/[\r\n]//g;
	next if $_=~/^\#/;
	@F=split(/\t/,$_);
	next if $F[2] ne "exon" && $F[2] ne "CDS";
	die "Error: \"$ARGV[1]\" does not appear to be a properly formatted GTF file\n" if !($_=~/transcript_id \"([^\"]+)\"/) || scalar(@F)!=9 || ($F[6] ne "+" && $F[6] ne "-");
	$trans=$1;
	$gene="";
	if($F[8]=~/gene_id \"([^\"]+)\"/) {
		$gene=$1;
	}
	$name="";
	if($F[8]=~/gene_name \"([^\"]+)\"/) {
		$name=$1;
	}
	$type="";
	if($F[8]=~/parent_type \"([^\"]+)\"/) {
		$type=$1;
	}
	if($F[8]=~/type \"([^\"]+)\"/) {
		$type=$1;
	}
	$genes{$F[0]}={} if !exists($genes{$F[0]});					#store annotation file in two-level hash, the first key being chromsomes, the second transcript
	$genes{$F[0]}->{$trans}=[$F[6],$gene,$name,$type,[],[],[9999999999,0,0,0,0]] if !exists($genes{$F[0]}->{$trans});	#each transcript entry consists of an anonymous array storing the strand, gene_id, gene_name, type, transcript start & end, and a list of exon start and ends
	push(@{$genes{$F[0]}->{$trans}->[5]},[$F[3],$F[4]]) if $F[2] eq "exon";				#add current exon start and end to list
	if($F[2] eq "CDS") {
		$genes{$F[0]}->{$trans}->[6]->[0]=$F[3] if $F[3]<$genes{$F[0]}->{$trans}->[6]->[0];
		$genes{$F[0]}->{$trans}->[6]->[1]=$F[4] if $F[4]>$genes{$F[0]}->{$trans}->[6]->[1];
	}
}
close(IN);
foreach $chr(keys(%genes)) {											#loop through transcripts and sort exons by start position
	foreach $trans(keys(%{$genes{$chr}})) {
		$genes{$chr}->{$trans}->[5]=[sort {$a->[0] <=> $b->[0]} @{$genes{$chr}->{$trans}->[5]}];
		if($genes{$chr}->{$trans}->[0] eq "-") {							#determine which end of CDS interval to place inside exons, calculate necessary offset to report accurate exon number for minus strand genes
			$CDS_offset=scalar(@{$genes{$chr}->{$trans}->[5]})+1;
			$CDS_index=1;
		}
		else {
			$CDS_offset=0;
			$CDS_index=0;
		}
		foreach(@{$genes{$chr}->{$trans}->[5]}) {							#loop through exons, determine which contains CDS start, apply offset to correct value for minus strand genes
			$genes{$chr}->{$trans}->[6]->[2]++;
			last if $genes{$chr}->{$trans}->[6]->[$CDS_index]>=$_->[0] && $genes{$chr}->{$trans}->[6]->[$CDS_index]<=$_->[1];
		}
		$genes{$chr}->{$trans}->[6]->[2]=abs($genes{$chr}->{$trans}->[6]->[2]-$CDS_offset);
		if(exists($obs{$trans})) {
			if($genes{$chr}->{$trans}->[0] eq "+") {
				$genes{$chr}->{$trans}->[4]->[0]=$obs{$trans};							#obsTSS is transcript start
				$genes{$chr}->{$trans}->[4]->[1]=$genes{$chr}->{$trans}->[5]->[-1]->[1];			#end of last exon is transcript end
			}
			else {
				$genes{$chr}->{$trans}->[4]->[0]=$genes{$chr}->{$trans}->[5]->[0]->[0];				#start of first exon is transcript start
				$genes{$chr}->{$trans}->[4]->[1]=$obs{$trans};							#obsTSS is transcript end
			}
		}
		else {
			$genes{$chr}->{$trans}->[4]->[0]=$genes{$chr}->{$trans}->[5]->[0]->[0];				#start of first exon is transcript start
			$genes{$chr}->{$trans}->[4]->[1]=$genes{$chr}->{$trans}->[5]->[-1]->[1];			#end of last exon is transcript end
		}
		if($genes{$chr}->{$trans}->[6]->[1]==0) {
			$genes{$chr}->{$trans}->[6]->[2]="X";
			$genes{$chr}->{$trans}->[6]->[3]="X";
			if($genes{$chr}->{$trans}->[0] eq "-") {
				$genes{$chr}->{$trans}->[6]->[4]=$genes{$chr}->{$trans}->[4]->[1]-$genes{$chr}->{$trans}->[5]->[-1]->[0]+1;
			}
			else {
				$genes{$chr}->{$trans}->[6]->[4]=$genes{$chr}->{$trans}->[5]->[0]->[1]-$genes{$chr}->{$trans}->[4]->[0]+1;
			}
		}
		else {
			if($genes{$chr}->{$trans}->[0] eq "-") {													#determine observed exon 1 length and 5' UTR length
				$genes{$chr}->{$trans}->[6]->[4]=$genes{$chr}->{$trans}->[4]->[1]-$genes{$chr}->{$trans}->[5]->[-1]->[0]+1;			#do not calculate UTR length if exon 1 length is less than 1
				if($genes{$chr}->{$trans}->[6]->[4]>0) {
					for($i=scalar(@{$genes{$chr}->{$trans}->[5]})-$genes{$chr}->{$trans}->[6]->[2];$i<scalar(@{$genes{$chr}->{$trans}->[5]});$i++) {			#loop through exons, starting at exon containing the CDS start
						if($i==(scalar(@{$genes{$chr}->{$trans}->[5]})-$genes{$chr}->{$trans}->[6]->[2])) {								#replace the annotated exon start and end with obsTSS and CDS start where appropriate
							$exon_start=$genes{$chr}->{$trans}->[6]->[1]+1;
						}
						else {
							$exon_start=$genes{$chr}->{$trans}->[5]->[$i]->[0];
						}
						if($i==(scalar(@{$genes{$chr}->{$trans}->[5]})-1)) {
							$exon_end=$genes{$chr}->{$trans}->[4]->[1];
						}
						else {
							$exon_end=$genes{$chr}->{$trans}->[5]->[$i]->[1];
						}
						$genes{$chr}->{$trans}->[6]->[3]+=($exon_end-$exon_start+1);
					}
				}
				else {
					$genes{$chr}->{$trans}->[6]->[3]="X";
				}
			}
			else {
				$genes{$chr}->{$trans}->[6]->[4]=$genes{$chr}->{$trans}->[5]->[0]->[1]-$genes{$chr}->{$trans}->[4]->[0]+1;
				if($genes{$chr}->{$trans}->[6]->[4]>0) {
					for($i=$genes{$chr}->{$trans}->[6]->[2]-1;$i>=0;$i--) {
						if($i==($genes{$chr}->{$trans}->[6]->[2]-1)) {
							$exon_end=$genes{$chr}->{$trans}->[6]->[0]-1;
						}
						else {
							$exon_end=$genes{$chr}->{$trans}->[5]->[$i]->[1];
						}
						if($i==0) {
							$exon_start=$genes{$chr}->{$trans}->[4]->[0];
						}
						else {
							$exon_start=$genes{$chr}->{$trans}->[5]->[$i]->[0];
						}
						$genes{$chr}->{$trans}->[6]->[3]+=($exon_end-$exon_start+1);
					}
				}
				else {
					$genes{$chr}->{$trans}->[6]->[3]="X";
				}
			}
		}
	}
}

#	%genes:
#		{ chromosome }:
#			-> { transcript id }:
#				-> [ ]:
#					-> strand
#					-> gene_id
#					-> gene_name
#					-> type
#					-> [ ]:
#						-> transcript start
#						-> transcript end
#					-> [ ]:
#						-> [ ]:
#							-> exon 1 start
#							-> exon 1 end
#						...
#						-> [ ]:
#							-> exon N start
#							-> exon N end
#					-> [ ]:
#						-> CDS start
#						-> CDS end
#						-> CDS start within exon
#						-> observed 5' UTR length
#						-> observed exon 1 length

for($i=0;$i<scalar(@table);$i++) {					#loop through @table from start to end, every obsTSS encountered becomes "previous" for subsequent uTSS's
	if($table[$i]->[0]=~/obsTSS/) {
		$prev=$i;
		$prevchr=$table[$i]->[5];
		$winners{$table[$i]->[0]}=0;
	}
	else {
		$prevs[$i]=$prev if $table[$i]->[5] eq $prevchr;
	}
	if($table[$i]->[0]=~/obsTSS|annoTSS/) {
		@translist=sort { $genes{$table[$i]->[5]}->{$a}->[6]->[3] <=> $genes{$table[$i]->[5]}->{$b}->[6]->[3] || length($a) <=> length($b) || $a cmp $b } split(/\;/,$table[$i]->[2]);			#select transcript with shorted 5' UTR, excluding any without annotated CDS and negative length UTR
		$table[$i]->[$col_count+16]=$translist[0];																			#for equal UTR length, select that with shortest transcript ID, then first in alphabetical sort
		$table[$i]->[$col_count+17]=$genes{$table[$i]->[5]}->{$translist[0]}->[6]->[4];
		$table[$i]->[$col_count+18]=$genes{$table[$i]->[5]}->{$translist[0]}->[6]->[3];
		$table[$i]->[$col_count+19]=$genes{$table[$i]->[5]}->{$translist[0]}->[6]->[2];
		if($table[$i]->[$col_count+18] eq "X" || $table[$i]->[$col_count+18]<0) {
			foreach(@translist) {
				if($genes{$table[$i]->[5]}->{$_}->[6]->[3] ne "X" && $genes{$table[$i]->[5]}->{$_}->[6]->[3]>=0) {
					$table[$i]->[$col_count+16]=$_;
					$table[$i]->[$col_count+17]=$genes{$table[$i]->[5]}->{$_}->[6]->[4];
					$table[$i]->[$col_count+18]=$genes{$table[$i]->[5]}->{$_}->[6]->[3];
					$table[$i]->[$col_count+19]=$genes{$table[$i]->[5]}->{$_}->[6]->[2];
					last;
				}
			}
		}
	}
}
for($i=scalar(@table)-1;$i>-1;$i--) {					#loop through @table from end to start, every obsTSS encountered becomes "next" for subsequen uTSS's
	if($table[$i]->[0]=~/obsTSS/) {
		$next=$i;
		$nextchr=$table[$i]->[5];
	}
	else {
		$nexts[$i]=$next if $table[$i]->[5] eq $nextchr;
	}
}
for($i=0;$i<scalar(@table);$i++) {					#loop through @table again, selecting closest obsTSS, "previous" or "next"
	next if $table[$i]->[0]=~/obsTSS/;
	next if $prevs[$i]==-1 && $nexts[$i]==-1;
	if($prevs[$i]!=-1) {
		$prevdist=$table[$i]->[6]-$table[$prevs[$i]]->[6];
	}
	else {
		$prevdist=999999999;
	}
	if($nexts[$i]!=-1) {
		$nextdist=$table[$i]->[6]-$table[$nexts[$i]]->[6];
	}
	else {
		$nextdist=999999999;
	}
	if(abs($prevdist)<abs($nextdist)) {
		$prevstr=$table[$prevs[$i]]->[4];
		$prevdist*=-1 if $prevstr eq "-";
		$table[$i]->[$col_count+8]=$table[$prevs[$i]]->[0];
		$table[$i]->[$col_count+9]=$table[$prevs[$i]]->[2];
		$table[$i]->[$col_count+10]=$table[$prevs[$i]]->[3];
		$table[$i]->[$col_count+11]=$table[$prevs[$i]]->[$name_field] if $name_field>=0;
		$table[$i]->[$col_count+12]=$table[$prevs[$i]]->[$type_field] if $type_field>=0;
		$table[$i]->[$col_count+13]=$prevdist;
		if($table[$i]->[4] eq $prevstr) {
			$table[$i]->[$col_count+14]="same";
		}
		else {
			$table[$i]->[$col_count+14]="opposite";
		}
		$winners{$table[$prevs[$i]]->[0]}++ if $table[$i]->[0]!~/annoTSS/;
	}
	elsif(abs($nextdist)<abs($prevdist)) {
		$nextstr=$table[$nexts[$i]]->[4];
		$nextdist*=-1 if $nextstr eq "-";
		$table[$i]->[$col_count+8]=$table[$nexts[$i]]->[0];
		$table[$i]->[$col_count+9]=$table[$nexts[$i]]->[2];
		$table[$i]->[$col_count+10]=$table[$nexts[$i]]->[3];
		$table[$i]->[$col_count+11]=$table[$nexts[$i]]->[$name_field] if $name_field>=0;
		$table[$i]->[$col_count+12]=$table[$nexts[$i]]->[$type_field] if $type_field>=0;
		$table[$i]->[$col_count+13]=$nextdist;
		if($table[$i]->[4] eq $nextstr) {
			$table[$i]->[$col_count+14]="same";
		}
		else {
			$table[$i]->[$col_count+14]="opposite";
		}
		$winners{$table[$nexts[$i]]->[0]}++ if $table[$i]->[0]!~/annoTSS/;
	}
	elsif($table[$i]->[4] eq "+") {				#if previous and next distances equal, select upstream relative to uTSS
		$prevstr=$table[$prevs[$i]]->[4];
		$prevdist*=-1 if $prevstr eq "-";
		$table[$i]->[$col_count+8]=$table[$prevs[$i]]->[0];
		$table[$i]->[$col_count+9]=$table[$prevs[$i]]->[2];
		$table[$i]->[$col_count+10]=$table[$prevs[$i]]->[3];
		$table[$i]->[$col_count+11]=$table[$prevs[$i]]->[$name_field] if $name_field>=0;
		$table[$i]->[$col_count+12]=$table[$prevs[$i]]->[$type_field] if $type_field>=0;
		$table[$i]->[$col_count+13]=$prevdist;
		if($table[$i]->[4] eq $prevstr) {
			$table[$i]->[$col_count+14]="same";
		}
		else {
			$table[$i]->[$col_count+14]="opposite";
		}
		$winners{$table[$prevs[$i]]->[0]}++ if $table[$i]->[0]!~/annoTSS/;
	}
	else {
		$nextstr=$table[$nexts[$i]]->[4];
		$nextdist*=-1 if $nextstr eq "-";
		$table[$i]->[$col_count+8]=$table[$nexts[$i]]->[0];
		$table[$i]->[$col_count+9]=$table[$nexts[$i]]->[2];
		$table[$i]->[$col_count+10]=$table[$nexts[$i]]->[3];
		$table[$i]->[$col_count+11]=$table[$nexts[$i]]->[$name_field] if $name_field>=0;
		$table[$i]->[$col_count+12]=$table[$nexts[$i]]->[$type_field] if $type_field>=0;
		$table[$i]->[$col_count+13]=$nextdist;
		if($table[$i]->[4] eq $nextstr) {
			$table[$i]->[$col_count+14]="same";
		}
		else {
			$table[$i]->[$col_count+14]="opposite";
		}
		$winners{$table[$nexts[$i]]->[0]}++ if $table[$i]->[0]!~/annoTSS/;
	}
}
print $head;
for($i=0;$i<scalar(@table);$i++) {									#final pass through @table
	$id=$table[$i]->[0];
	$chr=$table[$i]->[5];
	$str=$table[$i]->[4];
	$pos=$table[$i]->[6];
	@list_keys=();
	%list=();
	if($id=~/obsTSS/) {
		@list_keys=split(/;/,$table[$i]->[3]);
		foreach(@list_keys) {
			$list{$_}=1;
		}
	}
	foreach $trans(keys(%{$genes{$chr}})) {								#loop through all transcripts on current TSS chromosome
		next if $id=~/obsTSS/ && exists($list{$genes{$chr}->{$trans}->[1]});
		if($pos>=$genes{$chr}->{$trans}->[4]->[0] && $pos<=$genes{$chr}->{$trans}->[4]->[1]) {	#check first for overlap with whole transcript
			$table[$i]->[$col_count]=1;
			$table[$i]->[$col_count+1].="\;$trans";
			$table[$i]->[$col_count+2].="\;$genes{$chr}->{$trans}->[1]" if $genes{$chr}->{$trans}->[1] ne "";
			$table[$i]->[$col_count+3].="\;$genes{$chr}->{$trans}->[2]" if $genes{$chr}->{$trans}->[2] ne "";
			$table[$i]->[$col_count+4].="\;$genes{$chr}->{$trans}->[3]" if $genes{$chr}->{$trans}->[3] ne "";
			if($str eq $genes{$chr}->{$trans}->[0]) {
				$table[$i]->[$col_count+5]="sense" if $table[$i]->[$col_count+5] eq "X";
				$table[$i]->[$col_count+5]="both" if $table[$i]->[$col_count+5] eq "anti-sense";
			}
			else {
				$table[$i]->[$col_count+5]="anti-sense" if $table[$i]->[$col_count+5] eq "X";
				$table[$i]->[$col_count+5]="both" if $table[$i]->[$col_count+5] eq "sense";
			}
			foreach(@{$genes{$chr}->{$trans}->[5]}) {					#check for overlap with exons, quit looping after first overlap found
				last if $table[$i]->[$col_count+7]==1;
				if($pos>=$_->[0] && $pos<=$_->[1]) {
					$table[$i]->[$col_count+7]=1;
				}
			}
			if(exists($obs{$trans})) {							#determine distance from obsTSS, if a call was made; report value for transcript with shortest distance
				$disttemp=$pos-$obs{$trans};
				$disttemp*=-1 if $genes{$chr}->{$trans}->[0] eq "-";
				$table[$i]->[$col_count+6]=$disttemp if $table[$i]->[$col_count+6] eq "X" || ($disttemp<$table[$i]->[$col_count+6] && $disttemp>0);
			}
		}
	}
	$table[$i]->[$col_count+15]=$winners{$id} if $id=~/obsTSS/;
	@cols=();
	for($j=12;$j<$col_count;$j++) {
		push(@cols,$j);
	}
	for($j=($col_count+1);$j<=($col_count+4);$j++) {
		push(@cols,$j);
	}
	for($j=($col_count+8);$j<=($col_count+12);$j++) {
		push(@cols,$j);
	}
	foreach $k(@cols) {
		@vals=split(/\;/,$table[$i]->[$k]);
		next if scalar(@vals)==1;
		shift(@vals) if $vals[0] eq "X";
		%unique=();
		foreach(@vals) {
			$unique{$_}=0;
		}
		$table[$i]->[$k]=join("\;",keys(%unique));
	}
	$line=join("\t",@{$table[$i]});
	print "$line\n";
}
