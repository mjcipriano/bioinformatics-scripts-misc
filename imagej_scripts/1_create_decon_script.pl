#!/usr/bin/perl


use Getopt::Long;
use File::Copy;
use File::Path;
use File::Spec;


use strict;

# Requires as argument the full path to the files you would like to process

# ./1_decon_script.pl --indir=/home/mcipriano/Desktop/ij_scripts/test --outdir=/home/mccipriano/Desktop/ij_scripts/test/deconvolved 

# Rules: 1) Create directory and place all files in directory (called root directory or indir)
#        2) Create subdirectory named psf and place psf files in this directory named PSF_DAPI.h5, PSF_FITC.h5, and PSF_TRITC.h5
#        3) Create subdirectory named deconvolved (called outdir)

my $indir;
my $outdir;
my $celldir;
my $basename;

my $dx = "157.3";
my $dy = "157.3";
my $dz = "200";
my $dt = "1.000";

my $iterations = "100";
my $qual = "0.01";
my $timeOut = "10000";

my $bgMode = "auto";
my $brMode = "auto";
my $blMode = "auto";
my $pad = "auto";
my $mode = "fast";
my $bg = "0";


my $micType = "widefield";
my $riMedium = "1.415";
my $riLens = "1.515";
my $objQual = "good";
my $micNA = "1.40";

my $psfDir;
my $psfdapi;
my $psffitc;
my $psftritc;
my $copyDIC = 0;
my $taskID = 1;
my @taskIDS = ();

my $opts = GetOptions( 	"indir=s"	=>	\$indir,
			"outdir=s"	=>	\$outdir,
			"psfdir=s"	=>	\$psfDir,
			"basename=s"	=>	\$basename,
			"dx=f"		=>	\$dx,
			"dy=f"		=>	\$dy,
			"dt=f"		=>	\$dt,
			"iterations=i"	=>	\$iterations,
			"qual=f"	=>	\$qual,
			"timeout=i"	=>	\$timeOut,
			"bgMode=s"	=>	\$bgMode,
			"brMode=s"	=>	\$brMode,
			"blMode=s"	=>	\$blMode,
			"pad=s"		=>	\$pad,
			"mode=s"	=>	\$mode,
			"bg=i"		=>	\$bg,
			"micType=s"	=>	\$micType,
			"riMedium=f"	=>	\$riMedium,
			"riLens=f"	=>	\$riLens,
			"NA=f"		=>	\$micNA,
			"psfdapi=s"	=>	\$psfdapi,
			"psffitc=s"	=>	\$psffitc,
			"psftritc=s"	=>	\$psftritc,
			"copyDIC"	=>	\$copyDIC
			);

if(!defined($basename))
{
	$basename = 'unnamed';
}

if(!defined($outdir))
{
	$outdir = File::Spec->catdir($indir , "deconvolved");
	mkpath($outdir);

}

opendir(DIR, $indir);

my @files = readdir(DIR);

my $huygens_batch_name = File::Spec->catfile($indir, $basename . "_huygens_batch.hgsb");

open(BATCH, ">$huygens_batch_name");

my %filehash;

foreach my $file (@files)
{
	if(-f File::Spec->catfile($indir, $file) )
	{
		$filehash{$file} = 1;
		# print $file . "\n";
	
	}
}


my %roots;
my %roots_count;
my %filename2root;

my %excitationWL;
my %emissionWL;
my %snr;

$excitationWL{"dapi"} = 350;
$excitationWL{"fitc"} = 488;
$excitationWL{"tritc"} = 488;


$emissionWL{"dapi"} = 430;
$emissionWL{"fitc"} = 515;
$emissionWL{"tritc"} = 580;


$snr{"dapi"} = 33;
$snr{"fitc"} = 25;
$snr{"tritc"} = 6;

my %psfFiles;

if(!defined($psfDir))
{
	$psfDir         = File::Spec->catdir($indir, "psf");
}

if(!defined($psfdapi))
{
	$psfFiles{"dapi"}  = File::Spec->catfile($psfDir, "PSF_DAPI.h5");
}
if(!defined($psffitc))
{
	$psfFiles{"fitc"}  = File::Spec->catfile($psfDir, "PSF_FITC.h5");
}
if(!defined($psftritc))
{
	$psfFiles{"tritc"} = File::Spec->catfile($psfDir, "PSF_TRITC.h5");
}



my $macro_text = '';


# TODO: print decon header


while(my ($file, undef) = each(%filehash))
{
	# Find out what type of file it is
	if( $file =~ /^(.+_\d\d)_(w[1-6])(DAPI|TRITC|FITC|DIC)_(s\d+)(.+)$/)
	{
		# This matches a file with multiple wavelengths and multiple stage positions
		# rename it to be more friendly to parsing
		my $root = $1 . ":" . $4 . $5;
		$root =~ s/_cmle//g;
		$root = lc($root);
		$roots{$root} = 1;
		$filename2root{$file} = $root;
		$roots_count{$root}++;
		
		# print join("\t", "MWMS", $file, $rename, $root) . "\n";
	} elsif( $file =~ /^(.+_\d\d)_(w[1-6])(DAPI|TRITC|FITC|DIC)(.+)$/)
	{
		# This matches a file with multiple wavelengths
		my $root = $1 . ":" . $4;
		$root =~ s/_cmle//g;
		$root = lc($root);
		$roots{$root} = 2;
		$filename2root{$file} = $root;
		$roots_count{$root}++;
		# print "MW   " . $file . " " . $root . "\n";
	}

}

my $filetext = '';

while(my ($root, $root_type) = each(%roots) )
{
	my $this_dic;
	my $this_dapi;
	my $this_fitc;
	my $this_tritc;
	my @root_files = @{get_files_from_root($root)};
	foreach my $this_file (@root_files)
	{
		if($this_file =~ /w\dDAPI/)
		{
			if(!defined($this_dapi) || $this_file =~ /_cmle/)
			{
				$this_dapi = $this_file
			}
		}

		if($this_file =~ /w\dDIC/)
		{
			if(!defined($this_dic) || $this_file =~ /_cmle/)
			{
				$this_dic = $this_file
			}
		}

		if($this_file =~ /w\dFITC/)
		{
			if(!defined($this_fitc) || $this_file =~ /_cmle/)
			{
				$this_fitc = $this_file
			}
		}

		if($this_file =~ /w\dTRITC/)
		{
			if(!defined($this_tritc) || $this_file =~ /_cmle/)
			{
				$this_tritc = $this_file
			}
		}

	}
	my $macro_string = create_decon_string(-dir=>$indir, -dic=>$this_dic, -fitc=>$this_fitc, -tritc=>$this_tritc, -dapi=>$this_dapi, -root=>$root);

#	print join("\t", $root, "DIC:$this_dic", "FITC:$this_fitc", "TRITC:$this_tritc", "DAPI:$this_dapi" ) . "\n";
#	print join("\t", $root_type, $roots_count{$root}, scalar @root_files, $root, @root_files) . "\n";
	
	$filetext .= $macro_string ;
}


print BATCH getDeconHeader(-resultDir=>$outdir);
print BATCH $filetext;



sub get_files_from_root
{
	my $root = shift;

	my @this_files = ();
	while(my ($filename, $this_root) = each(%filename2root))
	{
		if($root eq $this_root)
		{
			push(@this_files, $filename);
		}
	}
	return \@this_files;
}


sub create_decon_string
{
        my @param = @_;
        my %param = @param;
        @param{ map { lc $_ } keys %param } = values %param; # lowercase key

	my $root = $param{-root};
	my $fitc = $param{-fitc};
	my $tritc = $param{-tritc};
	my $dapi = $param{-dapi};
	my $dic = $param{-dic};

	$root =~ s/\.tif//;
	$root =~ s/\:/_/g;
	if($root =~ /_$/)
	{
		# Do Nothing
	} else
	{
		$root .= "_";
	}

	my $text = '';

	if(defined($dapi))
	{
		$text .= getDeconString(-indir=>$indir, -file=>$dapi, -channel=>"dapi", -excitationWL=>$excitationWL{"dapi"}, -emissionWL=>$emissionWL{"dapi"}, -snr=>$snr{"dapi"} );
	}

	if(defined($fitc))
	{
		$text .= getDeconString(-indir=>$indir, -file=>$fitc, -channel=>"fitc", -excitationWL=>$excitationWL{"fitc"}, -emissionWL=>$emissionWL{"fitc"}, -snr=>$snr{"fitc"} );
	}

	if(defined($tritc))
	{
		$text .= getDeconString(-indir=>$indir, -file=>$tritc, -channel=>"tritc", -excitationWL=>$excitationWL{"tritc"}, -emissionWL=>$emissionWL{"tritc"}, -snr=>$snr{"tritc"} );
	}

	if(defined($dic))
	{
		if($copyDIC)
		{
			copy( File::Spec->catfile($indir, $dic), File::Spec->catfile($outdir,  $dic) );
		}
	}

	return $text;
}


sub varPrint
{
	my $varName = shift;
	my $varValue = shift;
	my $numPrint = shift;

	my $returnText = "";
	$returnText .= " " .  $varName . " {";

	for(my $i=0;$i<$numPrint;$i++)
	{
		$returnText .= $varValue . " " ;
	}

	$returnText .= "} ";

	$returnText .= "parState," . $varName . " {";
	for(my $i=0;$i<$numPrint;$i++)
	{
		$returnText .= "verified " ;
	}
	$returnText .= "} ";

	return $returnText;
}

sub varPrintShort
{
	my $varName = shift;
	my $varValue = shift;
	my $numPrint = shift;

	my $returnText = "";
	$returnText .= " "  . $varName . " {";
	for(my $i=0;$i<$numPrint;$i++)
	{
		$returnText .= $varValue . " " ;
	}

	$returnText .= "} ";


	return $returnText;
}


sub getDeconHeader
{

	# Variables to change

        my @param = @_;
        my %param = @param;
        @param{ map { lc $_ } keys %param } = values %param; # lowercase key

	

	my $thisText = '';

	# This is only done once per batch file
	$thisText .= "info {title {Batch processing template} version 2.2 templateName batch_auto date {Tue Oct 12 17:53:20 PDT 2010}}\n";
	$thisText .= "taskList {setEnv";
	foreach my $task (@taskIDS)
	{
		$thisText .= " taskID:$task";
	}
	$thisText .= "}\n";
	$thisText .= "setEnv {resultDir " . $outdir ." perJobThreadCnt auto concurrentJobCnt 1 exportFormat {type tiff16 multidir 1 cmode scale} timeOut " . $timeOut ."} ";
	$thisText .= "\n";

	return $thisText;
}

sub getDeconString
{
        my @param = @_;
        my %param = @param;
        @param{ map { lc $_ } keys %param } = values %param; # lowercase key

	my $filename = $param{-file};
	my $type = $param{-channel};
	my $excitationWL = $param{-excitationWL};
	my $emissionWL = $param{-emissionWL};
	my $signalToNoise = $param{-snr};



	my $rootDir = $param{-indir};

	my $inputfilename  = $rootDir . "/" . $filename;
	my $outputrootname = $filename;
	$outputrootname =~ s/\.tif$//;

	my $psfFile 	= $psfFiles{$type};



	my $dx_small = sprintf("%.4f", $dx / 1000);
	my $dy_small = sprintf("%.4f", $dy / 1000);
	my $dz_small = sprintf("%.4f", $dz / 1000);




	my $thisText = '';

	# Everything after this is done for each image
	$thisText .= "taskID:" . $taskID . " {info {state readyToRun tag {setp $type cmle Leica_$type} "; # This might not matter much
	$thisText .= "timeStartAbs 1286931200 timeOut " . $timeOut . "} taskList {imgOpen setp cmle imgSave} ";
	$thisText .= "imgOpen {path " . $inputfilename . " series auto index 0} imgSave {rootName " . $outputrootname ."} ";
	$thisText .= "setp {s { " . $dx_small . " " . $dy_small . " " . $dz_small . " " . $dt . "} parState,s verified ";
	$thisText .= "dx " . $dx . " parState,dx verified dy " . $dy . " parState,dy verified dz " . $dz . " parState,dz verified dt " . $dt . " parState,dt verified ";
	$thisText .= "iFacePrim 0.000 parState,iFacePrim verified iFaceScnd 1000000.000 parState,iFaceScnd verified ";


	$thisText .= varPrint("micr", $micType, 32);
	$thisText .= varPrint("na", $micNA, 32);
	$thisText .= varPrint("objQuality", $objQual, 32);
	$thisText .= varPrint("ri", $riMedium, 32);
	$thisText .= varPrint("ril", $riLens, 32);
	$thisText .= varPrint("ps", "2.530", 32);
	$thisText .= varPrint("pr", "280", 32);
	$thisText .= varPrint("ex", $excitationWL, 32); 
	$thisText .= varPrint("em", $emissionWL, 32); 
	$thisText .= varPrint("pcnt", "1", 32);
	$thisText .= varPrint("exBeamFill", "2.00", 32);
	$thisText .= varPrint("imagingDir", "upward", 32);

	$thisText .= " allVerified 1} ";


	$thisText .= " cmle {psfMode file psfPath " . $psfFile ." it " . $iterations . " q " . $qual ." bgMode " . $bgMode . " brMode " . $brMode . " blMode " . $blMode . " pad " . $pad ." mode " . $mode ." timeOut " . $timeOut ;
	$thisText .= varPrintShort("bg", "0.0", 32);
	$thisText .= varPrintShort("sn", "33.0", 32);

	$thisText .= "}} \n";

	push(@taskIDS, $taskID);
	$taskID++;
	return $thisText;
}



