#!/usr/bin/perl


use Getopt::Long;
use File::Copy;
use File::Path;
use File::Spec;

use strict;

# Requires as argument the full path to the files you would like to process

# ./process_ij.pl --indir=/home/mcipriano/Desktop/ij_scripts/test --outdir=/home/mcipriano/Desktop/ij_scripts/test/deconvolved --celldir=/home/mcipriano/Desktop/ij_scripts/test/deconvolved/cells


my $indir; 
my $outdir;
my $celldir;
my $basename;

my $opts = GetOptions( 	"indir=s"	=>	\$indir,
			"outdir=s"	=>	\$outdir,
			"celldir=s"	=>	\$celldir,
			"basename=s"	=>	\$basename
			);

if(!defined($indir))
{
	die("No --indir defined");
}

if(!defined($outdir))
{
	$outdir = File::Spec->catdir($indir, "processed");
	mkpath($outdir);
}
if(!defined($celldir))
{
	$celldir = File::Spec->catdir($outdir, "cells");
	mkpath($celldir);

}
if(!defined($basename))
{
	$basename = 'unnamed';
}

opendir(DIR, $indir);

my @files = readdir(DIR);

my $ij_pre_script_name = File::Spec->catfile($indir, $basename . "_ij_process.txt");

open(IJ, ">$ij_pre_script_name");

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

my $macro_text = '';

while(my ($file, undef) = each(%filehash))
{
	if( !($file =~ /\.tif+$/i))
	{
		next;
	}
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
#	print join("\t", $root, "DIC:$this_dic", "FITC:$this_fitc", "TRITC:$this_tritc", "DAPI:$this_dapi" ) . "\n";
#	print join("\t", $root_type, $roots_count{$root}, scalar @root_files, $root, @root_files) . "\n";
	my $macro_string = create_macro_string(-dir=>$indir, -dic=>$this_dic, -fitc=>$this_fitc, -tritc=>$this_tritc, -dapi=>$this_dapi, -root=>$root);
	print IJ $macro_string;
}





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


sub create_macro_string
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

	if(!defined($dic) || !defined($fitc) || !defined($tritc) || !defined($dapi) )
	{
		warn("Not all of  DIC/FITC/TRITC/DAPI images: ($root) $dic $fitc $tritc $dapi");
	}

	$text .= "path = \"$indir\";\n";
	$text .= "pathout = \"$outdir\";\n";
	$text .= "pathcell = \"$celldir\";\n";

	$text .= "fileBASE = \"$root\";\n";


	if(defined($fitc))
	{
		$text .= "fileFITC = \"$fitc\";\n";
	} else
	{
		$text .= "fileFITC = \"*None*\";\n";
	}

	if(defined($tritc))
	{
		$text .= "fileTRITC = \"$tritc\";\n";
	} else
	{
		$text .= "fileTRITC = \"*None*\";\n";
	}

	if(defined($dapi))
	{
		$text .= "fileDAPI = \"$dapi\";\n";

	} else
	{
		$text .= "fileDAPI = \"*None*\";\n";
	}

	if(defined($dic))
	{
		$text .= "fileDIC = \"$dic\";\n";
	} else
	{
		$text .= "fileDIC = \"*None*\";\n";
	}

	return $text . get_all_macro() ;

}

sub get_all_macro
{

	my $macro_string;
	$macro_string = <<'ENDEND';

// Start of Macro

function reSlice(start_slice, end_slice) 
{
    num_slices = nSlices;
    setSlice(1);
    i = 1;

    while(i < start_slice && nSlices > 1)
    {
        setSlice(1);
        run("Delete Slice");
        i++;
    }


    while(end_slice > 0 && nSlices > 1)
    {
        setSlice(nSlices);
        run("Delete Slice");
        end_slice--;
    }
    return 0;
}

function DAPI_split(dapiChannel, allFile, targetPath, targetBASE)
{

	selectWindow(dapiChannel);

	thisWidth = getWidth();
	thisHeight = getHeight();

	// Clear ROIs
	initROI = roiManager("count");
	for(i=0;i<initROI;i++)
	{
		roiManager("Select", 0);
		roiManager("Delete");
	}

	setAutoThreshold("Minimum dark");
	run("Convert to Mask");
	run("Options...", "iterations=8 edm=Overwrite count=1");
	run("Dilate");
	run("Options...", "iterations=25 edm=Overwrite count=1");
	run("Close-");
	run("Options...", "iterations=1 edm=Overwrite count=1");
	run("Erode");
	run("Create Selection");
	roiManager("Split");

	numROI = roiManager("count");
	thisX = newArray(numROI);
	thisY = newArray(numROI);
	thisW= newArray(numROI);
	thisH= newArray(numROI);


	addAmt = 85;

	for(i=0;i<numROI;i++)
	{
		roiManager("Select", 0);
		getSelectionBounds(thisX[i], thisY[i], thisW[i], thisH[i]);
		roiManager("Delete");

	}

	for(i=0;i<numROI;i++)
	{
		numthisX = thisX[i]- addAmt;
		numthisY = thisY[i]- addAmt;
		numthisdX = thisX[i] + addAmt;
		numthisdY = thisY[i] + addAmt;

		if(numthisX < 0)
		{
			numthisX = 0;
		}
		if(numthisY < 0)
		{
			numthisY = 0;
		}
		if(numthisdY >= thisHeight)
		{
			numthisdY = thisHeight-1;
		}
		if(numthisdX >= thisWidth)
		{
			numthisdX = thisWidth -1;
		}
		numthisW = numthisdX - numthisX;
		numthisH = numthisdY - numthisY;

		makeRectangle(numthisX,numthisY, numthisW, numthisH);
		roiManager("Add");
	}

	for(i=0;i<numROI;i++)
	{
		selectWindow(allFile);
		roiManager("Select", i);
		run("Duplicate...", "title=" + targetBASE + "cell" + i + " duplicate");
		saveAs("Tiff", targetPath + File.separator + targetBASE + "cell" + i + "_STACKS");
	}


}

setBatchMode(1);

startMISlice = 1;
stopMISlice = 0;
contrastValue = 0.2;

redMAX = "*None*";
blueMAX = "*None*";
greenMAX = "*None*";

numChannels = 0; // Does not count DIC channel

// Open Files
if(fileFITC != "*None*")
{
    open(path + File.separator + fileFITC);
    reSlice(startMISlice, stopMISlice);
    selectWindow(fileFITC);
    run("Z Project...", "start=1 stop="+ nSlices + " projection=[Max Intensity]");
    greenMAX = "MAX_" + fileFITC;
    numChannels++;
}
if(fileTRITC != "*None*")
{
    open(path + File.separator + fileTRITC);
    reSlice(startMISlice, stopMISlice);
    selectWindow(fileTRITC);
    run("Z Project...", "start=1 stop="+ nSlices + " projection=[Max Intensity]");
    redMAX = "MAX_" + fileTRITC;
    numChannels++;
}
if(fileDAPI != "*None*")
{
    open(path + File.separator + fileDAPI);
    reSlice(startMISlice, stopMISlice);
    selectWindow(fileDAPI);
    run("Z Project...", "start=1 stop="+ nSlices + " projection=[Max Intensity]");
    blueMAX = "MAX_" + fileDAPI;
    numChannels++;
}

if(fileDIC != "*None*")
{
    open(path + File.separator + fileDIC);
    reSlice(startMISlice, stopMISlice);
}

// Max Intensity Projection of stacks without DIC

run("Merge Channels...", "red=" + redMAX + " green=" + greenMAX + " blue=" + blueMAX + " gray=*None* create keep");

for(i=0;i<numChannels;i++)
{
    Stack.setChannel(i);
    run("Enhance Contrast", "saturated=" + contrastValue);
}

saveAs("Tiff", pathout + File.separator + fileBASE + "MAX_OVERLAY");

// Make Max intensity Montage

if(fileDIC != "*None*")
{
    selectWindow(fileDIC);
    if(nSlices() > 1)
    {
    	reSlice(round(nSlices()/2), round(nSlices()/2));
    }

    run("Merge Channels...", "red=" + redMAX + " green=" + greenMAX + " blue=" + blueMAX + " gray=" + fileDIC + " create keep");

    for(i=0;i<numChannels+1;i++)
    {
        Stack.setChannel(i);
        run("Enhance Contrast", "saturated=" + contrastValue);
    }



    run("Make Montage...", "columns=2 rows=2 scale=0.50 first=1 last=4 increment=1 border=1 font=12");
    saveAs("Tiff", pathout + File.separator + fileBASE + "MAX_MONTAGE");
} else
{
    run("Make Montage...", "columns=2 rows=2 scale=0.50 first=1 last=4 increment=1 border=1 font=12");
    saveAs("Tiff", pathout + File.separator + fileBASE + "MAX_MONTAGE");
}

// Create Hyperstack from all channels

run("Merge Channels...", "red="+fileTRITC+ " green="+ fileFITC + " blue="+ fileDAPI + " gray=*None* create keep");

setSlice(round(nSlices/2));

for(i=0;i<numChannels;i++)
{
    Stack.setChannel(i);
    run("Enhance Contrast", "saturated=" + contrastValue);
}

saveAs("Tiff", pathout + File.separator + fileBASE + "STACKS");

if(fileDAPI != "*None*")
{
    stacksFILE = fileBASE + "STACKS.tif";
    DAPI_split("MAX_" + fileDAPI, stacksFILE , pathcell, fileBASE);
}



while(nImages() > 0)
{
	close();
}

// End of Macro ////////////////////////////////////////////////

ENDEND

	return $macro_string;
}

