#!/usr/bin/perl


use Getopt::Long;

# Requires as argument the full path to the files you would like to process

# ./process_ij.pl --indir=/home/mcipriano/Desktop/ij_scripts/test --outdir=/home/mccipriano/Desktop/ij_scripts/test/deconvolved --celldir=/home/mcipriano/Desktop/ij_scripts/test/deconvolved/cells


my $dir = '';
my $outdir = '';
my $celldir = '';

my $opts = GetOptions( 	"indir=s"	=>	\$dir,
			"outdir=s"	=>	\$outdir,
			"celldir=s"	=>	\$celldir,
			"basename=s"	=>	\$basename
			);

if(!defined($basename))
{
	$basename = 'unnamed';
}
opendir(DIR, $dir);

my @files = readdir(DIR);

my $ij_pre_script_name = $basename . "_ij_preprocess.txt";
my $huygens_batch_name = $basename . "_huygens_batch.hgsb";
my $ij_post_script_name = $basename . "_ij_postprocess.txt";

open(IJ, ">$ij_pre_script_name");
open(IJPOST, ">$ij_post_script_name");
open(BATCH, ">$huygens_batch_name");

my %filehash;

foreach my $file (@files)
{
	if(-f $dir . "/" . $file)
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
	my $macro_string = create_macro_string(-dir=>$dir, -dic=>$this_dic, -fitc=>$this_fitc, -tritc=>$this_tritc, -dapi=>$this_dapi, -root=>$root);
#	print join("\t", $root, "DIC:$this_dic", "FITC:$this_fitc", "TRITC:$this_tritc", "DAPI:$this_dapi" ) . "\n";
#	print join("\t", $root_type, $roots_count{$root}, scalar @root_files, $root, @root_files) . "\n";
	my $macro_string = create_macro_string(-dic=>$this_dic, -fitc=>$this_fitc, -tritc=>$this_tritc, -dapi=>$this_dapi, -root=>$root);
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
		warn("Must have DIC/FITC/TRITC/DAPI images: $dic $fitc $tritc $dapi");
		return '';
	}
	$text .= "path = \"$dir\";\n";
	$text .= "pathcell = \"$celldir\";\n";
	$text .= "pathdecon = \"$outdir\";\n";
	$text .= "fileBASE = \"$root\";\n";


	if(defined($fitc))
	{
		$text .= "fileFITC = \"$fitc\";\n";
	}

	if(defined($tritc))
	{
		$text .= "fileTRITC = \"$tritc\";\n";
	}
	if(defined($dapi))
	{
		$text .= "fileDAPI = \"$dapi\";\n";
	}
	if(defined($dic))
	{
		$text .= "fileDIC = \"$dic\";\n";
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


	addAmt = 60;

	for(i=0;i<numROI;i++)
	{
		roiManager("Select", 0);
		getSelectionBounds(thisX[i], thisY[i], thisW[i], thisH[i]);
		roiManager("Delete");

	}

	for(i=0;i<numROI;i++)
	{
		numthisX = thisX[i]-addAmt;
		numthisY = thisY[i]-addAmt;
		numthisW = thisW[i]+2*addAmt;
		numthisH= thisH[i]+2*addAmt;
		if(numthisX < 0)
		{
			numthisX = 0;
		}
		if(numthisY < 0)
		{
			numthisY = 0;
		}

		makeRectangle(numthisX,numthisY, numthisW, numthisH);
		roiManager("Add");
	}

	for(i=0;i<numROI;i++)
	{
		selectWindow(allFile);
		roiManager("Select", i);
		run("Duplicate...", "title=" + targetBASE + "cell" + i + " duplicate");
		saveAs("Tiff", targetPath + File.separator + targetBASE+ "cell" + i + "_STACKS");
	}


}

setBatchMode(1);

startMISlice = 1;
stopMISlice = 0;
contrastValue = 0.2;


// Open Files
if(fileFITC != "*None*")
{
    open(path + File.separator + fileFITC);
    reSlice(startMISlice, stopMISlice);
}
if(fileFITC != "*None*")
{
    open(path + File.separator + fileTRITC);
    reSlice(startMISlice, stopMISlice);
}
if(fileFITC != "*None*")
{
    open(path + File.separator + fileDAPI);
    reSlice(startMISlice, stopMISlice);
}
if(fileFITC != "*None*")
{
    open(path + File.separator + fileDIC);
    reSlice(startMISlice, stopMISlice);
}

// Max Intensity Projection of stacks without DIC
selectWindow(fileFITC);
run("Z Project...", "start=1 stop="+ nSlices + " projection=[Max Intensity]");

selectWindow(fileTRITC);
run("Z Project...", "start=1 stop="+ nSlices + " projection=[Max Intensity]");

selectWindow(fileDAPI);
run("Z Project...", "start=1 stop="+ nSlices + " projection=[Max Intensity]");

run("Merge Channels...", "red=MAX_"+fileTRITC+ " green=MAX_"+ fileFITC + " blue=MAX_"+ fileDAPI + " gray=*None* create keep");

Stack.setChannel(1);
run("Enhance Contrast", "saturated=" + contrastValue);
Stack.setChannel(2);
run("Enhance Contrast", "saturated=" + contrastValue);
Stack.setChannel(3);
run("Enhance Contrast", "saturated=" + contrastValue);

saveAs("Tiff", path + File.separator + fileBASE + "MAX_SINGLE");

// Make Max intensity Montage

selectWindow(fileDIC);
if(nSlices() > 1)
{
	reSlice(round(nSlices()/2), round(nSlices()/2));
}

run("Merge Channels...", "red=MAX_"+fileTRITC+ " green=MAX_"+ fileFITC + " blue=MAX_"+ fileDAPI + " gray=" + fileDIC + " create keep");

Stack.setChannel(1);
run("Enhance Contrast", "saturated=" + contrastValue);
Stack.setChannel(2);
run("Enhance Contrast", "saturated=" + contrastValue);
Stack.setChannel(3);
run("Enhance Contrast", "saturated=" + contrastValue);
Stack.setChannel(4);
run("Enhance Contrast", "saturated=" + contrastValue);

run("Make Montage...", "columns=2 rows=2 scale=0.50 first=1 last=4 increment=1 border=1 font=12");
saveAs("Tiff", path + File.separator + fileBASE + "MAX_MONTAGE");

// Create Hyperstack from all channels

run("Merge Channels...", "red="+fileTRITC+ " green="+ fileFITC + " blue="+ fileDAPI + " gray=*None* create keep");

setSlice(round(nSlices/2));

Stack.setChannel(1);
run("Enhance Contrast", "saturated=" + contrastValue);

Stack.setChannel(2);
run("Enhance Contrast", "saturated=" + contrastValue);

Stack.setChannel(3);
run("Enhance Contrast", "saturated=" + contrastValue);

saveAs("Tiff", path + File.separator + fileBASE + "STACKS");



DAPI_split("MAX_" + fileDAPI, fileBASE + "STACKS.tif", path, fileBASE);



while(nImages() > 0)
{
	close();
}

// End of Macro

ENDEND

	return $macro_string;
}



sub getDeconHeader
{

	my $macro_string;
	$macro_string = <<'ENDEND';

function varPrint (varName, varValue, numPrint)
{
	returnText = "";
	returnText = returnText + " " + varName + " {";
	for(i=0;i<numPrint;i++)
	{
		returnText = returnText + varValue + " " ;
	}

	returnText = returnText + "} ";

	returnText = returnText + "parState," + varName + " {";
	for(i=0;i<numPrint;i++)
	{
		returnText = returnText + "verified " ;
	}
	returnText = returnText + "} ";

	return returnText;
}

function varPrintShort (varName, varValue, numPrint)
{
	returnText = "";
	returnText = returnText + " " + varName + " {";
	for(i=0;i<numPrint;i++)
	{
		returnText = returnText + varValue + " " ;
	}

	returnText = returnText + "} ";

	return returnText;
}


// Variables to change

rootDir = "/home/mcipriano/Desktop/imagej_scripts/test/";


inputfilename  = rootDir + "20101005-Nuf2-GFP_4_5hr_100x_01_w4DAPI_s11_cmle.tif";
outputrootname = rootDir + "20101005-Nuf2-GFP_4_5hr_100x_01_w4DAPI_s11_cmle";

outputfilename = rootDir + "HuygensBatch.hgsb";
psfDir         = rootDir + "psf/";
psfDAPI        = psfDir + "PSF_DAPI.h5";
psfFITC        = psfDir + "PSF_FITC.h5";
psfTRITC       = psfDir + "PSF_TRITC.h5";
psfFile        = psfDAPI;
resultDir      = rootDir + "deconvolved/";

iterations = "100";
qual = "0.01";

bgMode = "auto";
brMode = "auto";
blMode = "auto";
pad = "auto";
mode = "fast";
timeOut = "10000";
bg = "0";

dx = "157.3";
dy = "157.3";
dz = "200";
dTime = "1.000";

micType = "widefield";
riMedium = "1.415";
riLens = "1.515";
objQual = "good";
micNA = "1.40";
excitationWL = "350";
emissionWL = "430";

signalToNoise = "33";


dx_small = d2s( (parseFloat(dx) / 1000), 4); // TODO make 4decimal places 
dy_small = d2s( (parseFloat(dy) / 1000), 4);
dz_small = d2s( (parseFloat(dz) / 1000),4); 



outfile = File.open(outputfilename);

// This is only done once per batch file
print(outfile, "info {title {Batch processing template} version 2.2 templateName batch_2010-10-12 date {Tue Oct 12 17:53:20 PDT 2010}}");
print(outfile, "taskList {setEnv taskID:3}");
print(outfile, "setEnv {resultDir " + resultDir + " perJobThreadCnt auto concurrentJobCnt 1 exportFormat {type tiff16 multidir 1 cmode scale} timeOut " + timeOut + "} ");

ENDEND

	return $macro_string;
}

sub getDeconHeader
{

	my $macro_string;
	$macro_string = <<'ENDEND';


function varPrint (varName, varValue, numPrint)
{
	returnText = "";
	returnText = returnText + " " + varName + " {";
	for(i=0;i<numPrint;i++)
	{
		returnText = returnText + varValue + " " ;
	}

	returnText = returnText + "} ";

	returnText = returnText + "parState," + varName + " {";
	for(i=0;i<numPrint;i++)
	{
		returnText = returnText + "verified " ;
	}
	returnText = returnText + "} ";

	return returnText;
}

function varPrintShort (varName, varValue, numPrint)
{
	returnText = "";
	returnText = returnText + " " + varName + " {";
	for(i=0;i<numPrint;i++)
	{
		returnText = returnText + varValue + " " ;
	}

	returnText = returnText + "} ";

	return returnText;
}


// Variables to change

rootDir = "/home/mcipriano/Desktop/imagej_scripts/test/";


inputfilename  = rootDir + "20101005-Nuf2-GFP_4_5hr_100x_01_w4DAPI_s11_cmle.tif";
outputrootname = rootDir + "20101005-Nuf2-GFP_4_5hr_100x_01_w4DAPI_s11_cmle";

outputfilename = rootDir + "HuygensBatch.hgsb";
psfDir         = rootDir + "psf/";
psfDAPI        = psfDir + "PSF_DAPI.h5";
psfFITC        = psfDir + "PSF_FITC.h5";
psfTRITC       = psfDir + "PSF_TRITC.h5";
psfFile        = psfDAPI;
resultDir      = rootDir + "deconvolved/";

iterations = "100";
qual = "0.01";

bgMode = "auto";
brMode = "auto";
blMode = "auto";
pad = "auto";
mode = "fast";
timeOut = "10000";
bg = "0";

dx = "157.3";
dy = "157.3";
dz = "200";
dTime = "1.000";

micType = "widefield";
riMedium = "1.415";
riLens = "1.515";
objQual = "good";
micNA = "1.40";
excitationWL = "350";
emissionWL = "430";

signalToNoise = "33";


dx_small = d2s( (parseFloat(dx) / 1000), 4); // TODO make 4decimal places 
dy_small = d2s( (parseFloat(dy) / 1000), 4);
dz_small = d2s( (parseFloat(dz) / 1000),4); 



outfile = File.open(outputfilename);

// This is only done once per batch file
print(outfile, "info {title {Batch processing template} version 2.2 templateName batch_2010-10-12 date {Tue Oct 12 17:53:20 PDT 2010}}");
print(outfile, "taskList {setEnv taskID:3}");
print(outfile, "setEnv {resultDir " + resultDir + " perJobThreadCnt auto concurrentJobCnt 1 exportFormat {type tiff16 multidir 1 cmode scale} timeOut " + timeOut + "} ");


// Everything after this is done for each image
thisText = "taskID:3 {info {state readyToRun tag {setp DAPI_200nm_template cmle Leica_DAPI} "; // This might not matter much
thisText = thisText + "timeStartAbs 1286931200 timeOut " + timeOut + "} taskList {imgOpen setp cmle imgSave} ";
thisText = thisText + "imgOpen {path " + inputfilename + " series auto index 0} imgSave {rootName " + outputrootname + "} ";
thisText = thisText + "setp {s { " + dx_small + " " + dy_small + " " + dz_small + " " + dTime + "} parState,s verified ";
thisText = thisText + "dx " + dx + " parState,dx verified dy " + dy + " parState,dy verified dz " + dz + " parState,dz verified dt " + dTime + " parState,dt verified ";
thisText = thisText + "iFacePrim 0.000 parState,iFacePrim verified iFaceScnd 1000000.000 parState,iFaceScnd verified ";


thisText = thisText + varPrint("micr", micType, 32);
thisText = thisText + varPrint("na", micNA, 32);
thisText = thisText + varPrint("objQuality", objQual, 32);
thisText = thisText + varPrint("ri", riMedium, 32);
thisText = thisText + varPrint("ril", riLens, 32);
thisText = thisText + varPrint("ps", "2.530", 32);
thisText = thisText + varPrint("pr", "280", 32);
thisText = thisText + varPrint("ex", excitationWL, 32); // Excitation Wavelength -- Change this
thisText = thisText + varPrint("em", emissionWL, 32); // Emission Wavelength -- Change this
thisText = thisText + varPrint("pcnt", "1", 32);
thisText = thisText + varPrint("exBeamFill", "2.00", 32);
thisText = thisText + varPrint("imagingDir", "upward", 32);

thisText = thisText + " allVerified 1} ";


thisText = thisText + " cmle {psfMode file psfPath " + psfFile + " it " + iterations + " q " + qual + " bgMode " + bgMode + " brMode " + brMode + " blMode " + blMode + " pad " + pad + " mode " + mode + " timeOut " + timeOut ;
thisText = thisText + varPrintShort("bg", "0.0", 32);
thisText = thisText + varPrintShort("sn", "33.0", 32);

thisText = thisText + "}} \n";

print(outfile, thisText)

ENDEND

	return $macro_string;
}
