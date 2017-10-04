#!/usr/bin/python

#Should do ng/serialize and than ng/heatmap with export flag

import os
import time
import sys
import ConfigParser

def startStage(number):
	global watch
	watch = time.time()
	print "Running stage %s"%number
def endStage(number):
	global watch
	print "Stage %s done in %f seconds"%(number,time.time()-watch)
def execWithPushd(line,directory):
	os.system("pushd %s > /dev/null; %s; popd > /dev/null"%(directory,line))
def doStage(number,line,directory):
	startStage(number)
	execWithPushd(line,directory)
	endStage(number)
def doStageWithFlag(number,line,directory,flag):
	if flag:
		doStage(number,line,directory)
	else:
		print "Skipping stage %s"%number

def doMatlabStageWithFlag(number,matline,directory,dump,flag):
	line = "matlab -nodisplay -r \"try %s;catch E;fprintf('AAABBBCCCAAA');getReport(E),end;exit\" 1>%s 2>%s"%(matline,dump,dump)
	doStageWithFlag(number,line,directory,flag)
	if flag:
		with open(dump,'rb') as f:
			for ln in f.readlines():
				if 'AAABBBCCCAAA' in ln:
					print "Error occured in stage %s"%number
					sys.exit(-1)

def doMatlabStageWithFlagDisplay(number,matline,directory,dump,flag):
	line = "matlab -nodesktop -nosplash -r \"try %s;catch E;fprintf('AAABBBCCCAAA');getReport(E),end;exit\" 1>%s 2>%s"%(matline,dump,dump)
	doStageWithFlag(number,line,directory,flag)
	if flag:
		with open(dump,'rb') as f:
			for ln in f.readlines():
				if 'AAABBBCCCAAA' in ln:
					print "Error occured in stage %s"%number
					sys.exit(-1)

def chrsizeFromFile(chrname, chrfile):
	with open(chrfile, "rb") as f:
		for line in f.readlines():
			c_name, c_size = line.strip().split("\t",2)
			if c_name == chrname:
				return c_size
	raise Exception("Failed to find chr")

# Assumes input file is CSV: chrN\tSTART\tEND
def fixDomainsFile(inFile, outFile, chrsize):
	with open(inFile, "rb") as fIn, open(outFile, "wb") as fOut:
		badLineEncountered = False
		for line in fIn.readlines():
			ls = line.strip().split("\t")
			assert(not badLineEncountered)
			if len(ls) == 3:
				fOut.write("%s\t%s\n" % (ls[1], ls[2]))
			if len(ls) == 2:
				fOut.write("%s\t%s\n" % (ls[1], chrsize))
				badLineEncountered = True

def makeAbsFile(tmp):
	global output_dir,prefix,chrname
	fName = "%s/%s.%s.%s"%(output_dir,prefix,chrname,tmp)
	fName = os.path.abspath(fName)
	return fName

def config_get_with_default(c, sec, arg, default, boolean=False):
    if c.has_option(sec, arg):
        if boolean:
            return c.getboolean(sec, arg)
        else:
            return c.get(sec, arg)
    else:
        return default

def show_help():
    print
    print("PSYCHIC - finding enhancers using Hi-C data")
    print
    print("Usage: %s <ConfigFile> [<Display>]" % sys.argv[0])
    print
    print
    print("Config file's arguments:")
    print("\tres - resolution in bases")
    print("\twin - window size in bases")
    print("\tchrname - chromosome name")
    print("\tchrsize - path to a chrom.sizes file")
    print("\toutput_prefix - prepend to output files")
    print("\toutput_dir - output directory")
    print("\tinput_matrix - path to Hi-C matrix")
    print("\tgenes_file - path to genes bed file")
    print("Optional arguments:")
    print("\ttad_init - TAD initialization method ('DI', 'InsuationScore' or 'shuffle /path/to/tads.bed')")
    print("\tbilinear - whether to use bilinear regression, true (default) or false")
    print("\tskip_dixon - skip dixon TAD calling and use supplied domains")
    print("\tdomain_path - path to domains file")
    print("\tskip_hierarchy - skip domain merging")
    print
    print("For more information see: https://github.com/dhkron/PSYCHIC") 

#Load config from file
if len(sys.argv) < 2 or len(sys.argv) > 3 or \
   "-h" in sys.argv or "--help" in sys.argv:
    show_help()
    exit()
c = ConfigParser.ConfigParser()
c.read(sys.argv[1])
sec = c.sections()[0]

#Configurable stuff
res = c.getint(sec,'res')
win = c.getint(sec,'win')
chrname = c.get(sec,'chrname')
chrsize = c.get(sec,'chrsize')
prefix = c.get(sec,'output_prefix')
output_dir = c.get(sec,'output_dir')
input_matrix_path = c.get(sec,'input_matrix')
chrnum = chrname.replace('chr','')
fGenes = c.get(sec,'genes_file')
tad_init = config_get_with_default(c, sec, 'tad_init', 'di')
bilinear = config_get_with_default(c, sec, 'bilinear', 'true').lower()

# All paths should be absolute
chrsize = os.path.abspath(chrsize)
output_dir = os.path.abspath(output_dir)
input_matrix_path = os.path.abspath(input_matrix_path)
fGenes = os.path.abspath(fGenes)

skipdixon = config_get_with_default(c, sec, 'skip_dixon', False, True)
domainpath = config_get_with_default(c, sec,'domain_path', '')
skip_hierarchy = config_get_with_default(c, sec, 'skip_hierarchy', False, True)

flgDebugDixon = config_get_with_default(c, sec, 'debug_dixon_domains', False, True)
flgDebugCompare = config_get_with_default(c, sec, 'debug_compare', False, True)

if skip_hierarchy:
    skipdixon = True

#Params for Dixon subchain
hmm_min = 2
hmm_prob = 0.99
if len(sys.argv) >= 3:
	display = sys.argv[2]
else:
	display = False

#Inner parameters
path_to_dixon = "domaincall_software"
path_to_dixon_perl = path_to_dixon + "/perl_scripts"
path_to_matlab = "matlab"
path_to_insulation = "insulation"

#All files should be converted to absolute path, because I'm using pushd and popd
fMatrix = input_matrix_path
fMatrix2 = makeAbsFile("fixed_matrix")
fDI = makeAbsFile("DI")
fHMM = makeAbsFile("HMM")
f7col = makeAbsFile("7col")
fDomains = makeAbsFile("domains")
fMatrixDbg = makeAbsFile("matrix.txt")
fDomainsDbg = makeAbsFile("domains.txt")
fModel = makeAbsFile("model.mat")
fPrTMap = makeAbsFile("prob.tad.matrix.txt")
fPrBMap = makeAbsFile("prob.bg.matrix.txt")
fSupersum = makeAbsFile("supersum.txt")
fLLR = makeAbsFile("llr.txt")
fNewDomains = makeAbsFile("domains.new.txt")
fBed = makeAbsFile("hierarchy.bed")
fFig = makeAbsFile("hierarchy.png")
fBedModelEstimated = makeAbsFile("model.estimated.params.bed")
fMatrixModelEstimated = makeAbsFile("model.estimated.matrix.txt")
fHrrcDebug = makeAbsFile("hierarchy.small.png")
fEnh4 = makeAbsFile("enh_1e-4.bed")
fEnh3 = makeAbsFile("enh_1e-3.bed")
fEnh2 = makeAbsFile("enh_1e-2.bed")
fEnhR = makeAbsFile("enh_rand.bed")
fEnh52 = makeAbsFile("enh_5e-2.bed")

#Check file exist
flgShouldFixMatrix = not os.path.exists(fMatrix2) and not skipdixon
flgDI = not os.path.exists(fDI) and not skipdixon
flgHMM = not os.path.exists(fHMM) and not skipdixon
flg7col = not os.path.exists(f7col) and not skipdixon
flgDomains = not os.path.exists(fDomains) and not skipdixon
flgModel = not os.path.exists(fModel) and not skip_hierarchy
flgProbMatrixes = not (os.path.exists(fPrTMap) and \
                       os.path.exists(fPrBMap) and \
                       os.path.exists(fSupersum) and \
                       os.path.exists(fLLR)) and \
                  not skip_hierarchy
flgStage5 = not (os.path.exists(fMatrixDbg) and os.path.exists(fDomainsDbg)) and not skip_hierarchy
flgNewTads = not (os.path.exists(fMatrixDbg) and os.path.exists(fNewDomains)) and not skip_hierarchy
flgBed = not (os.path.exists(fBed)) and not skip_hierarchy
flgME = not (os.path.exists(fBedModelEstimated))
flgHrrcDebug = not (os.path.exists(fHrrcDebug))
flgEnh = not (os.path.exists(fEnh4) and os.path.exists(fEnh3) and os.path.exists(fEnh2) and os.path.exists(fEnhR) and os.path.exists(fEnh52)) 

#Matrix fixer
if flgShouldFixMatrix:
	with open(fMatrix,"rb") as f:
		try:
			line = f.readline()
			_chr,_start,_end,_rest = line.split("\t",3)
			if "chr" not in _chr:
				print "Creating a Dixon-compatible matrix..."
				#Fix matrix
				with open(fMatrix2,"wb") as f2:
					f2.write("%s\t%d\t%d\t%s"%(chrname,0,res,line))
					counter = 1
					for line in f.readlines():
						f2.write("%s\t%d\t%d\t%s"%(chrname,counter*res,(counter+1)*res,line))
						counter = counter+1
				#File write complete, switch pointers
				fMatrix = fMatrix2
				print "Done!"
		except Exception as e:
			print "Exception occured: %s"%e
			print "Probably this matrix is not in the right format"
			exit()
else:
	print "Not fixing matrix"


if tad_init.lower() == 'di':
    #Stage 1 - ./DI_from_matrix.pl matrix.chrN @res @win @chrsize > DI.chrN
    line_mat_to_di = "perl DI_from_matrix.pl %s %s %s %s > %s"
    line_mat_to_di = line_mat_to_di%(fMatrix,res,win,chrsize,fDI)

    doStageWithFlag("DixonDI",line_mat_to_di,path_to_dixon_perl,flgDI)

    #Stage 2 - matlab -nodisplay -r "HMM_calls DI.chrN HMM.chrN
    matlab_dump = output_dir + "/%s.%s.mdump1"%(prefix,chrname)
    matlab_dump = os.path.abspath(matlab_dump)
    #line_hmm = "matlab -nodisplay -r \"HMM_calls %s %s; exit\" 2>%s 1>%s"
    #line_hmm = line_hmm%(fDI,fHMM,matlab_dump,matlab_dump)
    line_hmm = "HMM_calls %s %s"%(fDI,fHMM)

    doMatlabStageWithFlag("DixonDomainHMMCalling",line_hmm,path_to_dixon,matlab_dump,flgHMM)

    #Stage 3 - perl file_ends_cleaner.pl HMM.chrN DI.chrN | perl converter_7col.pl > 7col.chrN
    stage_line = "perl file_ends_cleaner.pl %s %s | perl converter_7col.pl > %s"
    stage_line = stage_line%(fHMM,fDI,f7col)

    doStageWithFlag("DixonDomains1",stage_line,path_to_dixon_perl,flg7col)

    #Stage 4 - perl hmm_probability_correcter.pl 7col.chrN @hmm_min @hmm_prob @res | perl hmm-state_caller.pl @chrsize @chr | perl hmm-state_domains.pl > domains.chrN
    stage_line = "perl hmm_probablity_correcter.pl %s %s %s %s | perl hmm-state_caller.pl %s %s | perl hmm-state_domains.pl > %s"
    stage_line = stage_line%(f7col,hmm_min,hmm_prob,res,chrsize,chrname,fDomains)

    doStageWithFlag("DixonDomains2",stage_line,path_to_dixon_perl,flgDomains)

elif tad_init.lower() == 'insulationscore':
    # Replaces Stages 1 - 4
    # Crane et al. 2015. http://www.nature.com/nature/journal/v523/n7559/full/nature14450.html#methods
    stage_line = "python insulation_score_tads.py %s %s %d > %s" % (fMatrix, chrname, res, fDomains)
    doStageWithFlag("InsulationScoreDomains", stage_line, path_to_insulation, flgDomains)
elif tad_init.lower().startswith('shuffle'):
    # Replaces Stages 1 - 4
    boundaries_to_shuffle = tad_init[len('shuffle '):]
    if not os.path.exists(boundaries_to_shuffle):
        print('To use random TAD initialization, set in tad_int "shuffle /path/to/original/domains.txt"')
        exit(-1)
    boundaries_to_shuffle = os.path.abspath(boundaries_to_shuffle)
    stage_line = "python shuffle_tads.py %s %s > %s" % (boundaries_to_shuffle, chrname, fDomains)
    doStageWithFlag("RandomDomains", stage_line, path_to_insulation, flgDomains)
else:
    print('tad_init should be configured to either "DI" (default), "InsulationScore" or "shuffle /path/to/original/domains.txt"')
    exit(-1)

#Stage 5 - Generate useful easy loadable files
if not skipdixon:
	#Fix domains file
	fixDomainsFile(fDomains, fDomainsDbg, chrsizeFromFile(chrname, chrsize))
elif not domainpath == "":
	#If skipped dixon and using a specific filename, use it
	fDomainsDbg = domainpath

#Fix matrix file
stage_line2 = "cut -f4- %s > %s"%(fMatrix,fMatrixDbg)
doStageWithFlag("FixMatrixFormat","%s"%(stage_line2),".",flgStage5)

#Optional - DebugDixonDomains(fMatrixDbg, fDomainsDbg, res)
stage_dbg_line = "matlab -nosplash -nodesktop -display %s -r \"DebugDixonDomains %s %s %s\""
stage_dbg_line = stage_dbg_line%(display,fMatrixDbg,fDomainsDbg,res)
doStageWithFlag('domain_debug',stage_dbg_line,path_to_matlab,flgDebugDixon and bool(display))

#Stage 6 - GenerateModelFromTADs(fMatrix,fDomains,res,win,fModel)
matlab_dump = output_dir + "/%s.%s.mdump2"%(prefix,chrname)
matlab_dump = os.path.abspath(matlab_dump)
#	      GenerateModelFromTADs(fMatrix,fDomains,res,win,fOut)
stage_line = "GenerateModelFromTADs %s %s %s %s %s"%(fMatrixDbg,fDomainsDbg,res,win,fModel)
doMatlabStageWithFlag("GenerateProbModel",stage_line,path_to_matlab,matlab_dump,flgModel)

#Stage 7 - GenerateSupersumFromModel(fMatrix,fDomains,res,win,fOutT,foutB,foutR)
matlab_dump = output_dir + "/%s.%s.mdump3"%(prefix,chrname)
matlab_dump = os.path.abspath(matlab_dump)
#	      GenerateSupersumFromModel(fMatrix,fModel,res,win,fOutSupersum,fOutT,fOutB,fOutLLR)
stage_line = "GenerateSupersumFromModel %s %s %s %s %s %s %s %s"
stage_line = stage_line%(fMatrixDbg,fModel,res,win,fSupersum,fPrTMap,fPrBMap,fLLR)
doMatlabStageWithFlag("GenerateSupersum",stage_line,path_to_matlab,matlab_dump,flgProbMatrixes)

#Stage 8 - Create a bed file with TADs of new method
matlab_dump = output_dir + "/%s.%s.mdump4"%(prefix,chrname)
matlab_dump = os.path.abspath(matlab_dump)
#	      FindNewTads(fMat,fSupersum,res,box,fNewDomains)
stage_line = "FindNewTads %s %s %s %s %s"
stage_line = stage_line%(fMatrixDbg,fSupersum,res,"0",fNewDomains)
doMatlabStageWithFlag("FindTADsNewModel",stage_line,path_to_matlab,matlab_dump,flgNewTads)

#Stage 8 - Create a bed files out of all this stuff
# ** DEPRACATED: OLD MERGE METHOD **
# matlab_dump = output_dir + "/%s.%s.mdump4"%(prefix,chrname)
# matlab_dump = os.path.abspath(matlab_dump)
# stage_line = "CreateSingleHierarchyBed %s %s %s %s %s %s %s %s %s"
# stage_line = stage_line%(fMatrixDbg,fLLR,fSupersum,prefix,res,chrnum,"0",fBed,fFig)
# doMatlabStageWithFlag(8,stage_line,path_to_matlab,matlab_dump,flgBed)

#Stage 9 - Create a bed files out of all this stuff
matlab_dump = output_dir + "/%s.%s.mdump5"%(prefix,chrname)
matlab_dump = os.path.abspath(matlab_dump)
#	      CreateSingleHierarchyBedNewMethod(fMat,fNewDomains,fBgModel,prefix,res,chr,box,bedPath,figPath)
stage_line = "CreateSingleHierarchyBedNewMethod %s %s %s %s %s %s %s %s" # %s
stage_line = stage_line%(fMatrixDbg,fNewDomains,fModel,prefix,res,chrnum,"0",fBed) # fFig
doMatlabStageWithFlag("FindHierarchies",stage_line,path_to_matlab,matlab_dump,flgBed)

#Stage 10 - Create model estimated heatmap & BED describing the two power law
matlab_dump = output_dir + "/%s.%s.mdump6"%(prefix,chrname)
matlab_dump = os.path.abspath(matlab_dump)
#	      ModelEstimate(matPath,bedPath,bedPathOut,emMatOut,res)
stage_line = "ModelEstimate %s %s %s %s %d %s"
stage_line = stage_line%(fMatrixDbg,fBed,fBedModelEstimated,fMatrixModelEstimated,res,bilinear)
doMatlabStageWithFlag("ModelEstimate",stage_line,path_to_matlab,matlab_dump,flgME)

#Stage 11 - Create an image of the hierarchies 
matlab_dump = output_dir + "/%s.%s.mdump7"%(prefix,chrname)
matlab_dump = os.path.abspath(matlab_dump)
#	      PlotHrrcBed(fMatrix,fHrrcBed,s,e,res,filter,figPath)
stage_line = "PlotHrrcBed %s %s %d %d %d %d %s"
stage_line = stage_line%(fMatrixDbg,fBed,1,1000,res,0,fHrrcDebug)
# doMatlabStageWithFlagDisplay("DebugHierarchies",stage_line,path_to_matlab,matlab_dump,flgHrrcDebug)

#Stage 12 - Find potential enhancers - overrepresented areas interacting with promotor 
matlab_dump = output_dir + "/%s.%s.mdump8"%(prefix,chrname)
matlab_dump = os.path.abspath(matlab_dump)
#	      EnhancerPromoter(fnij,fme,fgenes,res,ch,out4,out3,out2,outRand)
stage_line = "EnhancerPromoter %s %s %s %d %s %s %s %s %s %s"
stage_line = stage_line%(fMatrixDbg,fMatrixModelEstimated,fGenes,res,chrnum,fEnh4,fEnh3,fEnh2,fEnhR,fEnh52) 
#print stage_line
doMatlabStageWithFlag("EnhancerPromoter",stage_line,path_to_matlab,matlab_dump,flgEnh)
