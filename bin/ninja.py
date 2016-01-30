#!/usr/bin/env python
import timeit
import argparse
import os
import subprocess
from subprocess import Popen, PIPE
import sys
import shutil
__version__ = "Ninja 1.3.0"

###
#   CLASSES
###

class Logger(object):
  """ A convenient logging object
      Prints output to a given log file and/or stdout
  """
  def __init__(self, logfp=None, use_std_out=True):
    # note: if logfp directories don't exist, make them.

    if logfp is not None:
      outdir = os.path.abspath(os.path.dirname(os.path.realpath(logfp)))

      # Checks for output directory. Makes it if necessary. 
      if not os.path.exists(outdir):
        os.makedirs(outdir)
    self.logfp = logfp
    logf = open(logfp,'w')
    logf.close()
    self.use_std_out = use_std_out

  def log(self, msg):
    if not msg.endswith('\n'):
      msg += '\n'        
    if self.logfp is not None:
      logf = open(self.logfp,'a')
      logf.write(msg)
      logf.close()
    if self.use_std_out:
      sys.stdout.write(msg)

###
#   UTILITY METHODS
###

# Gets arguments from command line using argparse, instead of the deprecated optparse
# Takes an argparse object as a parameter
def get_args(p):
    p.add_argument("-i", "--input",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input fasta file (use \"input1.fna,input2.fna\" for paired-end) [required].")
    p.add_argument("-o", "--output",
                   type = str,
                   default = '.',
                   metavar = '',
                   help = "Output folder (default current directory)")
    p.add_argument("-b", "--database",
                   type = str,
                   default = 'greengenes97',
                   metavar = '',
                   help = "Name of database folder in ninja top-level directory; folder must contain bowtie2 index with basename the same as the folder, a taxonomy file named basename.taxonomy, and a ninja db map (output from ninja_prep) named basename.db [default greengenes97]")
    p.add_argument("-t", "--trim",
                   type = int,
                   default = -1,
                   metavar = '',
                   help = "Trim sequences to a specified number of bp, cutting off from their ends (default no trim)")
    p.add_argument("-T", "--trim2",
                   type = int,
                   default = -1,
                   metavar = '',
                   help = "Trim reverse paired-end reads to a different length than forward reads. If trim2 is specified, trim must be too. (default no trim)")
    p.add_argument("-s", "--similarity",
                   type = float,
                   default = 0.97,
                   metavar = '',
                   help = "Minimum percent similarity - id - between query sequence and reference sequence (default 97%%)")
    p.add_argument("-I", "--insert",
                   type = int,
                   default = 1600,
                   metavar = '',
                   help = "Maximum total length for paired-end matches. Set this as small as possible (e.g. 400 for 515F-806R primers) (default 1600, for 16S)")
    p.add_argument("-p", "--threads",
                   type = int,
                   default = 4,
                   metavar = '',
                   help = "Number of threads/cores to run bowtie2 on (default 4)")
    p.add_argument("-m", "--mode",
                   type = str,
                   default = 'normal',
                   metavar = '',
                   help = "NINJA sensitivity mode: 'normal' (default), 'fast' (less sensitive), or 'max' (more sensitive, slower)")
    p.add_argument("-d", "--denoising",
                   type = int,
                   default = 1,
                   help = "Discards all reads that appear fewer than d times. No denoising = 1; " + \
                          " Moderate denoising = 2 (throws out all singleton reads;" + \
                          " Aggressive denoising = 6 (nearly guaranteed to eliminate all sequencing error - although not PCR error - in most data sets) (default 2)")
    p.add_argument("-F", "--full_output",
                   action = 'store_true',
                   help = "Output fasta files containing failed sequences and filtered sequences [default False]")
    p.add_argument("-l", "--legacy",
                   action = 'store_true',
                   help = "Output legacy (tab-delimited) QIIME OTU table [default False]")
    p.add_argument("-P", "--print_only",
                   action = 'store_true',
                   help = "Print commands only - do not run [default False]")
    p.add_argument("-S", "--suppress_stdout",
                   action = 'store_true',
                   help = "Suppress standard output [default False]")
    p.add_argument("-R", "--retain_intermediates",
                   action = 'store_true',
                   help = "Retain intermediate files [default False]")
    p.add_argument("-C", "--check_fasta",
                   action = 'store_true',
                   help = "Check fasta for correct formatting; otherwise assumes fasta is in QIIME-ready format [default False]")
    p.add_argument("-r", "--reverse_complement",
                   action = 'store_true',
                   help = "Flags sequences for reverse complementing (default no reverse complement)")

    args = p.parse_args()
    return args

# Checks if args work. Takes args and parser as input
def check_args(args, p):

    if args['input'] is None:
        p.print_help()
        sys.exit('\nPlease include an input sequences file in fasta format.')
    if args['trim'] == -1 and args['trim2'] != -1:
        sys.exit('\nIf trim2 is specified you must also specify trim1 (set to large number for no trimming).')
    

# Checks if an input fasta file is formatted correctly for QIIME. Returns boolean.
# Looks for data or titles with multiple lines and improper characters in seqs.
# Makes sure first line is a header. Doesn't allow spaces in seqs.
def check_fasta(f, logger):
    lineNumber = 1
    inData = True
    seqChars = ['A', 'T', 'G', 'C', 'N', '\n']
    for line in f:
        if line[0] == ">":
            if inData:
                inData = False
            else:
                logger.log("Warning: Multiline title in line " + str(lineNumber) + " of input fasta.")
                return False
        else:
            # Checks if line is actually a title
            for ch in line:
                c = ch.upper()
                if c not in seqChars:
                    logger.log('Warning: Unsupported sequence character "' + c + '" in line ' + str(lineNumber) + ' of input fasta.')
                    return False
            if not inData:
                inData = True
            else:
                logger.log("Warning: Multiline sequence in line " + str(lineNumber) + " of input fasta.")
                return False
        lineNumber += 1
    return True

# Generator for fasta files - returns [(name, seq)]
# Call using 'with open(file) as f'
def read_fasta(f):
    title = None
    data = None
    for line in f:
        if line.startswith(">"):
            if title != None:
                yield (title,data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if title is None:
      yield None
    else:
      yield (title,data)
    return

# Takes a list of [(title, seq), ...] tuples and prints as a FASTA file to the given filename
def write_fasta(listOfTuples, fileName):
    output = ""
    if fileName.find(".",0) == -1:
        fileName += ".fna"
    for data in listOfTuples:
        title = data[0].strip()
        if title[0] != ">":
            title = ">" + title
        output += (title + "\n")
        output += (data[1] + "\n")
    with open(fileName, 'wb') as outFile:
        outFile.write(output)

# Takes a dict[OTU, IDs...] and prints as an OTU map to the given filename
def write_map(otuMap, fileName):
    with open(fileName, 'wb') as outFile:
        for k in otuMap:
            outFile.write(k + "\t" + otuMap[k] + "\n")

 # Returns reverse complement of a DNA seq
def reverse_complement(seq): return complement(seq[::-1])

 # Returns complement of a DNA seq
def complement(seq):
    complement = []
    for c in seq:
        if c == 'A':
            complement.append('T')
        elif c == 'T':
            complement.append('A')
        elif c == 'C':
            complement.append('G')
        elif c == 'G':
            complement.append('C')
    return ''.join(complement)

###
#   MAIN METHODS
###

# Runs ninja_filter
# INPUT     inputSeqsFile:          input sequences in fasta format
#           inputSeqsFile2:         paired-end input sequences in fasta format
# OUTPUT    filteredSeqsFile:       input sequences filtered using ninja algorithm
#           seqsDB:                 utility file used in ninja_parse
# OPTIONAL  trim:                   trims sequences to <= uniform X bp (e.g. AGGC, GCG with trim 2 returns AG, GC)
#           RC:                     takes reverse complement of input sequences
#           denoising:              for argument x.y, discards all reads that appear less than x times and all kmers that
#                                   appear less than y times unless kmers in reads that appear more than x times
def ninja_filter(inputSeqsFile, inputSeqsFile2, file_prefix, trim, trim2, RC, denoising, logger, full_output=False,
                run_with_shell=True, print_only=False):

    # Converts optional arguments to args readable by ninja_filter
    argTrim = ''
    argRC = ''
    argDenoising = ''
    if trim != -1:
        argTrim = str(trim)
    if trim2 != -1:
        argTrim += "," + str(trim2)
        argTrim = '"' + argTrim + '"'
    if RC:
        argRC = "RC"

    # Sets the relevant binaries for mac and windows support.
    ninjaDirectory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    ninjaDirectory = os.path.abspath(os.path.join(ninjaDirectory, os.pardir))
    if sys.platform.startswith("darwin") or sys.platform.startswith("os"):      # Mac
      ninjaFilterFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_filter_mac"))
    elif sys.platform.startswith("win32") or sys.platform.startswith("cygwin"):   # Windows and cygwin
      ninjaFilterFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_filter.exe"))
    else:   # Linux
        ninjaFilterFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_filter_linux"))

    argDenoising = 'D ' + str(denoising)

    # Runs ninja_filter. Run in shell only on Mac
    cmd = ""
    cmd = '"' + ninjaFilterFile + '"'
    cmd += ' ' + '"' + inputSeqsFile + '"'
    if inputSeqsFile2 is not None:
      cmd += ' ' + 'PE "' + inputSeqsFile2 + '"'
    cmd += ' ' + '"' + file_prefix + '"'
    cmd += ' ' + argTrim
    cmd += ' ' + argRC
    cmd += ' ' + argDenoising
    if full_output:
      cmd += ' LOG'     
    logger.log(cmd)

    if not print_only:
      proc = Popen(cmd,shell=run_with_shell,universal_newlines=True,stdout=PIPE,stderr=PIPE)
      stdout, stderr = proc.communicate()
      logger.log(stdout)
      if proc.returncode != 0:
        logger.log(stderr + '\n')
        raise ValueError("ERROR: Filtering failed. One possible explanation is a problem with input FASTA formatting. Please rerun with '--check_fasta'. Exiting.")
    return cmd
      

# Runs bowtie2. Uses two presets for ninja normal and max
# INPUT     filteredSeqsFile:   filtered sequences output from ninja_filter
# OUTPUT    alignmentsFile:     main output of bowtie
# OPTIONAL  mode:               'ninja' or 'ninjaMax', for less and more sensitivity, respectively
#           threads:            number of threads/cores to run bowtie2 on
#           similarity:         minimum percent similarity between query sequence and reference sequence
def bowtie2(filteredSeqsFile, filteredSeqsFile2, alignmentsFile, bowtieDatabase, similarity, insert, threads, mode,
            logger, run_with_shell=True, print_only=False):

    # TODO: Automatically convert fasta file if formatted incorrectly

    # Checks if similarity is a percentage
    if similarity > 1 or similarity < 0:
        print("Similarity error. Enter similarity as a percentage between 0 and 100. Exiting.")
        sys.exit()

    similarity = 1 - similarity     # Converts to similarity readable by bowtie2

    # Switches between ninja normal and max according to user input. Only runs in shell on Mac/linux

    cmd = ['bowtie2-align-s','--no-head']
    cmd.append('-x ' + bowtieDatabase)
    cmd.append('-S ' + '"' + alignmentsFile + '"')
    cmd.append('--np 0')
    cmd.append('--mp "1,1"')
    cmd.append('--rdg "0,1"')
    cmd.append('--rfg "0,1"')
    cmd.append('--score-min "L,0,-' + str(similarity) + '"')
    cmd.append('--norc')
    if filteredSeqsFile2 is None:
      cmd.append('-f ' + '"'+ filteredSeqsFile + '"')
    else:
      cmd.append('-f -1 ' + '"'+ filteredSeqsFile + '"')
      cmd.append('-2 ' + '"'+ filteredSeqsFile2 + '"')
      cmd.append('--maxins ' + str(insert))
      cmd.append('--no-mixed --no-discordant')
    cmd.append('-p ' + str(threads))
    cmd.append('-k 1')

    if mode == 'fast':
      cmd.append('--fast')
    elif mode == 'max':
      cmd.append('--very-sensitive')

    # run command
    cmd = ' '.join(cmd)
    logger.log(cmd)
    if not print_only:
      proc = Popen(cmd,shell=run_with_shell,universal_newlines=True,stdout=PIPE,stderr=PIPE)
      stdout, stderr = proc.communicate()
      logger.log(stdout)
      if proc.returncode != 0:
        logger.log(stderr + '\n')
        raise ValueError("ERROR: Bowtie2 failed. Exiting.")
    return cmd

# Runs ninja_parse_filtered.
# INPUT     seqsDBFile:       db file output from ninja_filter
#           alignmentsFile:   alignment file output from bowtie2
#           masterDBFile:     master db file packaged with ninja
#           taxMapFile:       reference taxonomy map packaged with ninja
# OUTPUT    otuTableFile:     OTU table output that's really the point of all this
# AUTO      legacyTable:      OTU table in legacy format, output automatically
#           parseLog:         a utility file containing parsed sequences, used in post-processing and output automatically
def ninja_parse(file_prefix, alignmentsFile, masterDBFile, taxMapFile, full_output,
        logger, legacy_table=False, run_with_shell=True, print_only=False):
  # Sets the relevant binaries for mac and windows support.
  ninjaDirectory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
  ninjaDirectory = os.path.abspath(os.path.join(ninjaDirectory, os.pardir))
  if sys.platform.startswith("darwin") or sys.platform.startswith("os"):      # Mac
    ninjaParseFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_parse_filtered_mac"))
  elif sys.platform.startswith("win32") or sys.platform.startswith("cygwin"):   # Windows and cygwin
    ninjaParseFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_parse_filtered.exe"))
  else:   # Linux
      ninjaParseFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_parse_filtered_linux"))
            
  cmd = ['"' + ninjaParseFile + '"', '"' + file_prefix + '"', '"' + alignmentsFile + '"', '"' + masterDBFile + '"' , '"' + taxMapFile + '"' ]
  if legacy_table:
    cmd.append('--legacy')
  if full_output:
    cmd.append('LOG')
  cmd = ' '.join(cmd)
  logger.log(cmd)
  if not print_only:
    proc = Popen(cmd,shell=run_with_shell,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    stdout, stderr = proc.communicate()
    logger.log(stdout)
    if proc.returncode != 0:
      logger.log(stderr + '\n')
      raise ValueError("ERROR: Parsing failed. Exiting.")
  return cmd


# Performs housekeeping on files, deleting the intermediate ones listed below
# INPUT     inputSeqsFile:      original file of sequences passed to ninja_filter (.fna, .fasta)
#           filteredSeqsFile:   seqs with duplicates removed output from ninja_filter
#           filteredSeqsFile2:  paired-end reverse reads, or None
#           seqsDBFile:         db file output from ninja_filter
#           alignmentsFile:     main output of bowtie2
#           parseLogFile:       parse log generated from ninja_parse_filter
def clean(inputSeqsFile, filteredSeqsFile, filteredSeqsFile2, seqsDBFile, alignmentsFile, dupes_file=None, parseLogFile=None):
    try:
        os.remove(filteredSeqsFile)
        if filteredSeqsFile2 is not None:
          os.remove(filteredSeqsFile2)
        os.remove(seqsDBFile)
        os.remove(alignmentsFile)
        if dupes_file is not None:
          os.remove(dupes_file)
        if parseLogFile is not None:
            os.remove(parseLogFile)
            os.remove("map_seqid_reps.txt")

    except OSError as e:
        error(e, msg = "INTERNAL ERROR: Can't find all files marked for moving and/or deletion. Check working directory and output folder.")

# Runs ninja, bowtie2 and then processes output. All files output in specified output folder. 
# User must specify ninja's directory as an environment variable named 'NINJA_DIR'
def main(argparser):

    args = get_args(argparser)
    args = vars(args)

    # Opens logger to write to log and/or stdout
    # First stores original console location as a variable for error handling
    ninjaLog = os.path.join(args['output'], "ninja_log.txt")
    logger = Logger(ninjaLog, not args['suppress_stdout'])
    logger.log(__version__)

    check_args(args, argparser)

    # if paired end, store second file as a different parameter
    if("," in args['input']):
        files = args['input'].split(',')
        args['input'] = files[0]
        args['input2'] = files[1]
    else:
        args['input2'] = None

    # Checks if input sequences fasta is correctly formatted. Writes correct one if not
    if args['check_fasta']:
        fileName = args['input']
        if not check_fasta(open(fileName), logger):
            new_input_fasta = os.path.join(args['output'], "formatted_input_fasta.fna")
            logger.log("Warning: Input fasta formatted incorrectly for QIIME, e.g. sequences or title on multiple lines. Writing " + \
                  "corrected file to " + new_input_fasta)
            with open(fileName) as f:
                write_fasta(read_fasta(f), new_input_fasta)
                args['input'] = new_input_fasta
        if args['input2'] is not None:
          fileName = args['input2']
          if not check_fasta(open(fileName), logger):
              new_input_fasta = os.path.join(args['output'], "formatted_input2_fasta.fna")
              logger.log("Warning: Reverse input fasta formatted incorrectly for QIIME, e.g. sequences or title on multiple lines. Writing " + \
                    "corrected file to " + new_input_fasta)
              with open(fileName) as f:
                  write_fasta(read_fasta(f), new_input_fasta)
                  args['input2'] = new_input_fasta

    RC = args['reverse_complement']
    similarity = args['similarity']
    threads = args['threads']
    mode = args['mode']
    denoising = args['denoising']
    suppress_stdout = args['suppress_stdout']
    full_output = args['full_output']
    retain_intermediates = args['retain_intermediates']
    legacy_table = args['legacy']

    # Runs ninja pipeline

    # Gets ninja's directory relative to current working directory
    ninjaDirectory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    ninjaDirectory = os.path.abspath(os.path.join(ninjaDirectory, os.pardir))

    # Checks for output subdirectory of current working directory. Makes it if necessary. 
    # Edits global output folder variable
    outdir = os.path.join(os.getcwd(), args['output'])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outdir = os.path.abspath(outdir)
        
    # do not run commands with shell=TRUE on win32 or cygwin
    run_with_shell = not (sys.platform.startswith("win32") or sys.platform.startswith("cygwin"))

    try:
      bowtie2_cmd = "bowtie2-align-s"
      subprocess.check_call(bowtie2_cmd + " --version", shell=run_with_shell, stdout = sys.stdout)
    except OSError as e:
      try:
        bowtie2_cmd = os.path.join(ninjaDirectory,"bowtie2-align-s")
        subprocess.check_call(bowtie2_cmd + " --version", shell=run_with_shell, stdout = sys.stdout)
      except OSError as e:
        error(e = None, msg = "ERROR: Bowtie2 executable not found in system path or top-level NINJA package folder. Please install bowtie2 and add its accompanying executables to the system path or place bowtie2-align-s in the top-level ninja package folder (not a subfolder). " + \
                            "Check README.txt for additional instructions. Exiting.", exit = True)

    # Sets variables used in ninja calls. First, ninja_filter files
    file_prefix = os.path.join(outdir, "ninja")

    # Set paired-end file to None if this is not a paired-end run
    if(args['input2']) is not None:
      pe_file = file_prefix + "2_filt.fa"
    else:
      pe_file = None 

    # Bowtie2 files
    alignmentsFile = os.path.join(outdir, "alignments.txt")
    databasedir = os.path.join(ninjaDirectory, 'databases', args['database'])
    logger.log('Ninja database directory is is' + databasedir)

    masterDBFile = os.path.abspath(os.path.join(databasedir, args['database'] + ".db"))
    bowtieDatabase = os.path.abspath(os.path.join(databasedir, args['database']))
        # Ninja_parse files
    taxMapFile = os.path.abspath(os.path.join(databasedir, args['database'] + ".taxonomy"))
    otuTableFile = os.path.join(outdir, "otutable.biom")
        # Post-processing files
    seqOutFile = os.path.join(outdir, "failed_sequences.fna")
    mapOutFile = os.path.join(outdir, "otu_map.txt")

    # Runs ninja_filter, bowtie2 and ninja_parse. Processes ninja results, generating OTU map and a list of failed seqs
    logger.log("Running Ninja filter...")
    t1 = timeit.Timer(lambda:
      ninja_filter(args['input'], args['input2'], file_prefix, args['trim'], args['trim2'], RC, denoising, logger, full_output,
        run_with_shell=run_with_shell, print_only=args['print_only'])
    )
    logger.log("Ninja filter time: " + str(t1.timeit(1)))
    logger.log("Running Bowtie...")
    t2 = timeit.Timer(lambda:
      bowtie2(file_prefix + "_filt.fa", pe_file, alignmentsFile, bowtieDatabase, similarity, args['insert'], threads, mode,
        logger, run_with_shell=run_with_shell, print_only=args['print_only'])
    )
    logger.log("Bowtie time: " + str(t2.timeit(1)))
    logger.log("Running Ninja parse...")
    t3 = timeit.Timer(lambda:
      ninja_parse(file_prefix, alignmentsFile, masterDBFile, taxMapFile, full_output, 
          logger, legacy_table, run_with_shell=run_with_shell, print_only=args['print_only'])
    )
    logger.log("Ninja parse time: " + str(t3.timeit(1)) + "\n")

    if not retain_intermediates:
      clean(args['input'], file_prefix + "_filt.fa", pe_file, file_prefix + ".db", alignmentsFile)

# Wrapper for main function, called from command line
# Bare minimum args:
#   -i "seqs.fna" 
# Sample maximum args:
#   -i "seqs.fna" -o "output" -r -t 200 -mo 'max' -s 98 -d 1.005 -q
if __name__=='__main__':
    # Parses command line arguments
    p = argparse.ArgumentParser(description = "NINJA-OPS: NINJA Is Not Just Another OTU Picking Solution (v1.1)\n" + \
                                              "Knights Lab (www.ninja-ops.ninja)\n" + \
                                              "This program outputs an otu table and map from sequence reads in fasta format.", 
                                add_help = True, 
                                epilog ='NOTE: If one or more output files are empty, trying reverse complementing your input ' + \
                                        'sequences with -r')

    main(p)
