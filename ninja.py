import timeit
import argparse
import os
import subprocess
import sys
import shutil

###
#   GLOBAL VARIABLES
###

# Stores original console location and log file location 
console = None
ninjaLog = ""

# Stores output folder and ninja executables as global variables
out = ""
ninjaFilterFile = ""
bowtie2File = ""
ninjaParseFile = ""

# Stores booleans for subprocess shell args, depending on the OS, verbose output and input fasta correction 
shellBool = False
verbose = True
inputFastaCorrected = False


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
                   help="Input sequences originally passed to ninja_filter [required].")
    p.add_argument("-o", "--output",
                   type = str,
                   default = None,
                   metavar = '',
                   help = "Output folder (default /ninja_output)")
    p.add_argument("-t", "--trim",
                   type = int,
                   default = -1,
                   metavar = '',
                   help = "Trim sequences to a specified number of bp, cutting off from their ends (default no trim)")
    p.add_argument("-s", "--similarity",
                   type = float,
                   default = 0.97,
                   metavar = '',
                   help = "Minimum percent similarity - id - between query sequence and reference sequence (default 97%%)")
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
                   type = float,
                   default = 0.0001,
                   metavar = '',
                   help = "For argument X.Y, discards all reads that appear less than X times and all kmers that appear" + \
                          " less than Y times unless kmers in reads that appear more than X times (default 0.0001)")
    p.add_argument("-r", "--reverse_complement",
                   action = 'store_true',
                   help = "Flags sequences for reverse complementing (default no reverse complement)")
    p.add_argument("-q", "--quiet",
                   action = 'store_false',
                   help = "Flags silent output to disable writing to log file in output folder (default not quiet)")

    args = p.parse_args()
    return args

# Checks if args work. Takes args and parser as input
def check_args(args, p):
    # Yells at user if they don't specify an input file
    if args['input'] is None:
        p.print_help()
        sys.exit('\nPlease include an input sequences file in fasta format.')
    
    # Checks if input sequences fasta is correctly formatted. Writes correct one if not
    else:
        fileName = args['input']
        if not check_fasta(open(fileName)):
            print("ERROR: Input fasta formatted incorrectly for QIIME, e.g. sequences or title on multiple lines. Writing " + \
                  "corrected file in current working directory. File automatically moved to output folder after NINJA runs.")
            with open(fileName) as f:
                write_fasta(read_fasta(f), "formatted_input_fasta")
                args['input'] = "formatted_input_fasta.fna"
            global inputFastaCorrected
            inputFastaCorrected = True

    # Sets default output folder if user doesn't specify
    if args['output'] is None:
        args['output'] = "ninja_output"

# Checks if an input fasta file is formatted correctly for QIIME. Returns boolean.
# Looks for data or titles with multiple lines and improper characters in seqs.
# Makes sure first line is a header. Doesn't allow spaces in seqs.
def check_fasta(f):
    lineNumber = 1
    inData = True
    seqChars = ['A', 'T', 'G', 'C', 'N', '\n']
    for line in f:
        if line[0] == ">":
            if inData:
                inData = False
            else:
                print "ERROR: Multiline title in line " + str(lineNumber) + " of input fasta."
                return False
        else:
            # Checks if line is actually a title
            for ch in line:
                c = ch.upper()
                if c not in seqChars:
                    print 'ERROR: Unsupported sequence character "' + c + '" in line ' + str(lineNumber) + ' of input fasta.'
                    return False
            if not inData:
                inData = True
            else:
                print "ERROR: Multiline sequence in line " + str(lineNumber) + " of input fasta."
                return False
        lineNumber += 1
    return True

# Generator for fasta files - returns [(name, seq)]
# Call using 'with open(file) as f'
def read_fasta(f):
    title = None
    data = None
    for line in f:
        if line[0]==">":
            if title != None:
                yield (title,data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if title == None:
        yield None
    yield (title,data)

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

# Prints an error to the console, as opposed to the log. Exits program if exit flagged
def error(e, msg = "", exit = False):
    sys.stdout = console    # Switches output to console
    if e: print(e)
    print(msg)

    sys.stdout = open(ninjaLog, 'w')   # Resets output to log file
    if exit:
        print("CRITICAL ERROR. Check console for more information.")
        sys.exit(0)

###
#   MAIN METHODS
###

# Runs ninja_filter
# INPUT     inputSeqsFile:          input sequences in fasta format
# OUTPUT    filteredSeqsFile:       input sequences filtered using ninja algorithm
#           seqsDB:                 utility file used in ninja_parse
# OPTIONAL  trim:                   trims sequences to <= uniform X bp (e.g. AGGC, GCG with trim 2 returns AG, GC)
#           RC:                     takes reverse complement of input sequences
#           denoising:              for argument x.y, discards all reads that appear less than x times and all kmers that
#                                   appear less than y times unless kmers in reads that appear more than x times
def ninja_filter(inputSeqsFile, filteredSeqsFile, seqsDBFile, trim, RC, denoising):

    # TODO: Add global flag to return commands run from each function. (supress execution)

    # Converts optional arguments to args readable by ninja_filter
    argTrim = ''
    argRC = ''
    argDenoising = ''
    if trim != -1:
        argTrim = str(trim)
    if RC:
        argRC = "RC"

    # Runs ninja_filter. Run in shell only on Mac
    cmd = ""
    try:
        global shellBool
        cmd = ninjaFilterFile + ' ' + inputSeqsFile + ' ' + filteredSeqsFile + ' ' + seqsDBFile + ' ' + argTrim + \
              ' ' + argRC + ' ' + argDenoising
        subprocess.check_call(cmd, shell = shellBool, stdout = sys.stdout)
    except subprocess.CalledProcessError as e:
        error(e, msg = "ERROR: Filtering failed. Check input FASTA formatting and input file locations. Exiting.", exit = True)
    return cmd
    

# Runs bowtie2. Uses two presets for ninja normal and max
# INPUT     filteredSeqsFile:   filtered sequences output from ninja_filter
# OUTPUT    alignmentsFile:     main output of bowtie2
# OPTIONAL  mode:               'ninja' or 'ninjaMax', for less and more sensitivity, respectively
#           threads:            number of threads/cores to run bowtie2 on
#           similarity:         minimum percent similarity between query sequence and reference sequence
def bowtie2(filteredSeqsFile, alignmentsFile, bowtieDatabase, similarity, threads, mode):

    # TODO: Automatically convert fasta file if formatted incorrectly

    # Checks if similarity is a percentage
    if similarity > 1 or similarity < 0:
        print("Similarity error. Enter similarity as a percentage between 0 and 100. Exiting.")
        sys.exit()

    similarity = 1 - similarity     # Converts to similarity readable by bowtie2

    # Switches between ninja normal and max according to user input. Only runs in shell on Mac
    global shellBool
    cmd = ""
    if mode == 'normal':
        #if verbose: print("\nRunning normal bowtie2 to create alignment...\n")
        try:
            cmd = bowtie2File + ' --no-head --no-unal -x ' + bowtieDatabase +  ' -S ' + alignmentsFile + \
                ' --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-' + str(similarity) + \
                '" -k 1 --norc -f ' + filteredSeqsFile + ' -p ' + str(threads)
            subprocess.check_call(cmd, shell = shellBool, stderr = sys.stdout)
        except subprocess.CalledProcessError as e:
            print(e.cmd)
            error(e, msg = "ERROR: Bowtie2 failed. Exiting.", exit = True)
    elif mode == 'fast':
        #if verbose: print("\nRunning fast bowtie2 to create alignment...\n")
        try:
            cmd = bowtie2File + ' --no-head --no-unal -x ' + bowtieDatabase + ' -S ' + alignmentsFile + \
                ' --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-' + str(similarity) + \
                '" -k 1 --norc -f ' + filteredSeqsFile + ' -p ' + str(threads) + ' --fast'
            subprocess.check_call(cmd, shell = shellBool, stderr = sys.stdout)
        except subprocess.CalledProcessError as e:
            print(e.cmd)
            error(e, msg = "ERROR: Bowtie2 failed. Exiting.", exit = True)
    else:
        #if verbose: print("\nRunning very sensitive bowtie2 to create alignment...\n")
        try:
            cmd = bowtie2File + ' --no-head --no-unal -x ' + bowtieDatabase + ' -S ' + alignmentsFile + \
                ' --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-' + str(similarity) + \
                '" --norc -f ' + filteredSeqsFile + ' -p ' + str(threads) + \
                ' --very-sensitive'
            subprocess.check_call(cmd, shell = shellBool, stderr = sys.stdout)
        except subprocess.CalledProcessError as e:
            error(e, msg = "ERROR: Bowtie2 failed. Exiting.", exit = True)
    return cmd

# Runs ninja_parse_filtered.
# INPUT     seqsDBFile:       db file output from ninja_filter
#           alignmentsFile:   alignment file output from bowtie2
#           masterDBFile:     master db file packaged with ninja
#           taxMapFile:       reference taxonomy map packaged with ninja
# OUTPUT    otuTableFile:     OTU table output that's really the point of all this
# AUTO      legacyTable:      OTU table in legacy format, output automatically
#           parseLog:         a utility file containing parsed sequences, used in post-processing and output automatically
def ninja_parse(seqsDBFile, alignmentsFile, masterDBFile, taxMapFile, otuTableFile):
            
    global shellBool
    try:
        cmd = ninjaParseFile + ' ' + seqsDBFile + ' ' + alignmentsFile + ' ' + masterDBFile + ' ' + \
              taxMapFile + ' ' + otuTableFile + ' --legacy' 
        subprocess.check_call(cmd, shell = shellBool, stdout = sys.stdout)
        #if verbose: print("NINJA complete. Begin post-processing.\n")
    except subprocess.CalledProcessError as e:
        error(e, msg = "ERROR: Parsing failed. Exiting.", exit = True)

    # TODO (Gabe): Output legacy table and regular biom file 
    return cmd

# Processes NINJA output. Finds all rejected seqs and outputs them to a FASTA file. Also creates an OTU map using accepted seqs.
# INPUT     inputSeqsFile:      original file of sequences passed to ninja_filter (.fna, .fasta)
#           filteredSeqsFile:   seqs with duplicates removed output from ninja_filter
#           parseLogFile:       parse log generated from ninja_parse_filter
# OUTPUT    seqOutFile:         file to write failed sequences (.fna)
#           mapOutFile:         file to write OTU map (.txt)
# OPTIONAL  trim:               matches trim of input seqs to trim specified in ninja_filter (e.g. 200 bp)
#           RC:                 takes RC of input seqs if specified in ninja_filter
def process(inputSeqsFile, filteredSeqsFile, parseLogFile, seqOutFile, mapOutFile, trim, RC):

    # Opens filtered sequences file and stores in the dict "filteredSeqs" with format {seq:ninjaID}
    try:
        filteredSeqs = {}
        with open(filteredSeqsFile) as f:
            for header, seq in read_fasta(f):
                filteredSeqs[seq] = header
    except IOError as e:
        error(e, msg = "ERROR: Filtered sequences file not found. Exiting.", exit = True)

    # Opens parse log file and stores in the dict "parseLog" with format {ninjaID:OTU}
    try:
        parseLog = {}
        with open(parseLogFile) as f:
            for line in f:
                l = line.strip().split()
                parseLog[l[0]] = l[1]
                
    except IOError as e:
        error(e, msg= "ERROR: Parse log file not found. Exiting.", exit = True)

    # 1)  Opens input seq file and iterates through it
    # 1b) Changes seq to RC or trims it if user specified
    # 2)  Searches for each seq of input in filtered sequences, grabs its index
    # 3)  Searches for index in parse log
    # 4)  If index found, appends seqID to list of OTU:seqIds
    # 5)  If index not found, grabs header and sequence from inputSeqs
    try:
        # Temporary variables
        ninjaID = ""                        # Filtered seqs -> Parse Log
        seqID = ""                          # Parse log -> Reference Map
        otuID = ""                          # In reference map
        failuresOutput = []                 # Array of tuples (header, seq) written to failures file
        mapOutput = {}                      # Dict of {OTU ID: [seqID]} written to map file
        numberSequences = 0;                # Total number of sequences failed
        if verbose: print("Searching for failed sequences and generating OTU map...")
        for header, seq in read_fasta(open(inputSeqsFile)):       # Opens (large) input seq file
            if RC:                                 # Handles RC
                seq = reverse_complement(seq)
            if ((trim > 0) and (trim != -1)):      # Handles trim
                seq = seq[0:trim]
            if seq in filteredSeqs:                # Checks if input seq filtered  
                ninjaID = filteredSeqs[seq]        # Seq -> Ninja ID
                ninjaID = ninjaID.strip()
                if ninjaID in parseLog:
                    otuID = parseLog[ninjaID]      # Ninja ID -> OTU ID
                    if otuID not in mapOutput:           # Checks for existing OTU ID in OTU map
                        mapOutput[otuID] = header.split(' ', 1)[0] + "\t"     # Adds new line with OTU ID and one Seq ID to OTU map
                    else:
                        mapOutput[otuID] += header.split(' ', 1)[0] + "\t"    # Appends seq ID to list of seq ID's for a given OTU
                else:
                    failuresOutput.append((header, seq))
                    numberSequences += 1
  
  
        if verbose: print("OTU Map generated. A total of " + str(numberSequences) + " sequences were recovered.")
    
    except IOError as e:
        error(e, msg = "ERROR: Input seqs not found. Exiting.", exit = True)

    try:
        if verbose: print("Writing files...")
        write_fasta(failuresOutput, seqOutFile)
        write_map(mapOutput, mapOutFile)    # Sorts map in ascending order of OTUs
    except IOError as e:
        error(e, msg = "ERROR: Writing to files failed. Exiting.", exit = True)

# Performs housekeeping on files, deleting the intermediate ones listed below
# INPUT     inputSeqsFile:      original file of sequences passed to ninja_filter (.fna, .fasta)
#           filteredSeqsFile:   seqs with duplicates removed output from ninja_filter
#           seqsDBFile:         db file output from ninja_filter
#           alignmentsFile:     main output of bowtie2
#           parseLogFile:       parse log generated from ninja_parse_filter
def clean(inputSeqsFile, filteredSeqsFile, seqsDBFile, alignmentsFile, parseLogFile):
    if verbose: print("Cleaning up output...")
    try:
        global out
        os.remove(filteredSeqsFile)
        os.remove(seqsDBFile)
        os.remove(alignmentsFile)
        os.remove(parseLogFile)
        os.remove("map_seqid_reps.txt")

        # Moves corrected input sequences file to output if it was written
        global inputFastaCorrected
        if inputFastaCorrected:
            global out
            newInputSeqsFile = os.path.join(out, inputSeqsFile).replace("\\", "/") 
            os.rename(inputSeqsFile, newInputSeqsFile)

        # Deletes ninja log if quiet enabled
        if not verbose: os.remove(ninjaLog)
    except OSError as e:
        error(e, msg = "INTERNAL ERROR: Can't find all files marked for moving and/or deletion. Check working directory and output folder.")
    if verbose: print("Done.")

# Runs ninja, bowtie2 and then processes output. All files output in specified output folder. 
# User must specify ninja's directory as an environment variable named 'NINJA_DIR'
def main(inputSeqsFile, folder, trim, RC, similarity, threads, mode, denoising, verboseBool):

    # Gets ninja's directory relative to current working directory
    ninjaDirectory = os.path.relpath(os.path.dirname(os.path.realpath(__file__)), os.getcwd()).replace("/", "\\") 

    # Checks for output subdirectory of current working directory. Makes it if necessary. 
    # Edits global output folder variable
    global out
    out = os.path.join(os.getcwd(), folder)
    if not os.path.exists(out):
        os.makedirs(out)
    out = os.path.abspath(out)

    # Opens log file in output folder for printing and sets it as a global variable
    # First stores original console location as a variable for error handling
    global console
    global ninjaLog
    console = sys.stdout
    ninjaLog = os.path.join(out, "ninja_log.txt").replace("\\", "/") 
    sys.stdout = open(ninjaLog, 'w')
        
    global shellBool
    osName = sys.platform
    shellBool= not (osName.startswith("win32") or osName.startswith("cygwin"))

    # Checks for bowtie2 file if not in ninja_package.
    global bowtie2File
    bowtie2File = "bowtie2-align-s"
    # find out if bowtie2 is in the system path
    try:
        subprocess.check_call(bowtie2File + " --version", shell = shellBool, stdout = sys.stdout)
    except OSError as e:
        bowtie2File = os.path.abspath(os.path.join(ninjaDirectory, "bowtie2-align-s")).replace("\\", "/")
        if not os.path.exists(bowtie2File):
            error(e = None, msg = "ERROR: Bowtie2 executable not found in system path or top-level NINJA package folder. Please install bowtie2 and add its accompanying executables to the system path or place bowtie2-align-s in the top-level ninja package folder (not a subfolder). " + \
                              "Check README.txt for additional instructions. Exiting.", exit = True)

    # Sets the relevant globals and binaries for mac and windows support. Linux upcoming
    global ninjaFilterFile
    global ninjaParseFile
    global verbose
    verbose = verboseBool
    if osName.startswith("darwin") or osName.startswith("os"):			# Mac
    	ninjaFilterFile = os.path.join(ninjaDirectory, os.path.join("bin", "./ninja_filter_mac")).replace("\\", "/") 
    	ninjaParseFile = os.path.join(ninjaDirectory, os.path.join("bin", "./ninja_parse_filtered_mac")).replace("\\", "/") 
    elif osName.startswith("win32") or osName.startswith("cygwin"):		# Windows and cygwin
    	ninjaFilterFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_filter.exe")).replace("\\", "/") 
    	ninjaParseFile = os.path.join(ninjaDirectory, os.path.join("bin", "ninja_parse_filtered.exe")).replace("\\", "/")
        bowtie2File = bowtie2File + ".exe"
    else:   # Linux
        ninjaFilterFile = os.path.join(ninjaDirectory, os.path.join("bin", "./ninja_filter_linux")).replace("\\", "/") 
        ninjaParseFile = os.path.join(ninjaDirectory, os.path.join("bin", "./ninja_parse_filtered_linux")).replace("\\", "/")


        # Sets variables used in ninja calls. First, ninja_filter files
    filteredSeqsFile = os.path.join(out, "filtered_sequences.fa").replace("\\", "/") 
    seqsDBFile = os.path.join(out, "seqsDB.db").replace("\\", "/") 
        # Bowtie2 files
    alignmentsFile = os.path.join(out, "alignments.txt").replace("\\", "/") 
    masterDBFile = os.path.abspath(os.path.join(ninjaDirectory, os.path.join("bt2db", "masterDB97.map"))).replace("\\", "/") 
    bowtieDatabase = os.path.abspath(os.path.join(ninjaDirectory, os.path.join("bt2db", "Ninja97"))).replace("\\", "/") 
        # Ninja_parse files
    taxMapFile = os.path.abspath(os.path.join(ninjaDirectory, os.path.join("bt2db", "taxMap97.txt"))).replace("\\", "/") 
    parseLogFile = os.path.join(os.getcwd(), "parseLog.txt").replace("\\", "/") 
    otuTableFile = os.path.join(out, "otutable.biom").replace("\\", "/") 
        # Post-processing files
    seqOutFile = os.path.join(out, "failed_sequences.fna").replace("\\", "/") 
    mapOutFile = os.path.join(out, "otu_map.txt").replace("\\", "/") 

    # Runs ninja_filter, bowtie2 and ninja_parse. Processes ninja results, generating OTU map and a list of failed seqs
    print
    t1 = timeit.Timer(lambda: ninja_filter(inputSeqsFile, filteredSeqsFile, seqsDBFile, trim, RC, denoising))
    if verbose: 
        print("Ninja filter time: " + str(t1.timeit(1)))
    else:
        t1.timeit(1)
    t2 = timeit.Timer(lambda: bowtie2(filteredSeqsFile, alignmentsFile, bowtieDatabase, similarity, threads, mode))
    if verbose: 
        print("Bowtie time: " + str(t2.timeit(1)))
    else:
        t2.timeit(1)
    t3 = timeit.Timer(lambda: ninja_parse(seqsDBFile, alignmentsFile, masterDBFile, taxMapFile, otuTableFile))
    if verbose: 
        print("Ninja parse time: " + str(t3.timeit(1)) + "\n")
    else:
        t3.timeit(1)
    t4 = timeit.Timer(lambda: process(inputSeqsFile, filteredSeqsFile, parseLogFile, seqOutFile, mapOutFile, trim, RC))
    if verbose: 
        print("Post-processing time: " + str(t4.timeit(1)))
    else:
        t4.timeit(1)
    t5 = timeit.Timer(lambda: clean(inputSeqsFile, filteredSeqsFile, seqsDBFile, alignmentsFile, parseLogFile))
    if verbose: 
        print("Clean-up time: " + str(t5.timeit(1)))
    else:
        t5.timeit(1)

# Wrapper for main function, called from command line
# Bare minimum args:
#   -i "seqs.fna" 
# Sample maximum args:
#   -i "seqs.fna" -o "output" -r -t 200 -mo 'max' -s 98 -d 1.005 -q
if __name__=='__main__':
    # Parses command line arguments
    p = argparse.ArgumentParser(description = "NINJA OTU Picker: NINJA Is Not Just Another OTU Picker -- filter program\n" + \
                                              "Knights Lab (www.knightslab.org/ninja)\n" + \
                                              "This program outputs an otu table and map from sequence reads in fasta format.", 
                                add_help = True, 
                                epilog ='NOTE: If one or more output files are empty, trying reverse complementing your input ' + \
                                        'sequences with -r')
    args = get_args(p)
    args = vars(args)
    check_args(args, p)

    # Runs ninja pipeline
    t = timeit.Timer(lambda: main(args['input'], args['output'], args['trim'], args['reverse_complement'], args['similarity'], 
                                  args['threads'], args['mode'], args['denoising'], args['quiet']))
    if args['quiet']: 
        print("Total time: " + str(t.timeit(1)))
    else: 
        t.timeit(1)

