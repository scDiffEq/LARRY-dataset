############################################################
# create log
############################################################

def create_logger(logname="debug.log"):
    
    import os, logging
    
    logdir="/home/mvinyard/scripts/.process_logs"
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    
    logname = os.path.join(logdir, logname)
    
    if os.path.exists(logname):
        os.remove(logname)

    # now we will Create and configure logger
    # set the threshold of logger to INFO 
    logging.basicConfig(
        level="INFO",
        format="%(asctime)s - %(name)s - [ %(message)s ]",
        datefmt="%d-%b-%y | %H:%M:%S",
        force=True,
        handlers=[logging.FileHandler(logname), logging.StreamHandler()],
    )

    # Let us Create an object
    logger = logging.getLogger()

    return logger

logger = create_logger(logname="filter.fastqs.Weinreb.2020.debug.log")

logger.info("logger created.")
                       
############################################################
# import libraries
############################################################
logger.info("loading libraries.")

import scdiffeq as sdq

odeint, torch, np, pd, plt, nn, a, os, time, optim, sp, PCA, v = sdq.ut.devlibs()

logger.info("loading libraries complete.")
       
############################################################
# define functions
############################################################

def get_files_in_path(path):
    
    import glob

    if path.endswith("/"):
        mod_path = "".join([path, "*"])
    else:
        mod_path = "/".join([path, "*"])
        path = "".join([path, "/"])

    files_in_dir = glob.glob(mod_path)

    seq_libraries = []
    for file in files_in_dir:
        library = file.split(".")[-5].split("/")[-1]
        seq_libraries.append(library)        

    return seq_libraries, files_in_dir

def get_sequencing_library_final_barcodes(adata, library):

    unique_bcs = adata.obs.loc[adata.obs.Library == library]["Cell barcode"].unique()
    print(str(int(unique_bcs.shape[0])), "unique barcodes assocaited with:", library)

    return unique_bcs


def count_reads_in_file(filepath, lines_per_read=4, ret=False):

    print("Checking", filepath, "file size...")
    command_string = " ".join(["zcat", fq_path, "| wc -l"])
    num_lines = int(os.popen(command_string).read().split()[0])
    num_reads = int(num_lines / lines_per_read)

    if ret == True:
        return num_lines


def initiate_outFASTQ(sequencing_library, outdir="./"):

    output_path = "%s.LARRY.%s.fastq" % (sequencing_library, "filtered_on_final")
    full_outpath = os.path.join(outdir, output_path)
    outfile = open(full_outpath, "w")

    return outfile, full_outpath


def write_fastq_read(outfile, metadata, read, phredd):

    """"""

    outfile.writelines([metadata, "\n"])
    outfile.writelines([read, "\n"])
    outfile.writelines(["+", "\n"])
    outfile.writelines([phredd, "\n"])


def record_barcoded_read(line, input_FASTQ, filtered_FASTQ, line_counter):

    line_counter += 1
    read = input_FASTQ.readline().decode("utf-8").strip("\n")
    line_counter += 1
    plus = input_FASTQ.readline().decode("utf-8").strip("\n")
    line_counter += 1
    phredd = input_FASTQ.readline().decode("utf-8").strip("\n")
    # write
    write_fastq_read(outfile=filtered_FASTQ, metadata=line, read=read, phredd=phredd)

    return line_counter


def filter_FASTQ_on_barcodes(path, sequencing_library, unique_bcs, outpath="./"):

    """"""
    import gzip
    
    input_FASTQ = gzip.open(path)
    #     num_lines = count_reads_in_file(path, lines_per_read=4, ret=True)
    line = input_FASTQ.readline().decode("utf-8").strip("\n")

    filtered_FASTQ, filtered_outpath = initiate_outFASTQ(sequencing_library, outpath)
    line_counter, reads_in_barcodes, reads_not_in_barcodes = 0, 0, 0

    print(
        "Processing and filtering FASTQ reads based on final library cell barcodes...",
        "\n",
    )

    while not (line == ""):

        if line[0] == "@":
            cell_bc = line[1:].split(":")[0]

            if cell_bc in unique_bcs:
                reads_in_barcodes += 1
                line_counter = record_barcoded_read(
                    line, input_FASTQ, filtered_FASTQ, line_counter
                )
            else:
                reads_not_in_barcodes += 1

        #### update read processing status ####
        line_counter += 1
        if line_counter % (4e6) == 0:
            print(str(int(line_counter / 4e6)), "million reads processed...")
        line = input_FASTQ.readline().decode("utf-8").strip("\n")
    filtered_FASTQ.close()
                       
    return filtered_outpath
############################################################
# load data
############################################################   
logger.info("load data.")
                       
FASTQ_path = "/home/mvinyard/data/raw/weinreb_2020/fastqs/raw_downloaded/"
adata_path = "/home/mvinyard/data/preprocessed/weinreb_2020/author_preprocessed_download/LARRY.h5ad"

filterd_FASTQ_path = "/home/mvinyard/data/raw/weinreb_2020/fastqs/filtered/"

adata = a.read_h5ad(adata_path)
                       
logger.info(adata)
logger.info("load data complete.")
############################################################
# adjust names of libraries to ensure matching
############################################################

# get all libraries according to FASTQs on disk
lib, lib_paths = get_files_in_path(FASTQ_path)

print("The following library names will be changed in AnnData:")

# loop through adata.obs library values
for i in adata.obs.Library.unique():
    if i in lib:
        continue
    else:
        # loop through gsutil-stored library values
        for libval in lib:
            # I don't want to repalce the LSK cell libraries
            if i in libval and not ("LSK" in libval):
                print(libval, i)
                adata.obs.Library.replace(i, libval, inplace=True)
print("\n")
print("All AnnData Library values:")
# check library values stored in AnnData to make sure they match
for i in adata.obs.Library.unique():
    print(i)
    
############################################################
# do filtering
############################################################
    
for _lib, _lib_path in zip(lib, lib_paths):
    logger.info(" ".join(["Now processing library", _lib]))
    logger.info(" ".join(["Unfiltered library data path:", _lib_path]))
    print(_lib)
    print(_lib_path)
    library_bcs = get_sequencing_library_final_barcodes(adata, library=_lib)
    filtered_outpath = filter_FASTQ_on_barcodes(
        path=_lib_path, sequencing_library=_lib, unique_bcs=library_bcs, outpath=filterd_FASTQ_path
    )
    logger.info(" ".join(["Filtering complete for library", _lib, "Filtered data path:", filtered_outpath]))
