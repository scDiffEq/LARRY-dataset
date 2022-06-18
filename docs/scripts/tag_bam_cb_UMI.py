#!/usr/bin/python

import pysam, time


def _extract_cell_UMI_barcodes(read):

    """
    For a single read, fetch the cell barcode and the read barcode.

    Parameters:
    -----------
    read
        pysam <class 'pysam.libcalignedsegment.AlignedSegment'> generated using samfile.fetch

    Returns:
    --------
    cell_barcode, UMI_barcode

    Note:
    -----
    This protocol is designed for use with the LARRY dataset. Datasets have varying read structure which
    may prevent the further application of this barcode extraction function.
    """

    read_name = read.qname.split(":")

    cell_barcode = read_name[0]
    UMI_barcode = read_name[1]

    return cell_barcode, UMI_barcode


def _add_cellbarcode_UMI_to_bam_tags(inbam_path, out_path="./barcode_tagged.bam"):

    """
    For a bam with cell barcode and UMI in the read header, fetch the cell barcode and the read barcode and
    add these as tags to the .bam file.

    Parameters:
    -----------
    inbam_path
        path to bam file

    Returns:
    --------
    barcoded bam to out_path

    Note:
    -----
    This protocol is designed for use with the LARRY dataset. Datasets have varying read structure which
    may prevent the further application of this barcode extraction function.
    """

    untagged_bam = pysam.AlignmentFile(inbam_path, "rb")
    tagged_bam = pysam.AlignmentFile(out_path, "wb", template=untagged_bam)

    barcoded_read_count = 0
    start_time = time.time()
    for i, read in enumerate(untagged_bam.fetch(until_eof=True)):

        cell_bc, UMI_bc = _extract_cell_UMI_barcodes(read)

        tags = read.get_tags()
        tags.append(("CB", cell_bc, "Z"))
        tags.append(("UB", UMI_bc, "Z"))
        read.set_tags(tags)
        tagged_bam.write(read)

        

        if i % 10000000 == 0:
            if i != 0:
                current_time = time.time()
                elapsed_time = (current_time - start_time) / 60
                percent_barcoded = (barcoded_read_count / i) * 100
                print("Number of reads processed:", i)
                print("{:.2f}".format(elapsed_time), "minutes elapsed.")
                print(
                    "Number of reads barcoded:",
                    barcoded_read_count,
                    "Percentage barcoded:",
                    "{:.2f}".format(percent_barcoded),
                    "%",
                )
        barcoded_read_count += 1

    untagged_bam.close()
    tagged_bam.close()


LARRY_bam = "/home/mvinyard/data/raw/weinreb_2020/bams/LARRY.sorted.bam"
out_path = "/home/mvinyard/data/raw/weinreb_2020/bams/LARRY.sorted.tagged.bam"

print("Adding CB and UB tags to BAM file...")

_add_cellbarcode_UMI_to_bam_tags(LARRY_bam, out_path)

