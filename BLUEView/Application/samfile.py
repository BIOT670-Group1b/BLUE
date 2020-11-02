import pysam

from pyensembl import Genome
from pyensembl import fasta

def load_chromosome(file,type):
    global ChromData
    ChromData = fasta.parse_fasta_dictionary(file)[type]

# Loads BAM file and stores it in a global variable for class to use
def load_sam(file):
    global results

    sam = pysam.AlignmentFile(file, "rb")
    results = sam.fetch()

# Loads reference genome and stores it in a glabal variable
def load_genome(file):
    global data
    data = Genome(reference_name='GRCh38', annotation_name='ENSEMBL', gtf_path_or_url=file)
    data.index()

# Stores every read in list to be referenced by index
def store_all_reads():
    init_all_reads()

    for read in results:
        # Format contig to be numeric
        chromosome = read.reference_name
        chromosome = chromosome.replace('chr', '')
        chromosome = int(chromosome)

        all_queries.append(read.query_name)
        all_positions.append(read.pos + 1) # pysam counts the start position one behind the given position in the bam file
        all_lengths.append(read.query_length)
        all_chromosomes.append(chromosome)
        all_sequences.append(read.seq)
        all_cigars.append(read.cigar)

# Initializes global read lists
def init_all_reads():
    global all_queries
    global all_positions
    global all_lengths
    global all_chromosomes
    global all_sequences
    global all_cigars
    all_queries = []
    all_positions = []
    all_lengths = []
    all_chromosomes = []
    all_sequences = []
    all_cigars = []


# Get the unique reads
def get_unique_positions():
    unique_positions = []
    for pos in all_positions:
        if pos not in unique_positions:
            unique_positions.append(pos)

    return unique_positions

# Get the gene names associates with each read
def get_gene_names(Length, Position, Chromosome):
    start_position = Position
    end_position = start_position + Length
    chromosome = Chromosome

    genes = []

    gene_names = data.genes_at_locus(contig=chromosome, position=start_position, end=end_position)
    for gene in gene_names:
        if gene.biotype == 'protein_coding':
            genes.append(gene.gene_name)

    return genes

# Get the ID's of genes based on size, position and chromosome
def get_gene_ids(Length, Position, Chromosome):
    start_position = Position
    end_position = start_position + Length
    chromosome = Chromosome

    genes = []

    gene_names = data.genes_at_locus(contig=chromosome, position=start_position, end=end_position)
    for gene in gene_names:
        if gene.biotype == 'protein_coding':
            genes.append(gene.gene_id)

    return genes

# Gets all genes associated with a BAM file
def get_all_gene_names():
    read_gene_names = []

    for i in range(len(all_positions)):
        read_gene_names[i] = get_gene_names(all_lengths[i], all_positions[i], all_chromosomes[i])
    return read_gene_names

# Given a position and chromosome, return info about all genes located
def get_gene_info(start_position, chromosome):
    gene_names = data.genes_at_locus(contig=chromosome, position=start_position)
    genes = []
    for gene in gene_names:
        genes.append(gene)
    return genes

# Given a gene name, return it's ID
def get_id_by_name(gene_name):
    gene_id = data.gene_ids_of_gene_name(gene_name)
    return gene_id[0]

# Returns a string position range of a gene given a name (with size)
def get_positions_by_name(gene_name):
    positions = data.loci_of_gene_names(gene_name)
    start = positions[0].start
    end = positions[0].end
    return str(start), str(end)

# Returns the starting position of a gene by name
def get_start_position_by_name(gene_name):
    positions = data.loci_of_gene_names(gene_name)
    start_position = positions[0].start
    return start_position

# Returns the chromosome of a given gene name
def get_chromosome_by_name(gene_name):
    gene = data.genes_by_name(gene_name)
    chromsome = gene[0].contig
    return chromsome

# Get the transcript ids using the gene name
def get_transcript_ids_by_name(gene_name):
    transcript_ids = []
    transcripts = data.transcript_ids_of_gene_name(gene_name)
    for t in transcripts:
        transcript_ids.append(t)
    transcript_ids = sorted(transcript_ids)
    return transcript_ids

# Get the name of the transcript using the transcript id
def get_transcript_name_by_id(t_id):
    transcript = []
    t_name = data.transcript_name_of_transcript_id(t_id)
    transcript.append(t_name)
    return transcript[0]

# Get the start and end position of a transcript given an id
def get_transcript_positions_by_id(t_id):
    positions = data.locus_of_transcript_id(t_id)
    start = positions.start
    end = positions.end
    return str(start), str(end)

# Get all exon regions associated with a transcript id
def get_exon_regions_by_transcript_id(t_id):
    exon_starts = []
    exon_ends = []
    exon_ids = data.exon_ids_of_transcript_id(t_id)
    for exon_id in exon_ids:
        start = data.locus_of_exon_id(exon_id).start
        end = data.locus_of_exon_id(exon_id).end
        if start not in exon_starts:
            exon_starts.append(start)
            exon_ends.append(end)

    exon_starts = sorted(exon_starts)
    exon_end = sorted(exon_ends)
    return exon_starts, exon_ends

# Returns exon regions given a gene name
def get_exon_regions_by_name(gene_name):
    exon_regions = []
    exon_ids = data.exon_ids_of_gene_name(gene_name)
    for exon_id in exon_ids:
        start = data.locus_of_exon_id(exon_id).start
        end = data.locus_of_exon_id(exon_id).end
        exon_string = str(start) + " - " + str(end)
        if exon_string not in exon_regions:
            exon_regions.append(exon_string)

    exon_regions = sorted(exon_regions)
    return exon_regions


# Returns chromosome sequence
def get_chromosome_seq():
    return ChromData

# Return the position of a sequence in the loaded chromosome sequence
def get_chrom_seq_pos(seq):
    return ChromData.find(seq)


# Using the data stored, locate any fusion reads within the bam file and show it on screen
def map_fusion_reads():
    import GUIDesign

    global fusion_positions
    global fusion_lengths

    fusion_positions = []
    fusion_lengths = []

    positions = get_unique_positions()
    fusion_position = 0
    supp_position = 0
    position_headers = [] # Initialize positions for the reads
    paired_reads = []

    i = len(positions)
    total_fusions = 0 # Stores the number of fusion reads found
    isFirst = False # Used to create headers if fusion was found
    isSet = False # Used to capture the first index for comparisons against all subsequent ones
    start = 0 # Index incrememented after each comparison
    x = 0

    if i > 1:
        while x < len(all_positions):
            right_indx = x

            if not isSet:
                left_indx = x
                isSet = True

            x = x + 1
            # If we reach a new instance of positions, we want to test them to see if the reads fuse
            if all_positions[left_indx] != all_positions[right_indx] and left_indx not in paired_reads and right_indx not in paired_reads:
                fusion, fusion_indx, supplement, supp_indx = get_fusion_read(left_indx, right_indx)
                paired_reads.append(left_indx)
                paired_reads.append(right_indx)

                if fusion != "No Fusions Detected":
                    total_fusions = total_fusions + 1
                    # Store the fusion reads for larger scope reference
                    fusion_positions.append(all_positions[fusion_indx])
                    fusion_lengths.append(all_positions[fusion_indx])

                    fusion_length = all_lengths[fusion_indx]
                    fusion_position = all_positions[fusion_indx]
                    fusion_chr = all_chromosomes[fusion_indx]

                    supp_length = all_lengths[supp_indx]
                    supp_position = all_positions[supp_indx]
                    supp_chr = all_chromosomes[supp_indx]
                    # Get the difference of lengths to focus the fusion area
                    len_diff = fusion_length - supp_length
                    if total_fusions == 1:
                        isFirst = True
                        position_headers = get_position_headers(fusion_length, fusion_position, 50)

                    # If no reference genome has been loaded, we only can display the sequences, not names
                    try:
                        data
                    except NameError:
                        GUIDesign.insert_read(fusion, supplement, "", "", "", "", isFirst, position_headers,
                                              fusion_position, supp_position, fusion_length, supp_length)
                    else:
                        fusion_name = get_gene_names(fusion_length, fusion_position, fusion_chr)
                        fusion_id = get_gene_ids(fusion_length, fusion_position, fusion_chr)

                        supp_name = get_gene_names(supp_length, supp_position, supp_chr)
                        supp_id = get_gene_ids(supp_length, supp_position, supp_chr)

                        GUIDesign.insert_read(fusion, supplement, fusion_name, fusion_id[0], supp_name, supp_id[0],
                                              isFirst, position_headers, fusion_position, supp_position, fusion_length, supp_length)

                start = start + 1
                x = start
                isSet = False
                isFirst = False
    if total_fusions == 0:
        GUIDesign.insert_read(fusion, supplement, "", "", "", "", isFirst,
                              position_headers, fusion_position, supp_position,"","")


# Checks if two reads are a fusion based on unique indeces
def get_fusion_read(left_indx, right_indx):
    fusion = "No Fusions Detected"
    fusion_indx = -1
    supplement = "No Fusions Detected"
    supp_indx = -1
    fusion_pos = 0

    Read1Length = all_lengths[left_indx]
    Read1Seq = all_sequences[left_indx]

    Read2Length = all_lengths[right_indx]
    Read2Seq = all_sequences[right_indx]

    if Read1Length > Read2Length:
        if Read1Seq.endswith(Read2Seq) or Read1Seq.startswith(Read2Seq):
            fusion = Read1Seq
            fusion_indx = left_indx
            supplement = Read2Seq
            supp_indx = right_indx
    elif Read2Length > Read1Length:
        if Read2Seq.endswith(Read1Seq) or Read2Seq.startswith(Read1Seq):
            fusion = Read2Seq
            fusion_indx = right_indx
            supplement = Read1Seq
            supp_indx = left_indx

    return fusion, fusion_indx, supplement, supp_indx

# This takes a number and returns an evenly distributed list of positions to use as headers
def get_position_headers(length, start_position, pos_between):
    positions = []

    for x in range(length):
        if x == 0:
            positions.append(start_position)
        if x % pos_between == 0 and x != 0: # This adjustable number states the number of base pairs before a position header
            position = start_position + x
            positions.append(position)
        elif x == length-1:
            position = start_position + length
            positions.append(position)
    return positions

# Use this to reference the found fusions
def get_fusions():
    return fusion_positions, fusion_lengths

# Stores the cigar in two lists to be referenced by index
# 0 = match; 1 = insertion; 2 = deletion; 3 = skip
# 4 = soft clipping; 5 = hard clipping, 6 = padding
def parse_cigar(cigar):
    cigar_type = []
    cigar_length = []
    for (cigarType, cigarLength) in cigar:
            cigar_type.append(cigarType)
            cigar_length.append(cigar_length)
    return cigar_type, cigar_length