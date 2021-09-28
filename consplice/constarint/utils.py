from __future__ import print_function
import sys
import io
import os
import gzip
import yaml
from collections import defaultdict





def load_config(path):
    """
    load_config
    ===========
    This method is used to load the internal (or user supplied) config yaml file used for data processing.

    Parameters:
    -----------
    1) path: (str) The file path to the config file to load

    Returns:
    ++++++++
    1) The loaded config file as a dict.
    """
    
    print("\nLoading config file")
    assert os.path.exists(path), "\n!!ERROR!! The config.yml path does not exists. Please fix the problem and try again. Bad file path = '{}'".format(path)
    assert os.path.isfile(path), "\n!!ERROR!! The config.yml path does not exists. Please fix the problem and try again. Bad file path = '{}'".format(path)

    yaml_dict = dict()
    with open(path) as fh:
        yaml_dict = yaml.safe_load(fh)

    if check_config(yaml_dict):
        return(yaml_dict)

    else:
        return({})


def check_config(config_dict):
    """
    check_config
    ============
    This method is used to check that the config file is formatted correctly. 

    Parameters:
    -----------
    1) config_dict: (dict) The config yaml file loaded as a dictionary 

    Returns:
    ++++++++
    1) True (bool) If formated correctly, otherwise it raises an assertion error
    """
    
    print("\nChecking config file")

    ## GENOME BUILD section
    assert "GENOME_BUILD" in config_dict, "!!ERROR!! the 'GENOME_BUILD' key is missing from the config file."

    assert config_dict["GENOME_BUILD"] in ["GRCh37","GRCh38"], "!!ERROR!! only GRCh37 or GRCh38 are allowed for the 'GENOME_BUILD' section of the config file."


    ## Score bin section
    assert "SCORE_BINS" in config_dict, "!!ERROR!! the 'SCORE_BINS' key is missing from the config file"
    
    assert "sum_spliceai_score_bins" in config_dict["SCORE_BINS"], "!!ERROR!! the 'sum_spliceai_score_bins' key is missing from the SCORE_BINS section of the config file"
    assert "max_spliceai_score_bins" in config_dict["SCORE_BINS"], "!!ERROR!! the 'sum_spliceai_score_bins' key is missing from the SCORE_BINS section of the config file"

    assert isinstance(config_dict["SCORE_BINS"]["sum_spliceai_score_bins"], list), "!!ERROR!! the 'sum_spliceai_score_bins' in the 'SCORE_BINS' section should be formatted as a list. Please fix the error"
    assert isinstance(config_dict["SCORE_BINS"]["max_spliceai_score_bins"], list), "!!ERROR!! the 'max_spliceai_score_bins' in the 'SCORE_BINS' section should be formatted as a list. Please fix the error"

    ## PAR Section 
    assert "PAR" in config_dict, "!!ERROR!! the 'PAR' key is missing from the config file"

    for gb in ["GRCh37","GRCh38"]:
        assert gb in config_dict["PAR"], "!!ERROR!! The '{}' key is missing from the 'PAR' section of the config file".format(gb)

        for par_key in ["NOTE","X_PAR1","X_PAR2","X_NONPAR"]:
            assert par_key in config_dict["PAR"][gb], "!!ERROR!! The '{}' key is missing from the '{}' section of the 'PAR' section in the config file".format(par_key, gb)
    
            if par_key in ["X_PAR1", "X_PAR2", "X_NONPAR"]: 
                assert isinstance(config_dict["PAR"][gb][par_key]["start"], int), "!!ERROR!! The start position in the 'PAR' section needs to be an int"
                assert isinstance(config_dict["PAR"][gb][par_key]["end"], int), "!!ERROR!! The end position in the 'PAR' section needs to be an int"


    return(True)


def in_par(par_dict, pos):
    """
    in_par
    =====
    This method is used to check that a X chromosome position is in the pseudoautosomal region (PAR) or not. 
     It uses the PAR region information for the config file to identify if the position is in the PAR or not. 
     This method expects the 'PAR' dict from the config file. (Example: config_dict['PAR']['GRCh38'])

    Parameters:
    -----------
    1) par_dict: (dict) A genome build specific dictionary from the PAR section of the config file. 
    2) pos:      (int)  The X chromosome position to evaluate

    Returns:
    ++++++++
    1) True (bool) if the position is in the PAR, False if in the NonPAR region, ValueError otherwise
    """
    
    if pos > par_dict["X_PAR1"]["start"] and pos <= par_dict["X_PAR1"]["end"]:
        return(True)
    elif pos > par_dict["X_PAR2"]["start"] and pos <= par_dict["X_PAR2"]["end"]:
        return(True)
    elif pos > par_dict["X_NONPAR"]["start"] and pos <= par_dict["X_NONPAR"]["end"]:
        return(False)
    else:
        raise ValueError("Within gene position outside of X chromosome")


def correct_allele_by_strand(strand,allele):
    """
    correct_allele_by_strand
    ========================
    This method is used to get a corrected allele based on the strand a variant is on. 
     If the variant is on the positive strand the allele will stay the same. If the variant
     is on the negative strand the variant will be corrected. (A -> T, T => A, C -> G, G -> C)

    Paremters:
    ----------
    1) strand: (str) The strand the variant is on
    2) allele: (str) The allele to correct

    Returns:
    ++++++++
    1) (str) The corrected allele
    """
    
    corrected_allele = allele
    if strand == "-":
        
        if allele.upper() == "A":
            corrected_allele = "T"
        elif allele.upper() == "T":
            corrected_allele = "A"
        elif allele.upper() == "G":
            corrected_allele = "C"
        elif allele.upper() == "C":
            corrected_allele = "G"

    return(corrected_allele)


def get_alternative_gene_symbols(gene_info_file):
    """
    get_alternative_gene_symbols
    ============================
    This method is used to parse a file with accepted and alternative gene symbols and create a mapping between the pair. 
     This function is set up to work with the HGNC Database mapping file. http://www.genenames.org/. The HGNC has put together
     A large data file that contains mapping between genes across the different data providers and each of the unique ids that 
     each provider makes. Here, the function will use the accepted symbols and the altnerative symbols from this file 
     to create a mapping between the symbols. File url: ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

    Parameters:
    ----------
    1) gene_info_file: (str) The file path to the HGNC symbol mapping file. NOTE: The file needs a "#" at the begining of the 1st line in order to work. 

    Returns:
    +++++++
    1) (dictionary) A dict where keys represent a single gene symbol, and values represent a set of synonymous symbols. The keys will include the accepted and 
                    alternative gene symbols.

    NOTE: Any gene with no alternative symbols will be included as a key with the value as an empty set
    """

    ## Open the alt gene file 
    try:
        gene_info_fh = gzip.open(gene_info_file, "rt", encoding = "utf-8") if gene_info_file.endswith(".gz") else io.open(gene_info_file, "rt", encoding = "utf-8")
    except IOError as e:
        print("!!ERROR!! Unable to open the Alt Gene Symbols file: {}".format(gene_info_file))
        print(str(e))
        sys.exit(1)

    ## Main dictionary for alternative symbols
    alternative_gene_symbols = defaultdict(set)

    ## Get an index from the header 
    header_line = gene_info_fh.readline()
    assert "#" in header_line, "The alternative gene sybmol file does not have a header staring with '#'. A header is required. Please add a header and try again"
    header_index = header_line.strip().replace("#","").split("\t")

    ## Iterate over the file
    for line in gene_info_fh:
         
        ## Get a dictionary for the current line
        line_dict = dict(zip(header_index, line.strip().split("\t")))

        ## Get the offical gene symbol and the synonymous symbols
        gene_symbol = line_dict["symbol"]

        ## Get that alternative gene symbols
        synonymous_symbols = set(x for x in line_dict["alias_symbol"].strip().replace("\"","").split("|")) if line_dict["alias_symbol"] else set()
        ## Add any previous symbols 
        synonymous_symbols.update(set(x for x in line_dict["prev_symbol"].strip().replace("\"","").split("|")) if line_dict["prev_symbol"] else set())

        ## Add to dict 
        alternative_gene_symbols[gene_symbol].update(synonymous_symbols)

        ## iterate over each synonymous symbol
        for synonymous_symbol in synonymous_symbols:

            ## add synonymous symbol 
            alternative_gene_symbols[synonymous_symbol].update(set([x for x in synonymous_symbols if x != synonymous_symbol] + [gene_symbol]))

    ## Close fh
    gene_info_fh.close()

    return(alternative_gene_symbols)


def expand_gene_dict_with_alt_symbols(alt_symbol_dict, gene_strand_dict):
    """
    expand_gene_dict_with_alt_symbols
    =================================
    Method to update the alt symbol dict from the get_alternative_gene_symbols method with strand and gene_ids from 
     a dict with gene id and strand info by gene symbol

    Parmaeters:
    -----------
    1) alt_symbol_dict:  (dict) The alternative gene symbol dict from the get_alternative_gene_symbols method.
    2) gene_strand_dict: (dict) A dict that contains a gene id and strand per key, where each key is a gene symbol

    Returns:
    ++++++++
    1) (dict) An updated alt_symbol_dict with gene id and strand where available.
    """
    
    ## Iterate over all genes in the gene strand dict
    for key in list(gene_strand_dict.keys()) :
        
        ## check if the gene exists in the alt symobl dict
        if key in alt_symbol_dict:
            
            ## Check all alt symobls 
            for gene_name in alt_symbol_dict[key]:
                
                ## If alt symobl is not in the gene strand dict, add it
                if gene_name not in gene_strand_dict:
                    
                    gene_strand_dict[gene_name] = {"gene_name":gene_name, 
                                                   "gene_id":gene_strand_dict[key]["gene_id"], 
                                                   "strand": gene_strand_dict[key]["strand"]}

    return(gene_strand_dict)


def get_delta_score_bin(delta_score, d_score_bins):
    """
    get_delta_score_bin
    ===================
    Method to get the delta score bin a certain score resides in. 

    Parameteres:
    ------------
    1) delta_score:  (float) The delta score to get the bin for
    2) d_score_bins: (list) A list of delta score bin ranges. (example: ["0.0-0.4","0.4-1.0"])
                             used to determine which bin the delta score is in

    Returns:
    ++++++++
    1) (str) The delta score bin the delta score is in. (example: "0.3-0.5")
    """

    ## iterate over each delta score bin
    for d_score_bin in d_score_bins:

        ## get min and max delta scores for the delta score  bin
        min_d_score = float(d_score_bin.strip().split("-")[0])
        max_d_score = float(d_score_bin.strip().split("-")[1])

        ## If the max delta score is 1, add a small amount to it to account for the last bin
        max_d_score += 0.1 if max_d_score >= 1.0 else 0

        ## check the delta score againts bin values and get the correct bin
        if delta_score >= min_d_score and delta_score < max_d_score:
            delta_score_bin = d_score_bin
            break
    
    return(delta_score_bin)


def get_max_sum_delta_score(spliceai_ano_dict_by_symbol):
    """
    get_max_sum_delta_score
    ===================
    Method to get the MAX SUM of SpliceAI delta scores from a specific ref key in the spliceai_region_dict

    Parameters:
    -----------
    1) spliceai_ano_dict_by_symbol: (dict) A sub-dict of the spliceai_region_dict for a current pos/ref/alt 

    Returns:
    ++++++++
    1) (float) The max delta score
    """
    
    max_sum_score = 0.0 
    ## For each 
    for gene, spliceai_annotation in spliceai_ano_dict_by_symbol.items():

        ## Get the current max delta score
        sum_ds = sum([float(spliceai_annotation["DS_AG"]),
                      float(spliceai_annotation["DS_AL"]),
                      float(spliceai_annotation["DS_DG"]),
                      float(spliceai_annotation["DS_DL"])])

        if float(sum_ds) > float(max_sum_score):
            max_sum_score = float(sum_ds)

    return(max_sum_score)


def get_max_delta_score(spliceai_ano_dict_by_symbol):
    """
    get_max_delta_score
    ===================
    Method to get the max SpliceAI delta score from a specific ref/alt key in the spliceai_region_dict

    Parameters:
    -----------
    1) spliceai_ano_dict_by_symbol: (dict) A sub-dict of the spliceai_region_dict for a current pos/ref/alt 

    Returns:
    ++++++++
    1) (float) The max delta score
    """
    
    max_delta_score = 0.0 
    for gene, spliceai_annotation in spliceai_ano_dict_by_symbol.items():

        ## Get the current max delta score
        max_ds = max(spliceai_annotation["DS_AG"],spliceai_annotation["DS_AL"],spliceai_annotation["DS_DG"],spliceai_annotation["DS_DL"])

        if float(max_ds) > float(max_delta_score):
            max_delta_score = float(max_ds)

    return(max_delta_score)


def output_log(log_string,log_file, new_file = False):
    """  
    output_log
    ==========
    Method used to write info to one of the output logs

    Parameters:
    -----------
    1) log_string: (str): A string to write to the log file
    2) log_file:   (str): The log file to write to
    3) new_file:   (bool): Wether or not to create a new file. Default = False
    """
    
    out_type = "w"
    if not new_file:
        out_type = "a"

    with open(log_file, out_type) as out:
        final_log_string = log_string if log_string.endswith("\n") else log_string + "\n"
        out.write(final_log_string)


def create_interlap_from_ggd_pkg(file_path, zero_based = True, set_cutoff = False, cutoff_col = "NA", cutoff_score = 0.0, chrom_col = "chrom", start_col = "start", end_col = "end" ):
    """
    create_interlap_from_ggd_pkg
    ============================
    Method to create a per chromosome interlap object of genomic coordinates from a ggd pkg file. For example, if 
     you want to create an interlap object of segmental duplications from the 'grch38-segmental-dups-ucsc-v1' ggd
     pkg, you would provide the file path to the grch38-segmental-dups-ucsc-v1 bed file. This method will create an
     interlap object per chromsome for each seg dup that is in the provided file and return a per chromsome dictionary 
     of interlap objects representing the genomic coordinates of segdups.

    WHY? Interlap is a tree based data struture for fast interval look up. If you want to check if something is 
     overlaping a certain feature, like seg dups, you can create the seg dup interlap dict and preform a fast 
     interval overlap lookup without overhead. (It is much faster then a naive approach) 
    
    NOTE: This function is specifically for ggd pkgs that are based on genomic coordinates and it is assumed 
     that ggd pkg files has a header to describe the columns. This header is used to identify the chromosome, 
     start, and end positions of the feature. 

    The start position for each feature will be set to 1-based if it is zero-based, however, the zero_based 
     parameter needs to be set to True. This function does not know if the features are zero or one based unless
     it is specified by the "zero_based" parameter. This is done because vcf files are 1-based, which is commonly
     what these interlap objects are used to overlap with, therefore, it removes the off by 1 error that would 
     occure if not done. If you do not wish to change a zero_based file to one_based, set "zero_based" to False
     and the interlap object start positions will stay zero_based.

    Additionaly, this function can be given a score/match cutoff to use to filter results by. That is, a score can
     be used to set a cutoff of which features are added to the interlap dict. For example, to filter for segdups 
     with a match fraction score above 95%, one would set the cutoff_co parameter to the column for the match
     fraction score and set the cutoff_score to the score that represents 95% (0.95 for values that range from 0 to 1
     or 95 for values that range from 0 to 100). In order to do this the "set_cutoff" parameter must be set to True.

    Parameters:
    ----------
    1) file_path:     (str)  The full path and name of the ggd pkg file to create an interlap object for. This should
                              be a ggd data package representing genomic regions and have a header.
    2) zero_based:    (bool) Whether or no the data file genomic coordinates are zero based. If set to true, the start 
                              position of each region will be set to 1-based
    3) set_cutoff:    (bool) Whether or not to use a cutoff score to filter the features added to the interlap object
                               (Default = False)
    4) cutoff_col:    (str)   The name of the column that represents the cutoff score to use. (Default = "NA")
    5) cutoff_score:  (float) The value to use as the cutoff value. Any value in the cutoff_col >= this value 
                               will be kept, while anything < this value will be removed. (Default = 0.0)
    7) start_col:     (str)   The name of the column that represents the start position (default = "start")
    8) end_col:       (str)   The name of the column that represents the end position (default = "end")

    Returns:
    +++++++
    1) (dict) A dictionary of interlap objects. Keys are the chromsome while values are interlap object representing 
               the feature regions in the data file. If zero_based is set to true, the start position of each interlap 
               feature will be set to one_based by increment each start position by 1. 
    """

    from interlap import InterLap

    
    try:
        ggd_pkg_fh = gzip.open(file_path, "rt", encoding = "utf-8") if file_path.endswith(".gz") else io.open(file_path, "rt", encoding = "utf-8")
    except IOError as e:
        print("\n!!ERROR!! There was a problem reading the ggd pkg file: {}. ".format(file_path))
        print(str(e))
        sys.exit(1)

    assert isinstance(cutoff_score, float), "\n!!ERROR!! The cutoff score must be a float. Provided was: {} of type {}".format(cutoff_score, type(cutoff_score))

    if set_cutoff:
        print("\n\t\t Using a cutoff score of ({}) found in column '{}'".format(cutoff_score, cutoff_col))

    ## Iterate over the file
    interlap_dict = defaultdict(InterLap)
    header_index = []
    for line in ggd_pkg_fh:
        
        ## Get the header 
        if line[0] == "#":
            header_index = line.strip().replace("#","").split("\t")
            continue

        ## Load the next line in the file as a dict, with keys as head values and values as the items in that line
        line_dict = dict(zip(header_index, line.strip().split("\t")))

        ## If the set to check for cutoff
        if set_cutoff:
            
            ## If cutoff is below the cutoff score, skip it
            if float(line_dict[cutoff_col]) < float(cutoff_score):
                continue

        ## Set start to 1 based if zero based
        start = int(line_dict[start_col]) + 1 if zero_based else int(line_dict[start_col])

        ## Add start and end position of an interval to the chromsome specific interlap object
        interlap_dict[line_dict[chrom_col]].add((start, int(line_dict[end_col])))

    ## Close the fh 
    ggd_pkg_fh.close()

    return(interlap_dict)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
