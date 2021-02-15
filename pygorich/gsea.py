import json
import requests
import pandas as pnd
from scipy.stats import hypergeom, fisher_exact, binom_test
from statsmodels.stats.multitest import multipletests
import sys
sys.path += ["/home/biegm/projects/python_packages"]
import pyanno

class Enricher():
    '''
    Class for managing gene sets and performing GSEA.

    ...

    Attributes
    ----------
    __enrichr_url : string
        URL of EnrichR API

    __genesets : pandas.DataFrame
        Stores information about library id, geneset id, and geneset

    __library_ids : list<string>
        List of available libraries IDs.

    __annotated_foreground : pandas.DataFrame
        Stores Information about the foreground regions.

    __annotated_background : pandas.DataFrame
        Stores information about the background regions.

    __annotation_database : pandas.DataFrame
        Stores information about the gene bed file, and enhancer-promoter
        link file used for annotating the foreground and the background 
        genomic regions.

    Methods
    -------
    printLibraries()
        Method for printing available geneset libraries.

    loadAnnotationDatabase(genes_filename=None,
                           enhancer_link_filename=None,
                           max_distance_gene=1000000,
                           name_col_gene=6,
                           max_distance_enhancer=0,
                           name_col_enhancer=15)
        Loads annotation database used for annotating foregound and
        background regions.

    loadRegions(self,
               foreground_bed_filename = None,
               background_bed_filename = None)
        Load and annotate foreground and background genomic regions.

    loadLibrary(library_id,
                library_filename=None,
                from_enrichr=True)
        Load genesets from library.

    getAnnotationDatabase()
        Returns annotation database.

    getAnnotatedForground()
        Returns annotated foreground.

    getAnnotatedBackground()
        Returns annotated background.

    getGeneSets()
        Returns genesets
    '''

    ##############
    # Constructors
    def __init__(self):
        '''Standard constructor. Creates an empty Enricher object.
        '''
        # Define private Attributes
        self.__enrichr_url = "http://amp.pharm.mssm.edu/Enrichr"

        self.__genesets = pnd.DataFrame(columns = ["LIB.ID", 
                                                   "GENESET.ID", 
                                                   "GENE.LIST"])

        self.__library_ids = None

        self.__annotated_foreground = None

        self.__annotated_background = None
        self.__load_genesets_from_enrichr()

        self.__annotation_database = None

        self.__enrichment_results = pnd.DataFrame(columns = 
                                                  ["LIB.ID",
                                                   "GENESET.ID",
                                                   "n.FOREGROUND",
                                                   "n.BACKGROUND",
                                                   "n.FOREGROUND.IN.SET",
                                                   "n.BACKGROUND.IN.SET",
                                                   "p.FISHER",
                                                   "q.FISHER",
                                                   "odds.FISHER",
                                                   "p.HYPERGEOMETRIC",
                                                   "q.HYPERGEOMETRIC",
                                                   "p.BINOMIAL",
                                                   "q.BINOMIAL",
                                                   "REGION.GENES.PAIRS"])

    ################
    # Public Methods

    # Print Methods
    def printLibraries(self):
        '''Method for printing available geneset libraries.
        '''
        print("\n".join(self.__library_ids))

    # Getter Methods
    def getAnnotationDatabase(self):
        '''Returns annotation database.
        Returns
        -------
        annotation_database : pandas.DataFrame
            DataFrame storing information about the database used for
            annotation of foreground and background regions.
        '''
        return self.__annotation_database

    def getAnnotatedForeground(self):
        '''Returns annotated foreground.
        Returns
        -------
        annotated_foreground : pandas.DataFrame
            DataFrame storing the annotated foreground regions.
        '''
        return self.__annotated_foreground

    def getAnnotatedBackground(self):
        '''Returns annotated background.
        Returns
        -------
        annotated_background : pandas.DataFrame
            DataFrame storing the annotated background regions.
        '''
        return self.__annotated_background

    def getGeneSets(self):
        '''Returns genesets
        Returns
        -------
        genesets : pandas.DataFrame
            DataFrame storing loaded genesets
        '''
        return self.__genesets

    def getEnrichmentResults(self):
        '''Returns enrichment results.
        Returns
        -------
        enrichment_results : pandas.DataFrame
            DataFrame storing enrichment results.
        '''
        return self.__enrichment_results

    def getLibraryIDs(self):
        '''Return Library IDs
        Returns
        -------
        library_ids : list
            List of library IDs
        '''
        return self.__library_ids

    # Setter methods
    def resetEnrichmentResults(self):
        '''Reset enrichment results to empty DataFrame.
        '''
        self.__enrichment_results = pnd.DataFrame(columns = 
                                                  ["LIB.ID",
                                                   "GENESET.ID",
                                                   "n.FOREGROUND",
                                                   "n.BACKGROUND",
                                                   "n.FOREGROUND.IN.SET",
                                                   "n.BACKGROUND.IN.SET",
                                                   "p.FISHER",
                                                   "q.FISHER",
                                                   "odds.FISHER",
                                                   "p.HYPERGEOMETRIC",
                                                   "q.HYPERGEOMETRIC",
                                                   "p.BINOMIAL",
                                                   "q.BINOMIAL",
                                                   "REGION.GENES.PAIRS"])


    # Data load Methods
    def loadAnnotationDatabase(self,
                               genes_filename=None,
                               enhancer_link_filename=None,
                               max_distance_gene=1000000,
                               name_col_gene=6,
                               max_distance_enhancer=0,
                               name_col_enhancer=15):
        '''Load Annotation Database for foreground and background regions.
        Parameters
        ----------
        genes_filename : string
            Path to bed file containing gene regions.
        enhancer_link_filename : string
            Path to bed file containing enhancer promoter interactions.
        max_distance_gene : int
            Maximal distance in base pairs to TSS for annotating a region
            to a gene.
        name_col_gene : int
            Column of gene name in genes_filename.
        max_distance_enhancer : int
            Maximal distance of region to enhancer for annotating a region via
            an enhancer to a gene.
        name_col_enhancer : int
            Column of gene name in enhancer promoter link file.
        '''

        self.__annotation_database = pnd.DataFrame(columns = ["FILENAME",
                                                              "REGION.TYPE",
                                                              "SOURCE",
                                                              "ANNOTATION.BY",
                                                              "MAX.DISTANCE",
                                                              "DISTANCE.TO",
                                                              "N.HITS",
                                                              "NAME.COL"])
        if(genes_filename is not None):
            self.__annotation_database.loc["GENES", :] = [genes_filename,
                                                          "genes",
                                                          "genes",
                                                          "NAME",
                                                          max_distance_gene,
                                                          "START",
                                                          2,
                                                          name_col_gene]
        if(enhancer_link_filename is not None):
            self.__annotation_database.loc["ENHANCER",
                                           :] = [enhancer_link_filename,
                                                 "enhancer",
                                                 "enhancer",
                                                 "NAME",
                                                 max_distance_enhancer,
                                                 "REGION",
                                                 1,
                                                 name_col_enhancer]

    def loadRegions(self,
                    foreground_bed_filename = None,
                    background_bed_filename = None):
        '''Load and annotate foreground and background genomic regions.
        Parameters
        ----------
        foreground_bed_filename : string
            Path to bed file containing foreground regions. Important: Must
            contain header that starts with "#chrom start end".
        background_bed_filename : string
            Path to bed file containing background regions. Important: Must
            contain header that starts with "#chrom start end".
        '''
        # Load and annotated foreground regions
        # Create a new GenomicRegionAnnotator instance
        gra = pyanno.Annotator.GenomicRegionAnnotator()

        # load base
        gra.load_base_from_file(foreground_bed_filename)

        # load database
        gra.load_database_from_dataframe(self.__annotation_database)

        # Annotate base against all database genomic region files
        gra.annotate()

        # Retrieve annotated base intervals as pandas.DataFrame instance
        self.__annotated_foreground = gra.get_base()

        self.__annotated_foreground.loc[:, "GENESET"] = [
        ";".join([ gene.split("(")[0] for gene in row["enhancer"].split(";") ]) 
        if row["enhancer"] != "NA" 
        else row["genes"].split(";")[0].split("(")[0] 
        for i, row in self.__annotated_foreground.iterrows()]

        # Load and annotated background regions
        # Create a new GenomicRegionAnnotator instance
        gra = pyanno.Annotator.GenomicRegionAnnotator()

        # load base
        gra.load_base_from_file(background_bed_filename)

        # load database
        gra.load_database_from_dataframe(self.__annotation_database)

        # Annotate base against all database genomic region files
        gra.annotate()

        # Retrieve annotated base intervals as pandas.DataFrame instance
        self.__annotated_background = gra.get_base()

        self.__annotated_background.loc[:, "GENESET"] = [
        ";".join([ gene.split("(")[0] for gene in row["enhancer"].split(";") ])
        if row["enhancer"] != "NA" 
        else row["genes"].split(";")[0].split("(")[0] 
        for i, row in self.__annotated_background.iterrows()]

    def loadLibrary(self,
                    library_id,
                    library_filename=None,
                    from_enrichr=True):
        '''Load genesets from library
        Parameters
        ----------
        library_id : string
            Identifier of geneset library.
        library_filename : string
            Path to geneset library file. Library file must contain one header
            line! The columns are: 1. Library ID, 2. Geneset ID, 3. Semicolon
            separated list of gene ids.
        from_enrichr : boolean
            Whether geneset library shall be loaded from enrichr database or
            not.
        '''
        if(from_enrichr):
            # Retrieve genesets for each library and store in self.__genesets
            query_string="/geneSetLibrary?mode=text&libraryName=%s"
            response = requests.get(self.__enrichr_url + 
                                    query_string % library_id)
            if not response.ok:
                    raise Exception('Error searching for terms')
            for geneset_string in response.text.split("\n"):
                geneset_list = geneset_string.split("\t")
                geneset_id = geneset_list[0]
                geneset = ";".join([gene.split(",")[0] for 
                                    gene in geneset_list[2:]])
                self.__genesets.loc[library_id+
                                    "@"+
                                    geneset_id, :] = [library_id,
                                                      geneset_id,
                                                      geneset]
        else:
            if(not library_filename is None):
                library_file = open(library_filename, "r")
                c = 0
                for line in library_file:
                    if(c == 0):
                        c += 1
                        continue
                    split_line = line.rstrip().split("\t")
                    if(len(split_line) < 3):
                        continue
                    lib_id = split_line[0]
                    if(not(lib_id == library_id)):
                        continue
                    geneset_id = split_line[1]
                    genelist = split_line[2]
                    self.__genesets.loc[lib_id+
                                        "@"+
                                        geneset_id, :] = [lib_id,
                                                          geneset_id,
                                                          genelist]


    # Enrichment Methods
    def enrich(self,
               method="all"):
        '''Perform Enrichment Analysis.
        Parameters
        ----------
        method : string
            Statistical method used for enrichment analysis. Can be either of
            fisher, hypergeometric, binomial, all.
        '''
        geneset_list_foreground_regions = [ set(geneset.split(";")) for
                                        geneset in 
                                        self.__annotated_foreground["GENESET"]]
        foreground_region_ids = list(self.__annotated_foreground.index)
        geneset_list_background_regions = [ set(geneset.split(";")) for
                                        geneset in 
                                        self.__annotated_background["GENESET"]]

        n_foreground = len(self.__annotated_foreground.index)
        n_background = len(self.__annotated_background.index)

        for library_id in list(set(self.__genesets.loc[:, "LIB.ID"])):
            print("Calculate Enrichment for: "+library_id)
            for geneset_id in self.__genesets[self.__genesets["LIB.ID"] == 
                                          library_id].loc[:, "GENESET.ID"]:
                if(not library_id+"@"+geneset_id in 
                    set(self.__enrichment_results.index)):
#                    print("\t"+geneset_id)
                    geneset = set(self.__genesets.loc[library_id+
                                                      "@"+
                                                      geneset_id,
                                                      "GENE.LIST"].split(";"))
                    (n_foreground_in_geneset, 
                     foreground_region_genes_pairs) = self.__calculateGenesetOverlaps(
                                            geneset_list_foreground_regions,
                                            geneset,
                                            region_id_list = foreground_region_ids)
                    (n_background_in_geneset,
                     background_region_genes_pairs) = self.__calculateGenesetOverlaps(
                                            geneset_list_background_regions,
                                            geneset)

                    p_val_fisher = None
                    odds_fisher = None
                    p_val_hyper = None
                    p_val_binom = None
                    if(method == "fisher" or method == "all"):
                        n_foreground_not_in_geneset = (n_foreground-
                                                       n_foreground_in_geneset)
                        n_background_not_in_geneset = (n_background-
                                                       n_background_in_geneset)
                        ct = [ [n_foreground_in_geneset,
                                n_foreground_not_in_geneset], 
                               [n_background_in_geneset,
                                n_background_not_in_geneset] ]
                        odds_fisher, p_val_fisher = fisher_exact(ct, 
                                                    alternative = "greater")
                    if(method == "hypergeometric" or method == "all"):
                        M = n_background
                        n = n_background_in_geneset
                        N = n_foreground
                        k = n_foreground_in_geneset
                        p_val_hyper = 1.
                        if(n > 0 and k > 0):
                            p_val_hyper = 1.-hypergeom.cdf(k, M, n, N)
                    if(method == "binomial" or method == "all"):
                        p = float(n_background_in_geneset)/float(n_background)
                        x = n_foreground_in_geneset
                        n = n_foreground
                        p_val_binom = binom_test(x, 
                                             n = n, 
                                             p = p, 
                                             alternative="greater")

                results = [library_id,
                           geneset_id,
                           n_foreground,
                           n_background,
                           n_foreground_in_geneset,
                           n_background_in_geneset,
                           p_val_fisher,
                           None,
                           odds_fisher,
                           p_val_hyper,
                           None,
                           p_val_binom,
                           None,
                           ";".join(foreground_region_genes_pairs)]
                self.__enrichment_results.loc[
                    library_id+"@"+geneset_id, :] = results

            # Perform multiple testing correction
            indices = self.__enrichment_results[
                self.__enrichment_results["LIB.ID"] == library_id].index

            if(method == "fisher" or method == "all"):
                p_values = self.__enrichment_results.loc[indices, "p.FISHER"]
                r, q_values, a_s, a_b = multipletests(p_values, 
                                                      method = "fdr_bh")
                self.__enrichment_results.loc[indices, "q.FISHER"] = q_values
            if(method == "hypergeometric" or method == "all"):
                p_values = self.__enrichment_results.loc[indices, 
                                                         "p.HYPERGEOMETRIC"]
                r, q_values, a_s, a_b = multipletests(p_values, 
                                                      method = "fdr_bh")
                self.__enrichment_results.loc[indices, 
                                              "q.HYPERGEOMETRIC"] = q_values
            if(method == "binomial" or method == "all"):
                p_values = self.__enrichment_results.loc[indices, "p.BINOMIAL"]
                r, q_values, a_s, a_b = multipletests(p_values, 
                                                      method = "fdr_bh")
                self.__enrichment_results.loc[indices, "q.BINOMIAL"] = q_values

    #################
    # Private Methods
    def __load_genesets_from_enrichr(self):
        '''Load available gene set libraries from EnrichR.
        '''
        # Retrieve set of libraries
        lib_url = self.__enrichr_url+"/datasetStatistics"
        lib_json = json.loads(requests.get(lib_url).text)
        self.__library_ids = [lib["libraryName"] for 
                              lib in lib_json["statistics"]]

    def __calculateGenesetOverlaps(self,
                                   region_associated_geneset_list,
                                   geneset,
                                   region_id_list = None):
        '''Calculate number of regions in geneset
        Parameters
        ----------
        region_associated_geneset_list : list<set>
            List of genesets the regions are associated with
        geneset : set
            Set of genes

        Returns
        -------
        overlaps : int
            Number of regions in geneset
        '''
        overlaps = 0
        i = 0
        region_genes_pairs = []
        for region_associated_geneset in region_associated_geneset_list:
            overlap_genes = list(region_associated_geneset & geneset)
            if(len(overlap_genes) > 0):
                overlaps += 1
                if(not region_id_list is None):
                    region_id = region_id_list[i]
                    region_genes_pairs += [region_id+
                                           "="+
                                           ",".join(overlap_genes)]
            i += 1
        return overlaps, region_genes_pairs

