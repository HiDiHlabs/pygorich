import json
import requests
import pandas as pnd
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
        '''Load and annotate foreground and background genomic regions
        Parameters
        ----------
        foreground_bed_filename : string
            Path to bed file containing foreground regions. Important: Must
            contain header that starts with "#chrom start end".
        background_bed_filename : string
            Path to bed file containing background regions. Important: Must
            contain header that starts with "#chrom start end"
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

        # Retrieve genesets for each library and store in self.__genesets
        for library_id in self.__library_ids[:3]:
            print(library_id)
            query_string="/geneSetLibrary?mode=text&libraryName=%s"
            response = requests.get(self.__enrichr_url + 
                                    query_string % library_id)
            if not response.ok:
                    raise Exception('Error searching for terms')
            for geneset_string in response.text.split("\n"):
                geneset_list = geneset_string.split("\t")
                geneset_id = geneset_list[0]
                geneset = ";".join(geneset_list[2:])
                self.__genesets.loc[library_id+
                                    "@"+
                                    geneset_id, :] = [library_id,
                                                      geneset_id,
                                                      geneset]
