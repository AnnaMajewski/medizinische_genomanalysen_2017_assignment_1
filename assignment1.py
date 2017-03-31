import mysql.connector
import pysam
import pybedtools

__author__ = 'Anna Majewski'


class Assignment1:
    
    def __init__(self):
        ## Your gene of interest
        self.gene = "SMPD1"
        ## Zuerst habe ich die Information in einer Zeile ausgelesen.
        # Da ich die Infos aber weiter verwenden will, schreibe ich sie in Variablen.
        # Die Namen der Spalten entnehme ich der SQL-query.
        # Diese werden beim Ausfuehren des Query aufgefuellt.
        self.info = []
        self.name = ""
        self.chrom = ""
        self.txStart = ""
        self.txEnd = ""
        self.strand = ""
        self.exonCount = ""
        self.exonStarts = ""
        self.exonEnds = ""

    def fetch_gene_coordinates(self, genome_reference, file_name):
        
        print ("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        ## Execute query
        ergebnis = cursor.execute(query)
        
        ## Write to file
        ## TODO this may need some work 
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
                ## Wenn mein Gen erkannt wird, werden die Informationen in ein Tupel "row" geschrieben
                if row[0] == self.gene:
                    info = row
                    ## Dieses tupel loese ich so auf, dass jeder Teil der Info in einer Variable gespeichert wird.
                    self.gene = info[0]
                    self.name = info[1]
                    self.chrom = info[2]
                    self.txStart = info[3]
                    self.txEnd = info[4]
                    self.strand = info[5]
                    self.exonCount = info[6]
                    self.exonStarts = info[7] # Das ist ein bytestring
                    self.exonEnds = info[8] # Das ist ein bytestring

        #print("Genname:", self.gene)
        #print("Gennname2:", self.name)
        #print("Chromosom:", self.chrom)
        #print("txStart:", self.txStart)
        #print("txEnd:", self.txEnd)
        #print("strand:", self.strand)
        #print("Exone:", self.exonCount)
        #print("Exons:", self.exonStarts)
        #print("Exons:", self.exonEnds)

        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print ("Done fetching data")
                
    def get_sam_header(self):
        # Um ein BAM file lesen zu koennen muss mna ein alignmentfile object erstellen
        # zuerst gebe ich fix hier ein pathfile an, wird spaeter ersetzt
        path_to_file = '/home/taka/medgen/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam'
        samfile = pysam.AlignmentFile(path_to_file, "rb")
        header=[]
        # laut manual enthaelt das first level "HD" oder "SQ".
        for i in samfile.header['SQ']:
            header.append(i)
            print("Header: {}".format(i))

        samfile.close()
        
    def get_properly_paired_reads_of_gene(self):
        print("todo")
        
    def get_gene_reads_with_indels(self):
        print("todo")
        
    def calculate_total_average_coverage(self):
        print("todo")
        
    def calculate_gene_average_coverage(self):
        print("todo")
        
    def get_number_mapped_reads(self):
        print("todo")
        
    def get_gene_symbol(self):
        print("todo")
        
    def get_region_of_gene(self):
        print("todo")
        
    def get_number_of_exons(self):
        print("todo")
    
    def print_summary(self):
        print("Print all results here")
    
        
if __name__ == '__main__':
    print ("Assignment 1")
    assignment1 = Assignment1()
    assignment1.get_sam_header()
    assignment1.print_summary()
    
    assignment1.fetch_gene_coordinates("hg19", "MYFILE.TXT") ## TODO change filename
    
    

