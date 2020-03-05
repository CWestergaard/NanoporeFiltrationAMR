"""
Program: NanoporeFiltrationAMR
Description:Run KMA on input fastq.gz file and group resulting gene-hits into
            gene-groups. Apply filtrations on gene-groups based on highest 
            depth gene in each gene-group, as well as a minimum depth per 
            gene-group relative to the average genome coverage of the input file.
            To only include genes that are assumed to be functional, a cutoff
            for identity and coverage is applied as well (Default 90).
Version: 1.0
Author: Casper Westergaard
"""

# Import libraries
import gzip
import sys
import os
import subprocess
import csv
import argparse

###########################################################################
# FUNCTIONS
########################################################################### 

class Gene:
    """
        Gene class - Save information about Antimicrobial class, predicted phenotype and group
    """
    def __init__(self, geneClass, genePhenotype, geneGroup):
        self.geneClass = geneClass
        self.genePhenotype = genePhenotype
        self.geneGroup = geneGroup
        

def genomeSize_convert(genomeSize_string):
    """
        Takes a genome size as string and converts it to int.
        
        parameters:
            genomeSize_string = Size of genome as bp, Kbp, Mbp or Gbp.
        returns:
            genomeSize = Size of genome as int.
    """    
    # Convert input genome size to int
    if genomeSize_string[-1].upper() == 'K':
        genomeSize = int(genomeSize_string[0:-1]) * 1000
    elif genomeSize_string[-1].upper() == 'M':
        genomeSize = int(genomeSize_string[0:-1]) * 1000000
    elif genomeSize_string[-1].upper() == 'G':
        genomeSize = int(genomeSize_string[0:-1]) * 1000000000
    else:
        genomeSize = int(genomeSize_string)
    
    return genomeSize
    
def BaseCount(input_file):
    """"
        Counts the number of bases and reads in the input file.
        
        parameters:
            input_file: path to fastq.gz input file.
        returns:
            readCount: Count of reads in input file.
            baseCount: Count of bases in input file.
    """
    readCount = 0
    baseCount = 0
    with gzip.open(input_file, 'r') as reads:
        
        for header in reads:
            readCount += 1
            
            if header.startswith(b'@') and b'start_time' in header:  
                seq = next(reads)
                baseCount += len(seq.strip())
                next(reads)
                next(reads)               
            else:
                sys.exit('Reads should include 4 lines, starting with a header.'
                         'Read header of read {0} does not look as expected.\n{1}'
                         .format(readCount,header))   
                
    return baseCount,readCount

def runKMA(input_file, output_file, refdb):
    """
        Run KMA on the input file using the supplied reference database.
        
        parameters:
            input_file: path to fastq.gz input file to run KMA on.
            output_file: path and filename for KMA output.
            refdb: path to KMA indexed reference database.
    """
    os.makedirs(output_folder, exist_ok=True)
    cmd_kma = ['kma', '-i', input_file, '-o', output_file,
               '-t_db', refdb, '-mem_mode', '-mp', '20', '-bcNano']
    subprocess.run(cmd_kma)   

def createPhenotypeDB(phenotypesFile):
    """
        Create dict of gene classes with antimicrobial class, phenotype and gene group as values.
        
        Parameters:
            phenotypesFile = phenotypes.txt file from ResFinder database.
        return:
            genePhenotypes_db = dict of gene classes with gene name as key and
                                antimicrobial class, phenotype and gene group as values.
    """
    genePhenotypes_db = dict()
    with open(phenotypesFile, 'r') as genePhenotypes:
        genePhenotypes.readline()    # Skip header
        for line in genePhenotypes:
            line = line.split('\t')
            geneName = line[0].strip()
            geneClass = line[1].strip().upper()
            genePhenotype = line[2].strip()
            geneGroup = getGeneGroup(geneName)
            genePhenotypes_db[geneName] = Gene(geneClass, genePhenotype, geneGroup)
    
    return genePhenotypes_db

def getGenes(csvfile, coverageCutoff, identityCutoff):
    """ 
        Object stores a dict of genes that is above the given identity and
        coverage thresholds.
        Each gene is a string created from gene name and accession no.
        
        parameters:
            csvfile = Input csv file generated from KMA.
            coverageCutoff = Minimum coverage value to include genes.
            identityCutoff = Minimum identity value to include genes
        returns:
            genes = Dictionary with gene name and accession no. as key and
                    a list of alignment statistics as value.
    """
    genes = dict()

    with open(csvfile, 'r') as results_file:

        reader = csv.reader(results_file, delimiter='\t')
        # ignore first row
        next(reader)

        for row in reader:
            templateIdentity = float(row[4])
            templateCoverage = float(row[5])
            queryIdentity = float(row[6])
            queryCoverage = float(row[7])
            depth = float(row[8])
            if '~~~' in row[0]:
                geneName = row[0].split('~~~')[1]
                accesion = row[0].split('~~~')[2]
                gene = '{}_{}'.format(geneName, accesion)
            else:
                gene = row[0]

            if (templateIdentity >= identityCutoff and templateCoverage >= coverageCutoff \
            and queryIdentity >= identityCutoff and queryCoverage >= coverageCutoff):
                genes[gene] = [templateIdentity,templateCoverage,queryIdentity,queryCoverage,depth]

    return genes
        
def getGeneGroup(geneName):
    """ 
        Extract the gene-group from the gene name and returns it.
        
        parameters:
            geneName: Name of gene.
        returns:
            geneGroup: Gene-group based on name of gene.
    """
    geneGroup = geneName[0:3]

    if geneGroup == 'bla':
        if len(geneName.split('-')) == 1:

            if geneName[0:6] == 'blaBEL':
                geneGroup = 'blaBEL'
            else:
                geneGroup = geneName.split('_')[0]

        else:
            geneGroup = geneName.split('-')[0]

    return geneGroup

def addToDict(geneDict, geneName, geneStats, geneClass, genePhenotype, geneGroup):
    """
        Adds genes to input dictionary.
        geneClass{}: genePhenotype{}: geneGroup{}: [geneName, geneStats]
        
        parameters:
            geneDict = Dictionary to add to.
            geneName = Name of gene.
            geneStats = Alignment statistics of gene from KMA.
            geneClass = Predicted antimicrobial class that gene confers resistance against.
            genePhenotype = Predicted phenotype of gene.
            geneGroup = Gene-group of gene.
        returns:
            geneDict = Updated input dict.        
    """
    geneDict.setdefault(geneClass, {})
    geneDict[geneClass].setdefault(genePhenotype, {})
    geneDict[geneClass][genePhenotype].setdefault(geneGroup, []).append([geneName, geneStats])
    
    return geneDict

###########################################################################
# GET INPUT
###########################################################################
    
# Input from commandline
# Required input
parser = argparse.ArgumentParser(description='Apply filtrations to KMA results and report results.'
                                 'Filtrations are meant to reduce hits from contamination and sequencing errors in Nanopore data.'
                                 'Genes remaining after the filtrations are listed with their predicted phenotype in the output.')

parser.add_argument('-o', type=str, dest='output_folder',
                    help='Path to output folder', required=True)

parser.add_argument('-i', type=str, dest='input_file',
                    help='Path to input file', required=True)

parser.add_argument('-pf', type=str, dest='phenotypesFile',
                    help='Paht to input file containing name and phenotype of all genes in ResFinder database', required=True)

parser.add_argument('-gs', type=str, dest='genomeSize_string', 
                    help='Size of genome', required=True)

parser.add_argument('-refdb', type=str, dest='refdb',
                    help='Path to KMA indexed ResFinder database', required=True)

# Optional input
parser.add_argument('-gd', type=float, dest='groupDepth', default='0.45',
                    help='Set depth filtration, default 0.45, input value between 0 and 1', required=False)

parser.add_argument('-md', type=float, dest='minDepth', default='0.2',
                    help='Minimum depth for gene-groups to be included, default 0.2, input value between 0 and 1', required=False)

parser.add_argument('-id', type=float, dest='identityCutoff', default='90',
                    help='Minimum identity for genes to be included, default 90, input value between 0 and 100', required=False)

parser.add_argument('-cov', type=float, dest='coverageCutoff', default='90',
                    help='Minimum coverage for genes to be included, default 90, input value between 0 and 100', required=False)

args = parser.parse_args()

###########################################################################
# CHECK INPUT
###########################################################################

input_file = args.input_file
output_folder = args.output_folder
refdb = args.refdb
phenotypesFile = args.phenotypesFile
genomeSize_string = args.genomeSize_string
minDepth = args.minDepth
groupDepth = args.groupDepth
coverageCutoff = args.coverageCutoff
identityCutoff = args.identityCutoff
filename = os.path.basename(input_file)

# Create output folder
os.makedirs('{}'.format(output_folder), exist_ok=True)

# Convert input genome size to int
genomeSize = genomeSize_convert(genomeSize_string)

###########################################################################
# COUNT BASES
###########################################################################

# Count bases and generate minimum depth cutoff based on average genome coverage
baseCount,readCount = BaseCount(input_file)
coverage = baseCount/genomeSize
minDepthCutoff = coverage*minDepth

###########################################################################
# CREATE PHENOTYPIC GENE DATABASE
###########################################################################

# Create gene class containing Antibiotic class, Phenotype and Gene-group for each gene
# Based on phenotypes.txt from the ResFinder database
genePhenotypes_db = createPhenotypeDB(phenotypesFile)

###########################################################################
# RUN KMA THEN SAVE AND GROUP VALID RESULTS
###########################################################################

# Run KMA
kma_filename = output_folder+'/kma_resfinder.'+filename
runKMA(input_file, kma_filename, refdb)

# Get results from KMA above coverage and identity thresholds
validGenes = getGenes(kma_filename+'.res', coverageCutoff, identityCutoff)
    
# Gather the valid genes in geneGroups based on their geneGroup from genePhenotypes_db
geneGroups = dict()
for geneName in validGenes:
    if geneName in genePhenotypes_db:
        geneInfo = genePhenotypes_db[geneName]
        geneGroups.setdefault(geneInfo.geneGroup, []).append([geneName, validGenes[geneName]])        
    else:
        print('{0} is found in KMA results, but is not found in the phenotype'
              'translation file. The gene has been discarded from the analysis.'
              .format(geneName))
      
###########################################################################
# DEPTH FILTRATION ON NANOPORE DATA
########################################################################### 
            
# Sort genes in geneGroups based on depth, and save total depth of each gene-group in geneGroup_depths
geneGroup_depths = dict() 
for geneGroup in geneGroups:
    geneGroups[geneGroup].sort(key=lambda x: x[1][4],reverse = True)    # Sort genes based on depth
    geneGroup_depths[geneGroup] = sum([geneGroups[geneGroup][x][1][4] for x in range(len(geneGroups[geneGroup]))])

# Depth filtration    
filteredGenes = dict()
removedGenes = dict()
for geneGroup in geneGroups:
    
    # Minimum depth
    if geneGroup_depths[geneGroup] < minDepthCutoff:
        for gene in geneGroups[geneGroup]:
            geneName = gene[0]
            geneStats = gene[1]
            geneInfo = genePhenotypes_db[geneName]
            removedGenes = addToDict(removedGenes, geneName, geneStats, geneInfo.geneClass, geneInfo.genePhenotype, geneInfo.geneGroup)
            
    # Gene-group depth
    else:
        geneGroup_highestDepth = geneGroups[geneGroup][0][1][4]
        geneGroup_depthCutoff = groupDepth*geneGroup_highestDepth
        
        for gene in geneGroups[geneGroup]:
            geneName = gene[0]
            geneStats = gene[1]
            geneDepth = geneStats[4] 
            geneInfo = genePhenotypes_db[geneName]
            
            if geneDepth < geneGroup_depthCutoff:
                removedGenes = addToDict(removedGenes, geneName, geneStats, geneInfo.geneClass, geneInfo.genePhenotype, geneInfo.geneGroup)
                
            else:
                filteredGenes = addToDict(filteredGenes, geneName, geneStats, geneInfo.geneClass, geneInfo.genePhenotype, geneInfo.geneGroup)
     
###########################################################################
# OUTPUT RESULTS
###########################################################################    

cutoffString = ('Input file: {0}\nCutoffs for including genes are:\nIdentity >= {1}\n'
                'Coverage >= {2}\nDepth(Relative to highest depth gene from that gene-group) >= {3}%'
                '\nMinimum gene-group depth (Relative to avg. genome coverage) >= {4}%\n'
                .format(filename,identityCutoff,coverageCutoff,100*groupDepth,100*minDepth))

# Output effect of depth filtration
filtrationFilename = '{0}/{1}_Depth={2}_Min_Depth={3}_filtration.txt'.format(output_folder,filename,groupDepth,minDepth)
with open(filtrationFilename, 'w') as outfile:
    # Output filtration parameters
    outfile.write(cutoffString)

    # Output results of filtration
    for geneGroup in geneGroups:
        geneGroupDepth = geneGroup_depths[geneGroup]
        outfile.write('\n# {0} total depth: {1}'
                      .format(geneGroup,round(geneGroupDepth,2)))
        
        # Discarded gene-groups due to minDepth
        if geneGroupDepth < minDepthCutoff:
            outfile.write(' - Removed due to minimum depth cutoff.\n')
            outfile.write('Discarded {0} genes:\n'.format(geneGroup))
            for gene in geneGroups[geneGroup]:
                geneName = gene[0]
                geneDepth = gene[1][4]
                outfile.write('{0}: {1}\n'.format(geneName,geneDepth))   

        # Included gene-groups
        else:
            geneGroup_highestDepth_name = geneGroups[geneGroup][0][0]
            geneGroup_highestDepth = geneGroups[geneGroup][0][1][4]
            geneGroup_depthCutoff = groupDepth*geneGroup_highestDepth
            outfile.write('\nHighest depth from gene-group {0} is: {1} ({2})'
                          .format(geneGroup,geneGroup_highestDepth,geneGroup_highestDepth_name))          
            outfile.write('\nDepth threshold for including genes from gene-group {0} is: {1}'
                          .format(geneGroup,round(geneGroup_depthCutoff,2)))
            
            # Genes included by gene-group depth filtration
            outfile.write('\n\nIncluded {0}-genes:\n'.format(geneGroup))
            for gene in geneGroups[geneGroup]:
                geneName = gene[0]
                geneDepth = gene[1][4]
                if geneDepth >= geneGroup_depthCutoff:
                    outfile.write('{0}: {1}\n'.format(geneName,geneDepth)) 
                else:
                    break
                
            # Genes discarded by gene-group depth filtration
            outfile.write('\nDiscarded {0}-genes:\n'.format(geneGroup))
            for gene in geneGroups[geneGroup]:
                geneName = gene[0]
                geneDepth = gene[1][4]
                if geneDepth < geneGroup_depthCutoff:
                    outfile.write('{0}: {1}\n'.format(geneName,geneDepth)) 

                     
# Output results after filtration
resultsFilename = '{0}/{1}_Depth={2}_Min_Depth={3}_results.txt'.format(output_folder,filename,groupDepth,minDepth)
with open(resultsFilename, 'w') as outfile:
    outfile.write(cutoffString)
    outfile.write('\nAverage genome coverage: {0} (Based on a genome size of: {1})\n'.format(round(coverage,2),genomeSize_string))
    if coverage < 40:
        outfile.write('Coverage is below 40 - Results might not be accurate.')
    else:
        outfile.write('Coverage is above 40 - Results are expected to be accurate.')
    outfile.write('\n\nGene information:\nGene-template - Template-Identity\tTemplate-Coverage\tQuery-Identity\tQuery-Coverage\tDepth\n')
    
    for geneClass in filteredGenes:
        outfile.write('\n###Antimicrobial class - {0}'.format(geneClass))
                      
        for genePhenotype in filteredGenes[geneClass]:
            outfile.write('\n##Phenotype - {0}'.format(genePhenotype))
                          
            for geneGroup in filteredGenes[geneClass][genePhenotype]:
                outfile.write('\n#Gene-group - {0} - Total depth before filtration, across all phenotypes: {1}\n'
                              .format(geneGroup,round(geneGroup_depths[geneGroup],2)))
                
                for gene in filteredGenes[geneClass][genePhenotype][geneGroup]:
                    geneTemplate = gene[0]
                    templateIdentity = gene[1][0]
                    templateCoverage = gene[1][1]
                    queryIdentity = gene[1][2]
                    queryCoverage = gene[1][3]
                    depth = gene[1][4]
                    outfile.write('{0} - {1}\t{2}\t{3}\t{4}\t{5}\n'
                                  .format(geneTemplate,templateIdentity,templateCoverage,queryIdentity,queryCoverage,depth))
