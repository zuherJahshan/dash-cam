import glob
import os
import random
from typing import List

class Reads(object):
    def __init__(self,
                 datadir,
                 kmer_size):
        self.kmer_size = kmer_size
        self.datadir = datadir
        self.read_dir = datadir + "/reads/"
        self.platforms = {
            "illumina": self._getIlluminaReads,
            "pacbio": self._getPBReads,
            "roche": self._get454Reads
        }
        os.makedirs(self.read_dir, exist_ok=True)

    def getPlatforms(self):
        return list(self.platforms.keys())

    def getReads(
        self,
        platform: str,
        genome_filepath,
        reads_num = 10000) -> List[str]:
            if not platform in self.platforms:
                print(f"There is no existing platform named {platform}. Please run Reads.getPlatforms() to get available platforms.")
                return []
            else:
                return self.platforms[platform](genome_filepath, reads_num)

    def _getKmers(
        self,
        s,
        kmer_size
    ):
        if kmer_size > len(s):
            return []

        subsequences = [s[i:i+kmer_size] for i in range(len(s) - kmer_size + 1)]
        return subsequences

    def _getReadsFromFastq(
            self,
            fastq_name
    ):
        with open(fastq_name, 'r') as f:
            return [line for idx, line in enumerate(f.read().split("\n")) if idx % 4 == 1]

    def _getIlluminaReads(
        self,
        genome_file,
        sample_size
    ):
        # Run to generate reads file
        reads_filepath = self.read_dir + "/illumina-" + genome_file.split("/")[-1].split(".")[0]
        if not glob.glob(reads_filepath + "*.fq"):
            os.system(f"art_illumina -sam -i {genome_file} -l {self.kmer_size} -ss HS25 -f 10 -o {reads_filepath} > /dev/null")
        
        read_files = glob.glob(reads_filepath + "*.fq")

        # extract reads from file
        reads = []
        for read_file in read_files:
            reads += self._getReadsFromFastq(read_file)
        return random.sample(reads, min(sample_size, len(reads)))

    def _getPBReads(
        self,
        genome_file,
        sample_size,
        accuracy = 0.9
    ):

        reads_filepath = self.read_dir + "/pacbio-" + genome_file.split("/")[-1].split(".")[0]
        if not glob.glob(reads_filepath + "*.fastq"):
            # Run to generate reads file
            os.system(f"pbsim {genome_file} --hmm_model pbsim2/data/P6C4.model --prefix {reads_filepath} --depth 5 --length-mean {self.kmer_size} --length-sd 0 --accuracy-mean {accuracy} > /dev/null")
        
        # extract reads from file
        reads = []
        for file in glob.glob(reads_filepath + "*.fastq"):
            reads += self._getReadsFromFastq(file)
        return random.sample(reads, min(sample_size, len(reads)))

    def _get454Reads(
        self,
        genome_file,
        sample_size
    ):
        reads_filepath = self.read_dir + "/roche-" + genome_file.split("/")[-1].split(".")[0]
        if not glob.glob(reads_filepath + "*.fq"):
            # Run to generate reads file
            os.system(f"art_454 -s -t {genome_file} {reads_filepath} 10 101 0 > /dev/null")

        reads = []
        for file in glob.glob(reads_filepath + "*.fq"):
            reads += self._getReadsFromFastq(file)
        kmers = []
        for read in reads:
            kmers += self._getKmers(read, kmer_size=self.kmer_size)
        return random.sample(kmers, min(sample_size, len(kmers)))