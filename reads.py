import glob
import os
import random
from typing import List

class Reads(object):
    def __init__(self,
                 datadir,
                 read_length):
        self.read_length = read_length
        self.datadir = datadir
        self.read_dir = datadir + "/reads/"
        self.platforms = {
            "illumina": {
                "sequencer_command": lambda genome_filepath, read_len, reads_filepath, accuracy, coverage:\
                    f"art_illumina -sam -i {genome_filepath} -l {read_len} -ss HS25 -f {coverage} -o {reads_filepath} > /dev/null",
                "read_file_suffix": ".fq"
            },
            "pacbio": {
                "sequencer_command": lambda genome_filepath, read_len, reads_filepath, accuracy, coverage:\
                    f"pbsim {genome_filepath} --hmm_model pbsim2/data/P6C4.model --prefix {reads_filepath} \
                    --depth {coverage} --length-mean {read_len} --length-sd 0 --accuracy-mean {accuracy} > /dev/null",
                "read_file_suffix": ".fastq"
            },
            "roche": {
                "sequencer_command": lambda genome_filepath, read_len, reads_filepath, accuracy, coverage:\
                    f"art_454 -s -t {genome_filepath} {reads_filepath} {coverage} {max(101, read_len)} 0 > /dev/null",
                "read_file_suffix": ".fq"
            },
        }
        os.makedirs(self.read_dir, exist_ok=True)

    
    def getPlatforms(self):
        return list(self.platforms.keys())

    def getReadFiles(self,
                     platform,
                     genome_filepath) -> List[str]:
        if not self._isPlatformValid(platform):
            return []
        reads_filepath = f"{self.read_dir}/{platform}-{genome_filepath.split('/')[-1].split('.')[0]}"
        if not glob.glob(f"{reads_filepath}*{self.platforms[platform]['read_file_suffix']}"):
            os.system(self.platforms[platform]['sequencer_command'](genome_filepath, self.read_length, reads_filepath, 0.9, 10))
        return glob.glob(f"{reads_filepath}*{self.platforms[platform]['read_file_suffix']}")

    def getReads(self,
                 platform,
                 genome_filepath) -> List[str]:
        if not self._isPlatformValid(platform):
            return []
        else:
            read_filepaths = self.getReadFiles(platform, genome_filepath)
            reads = []
            for file in read_filepaths:
                reads += self._getReadsFromFastq(file)
            return reads
        
    def getKmers(
        self,
        platform: str,
        genome_filepath,
        kmer_size,
        kmers_num = 10000) -> List[str]:
            kmers = []
            for read in self.getReads(platform=platform, genome_filepath=genome_filepath):
                kmers += self._getKmers(read, kmer_size=kmer_size)
            return random.sample(kmers, min(kmers_num, len(kmers)))

    def _isPlatformValid(self, platform):
        if not platform in self.platforms:
            print(f"There is no existing platform named {platform}. Please run Reads.getPlatforms() to get available platforms.")
            return False
        return True

    
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