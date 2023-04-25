import copy
import cupy as cp
import csv
import os
import random
from typing import List

class ParallelSearch:
    def __init__(self,
                 kmer_size = 32):
        self.blocks = {}
        self.base_mapping = {'N': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
        self.kmer_size = kmer_size
        self.results = {}
        self.record_results = False

    def seq_to_numeric(self, seq):
        return [self.base_mapping.get(s, 0) for s in seq]

    def buildBlock(self,
                   block_name: str,
                   genome_file: str):
        
        def reverseComplement(sequence):
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            return ''.join(complement[base] for base in reversed(sequence))
        
        if self.record_results:
            return
        # Extract genome from file
        f = open(genome_file, 'r')
        genome = ["".join(fragment.split('\n')[1:]) for fragment in f.read().split('>')[1:]]
        genome = genome + [reverseComplement(fragment) for fragment in genome]
        f.close()

        # Extract kmers from genome
        kmers_in_genome = []
        for fragment in genome:
            kmers_in_fragment_num = len(fragment) - self.kmer_size + 1
            for kmer_idx in range(kmers_in_fragment_num):
                kmers_in_genome.append(fragment[kmer_idx: kmer_idx + self.kmer_size])

        self._buildBlock(block_name=block_name,
                         sequences=kmers_in_genome)

    def search(self,
               patterns,
               threshold = 0,
               time = 0,
               discharge_rate = 0,
               true_genome = 'Nan'):
                
        self._perturbateBlocks(
            time=time,
            discharge_rate=discharge_rate
        )

        # Build results
        results = None
        if (threshold, time, discharge_rate) in self.results:
            results = self.results[(threshold, time, discharge_rate)]
        else:
            results = {
                "threshold": threshold,
                "time": time,
                "discharge rate": discharge_rate
            }
            for block_name in self.blocks:
                results.update({
                        block_name + "-tp": 0,
                        block_name + "-fp": 0,
                        block_name + "-fn": 0,
                        block_name + "-tn": 0
                })
            self.results[(threshold, time, discharge_rate)] = results
        for pattern in patterns:
            if len(pattern) != self.kmer_size:
                continue
            ########### Search Pattern ###########
            pattern_array = cp.array(self.seq_to_numeric(pattern), dtype='int32')
            matching_blocks = []

            for block_name, block in self.blocks.items():
                def baseNeq(b1, b2):
                    return not (b1 == b2 or b1 == 0 or b2 == 0)
                
                patternNeq = cp.vectorize(baseNeq)
                distances = cp.sum(patternNeq(block, pattern_array), axis=1)
                if cp.any(distances <= threshold):
                    matching_blocks.append(block_name)
            ######################################

            ########### Update Results ###########
            for matching_block in matching_blocks:
                if true_genome in matching_block:
                    self.results[(threshold, time, discharge_rate)][matching_block + "-tp"] += 1
                else:
                    self.results[(threshold, time, discharge_rate)][matching_block + "-fp"] += 1
            if not true_genome in matching_blocks and not true_genome == "Nan":
                self.results[(threshold, time, discharge_rate)][true_genome + "-fn"] += 1

            for block_name in self.blocks:
                if block_name != true_genome and not block_name in matching_blocks:
                    self.results[(threshold, time, discharge_rate)][block_name + "-tn"] += 1
            ######################################
        self._resetBlocks()
    
    def recordResults(self,
                      results_file: str):
        os.makedirs("/".join(results_file.split("/")[:-1]), exist_ok=True)
        self.record_results = True
        self.results_file = results_file
        self.results_headers = ["threshold", "time", "discharge rate"]
        for block_name in self.blocks:
            self.results_headers.append(block_name + "-tp")
            self.results_headers.append(block_name + "-fp")
            self.results_headers.append(block_name + "-fn")
            self.results_headers.append(block_name + "-tn")
        self.results = {}


    def stopRecording(self):
        with open(self.results_file, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(self.results_headers)
            for result in self.results.values():
                writer.writerow([result[header] for header in self.results_headers])
        self.record_results = False
        self.results = {}

    def _buildBlock(self, block_name, sequences):
        seq_matrix = cp.array([self.seq_to_numeric(s) for s in sequences], dtype='int32')
        self.blocks[block_name] = seq_matrix

    def _perturbateBlocks(self, 
                          time,
                          discharge_rate,
                          mu = 100,
                          sigma = 16.6,
    ):
        self.orig_blocks = copy.deepcopy(self.blocks)
        if discharge_rate > 0:
            self._readBaseFlips(discharge_rate)
        if time > 0:
            self._timeBaseFlips(time,
                                mu=mu,
                                sigma=sigma)
            
    def _resetBlocks(self):
        self.blocks = copy.deepcopy(self.orig_blocks)

    def _readBaseFlips(
        self,
        discharge_rate):
        for block_name, block in self.blocks.items():
            self.blocks[block_name] = cp.where(
                cp.random.rand(block.shape[0], block.shape[1]) < discharge_rate,
                self.base_mapping["N"],
                block
            )

    def _timeBaseFlips(
        self,
        time,
        mu,
        sigma):
        for block_name, block in self.blocks.items():
            timingOfFlips = cp.random.normal(mu, sigma, block.shape)
            print(time)
            print(timingOfFlips)
            print(cp.where(timingOfFlips < time, self.base_mapping["N"], block))
            self.blocks[block_name] = cp.where(timingOfFlips < time, self.base_mapping["N"], block)
