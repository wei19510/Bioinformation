import argparse
import pandas as pd
import configparser
import logging
import os
import subprocess
from subprocess import PIPE
import sys
import pysam
import gzip
from Bio import SeqIO
import shutil


class SampleCodeClass:
    def __init__(
            self, indir, q, fada, rada,
            F_remove, R_remove, t, min_len, outdir,
            m, M, N, r, pop_opts):
        self.args = None
        self.indir = indir
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.log = os.path.join(outdir, 'log.txt')
        self.setup_logger()
        self.q = q
        self.F_remove = F_remove
        self.R_remove = R_remove
        self.fada = fada
        self.rada = rada
        self.t = t
        self.min_len = min_len
        self.path_to_bbmap = "/usr/local/src/bbmap"
        self.path_to_bwa = "bwa"
        self.path_to_samtools = "samtools"
        self.path_to_bin = "/usr/local/bin/"
        self.m = m
        self.M = M
        self.N = N
        self.r = r
        self.pop_opts = pop_opts
        self.path_to_IQtree= "iqtree2"

    def make_sample_ini(self):
        name_list = []
        ini_list = os.listdir(self.indir)
        for name in ini_list:
            if name.endswith("R1_001.fastq.gz"):
                sample_name = name.split('_L001_R1_')[0]
                r1_file = name
                r2_file = sample_name + '_L001_R2_' + name.split('_L001_R1_')[1]
                name_list.append([sample_name, r1_file, r2_file])
        df = pd.DataFrame(name_list, columns=["#Sample_Name", "R1_File", "R2_File"])
        csv = "samples.ini"
        df.to_csv(csv, index=False)
        self.ini = csv
        print(f"csv file'{self.ini}'created.")

    def setup_logger(self, name=None):
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # creates a file handler that logs messages above DEBUG level
        fh = logging.FileHandler(self.log)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s %(funcName)s :\n%(message)s')
        fh.setFormatter(fh_formatter)
        # creates a file handler that logs messages above INFO level
        sh = logging.StreamHandler()
        sh.setLevel(logging.INFO)
        sh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s : %(message)s',
            '%Y-%m-%d %H:%M:%S')
        sh.setFormatter(sh_formatter)
        # add the handlers to logger
        logger.addHandler(fh)
        logger.addHandler(sh)
        self.logger = logger

    def execute_cmd(self, cmd, shell=True):
        self.logger.info('[cmd] {}'.format(cmd))
        proc = subprocess.run(
            cmd, shell=shell, stdout=PIPE,
            stderr=PIPE, text=True)
        stdout = proc.stdout
        stderr = proc.stderr
        if len(stdout) > 1:
            self.logger.info(stdout)
        if len(stderr) > 1:
            self.logger.info(stderr)
        return stdout, stderr

    def exec_cutadapt_help(self):
        cmd = 'cutadapt --help'
        self.execute_cmd(cmd)

    def load_ini(self):
        with open(self.ini, 'r') as ini:
            lines = ini.readlines()
        self.samples = {}
        for line in lines:
            if line.startswith('#'):
                continue
            data = [i.strip() for i in line.split(sep=',')]
            fq1 = os.path.join(self.indir, data[1])
            fq2 = os.path.join(self.indir, data[2])
            self.samples[data[0]] = [fq1, fq2]

    def fastq_quality_filter(self):
        fqfd = os.path.join(self.outdir, 'fastq_quality_filter')
        ftd = os.path.join(self.outdir, 'fastq_quality_filter')
        os.makedirs(fqfd, exist_ok=True)
        for sample in self.samples.keys():
            # forward read
            infq1 = self.samples[sample][0]
            outfq1 = os.path.join(fqfd, f'{sample}_R1.fastq.gz')
            if int(self.F_remove) == 0:
                cmd_fr = f'gunzip -c {infq1} | fastq_quality_filter -v -Q 33 -q {self.q} -p 40 -i - | gzip -c > {outfq1}'
            else:
                cmd_fr = f'gunzip -c {infq1} |  fastx_trimmer -Q 33 -f {self.F_remove} -i - |fastq_quality_filter -v -Q 33 -q {self.q} -p 40 -i - | gzip -c > {outfq1}'
            #execute one by one
            self.execute_cmd(cmd_fr)
            #reverse read
            infq2 = self.samples[sample][1]
            outfq2 = os.path.join(fqfd, f'{sample}_R2.fastq.gz')
            if int(self.R_remove) == 0:
                cmd_rr = f'gunzip -c {infq2} | fastq_quality_filter -v -Q 33 -q {self.q} -p 40 -i - | gzip -c > {outfq2}'
            else:
                cmd_rr = f'gunzip -c {infq2} |  fastx_trimmer -Q 33 -f {self.R_remove} -i - |fastq_quality_filter -v -Q 33 -q {self.q} -p 40 -i - | gzip -c > {outfq2}'
            self.execute_cmd(cmd_rr)
            
    def cutadapt(self):
        ca = os.path.join(self.outdir, 'cutadapt')
        os.makedirs(ca, exist_ok=True)
        for sample in self.samples.keys():
            # forward read
            infq1 = self.samples[sample][0]
            outfq1 = os.path.join(ca, f'{sample}_R1.fastq.gz')
            cmd_fr = f'cutadapt -e 0.05 -b {self.fada} -m {self.min_len} -o {outfq1} {infq1}'
            self.execute_cmd(cmd_fr)
            # reverse read
            infq2 = self.samples[sample][1]
            outfq2 = os.path.join(ca, f'{sample}_R2.fastq.gz')
            cmd_rr = f'cutadapt -e 0.05 -b {self.rada} -m {self.min_len} -o {outfq2} {infq2}'
            self.execute_cmd(cmd_rr)
            self.samples[sample] = [outfq1, outfq2]

    def re_pair_fastq(self):
        rp = os.path.join(self.outdir, 'repair/')
        os.makedirs(rp, exist_ok=True)
        for sample in self.samples.keys():
            # inf and inr from the dic
            inf = self.samples[sample][0]
            inr = self.samples[sample][1]
            # out1 and out2
            rpf = os.path.join(rp, f'{sample}_rp_R1.fastq')
            rpr = os.path.join(rp, f'{sample}_rp_R2.fastq')
            # repair; infiles are from the previously replaced dic
            cmd_rp = f"{os.path.join(self.path_to_bbmap,'repair.sh')} -Xmx50g in={inf} in2={inr} out1={rpf} out2={rpr} overwrite=true"
            # run in shell
            subprocess.run(cmd_rp, shell=True, check=True)
            # replace values in the dictionary
            self.samples[sample] = [rpf, rpr]
            
    def convert_fastq_to_fasta(self):
        fa = os.path.join(self.outdir, 'fa')
        os.makedirs(fa, exist_ok=True)
        print(f"Directory {fa} created")
        for sample in self.samples.keys():
            # inf and inr from the dic
            fastq_file_r1 = self.samples[sample][0]
            fastq_file_r2 = self.samples[sample][1]
            # output file
            fasta_file_r1 = os.path.join(fa, f'{sample}.1.fasta')
            fasta_file_r2 = os.path.join(fa, f'{sample}.2.fasta')
            # conversion
            SeqIO.convert(fastq_file_r1, "fastq", fasta_file_r1, "fasta")
            SeqIO.convert(fastq_file_r2, "fastq", fasta_file_r2, "fasta")
            self.samples[sample]=[fasta_file_r1, fasta_file_r2]
        self.inputfasta = [file for sublist in self.samples.values() for file in sublist]

    def rename(self):
        # Define directory for renamed files
        renamed = os.path.join(self.outdir, 'renamed')
        os.makedirs(renamed, exist_ok=True)

        for sample, file_paths in self.samples.items():
            for i, file_path in enumerate(file_paths):
            # Create new filename for the renamed file
                renamed_filename = f'{sample}.{i+1}.fasta.gz'
                renamed_file_path = os.path.join(renamed, renamed_filename)
            
            # Open the original fasta file
                with pysam.FastxFile(file_path) as fh_in, gzip.open(renamed_file_path, 'wt') as fh_out:
                    for entry in fh_in:
                        # Replace the name of each sequence
                        entry.name = f'{sample}:{entry.name}/{i+1}'
                        print('>' + entry.name, entry.sequence, sep='\n', file=fh_out)
                    
            # Update file paths in samples dict
            self.samples[sample][i] = renamed_file_path
            
    def pop_map_out(self, popmap=None):
        if popmap is not None:
            self.popmap = popmap
        else:
            names = list(self.samples.keys())
            names.sort()
            with open(os.path.join(self.outdir, "popmap.tsv"), "w") as f:
                for name in names:
                    f.write(f"{name}\t{name}\n")
            self.popmap = os.path.join(self.outdir, "popmap.tsv")
           
    def pl_stacks(self, from_fa=None):
        renamed = os.path.join(self.outdir, 'renamed')
        pl = os.path.join(self.outdir, 'pl')
        os.makedirs(pl, exist_ok=True)
        if from_fa is not None:
            renamed = from_fa
        cmd = f'denovo_map.pl -M 2 -T {self.t} -o {pl} --popmap {self.popmap} --samples {renamed} --paired -N 4 -X "ustacks:-m 3 --force-diff-len" -X "populations: -M {self.popmap} -R {self.r}  --max-obs-het 0.99 --write-single-snp --min-maf 0.01 --vcf --structure --plink --treemix --phylip -t {self.t} {self.pop_opts}"'
        self.execute_cmd(cmd)
        # ustacks options
        #-f  input file path.
        #-i  a unique integer ID for this sample.
        #-o  output path to write results.
        #-M  Maximum distance (in nucleotides) allowed between stacks (default 2).
        #-m  Minimum depth of coverage required to create a stack (default 3).
        #-N  Maximum distance allowed to align secondary reads to primary stacks (default: M + 2).
        #-p  enable parallel execution with num_threads threads.
        #-t  input file type. Supported types: fasta, fastq, gzfasta, or gzfastq (default: guess).
        #--name  a name for the sample (default: input file name minus the suffix).
        #-R  retain unused reads.
        #-H  disable calling haplotypes from secondary reads

    def md_phylip(self):
        # Create a md file and open it for writing
        infile_path = os.path.join(self.outdir,"pl", "populations.fixed.phylip") 
        outfile_path = os.path.join(self.outdir,"pl","pop.phy")
        with open(infile_path,"r") as infile, open(outfile_path, "w") as outfile:
            for line in infile:
                if not line.strip().startswith("#"):
                    outfile.write(line)

    def iqtree(self):
        phy_file = os.path.join(self.outdir, "pl", 'pop.phy')
        outpre = os.path.join(self.outdir, 'iq')
        cmd = f"{self.path_to_IQtree} -s {phy_file} -m MFP+ASC -bb 1000 --alrt 1000 -nt AUTO -pre {outpre}"
        self.execute_cmd(cmd)

    def run(self, popmap=None, from_fa=None):
        if from_fa is not None:
            print('From fa')
            if popmap is None:
                print('Please input popmap file')
                sys.exit(1)
            self.popmap = popmap
            self.pl_stacks(from_fa=from_fa)
            self.md_phylip()
            self.iqtree()
            return 0
        self.make_sample_ini()
        self.load_ini()
        self.fastq_quality_filter()
        self.cutadapt()
        self.re_pair_fastq()
        self.convert_fastq_to_fasta()
        self.rename()
        self.pop_map_out()
        self.pl_stacks()
        self.md_phylip()
        self.iqtree()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=str, default="rawdata_digitata", help="Input directory")
    # parser.add_argument("--ref_genome", type=str, default="./refgenome/RedCoralContigSpades.fasta", help="Reference genome")
    parser.add_argument("--q", type=int, default= 30, help="quality_value")
    parser.add_argument("--F_remove", type=int, default= 0, help="forward_remove_bp")
    parser.add_argument("--R_remove", type=int, default= 0, help="reverse_remove_bp")
    parser.add_argument("--outdir", type=str, default="outdir", help="Output directory")
    parser.add_argument("-t", type=int, default=4, help="Number of threads")
    parser.add_argument("--fada", type=str, default="GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", help="fada value")
    parser.add_argument("--rada", type=str, default="CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC", help="rada value")
    parser.add_argument("--min_len", type=int, default=80, help="Minimum length")
    parser.add_argument('--m', type=int, default=5, help='Minimum depth of coverage required to create a stack (default: 5)')
    parser.add_argument('--M', type=int, default=2, help='Maximum distance allowed between stacks (default: 2)')
    parser.add_argument('--N', type=int, default=1, help='Maximum distance allowed to align secondary reads to primary stacks (default: 1)')
    parser.add_argument("-r", type=float, default=0.1, help="r value")
    parser.add_argument("--popmap", type=str, default=None, help="Population map file")
    parser.add_argument("--from_fa", type=str, default=None, help="fasta file directory")
    parser.add_argument(
        "-pop_opts", '--populations_options', dest='popo',
        type=str,
        default='',
        help='populations command additional options')
    args = parser.parse_args()

    SCC = SampleCodeClass(
        indir=args.indir,
        q=args.q,
        F_remove=args.F_remove,
        R_remove=args.R_remove,
        outdir=args.outdir,
        t=args.t,
        fada=args.fada,
        rada=args.rada,
        min_len=args.min_len,
        m=args.m,
        M=args.M,
        N=args.N,
        r=args.r,
        pop_opts=args.popo
    )
    SCC.run(popmap=args.popmap, from_fa=args.from_fa)


if __name__ == '__main__':
    main()