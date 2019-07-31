import subprocess
import re
import os
import tempfile

#overall structure
# Type F with UMI:
# NNNNNN (umi) cttgTTG (sgRNA) GTTTAAGAGCTAAGCTG
# Type F without UMI:
# (sgRNA) GTTTAAGAGCTAAGCTGGAAACAGCATAGC
# Type R without UMI:
# GCTCTTAAAC revcomp(sgRNA) CAACAAGGTGGTTCTCCAAG
# Type R with UMI:
# GCTCTTAAAC revcomp(sgRNA) CAACAAG NNNNNN (umi) GTGGTTC

class process:
  def __init__(self, metadata, iAnalyzeRdir, reffile, log, hasUMI=True, aligner="BWA", runMode=True):
    self.metadata = metadata
    self.iAnalyzeRdir = iAnalyzeRdir
    self.reffile = reffile
    self.log = log
    self.umi = hasUMI
    self.aligner = aligner
    self.rm = runMode #print mode if set to false                                                                                       

  def collpase (self):
    if (not self.umi): return True
    
    meta = self.metadata
    log = self.log
    log.write("####Collapsing identical reads#####\n")

    for sample in meta.index:
      print(sample)
      fq = meta.loc[sample, 'fastq']
      tmpdir = tempfile.mkdtemp()
      cmd = "gunzip -c %s | awk \'{ if(NR%%4==1) {print $1} } \' > %s/id" % (fq, tmpdir)
      print(cmd)
      log.write("%s\n" % (cmd))
      if self.rm: subprocess.call(cmd,shell=True)
      cmd = "gunzip -c %s | awk \'{ if(NR%%4==2) {print $1} } \' > %s/seq" % (fq, tmpdir)
      print(cmd)
      log.write("%s\n" % (cmd))
      if self.rm: subprocess.call(cmd,shell=True)
      cmd = "gunzip -c %s | awk \'{ if(NR%%4==0) {print $1} } \' > %s/qual" % (fq, tmpdir)
      print(cmd)
      log.write("%s\n" % (cmd))
      if self.rm: subprocess.call(cmd,shell=True)

      out = re.sub(".fastq.gz", ".c.fastq.gz", fq)
      meta.loc[sample, 'fastq'] = out
      cmd = "paste %s/id %s/qual %s/seq | sort -k 3 | uniq -f 2 -c | awk '{print $2 \"#\" $1 \"\\n\" $4 \"\\n+\\n\" $3 }' | gzip -c > %s"  %(tmpdir,tmpdir,tmpdir,out)
      print(cmd)
      log.write("%s\n" % (cmd))
      if self.rm: subprocess.call(cmd,shell=True)

  def umi_extract (self):
    if (not self.umi): return True

    meta = self.metadata
    log = self.log
    log.write("####Extracting UMIs#####\n")
    for sample in meta.index:
      print(sample)
      fq = meta.loc[sample, 'fastq']
      out = re.sub(".fastq.gz", ".bc.fastq.gz", fq)
      
      meta.loc[sample, 'fastq'] = out
      lib_type = meta.loc[sample, 'lib.type']
      if (lib_type=="F"):
        cmd = 'umi_tools extract --stdin={} --bc-pattern=NNNNNNXXXXXXX -L log.out --stdout={}'.format(fq, out)
        print(cmd)
        log.write("%s\n" % (cmd))
        if self.rm: subprocess.call(cmd,shell=True)
      elif (lib_type=="R"):
        cmd = 'umi_tools extract --stdin={} --3prime --bc-pattern=NNNNNNXXXXXXX -L log.out --stdout={}'.format(fq, out)
        print(cmd)
        log.write("%s\n" % (cmd))
        if self.rm: subprocess.call(cmd,shell=True)
    log.write("\n\n")

  def trim(self):
    meta = self.metadata
    log = self.log
    log.write("####Removing adaptors#####\n")

    for sample in meta.index:
      print(sample)
      fq = meta.loc[sample, 'fastq']
      out = re.sub(".fastq.gz", ".trim.fastq.gz", fq)
      meta.loc[sample, 'fastq'] = out
      lib_type = meta.loc[sample, 'lib.type']
      if (lib_type=="F"):
        cmd = 'cutadapt -j 8 -m 19 -M 21 -q 15 -a "^CTTGTTG;optional...GTTTAAGAGCTAAGCTG" -o {} {}'.format(out,fq)
        print(cmd)
        log.write("%s\n" % (cmd))
        if self.rm: subprocess.call(cmd,shell=True)
      elif (lib_type=="R"):
        cmd = 'cutadapt -j 8 -m 19 -M 21 -q 15 -a ^GCTCTTAAAC...CAACAAGGTGGTTCTCCAAG -o {} {}'.format(out,fq)
        print(cmd)
        log.write("%s\n" % (cmd))
        if self.rm: subprocess.call(cmd,shell=True)
  
  def align(self):
    meta = self.metadata
    log = self.log
    for sample in meta.index:
      print(sample)
      fq = meta.loc[sample, 'fastq']
      lib_type = meta.loc[sample, 'lib.type']
      if (self.aligner=="bowtie2"):
        cmd = 'bowtie2 --sensitive --end-to-end -N 1 -p 6 -x {} -U {} | samtools view - -Sb  -h -t {}.fa.fai -o {}.bam'.format(self.reffile,fq,self.reffile,sample)
      if (lib_type=="F"):
        cmd = cmd + '\nsamtools view -F 0x14 -q 10 -Sb {}.bam > {}.f.bam; samtools sort -@ 6 -o {}.srt.bam {}.f.bam; samtools index {}.srt.bam'.format(sample,sample,sample,sample,sample)
      elif (lib_type=="R"):
        cmd = cmd + '\nsamtools view -f 0x10 -q 10 -Sb {}.bam > {}.f.bam; samtools sort -@ 6 -o {}.srt.bam {}.f.bam; samtools index {}.srt.bam'.format(sample,sample,sample,sample,sample)
      print(cmd)
      log.write("%s\n" % (cmd))
      if self.rm: subprocess.call(cmd,shell=True)
    log.write("\n\n")

  def dedup(self):
    if (not self.umi): return True

    meta = self.metadata
    log = self.log
    for sample in meta.index:
      cmd = 'umi_tools dedup -I {}.srt.bam --output-stats=deduplicated -S {}.dd.bam'.format(sample,sample)
      print(cmd)
      log.write("%s\n" % (cmd))
      if self.rm: subprocess.call(cmd,shell=True)
    log.write("\n\n")

  def count(self):
    meta = self.metadata
    log = self.log
    for sample in meta.index:
      if (self.umi):
        cmd = 'samtools view %s.dd.bam | cut -f3 | sort | uniq -c | awk \'{ print $2 "\t" $1}\' > %s.cnt' % (sample,sample)
      else:
        cmd = 'samtools view %s.srt.bam | cut -f3 | sort | uniq -c | awk \'{ print $2 "\t" $1}\' > %s.cnt' % (sample,sample)
      print(cmd)
      log.write("%s\n" % (cmd))
      if self.rm: subprocess.call(cmd,shell=True)
    log.write("\n\n")

  def univariate(self, metadata, design, ref, weighting, outfile):
    log = self.log
    log.write("####Performing univariate analysis using DESeq2#####\n")
    cmd = 'Rscript {}/univariate_analysis.R {} {} {} {} {}'.format(self.iAnalyzeRdir, metadata, design, ref, weighting, outfile)
    print(cmd)
    subprocess.call(cmd,shell=True)
    log.write("\n\n")
