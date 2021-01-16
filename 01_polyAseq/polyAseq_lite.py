#!/usr/bin/env python
'''
COPYRIGHT (C) 2017 changjin hong
author: hongc2@ccf.org
'''
import os, sys, argparse
import inspect
from lib.util.lib_utils import run_syscmd_wrapper, get_abspath_or_error, ensure_dir, check_if_file_valid, msgout
from lib.util.lib_par import par_shell_runs, find_error_index_from_retcodes
from lib.file_configs import ProgramConfig

PIPELINE_NAME = 'polyAseq'
VERSION = '0.1'
author_email = 'hongc2@ccf.org'

def split_barcodes_core(fastq_dir,
												prog_path,
												out_dir,
												log_prefix,
												barcode_fn,
												max_mismatch=2,
												ncpu = 2,
												debug = False):
	'''
	Note: this def can be generalized to dispatch multiple commands
	:param fastq_dir:
	:param prog_path:
	:param exp_dir:
	:param max_mismatch:
	:param ncpu:
	:return:
	'''

	fns = os.listdir(fastq_dir)
	arg_sets = []

	for fn in fns:
		msgout("notice","checking if [%s] is FASTQ file ..."%fn)
		if fn.endswith('.fastq.gz') or fn.endswith('.fastq'):

			cmd = prog_path
			cmd += " --outdir %s" % out_dir
			cmd += " --mm %d" % max_mismatch
			cmd += " --barcode_fn %s" % barcode_fn
			cmd += " --inpath %s"%os.path.join(fastq_dir,fn)
			stdout_fn = "%s_%s.out" % (log_prefix,fn)
			stderr_fn = "%s_%s.err" % (log_prefix,fn)
			msgout("notice",cmd)
			arg_sets.append((cmd,stdout_fn,stderr_fn,False,))

	if not debug: #debug
		retcodes = par_shell_runs(arg_sets,ncpu=ncpu,job_name=inspect.stack()[0][3])
		error_idx = find_error_index_from_retcodes(retcodes)
		if error_idx:
			for j in error_idx:
				print('an error in the command [%s]' % arg_sets[j][0])
			raise RuntimeError('there are some errors in split_barcdes_core()')

	return arg_sets

def bam_index(mod_bam_dir,debug=False):

	bams = os.listdir(mod_bam_dir)
	for bam in bams:
		if bam.endswith('.bam'):
			bai = os.path.join(mod_bam_dir,bam+'.bai')
			cmd = "samtools index %s"%(os.path.join(mod_bam_dir,bam))
			run_syscmd_wrapper(cmd,debug=debug)

class mStep:
	def __init__(self,name1,cnt1):
		self.name = name1
		self.cnt = cnt1
		self.substeps = []
		self.out_dir = None

	def get_out_dir(self,work_dir):
		if self.out_dir:
			return self.out_dir
		else:
			self.out_dir = ensure_dir(os.path.join(work_dir,'%02d_%s'%(self.cnt,self.name)))
			return self.out_dir

class cPolyAseq:
	def __init__(self, args):

		self.fastq_dir = get_abspath_or_error(args.fastq_dir)
		self.work_dir = ensure_dir(args.work_dir)
		self.comp_group_fn = check_if_file_valid(args.comp_group_fn,error_flag=True)
		self.ctrl_tag = args.ctrl_tag
		self.exp_tag = args.exp_tag
		self.barcode_fn = args.barcode_fn
		self.fq_index_fn = check_if_file_valid(
			os.path.join(os.path.dirname(self.comp_group_fn), 'fastq_index.csv'),error_flag=True)

		self.bam_dir = None
		self.gfa_bam_dir = None
		self.analy_dir = None
		self.prepropa_rd = None

		self.pconf = self.load_config(args)
		self.step = 0 #categorical step
		self.sstep = 0 #serial unique step
		self.mSteps = []

	def get_log_dir(self):
		log_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'log')
		if not os.path.exists(log_dir):
			os.makedirs(log_dir)
		return log_dir

	def launch_job(self):

		self.split_barcodes()
		self.make_filter_runs()
		self.prepropa()
		self.mostupdown()

	def print_msg(self,mS,sstep_name):
		msg = "[%02d-%s][%02d-%s]"%(mS.cnt,mS.name,self.sstep,sstep_name)
		print('%s in progress'%msg)

	def split_barcodes(self):
		self.step += 1 #1
		mS = mStep('ProcessSeq', self.step)
		self.mSteps.append(mS)
		sstep_name = 'split_barcodes'
		self.sstep += 1
		self.print_msg(mS,sstep_name)

		max_mismatch = 2
		out_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output', '%s_mm_%d' % (sstep_name,max_mismatch))
		ensure_dir(out_dir)

		log_prefix = os.path.join(self.get_log_dir(), '%02d_%s_mm'%(self.sstep,sstep_name))

		if self.sstep >= self.pconf.resume:

			_ = split_barcodes_core(fastq_dir = self.fastq_dir,
															prog_path = self.pconf.sections[sstep_name]['bin_path'],
															out_dir = out_dir,
															log_prefix=log_prefix,
															barcode_fn=self.barcode_fn,
															max_mismatch = max_mismatch,
															ncpu = self.pconf.ncpu,
															debug=self.pconf.debug)

		self.plot_barcodes(out_dir)

	def plot_barcodes(self,inpath):

		sstep_name = 'plot_barcodes'
		self.sstep += 1 #2
		self.print_msg(self.mSteps[-1], sstep_name)
		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']
			log_prefix = os.path.join(self.get_log_dir(),'%02d_%s'%(self.sstep,sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s"%inpath,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.a_run_hist(inpath,mm=2)

	def a_run_hist(self,demux_fastq_dir,mm=2):

		sstep_name = 'a_run_hist'
		self.sstep += 1 #3
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dirs = []
		prog_path = self.pconf.sections[sstep_name]['bin_path']

		outbase_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')
		ensure_dir(outbase_dir)
		suffix = '%02d_%s' % (self.sstep, sstep_name)

		log_prefix = os.path.join(self.get_log_dir(),suffix)
		for m in range(mm+1):
			suffix2 = "%s_mm%d"%(suffix,m)
			out_dir = os.path.join(outbase_dir,suffix2)
			ensure_dir(out_dir)

			cmd = "ls %s/*.fastq.gz | parallel -j %d %s --mm %d --outdir %s {}"%\
						(demux_fastq_dir,self.pconf.ncpu,prog_path,m,out_dir)

			if self.sstep >= self.pconf.resume:
				run_syscmd_wrapper(cmd,
													 stdout_fn="%s_%02d.out" % (log_prefix,m),
													 stderr_fn="%s_%02d.err" % (log_prefix,m),
													 debug = self.pconf.debug)

			out_dirs.append(out_dir)

		self.a_run_hist_plot(out_dirs[1],demux_fastq_dir)

	def a_run_hist_plot(self,inpath,demux_fastq_dir):
		sstep_name = 'a_run_hist_plot'
		self.sstep += 1 #4
		self.print_msg(self.mSteps[-1], sstep_name)
		outpath = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')
		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']
			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s'%(self.sstep,sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % inpath,
												 outopt="-o %s" % outpath,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.chop_a_tails(demux_fastq_dir)

	def chop_a_tails(self,demux_fastq_dir):
		sstep_name = 'chop_a_tails'
		self.sstep += 1 #5
		prog_path = self.pconf.sections[sstep_name]['bin_path']
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output', sstep_name)

		if self.sstep >= self.pconf.resume:
			if not os.path.exists(out_dir):
				os.makedirs(out_dir)

			cmd = "ls %s/*.fastq.gz | parallel -j %d %s --outdir %s {}" % \
						(demux_fastq_dir, self.pconf.ncpu, prog_path, out_dir)

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(cmd,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.a_chop_plot(out_dir)

	def a_chop_plot(self,inpath):
		sstep_name = 'a_chop_plot'
		self.sstep += 1 #6
		self.print_msg(self.mSteps[-1], sstep_name)
		outpath = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']
			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s'%(self.sstep,sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % inpath,
												 outopt="-o %s" % outpath,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.make_bowtie_runs(inpath)

	def make_bowtie_runs(self,inpath):
		sstep_name = 'make_bowtie_runs'
		self.sstep += 1 #7
		self.print_msg(self.mSteps[-1], sstep_name)
		self.bam_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'bowtie_alignments')
		ensure_dir(self.bam_dir)
		script_dir = os.path.join(self.bam_dir,'run_scripts')
		ensure_dir(script_dir)

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			runopt = "-f %s" % self.fq_index_fn

			runopt += " -n %d" % self.pconf.ncpu

			bowtie_index = self.pconf.sections['bowtie2-build']['ref_hg19_prefix']
			runopt += " -r %s" % bowtie_index

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % inpath,
												 outopt="-o %s" % script_dir,
												 runopt=runopt,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		script_fpaths = [os.path.join(script_dir,f) for f in os.listdir(script_dir) if f.endswith('.sh')]

		self.run_bowtie_ser(script_fpaths)

	def run_bowtie_ser(self,bowtie_scripts):
		sstep_name = 'run_bowtie_ser'
		self.sstep += 1 #8
		self.print_msg(self.mSteps[-1], sstep_name)
		if self.sstep >= self.pconf.resume:

			for c,bowtie_script in enumerate(bowtie_scripts):
				cmd = "bash %s" % bowtie_script

				log_prefix = os.path.join(self.get_log_dir(), '%02d_%s_%02d' % (self.sstep, sstep_name, c))

				run_syscmd_wrapper(cmd,
													 stdout_fn="%s.out" % log_prefix,
													 stderr_fn="%s.err" % log_prefix,
													 debug=self.pconf.debug)

	def make_filter_runs(self):
		self.step += 1
		mS = mStep('ProcessBams', self.step)
		self.mSteps.append(mS)

		sstep_name = 'make_filter_runs'
		self.sstep += 1 #9
		self.print_msg(self.mSteps[-1], sstep_name)
		mod_bam_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(),'%02d_%s' % (self.sstep, sstep_name))
			runopt = "-n %d"%self.pconf.ncpu
			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % self.bam_dir,
												 outopt="-o %s" % mod_bam_dir,
												 runopt=runopt,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

			sort_script_dir = os.path.join(mod_bam_dir, 'run_scripts_sort')
			filt_script_dir = os.path.join(mod_bam_dir, 'run_scripts_filter')

			so_script_fns = [os.path.join(sort_script_dir,f) for f in os.listdir(sort_script_dir) if f.endswith('.sh')]

			self.par_sort(so_script_fns)

			filt_script_fns = [os.path.join(filt_script_dir,f) for f in os.listdir(filt_script_dir) if f.endswith('.sh')]

			self.par_filter(filt_script_fns)

			bam_index(mod_bam_dir)

			self.rm_pc_dups(mod_bam_dir)
		else:
			self.gfa_bam_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output', 'rm_int_a_7of10or6')
			self.sstep += 7

	def par_sort(self, scripts):
		sstep_name = 'par_sort'
		self.sstep += 1 #10

		if self.sstep >= self.pconf.resume:
			for c, script in enumerate(scripts):
				log_prefix = os.path.join(self.get_log_dir(), '%02d_%s_%02d' % (self.sstep, sstep_name, c))
				run_syscmd_wrapper("bash %s" % script,
													 stdout_fn="%s.out" % log_prefix,
													 stderr_fn="%s.err" % log_prefix,
													 debug=self.pconf.debug)

	def par_filter(self, scripts):
		sstep_name = 'par_filter'
		self.sstep += 1 #11

		if self.sstep >= self.pconf.resume:
			for c, script in enumerate(scripts):
				log_prefix = os.path.join(self.get_log_dir(), '%02d_%s_%02d' % (self.sstep, sstep_name, c))
				run_syscmd_wrapper("bash %s" % script,
													 stdout_fn="%s.out" % log_prefix,
													 stderr_fn="%s.err" % log_prefix,
													 debug=self.pconf.debug)

	def rm_pc_dups(self,mod_bam_dir):
		sstep_name = 'rm_pc_dups'
		self.sstep += 1 #12
		prog_path = self.pconf.sections[sstep_name]['bin_path']
		subd = '%02d_%s' % (self.sstep, sstep_name)
		dedup_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), subd)
		_ = ensure_dir(dedup_dir)
		if self.sstep >= self.pconf.resume:

			cmd = "ls %s/*.02_mapq*.bam | parallel -j %d %s --log --outdir %s {}" % \
						(mod_bam_dir, self.pconf.ncpu, prog_path, dedup_dir)

			log_prefix = os.path.join(self.get_log_dir(),subd)

			run_syscmd_wrapper(cmd,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

			bam_index(dedup_dir)

		self.plot_dups(dedup_dir)
	
	def plot_dups(self,dedup_dir):
		#work from here
		sstep_name = 'plot_dups'
		self.sstep += 1 #13

		out_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')
	
		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']
	
			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))
			runopt = "-f %s"%self.fq_index_fn
			run_syscmd_wrapper(prog_path,
												 inopt = "-i %s" % dedup_dir,
												 runopt = runopt,
												 outopt = "-o %s"% out_dir,
												 stdout_fn = "%s.out" % log_prefix,
												 stderr_fn = "%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.plot_alignment_stats(dedup_dir)

	def plot_alignment_stats(self,dedup_dir):
		sstep_name = 'plot_alignment_stats'
		self.sstep += 1 #14
		mod_bam_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']
	
			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))
			runopt = "-n %d"%self.pconf.ncpu
			run_syscmd_wrapper(prog_path,
												 inopt = "-i %s" % dedup_dir,
												 outopt = "-o %s" % mod_bam_dir,
												 runopt = runopt,
												 stdout_fn = "%s.out" % log_prefix,
												 stderr_fn = "%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.filter_genomic_a(dedup_dir)

	def filter_genomic_a(self,dedup_dir):
		sstep_name = 'filter_genomic_a'
		self.sstep += 1 #15
		self.print_msg(self.mSteps[-1], sstep_name)
		self.gfa_bam_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output','rm_int_a_7of10or6')
		
		print(self.gfa_bam_dir)
		
		prog_path = self.pconf.sections[sstep_name]['bin_path']

		ref_genome_fa = self.pconf.sections['bowtie2-build']['ref_hg19_path']

		if self.sstep >= self.pconf.resume:
			if not os.path.exists(self.gfa_bam_dir):
				os.makedirs(self.gfa_bam_dir)

			cmd = "ls %s/*.bam | parallel -j%d %s --log --ref_genome_fa %s --outdir %s --win 10 --mm 3 --print --orfirst 6 {}" % \
						(dedup_dir, self.pconf.ncpu, prog_path, ref_genome_fa, self.gfa_bam_dir)

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(cmd,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

			bam_index(self.gfa_bam_dir)
			
		self.plot_check_genomic_a()

	def plot_check_genomic_a(self):
		sstep_name = 'plot_check_genomic_a'
		self.sstep += 1 #16
		out_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')
		self.print_msg(self.mSteps[-1], sstep_name)
		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % self.gfa_bam_dir,
												 outopt="-o %s" % out_dir,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)
			
			
	def prepropa(self):
		self.step += 1 #17
		mS = mStep('CallApa', self.step)
		self.mSteps.append(mS)

		sstep_name = 'prepropa'
		self.sstep += 1
		self.print_msg(self.mSteps[-1], sstep_name)
		work_dir = self.mSteps[-1].get_out_dir(self.work_dir)
		self.prepropa_rd = os.path.join(work_dir, 'output', 'prepropa.rd')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			track_dir = os.path.join(work_dir, 'output', 'tracks')
			_ = ensure_dir(track_dir)

			runopt = "-t %s" % track_dir
			runopt += " -B %s" % self.gfa_bam_dir
			runopt += " -f %s" % self.pconf.sections['ftp_server_url']['value']
			runopt += " -b %s" % self.pconf.sections['ensembl']['biotype_table']
			runopt += " -g %s" % self.pconf.sections['ensembl']['gtf_file']
			runopt += " -T %s" % self.pconf.sections['cache_dir']['file_path']
			runopt += " -n %d" % self.pconf.ncpu

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % self.fq_index_fn,
												 outopt="-o %s" % self.prepropa_rd,
												 runopt=runopt,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.testbatchadj()

	
	def testbatchadj(self):

		sstep_name = 'testbatchadj'
		self.sstep += 1 #18
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir = self.mSteps[-1].get_out_dir(self.work_dir)
		out_fn = os.path.join(out_dir, 'output', 'apa.rd')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			runopt = "-n %d" % self.pconf.ncpu
			runopt += " -c %s" % self.comp_group_fn

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % self.prepropa_rd,
												 runopt=runopt,
												 outopt="-o %s" % out_fn,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.alttests(out_fn)

	def alttests(self,apa_rd):
		sstep_name = 'alttests'
		self.sstep += 1 #19
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir = self.mSteps[-1].get_out_dir(self.work_dir)
		out_fn = os.path.join(out_dir, 'output', 'apa.alt.rd')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % apa_rd,
												 runopt="-w %s" % self.prepropa_rd,
												 outopt="-o %s" % out_fn,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.callsig(out_fn)

	def callsig(self,apa_alt_rd):
		sstep_name = 'callsig'
		self.sstep += 1 #20
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir = self.mSteps[-1].get_out_dir(self.work_dir)

		out_fn = os.path.join(out_dir, 'output', 'apa.sig.rd')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % apa_alt_rd,
												 runopt="-w %s" % self.prepropa_rd,
												 outopt="-o %s" % out_fn,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix)

		self.reporting(out_fn)

	def reporting(self,apa_sig_rd):
		sstep_name = 'reporting'
		self.sstep += 1 #21
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir= ensure_dir(os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output'))

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			runopt = "-w %s" % self.prepropa_rd
			runopt +=" -s %s" % self.pconf.sections['bowtie2-build']['chrom_size']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % apa_sig_rd,
												 outopt="-o %s" % out_dir,
												 runopt=runopt,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)

		self.intersect_gainedlost(apa_sig_rd)

	def intersect_gainedlost(self,apa_sig_rd):
		sstep_name = 'intersect_gainedlost'
		self.sstep += 1 #22
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir= ensure_dir(os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output'))

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			# runopt = "-w %s" % self.prepropa_rd
			runopt = " -c %s" % self.ctrl_tag
			runopt += " -e %s" % self.exp_tag

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % apa_sig_rd,
												 outopt="-o %s" % out_dir,
												 runopt=runopt,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)
		
		self.plot_hexamer()

	def plot_hexamer(self):
		sstep_name = 'plot_hexamer'
		self.sstep += 1 #23
		out_dir = os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output')

		self.print_msg(self.mSteps[-1], sstep_name)
		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s -d %s" % (self.prepropa_rd,self.pconf.sections['polya_db2']['bed_path']),
												 outopt="-o %s" % out_dir,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)
			
		self.prepropa_to_excel()
	
	def prepropa_to_excel(self):
		sstep_name = 'prepropa_to_excel'
		self.sstep += 1 #24
		in_dir = self.mSteps[-1].get_out_dir(self.work_dir)

		self.print_msg(self.mSteps[-1], sstep_name)
		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % in_dir,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)
		
		self.rep_cor()
	
	def rep_cor(self):
		sstep_name = 'rep_cor'
		self.sstep += 1 #25
		in_dir = self.mSteps[-1].get_out_dir(self.work_dir)
		apa_sig_rd = os.path.join(in_dir, 'output', 'apa.sig.rd')

		self.print_msg(self.mSteps[-1], sstep_name)
		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i1 %s -i2 %s" % (self.prepropa_rd,apa_sig_rd),
												 outopt="-o %s" % os.path.join(in_dir, 'output'),
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)


	def mostupdown(self):

		prev_out_dir = self.mSteps[-1].get_out_dir(self.work_dir)
		apa_sig_rd = os.path.join(prev_out_dir, 'output', 'apa.sig.rd')

		self.step += 1 #26
		mS = mStep('AnnoApa', self.step)
		self.mSteps.append(mS)
		sstep_name = 'mostupdown'
		self.sstep += 1
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir= ensure_dir(os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output'))
		out_fn = os.path.join(out_dir, 'apa.ann.rd')

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % apa_sig_rd,
												 outopt="-o %s" % out_fn,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)
		
		_ = self.mostupdown_all_samples(apa_sig_rd)
	
	
	def mostupdown_all_samples(self,apa_sig_rd):

		self.step += 1 #27
		mS = mStep('AnnoApa', self.step)
		self.mSteps.append(mS)
		sstep_name = 'mostupdown_all_samples'
		self.sstep += 1
		self.print_msg(self.mSteps[-1], sstep_name)
		out_dir= ensure_dir(os.path.join(self.mSteps[-1].get_out_dir(self.work_dir), 'output'))

		if self.sstep >= self.pconf.resume:
			prog_path = self.pconf.sections[sstep_name]['bin_path']

			log_prefix = os.path.join(self.get_log_dir(), '%02d_%s' % (self.sstep, sstep_name))

			run_syscmd_wrapper(prog_path,
												 inopt="-i %s" % apa_sig_rd,
												 outopt="-o %s" % out_dir,
												 stdout_fn="%s.out" % log_prefix,
												 stderr_fn="%s.err" % log_prefix,
												 debug=self.pconf.debug)
		
		return out_dir
		
	def load_config(self,args):
		"""
		objective: to load a program config file that contains subprogram paths, tools, run options, and databases
		:return:
		"""
		cConf = ProgramConfig(config_file=args.config_file,
													nnode=args.nnode,
													ncpu=args.ncpu,
													resume=args.resume,
													debug=args.debug)
		cConf.load_items()

		return cConf

def main():
	parser = argparse.ArgumentParser(description="%s_v%s [author:%s]" % (PIPELINE_NAME,VERSION, author_email))
	parser.add_argument('-C', '--config_file', dest='config_file', action='store', required=False, default=None,
											help='a program configuration file [configs/program.conf]')
	parser.add_argument('-i', '--fastq_dir', dest='fastq_dir', action='store', required=True, help='fastq directory.')

	parser.add_argument('-b', '--barcode_fn', dest='barcode_fn', action='store', required=True, help='barcode_fn')

	parser.add_argument('-e', '--exp_tag', dest='exp_tag', action='store', required=True, help='experimental group tag')
	parser.add_argument('-c', '--ctrl_tag', dest='ctrl_tag', action='store', required=True, help='control group tag')
	parser.add_argument('-d', '--comp_group_fn', dest='comp_group_fn', action='store', required=True,
											help='comp_group_fn')

	parser.add_argument('-n', '--nnode', action='store', dest='nnode', required=False, default=1, type=int,
											help='the number of nodes to use in MPI mode [1];under development')
	parser.add_argument('-p', '--ncpu', action='store', dest='ncpu', required=False, default=4, type=int,
											help='the number of cpus to utilize [4]')
	parser.add_argument('-o', '--work_dir', action='store', dest='work_dir', required=True,
											help='output directory without white space.')
	parser.add_argument('-r', '--resume', action='store', dest='resume', required=False, default=-1, type=int,
											help='specify the step number you want to resume if the successful result is available [-1]')
	parser.add_argument('-s', '--sample_sheet', dest='sample_sheet', required=False, default=None, action='store',
											help='Sample sheet [None];this feature is under development.')
	parser.add_argument('--debug', action='store_const', dest='debug', required=False, default=False, const=True,
											help='print commandline only in debug mode [False]')

	args = parser.parse_args()

	cP = cPolyAseq(args)
	cP.launch_job()

if __name__ == '__main__':
	main()
