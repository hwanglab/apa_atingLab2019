#!/usr/bin/env python
# NOTE: separate anything related FASTA or FASTQ file and put them into lib_fastq.py

import os, sys, re, inspect, shutil, zlib, gzip, random
import subprocess as sp
import numpy as np
from datetime import datetime

def open2(fn,open_mode):
	fp = None
	if open_mode == 'r':
		if not check_if_file_valid(fn):
			return fp
		if fn.endswith('.gz'):
			try:
				fp=gzip.open(fn,'r')
			except IOError:
				try:
					fp=gzip.open(fn,'rb')
				except IOError:
					raise IOError('cannot open the file[%s] to read'%fn)
		else:
			fp=open(fn,'r')
	elif open_mode == 'w':
		if fn.endswith('.gz'):
			try:
				fp=gzip.open(fn,'w')
			except IOError:
				try:
					fp=gzip.open(fn,'wb')
				except IOError:
					raise IOError('cannot open the file[%s] to write'%fn)
		else:
			fp=open(fn,'w')
	elif open_mode == 'a':
		if fn.endswith('.gz'):
			try:
				fp=gzip.open(file,'a')
			except IOError:
				try:
					fp=gzip.open(file,'ab')
				except IOError:
					raise IOError('cannot open the file[%s] to append'%file)
		else:
			fp=open(fn,'a')
	else:
		raise IOError('check file[%s] and open_mode[%s]'%(fn,open_mode))
	return fp

def check_if_file_valid(fpath,error_flag=True):
	ret_fpath=None
	if os.path.exists(fpath):
		if os.path.isdir(fpath):
			ret_fpath = fpath
		elif os.stat(fpath).st_size>0:
			ret_fpath = fpath

	if error_flag and not ret_fpath:
		raise IOError('The file [%s] does not exist!'%fpath)
	return ret_fpath

def count_num_lines(file2):
	
	if file2.endswith('gz'):
		cmd='zcat %s | wc -l | sed -e \'s/^[ \\t]*//\' | awk -F\" *\" \'{print $1}\'' % file2
	else:
		cmd='wc -l %s | sed -e \'s/^[ \\t]*//\' | awk -F\" *\" \'{print $1}\'' % file2
		
	#print cmd #debug
	
	proc = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
	(msg, err) = proc.communicate()
	if msg:
		try: num_lines=int(msg.strip())
		except: num_lines = 0
	else:
		num_lines=0
		print 'num_lines=0'
	return num_lines

def intersect(a, b):
	# return the intersection of two lists
	return list(set(a) & set(b))

def union(a, b):
	# the union of two lists
	return list(set(a) | set(b))

def difference(a, b):
	# show whats in list b which isn't in list a
	return list(set(b).difference(set(a)))

def joined(things, delimiter):
	return delimiter.join(map(str, things))

def unique(a):
	# return the list with duplicate elements removed
	return list(set(a))

def ensure_dir(d):
	ret_d = None
	if os.path.exists(d):
		ret_d = d
	else:
		try:
			os.makedirs(d)
			ret_d = d
		except IOError:
			raise IOError('cannot create the directory [%s]'%d)
	return ret_d

def get_abspath_or_error(fpath):
	try :
		fpath = os.path.abspath(fpath)
		if os.path.exists(fpath):
			return fpath
		else:
			raise IOError('the file path [%s] does not exist ' %fpath)
	except IOError:
		raise IOError('the file path [%s] does not exist' % fpath)

def file2list(fn,upper=False):
	
	if not os.path.exists(fn):
		raise IOError('check input arguments %s' % fn)
		
	myList=[]
	fp = open(fn,'r')
	for i in fp:
		if i.startswith('#'):continue
		item=i.strip()
		if upper:
			item=item.upper()
		myList.append(item)
	fp.close()
	return myList

def file_tag2(filename,tag,newExt):
	D,_,fBase,fExt = separateDirFn2(filename)
	if not newExt:
		newExt=fExt
	if not tag:
		taggedFn='%s/%s' % (D,fBase)
	else:
		taggedFn='%s/%s_%s' % (D,fBase,tag)

	if newExt:
		taggedFn='%s.%s' % (taggedFn,newExt)
	if filename == taggedFn:
		msgout('warning','tag file name[%s] is same as the original input [%s]' % (taggedFn,filename))
	return taggedFn

def get_stat_dic(myDic,mode):
	
	if mode == 'min':
		myStat = 1e7
		for val in myDic.itervalues():
			if val<myStat:
				myStat = val
	return myStat

def list_to_dic(myList1,initVal):
	'''
	Return a dictionary containing key from the list.  The dictionary values will be set to initVal
	'''
	myDic = {}
	if myList1:
		for i in myList1:
			myDic[i]=initVal
	return myDic


def merge_multiple_gz_in_order(fileList,sortList,bigCatGz):
	
	if sortList:
		fileList = fileList.sort()
	
	#check if zcat is available (TODO)
	inGzStr=''
	for f in fileList:
		#check if f is gz file
		if not f.endswith('gz'):
			raise IOError('the input file[%s] should be gzipped'%f)
		inGzStr += ' %s' % f

	if bigCatGz.endswith('gz'):
		cmd = 'zcat %s | gzip -cf > %s' % (inGzStr,bigCatGz)
	else:
		cmd = 'zcat %s > %s' % (inGzStr,bigCatGz)
	runcmd(cmd,'merge_multiple_gz_in_order')
	
def runcmd(cmd,call_from=inspect.stack()[1][3],\
           stdout2=sp.PIPE,stderr2=sp.PIPE):

		subproc = sp.Popen(cmd, stdout=stdout2, stderr=stderr2, shell=True)
		output,error = subproc.communicate()
		exit_code = subproc.wait()

		if exit_code>0:
			raise RuntimeError('Error[%s] occurs at running [%s]; call from [%s]'%\
				(error,cmd,call_from))

def run_syscmd_old(cmd_entries, stdout_fn=None, stderr_fn=None, handle_err=False, pipe=None):

	if pipe:
		cmd_entries.append('| %s'%pipe)

	if stdout_fn:
		stdout_fp = open(stdout_fn,'w')
	else:
		stdout_fp = sp.PIPE

	if stderr_fn:
		stderr_fp = open(stderr_fn,'w')
	else:
		stderr_fp = sp.PIPE

	err_flag = False
	try:
		proc = sp.Popen(cmd_entries, stdout=stdout_fp, stderr=stderr_fp, shell=True, close_fds=True)
		output,retcode = proc.communicate()
	except RuntimeError:
		err_flag = True

	if stdout_fp: stdout_fp.close()
	if stderr_fp: stderr_fp.close()

	if err_flag > 0:
		print 'retcode:', retcode
		if handle_err:
			return retcode, cmd_str
		else:
			raise RuntimeError('[%s] failed' % cmd_str)
	else:
		return retcode, cmd_str

def run_syscmd(cmd_str, stdout_fn=None, stderr_fn=None, handle_err=False):
	'''
	  p = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    return(stdout, stderr, p.returncode)
	'''
	print(cmd_str)
	# return (0,cmd_str)

	if stdout_fn:
		stdout_fp = open(stdout_fn, 'w')
	else:
		stdout_fp = sp.PIPE

	if stderr_fn:
		stderr_fp = open(stderr_fn, 'w')
	else:
		stderr_fp = sp.PIPE

	print('running [%s]'%cmd_str)
	err_flag = False
	try:
		proc = sp.Popen(cmd_str, stdout=stdout_fp, stderr=stderr_fp, shell=True, close_fds=True)
		_, _ = proc.communicate()
	except RuntimeError:
		err_flag = True

	if stdout_fn: stdout_fp.close()
	if stderr_fn: stderr_fp.close()

	retcode = int(proc.returncode)
	if err_flag or retcode != 0:
		print 'retcode:', retcode
		if handle_err:
			return retcode, cmd_str
		else:
			raise RuntimeError('[%s] failed' % cmd_str)
	else:
		return retcode, cmd_str

def get_process_msg_handler(log_dir, step_name):
	
	stdout_fn = os.path.join(log_dir, '%s.out'%step_name)
	stdofp = open(stdout_fn, 'w')
	stderr_fn = os.path.join(log_dir, '%s.err'%step_name)
	stdefp = open(stderr_fn, 'w')
	return stdofp, stdefp
	
def msgout(mode,msg,call_from=inspect.stack()[1][3]):
	if mode in ['out','notice']:
		print '%s [INFO:%s] %s\n'%(str(datetime.now()), call_from, msg)
	elif mode == 'warning':
		print '%s [WARNING:%s] %s'%(str(datetime.now()), call_from, msg)
	elif mode == 'banner':
		print '############################################'
		print '%s [PART:%s] %s\n'%(str(datetime.now()), call_from, msg)
		print '############################################'

def normalize_dic(myDic,mode='sum'):
	
	if mode == 'max':
		denom = -1e7
		for val in myDic.itervalues():
			if val>denom:
				denom = val
	elif mode == 'sum':
		denom = 0.
		for val in myDic.itervalues():
			denom += val
	else:
		print 'mode[%s] is not registered'%mode
		raise ValueError
	for key in myDic.keys():
		myDic[key] /= denom
		
	return myDic


def separateDirFn2(fullPath): #fullPath can file w|w/o extension or directory
	
	if os.path.exists(fullPath):
		myDir=os.path.dirname(os.path.abspath(fullPath))
		fname=os.path.basename(fullPath).strip()
	else:
		loc=fullPath.rfind(os.path.sep)
		myDir=fullPath[:loc]
		fname=fullPath[loc+1:]
	
	if os.path.isdir(fullPath):
		fileBase = fname
		fileExt = ''
	else: #handle both actual file or preserved(not existing) filename
		#orgPath = myDir+'/'+fname
		gzExt = False
		if fname.endswith('.gz'):
			fname=fname[:-3]
			gzExt = True
		loc=fname.rfind('.')
		if loc:
			fileBase=fname[:loc]
			fileExt=fname[loc+1:]
		else:
			fileBase=fname
			fileExt=''
		if gzExt:
			fileExt+='.gz'

	return (myDir,fname,fileBase,fileExt)


def sort_tsv_by_col2(fn,Col2sort,Modes,uniqF,outfile,temp_dir=None):
	#check if file exists
	if not check_if_file_valid(fn):
		raise RuntimeError('error',fn+' is not valid','lib_utils')
	
	files_to_delete = []
	if fn.endswith('.gz'):
		file2 = fn[:-3]
		archive_file(fn,file2,'gunzip')
		files_to_delete.append(file2)
	else:
		file2 = fn

	outfile2 = outfile
	if outfile.endswith('.gz'):
		outfile2 = outfile[:-3]
	
	#count the number of header
	header2skip = count_heads(file2,'#')
	
	if header2skip>0:
		cmd1 = 'head -n%d %s > %s' % (header2skip,file2,outfile2)
		runcmd(cmd1,'sort_tsv_by_col2')
	
	if header2skip>0:
		cmd1 = 'tail -n +%d %s' % (header2skip+1,file2)
	else:
		cmd1 = 'cat %s' % file2
		
	C = len(Col2sort)
	argSort='sort'
	if temp_dir is not None:
		argSort += ' -T %s'%temp_dir
	for i in range(C):
		argSort += ' -k%d,%d%s' % (Col2sort[i],Col2sort[i],Modes[i])
	cmd1 += ' | %s >> %s' % (argSort,outfile2)
	runcmd(cmd1,'sort_tsv_by_col2')
	
	if outfile.endswith('.gz'):
		archive_file(outfile2,outfile,'gzip')
	
	unlink_fns(files_to_delete)
	
def archive_file(file2,out,fileOp):
	if fileOp=='delete':
		os.remove(file2)
	elif fileOp=='gzip':
		runcmd('gzip -f '+file2)
	elif fileOp=='gunzip':
		if out:
			runcmd('gunzip %s -fc > %s' % (file2,out))
		else:
			runcmd('gunzip '+file2)
	elif fileOp=='copy':
		shutil.copyfile(file2,out)

def unlink_fns(files):
	for fn in files:
		if os.path.exists(fn):
			os.unlink(fn)

def count_heads(fn,start_with_this):
	fp=open2(fn,'r')
	H=0
	in_head_section = False
	for i in fp:
		if in_head_section:
			if i.startswith('%s'%start_with_this): H+=1
			else: break
		elif i.startswith('%s'%start_with_this):
			in_head_section = True
			H+=1
	fp.close()
	return H

def get_work_directory(prefix='divine'):
  new_created = False
  attempt = 0
  max_attempt = 10
  while (not new_created) and (attempt<max_attempt):
    workD = '/tmp/tmp_%s_%d'%(prefix,random.randint(1,1e7))
    attempt += 1 
    if not os.path.exists(workD):
      os.makedirs(workD)
      new_created = True
  return workD

def file_tag(filename,tag,newExt):
	D,_,fBase,fExt = separateDirFn2(filename)
	if not newExt:
		newExt=fExt

	if tag:
		taggedFn='%s/%s_%s' % (D,fBase,tag)
	else:
		taggedFn='%s/%s' % (D,fBase)

	if newExt:
		taggedFn='%s.%s' % (taggedFn,newExt)
		
	if filename == taggedFn:
		msgout('warning','tag file name[%s] is same as the original input [%s]' %\
								(taggedFn,filename),'file_tag2')
	return taggedFn

def segments_intersect(min1, max1, min2, max2):
	return max(0, min(max1, max2) - max(min1, min2) + 1)

def month_to_num(month):
	months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
	month_in_num = months.index(month.lower())
	if month_in_num:
		return '%02d' % (month_in_num+1)
	else:
		return None

def get_file_base(file_path,file_ext):
	file_base = None
	if file_path.endswith(file_ext):
		file_base = file_path.split(file_ext)[0]

	print 'file_base:',file_base
	return file_base

class py_struct(object):
	def __init__(self, **kwargs):

		self.set_attr(**kwargs)

		if 'source' not in self.__dict__.keys():
			self.source = 'py_struct'

	def tup2dic(self, dbrec):
		'''
		take tuple records from db and store into class member variables (list)
		:param dbrec:
		:return:
		'''
		for key in dbrec._fields:
			if key not in self.__dict__:
				self.__dict__[key] = [] # variables, indicator whether it prints to VCF info
			self.__dict__[key].append(getattr(dbrec.key))

	# def set_attr(self, **kwargs):
	# 	self.__dict__.update(kwargs)

	def get_attr(self, name):
		try:
			return self[name]
		except KeyError:
			raise AttributeError(name)

	def set_attr(self, **kwargs):
		'''
		:param kwargs:
		:return:
		'''
		for key, value in kwargs.iteritems():

			if key in self.__dict__.keys():
				if not value or (type(self.__dict__[key]) != list):
					self.__dict__[key] = value
				else:
					if type(value) is list:
						self.__dict__[key].extend(value)
					else:
						self.__dict__[key].append(value)
			else:
				self.__dict__[key] = value

	def dedup_attr(self):
		for key, value in self.__dict__.iteritems():
			if type(value) is list:
				self.__dict__[key] = list(set(value))

	def rm_attr(self):
		for key in self.__dict__.keys():
			self.__dict__.pop(key, None)

	def reset_attr(self,excl_keys=['name']):
		for key, value in self.__dict__.iteritems():
			if key not in excl_keys:
				value_type = type(value)
				if value_type is list:
					self.__dict__[key] = []
				elif value_type is bool:
					self.__dict__[key] = False
				elif value_type is dict:
					self.__dict__[key] = {}
				else:
					self.__dict__[key] = None

def add_new_column_to_dataframe(df, new_cols):
	new_cols_upd = []
	for new_col in new_cols:
		if new_col not in df.columns:
			new_cols_upd.append(new_col)

	if new_cols_upd:
		df = df.reindex(columns=np.append(df.columns.values, new_cols_upd))
	return df

def cmd_exists(cmd):
	return any(
			os.access(os.path.join(path, cmd), os.X_OK)
			for path in os.environ["PATH"].split(os.pathsep)
	)

def count_existing_paths(fpaths):
	fcnt = 0
	if fpaths:
		for fpath in fpaths:
			if os.path.isdir(fpath):
				if os.listdir(fpath):
					fcnt+=1
			elif os.path.exists(fpath):
				fcnt+=1
			else:
				print('the file path[%s] does not exist'%fpath)
	else:
		print('input lists (fpaths) are empty')

	return fcnt

def run_syscmd_wrapper(prog_cmd,
											 inopt=None,
											 outopt=None,
											 runopt=None,
											 stdout_fn=None,
											 stderr_fn=None,
											 debug=False):

	cmd = prog_cmd
	if inopt:
		cmd += " %s" % inopt
	if outopt:
		cmd += " %s" % outopt
	if runopt:
		cmd += " %s" % runopt

	msgout("notice", cmd)
	if not debug:
		run_syscmd(cmd, stdout_fn=stdout_fn, stderr_fn=stderr_fn)

def parse_file_name(fpath):

	fbase = os.path.basename(fpath)
	fdir = os.path.dirname(fpath)
	if fbase.endswith('.gz') or fbase.endswith('bz2'):
		fitems = fbase.split('.')
		if len(fitems) > 2:
			fname = '.'.join(fitems[:-2])
			fprefix = os.path.join(fdir, fname)
			fext = '.'.join(fitems[-2:])
		else:
			fprefix, fext = os.path.splitext(fpath)
	else:
		fprefix, fext = os.path.splitext(fpath)

	fname = fprefix.split(os.sep)[-1]
	file = py_struct(fpath=fpath,
									 fdir=fdir,
									 sample_name=fname,
									 fprefix=fprefix,
									 fext=fext)

	return file
