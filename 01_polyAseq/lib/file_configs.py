import os, sys
from lib.util.lib_utils import cmd_exists
from ConfigParser import SafeConfigParser

class ProgramConfig:
	def __init__(self,
							 config_file,
							 nnode=1,
							 ncpu=1,
							 verbose=False,
							 resume=-1,
							 debug=False):

		self.fpath = None
		self.prog_rootdir = None
		self.sections = {}
		self.set_config_file(config_file)

		self.nnode = nnode
		self.ncpu = ncpu
		self.verbose = verbose
		self.resume = resume
		self.debug = debug

	def set_config_file(self, config_file):

		if not os.path.exists(config_file):
			raise IOError('the config file [%s] does not exist! refer to README file' % config_file)

		self.prog_rootdir = os.path.dirname(os.path.dirname(os.path.abspath(config_file)))
		self.fpath = config_file

	def complete_filepath(self,fpath):
		fpath2 = fpath.split()
		if len(fpath2)>2:
			raise RuntimeError('check the file configuration [%s] at [%s]'%(fpath,self.fapth))
		fpath1 = fpath2[-1]
		if os.path.sep in fpath1:
			if not fpath1.startswith(os.path.sep) or not fpath1.startswith("~"):
				fpath1 = os.path.join(self.prog_rootdir, fpath1)
		if len(fpath2)==2:
			fpath1 = "%s %s"%(fpath2[0],fpath1)
		return fpath1

	def load_items(self):
		if not self.fpath:
			self.set_config_file()

		parser = SafeConfigParser()
		parser.read(self.fpath)

		self.sections = {}

		for section_name in parser.sections():
			self.sections[section_name] = {}
			for name, value in parser.items(section_name):
				value = value.strip()
				if name.startswith('runopt'):
					if 'runopts' not in self.sections[section_name]:
						self.sections[section_name]['runopts']=[]

					option_vars = [var.strip() for var in value.split(' ')]
					L = len(option_vars)
					if L > 2:
						raise AttributeError('define a single run option at each line!')

					if name.endswith('_path'):
						if L != 2:
							raise AttributeError('this line contains runopt and path but the # of items after split is not 2!')
						else:
							option_vars[1] = self.complete_filepath(option_vars[1])

					self.sections[section_name]['runopts'].append(option_vars)
				elif name.endswith('_path'):
					self.sections[section_name][name] = self.complete_filepath(value)
				else:
					self.sections[section_name][name] = value

				if self.verbose:
					print('section_name[%s],name[%s],value[%s]'%(section_name,name,value))

	def required(self,required_sections):
		for section in required_sections:
			if section in self.sections:
				bin_path_found = False
				for name, value in self.sections[section].iteritems():
					path_found = False
					if name.endswith('_prefix'):
						path_found = False
						index_name = value.split(os.path.sep)[-1]
						for fn in os.listdir(os.path.dirname(value)):
							if fn.startswith(index_name):
								path_found = True
								break
					elif name.endswith('_path'):
						if os.path.exists(value):
							path_found = True
						elif cmd_exists(value):
							path_found = True
						else:
							raise EnvironmentError('check section[%s],name[%s],value[%s] in %s'%(section,name,value,self.fpath))
					if name == 'bin_path':
						if path_found:
							bin_path_found = True
				if bin_path_found:
					if self.verbose:
						print('found a valid bin_path at section[%s]'%(section))
				else:
					raise EnvironmentError('name[bin_path] in section[%s] is not valid. Check [%s].'%(section, self.fpath))
			else:
				raise EnvironmentError('section[%s] is required in [%s]' % (section,self.fpath))

if __name__ == '__main__':
	#testing section
	print('testing ...')
	cP = ProgramConfig(config_file = '../configs/program.conf',verbose=True)
	cP.load_items()
	cP.required(['bowtie2','bowtie2-build','bwa','demux','star'])
	print('Done.')
