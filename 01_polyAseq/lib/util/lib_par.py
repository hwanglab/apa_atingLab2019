from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
import time
from lib.util.lib_utils import run_syscmd

def par_shell_runs(arg_sets,
									 ncpu=1,
									 job_name='bash_shell'):
	if ncpu < 1:
		raise RuntimeError("Number of processors to use must be greater than 0")

	C = int(cpu_count() * .75) #check if the number of cpu requested is larger than the max in the node
	if ncpu > C:
		ncpu = C

	if ncpu < 1:
		ncpu = 1

	# Start my pool
	pool = ThreadPool(ncpu)
	print("Running %d [%s] jobs using %d processors..." % (len(arg_sets), job_name, ncpu))

	start_time = time.time()
	# Run tasks
	results = []
	for arg_tuple in arg_sets:
		results.append(pool.apply_async(run_syscmd, arg_tuple))

	print "Waiting for result..."

	# Close the pool and wait for each running task to complete
	pool.close()
	pool.join()

	# Process results
	retcodes = []
	for result in results:
		retcode, out = result.get()
		print("out: {} err: {}".format(out, retcode))
		retcodes.append(retcode)

	print "Elapsed time: %.2f secs" % (time.time() - start_time)
	return retcodes

def find_error_index_from_retcodes(retcodes):
	error_idx = []
	for i,retcode in enumerate(retcodes):
		if retcode != 0:
			error_idx.append(i)
	return error_idx