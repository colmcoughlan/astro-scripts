#                                                       LOFAR PIPELINE FRAMEWORK
#
#                                                           Run a remote command
#                                                            John Swinbank, 2010
#                                                      swinbank@transientskp.org
# ------------------------------------------------------------------------------

from __future__ import with_statement
from collections import defaultdict
from threading import BoundedSemaphore

import re
import os
import signal
import threading
import time
import xml.dom.minidom as xml

from lofarpipe.support.pipelinelogging import log_process_output
from lofarpipe.support.utilities import spawn_process
from lofarpipe.support.lofarexceptions import PipelineQuit
from lofarpipe.support.jobserver import job_server
import lofarpipe.support.lofaringredient as ingredient
from lofarpipe.support.xmllogging import add_child

# By default, Linux allocates lots more memory than we need(?) for a new stack
# frame. When multiplexing lots of threads, that will cause memory issues.
threading.stack_size(1048576)

class ParamikoWrapper(object):
    """
    Sends an SSH command to a host using paramiko, then emulates a Popen-like
    interface so that we can pass it back to pipeline recipes.
    """
    def __init__(self, paramiko_client, command):
        self.returncode = None
        self.client = paramiko_client
        self.chan = paramiko_client.get_transport().open_session()
        self.chan.get_pty()
        self.chan.exec_command(command)
        self.stdout = self.chan.makefile('rb', -1)
        self.stderr = self.chan.makefile_stderr('rb', -1)

    def communicate(self):
        if not self.returncode:
            self.returncode = self.chan.recv_exit_status()
        stdout = "\n".join(line.strip() for line in self.stdout.readlines()) + "\n"
        stderr = "\n".join(line.strip() for line in self.stdout.readlines()) + "\n"
        return stdout, stderr

    def poll(self):
        if not self.returncode and self.chan.exit_status_ready():
            self.returncode = self.chan.recv_exit_status()
        return self.returncode

    def wait(self):
        if not self.returncode:
            self.returncode = self.chan.recv_exit_status()
        return self.returncode

    def kill(self):
        self.chan.close()

def run_remote_command(config, logger, host, command, env, arguments = None):
    """
    Run command on host, passing it arguments from the arguments list and
    exporting key/value pairs from env(a dictionary).

    Returns an object with poll() and communicate() methods, similar to
    subprocess.Popen.

    This is a generic interface to potentially multiple ways of running
    commands (SSH, mpirun, etc). The appropriate method is chosen from the
    config block supplied (with SSH as a fallback).
    """
    try:
        method = config.get('remote', 'method')
    except:
        method = None

    logger.info("********************** Remote method is %s" % method)

    if method == "paramiko":
        try:
            key_filename = config.get('remote', 'key_filename')
        except:
            key_filename = None
        return run_via_paramiko(logger, host, command, env, arguments, key_filename)
    elif method == "mpirun":
        return run_via_mpirun(logger, host, command, env, arguments)
    elif method == "mpirun_fionn":
        return run_via_mpirun_fionn(logger, host, command, env, arguments)
    elif method == "local":
        return run_via_local(logger, command, arguments)
    elif method == "juropa_mpi":
        return run_via_mpiexec(logger, command, arguments, host)
    elif method == "cep_mpi":
        return run_via_mpiexec_cep(logger, command, arguments, host)
    elif method == "slurm_srun_cep3":
        return run_via_slurm_srun_cep3(logger, command, arguments, host)
    elif method == "custom_cmdline":
        return run_via_custom_cmdline(logger, host, command, env, arguments, config)
    else:
        return run_via_ssh(logger, host, command, env, arguments)

def run_via_slurm_srun_cep3(logger, command, arguments, host):
    logger.debug("Dispatching command to %s with srun" % host)
    for arg in arguments:
        command = command + " " + str(arg)
    commandstring = ["srun","-N 1","-n 1","-w",host, "/bin/sh", "-c", "hostname && " + command]
    #commandstring = ["srun","-N 1","--cpu_bind=map_cpu:none","-w",host, "/bin/sh", "-c", "hostname && " + command]
    # we have a bug that crashes jobs when too many get startet at the same time
    # temporary NOT 100% reliable workaround
    #from random import randint
    #time.sleep(randint(0,10))
    ##########################
    process = spawn_process(commandstring, logger)
    process.kill = lambda : os.kill(process.pid, signal.SIGKILL)
    return process

def run_via_mpirun_fionn(logger, host, command, environment, arguments):
    """
    Dispatch a remote command via mpirun.

    Return a Popen object pointing at the MPI command, to which we add a kill
    method for shutting down the connection if required.
    """
    logger.debug("Dispatching command to %s with mpirun" % host)
    mpi_cmd = ["mpirun", "-hosts", host, "-np", "1"]
    envlst = ''
    for key in environment.keys():
        envlst = envlst + ',' + str(key)# remember to remove first comma)
    mpi_cmd.extend(['-envlist',envlst[1:]])
    mpi_cmd.extend(command.split())  # command is split into (python, script)
    mpi_cmd.extend(str(arg) for arg in arguments)
    # print("MPI command NEW = "+str(mpi_cmd))
    env = os.environ
    env.update(environment)
    process = spawn_process(mpi_cmd, logger, env = env)
    # mpirun should be killed with a SIGTERM to enable it to shut down the
    # remote command.
    process.kill = lambda : os.kill(process.pid, signal.SIGTERM)
    return process

def run_via_mpirun(logger, host, command, environment, arguments):
    """
    Dispatch a remote command via mpirun.

    Return a Popen object pointing at the MPI command, to which we add a kill
    method for shutting down the connection if required.
    """
    logger.debug("Dispatching command to %s with mpirun" % host)
    mpi_cmd = ["/usr/bin/mpirun", "-host", host]
    for key in environment.keys():
        mpi_cmd.extend(["-x", key])
    mpi_cmd.append("--")
    mpi_cmd.extend(command.split())  # command is split into (python, script)
    mpi_cmd.extend(str(arg) for arg in arguments)
    env = os.environ
    env.update(environment)
    process = spawn_process(mpi_cmd, logger, env = env)
    # mpirun should be killed with a SIGTERM to enable it to shut down the
    # remote command.
    process.kill = lambda : os.kill(process.pid, signal.SIGTERM)
    return process

# let the mpi demon manage free resources to start jobs
def run_via_mpiexec(logger, command, arguments, host):
    for arg in arguments:
        command = command + " " + str(arg)
    commandstring = ["mpiexec", "-x", "-np=1", "/bin/sh", "-c", "hostname && " + command]
    process = spawn_process(commandstring, logger)
    process.kill = lambda : os.kill(process.pid, signal.SIGKILL)
    return process

# start mpi run on cep
# TODO: rsync fails on missing ssh key??
def run_via_mpiexec_cep(logger, command, arguments, host):
    for arg in arguments:
        command = command + " " + str(arg)
    commandstring = ["mpiexec", "-x", "PYTHONPATH", "-x", "LD_LIBRARY_PATH", "-x", "PATH", "-H", host, "/bin/sh", "-c", "hostname ; " + command]
    process = spawn_process(commandstring, logger)
    process.kill = lambda : os.kill(process.pid, signal.SIGKILL)
    return process


def run_via_local(logger, command, arguments):
    commandstring = ["/bin/sh", "-c"]
    for arg in arguments:
        command = command + " " + str(arg)
    commandstring.append(command)
    process = spawn_process(commandstring, logger)
    process.kill = lambda : os.kill(process.pid, signal.SIGKILL)
    return process

def run_via_ssh(logger, host, command, environment, arguments):
    """
    Dispatch a remote command via SSH.

    We return a Popen object pointing at the SSH session, to which we add a
    kill method for shutting down the connection if required.
    """
    logger.debug("Dispatching command to %s with ssh" % host)
    ssh_cmd = ["ssh", "-n", "-tt", "-x", host, "--", "/bin/sh", "-c"]

    commandstring = ["%s=%s" % (key, value) for key, value in environment.items()]
    commandstring.append(command)
    commandstring.extend(re.escape(str(arg)) for arg in arguments)
    ssh_cmd.append('"' + " ".join(commandstring) + '"')
    process = spawn_process(ssh_cmd, logger)
    process.kill = lambda : os.kill(process.pid, signal.SIGKILL)
    return process

def run_via_custom_cmdline(logger, host, command, environment, arguments, config):
    """
    Dispatch a remote command via a customisable command line

    We return a Popen object pointing at the running executable, to which we add a
    kill method for shutting down the connection if required.

    The command line is taken from the "remote.cmdline" configuration option,
    with the following strings replaced:

      {host}         := host to execute command on
      {command}      := bash command line to be executed
      {uid}          := uid of the calling user

      {image}        := docker.image configuration option
      {slurm_job_id} := the SLURM job id to allocate resources in

    """
    commandArray = ["%s=%s" % (key, value) for key, value in environment.items()]
    commandArray.append(command)
    commandArray.extend(re.escape(str(arg)) for arg in arguments)
    commandStr = " ".join(commandArray)

    try:
        image = config.get('docker', 'image')
    except:
        image = "lofar"

    # Construct the full command line, except for {command}, as that itself
    # can contain spaces which we don't want to split on.
    full_command_line = config.get('remote', 'cmdline').format(
      uid          = os.geteuid(),
      slurm_job_id = os.environ.get("SLURM_JOB_ID"),
      docker_image = image,
      host         = host,
      command      = "{command}"
    ).split(' ')

    # Fill in {command} somewhere
    full_command_line = [x.format(command = commandStr) for x in full_command_line]

    logger.debug("Dispatching command to %s with custom_cmdline: %s" % (host, full_command_line))

    process = spawn_process(full_command_line, logger)
    process.kill = lambda : os.kill(process.pid, signal.SIGKILL)
    return process

def run_via_paramiko(logger, host, command, environment, arguments, key_filename):
    """
    Dispatch a remote command via paramiko.

    We return an instance of ParamikoWrapper.
    """
    logger.debug("Dispatching command to %s with paramiko" % host)
    import paramiko
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host, key_filename = key_filename)
    commandstring = ["%s=%s" % (key, value) for key, value in environment.items()]
    commandstring.append(command)
    commandstring.extend(re.escape(str(arg)) for arg in arguments)
    return ParamikoWrapper(client, " ".join(commandstring))

class ProcessLimiter(defaultdict):
    """
    Provide a dictionary-like structure of bounded semaphores with arbitrary
    keys.

    This gives a convenient way to keep tabs on the number of simultaneous
    jobs running on a given host.

    :param nproc: Bound value for semaphore (ie, maximum number of jobs)
    :type nproc: integer or none
    """
    def __init__(self, nproc = None):
        if nproc:
            super(ProcessLimiter, self).__init__(
                lambda: BoundedSemaphore(int(nproc))
            )
        else:
            class Unlimited(object):
                """
                Dummy semaphore for the unlimited case.
                Acquire and release always succeed.
                """
                def acquire(self):
                    return True
                def release(self):
                    return True
            super(ProcessLimiter, self).__init__(Unlimited)

class ComputeJob(object):
    """
    Container for information about a job to be dispatched to a compute node.

    :param host: Target host for job
    :param command: Full path to command to be run on target host
    :param arguments: List of arguments which will be passed to command
    """
    def __init__(self, host, command, arguments = []):
        self.host = host
        self.command = command
        self.arguments = arguments
        self.results = {}
        self.results['returncode'] = 123456  # Default to obscure code to allow
        # test of failing ssh connections 
        # TODO: This could be done nicer!!! THis error is shown in the logfile
        # and tends to confuse users

    def dispatch(self, logger, config, limiter, id, jobhost, jobport,
                  error, killswitch):

        """
        Dispatch this job to the relevant compute node.

        Note that error is an instance of threading.Event, which will be set
        if the remote job fails for some reason.
        """
        self.id = id
        limiter[self.host].acquire()
        # Start the time after we aquire the lock!
        time_info_start = time.time()
        try:
            if killswitch.isSet():
                logger.debug("Shutdown in progress: not starting remote job")
                self.results['returncode'] = 1
                error.set()
                return 1
            process = run_remote_command(
                config,
                logger,
                self.host,
                self.command,
                {
                    "PATH": os.environ.get('PATH'),
                    "PYTHONPATH": os.environ.get('PYTHONPATH'),
                    "LD_LIBRARY_PATH": os.environ.get('LD_LIBRARY_PATH'),
                    "LOFARROOT" : os.environ.get('LOFARROOT'),
                    "QUEUE_PREFIX" : os.environ.get('QUEUE_PREFIX','')
                },
                arguments = [id, jobhost, jobport]
            )
            # Wait for process to finish. In the meantime, if the killswitch
            # is set (by an exception in the main thread), forcibly kill our
            # job off.
            while process.poll() == None:
                if killswitch.isSet():
                    process.kill()
                else:
                    time.sleep(1)
            sout, serr = process.communicate()

            serr = serr.replace("Connection to %s closed.\r\n" % self.host, "")
            log_process_output("Remote command", sout, serr, logger)
        except Exception, e:
            logger.exception("Failed to run remote process %s (%s)" % (self.command, str(e)))
            self.results['returncode'] = 1
            error.set()
            return 1
        finally:
            limiter[self.host].release()

        if process.returncode != 0:
            logger.error(
                "Remote process %s %s failed on %s (status: %d)" % \
                (self.command, self.arguments, self.host, process.returncode)
            )
            error.set()

        # after node returned.
        # add the duration of
        time_info_end = time.time()
        self.results["job_duration"] = str(time_info_end - time_info_start)
        self.results['returncode'] = process.returncode

        logger.debug(
            "compute.dispatch results job {0}: {1}: {2}, {3}: {4} ".format(
              self.id, "job_duration", self.results["job_duration"],
                     "returncode", self.results["returncode"] ))
        return process.returncode


def threadwatcher(threadpool, logger, killswitch):
    """
    Start and watch a pool of threads. If an exception is thrown during
    processing, set the killswitch so that all threads can shut down cleanly,
    then join all the threads to wait for them to finish.

    :param threadpool: Pool of threads to handle
    :param logger: Logger
    :type logger: logging.Logger or descendant
    :param killswitch: Indication for threads to abort
    :type killswitch: threading.Event
    """
    # If we receive a SIGTERM, shut down processing.
    signal.signal(signal.SIGTERM, killswitch.set)
    try:
        # Start all the threads, but don't just join them, as that
        # blocks all exceptions in the main thread. Instead, we wake
        # up every second to handle exceptions.
        [thread.start() for thread in threadpool]
        logger.info("Waiting for compute threads...")

        while True in [thread.isAlive() for thread in threadpool]:
            time.sleep(1)
    except:
        # If something throws an exception (normally a
        # KeyboardException, ctrl-c) set the kill switch to tell the
        # comput threads to terminate, then wait for them.
        logger.warn("Processing interrupted: shutting down")
        killswitch.set()
    finally:
        # Always make sure everything has finished. Note that if an exception
        # is thrown before all the threads have started, they will not all be
        # alive (and hence not join()-able).
        [thread.join() for thread in threadpool if thread.isAlive()]


class RemoteCommandRecipeMixIn(object):
    """
    Mix-in for recipes to dispatch jobs using the remote command mechanism.
    """
    def _schedule_jobs(self, jobs, max_per_node = None):
        """
        Schedule a series of compute jobs. Blocks until completion.

        :param jobs: iterable of :class:`~lofarpipe.support.remotecommand.ComputeJob` to be scheduled
        :param max_per_node: maximum number of simultaneous jobs on any given node
        :type max_per_node: integer or none
        :rtype: dict mapping integer job id to :class:`~lofarpipe.support.remotecommand.ComputeJob`
        """
        threadpool = []
        jobpool = {}
        if not max_per_node and self.config.has_option('remote', 'max_per_node'):
            max_per_node = self.config.getint('remote', 'max_per_node')
        limiter = ProcessLimiter(max_per_node)
        killswitch = threading.Event()

        if max_per_node:
            self.logger.info("Limiting to %d simultaneous jobs/node" % max_per_node)
        '''
        print('Scheduling jobs:')
        for job in jobs:
                print('job command = '+str(job.command))
                print('job host = '+str(job.host))
                print('job arguements = '+str(job.arguments))
        '''
        with job_server(self.logger, jobpool, self.error) as (jobhost, jobport):
            self.logger.debug("Job dispatcher at %s:%d" % (jobhost, jobport))
            for job_id, job in enumerate(jobs):
                jobpool[job_id] = job
                threadpool.append(
                    threading.Thread(
                        target = job.dispatch,
                        args = (
                            self.logger, self.config, limiter, job_id,
                            jobhost, jobport, self.error, killswitch
                        )
                    )
                )
            threadwatcher(threadpool, self.logger, killswitch)

        if killswitch.isSet():
            raise PipelineQuit()

        # Add information regarding specific nodes to an xml node.
        self.logger.debug("Adding node_logging_information")
        local_document = xml.Document()
        node_durations = local_document.createElement("nodes")
        for job_id, job in enumerate(jobs):
            # Test if the duration is there
            # fixme the name of node_durations is not logical
            if "job_duration" in job.results:
                child_node_duration = add_child(node_durations, "job")
                child_node_duration.setAttribute("job_id", str(job_id))
                child_node_duration.setAttribute("job_host", str(job.host))
                child_node_duration.setAttribute("duration",
                     str(job.results["job_duration"]))

                # return code if present (Not there on error)
                if "returncode" in job.results:
                    child_node_duration.setAttribute(
                        "returncode", str(job.results['returncode']))
                else:
                    child_node_duration.setAttribute(
                        "returncode", str(-1))

                ## If there is 'node level' resource logging available
                if "monitor_stats" in job.results:
                      return_node = xml.parseString(
                          job.results['monitor_stats']).documentElement

                      child_node_duration.appendChild(return_node)


        # manually add the result xml as an ingredient output.
        # this allows backward compatible logging: If not read an additional
        # output does not matter
        self.outputs._fields["return_xml"] = ingredient.StringField(
                                                help = "XML return data.")
        self.outputs["return_xml"] = node_durations.toxml(encoding = "ascii")

        return jobpool
