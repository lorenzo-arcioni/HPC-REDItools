import sys
import os
from optparse import OptionParser
import subprocess as sp

inputpath  = None
binary     = None
reference  = None
wlm        = None
threads    = 1

#Initialize all options and give a velue for all global variables
def init():
    """
    Initializes global variables based on command line input arguments.

    Global Variables:
        threads (int): The number of threads to use.
        inputpath (string): The input BAM file path.
        binary (string): The binary file path.
        reference (string): The reference file path.
        wlm (string): The WLM path.

    Returns:
        None
    """
    
    # Declare global variables
    global threads, inputpath, binary, reference, wlm

    # Define the options
    parser = OptionParser()

    # Add options to the parser
    parser.add_option(
        "-i", 
        "--input-bam", 
        dest="input_path", 
        help="Input BAM file directory path"
    )
    parser.add_option(
        "-r", 
        "--reference", 
        dest="rf", 
        help="Reference file path"
    )
    parser.add_option(
        "-w", 
        "--wlm", 
        dest="wlm", 
        help="WLM path"
    )
    parser.add_option(
        "-t", 
        "--threads", 
        dest="threads", 
        help="Number of threads to use"
    )

    # Parse the command line arguments
    (options, args) = parser.parse_args()

    # Set global variables based on command line arguments
    inputpath = options.input_path or inputpath
    reference = options.rf or reference
    wlm = options.wlm or wlm
    threads = options.threads or threads

def fill_startbase(threads, inputpath, reference):
    """
    Writes a start script that includes configuration information about the 
    input job and environment

    Returns:
        None
    """

    # Create the start script file
    with open("./start.sh", "w") as start_file:

        # Write the bash script header
        start_file.write("#!/bin/bash\n\n")

        # Write the appropriate start script template based on the workload manager (WLM)
        if wlm == 'htcondor':
            pass
        elif wlm == 'slurm':
            # Load the slurm start script template
            with open("./bases/slurm_start_base.txt", "r") as slurm_file:
                start_file.write(slurm_file.read())
                slurm_file.close()
        else:
            pass

        # Load the main start script template
        with open("./bases/start_base.txt", "r") as base_file:
            base_template = base_file.read()
            # Fill the template with the appropriate values
            base_template = base_template.format(inputpath, threads, reference)
            # Write the filled template to the start script file
            start_file.write('\n' + base_template)
            base_file.close()

    start_file.close()

def fill_readbase(threads, wlm, reference, inputpath):
    """Generates a file with the necessary parameters to execute the read.py script."""

    # Open read.py file to write
    with open("./read.py", "w") as read:
        
        # Open reditools_additional_options.txt file to read and remove new line characters
        with open("./bases/reditools_additional_options.txt", "r") as f:
            rtao = f.read().replace("\n", "")
            f.close()
        
        # Check the workload manager used
        if wlm == 'slurm':
            # Open slurm_partial_script_base.txt file to read and replace new line characters with escaped characters
            with open("./bases/slurm_partial_script_base.txt", "r") as f:
                header = f.read().replace("\n", "\\n")
                f.close()
        elif wlm == 'htcondor':
             pass
        # If not using slurm, set header as empty string
        else:
            header = ""
        
        # Open read_base.txt file to read
        with open("./bases/read_base.txt", "r") as b:
            # Read the file contents and format it with the necessary parameters
            base = b.read().format(inputpath, reference, threads, wlm, header, rtao)
            # Write the formatted contents to the read.py file
            read.write(base)

        # Close the read.py file
        read.close()

def fill_controlscriptbase(wlm):
	"""
	Fills the control script base according to the Workload Manager (WLM) the user chose.
	"""

	# Writes the shebang to the control script file
	with open("./control_script.sh", "w") as control:
		control.write("#!/bin/bash\n\n")

		# If the WLM is htcondor, do nothing
		if wlm == 'htcondor':
			pass

		# If the WLM is slurm, fill the control script base with the slurm base and sbatch
		elif wlm == 'slurm':

			# Reads the slurm control script base file and writes its content to the control script file
			with open("./bases/slurm_controlscript_base.txt", "r") as f:
				control.write(f.read())
				f.close()

			# Reads the control script base file and replaces the placeholder with sbatch, then writes to the control script file
			with open("./bases/controlscript_base.txt", "r") as f:
				base = f.read().format("sbatch")
				control.write(base)
				f.close()

		# If the WLM is not htcondor or slurm, fill the control script base with bash
		else:

			# Reads the control script base file and replaces the placeholder with bash, then writes to the control script file
			with open("./bases/controlscript_base.txt", "r") as f:
				base = f.read().format("bash")
				control.write(base)
				f.close()

		control.close()  # Closes the control script file

def fill_monitor(wlm):
    """
    Generates a monitor script based on the workload manager (WLM) being used.
    If no WLM is being used, the monitor script is deleted.
    """

    with open("monitor.sh", "w") as monitor:
        monitor.write("#!/bin/bash\n")

        # If using the SLURM workload manager, add a header and print job information
        if wlm == "slurm":
            monitor.write("# Print job information\n")
            monitor.write("echo \"             JOBID            PARTITION                           NAME     USER    STATE       TIME TIME_LIMI  NODES NODELIST(REASON)\"\n")
            monitor.write("squeue --format=\"%.18i %.20P %.30j %.8u %.8T %.10M %.9l %.6D %R\" --me | grep \"PA_proc-\"\n")

        # If using the HTCondor workload manager, do nothing
        elif wlm == "htcondor":
            pass

        monitor.close()

    # If no WLM is being used, remove the monitor script
    if wlm == "none":
        os.remove("monitor.sh")

def fill_cancel(wlm: str) -> None:
    """
    Create a cancel script to cancel all processes of a computation.

    Args:
        wlm: A string representing the workload manager being used.
            Valid options are 'slurm' or 'htcondor', or 'none' if not using a workload manager.
    """
    
    with open("cancel.sh", "w") as cancel:
        cancel.write("#!/bin/bash\n")
        cancel.write("#Cancel all Processes of the computation\n")

        if wlm == "slurm":
            cancel.write("squeue --me | grep \"PA_proc-\" | awk '{print $1}' | xargs -n 1 scancel\n")
            cancel.write("echo Computation aborted >> ./general.log\n")
            cancel.write("echo '-------------------------------------------------' >> ./general.log")

        elif wlm == "htcondor":
            # This section is empty as there are no commands to cancel processes in HTCondor.
            pass

        cancel.close()

    if wlm == "none":
        os.remove("cancel.sh")

def main():
    """
    Main function that fills the bases and sets the necessary file permissions.
    """
    
    global threads   # number of threads
    global inputpath # input file
    global wlm       # Workload Manager
    global reference # reference file
    
    fill_startbase(threads, inputpath, reference)
    fill_readbase(threads, wlm, reference, inputpath)
    fill_controlscriptbase(wlm)
    # fill_monitor()
    # fill_cancel()

    sp.call("chmod 777 start.sh", shell=True)
    sp.call("chmod 777 control_script.sh", shell=True)
    sp.call("chmod 777 read.py", shell=True)
    
if __name__ == '__main__':
    init()
    sys.exit(main())
