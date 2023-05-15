#!/usr/bin/env python
import subprocess

# Specify the directory path
directory_path = './temp'

# Run MultiQC command on the directory
multiqc_cmd = f"multiqc {directory_path}"
subprocess.run(multiqc_cmd, shell=True)
