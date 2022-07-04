"""

Author: Janith Weerman
Date: 08/06/2022

Wrapper script for carrying out network propagation

"""

import os


# Create directories to hold output files
if not os.path.exists("/outputs"):
    os.makedirs('/outputs')
if not os.path.exists("outputs/graphs"):
    os.makedirs('outputs/graphs')
if not os.path.exists("outputs/results"):
    os.makedirs('outputs/results')

