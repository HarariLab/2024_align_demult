import pandas as pd
import os

POOLS = pd.read_csv(config['libraries'])
SAMPLES = pd.read_csv(config['pools_samples'])