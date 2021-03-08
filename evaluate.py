'''python evaluate.py  -o data/test_ascot_performance.h5
'''
from mmsplice import MTSplice
from os.path import isfile
from kipoi.writers import HDF5BatchWriter
from mtsplice.utils.functional import logit
from scipy.special import expit
import h5py
import argparse

from keras.models import load_model
import numpy as np
from mtsplice.data.ascot_dl import Ascot, tissues
from mtsplice import get_data_dir
DATA = get_data_dir()

# command line args
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output",
                    dest='output',
                    help="the output file")
args = parser.parse_args()

if isfile(args.output):
    print("File exist:")

spline_dl = Ascot(ascot=DATA+"ascot_gtex_test.csv.gz",
                  fasta_file=DATA+"hg38.fa",
                  pad_trim_same_l=False, mean_inpute=False,
                  region_anno=False, length=400, flanking=300,
                  seq_align='both', encode=True, flanking_exons=False)

model = MTSplice()  # predict on logit scale
diter = spline_dl.batch_pred_iter(shuffle=False)

pred = []
measured = []
mean = []
x_hat = []

for x, y in diter:
    batch = {
        "acceptor": x[0],
        "donor": x[1]
    }
    _pred = model.predict_on_batch(batch)
    x_hat.append(_pred)
    pred.append(expit(_pred + x[2]))
    measured.append(y)
    mean.append(x[2])

x_hat = np.concatenate(x_hat)
pred = np.concatenate(pred)
measured = np.concatenate(measured)
mean = np.concatenate(mean)
x = logit(measured) - np.nanmean(logit(measured), axis=-1, keepdims=True)

HDF5BatchWriter.dump(args.output,
                     {'measured': measured,
                      'x_hat': x_hat,
                      'x': x,
                      'pred': pred})
